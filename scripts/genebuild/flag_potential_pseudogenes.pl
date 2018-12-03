use warnings;
use strict;
use feature 'say';

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);
use Getopt::Long qw(:config no_ignore_case);

my $dbname = '';
my $user_r = 'ensro';
my $host   = '';
my $port   = '';

my $dna_dbname = '';
my $dna_host   = '';
my $dna_port   = '';

my $user_w = '';
my $pass   = '';

my $options = GetOptions ("user_r=s"        => \$user_r,
                          "user_w=s"        => \$user_w,
                          "host|dbhost|h=s" => \$host,
                          "port|dbport|P=i" => \$port,
                          "dbname|db|D=s"   => \$dbname,
                          "dbpass|pass|p=s" => \$pass,
                          "dna_dbname=s"    => \$dna_dbname,
                          "dna_host=s"      => \$dna_host,
                          "dna_port=i"      => \$dna_port);

my $dna_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $dna_port,
  -user    => $user_r,
  -host    => $dna_host,
  -dbname  => $dna_dbname);

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $port,
  -user    => $user_w,
  -pass    => $pass,
  -host    => $host,
  -dbname  => $dbname);

$db->dnadb($dna_db);

my $genes = $db->get_GeneAdaptor()->fetch_all;
my $potential_pseudogenes = flag_potential_pseudogenes($genes);
write_output($potential_pseudogenes,$db);

exit;


sub write_output {
  my ($output_genes,$output_db) = @_;

  my $gene_adaptor = $output_db->get_GeneAdaptor();
  foreach my $output_gene (@{$output_genes}) {
    $gene_adaptor->update($output_gene);
  }
}


sub flag_potential_pseudogenes {
  my ($genes) = @_;
  my $output_genes = [];

  my $frameshift_intron_length = 75;
  my $max_frameshifts = 2;
  my $max_non_canonical_ratio = 0.1;

  # First pass implementation, just takes the most supported transcripts in each position and then
  # looks to see if there are small introns. These should be passed to the pseudogene module
  foreach my $gene (@$genes) {
    my $transcripts = $gene->get_all_Transcripts;
    if(scalar(@$transcripts) > 1) {
      throw("More than one transcript per gene is not supported. Despite what the loop below this code would suggest");
    }


    my $flag = 0;
    my $non_canonical_flag = 0;
    foreach my $transcript (@$transcripts) {
      say "Transcript: ".$transcript->dbID;
      my $introns = $transcript->get_all_Introns;
      unless(scalar(@$introns) >= 1) {
        next;
      }

      my $small_intron_counter = 0;
      my $non_canonical_count = 0;

      # Add a penalty for not starting with a met
      unless($transcript->translation->seq =~ /^M/) {
        $small_intron_counter++;
      }

      my $total_length = 0;
      foreach my $intron (@$introns) {
        my $intron_length = $intron->length;
        if($intron_length < $frameshift_intron_length) {
          $small_intron_counter++;
        } elsif(!$intron->is_splice_canonical) {
          $non_canonical_count++;
        }
        $total_length += $intron_length;
      }

      my $avg_length = $total_length / scalar(@$introns);
      if(($avg_length < $frameshift_intron_length) ||
         (scalar(@$introns) == 1 && $small_intron_counter == 1) ||
         ($small_intron_counter > $max_frameshifts)) {
         $flag = 1;
         $transcript->biotype($transcript->biotype()."_pseudo");
      } elsif(($non_canonical_count / scalar(@$introns)) > $max_non_canonical_ratio) {
        $non_canonical_flag = 1;
        $transcript->biotype($transcript->biotype()."_noncanon");
      }
    }

    if($flag) {
      say "FLAGGING FRAMESHIFT: ".$gene->dbID;
      $gene->biotype($gene->biotype()."_pseudo");
      push(@{$output_genes},$gene);
    } elsif($non_canonical_flag) {
      say "FLAGGING NON-CANON: ".$gene->dbID;
      $gene->biotype($gene->biotype()."_noncanon");
      push(@{$output_genes},$gene);
    }
  }

  return($output_genes);
}
