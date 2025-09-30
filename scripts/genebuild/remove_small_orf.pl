use warnings;
use strict;
use feature 'say';


use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long qw(:config no_ignore_case);

my $dbname = '';
my $user_w   = '';
my $user_r   = '';
my $host   = '';
my $port   = '';
my $pass   = '';


my $dna_dbname = '';
my $dna_host   = '';
my $dna_port   = '';

my $small_orf_cutoff = 100;
my $small_intron = 75;
my $flag;
my $help;

my $options = GetOptions ("user_w|dbuser|u=s"  => \$user_w,
                          "host|dbhost|h=s"  => \$host,
                          "port|dbport|P=i"  => \$port,
                          "dbname|db|D=s"    => \$dbname,
                          "dbpass|pass|p=s"  => \$pass,
                          "dna_dbname=s"      => \$dna_dbname,
			  "user_r|u_r=s"  => \$user_r,
                          "dna_host|dna_h=s"  => \$dna_host,
                          "dna_port|dna_P=i"  => \$dna_port,
                          "orf_cutoff|o=i"   => \$small_orf_cutoff,
                          "intron_cutoff|i=i"   => \$small_intron,
                          "flag|f"           => \$flag,
                          "help|h"           => \$help,);

if ($help){
    print("\nRemove or Flag small ORFs\n\nUSAGE:\nperl flag_small_orf.pl -user ensadmin -host <host> -port <port> -dbname <db_name> -dnahost <dna_host> -dnaport <dna_port> -dnadbname <dnadb_name> [options]\n\nOPTIONS:\n-f\tFlag small ORFs (do not remove them, assign 'small_orf' biotype)\n-o n\tSmall ORF cutoff (default=100)\n-i n\tSmall intron cutoff (default=75)\n-h\tShow this help and exit\n\n");
    exit;
  }

my $dna_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $dna_port,
  -user    => $user_r,
  -host    => $dna_host,
  -dbname  => $dna_dbname);

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $port,
  -user    => $user_w,
  -host    => $host,
  -dbname  => $dbname,
  -pass    => $pass);

$db->dnadb($dna_db);

my $gene_adaptor = $db->get_GeneAdaptor();
my $genes = $gene_adaptor->fetch_all();

#say "Looping through genes...";
foreach my $gene (@$genes) {
  unless($gene->biotype eq 'protein_coding') {
    next;
  }

  my $single_exon_only = 1;
  my $transcripts = $gene->get_all_Transcripts;

  foreach my $transcript (@{$transcripts}) {
    my $introns = $transcript->get_all_Introns;
    if(scalar(@$introns) > 1) {
      $single_exon_only = 0;
      last;
    }

    my $intron = shift(@$introns);
    if($intron && $intron->length >= $small_intron) {
      $single_exon_only = 0;
      last;
    }

    if($transcript->translation->length >= $small_orf_cutoff) {
      $single_exon_only = 0;
      last;
    }
  }

  if($single_exon_only) {
    if ($flag){
      say "Flagging gene ".$gene->dbID();
      $gene->set_Biotype('small_orf');
      $gene_adaptor->update($gene);
    }
    else {
      say "Removing gene ".$gene->dbID();
      $gene_adaptor->remove($gene);
    }
  }
}

exit;
