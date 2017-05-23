use warnings;
use strict;
use feature 'say';


use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::Utilities;
use Bio::EnsEMBL::Utils::Exception;

my $host;
my $port=3306;
my $dbname;
my $user;
my $pass;
my $dnahost;
my $dnaport;
my $dnadbname;
my $dnauser;

GetOptions( 'dbhost:s'        => \$host,
            'dbport:n'        => \$port,
            'dbname:s'        => \$dbname,
            'dbuser:s'        => \$user,
            'dbpass:s'        => \$pass,
            'dnadbhost:s'     => \$dnahost,
            'dnadbport:s'     => \$dnaport,
            'dnadbuser:s'     => \$dnauser,
            'dnadbname:s'     => \$dnadbname);

my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $host,
                                              -user   => $user,
                                              -port   => $port,
                                              -dbname => $dbname,
                                              -pass   => $pass, );



my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                                -port    => $dnaport,
                                                -user    => $dnauser,
                                                -host    => $dnahost,
                                                -dbname  => $dnadbname);

$dba->dnadb($dnadb);

my $gene_adaptor = $dba->get_GeneAdaptor;
my $slice_adaptor = $dba->get_SliceAdaptor;
my $slices = $slice_adaptor->fetch_all('toplevel');

foreach my $slice (@{$slices}) {
  my $genes = $gene_adaptor->fetch_all_by_Slice($slice);
  my $gene_strings;

  foreach my $gene (@{$genes}) {
    my $gene_string = "";# $gene->start.":".$gene->end.":".$gene->seq_region_name;
    my $transcripts = $gene->get_all_Transcripts();
    my $transcript_strings;
    foreach my $transcript (@{$transcripts}) {
      my $transcript_string = $transcript->start.":".$transcript->end.":".$transcript->seq_region_name;
      my $exons = $transcript->get_all_Exons();
      my $exon_string = generate_exon_string($exons);
      $transcript_string .= ":".$exon_string;
      if($transcript->translation) {
        $transcript_string .= ":".$transcript->translation->seq;
      }

      if($transcript_strings->{$transcript_string}) {
        say "Found a duplicate transcript within a gene:";
        say "Duplicate transcript id: ".$transcript->dbID;
        say "Duplicate transcript start: ".$transcript->start;
        say "Duplicate transcript end: ".$transcript->end;
        say "Duplicate strand: ".$transcript->strand;
        say "Duplicate name: ".$transcript->seq_region_name;
      } else {
        $transcript_strings->{$transcript_string} = 1;
      }

      $gene_string .= $transcript_string
    } # end foreach my $transcript

    if($gene_strings->{$gene_string}) {
      say "Found a duplicate gene:";
      say "Duplicate id: ".$gene->dbID;
      say "Duplicate start: ".$gene->start;
      say "Duplicate end: ".$gene->end;
      say "Duplicate strand: ".$gene->strand;
      say "Duplicate name: ".$gene->seq_region_name;

      eval {
        say "Deleting duplicate: ".$gene->dbID();
        $gene_adaptor->remove($gene);
      }; if($@){
        say "Couldn't remove gene ".$gene->dbID.":\n($@)\n";
      }

    } else {
      $gene_strings->{$gene_string} = 1;
    }

  } # end foreach my $gene

} # end foreach my $slice

exit;


sub generate_exon_string {
  my ($exon_array) = @_;

  my $exon_string = "";
  foreach my $exon (@{$exon_array}) {
    my $start = $exon->start();
    my $end = $exon->end();
    $exon_string .= $start."..".$end.":";
  }

  return($exon_string);
}

