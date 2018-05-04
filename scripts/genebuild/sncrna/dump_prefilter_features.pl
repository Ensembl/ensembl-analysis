use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

my ($dbname, $dbhost, $dbport, $dbuser, $working_dir, $logic_name) = @ARGV;

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-DBNAME => $dbname,
  	-HOST => $dbhost,
  	-PORT => $dbport,
  	-USER => $dbuser,
	-DRIVER => 'mysql',
);

my $daf_adaptor = $db->get_DnaAlignFeatureAdaptor();

my @dafs = @{ $daf_adaptor->fetch_all_by_logic_name($logic_name)};

my $fn = $working_dir . "/" . $logic_name . "_dafs.bed";

open(FH, '>', $fn) or die "Could not write to $fn";

foreach my $daf (@dafs){
	my $strand = $daf->strand() > 0 ? "+" : "-";


	print FH "chr" . $daf->seq_region_name(), "\t",
		$daf->seq_region_start(), "\t",
		$daf->seq_region_end(), "\tchr",
		$daf->seq_region_name(), ":",
		$daf->seq_region_start(), "-",
		$daf->seq_region_end(), "\t",
		$daf->score(), "\t",
		$strand, "\t",
		$daf->hseqname(), "\t",
		$daf->p_value(), "\t",
		$daf->percent_id(), "\t",
		$daf->cigar_string(),  "\n";

}

close(FH);

# dump repeat features
my $rfa = $db->get_RepeatFeatureAdaptor();
$fn = $working_dir . "/repeats.bed";
open(FH, '>', $fn) or die "Could not write to $fn";

my @repeats = @{ $rfa->fetch_all() };

foreach my $repeat (@repeats){
  my $strand = $repeat->strand() > 0 ? "+" : "-";
  print FH "chr" . $repeat->seq_region_name(), "\t",
    $repeat->seq_region_start(), "\t",
    $repeat->seq_region_end(), "\t",
    $strand, "\n";
}

close(FH);

# dump putative stem-loops
my $gene_adaptor = $db->get_GeneAdaptor();
my @genes = @{ $gene_adaptor->fetch_all_by_biotype('miRNA')};

$fn = $working_dir . "/identified_mirnas.bed";

open(FH, '>', $fn) or die "Could not write to $fn";

foreach my $gene (@genes){
    my $strand = $gene->strand() > 0 ? "+" : "-";


      print FH "chr" . $gene->seq_region_name(), "\t",
          $gene->seq_region_start(), "\t",
          $gene->seq_region_end(), "\tchr",
          $gene->seq_region_name(), ":",
          $gene->seq_region_start(), "-",
          $gene->seq_region_end(), "\t0\t",
          $strand, "\t",
          $gene->dbID(), "\n";

}

close(FH);

