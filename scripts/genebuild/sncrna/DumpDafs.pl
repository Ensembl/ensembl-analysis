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


# print "#chr\tstart\tstop\tcoord\tscore\tstrand\thitname\tevalue\tpid\tcigar\n";

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
