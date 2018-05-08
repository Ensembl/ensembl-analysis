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

