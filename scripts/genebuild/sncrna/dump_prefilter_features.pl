#!/usr/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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


	print FH $daf->seq_region_name(), "\t",
		$daf->seq_region_start(), "\t",
		$daf->seq_region_end(), "\t",
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

## dump repeat features
#my $rfa = $db->get_RepeatFeatureAdaptor();
#$fn = $working_dir . "/repeats.bed";
#open(FH, '>', $fn) or die "Could not write to $fn";
#
#my @repeats = @{ $rfa->fetch_all() };
#
#foreach my $repeat (@repeats){
#  my $strand = $repeat->strand() > 0 ? "+" : "-";
#  print FH $repeat->seq_region_name(), "\t",
#    $repeat->seq_region_start(), "\t",
#    $repeat->seq_region_end(), "\t",
#    $strand, "\n";
#}
#
#close(FH);

# dump putative stem-loops
my $gene_adaptor = $db->get_GeneAdaptor();
my @genes = @{ $gene_adaptor->fetch_all_by_biotype('miRNA')};

$fn = $working_dir . "/identified_mirnas.bed";

open(FH, '>', $fn) or die "Could not write to $fn";

foreach my $gene (@genes){
    my $strand = $gene->strand() > 0 ? "+" : "-";


      print FH $gene->seq_region_name(), "\t",
          $gene->seq_region_start(), "\t",
          $gene->seq_region_end(), "\t",
          $gene->seq_region_name(), ":",
          $gene->seq_region_start(), "-",
          $gene->seq_region_end(), "\t0\t",
          $strand, "\t",
          $gene->dbID(), "\n";

}

close(FH);

