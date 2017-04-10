#!/usr/bin/env perl

# Copyright [2017] EMBL-European Bioinformatics Institute
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

=head1 NAME

  create_gencode_web_data_patch.pl

=head1 DESCRIPTION

A script to generate the production db patch for the GENCODE web_data/analysis_web data

=head1 OPTIONS

  --prod_host           The host the production db resides on

  --prod_port           The port for the production db

  --prod_dbname         The name of the production db

  --prod_user           The user for the production db

  --species             The species to generate the patch for

  --gencode_version     The new GENCODE version

=head1 EXAMPLE

  perl create_gencode_web_data_patch.pl --species mus_musculus --gencode_version 14

=cut

use warnings;
use strict;
use feature 'say';
use Getopt::Long;

use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $prod_host = 'mysql-ens-sta-1';
my $prod_port = 4519;
my $prod_user = 'ensro';
my $prod_dbname = 'ensembl_production';
my $species;
my $gencode_version;
my $previous_gencode_version;

GetOptions( 'prod_host:s'       => \$prod_host,
            'prod_port:n'       => \$prod_port,
            'prod_user:s'       => \$prod_user,
            'prod_dbname:s'     => \$prod_dbname,
            'species:s'         => \$species,
            'gencode_version:n' => \$gencode_version);

unless($species eq 'mus_musculus' || $species eq 'homo_sapiens') {
  throw("You must enter either 'homo_sapiens' or 'mus_musculus' at as the species, e.g. -species mus_musculus");
}

unless($gencode_version && $gencode_version !~ /\D/) {
  throw("You must specify the new GENCODE version number (do not include the M in the mouse version), e.g. -gencode_version 14");
}

$previous_gencode_version = $gencode_version - 1;
unless($previous_gencode_version > 0) {
  throw("The previous GENCODE version is calculated as being <= 0, this is wrong. New GENCODE version passed in to script: ".$gencode_version);
}

my $production_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                    -host   => $prod_host,
                    -user   => $prod_user,
                    -port   => $prod_port,
                    -dbname => $prod_dbname);


my $select_species_id = "SELECT species_id FROM analysis_web_data WHERE species_id=(select species_id from species where production_name=?)";

my $sth = $production_db->dbc->prepare($select_species_id);
$sth->bind_param(1, $species);
$sth->execute();
my ($result_row) = $sth->fetchrow_array();
unless($result_row) {
  throw("Failed to find the species id from the species table in the production db. The table was queried with the following:\n".
        "SELECT species_id FROM species WHERE production_name = '".$species."';");
}

my $species_id = $result_row;
my $select_previous_web_data_id;
my $new_web_data;
my $new_comment;
if($species eq 'mus_musculus') {
  $new_web_data = '{"caption" => "Genes (Comprehensive set from GENCODE M'.$gencode_version.')","colour_key" => "[biotype]","default" => '.
                  '{"MultiBottom" => "collapsed_label","MultiTop" => "gene_label","alignsliceviewbottom" => "as_collapsed_label"'.
                  ',"contigviewbottom" => "transcript_label","contigviewtop" => "gene_label","cytoview" => "gene_label"},'.
                  '"key" => "ensembl","label_key" => "[biotype]","multi_name" => "GENCODE M'.$gencode_version.' Comprehensive gene set",'.
                  '"name" => "Comprehensive Gene Annotations from GENCODE M'.$gencode_version.'"}';
  $new_comment = 'Mouse GENCODE M'.$gencode_version;
  $select_previous_web_data_id = 'SELECT web_data_id FROM web_data WHERE data LIKE "%GENCODE M'.$previous_gencode_version.' Comprehensive gene set%"';
} else {
  $new_web_data = '{"caption" => "Genes (Comprehensive set from GENCODE '.$gencode_version.')","colour_key" => "[biotype]","default" =>'.
                  ' {"MultiBottom" => "collapsed_label","MultiTop" => "gene_label","alignsliceviewbottom" => "as_collapsed_label",'.
                  '"contigviewbottom" => "transcript_label","contigviewtop" => "gene_label","cytoview" => "gene_label"},"key" => '.
                  '"ensembl","label_key" => "[biotype]","multi_name" => "GENCODE '.$gencode_version.' Comprehensive gene set",'.
                  '"name" => "Comprehensive Gene Annotations from GENCODE '.$gencode_version.'"}';
  $new_comment = 'Use this for human GENCODE';
  $select_previous_web_data_id = 'SELECT web_data_id FROM web_data WHERE data LIKE "%GENCODE '.$previous_gencode_version.' Comprehensive gene set%"';
}


$sth = $production_db->dbc->prepare($select_previous_web_data_id);
$sth->execute();
($result_row) = $sth->fetchrow_array();
unless($result_row) {
  throw("Failed to find the web data id from the web table in the production db for the previous GENCODE version.".
        " The table was queried with the following:\n".$select_previous_web_data_id);
}

my $previous_web_data_id = $result_row;
my $insert_new_web_data = 'INSERT into web_data (data,comment,created_at) values(\''.$new_web_data.'\',\''.$new_comment.'\',now());';
my $update_web_data_id;
my $update_modified_at;
if($species eq 'mus_musculus') {
  $update_web_data_id = 'UPDATE analysis_web_data set web_data_id=(select web_data_id from web_data where data like "%GENCODE M'.$gencode_version.
                        ' Comprehensive gene set%") where web_data_id='.$previous_web_data_id.';';
  $update_modified_at = 'UPDATE analysis_web_data set modified_at=now() where web_data_id=(select web_data_id from web_data where data like "%GENCODE M'.
                         $gencode_version.' Comprehensive gene set%");';
} else {
  $update_web_data_id = 'UPDATE analysis_web_data set web_data_id=(select web_data_id from web_data where data like "%GENCODE '.$gencode_version.
                        ' Comprehensive gene set%") where web_data_id='.$previous_web_data_id.';';
  $update_modified_at = 'UPDATE analysis_web_data set modified_at=now() where web_data_id=(select web_data_id from web_data where data like "%GENCODE '.
                         $gencode_version.' Comprehensive gene set%");';
}

say $insert_new_web_data;
say $update_web_data_id;
say $update_modified_at;

exit;
