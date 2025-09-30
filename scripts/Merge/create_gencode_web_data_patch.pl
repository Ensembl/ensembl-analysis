#!/usr/bin/env perl

# Copyright [2017-2024] EMBL-European Bioinformatics Institute
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

my $host = 'mysql-ens-vertannot-staging';
my $port = 4573;
my $user = 'ensro';
my $dbname;
my $prod_host = 'mysql-ens-sta-1';
my $prod_port = 4519;
my $prod_user = 'ensro';
my $prod_dbname = 'ensembl_production';
my $species;
my $gencode_version;
my $user_prod_id;

GetOptions( 'host:s'            => \$host,
            'port:n'            => \$port,
            'user:s'            => \$user,
            'dbname:s'          => \$dbname,
            'prod_host:s'       => \$prod_host,
            'prod_port:n'       => \$prod_port,
            'prod_user:s'       => \$prod_user,
            'prod_dbname:s'     => \$prod_dbname,
            'species:s'         => \$species,
            'user_prod_id:s'    => \$user_prod_id,
            'gencode_version:n' => \$gencode_version);

unless($species eq 'mus_musculus' || $species eq 'homo_sapiens') {
  throw("You must enter either 'homo_sapiens' or 'mus_musculus' at as the species, e.g. -species mus_musculus");
}

if ($gencode_version && $gencode_version =~ /\D/) {
  throw("You must specify the new GENCODE version number (do not include the M in the mouse version), e.g. -gencode_version 14");
}

throw("You must specify your web data id so we can track who makes changes") unless ($user_prod_id && $user_prod_id !~ /\D/);

my $production_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                    -host   => $prod_host,
                    -user   => $prod_user,
                    -port   => $prod_port,
                    -dbname => $prod_dbname);


my $query_data = 'SELECT wd.web_data_id, wd.data FROM species s LEFT JOIN analysis_web_data awd ON s.species_id = awd.species_id LEFT JOIN analysis_description ad ON awd.analysis_description_id = ad.analysis_description_id LEFT JOIN web_data wd ON wd.web_data_id = awd.web_data_id WHERE ad.logic_name = "ensembl" AND awd.db_type = "core" AND s.production_name = ?';

my $sth = $production_db->dbc->prepare($query_data);
$sth->bind_param(1, $species);
$sth->execute();
my @data = $sth->fetchrow_array();
my ($web_data_id, $data) = @data;
@data = $sth->fetchrow_array();
throw('More than one result for your query') unless (scalar(@data) == 0);
throw("Failed to find the species id from the species table in the production db. The table was queried with the following:\n".
      "SELECT species_id FROM species WHERE production_name = '".$species."';")
  unless($web_data_id);

#if($species eq 'mus_musculus') {
#  throw('This does not look like the web_data wanted: '.$web_data_id."\n$data\n") unless ($data =~ /GENCODE M(\d+)/);
#  $update_web_data = "GENCODE M%d', 'GENCODE M%d";
#} else {
#  throw('This does not look like the web_data wanted: '.$web_data_id."\n$data\n") unless ($data =~ /GENCODE (\d+)/);
#  $update_web_data = "GENCODE %d', 'GENCODE %d";
#}
throw('This does not look like the web_data wanted: '.$web_data_id."\n$data\n") unless ($data =~ /GENCODE([ M]+)(\d+)/);
my $update_web_data = "UPDATE web_data SET data = REPLACE(data, 'GENCODE$1%d', 'GENCODE$1%d'), modified_at = NOW(), modified_by = %d WHERE web_data_id = %d;\n";
my $previous_gencode_version = $2;
if (!$gencode_version) {
  if (!$dbname) {
    $sth = $production_db->dbc->prepare("SELECT db_type, db_release, db_assembly FROM db JOIN species USING(species_id) WHERE db.db_type = 'core' AND db.is_current = 1 AND production_name = '$species'");
    $sth->execute;
    @data = $sth->fetchrow_array;
    $dbname = join('_', $species, @data);
    @data = $sth->fetchrow_array;
    throw('Something is wrong with meta_key "gencode.version"') unless (scalar(@data) == 0);
  }
  my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                    -host   => $host,
                    -user   => $user,
                    -port   => $port,
                    -dbname => $dbname);
  $sth = $db->dbc->prepare('SELECT meta_value FROM meta WHERE meta_key = "gencode.version"');
  $sth->execute;
  @data = $sth->fetchrow_array;
  ($gencode_version) = $data[0] =~ /(\d+)/;;
  @data = $sth->fetchrow_array;
  throw('Something is wrong with meta_key "gencode.version"') unless (scalar(@data) == 0);
}
$sth->finish;
throw("Versions do not follow $previous_gencode_version -> $gencode_version") unless ($previous_gencode_version+1 == $gencode_version);

say "use $prod_dbname;";
printf $update_web_data, $previous_gencode_version, $gencode_version, $user_prod_id, $web_data_id;
