# Copyright [1999-2016] EMBL-European Bioinformatics Institute
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

use Test::More;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::RunPipeline;

ok( 1, 'Startup test' );

my $options = "" ;

ok( 1, 'Startup test' );

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new('hive');
my $hive_dba = $multi_db->get_DBAdaptor('hive') or BAIL_OUT 'Cannot get HIVE DB. Stopping.';

#TO DO - read these values out from the config file directly (or perhaps via versa)
my $db_name = 'eat_rat_core_db' ;
my $db_server = 'localhost' ;
my $db_user   = '_EAT_DB_USER_' ;
my $db_password = '_EAT_DB_PASS_' ;
my $db_port = '3306' ;

{
  my $module = 'download_assembly_init_conf' ;
  my $pipeline = Bio::EnsEMBL::Test::RunPipeline->new( $module, $options );
  $pipeline->run();

  #open database adaptor for checks
  my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -user   => $db_user,
    -pass   => $db_password,
    -dbname => $db_name,
    -host   => $db_server,
    -driver => 'mysql'
  );

  #repeat tests - don't make them fixed as programs and databases may change on the farm
  # outside our control
  # TODO update this example
  my $sth_select = $dba->dbc->prepare('SELECT COUNT(*) FROM dna');
  $sth_select->execute();
  my ($checkval) = $sth_select->fetchrow_array;
  is($checkval, 443, "num dna sequences loaded") ;
  $sth_select = $dba->dbc->prepare('SELECT COUNT(*) FROM seq_region');
  $sth_select->execute();
  my ($checkval) = $sth_select->fetchrow_array;
  is($checkval, 463, "num seq regions loaded") ;


}


done_testing();