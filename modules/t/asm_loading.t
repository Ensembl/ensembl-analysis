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

use File::Copy;

use Test::More;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::RunPipeline;

use Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveDownloadContigs;

my $options = "" ;

#corresponds to 'output_path' in the pipeline config files
my $pipe_output_path = "/home/rishi/eat_files" ;

ok( 1, 'Startup test' );

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new('hive');
my $hive_dba = $multi_db->get_DBAdaptor('hive') or BAIL_OUT 'Cannot get HIVE DB. Stopping.';

#TO DO - read these values out from the config file directly (or perhaps via versa)
my $db_name = 'eat_rat_core_db' ;
my $db_server = 'localhost' ;
my $db_user   = 'eat_tester' ;
my $db_password = 'apwd1234' ;
my $db_port = '3306' ;

{
  my $module = 'download_assembly_init_conf' ;
  my $pipeline = Bio::EnsEMBL::Test::RunPipeline->new( $module, $options );
  $pipeline->run();
  my $dir = $pipe_output_path.'/Primary_Assembly' ;
  ok(-d $dir, "Primary Assembly directory exists") or done_testing, exit;
  ok(-e $dir."/assembly_report.txt", "Assembly report exists") ;
  ok(-e $dir."/AGP/chr20.agp", "chr20.agp exists") or done_testing, exit ;
  ok(-e $dir."/contigs/contigs.fa", "contigs.fa exists") or done_testing, exit ;

  $module = 'download_assembly_continue_conf' ;
  $pipeline = Bio::EnsEMBL::Test::RunPipeline->new( $module, $options );
  $pipeline->run();
  #open database adaptor, do checks and close
  my $dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -user   => $db_user,
    -pass   => $db_password,
    -dbname => $db_name,
    -host   => $db_server,
    -driver => 'mysql'
  );

  #load contig value checks
  my $sth_select = $dba->dbc->prepare('SELECT COUNT(*) FROM dna');
  $sth_select->execute();
  my ($checkval) = $sth_select->fetchrow_array;
  is($checkval, 1570, "num dna sequences loaded") ;
  $sth_select = $dba->dbc->prepare('SELECT COUNT(*) FROM seq_region');
  $sth_select->execute();
  my ($checkval) = $sth_select->fetchrow_array;
  is($checkval, 1594, "num seq regions loaded") ;

  #load_assembly_info
  my $sth_select = $dba->dbc->prepare('SELECT COUNT(*) FROM assembly');
  $sth_select->execute();
  my ($checkval) = $sth_select->fetchrow_array;
  is($checkval, 3125, "assembly info loaded") ;

}

done_testing();
