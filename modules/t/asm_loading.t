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

use Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveDownloadContigs;

my $options;
if (@ARGV) {
  $options = join( q{ }, @ARGV );
}
else {
  $options = '-run_all 1';
}

#corresponds to 'output_path' in the pipeline config files
my $pipe_output_path = "/home/rishi/eat_files" ;

ok( 1, 'Startup test' );

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new('hive');
my $hive_dba = $multi_db->get_DBAdaptor('hive') or BAIL_OUT 'Cannot get HIVE DB. Stopping.';

my $rat_multi_db = Bio::EnsEMBL::Test::MultiTestDB->new('rattus_norvegicus');
my $rat_empty_dba = $rat_multi_db->get_DBAdaptor('empty') or BAIL_OUT 'Cannot setup Empty core DB. Stopping.';


{
  my $module = 'download_assembly_init_conf' ;
  my $pipeline = Bio::EnsEMBL::Test::RunPipeline->new( $module, $options );
  $pipeline->run();
  my $dir = $pipe_output_path.'/Primary_Assembly' ;
  ok(-d $dir, "Primary Assembly directory exists") or done_testing, exit;
  ok(-e $dir."/assembly_report.txt", "Assembly report exists") ;
  ok(-e $dir."/chr20.agp", "chr20.agp exists") or done_testing, exit ;
  ok(-e $dir."/contigs.fa", "contigs.fa exists") or done_testing, exit ;

  $module = 'download_assembly_continue_conf' ;
  $pipeline = Bio::EnsEMBL::Test::RunPipeline->new( $module, $options );
  $pipeline->run();
  #TODO database value existance checks
}



done_testing();
