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

ok( 1, 'Startup test' );

my $human = Bio::EnsEMBL::Test::MultiTestDB->new('mus_musculus');
my $human_dba = $human->get_DBAdaptor('empty');
ok( 2, 'Emptyness is available') or BAIL_OUT 'Cannot setup Empty core DB. Do not continue.';

my $multi_db = Bio::EnsEMBL::Test::MultiTestDB->new('hive');
my $hive_dba = $multi_db->get_DBAdaptor('hive') or BAIL_OUT 'Cannot get HIVE DB. Do not continue';

{
  my $module = 'download_assembly_init_conf' ;
  my $pipeline = Bio::EnsEMBL::Test::RunPipeline->new( $module, $options );
  $pipeline->run();
  #TODO file existance content checks

  my $module = 'download_assembly_continue_conf' ;
  my $pipeline = Bio::EnsEMBL::Test::RunPipeline->new( $module, $options );
  $pipeline->run();

}


done_testing();
