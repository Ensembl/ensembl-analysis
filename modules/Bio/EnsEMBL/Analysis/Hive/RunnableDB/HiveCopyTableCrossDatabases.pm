=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#Copyright [2016-2022] EMBL-European Bioinformatics Institute
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


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCopyTableCrossDatabases - 

=head1 SYNOPSIS


=head1 DESCRIPTION

Copy one table's data to another

=head1 METHODS

=cut
package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCopyTableCrossDatabases;


use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(run_command);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');



=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCopyTableCrossDatabases
  Function  :  
  file
  Returntype: 
  Exceptions: 
  Example   : 

=cut



use Bio::EnsEMBL::DBSQL::DBAdaptor;


sub param_defaults {
    return {
      dump_translations_path => '~/enscode/ensembl-analysis/scripts/protein/',
      dump_translations_name => 'dump_translations.pl',

      # dump_translations.pl script parameters
      # logic => '', # to set the genes and transcripts analysis (logic names) (optional)
      sourcehost => undef,
      sourceuser => undef,
      sourceport => '3306',
      sourcepass => undef,
      sourcedbname => undef,
      targethost => undef,
      targetdbname => undef,
      targetuser => undef,
      targetport => '3306',
      file => undef
    }
}

sub fetch_input {
  my $self = shift;

  return 1;
}

sub run {

  my $self = shift;
  # lincRNA_output_db genebuild12
  # kb15_papio_anubis_core_77_2_linc genebuild11
  $self->param_required('sourcehost');
  $self->param_required('sourceuser');
  $self->param_required('sourceport');
  $self->param_required('sourcepass');
  $self->param_required('sourcedbname');
  $self->param_required('targethost');
  $self->param_required('targetdbname');
  $self->param_required('targetuser');
  $self->param_required('targetpass');  
  $self->param_required('targetport');


  my $cmd_analysis_id = "mysql --host " . $self->param('sourcehost') . " --port 3306 " . 
                " --user " . $self->param('sourceuser')  . 
                " --pass=" . $self->param('sourcepass')  .
                "  " . $self->param('sourcedbname')  . 
                "  -NB -e  \"select analysis_id from analysis where logic_name = \'pfam\';\" ";
  print "xxxxxxxx \n";
  my $analysis_id = `$cmd_analysis_id`; 

  chomp($analysis_id); 
  my $cmd_cp_table   =  "mysqldump --opt -uensro -h " .  $self->param('sourcehost') . "  " . $self->param('sourcedbname') . 
                        "  --tables protein_feature    --where=\" analysis_id = " .  $analysis_id . "\" --no-create-info  |  mysql -h " . 
                        $self->param_required('targethost') . " -P 3306 -u " .  $self->param('targetuser')  . " --password=" . 
                        $self->param_required('targetpass') . " -D " .  $self->param('targetdbname')   ;  
  print "-->>>Command2: " . $cmd_cp_table . "\n"; 

  print run_command($cmd_cp_table,"Copying tables...");

  # throw("test stop!!\n"); 

  return 1;
}

sub write_output {
  my $self = shift;
  return 1;
}

1;



#################### how the config should look like: ####################
 #      {
 #       -logic_name => 'HiveCopyTableCross_test_pi',
 #       -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCopyTableCrossDatabases',
 #       -parameters => {
 #                              # dump_translations.pl script parameters
 #                              sourcehost => $self->o('reference_db','-host'),
 #                              sourceuser => $self->o('reference_db','-user'),
 #                              sourceport => $self->o('reference_db','-port'),
 #                              sourcepass => $self->o('reference_db','-pass'),
 #                              sourcedbname => $self->o('reference_db','-dbname'),
 #                              targethost => $self->o('lincRNA_output_db','-host'),
 #                              targetdbname => $self->o('lincRNA_output_db','-dbname'),
 #                              targetuser => $self->o('lincRNA_output_db','-user'),
 #                              targetpass => $self->o('lincRNA_output_db','-pass'),
 #                              targetport => $self->o('lincRNA_output_db','-port'),
 #                            },
 #       -rc_name    => 'default',
 #       -input_ids => [{}],
 #       -wait_for   => ['StoreFeatures_test_pi'], 

        # -flow_into     => ['SplitDumpFiles_test_pi'],
 #      },



