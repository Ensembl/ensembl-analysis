=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMetaPipelineRun

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMetaPipelineRun;

use strict;
use warnings;

use File::Basename;
use File::Spec::Functions qw(catfile);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(execute_with_wait);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    action => 'loop', # this can be loop, run, 
    _input_id_name => 'target_db',
    perl_path => 'perl',
    beekeeper_script => 'beekeeper.pl',
    tweak_script => 'tweak_pipeline.pl',
  }
}

sub fetch_input {
  my ($self) = @_;

  my $url = $self->param_required('ehive_url');
  my $base_cmd = $self->param_required('beekeeper_script').' -url '.$url;
  if (-e $self->param('beekeeper_script')) {
    $base_cmd = $self->param('perl_path').' '.$base_cmd;
  }
  execute_with_wait($base_cmd.' -sync');
  $self->param('base_cmd', $base_cmd);
}

sub run {
  my ($self) = @_;

  my $cmd = $self->param('base_cmd').' -'.$self->param('action');
  $cmd .= ' '.$self->param('commandline_params')
    if ($self->param_is_defined('commandline_params') and $self->param('commandline_params'));
  execute_with_wait($cmd);
  execute_with_wait($self->param('base_cmd').' -sync');
  execute_with_wait($cmd);
}


sub write_output {
  my ($self) = @_;

  $self->dataflow_output_id({$self->param('_input_id_name') => $self->param($self->param('_input_id_name'))}, $self->param('_branch_to_flow_to'));
}

sub pre_cleanup {
  my ($self) = @_;

# I need to execute fetch_input to set base_cmd and it's also doing some checks
  $self->fetch_input;
  if ($self->param_is_defined('reset_excluded')) {
    my $reset_cmd = $self->param('tweak_script');
    if (-e $self->param('tweak_script')) {
      $reset_cmd = $self->param('perl_path')." $reset_cmd";
    }
    else {
      if (-e $self->param('beekeeper_script')) {
        my $directory = dirname($self->param('beekeeper_script'));
        $reset_cmd = $self->param('perl_path').' '.catfile($directory, $reset_cmd);
      }
    }
    execute_with_wait("$reset_cmd -url ".$self->param_required('ehive_url')." -SET 'analysis[".$self->param_required('reset_excluded')."].is_excluded=0'");
  }
  execute_with_wait($self->param('base_cmd').' -unkwn');
  execute_with_wait($self->param('base_cmd').' -reset_failed_jobs');
}

1;
