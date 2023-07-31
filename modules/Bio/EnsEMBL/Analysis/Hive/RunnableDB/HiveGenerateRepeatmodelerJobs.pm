=head1 LICENSE

Copyright [2018-2022] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadIMGT

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenerateRepeatmodelerJobs;

use strict;
use warnings;
use feature 'say';

use Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor;
use Bio::EnsEMBL::Hive::Utils qw(destringify);
use POSIX;
use File::Spec::Functions;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults()},
    repeatmodeler_job_count => 5,
  }
}


sub fetch_input {
  my $self = shift;

}


sub run {
  my $self = shift;

  my $gca = $self->param_required('iid');
  my $path_to_genomic_fasta = $self->param_required('path_to_genomic_fasta');
  my $run_count = $self->param('run_count');
  $self->create_output_subdirs($path_to_genomic_fasta,$run_count);
}


sub write_output {
  my $self = shift;

  foreach my $run_dir (@{$self->output()}) {
    my $job_params = destringify($self->input_job->input_id());
    $job_params->{'repeat_run_dir'} = $run_dir;
    $self->dataflow_output_id($job_params,2);
  }
}

sub create_output_subdirs {
  my ($self,$path_to_genomic_fasta,$run_count) = @_;

  my $master_output = catfile($path_to_genomic_fasta,'output');
  if(system('mkdir -p '.$master_output)) {
    $self->throw("Could not make main output dir on path: ".$master_output);
  }

  for(my $i=0; $i<$run_count; $i++) {
    my $run_subdir = catfile($master_output,$i);
    if(system('mkdir -p '.$run_subdir)) {
      $self->throw("Could not make main output dir on path: ".$master_output);
    }
    push(@{$self->output},$run_subdir);
  }
}

sub path_to_genomic_fasta {
  my ($self,$path) = @_;
  if($path) {
    $self->param('_path_to_genomic_fasta',$path);
  }

  return($self->param('_path_to_genomic_fasta'));
}

1;
