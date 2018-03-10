=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSeedRepeatPipeline;

use strict;
use warnings;
use feature 'say';

use POSIX;
use File::Spec::Functions;
use Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;
  my $assembly_registry_dba = new Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor(%{$self->param_required('assembly_registry_db')});
  $self->hrdb_set_con($assembly_registry_dba,'assembly_registry_db');
}


sub run {
  my $self = shift;

  my $assembly_registry_dba = $self->hrdb_get_con('assembly_registry_db');
  my $min_contig_n50        = $self->param('min_contig_n50');
  my $min_scaffold_n50      = $self->param('min_scaffold_n50');
  my $min_total_length      = $self->param('min_total_length');
  my $max_version_only      = $self->param('max_version_only');
  say "Fetching GCAs by constraints...";
  my $inital_gca_list       = $assembly_registry_dba->fetch_gca_by_constraints($min_contig_n50,$min_scaffold_n50,$min_total_length,undef,$max_version_only);
  say "Processing GCAs to determine what to process...";
  my $gcas_to_process       = $self->find_gcas_to_process($inital_gca_list);
  $self->gcas_to_process($gcas_to_process);
}


sub write_output {
  my $self = shift;

  my $gcas_to_process = $self->gcas_to_process();

  say "Preparing to flow output ids for downstream processing...";
  # First output the new jobs
  foreach my $gca (@{$gcas_to_process}) {
    my $processing_hash->{'iid'} = $gca;
    $self->dataflow_output_id($processing_hash,2);
  }

}


sub post_cleanup {
  my ($self) = @_;

  # Sleep and output a new job for this analysis
  my $sleep_length_hours = $self->param('sleep_length_hours');
  unless($sleep_length_hours) {
    $sleep_length_hours = 1;
  }

  ceil($sleep_length_hours);
  say "Going to sleep for ".$sleep_length_hours." hours...";
  sleep($sleep_length_hours * 600);

  say "Seeding new run id for this analysis...";
  my $run_id = $self->param('iid');
  $run_id++;

  my $repeat_job_hash->{'iid'} = $run_id;
  $self->dataflow_output_id($repeat_job_hash,3);
}


sub find_gcas_to_process {
  my ($self,$inital_gca_list) = @_;

  my $gcas_to_process = [];
  my $base_repeat_dir = $self->param_required('base_repeat_dir');
  my $gca_dir = catfile($base_repeat_dir,'GCA');
  my $species_dir = catfile($base_repeat_dir,'species');
  foreach my $gca (@{$inital_gca_list}) {
    $gca =~ /^(GCA\_\d{9})\.(\d+)$/;
    my $chain = $1;
    my $version = $2;

    my $full_gca_path = catfile($gca_dir,$chain,$version);
    unless(-e $full_gca_path) {
      my $result = system('mkdir -p '.$full_gca_path);
      if($result) {
        $self->throw("Could not write output dir to path: ".$full_gca_path);
      }
      push(@{$gcas_to_process},$gca);
    }
  } # end foreach my $gca

  return($gcas_to_process);
}


sub gcas_to_process {
  my ($self,$gcas_to_process) = @_;

  if($gcas_to_process) {
    $self->param('_gcas_to_process',$gcas_to_process);
  }

  return($self->param('_gcas_to_process'));
}

1;
