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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMinimapAlign

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMinimapAlign;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable::Minimap;

=head2 fetch_input

 Arg [1]    : None
 Description: Create the Minimap runnable to align the isoseq sequences in 'isoseq_file'
              to the genome in 'genome_file' using minimap2. It will create a sorted BAM
              file which name will be "'isoseq_file'.bam"
 Returntype : None
 Exceptions : Throws if 'genome_file' does not exist
              Throws if 'isoseq_file' does not exist

=cut

sub fetch_input {
  my ($self) = @_;

  my $genome_file = $self->param_required('genome_file');
  $self->throw("Could not find the genome file $genome_file") unless (-e $genome_file);
  my $isoseq_file = $self->param_required('isoseq_file');
  $self->throw("Could not find the genome file $isoseq_file") unless (-e $isoseq_file);
  $self->create_analysis;
  my $options;
  if ($self->param_is_defined('commandline_params')) {
    $options = $self->param('commandline_params');
  }
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Minimap->new(
    -analysis => $self->analysis,
    -targetfile => $genome_file,
  );
  $runnable->queryfile($isoseq_file);
  if ($self->param_is_defined('samtools')) {
    $runnable->samtools($self->param('samtools'));
  }
  $self->runnable($runnable);
}

=head2 write_output

 Arg [1]    : None
 Description: The BAM file is written during run, we only dataflow the name of the file on branch
              '_branch_to_flow_to'
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;

  $self->dataflow_output_id({filename => $self->output->[0]}, $self->param('_branch_to_flow_to'));
}

1;
