=head1 LICENSE

# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Runnable::SamtoolsMerge -

=head1 SYNOPSIS

  Do NOT instantiate this class directly: must be instantiated
  from a subclass (see ExonerateTranscript, for instance).

=head1 DESCRIPTION


=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Analysis::Runnable::SamtoolsMerge;

use warnings;
use strict;

use parent ('Bio::EnsEMBL::Analysis::Runnable::BaseBamMerge');


sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    $self->samtools($self->program);

    return $self;
}



############################################################
#
# Analysis methods
#
############################################################

=head2 run

Usage   :   $obj->run($workdir, $args)
Function:   Runs exonerate script and puts the results into the file $self->results
            It calls $self->parse_results, and results are stored in $self->output
=cut

sub run {
  my ($self) = @_;

  my $input_files = $self->input_files;
  $input_files = $input_files->[0] if (scalar(@{$input_files}) == 1);

  $self->samtools->merge($self->options, $self->output_file, $input_files);
  $self->check_output_file;

  return 1;
}


1;
