=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::MergeBamFiles -

=head1 SYNOPSIS

  [mergebamfiles]
  module=MergeBamFiles

=head1 DESCRIPTION

  Merge the bam files to create the "super" bam file needed by
  Bam2Genes, the rough model step.
  It uses samtools to check but it uses picard for merging,
  no config file as everything is fetched from the analysis table
  The module is looking for files named: *_sorted.bam

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMergeBamFiles;

use warnings ;
use strict;

use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub fetch_input {
    my ($self) = @_;

    $self->create_analysis;
    if (scalar(@{$self->param('input_files')}) == 0) {
        $self->throw('You did not specify input files for '.$self->analysis->logic_name);
    }
    elsif (scalar(@{$self->param('input_files')}) == 1 and $self->param('options') !~ /-b /) {
        # In samtools merge you can specify a file with a list of files using -b
        $self->input_is_void(1);
    }
    if (defined $self->param('picard_lib')) {
        $self->require_module('Bio::EnsEMBL::Analysis::Runnable::PicardMerge');
        $self->runnable(Bio::EnsEMBL::Analysis::Runnable::PicardMerge->new(
            -program => $self->param('java') || 'java',
            -java_options => $self->param('java_options'),
            -lib => $self->param('picard_lib'),
            -options => $self->param('options'),
            -analysis => $self->analysis,
            -output_file => $self->param('output_file'),
            -input_files => $self->param('input_files'),
            -use_threading => $self->param('use_threading'),
            -samtools => $self->param('samtools') || 'samtools',
            ));
    }
    else {
        $self->require_module('Bio::EnsEMBL::Analysis::Runnable::SamtoolsMerge');
        $self->runnable(Bio::EnsEMBL::Analysis::Runnable::SamtoolsMerge->new(
            -program => $self->param('samtools') || 'samtools',
            -options => $self->param('options'),
            -analysis => $self->analysis,
            -output_file => $self->param('output_file'),
            -input_files => $self->param('input_files'),
            -use_threading => $self->param('use_threading'),
            ));
    }
}

sub run {
    my ($self) = @_;

    foreach my $runnable (@{$self->runnable}) {
        $runnable->run;
        $runnable->check_output_file;
    }
    return 1;
}

sub write_output {
    my ($self) = @_;

    return 1;
}

1;
