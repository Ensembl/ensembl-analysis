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

Bio::EnsEMBL::Analysis::Runnable::PicardMerge -

=head1 SYNOPSIS

  Do NOT instantiate this class directly: must be instantiated
  from a subclass (see ExonerateTranscript, for instance).

=head1 DESCRIPTION

This is an abstract superclass to handle the common functionality for
Exonerate runnables: namely it provides
- a consistent external interface to drive exonerate regardless of
  what features youre finally producing, and
- a common process to stop people duplicating function (eg how to
  arrange command-line arguments, what ryo-string to use etc).

It does NOT provide the parser to convert the exonerate output
into Transcripts or AffyFeatures etc. That is the job of the
subclasses, which MUST implement the parse_results method.

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Analysis::Runnable::PicardMerge;

use warnings;
use strict;
use vars qw(@ISA);

use File::Copy;
use Bio::EnsEMBL::Analysis::Runnable::BaseBamMerge;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_verbosity logger_info);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::BaseBamMerge);


sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    my ($samtools, $lib, $java_options) = rearrange( [ qw( SAMTOOLS LIB JAVA_OPTIONS) ], @args);

    logger_verbosity($verbose) if $verbose;

    $self->samtools($samtools);
    $self->picard_lib($lib);
    $self->java_options($java_options);

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

    my $cmd = join(' ', $self->program, $self->java_options, '-jar', $self->picard_lib, $self->options, 'OUTPUT='.$self->output_file, join(' INPUT=', @{$self->input_files}));
    $cmd .= ' USE_THREADING=true' if ($self->use_threading);
    throw('Could not execute picard command '.$cmd) unless (system($cmd));
    my $index_file = $self->output_file;
    $index_file =~ s/bam$/bai/;
    # It's needed because Bio::DB::Sam is looking for *.bam.bai and picard create *.bai
    # Return 1 on success
    throw("Failed to move index file $index_file to $outfile.bai") unless (move($index_file, $outfile.'.bai'));
    # It's needed because the index can be younger that the bam file but Bio::DB::Sam doesn't like it
    # Return 0 on success
    throw("Failed to update the timestamp for $outfile.bai") if (system("touch $outfile.bai"));
    $self->check_output_file;

    return 1;
}

sub java_options {
    my ($self, $files) = @_;
    if ($files){
        $self->{_java_options} = $files;
    }
    return $self->{_java_options};
}


sub picard_lib {
    my ($self, $files) = @_;
    if ($files){
        $self->{_picard_lib} = $files;
    }
    return $self->{_picard_lib};
}


1;
