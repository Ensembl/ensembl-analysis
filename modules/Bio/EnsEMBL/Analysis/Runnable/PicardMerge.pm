=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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


=head1 DESCRIPTION

Merge BAM files using Picard

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Analysis::Runnable::PicardMerge;

use warnings;
use strict;

use File::Copy;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);
use Bio::EnsEMBL::Utils::Exception qw(throw);

use parent qw(Bio::EnsEMBL::Analysis::Runnable::BaseBamMerge);


=head2 new

 Arg [SAMTOOLS]    : String
 Arg [LIB]         : String
 Arg [JAVA_OPTIONS]: String
 Description: Create a Bio::EnsEMBL::Analysis::Runnable::PicardMerge object to merge BAM files using picard
 Returntype : Bio::EnsEMBL::Analysis::Runnable::PicardMerge
 Exceptions : None

=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    my ($samtools, $lib, $java_options) = rearrange( [ qw( SAMTOOLS LIB JAVA_OPTIONS) ], @args);

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

 Arg [1]    : None
 Description: Merge the BAM files using Picard and rename the index file to be
              useable by samtools
 Returntype : Integer, 1
 Exceptions : Throws if it cannot execute the command
              Throws if it fails to rename the file
              Throws if it fails to update the timestamp on the index file

=cut

sub run {
    my ($self) = @_;

    my $picard_command = $self->{_picard_version} ? 'MergeSamFiles' : '';
    my $cmd = join(' ', $self->program, $self->java_options, '-jar', $self->picard_lib, $picard_command, $self->options, 'OUTPUT='.$self->output_file, 'INPUT='.join(' INPUT=', @{$self->input_files}));
    $cmd .= ' USE_THREADING=true' if ($self->use_threading);
    logger_info($cmd);
    throw('Could not execute picard command '.$cmd) if (system($cmd));
    my $index_file = $self->output_file;
    $index_file =~ s/bam$/bai/;
    # It's needed because Bio::DB::Sam is looking for *.bam.bai and picard create *.bai
    # Return 1 on success
    throw("Failed to move index file $index_file to ".$self->output_file.'.bai') unless (move($index_file, $self->output_file.'.bai'));
    # It's needed because the index can be younger that the bam file but Bio::DB::Sam doesn't like it
    utime(undef, undef, $self->output_file.'.bai') || throw('Failed to update the timestamp for '.$self->output_file.'.bai');

    return 1;
}


=head2 java_options

 Arg [1]    : (optional) String options
 Description: Getter/setter for the java options like Xmx
 Returntype : String
 Exceptions : None

=cut

sub java_options {
    my ($self, $files) = @_;
    if ($files){
        $self->{_java_options} = $files;
    }
    return $self->{_java_options};
}


=head2 picard_lib

 Arg [1]    : (optional) String path
 Description: Getter/setter for the path to the picard library
 Returntype : String
 Exceptions : None

=cut

sub picard_lib {
    my ($self, $files) = @_;
    if ($files){
        $self->{_picard_lib} = $files;
# If the picard lib is picard.jar, we have the "new" version which needs
# a command
        $self->{_picard_version} = $files =~ /picard\.jar$/ ? 1 : 0;
    }
    return $self->{_picard_lib};
}

1;
