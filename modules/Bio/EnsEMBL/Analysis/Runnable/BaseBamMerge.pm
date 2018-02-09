=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Runnable::BaseBamMerge -

=head1 SYNOPSIS

  Do NOT instantiate this class directly: must be instantiated
  from a subclass (see ExonerateTranscript, for instance).

=head1 DESCRIPTION

This is an abstract superclass to merge BAM files using either
samtools or picard

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Analysis::Runnable::BaseBamMerge;

use warnings;
use strict;

use Bio::EnsEMBL::Analysis::Runnable::Samtools;
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use parent qw(Bio::EnsEMBL::Analysis::Runnable);


=head2 new

 Arg [OUTPUT_FILE]  : String
 Arg [INPUT_FILES]  : Arrayref String
 Arg [USE_THREADING]: Integer
 Arg [VERBOSE]      : Integer
 Description        : Creates an object to merge BAM files
 Returntype         : Bio::EnsEMBL::Analysis::Runnable::BaseBamMerge
 Exceptions         : Throws if INPUT_FILES is not defined or does not exists
                      Throws if OUTPUT_FILE is not defined

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($output_file, $input_files, $use_threading, $verbose) = rearrange( [ qw( OUTPUT_FILE INPUT_FILES USE_THREADING VERBOSE) ], @args);

  logger_verbosity($verbose) if $verbose;
  $self->use_threading($use_threading);

  if (defined($input_files)) {
    my $missing = 0;
    if(ref($input_files) ne "ARRAY"){
      throw("You must supply an array reference with -input_files");
    }
    else {
        foreach my $file (@$input_files) {
            if (! -e $file) {
                warning($file.' is missing');
                $missing++;
            }
        }
    }
    throw('At least one of the input file does not exist!') if ($missing);
    $self->input_files($input_files);
  }
  else {
      throw('You need to give a array of input files');
  }
  if (defined $output_file) {
    warning('The given output file '.$output_file.' exists') if ( -e $output_file);
    $self->output_file($output_file);
  }
  else {
      throw('You need to specify an output file');
  }

  return $self;
}



############################################################
#
# Analysis methods
#
############################################################

=head2 run

 Arg [1]    : None
 Description: This method should be implemented by the inheriting module
 Returntype : None
 Exceptions : Throws if not implemented

=cut

sub run {
  my ($self) = @_;

  throw('Need to be implemented by the inheriting object');
}


=head2 check_output_file

 Arg [1]    : None
 Description: Check if the output file has been successfully created using
              flagstat from Bio::EnsEMBL::Analysis::Runnable::Samtools
              Store the values in $self->output
 Returntype : None
 Exceptions : None

=cut

sub check_output_file {
    my ($self) = @_;

    $self->output($self->samtools->flagstat($self->output_file, 1));
}


=head2 input_files

 Arg [1]    : Arrayref String files
 Description: Getter/Setter for the input files
 Returntype : Arrayref
 Exceptions : None

=cut

sub input_files {
    my ($self, $files) = @_;
    if ($files){
        $self->{_input_files} = $files;
    }
    return $self->{_input_files};
}


=head2 output_file

 Arg [1]    : String file
 Description: Getter/Setter for the output files
 Returntype : String
 Exceptions : None

=cut

sub output_file {
    my ($self, $file) = @_;

    if ($file) {
        $self->{_output_file} = $file;
    }
    return $self->{_output_file};
}


=head2 use_threading

 Arg [1]    : Boolean, 0 or 1
 Example    : $sef->use_threading(1);
 Description: Getter/Setter, set this to 1 if you want picard to use threads, is ~20% faster but use ~20% more memory
 Returntype : Boolean
 Exceptions : None

=cut

sub use_threading {
    my ($self, $file) = @_;

    if ($file) {
        $self->{_use_threading} = $file;
    }
    return $self->{_use_threading};
}


=head2 samtools

 Arg [1]    : String, '/usr/bin/samtools'
 Example    : $self->samtools('/usr/bin/samtools');
 Description: Creates a object to run samtools command
 Returntype : Bio::EnsEMBL::Analysis::Runnable::Samtools
 Exceptions : None

=cut

sub samtools {
    my ($self, $samtools, $verbose) = @_;

    if ($samtools) {
        $self->{_samtools} = Bio::EnsEMBL::Analysis::Runnable::Samtools->new(
                              -program => $samtools,
                              -verbose => $verbose,
                              -use_threading => $self->use_threading
                              );
    }
    return $self->{_samtools};
}

1;
