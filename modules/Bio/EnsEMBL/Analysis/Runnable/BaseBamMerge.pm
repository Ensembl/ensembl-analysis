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

Bio::EnsEMBL::Analysis::Runnable::SamtoolsMerge -

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

package Bio::EnsEMBL::Analysis::Runnable::SamtoolsMerge;

use warnings;
use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Analysis::Runnable::Samtools;
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_verbosity logger_info);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);


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
        foreach my $file (@{$self->INPUT_FILES}) {
            if (! -e $file) {
                warning($file.' is missing');
                $missing++;
            }
        }
    }
    throw('At least one of the input file does not exist!') if ($missing);
    $self->input_files($input_files);
  }
  elsif ($self->options =~ /-b\s+(\S+)/){
      throw('Could not access '.$1) unless (-e $1);
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

Usage   :   $obj->run($workdir, $args)
Function:   Runs exonerate script and puts the results into the file $self->results
            It calls $self->parse_results, and results are stored in $self->output
=cut

sub run {
  my ($self) = @_;

  throw('Need to be implemented by the inheriting object');
}

=head2 parse_results

  Arg [1]   :
  Arg [2]   :
  Function  :
  Returntype:
  Example   :

=cut

sub check_output_file {
    my ($self) = @_;

    $self->output($self->samtools->flagstat($self->output_file, 1));
}



sub input_files {
    my ($self, $files) = @_;
    if ($files){
        $self->{_input_files} = $files;
    }
    return $self->{_input_files};
}


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
        $self->{_samtools} = Bio::EnsEMBL::Analysis::Runnable::Samtools->new(-program => $samtools, -verbose => $verbose, -use_threading => $self->use_threading);
    }
    return $self->{_samtools};
}


1;
