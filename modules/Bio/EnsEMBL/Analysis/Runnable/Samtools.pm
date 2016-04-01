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

Bio::EnsEMBL::Analysis::Runnable::Samtools -

=head1 SYNOPSIS

  Do NOT instantiate this class directly: must be instantiated
  from a subclass (see ExonerateTranscript, for instance).

=head1 DESCRIPTION


=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Analysis::Runnable::Samtools;

use warnings;
use strict;

use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_verbosity logger_info);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );


sub new {
  my ($class,@args) = @_;
  my $self = bless {},$class;
  my ($program, $use_threading, $verbose) = rearrange( [ qw( PROGRAM USE_THREADING VERBOSE) ], @args);


  $program = 'samtools' unless ($program);
  throw($program.'samtools cannot be executed ') unless (-x $program);
  $self->program($program);
  $self->use_threading($use_threading);
  logger_verbosity($verbose) if ($verbose);

  return $self;
}



############################################################
#
# Analysis methods
#
############################################################

=head2 run_command

Usage   :   $obj->run_command($command, $options, $output_file, $input_file)
Function:   BE VERY CAREFUL: this run any samtools command, if you specify an output and input files
            it will run the command like this: samtools <command> <options> <outfile> <infile>
            This works for most of the command but read the samtools man page first
=cut

sub run_command {
  my ($self, $command, $options, $output_file, $input_file) = @_;

  my $cmd = join(' ', $self->program, $command, $options);
  $cmd .= ' -@ '.$self->use_threading if ($self->use_threading);
  if ($command eq 'view') {
      $cmd .= join(' ', '-o', $output_file, $input_file);
  }
  else {
      $cmd .= join(' ', $output_file, $input_file);
  }
  throw('Command '.$cmd.' failed') if (system($cmd));
}


sub merge {
  my ($self, $options, $output_file, $input_files) = @_;

  my $cmd = join(' ', $self->program, 'merge', $options);
  $cmd .= ' -@ '.$self->use_threading if ($self->use_threading);
  if (ref($input_files) ne 'ARRAY') {
       throw('input file does not exist: '.$input_files) unless (-e $input_files);
       $cmd .= join(' ', '-b', $input_files, $output_file);
  }
  else {
      $cmd .= join(' ', $output_file, @{$self->input_files});
  }
  logger_info($cmd);
  throw($output_file.' merge failed') if (system($cmd));
  $self->index($output_file);
}

sub index {
  my ($self, $file) = @_;

  my $cmd = join(' ', $self->program, 'index', $file);
  logger_info($cmd);
  throw($file.' indexing failed') if (system($cmd));
}


sub flagstat {
  my ($self, $file, $stat) = @_;

  my $cmd = join(' ', $self->program, 'flagstat', $file);
  my @output;
  open(CMD, $cmd.' 2>&1 | ') || throw("Could not open command $cmd");
  while(<CMD>) {
      throw($file.' is truncated') if (/truncated/);
      if ($stat and /^\s*(\d+).*in total/) {
          $output[0] = $1;
      }
      elsif ($stat and /^\s*(\d+).*mapped\s+\(/) {
          $output[1] = $1;
      }
      elsif ($stat and /\s*(\d+).*properly paired \(/) {
          $output[2] = $1;
      }
      else {
          push(@output, $_);
      }
  }
  close(CMD) || throw("Could not close $cmd");
  return \@output;
}


sub program {
    my ($self, $file) = @_;

    if ($file) {
        $self->{_program} = $file;
    }
    return $self->{_program};
}


sub use_threading {
    my ($self, $file) = @_;

    if ($file) {
        $self->{_use_threading} = $file;
    }
    return $self->{_use_threading};
}


1;
