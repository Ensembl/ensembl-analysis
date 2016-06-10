=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016] EMBL-European Bioinformatics Institute
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


=head1 DESCRIPTION

Run samtools in a Perl OO manner

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


=head2 new

 Arg [PROGRAM]      : String
 Arg [USE_THREADING]: Integer
 Arg [VERBOSE]      : Integer
 Description: Creates a Bio::EnsEMBL::Analysis::Runnable::Samtools object to run samtools command
 Returntype : Bio::EnsEMBL::Analysis::Runnable::Samtools
 Exceptions : None

=cut

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

 Arg [1]    : String a samtools command like view, index,...
 Arg [2]    : String samtools options
 Arg [3]    : String output file
 Arg [4]    : String input file
 Example    : $samtools->run_command('view', '-b', $output_file, $input_file);
 Description: Run a samtools command. If use_threading is defined it will ask
              for the number of threads in use_threading
 Returntype : None
 Exceptions : Throws if samtools failed

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


=head2 merge

 Arg [1]    : String options
 Arg [2]    : String output file
 Arg [3]    : Arrayref input files or the filename is -b is used
 Example    : $samtools->merge(undef, $outfile, \@infiles;
 Description: Merge the BAM files and index the output file
 Returntype : None
 Exceptions : Throws if you don't give an array unless you use -b
              and an existing file
              Throws if the merge fails

=cut

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


=head2 index

 Arg [1]    : String filename
 Example    : $samtools->index($filename);
 Description: Index a BAM file
 Returntype : None
 Exceptions : Throws if the command failed

=cut

sub index {
  my ($self, $file) = @_;

  my $cmd = join(' ', $self->program, 'index', $file);
  logger_info($cmd);
  throw($file.' indexing failed') if (system($cmd));
}


=head2 flagstat

 Arg [1]    : String filename
 Arg [2]    : Boolean mapping_precentage
 Example    : $samtools->flagstat($file, 1);
 Description: Run flagstat. If you set Arg[2] to 1, it will only return the
              total number of reads, the percentage of mapped reads and the
              percentage of properly paired reads. Otherwise it will return
              the full ouptput from flagstat.
 Returntype : Arrayref
 Exceptions : Throws if the BAM file is truncated
              Throws if the command failed

=cut

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


=head2 program

 Arg [1]    : (optional) String
 Description: Getter/setter for the binary
 Returntype : String
 Exceptions : None

=cut

sub program {
    my ($self, $file) = @_;

    if ($file) {
        $self->{_program} = $file;
    }
    return $self->{_program};
}


=head2 use_threading

 Arg [1]    : (optional) Integer
 Description: Getter/setter for the number of threads to use
 Returntype : Integer
 Exceptions : None

=cut

sub use_threading {
    my ($self, $file) = @_;

    if ($file) {
        $self->{_use_threading} = $file;
    }
    return $self->{_use_threading};
}


1;
