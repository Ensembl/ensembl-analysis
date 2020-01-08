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
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(locate_executable execute_with_wait);


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

  my $cmd = $self->make_commandline($command, $options, $output_file, $input_file);
  execute_with_wait($cmd, 'Command '.$cmd.' failed');
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

  my $cmd = $self->_base_command('merge', $options);
  if (ref($input_files) ne 'ARRAY') {
       throw('input file does not exist: '.$input_files) unless (-e $input_files);
       $cmd .= join(' ', '-b', $input_files, $output_file);
  }
  else {
      $cmd .= join(' ', $output_file, @{$input_files});
  }
  logger_info($cmd);
  execute_with_wait($cmd, $output_file.' merge failed');
  $self->index($output_file);
}

=head2 _base_command

 Arg [1]    : String $command
 Arg [2]    : String (optional) $options
 Description: Return the base of a samtools command. It adds the number of threads
              if you use use_threading.
 Returntype : String command
 Exceptions : None

=cut

sub _base_command {
  my ($self, $command, @options) = @_;

  my $cmd = $self->program.' '.$command;
  $cmd .= ' '.join(' ', @options) if (defined $options[0]);
  if ($self->version) {
    $cmd .= ' -@ '.$self->use_threading if ($self->use_threading);
  }
  return $cmd.' ';
}

=head2 sort

 Arg [1]    : String output_file
 Arg [2]    : String input_file
 Arg [3]    : String (optional) temp_prefix
 Arg [4]    : String (optional) options
 Description: Sort a BAM file, it should be smart enough to do it the right way
 Returntype : None
 Exceptions : Throws if it cannot execute the command

=cut

sub sort {
  my ($self, $output_file, $input_file, $temp_prefix, $options) = @_;

  my $cmd = $self->_base_command('sort', $options);
  if($self->version) {
    $temp_prefix = $output_file.$$ unless ($temp_prefix);
    $cmd .= ' -T '.$temp_prefix.' -o '.$output_file.'.bam '.$input_file;
  }
  else {
    $cmd .= $input_file.' '.$output_file;
  }
  execute_with_wait($cmd, 'Sort command: '.$cmd.' failed');
}

=head2 index

 Arg [1]    : String filename
 Arg [2]    : String (optional) filename
 Example    : $samtools->index($filename);
 Description: Index a BAM file
 Returntype : None
 Exceptions : Throws if the command failed

=cut

sub index {
  my ($self, $file, $options) = @_;

  $options = '' unless (defined $options);
  my $cmd = join(' ', $self->_base_command('index', $options), $file);
  logger_info($cmd);
  execute_with_wait($cmd, $file.' indexing failed');
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

  my $cmd = join(' ', $self->_base_command( 'flagstat'), $file);
  my @output;
  open(CMD, $cmd.' 2>&1 | ') || throw("Could not open command $cmd");
  while(<CMD>) {
      throw($file.' is truncated, something went wrong: '.$_)
        if (/truncated/ or /EOF marker is absent/ or /invalid BAM binary header/);
      if ($stat and /^\s*(\d+).*in total/) {
          $output[0] = $1;
      }
      elsif ($stat and /^\s*\d+.*mapped\s+\(\s*([0-9\.]+)\s*%/) {
          $output[1] = $1;
      }
      elsif ($stat and /\s*\d+.*properly paired \(\s*([0-9\.]+)\s*%/) {
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

 Arg [1]    : String (optional) $program
 Description: Getter/setter for the binary. It tries to locate properly the executable
 Returntype : String
 Exceptions : Throw if it can't find the executable Arg[1]
              Throw if it has problem guessing the version

=cut

sub program {
    my ($self, $file) = @_;

    if ($file) {
        $self->{_program} = locate_executable($file);
        $self->{_version} = 1;
        open(FH, $self->{_program}.' 2>&1 |') || throw('Could not run samtools');
        while (<FH>) {
          if (/Version:\s+(\d)/) {
            $self->{_version} = $1;
            last;
          }
        }
        close(FH); # samtools exits with 1 so we should never check the result here, this is just for checking which version so it's OK
        warning('It seems you are on a 0.x.x version, we recommend you to update samtools') unless $self->version;
    }
    return $self->{_program};
}

=head2 version

 Arg [1]    : None
 Description: Getter for the version of samtools to know how to construct commands
 Returntype : String version
 Exceptions : None

=cut

sub version {
  my ($self) = @_;

  return $self->{_version};
}

=head2 use_threading

 Arg [1]    : Integer (optional)
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


=head2 make_commandline

 Arg [1]    : String, a samtools command like view, index,...
 Arg [2]    : String, samtools options
 Arg [3]    : String, output file
 Arg [4]    : String, input file
 Example    : $samtools->make_commandline('view', '-b', $output_file, $input_file);
 Description: Prepare a samtools command. If use_threading is defined it will use
              the value in use_threading for the number of threads
 Returntype : String $command
 Exceptions : None

=cut

sub make_commandline {
  my ($self, $command, $options, $output_file, $input_file) = @_;

  my $cmd = $self->_base_command($command, $options);
  if ($command eq 'view' or $self->version) {
      $cmd .= join(' ', '-o', $output_file, $input_file);
  }
  else {
      $cmd .= join(' ', $output_file, $input_file);
  }
  return $cmd;
}


=head2 index_genome

 Arg [1]    : String, path to a fasta file containting the genome
 Description: Index the genome using faidx
 Returntype : None
 Exceptions : Throws if Arg[1] does not exist
              Throws if command fails

=cut

sub index_genome {
  my ($self, $genome) = @_;

  throw("Could not find genome: $genome") unless (-e $genome);
  execute_with_wait($self->program." faidx $genome");
}


=head2 reheader

 Arg [1]    : String, path to header file
 Arg [2]    : String, source BAM file
 Arg [3]    : String, path to target BAM file
 Description: Use the reheader command from samtools to change the
              header of a BAM file
 Returntype : None
 Exceptions : Throws if Arg[1] does not exist
              Throws if Arg[2] does not exist

=cut

sub reheader {
  my ($self, $header, $source_bam, $target_bam) = @_;
  throw("Could not find $source_bam") unless (-e $source_bam);
  throw("Could not find $header") unless (-e $header);
  warning("Will overwrite $target_bam") if (-e $target_bam);
  execute_with_wait($self->program." reheader $header $source_bam > $target_bam");
}

1;
