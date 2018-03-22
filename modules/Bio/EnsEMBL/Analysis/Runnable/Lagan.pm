=head1 LICENSE

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Runnable::Lagan - 

=head1 SYNOPSIS

  my $runnable = new Bio::EnsEMBL::Analysis::Runnable::Lagan
     (-workdir => $workdir,
      -fasta_files => $fasta_files,
      -program => "/path/to/program");
  $runnable->run;
  my @output = @{$runnable->output};

=head1 DESCRIPTION

Lagan expects to run the program lagan.pl, a global pairwise aligner,
using a couple of fasta files. The output (multiple alignment) is parsed
and returned as a Bio::EnsEMBL::Compara::GenomicAlignBlock object.
=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::Runnable::Lagan;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception;
use Bio::EnsEMBL::Utils::Argument;
use Bio::EnsEMBL::Compara::GenomicAlign;
use Bio::EnsEMBL::Compara::GenomicAlignBlock;

use Bio::EnsEMBL::Analysis::Runnable;
our @ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

$ENV{'LAGAN_DIR'} = "/usr/local/ensembl/lagan-1.21/";
my $BIN_DIR = "/usr/local/ensembl/lagan-1.21/";

=head2 new

  Arg [1]   : -workdir => "/path/to/working/directory"
  Arg [2]   : -fasta_files => "/path/to/fasta/file"
  Arg [3]   : -tree_string => "/path/to/tree/file" (optional)
  Arg [4]   : -parameters => "parameter" (optional)

  Function  : contruct a new Bio::EnsEMBL::Analysis::Runnable::Mlagan
  runnable
  Returntype: Bio::EnsEMBL::Analysis::Runnable::Mlagan
  Exceptions: none
  Example   :

=cut


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($workdir, $fasta_files, $tree_string, $parameters) = rearrange(['WORKDIR', 'FASTA_FILES', 'CONSTRAINTS_FILE'], @args);

  unless (defined $self->program) {
    $self->program("$BIN_DIR/lagan.pl");
  }
  chdir $self->workdir;
  $self->fasta_files($fasta_files) if (defined $fasta_files);

  return $self;
}

sub fasta_files {
  my $self = shift;
  $self->{'_fasta_files'} = shift if(@_);
  foreach my $this_file (@{$self->{'_fasta_files'}}) {
#     $self->files_to_delete($this_file);
  }
  return $self->{'_fasta_files'};
}

sub parameters {
  my $self = shift;
  $self->{'_parameters'} = shift if(@_);
  return $self->{'_parameters'};
}

=head2 run_analysis

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Mlagan
  Arg [2]   : string, program name
  Function  : create and open a commandline for the program trf
  Returntype: none
  Exceptions: throws if the program in not executable or if the results
  file doesnt exist
  Example   : 

=cut

sub run_analysis {
  my ($self, $program) = @_;

  $self->run_lagan;

  $self->parse_results;

  return 1;
}

sub run_lagan {
  my $self = shift;

  chdir $self->workdir;

  throw($self->program . " is not executable Lagan::run_analysis ")
    unless ($self->program && -x $self->program);

  my $command = $self->program;
  foreach my $fasta_file (@{$self->fasta_files}) {
    $command .= " $fasta_file";
  }
  if ($self->parameters) {
    $command .= " " . $self->parameters;
  }
  $command .= " -mfa -out lagan.mfa 1>/dev/null 2>&1";
  print "Running mlagan " . $command . "\n";
  unless (system($command) == 0) {
    throw("lagan execution failed\n");
  }
}

=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Mlagan
  Function  : parse the specifed file and produce RepeatFeatures
  Returntype: nine
  Exceptions: throws if fails to open or close the results file
  Example   : 

=cut


sub parse_results{
  my ($self, $run_number) = @_;

  my $alignment_file = $self->workdir . "/lagan.mfa";
  my $gab = new Bio::EnsEMBL::Compara::GenomicAlignBlock;

  open F, $alignment_file || throw("Could not open $alignment_file");
  my $id;
  my $seq = "";
  while (<F>) {
    next if (/^\s*$/);
    chomp;
    ## The "DnaFrag" and "."
    ## are used to avoid having two sequence names where one of them is a
    ## substring of the other (this confuses MLAGAN) and maybe Lagan as well... :-)
    if (/^>DnaFrag(\d+)\./) {
      my $new_id = $1;
      if (defined $id && defined $seq) {
        my $ga = new Bio::EnsEMBL::Compara::GenomicAlign;
        $ga->dnafrag_id($id);
        $ga->aligned_sequence($seq);
        $gab->add_GenomicAlign($ga);
        $id = undef;
        $seq = "";
      }
      $id = $new_id;
    } else {
      $seq .= $_;
    }
  }
  close F;
  my $ga = new Bio::EnsEMBL::Compara::GenomicAlign;
  $ga->dnafrag_id($id);
  $ga->aligned_sequence($seq);
  $gab->add_GenomicAlign($ga);
  
  $self->output([$gab]);
}


1;
