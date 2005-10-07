# Ensembl module for Bio::EnsEMBL::Analysis::Runnable::Mavid
#
# Copyright (c) 2005 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Mlagan

=head1 SYNOPSIS

  my $runnable = new Bio::EnsEMBL::Analysis::Runnable::Mlagan
     (-workdir => $workdir,
      -fasta_files => $fasta_files,
      -tree_string => $tree_string,
      -program => "/path/to/program");
  $runnable->run;
  my @output = @{$runnable->output};

=head1 DESCRIPTION

Mavid expects to run the program mavid, a global multiple aligner for large genomic sequences,
using a fasta file and a tree file (Newick format), and eventually a constraints file.
The output (multiple alignment) is parsed and return as a Bio::EnsEMBL::Compara::GenomicAlignBlock object.

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut


package Bio::EnsEMBL::Analysis::Runnable::Mlagan;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception;
use Bio::EnsEMBL::Utils::Argument;
use Bio::EnsEMBL::Compara::DBSQL::ProteinTreeAdaptor;
use Bio::EnsEMBL::Compara::GenomicAlign;
use Bio::EnsEMBL::Compara::GenomicAlignBlock;

use Bio::EnsEMBL::Analysis::Runnable;
our @ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

#$ENV{'LAGAN_DIR'} = "/nfs/acari/abel/ftp/lagan12_i386";
#my $BIN_DIR = "/nfs/acari/abel/ftp/lagan12_i386";
$ENV{'LAGAN_DIR'} = "/nfs/acari/abel/ftp/lagan12_turing";
my $BIN_DIR = "/nfs/acari/abel/ftp/lagan12_turing";

#if (! -e "/proc/version") {
#  $ENV{'LAGAN_DIR'} = "/nfs/acari/abel/ftp/lagan12_alpha";
#  $BIN_DIR = "/nfs/acari/abel/ftp/lagan12_alpha";
#}

print "BIN_DIR $BIN_DIR\n";
print $ENV{'MACHTYPE'},"\n";

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
  my ($workdir, $fasta_files, $tree_string, $parameters) = rearrange(['WORKDIR', 'FASTA_FILES', 'TREE_STRING','CONSTRAINTS_FILE'], @args);

  unless (defined $self->program) {
    $self->program("$BIN_DIR/mlagan");
  }
  chdir $self->workdir;
  $self->fasta_files($fasta_files) if (defined $fasta_files);
  if (defined $tree_string) {
    $self->tree_string($tree_string);
  }

  return $self;
}

sub fasta_files {
  my $self = shift;
  $self->{'_fasta_files'} = shift if(@_);
  return $self->{'_fasta_files'};
}

sub tree_string {
  my $self = shift;
  $self->{'_tree_string'} = shift if(@_);
  return $self->{'_tree_string'};
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

  $self->run_mlagan;

  $self->parse_results;
#  unlink(<$dir/*>);
#  rmdir($dir);

  return 1;
}

sub run_mlagan {
  my $self = shift;

  chdir $self->workdir;

  throw($self->program . " is not executable Mlagan::run_analysis ")
    unless ($self->program && -x $self->program);

  my $command = $self->program;
  foreach my $fasta_file (@{$self->fasta_files}) {
    $command .= " $fasta_file";
  }
  if ($self->tree_string) {
    $command .= " -tree \"" . $self->tree_string . "\"";
  }
  if ($self->parameters) {
    $command .= " " . $self->parameters;
  }
  $command .= " -out mlagan.mfa";
  print "Running mlagan " . $command . "\n";
  unless (system($command) == 0) {
    throw("mlagan execution failed\n");
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

  my $alignment_file = $self->workdir . "/mlagan.mfa";
  my $gab = new Bio::EnsEMBL::Compara::GenomicAlignBlock;
#  $gab->method_link_species_set($self->{'method_link_species_set'});
#  $gab->genomic_align_array([$genomic_align1, $genomic_align2]);
#  $gab->score($fp->score);
#  $gab->perc_id($fp->percent_id);
#  $gab->length($fp->alignment_length);


  open F, $alignment_file || throw("Could not open $alignment_file");
  my $id;
  my $seq = "";
  while (<F>) {
    next if (/^\s*$/);
    chomp;
    if (/^>(\S+)\s*.*$/) {
      my $new_id = $1;
      $new_id =~ s/_aligned//;
      if (defined $id && defined $seq) {
        my $ga = new Bio::EnsEMBL::Compara::GenomicAlign;
#        $ga->method_link_species_set($self->{'method_link_species_set'});
        $ga->dnafrag_id($id);
#        $ga->dnafrag_start($qyChunk->seq_start + $fp->start -1);
#        $ga->dnafrag_end($qyChunk->seq_start + $fp->end -1);
#        $ga->dnafrag_strand($fp->strand);
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
  #        $ga->method_link_species_set($self->{'method_link_species_set'});
  $ga->dnafrag_id($id);
  #        $ga->dnafrag_start($qyChunk->seq_start + $fp->start -1);
  #        $ga->dnafrag_end($qyChunk->seq_start + $fp->end -1);
  #        $ga->dnafrag_strand($fp->strand);
  $ga->aligned_sequence($seq);
  $gab->add_GenomicAlign($ga);
  
  $self->output([$gab]);
}


1;
