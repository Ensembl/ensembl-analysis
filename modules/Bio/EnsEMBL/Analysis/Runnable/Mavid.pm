=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Runnable::Mavid - 

=head1 SYNOPSIS

  my $runnable = new Bio::EnsEMBL::Analysis::Runnable::Mavid
     (-input_dir => $input_dir,
      -output_dir => $output_dir,
      -genome_names => ["homo_sapiens","mus_musculus","rattus_norvegicus"],
      -program => "/path/to/program");
  $runnable->run;
  my @output = @{$runnable->output};

=head1 DESCRIPTION

Mavid expects to run the program mavid, a global multiple aligner for large genomic sequences,
using a fasta file and a tree file (Newick format), and eventually a constraints file.
The output (multiple alignment) is parsed and return as a Bio::EnsEMBL::Compara::GenomicAlignBlock object.

=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::Runnable::Mavid;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception;
use Bio::EnsEMBL::Utils::Argument;
use Bio::EnsEMBL::Compara::Graph::NewickParser;
use Bio::EnsEMBL::Compara::GenomicAlign;
use Bio::EnsEMBL::Compara::GenomicAlignBlock;

use Bio::EnsEMBL::Analysis::Runnable;
our @ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

my $BIN_DIR = "/nfs/acari/abel/ftp/mavid_mercator/mavid0.9new/";

=head2 new

  Arg [1]   : -workdir => "/path/to/working/directory"
  Arg [2]   : -fasta_file => "/path/to/fasta/file"
  Arg [3]   : -tree_file => "/path/to/tree/file" (optional)
  Arg [4]   : -constraints_file => "/path/to/constraints/file" (optional)

  Function  : contruct a new Bio::EnsEMBL::Analysis::Runnable::Mavid
  runnable
  Returntype: Bio::EnsEMBL::Analysis::Runnable::Mavid
  Exceptions: none
  Example   :

=cut


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($workdir, $fasta_file, $tree_file, $constraints_file) = rearrange(['WORKDIR', 'FASTA_FILE', 'TREE_FILE','CONSTRAINTS_FILE'], @args);

  unless (defined $self->program) {
    $self->program("$BIN_DIR/mavid");
  }
  chdir $self->workdir;
  $self->fasta_file($fasta_file) if (defined $fasta_file);
  if (defined $tree_file) {
    $self->tree_file($tree_file);
  } else {
    $self->random_tree(1);
    $self->improve_tree_with_phyml(1);
  }
  $self->constraints_file($constraints_file) if (defined $constraints_file);

  return $self;
}

sub fasta_file {
  my $self = shift;
  $self->{'_fasta_file'} = shift if(@_);
  return $self->{'_fasta_file'};
}

sub tree_file {
  my $self = shift;
  $self->{'_tree_file'} = shift if(@_);
  return $self->{'_tree_file'};
}

sub constraints_file {
  my $self = shift;
  $self->{'_constraints_file'} = shift if(@_);
  return $self->{'_constraints_file'};
}

sub run_number {
  my $self = shift;
  $self->{'_run_number'} = shift if(@_);
  return $self->{'_run_number'};
}

sub random_tree {
  my $self = shift;
  $self->{'_random_tree'} = shift if(@_);
  return $self->{'_random_tree'};
}

sub improve_tree_with_phyml {
  my $self = shift;
  $self->{'_improve_tree_with_phyml'} = shift if(@_);
  return $self->{'_improve_tree_with_phyml'};
}

=head2 run_analysis

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Mavid
  Arg [2]   : string, program name
  Function  : create and open a commandline for the program trf
  Returntype: none
  Exceptions: throws if the program in not executable or if the results
  file doesnt exist
  Example   : 

=cut

sub run_analysis {
  my ($self, $program) = @_;
  $self->run_number(1);

  $self->run_mavid($self->run_number);

  while ($self->improve_tree_with_phyml) {
    $self->run_phyml($self->run_number);
    last unless ($self->improve_tree_with_phyml);
    $self->run_mavid($self->run_number);
  }
  $self->parse_results($self->run_number - 1);
  for (my $i = 1; $i <= $self->run_number; $i++) {
    my $dir = $self->workdir . "/run_number_" . $i;
    unlink(<$dir/*>);
    rmdir($dir);
  }
  return 1;
}

=head2 parse_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::Mavid
  Arg [2]   : string, results filename
  Function  : parse the specifed file and produce RepeatFeatures
  Returntype: nine
  Exceptions: throws if fails to open or close the results file
  Example   : 

=cut


sub parse_results{
  my ($self, $run_number) = @_;

  $run_number = 1 unless ($run_number);

  my $alignment_file = $self->workdir . "/run_number_" . $run_number . "/mavid.mfa";
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

sub run_mavid {
  my ($self, $run_number) = @_;

  my $run_dir = $self->workdir . "/run_number_$run_number";
  if (! -e $run_dir) {
    mkdir($run_dir, 0777);
  }
  chdir $run_dir;

  throw($self->program." is not executable Mavid::run_analysis ") 
    unless($self->program && -x $self->program);

  if ($self->random_tree) {
    my $command = "$BIN_DIR/randtree " . $self->fasta_file . " > " . $run_dir . "/treefile";
    unless (system($command) == 0) {
      throw("Couldn't generate random tree, $!\n");
    }
    $self->tree_file($run_dir . "/treefile");
  }

  my $command = $self->program . " ";
  if ($self->constraints_file) {
    $command .= "-c " . $self->contraints_file . " ";
  }
  $command .= $self->tree_file . " " . $self->fasta_file;
  print "Running mavid ".$command."\n";
  unless (system($command) == 0) {
    throw("mavid execution failed\n");
  }
  $self->random_tree(0);
  $self->run_number($self->run_number + 1);
  chdir $self->workdir;
}

sub run_phyml {
  my ($self, $run_number) = @_;

  my $run_dir_minus_1 = $self->workdir . "/run_number_". ($run_number - 1);
  my $run_dir_minus_2 = $self->workdir . "/run_number_". ($run_number - 2);
  my $program = "phyml";
  my $command = "$program ". $run_dir_minus_1 . "/mavid.phy 0 i 3 0 JC69 e e 4 e BIONJ y y";
  print "Running phyml: $command\n";
  unless (system($command) == 0) {
    throw("phyml execution failed\n");
  }
  $self->tree_file($self->binary_tree_file($run_dir_minus_1 . "/mavid.phy_phyml_tree.txt"));
  if ($run_number >= 3) {
    unless ($self->better_likelihood($run_dir_minus_2 . "/mavid.phy_phyml_lk.txt",$run_dir_minus_1 . "/mavid.phy_phyml_lk.txt")) {
      $self->improve_tree_with_phyml(0);
    }
  }
}

sub binary_tree_file {
  my $self = shift;
  my $tree_file = shift;
  open F, $tree_file || throw("Can not open $tree_file");
  my $newick = "";
  while (<F>) {
    chomp;
    if (/^\s*(.*)\s*$/) {
      $newick .= $1;
    }
  }
  close F;
  my $tree = Bio::EnsEMBL::Compara::Graph::NewickParser::parse_newick_into_tree($newick);
  $self->binary_tree($tree);
  my $binary_tree_file = $tree_file . ".bin";
  open F, ">$binary_tree_file";
  print F $tree->newick_format('simple'),"\n";
  close F;
  $tree->release;
  return $binary_tree_file;
}

sub binary_tree {
  my $self = shift;
  my $tree = shift;

  if ($tree->get_child_count > 2) {
    my $children = $tree->children;
    my %distance_between_two_children;
    while (my $c1 = shift @{$children}) {
      foreach my $c2 (@{$children}) {
        my $ancestor = $c1->find_first_shared_ancestor($c2);
        my $dist = $c1->distance_to_ancestor($ancestor) + $c2->distance_to_ancestor($ancestor);
        $distance_between_two_children{$dist} = [$c1, $c2];
      }
    }
    my ($two_children) = map($distance_between_two_children{$_}, sort {$a <=> $b} keys %distance_between_two_children);
    my $new_node = new Bio::EnsEMBL::Compara::NestedSet;
    $new_node->distance_to_parent(0);
    $tree->add_child($new_node);
    map($new_node->add_child($_), @{$two_children});
    $self->binary_tree($tree);
  } else {
    foreach my $c (@{$tree->children}) {
      $self->binary_tree($c);
    }
  }
}

sub better_likelihood {
  my $self = shift;
  my $lk_file1 = shift;
  my $lk_file2 = shift;
  my ($lk1, $lk2);
  if (! -e $lk_file1) {
    throw("$lk_file1 does not exist");
  }
  open F, $lk_file1 || throw("Can not open $lk_file1");
  while (<F>) {
    if (/^([+-]?\d+\.\d+)$/) {
      $lk1 = $1;
      last;
    }
  }
  close F;
  if (! -e $lk_file2) {
    throw("$lk_file2 does not exist");
  }
  open F, $lk_file2 || throw("Can not open $lk_file2");
  while (<F>) {
    if (/^([+-]?\d+\.\d+)$/) {
      $lk2 = $1;
      last;
    }
  }
  close F;
  if ($lk2 > $lk1) {
    return 1;
  } else {
    return 0;
  }
}

1;
