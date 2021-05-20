=head1 LICENSE

 Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
 Copyright [2016-2020] EMBL-European Bioinformatics Institute

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

=head1 CONTACT

 Please email comments or questions to the public Ensembl
 developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

 Questions may also be sent to the Ensembl help desk at
 <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Star

=head1 SYNOPSIS

  my $runnable =
    Bio::EnsEMBL::Analysis::Runnable::Star->new();

 $runnable->run;
 my @results = $runnable->output;

=head1 DESCRIPTION

This module uses Star to align fastq to a genomic sequence. Star is a splice aware
aligner. It creates output files with the reads overlapping splice sites and the reads
aligning on the exons. Some reads are aligned multiple times in the genome.

=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::Runnable::Star;

use warnings;
use strict;

use File::Spec::Functions qw(catdir catfile);
use File::Path qw(remove_tree);

use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use parent ('Bio::EnsEMBL::Analysis::Runnable::BaseShortReadAligner');


=head2 new

 Arg [DECOMPRESS]  : String, as a command like 'gzip -c -'
 Arg [SAMPLE_NAME] : String, specify the sample name, used for the temporary files and directories
 Arg [GENOME_DIR]  : String, path to a directory containing the index SA
 Arg [THREADS]     : Int, number of threads to use
 Description       : Creates a  object to align reads to a genome using STAR
 Returntype        : Bio::EnsEMBL::Analysis::Runnable::Star
 Exceptions        : Throws if the genome has not been indexed

=cut

sub new {
  my ( $class, @args ) = @_;

  my $self = $class->SUPER::new(@args);
  my ($decompress, $sample_id, $genome_dir, $threads) = rearrange([qw (DECOMPRESS SAMPLE_NAME GENOME_DIR THREADS)],@args);
  $self->throw("Genome file must be indexed, '$genome_dir/SA' does not exist\n") unless (-e $genome_dir.'/SA');
  $self->genome($genome_dir);
  $self->sample_id($sample_id);
  $self->threads($threads);
  $self->is_file_compressed($decompress);
  return $self;
}

=head2 run

 Arg [1]    : None
 Description: Run Star to align reads to an indexed genome. The resulting output file will be stored in $self->output
              It will delete the previous temporary files if they exist
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;

  my $sample_id = $self->sample_id();
  my $fastq = $self->fastq;
  my $fastqpair = $self->fastqpair;
  my $options = $self->options;
  my $out_dir = catfile($self->outdir, $sample_id.'_');
  my $tmp_dir = catdir($self->outdir, $sample_id.'_tmp');
  if (-d $tmp_dir) {
    remove_tree($tmp_dir, {safe => 1});
  }
  $options .= ' --readFilesCommand '.$self->is_file_compressed if ($self->is_file_compressed);
  $options .= ' --outSAMattrRGline "'.$self->rg_lines.'"' if ($self->rg_lines);
  $options .= ' --outSAMattributes "'.$self->sam_attributes.'"' if ($self->sam_attributes);

  # run STAR
  my $command = $self->program." --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMstrandField intronMotif --runThreadN ".$self->threads." --twopassMode Basic --runMode alignReads --genomeDir ".$self->genome." --readFilesIn $fastq $fastqpair --outFileNamePrefix $out_dir $options --outTmpDir $tmp_dir --outSAMtype BAM SortedByCoordinate";
  $self->warning("Command: $command\n");
  if (system($command)) {
      $self->throw("Error aligning $fastq $fastqpair\nCommandline used: $command\nError code: $?\n");
  }
  $self->output([$out_dir.'Aligned.sortedByCoord.out.bam']);
}




#Containers
#=================================================================


=head2 is_file_compressed

 Arg [1]    : (optional) String
 Description: Getter/setter for the command to execute with --readFilesCommand
 Returntype : String
 Exceptions : None

=cut

sub is_file_compressed {
    my ($self, $value) = @_;
    if (defined $value) {
      $self->{'_decompress'} = $value;
    }
    return $self->{'_decompress'};
}


=head2 rg_lines

 Arg [1]    : (optional) String
 Example    : $self->rg_lines("ID:blood\tCN:ENA\tPL:ILLUMINA");
 Description: getter/setter for the attributes expected in the RG line
 Returntype : String
 Exceptions : None

=cut

sub rg_lines {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_rg_lines'} = $value;
  }

  if (exists($self->{'_rg_lines'})) {
    return $self->{'_rg_lines'};
  } else {
    return;
  }
}


=head2 sam_attributes

 Arg [1]    : (optional) String
 Example    : $self->sam_attributes('NH HI AS NM MD');
 Description: getter/setter for the attributes expected in the SAM alignment
 Returntype : String
 Exceptions : None

=cut

sub sam_attributes {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_sam_attributes'} = $value;
  }

  if (exists($self->{'_sam_attributes'})) {
    return $self->{'_sam_attributes'};
  } else {
    return;
  }
}


=head2 sample_id

 Arg [1]    : (optional) String
 Description: Getter/setter for the sample name, which becomes the file name for the output
 Returntype : String
 Exceptions : None

=cut

sub sample_id {
  my ($self, $value) = @_;
  if (defined $value) {
    $self->{'_sample_id'} = $value;
  }
  return $self->{'_sample_id'};
}


=head2 threads

 Arg [1]    : (optional) int
 Description: Getter/setter for the number of threads to use
 Returntype : String
 Exceptions : None

=cut

sub threads {
  my ($self, $value) = @_;
  if (defined $value) {
    $self->{'_threads'} = $value;
  }
  return $self->{'_threads'};
}


1;
