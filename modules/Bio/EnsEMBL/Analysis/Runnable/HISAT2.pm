=head1 LICENSE

 Copyright [2020] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Runnable::HISAT2

=head1 SYNOPSIS

  my $runnable =
    Bio::EnsEMBL::Analysis::Runnable::HISAT2->new();

 $runnable->run;
 my @results = $runnable->output;

=head1 DESCRIPTION

This module uses HISAT2 to align fastq to a genomic sequence. HISAT2 is a splice aware
aligner. It creates output files with the reads overlapping splice sites and the reads
aligning on the exons. Some reads are aligned multiple times in the genome.

=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::Runnable::HISAT2;

use warnings;
use strict;

use File::Spec::Functions;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use parent ('Bio::EnsEMBL::Analysis::Runnable::BaseShortReadAligner');


=head2 new

 Arg [DECOMPRESS]           : String as a command like 'gzip -c -'
 Arg [EXPECTED_ATTRIBUTES]  : String specify the attribute expected for the output, see STAR manual
 Description                : Creates a  object to align reads to a genome using STAR
 Returntype                 : 
 Exceptions                 : Throws if WORKDIR does not exist
                              Throws if the genome has not been indexed

=cut

sub new {
  my ( $class, @args ) = @_;

  my $self = $class->SUPER::new(@args);
  my ($decompress, $sam_attributes, $sample_id, $genome_index, $genome, $threads) = rearrange([qw (DECOMPRESS RG_LINES SAMPLE_NAME GENOME_INDEX GENOME THREADS)],@args);
#  $self->throw("Genome file must be indexed, index specified does not exist. Index specified:\n".$genome) unless (-e $genome.'.1.ht2');
  $self->genome($genome);
  $self->genome_index($genome_index);
  $self->sample_id($sample_id);
  $self->threads($threads);
  return $self;
}

=head2 run

 Arg [1]    : None
 Description: Run HISAT2 to align reads to an indexed genome. The resulting output file will be stored in $self->output
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;

  my $sample_id = $self->sample_id();
  my $samtools_path = $self->samtools_path();
  my $fastq = $self->fastq;
  my $fastqpair = $self->fastqpair;
  my $options = $self->options;
  my $threads = $self->threads();
  my $out_dir = catfile($self->outdir, $sample_id);
  my $tmp_dir = catfile($self->outdir, $sample_id."_tmp");

  my $hisat2_path = $self->program();
  unless(-e $hisat2_path) {
    $self->warning("HISAT2 not found on the specified path, or path not specified. Defaulting to 'hisat2'");
    $hisat2_path = "hisat2";
  }

  # run HISAT2
  my $command = "";
  my $sam_file = $out_dir.".sam";
  if($fastqpair) {
    $command =  $hisat2_path." -q -x ".$self->genome_index." -1 ".$fastq." -2 ".$fastqpair." -S ".$sam_file;
  } else {
    $command =  $hisat2_path." -q -x ".$self->genome_index." -U ".$fastq." -S ".$sam_file;
  }

  $self->warning("Command: $command\n");
  if (system($command)) {
    $self->throw("Error aligning $fastq $fastqpair\nCommandline used: $command\nError code: $?\n");
  }

  # Create sorted bam file
  my $samtools_command = "";
  unless(-e $samtools_path) {
    $self->warning("Samtools not found on the specified path, or path not specified. Defaulting to 'samtools'");
    $samtools_path = "samtools";
  }

  my $bam_file = $out_dir.'.bam';
  $samtools_command = $samtools_path." sort -o ".$bam_file." ".$sam_file;
  if(system($samtools_command)) {
    $self->throw("Error creating sorted bam file\nCommandline used: $samtools_command\nError code: $?\n");
  }

  $self->output([$bam_file]);
}




#Containers
#=================================================================

=head2 samtools_path

 Arg [1]    : (optional) String
 Description: Getter/setter for the path to samtools
 Returntype : String
 Exceptions : None

=cut

sub samtools_path {
    my ($self, $value) = @_;
    if (defined $value) {
      $self->{'_samtools_path'} = $value;
    }
    return $self->{'_samtools_path'};
}


=head2 genome_index

 Arg [1]    : (optional) String
 Description: Getter/setter for the command to genome index. This is used instead of the parent
              genome subroutine as the hisat2 index consists of multiple files with different extensions
 Returntype : String
 Exceptions : None

=cut

sub genome_index {
    my ($self, $value) = @_;
    if (defined $value) {
      $self->{'_genome_index'} = $value;
    }
    return $self->{'_genome_index'};
}


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
