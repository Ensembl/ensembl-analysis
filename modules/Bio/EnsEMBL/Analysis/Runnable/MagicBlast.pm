=head1 LICENSE

 Copyright [2019] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Runnable::MagicBlast

=head1 SYNOPSIS

  my $runnable =
    Bio::EnsEMBL::Analysis::Runnable::MagicBlast->new();

 $runnable->run;
 my @results = $runnable->output;

=head1 DESCRIPTION

This module uses Star to align fastq to a genomic sequence. Star is a splice aware
aligner. It creates output files with the reads overlapping splice sites and the reads
aligning on the exons. Some reads are aligned multiple times in the genome.

=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::Runnable::MagicBlast;

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
  my ($decompress, $sam_attributes, $genome, $sample_id, $threads) = rearrange([qw (DECOMPRESS RG_LINES GENOME SAMPLE_NAME THREADS)],@args);
  $self->genome($genome);
  $self->sample_id($sample_id);
  $self->threads($threads);
  return $self;
}

=head2 run

 Arg [1]    : None
 Description: Run MagicBlast to align reads to an indexed genome. The resulting output file will be stored in $self->output
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;

  my $sample_id = $self->sample_id();
  my $fastq = $self->fastq;
  my $fastqpair = $self->fastqpair;
  my $options = $self->options;
  my $threads = $self->threads();
  my $out_file = catfile($self->outdir, $sample_id.".sam");

  # run magicblast
  my $command = $self->program." -num_threads ".$threads." -paired -infmt fastq -db ".$self->genome." -outfmt sam -query ".$fastq." -query_mate ".$fastqpair." -out ".$out_file;

  $self->warning("Command: $command\n");
  if (system($command)) {
      $self->throw("Error aligning $fastq $fastqpair\nCommandline used: $command\nError code: $?\n");
  }
  $self->output([$out_file]);
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
