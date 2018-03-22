=head1 LICENSE

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
 Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

use File::Spec;

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
  my ($decompress, $sam_attributes) = rearrange([qw (DECOMPRESS RG_LINES)],@args);
  $self->throw('Your working directory '.$self->workdir." does not exist!\n") unless (-d $self->workdir);
  my (undef,$genome_dir,) = File::Spec->splitpath($self->genome);
  $self->throw("Genome file must be indexed, '$genome_dir/SA' does not exist\n") unless (-e $genome_dir.'/SA');
  $self->genome($genome_dir);
  $self->sam_attributes($sam_attributes);
  $self->is_file_compressed($decompress || '');
  return $self;
}

=head2 run

 Arg [1]    : None
 Description: Run Star to align reads to an indexed genome. The resulting output file will be stored in $self->output
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;

  my $fastq = $self->fastq;
  my $fastqpair = $self->fastqpair;
  my $options = $self->options;
  $options .= ' --readFilesCommand '.$self->is_file_compressed if ($self->is_file_compressed);
  $options .= ' --outSAMattrRGline '.$self->sam_attributes;
  my (undef,undef,$filename) = File::Spec->splitpath($fastq);
  my $outdir = File::Spec->catfile($self->outdir, $filename);

  # run STAR
  my $command = $self->program." --twopassMode Basic --runMode alignReads --genomeDir ".$self->genome." --readFilesIn $fastq $fastqpair --outFileNamePrefix $outdir $options --outTmpDir ".$self->workdir;
  $self->warning("Command: $command\n");
  if (system($command)) {
      $self->throw("Error aligning $fastq $fastqpair\nError code: $?\n");
  }
  $self->output([$outdir.'.bam']);
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

1;
