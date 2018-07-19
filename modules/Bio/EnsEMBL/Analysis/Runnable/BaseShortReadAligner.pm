=head1 LICENSE

# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Runnable::BaseShortReadAligner

=head1 SYNOPSIS


=head1 DESCRIPTION

Modules which will align short reads to a genome file using software such as Star,
BWA, should inherits from this module. It defines the base options; fastq, fastqpair,
outdir and genome.

=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::Runnable::BaseShortReadAligner;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use parent ('Bio::EnsEMBL::Analysis::Runnable');


=head2 new

  Arg [FASTQ]     : String, absolute path of the FASTQ file containing mates 1
  Arg [FASTQPAIR] : String, (optional) absolute path of the FASTQ file containing mates 2
  Arg [OUTDIR]    : String, directory where you will store the results
  Arg [GENOME]    : String, absolute path of the genome FASTA file
  Description     : Constructor, fastqpair is not checked as you may use single end reads
  Returntype      : Throws if FASTQ, OUTDIR or GENOME are not defined

=cut

sub new {
  my ( $class, @args ) = @_;

  my $self = $class->SUPER::new(@args);
  my ($fastq, $fastqpair, $outdir, $genome ) = rearrange([qw (FASTQ FASTQPAIR OUTDIR GENOME)],@args);
  $self->throw("Your fastq file '$fastq' does not exist\n") unless (-e $fastq);
  $self->fastq($fastq);
  $self->fastqpair($fastqpair);
  $self->throw("Your output directory '$outdir' does not exist!\n") unless (-d $outdir);
  $self->outdir($outdir);
  $self->throw("Your genome file '$genome' does not exist!\n") unless (-e $genome);
  $self->genome($genome);
  return $self;
}


=head2 

 Arg [1]    : None
 Description: Run method that need to be implemented in the module inheriting from
              Bio::EnsEMBL::Analysis::Runnable::BaseShortReadAligner. It is encourage to store
              the resulting file(s) into $self->output
 Returntype : None
 Exceptions : Throws if the method has not been implemented in the inheriting module

=cut

sub run {
  my ($self) = @_;

  $self->throw('You should override the run method if you want to use Bio::EnsEMBL::Analysis::Runnable::BaseShortReadAligner');
}




#Containers
#=================================================================

=head2 fastq

 Arg [1]    : (optional) String
 Description: Getter/setter for the fastq file for the first mate
 Returntype : String
 Exceptions : None

=cut

sub fastq {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_fastq'} = $value;
  }

  if (exists($self->{'_fastq'})) {
    return $self->{'_fastq'};
  } else {
    return;
  }
}


=head2 fastqpair

 Arg [1]    : (optional) String
 Description: Getter/setter for the fastq file for the second mate
 Returntype : String
 Exceptions : None

=cut

sub fastqpair {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_fastqpair'} = $value;
  }

  if (exists($self->{'_fastqpair'})) {
    return $self->{'_fastqpair'};
  } else {
    return;
  }
}


=head2 outdir

 Arg [1]    : (optional) String
 Description: Getter/setter for the output directory
 Returntype : String
 Exceptions : None

=cut

sub outdir {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_outdir'} = $value;
  }

  if (exists($self->{'_outdir'})) {
    return $self->{'_outdir'};
  } else {
    return;
  }
}


=head2 genome

 Arg [1]    : (optional) String
 Description: Getter/setter for the genome file
 Returntype : String
 Exceptions : None

=cut

sub genome {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_genome'} = $value;
  }

  if (exists($self->{'_genome'})) {
    return $self->{'_genome'};
  } else {
    return;
  }
}

1;
