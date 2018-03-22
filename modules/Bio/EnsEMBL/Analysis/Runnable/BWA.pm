=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Runnable::BWA

=head1 SYNOPSIS

  my $runnable = 
    Bio::EnsEMBL::Analysis::Runnable::BWA->new();

 $runnable->run;
 my @results = $runnable->output;
 
=head1 DESCRIPTION

This module uses BWA to align fastq to a genomic sequence

=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::Runnable::BWA;

use warnings ;
use strict;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(execute_with_wait);

use parent ('Bio::EnsEMBL::Analysis::Runnable::BaseShortReadAligner');

=head2 new

 Arg [OPTIONS]: String
 Arg [FASTQ]  : String
 Arg [OUTDIR] : String
 Arg [GENOME] : String
 Description  : Creates a new Bio::EnsEMBL::Analysis::Runnable::BWA object to align short reads on a genome using BWA
 Returntype   : Bio::EnsEMBL::Analysis::Runnable::BWA
 Exceptions   : Throws if FASTQ, OPTIONS, OUTDIR or GENOME are not set
                Throws if GENOME has not been indexed

=cut

sub new {
  my ( $class, @args ) = @_;

  my $self = $class->SUPER::new(@args);
  $self->throw("Genome file must be indexed \ntry ".$self->program.' index '.$self->genome."\n") unless (-e $self->genome.'.ann');
  return $self;
}

=head2 run

  Args       : none
  Description: Run BWA to align reads to an indexed genome
  Returntype : none

=cut 

sub run {
  my ($self) = @_;

  my $fastq = $self->fastq;
  my $options = $self->options;
  my $outdir = $self->outdir;
  my $program = $self->program;
  my @tmp = split(/\//,$fastq);
  my $filename = pop @tmp;
  # run bwa
  my $command = "$program aln $options -f $outdir/$filename.sai ".$self->genome." $fastq";
  print STDERR "Command: $command\n";
  $self->warning("Command: $command\n");
  execute_with_wait($command);
}




#Containers
#=================================================================

=head2 fastq

 Arg [1]    : (optional) String
 Description: Getter/setter
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


=head2 options

 Arg [1]    : (optional) String
 Description: Getter/setter
 Returntype : String
 Exceptions : None

=cut

sub options {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_options'} = $value;
  }
  
  if (exists($self->{'_options'})) {
    return $self->{'_options'};
  } else {
    return;
  }
}


=head2 outdir

 Arg [1]    : (optional) String
 Description: Getter/setter
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
 Description: Getter/setter
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
