# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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
=pod

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::BWA




=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::RunnableDB::BWA->new( );

$runnableDB->fetch_input();
$runnableDB->run();

=head1 DESCRIPTION

This module uses BWA to align fastq to a genomic sequence

=head1 CONTACT

Post general queries to B<http://lists.ensembl.org/mailman/listinfo/dev>

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::BWA2BAM;

use warnings ;
use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Runnable::BWA2BAM;
use Bio::EnsEMBL::Analysis::RunnableDB::BWA;
use vars qw(@ISA);

@ISA =  ("Bio::EnsEMBL::Analysis::RunnableDB::BWA");


sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
  return $self;
}

sub fetch_input {
  my ($self) = @_;
  my %parameters = %{$self->parameters_hash};
  my $program = $self->analysis->program_file;
  $self->throw("BWA program not defined in analysis \n")
    if not defined $program;
  my $filename =  $self->INDIR ."/" .$self->input_id;
  my $fastqpair;

  my $method;
  if ( $self->PAIRED ) {
    $method = " sampe " . $self->SAMPE_OPTIONS;
    # need to seaprate the 2 file names in the input id
    my @files = split(/:/,$self->input_id);
    $self->throw("Cannot parse 2 file names out of $self->input_id should be separated by :\n")
      unless scalar(@files) == 2;
    $filename  =  $self->INDIR ."/" .$files[0];
    $fastqpair =  $self->INDIR ."/" .$files[1];
  } else {
    $method = " samse " . $self->SAMSE_OPTIONS;
  }
  $self->throw("Fastq file  $filename not found\n")
    unless ( -e $filename );
  if ( $self->PAIRED) {
    $self->throw("Fastq pair file  $fastqpair not found\n")
      unless ( -e $fastqpair );
  }
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::BWA2BAM->new
    (
     -analysis   => $self->analysis,
     -program    => $program,
     -fastq      => $filename,
     -fastqpair  => $fastqpair,
     -options    => $self->OPTIONS,
     -outdir     => $self->OUTDIR,
     -genome     => $self->GENOMEFILE,
     -method     => $method,
     -samtools   => $self->SAMTOOLS_PATH,
     -header     => $self->HEADER,
     -min_mapped => $self->MIN_MAPPED,
     -min_paired => $self->MIN_PAIRED,
     %parameters,
    ); 
  $self->runnable($runnable);
}
