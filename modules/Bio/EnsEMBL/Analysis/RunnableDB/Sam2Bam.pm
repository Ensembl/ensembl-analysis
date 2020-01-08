# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::Sam2Bam




=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::RunnableDB::Sam2Bam->new( );

$runnableDB->fetch_input();
$runnableDB->run();

=head1 DESCRIPTION

This module uses samtools to convert a directory containing SAM
files into a single sorted indexed merged BAM file

=head1 CONTACT

Post general queries to B<http://lists.ensembl.org/mailman/listinfo/dev>

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Sam2Bam;

use warnings ;
use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Config::GeneBuild::Sam2Bam;
use Bio::EnsEMBL::Analysis::Runnable::Sam2Bam;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use vars qw(@ISA);

@ISA =  ("Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild");


sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
  $self->read_and_check_config($SAM2BAM_CONFIG_BY_LOGIC); 
  return $self;
}

sub fetch_input {
  my ($self) = @_;
  my %parameters = %{$self->parameters_hash};
  my $program = $self->analysis->program_file;
  $self->throw("Samtools program not defined in analysis \n")
    if not defined $program;
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Sam2Bam->new
    (
     -analysis => $self->analysis,
     -header   => $self->HEADERFILE,
     -program  => $program,
     -regex    => $self->REGEX,
     -samdir   => $self->SAM_DIR,
     -bamfile  => $self->BAMFILE,
     -genome   => $self->GENOMEFILE,
     %parameters,
    ); 
  $self->runnable($runnable);
}


sub run {
  my ($self) = @_;
  $self->throw("Can't run - no runnable objects") unless ( $self->runnable );
  my ($runnable) = @{$self->runnable};
  $runnable->run;
}

# override write output as we have nothing for the db
sub write_output {
  my ($self) = @_;
}

#Containers
#=================================================================

sub SAM_DIR {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_SAM_DIR'} = $value;
  }
  
  if (exists($self->{'_SAM_DIR'})) {
    return $self->{'_SAM_DIR'};
  } else {
    return undef;
  }
}

sub REGEX {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_REGEX'} = $value;
  }
  
  if (exists($self->{'_REGEX'})) {
    return $self->{'_REGEX'};
  } else {
    return undef;
  }
}

sub BAMFILE {
  my ($self,$value) = @_;

  if (defined $value) {
    if ( $value =~ /(.+)\..+/){
      $value = $1;
    }
    $self->{'_BAMFILE'} = $value;
  }
  
  if (exists($self->{'_BAMFILE'})) {
    return $self->{'_BAMFILE'};
  } else {
    return undef;
  }
}

sub GENOMEFILE {
  my ($self,$value) = @_;

  if (defined $value) {
    if ( $value =~ /(.+)\..+/){
    }
    $self->{'_GENOMEFILE'} = $value;
  }
  
  if (exists($self->{'_GENOMEFILE'})) {
    return $self->{'_GENOMEFILE'};
  } else {
    return undef;
  }
}

sub HEADERFILE {
  my ($self,$value) = @_;

  if (defined $value) {
    if ( $value =~ /(.+)\..+/){
    }
    $self->{'_HEADERFILE'} = $value;
  }
  
  if (exists($self->{'_HEADERFILE'})) {
    return $self->{'_HEADERFILE'};
  } else {
    return undef;
  }
}
