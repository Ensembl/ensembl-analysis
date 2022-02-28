
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

Bio::EnsEMBL::Analysis::RunnableDB::Gsnap

=head1 SYNOPSIS

my $runnableDB =  Bio::EnsEMBL::Analysis::RunnableDB::Gsnap->new( );

$runnableDB->fetch_input();
$runnableDB->run();

=head1 DESCRIPTION

This module uses Gsnap to align fastq to a genomic sequence

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Gsnap;

use warnings ;
use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Config::GeneBuild::Gsnap;
use Bio::EnsEMBL::Analysis::Runnable::Gsnap;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use vars qw(@ISA);

@ISA =  ("Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild");


sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
  $self->read_and_check_config($GSNAP_CONFIG_BY_LOGIC); 
  return $self;
}

sub fetch_input {
  my ($self) = @_;
  my %parameters = %{$self->parameters_hash};
  my $program = $self->analysis->program_file;
  my $filename =  $self->INDIR ."/" .$self->input_id;

  $self->throw("Gsnap program not defined in analysis \n")
    if not defined $program;
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Gsnap->new
    (
     -analysis => $self->analysis,
     -program  => $program,
     -options  => $self->OPTIONS,
     -indir    => $self->INDIR,
     -outdir   => $self->OUTDIR,
     -genome   => $self->GENOMEDIR,
     -genomename  => $self->GENOMENAME,
     -fastq    => $filename,
     -paired   => $self->PAIRED,
     -samtools => $self->SAMTOOLS_PATH,
     -header   => $self->HEADER,
     %parameters,
    );
  $self->runnable($runnable);
}


sub run {
  my ($self) = @_;
  $self->throw("Can't run - no runnable objects") unless ( $self->runnable );
  my ($runnable) = @{$self->runnable};
  eval {
    $runnable->run;
    }; 
    if(my $err = $@){
      chomp $err;
	$self->throw("ERROR $err \n");
      }
}

# override write output as we have nothing for the db
sub write_output {
  my ($self) = @_;
}

#Containers
#=================================================================


sub OUTDIR {
  my $total_reads = 0;
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_OUTDIR'} = $value;
  }
  
  if (exists($self->{'_OUTDIR'})) {
    return $self->{'_OUTDIR'};
  } else {
    return undef;
  }
}


sub INDIR {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_INDIR'} = $value;
  }
  
  if (exists($self->{'_INDIR'})) {
    return $self->{'_INDIR'};
  } else {
    return undef;
  }
}


sub GENOMEDIR {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_GENOMEDIR'} = $value;
  }
  
  if (exists($self->{'_GENOMEDIR'})) {
    return $self->{'_GENOMEDIR'};
  } else {
    return undef;
  }
}

sub GENOMENAME {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_GENOMENAME'} = $value;
  }
  
  if (exists($self->{'_GENOMENAME'})) {
    return $self->{'_GENOMENAME'};
  } else {
    return undef;
  }
}

sub OPTIONS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_OPTIONS'} = $value;
  }
  
  if (exists($self->{'_OPTIONS'})) {
    return $self->{'_OPTIONS'};
  } else {
    return undef;
  }
}


sub PAIRED {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_PAIRED'} = $value;
  }
  
  if (exists($self->{'_PAIRED'})) {
    return $self->{'_PAIRED'};
  } else {
    return undef;
  }
}


sub SAMTOOLS_PATH {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_SAMTOOLS_PATH'} = $value;
  }
  
  if (exists($self->{'_SAMTOOLS_PATH'})) {
    return $self->{'_SAMTOOLS_PATH'};
  } else {
    return undef;
  }
}

sub HEADER {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_HEADER'} = $value;
  }
  
  if (exists($self->{'_HEADER'})) {
    return $self->{'_HEADER'};
  } else {
    return undef;
  }
}
