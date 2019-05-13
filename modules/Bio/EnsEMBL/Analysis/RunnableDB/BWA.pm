# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::RunnableDB::BWA;

use warnings ;
use strict;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Config::GeneBuild::BWA;
use Bio::EnsEMBL::Analysis::Runnable::BWA;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use vars qw(@ISA);

@ISA =  ("Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild");


sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
  $self->read_and_check_config($BWA_CONFIG_BY_LOGIC); 
  return $self;
}

sub fetch_input {
  my ($self) = @_;
  my %parameters = %{$self->parameters_hash};
  my $program = $self->analysis->program_file;
  my $filename =  $self->INDIR ."/" .$self->input_id;
  $self->throw("Fastq file  $filename not found\n")
    unless ( -e $filename );
  $self->throw("BWA program not defined in analysis \n")
    if not defined $program;
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::BWA->new
    (
     -analysis => $self->analysis,
     -program  => $program,
     -options  => $self->OPTIONS,
     -outdir   => $self->OUTDIR,
     -genome   => $self->GENOMEFILE,
     -fastq    => $filename,
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


sub OUTDIR {
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



sub GENOMEFILE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_GENOMEFILE'} = $value;
  }
  
  if (exists($self->{'_GENOMEFILE'})) {
    return $self->{'_GENOMEFILE'};
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

sub SAMPE_OPTIONS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_SAMPE_OPTIONS'} = $value;
  }
  
  if (exists($self->{'_SAMPE_OPTIONS'})) {
    return $self->{'_SAMPE_OPTIONS'};
  } else {
    return undef;
  }
}


sub SAMSE_OPTIONS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_SAMSE_OPTIONS'} = $value;
  }
  
  if (exists($self->{'_SAMSE_OPTIONS'})) {
    return $self->{'_SAMSE_OPTIONS'};
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

sub MIN_MAPPED {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_MIN_MAPPED'} = $value;
  }

  if (exists($self->{'_MIN_MAPPED'})) {
    return $self->{'_MIN_MAPPED'};
  } else {
    return undef;
  }
}

sub MIN_PAIRED {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_MIN_PAIRED'} = $value;
  }

  if (exists($self->{'_MIN_PAIRED'})) {
    return $self->{'_MIN_PAIRED'};
  } else {
    return undef;
  }
}

1;
