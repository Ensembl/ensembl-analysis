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

=head1 NAME

Bio::EnsEMBL::Analysis::PmatchFeature

=head1 SYNOPSIS

    

=head1 DESCRIPTION


=head1 CONTACT

http://lists.ensembl.org/mailman/listinfo/dev

=head1 APPENDIX


=cut


# Let the code begin...

package Bio::EnsEMBL::Pipeline::Tools::Pmatch::PmatchFeature;

use warnings ;
use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

@ISA = qw();

sub new {
  my($class,@args) = @_;

  my $self = bless {},$class;
  my ($protein_id,$cdna_id,$chr_name,$start,$end,$coverage,$analysis) = rearrange([qw(PROTEIN_ID 
                       CDNA_ID                                                               CHR_NAME
										 START
										 END
										 COVERAGE
										 ANALYSIS	      
										)],@args);
  
  
  $self->protein_id($protein_id);
  $self->chr_name($chr_name);
  $self->start($start);
  $self->end($end);
  $self->cdna_id($cdna_id);
  $self->coverage($coverage);
  $self->analysis($analysis);
  return $self;
}

sub  cdna_id {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_cdna_id} = $arg;
  }

  return $self->{_cdna_id};
}

sub coverage {
  my ($self,$arg) = @_;

  if (defined($arg)) {

    if ($arg < 0 || $arg > 100) {
      throw("Coverage must be betwee 0 and 100.  Trying to set it to $arg");
    }
    
    $self->{_coverage} = $arg;
  }

  return $self->{_coverage};
}

sub score { 
  my ($self) = @_ ; 
  return $self->{_coverage};
}


sub protein_id {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_protein_id} = $arg;
  }

  return $self->{_protein_id};
}

sub chr_name{
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_chr_name} = $arg;
  }
  return $self->{_chr_name};
}

sub start {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_start} = $arg;
  }
  return $self->{_start};
}

sub end {
  my ($self,$arg) = @_;

  if (defined($arg)) {
    $self->{_end} = $arg;
  }
  return $self->{_end};
}

sub analysis {
   my ($self,$value) = @_;

   if(!$self->{_analysis}){
     $self->{_analysis} = undef;
   }
   if ($value) {
     unless($value->isa('Bio::EnsEMBL::Analysis')) {
       throw("Analysis is not a Bio::EnsEMBL::Analysis object "
		    . "but a $value object");
     }

     $self->{_analysis} = $value;
   } 

   return $self->{_analysis};
}


1;

