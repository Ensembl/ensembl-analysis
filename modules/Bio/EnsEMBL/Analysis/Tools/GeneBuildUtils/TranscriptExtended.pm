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

Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::TranscriptExtended; 

=head1 SYNOPSIS 


  Creation by re-blessing:

  for my $t ( @transcripts ) { 
    bless $t,"Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended";

    $t->ev_set('est') ; 
  }  

  Creation without re-blessing : 

  my $tran = new Bio::EnsEMBL::Transcript(-EXONS => \@exons);


=head1 DESCRIPTION

This module extends a Bio::EnsEMBL::Transcript with smoe other
methods. You have to re-bless your Bio::EnsEMBL::Transcript objects 
to use it. 

=head1 CONTACT

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

=cut



package Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended; 

use warnings ;
use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Feature ; 
use Bio::EnsEMBL::Transcript ; 
use Bio::EnsEMBL::Utils::Exception qw(throw warning );

@ISA = qw(Bio::EnsEMBL::Transcript Bio::EnsEMBL::Feature);


sub new  {
  my($class) = shift;
  if( ref $class ) { 
      $class = ref $class;
  }
  my $self = $class->SUPER::new(@_);
  return $self ;
   
  $self->{different_est_support} = {} ;
} 

=head

   Name : ev_set
   Arg  : String
   Func : getter/setter for the evidence-set 
   Return : String describing ev_set (defined in GeneBuild/TrancriptCoalescer.pm)

=cut

sub ev_set {
  my ($self,$value) = @_;
  if (defined($value)) {
    $self->{'set'} = $value;
  }
  return $self->{'set'};
}


=head

   Name : remove_Exon 
   Arg  : Bio::EnsEMBL::Exon 
   Func : removes exon from Bio::EnsEMBL::Transcript 
   Return : none 
   Example : $transcript->remove_Exon($exon) ; 
=cut



sub remove_Exon { 
  my ($self,$ex_to_remove ) = @_;
  unless ($ex_to_remove){ 
    throw ("you have to supply an Bio::EnsEMBL::Exon object ".
           "which will be removed from the transcript\n");
  }
  my $clone=[] ; 
  my $ea = $self->{'_trans_exon_array'}; 
  my $nr_before = scalar(@$ea); 
  for my $e (@$ea) {
    if  ($e ne $ex_to_remove ) { 
       push @$clone, $e; 
    }
   } 
  $self->{'_trans_exon_array'} = $clone;
  if ($nr_before == scalar(@$clone) ) { 
    warning( "No matching exon found, exon could not be removed\n" );
   }
  return ; 
}




sub has_3prim_support {
  my ($self,$value) = @_;
  if (defined($value)) {
    $self->{'has_3prim_sup'} = $value;
  }
  return $self->{'has_3prim_sup'};
}


sub has_5prim_support {
  my ($self,$value) = @_;
  if (defined($value)) {
    $self->{'has_5prim_sup'} = $value;
  }
  return $self->{'has_5prim_sup'};
}




sub extend_3prim_end {
  my ($self,$value) = @_;
  if (defined($value)) {
    $self->{'extend_3prim_end'} = $value;
  }
  return $self->{'extend_3prim_end'};
}



sub extend_5prim_end {
  my ($self,$value) = @_;
  if (defined($value)) {
    $self->{'extend_5prim_end'} = $value;
  }
  return $self->{'extend_5prim_end'};
}



sub nr_exons_overlapped_by_est {
  my ($self,$value) = @_;
  if (defined($value)) {
    $self->{'nr_exons_overlapped_est'}+= $value;
  }
  return $self->{'nr_exons_overlapped_est'};
} 

=head  

Name : different_est_support($transcript) 
Arg  : Bio::EnsEMBL::Transcript 
Func : checks if $transcript is not overlapping the other transcripts,
       if there is no overlap between already stored transcripts this transcript is added 

=cut

 
sub different_est_support {
  my ($self,$new_tr) = @_;

  my $overlap ;  
  if ($new_tr) {  
    for my $key  ( keys %{ $self->{different_est_support} } ) {  
      my $stored_tr = $self->{different_est_support}{$key} ;  
  
       if ($stored_tr->seq_region_start <= $new_tr->seq_region_end  && 
           $new_tr->seq_region_start <= $stored_tr->seq_region_end  ) { 
           # new_tr overlaps stored_tr 
  
         # check which transcript has more exons 
         my $new_tr_exons = scalar ( @{  $new_tr->get_all_Exons } )  ; 
         my $stored_tr_exons = scalar(@{$self->{different_est_support}{$stored_tr}->get_all_Exons }); 
         if ($new_tr_exons > $stored_tr_exons ) {  
            delete $self->{different_est_support}{$stored_tr} ;
            if (exists $self->{different_est_support}{$stored_tr}) { 
              throw(" key wasn't delted\n" ) ;  
            }else { 
              warn(" key is deleted\n" ) ;
            } 
            $self->{different_est_support}{$new_tr} = $new_tr ;
         }  
       }else {
         ${$self->{different_est_support}}{$new_tr}=$new_tr ; 
       }
    }
  } 
  return $self->{different_est_support} ;   
} 





=head  

Name   : exchange_exon
Arg[1] : Bio::EnsEMBL::Exon
Arg[2] : Bio::EnsEMBL::Exon
Func   : splices a new exon in a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended-object 
         and cut's ond the old Bio::EnsEMBL::Exon object 
Returnval : Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended()

=cut



sub exchange_exon {
  my ($self, $old_exon, $new_exon  ) = @_ ; 

  my $biotype = $self->biotype ;  
  my @exons = @{ $self->get_all_Exons } ;  
  my $new_tr = Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended->new(); 

  for my $e  (@exons) {
    if ($e eq $old_exon) {
     $e = $new_exon ;
    }
    $new_tr->add_Exon($e) ; 
  }
  $new_tr->biotype($biotype) ; 
  return $new_tr ; 
}

=head  

Name      : score
Arg[1]    : Scalar score
Func      : Get / Set a score value for the transcript
Returnval : Scalar score

=cut

sub score {
  my ($self,$value) = @_;
  if (defined($value)) {
    $self->{'score'} = $value;
  }
  return $self->{'score'};
}

=head  

Name      : similarity
Arg[1]    : Bio::EnsEMBL::Transcript
Func      : Get / Set a similarity transcript associated with this transcript
Returnval : Bio::EnsEMBL::Transcript

=cut

sub similarity {
  my ($self,$value) = @_;
  if (defined($value)) {
    $self->{'similarity'} = $value;
  }
  return $self->{'similarity'};
}



1;
