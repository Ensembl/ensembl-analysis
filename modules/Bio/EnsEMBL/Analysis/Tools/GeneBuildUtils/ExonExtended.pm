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
#

=head1 NAME

Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended; 

=head1 SYNOPSIS 


  Creation by re-blessing:

  for my $t ( @transcripts ) { 
    bless $t,"Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended";

    $t->ev_set('est') ; 
  }  

  Creation without re-blessing : 

  my $tran = new Bio::EnsEMBL::Transcript(-EXONS => \@exons);


=head1 SYNOPSIS 

  Creation by re-blessing:

  for my $t ( @transcripts ) { 
    bless $t,"Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended";
    $t->ev_set('est') ; 
  }  

  Creation without re-blessing : 

  my $exon= new Bio::EnsEMBL::Exon(
                                  -START  =>   100 , 
                                  -END    =>   300,
                                  -STRAND => '-1', 
                                  -SLICE  => $slice,
                                  -ANALYSIS => $analysis, 
                                  ); 

=head1 DESCRIPTION

This module extends a Bio::EnsEMBL::Transcript with smoe other
methods. You have to re-bless your Bio::EnsEMBL::Transcript objects 
to use it. 

=head1 CONTACT

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

=cut



package Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended; 

use warnings ;
use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Feature ; 
use Bio::EnsEMBL::Exon ; 
use Bio::EnsEMBL::Utils::Exception qw(throw warning );

@ISA = qw(Bio::EnsEMBL::Exon Bio::EnsEMBL::Feature);



sub new  {
  my($class) = shift;
  if( ref $class ) {
      $class = ref $class;
  }
  my $self = $class->SUPER::new(@_);
  return $self ;
}


=head2  biotype 

   Name       : biotype 
   Arg        : String
   Func       : getter/setter for the biotype to which this Exon belongs to 
   Returntype : String describing biotype 

=cut 


sub biotype {
  my ($self,$value) = @_;
  if (defined($value)) {
    $self->{'type'} = $value;
  }
  return $self->{'type'};
}


=head2 ev_set

   Name       : ev_set
   Arg        : String
   Func       : getter/setter for the evidence-set 
   Returntype : String describing ev_set (defined in GeneBuild/TrancriptCoalescer.pm)

=cut 


sub ev_set {
  my ($self,$value) = @_;
  if (defined($value)) {
    $self->{'set'} = $value;
  }
  return $self->{'set'};
}


=head2 transcript

   Name       : transcript 
   Arg        : Bio::EnsEMBL::Transcript
   Func       : getter/setter for the Bio::EnsEMBL::Transcript the Exon belongs to 
   Returntype : Bio::EnsEMBL::Transcript

=cut 


sub transcript {
  my ($self,$value) = @_;
  if (defined($value)) {
    $self->{'transcript'} = $value;
  }
  return $self->{'transcript'};
}

=head1

   Name       : number_exons 
   Arg        : int
   Func       : getter/setter number of exons in Transcript
   Returntype : int

=cut 


sub number_exons {
 my ($self) = @_ ; 
  if (defined($self->transcript)){ 
    return scalar( @{ $self->transcript->get_all_Exons} ) ; 
  }else {
    warning("Exon has no Bio::EnsEMBL::Transript-object attached - can't get number of exons\n" ) ; 
  }
  return undef ; 
}

=head2 prev_exon

   Name       : prev_exon 
   Arg        : Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended; 
   Func       : points to previous exon (5'prim) 
   Returntype : Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended; 

=cut 

sub prev_exon {
  my ($self,$value) = @_;
  if (defined($value)) {
    $self->{'prev_exon'} = $value;
  }
  return $self->{'prev_exon'};
}


=head2 next_exon

   Name       : next_exon  
   Arg        : Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended; 
   Func       : points to next exon (3'prim) 
   Returntype : Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended; 

=cut 


sub next_exon {
  my ($self,$value) = @_;
  if (defined($value)) {
    $self->{'next_exon'} = $value;
  }
  return $self->{'next_exon'};
}

=head2 is_terminal_exon

   Name       : is_terminal_exon 
   Arg        : none 
   Func       : returns 1 if exon is terminal and 0 othewise 
   Returntype : int (0 || 1) 

=cut 

sub is_terminal_exon {
  my ($self) = @_;
  if ($self->next_exon && $self->prev_exon) {
    return 0 ;
  }
  return  1 ;
}

=head2 is_3prim_exon

   Name       : is_3prim_exon 
   Arg        : none
   Func       : returns 1 if exon is at 3prim end of transcript 
   Returntype : int (0 || 1) 

=cut 

sub is_3prim_exon {
  my ($self) = @_;
  return 1 if ($self->prev_exon && !$self->next_exon) ; 
  return 0 ;
}

=head2 is_5prim_exon

   Name       : is_5prim_exon 
   Arg        : none 
   Func       : returns 1 if exon is at 5prim end of transcript 
   Returntype : int (0 || 1) 

=cut 

sub is_5prim_exon {
  my ($self) = @_;
  return 1 if ($self->next_exon && !$self->prev_exon) ; 
  return 0 ;   
}

=head2 cluster

   Name       : cluster
   Arg        : Bio::EnsEMBL::Transcript
   Func       : getter/setter for the Bio::EnsEMBL::Transcript the Exon belongs to 
   Returntype : Bio::EnsEMBL::Transcript

=cut 

sub cluster {
  my ($self,$value) = @_;
  if (defined($value)) {
    $self->{'cluster'} = $value;
  }
  return $self->{'cluster'};
}


=head2 visited ( $val ) 

  Function   : marks if exon has been visited in recursion procedure or not 
  Arg        : integer (true / false ) 
  Returntype : integer (true || false ) 
  Caller     : recursion procedure in Condense_EST.pm

=cut
 
sub visited  {
  my ($self,$value) = @_;
  if (defined($value)) {
    $self->{'visited'} = $value;
  }
  return $self->{'visited'};
}

  
=head2 get_percentage_exon_conversation_in_exon_cluster() 

  Example    : $self->get_percentage_exon_conversation_in_exon_cluster
  Arg        : none
  Function   : returns the precentage of other Exons in the same cluster which 
               share the same Exon-boundaries. 
  Returntype : float 

=cut
 


sub get_percentage_exon_conversation_in_exon_cluster {
  my ( $self ) = @_ ;

  my @ex_clust = @{ $self->cluster->get_all_Exons_in_ExonCluster } ;

  # uniquify exon acc. to their hashkey 
  my %uniq_exons ;
  for (@ex_clust) {
    unless ($_->is_terminal_exon) {
      $uniq_exons { $_->hashkey }++ ;
    }
  }

  # get maximum value (most consered boundaries) out of hash
  my @tmp = reverse sort (values %uniq_exons ) ;

  my $max_cons = -1;
  if (scalar(@tmp) > 0) {
    $max_cons = shift @tmp ;
  }

  my $percentage_exon_conservation = 0 ;
  if  ($max_cons > 0  ) {
    if ($uniq_exons{$self->hashkey} ){
     $percentage_exon_conservation = $uniq_exons{$self->hashkey} / $max_cons ;
    }
  } else {
   $percentage_exon_conservation = 0 ;
  }
  return $percentage_exon_conservation ;
}



