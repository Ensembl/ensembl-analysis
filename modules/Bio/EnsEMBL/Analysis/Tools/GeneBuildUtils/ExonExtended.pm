
package Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended; 

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Feature ; 
use Bio::EnsEMBL::Exon ; 
use Bio::EnsEMBL::Utils::Exception qw(throw warning );

@ISA = qw(Bio::EnsEMBL::Exon Bio::EnsEMBL::Feature);




=head

   Name : biotype 
   Arg  : String
   Func : getter/setter for the biotype to which this Exon belongs to 
   Return : String describing biotype 

=cut 


sub biotype {
  my ($self,$value) = @_;
  if (defined($value)) {
    $self->{'type'} = $value;
  }
  return $self->{'type'};
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

#=head
#
#   Name : source_transcript
#   Func : returns the orignial transcript to which this exon belongs to 
#   Return : Bio::EnsEMBL::Transcript object 
#
#=cut 
#
#
#sub source_transcript {
#  my ($self,$value) = @_;
#  if (defined($value)) {
#    $self->{'source_transcript'} = $value;
#  }
#  return $self->{'source_transcript'};
#  
#}


=head

   Name : transcript 
   Arg  : Bio::EnsEMBL::Transcript
   Func : getter/setter for the Bio::EnsEMBL::Transcript the Exon belongs to 
   Return : Bio::EnsEMBL::Transcript

=cut 


sub transcript {
  my ($self,$value) = @_;
  if (defined($value)) {
    $self->{'transcript'} = $value;
  }
  return $self->{'transcript'};
}


#=head
#
#   Name : exon_is_overlapped 
#   Arg  : Bio::EnsEMBL::Transcript
#   Func : getter/setter for the Bio::EnsEMBL::Transcript the Exon belongs to 
#   Return : Bio::EnsEMBL::Transcript
#
#=cut 
#
#
#sub exon_is_overlapped {
#  my ($self,$value) = @_;
#  if (defined($value)) {
#    $self->{'exon_is_overlapped'} = $value;
#  }
#  return $self->{'exon_is_overlapped'};
#}
#
#=head
#
#   Name : number_exons 
#   Arg  : int
#   Func : getter/setter number of exons in Transcript
#   Return : int
#
#=cut 
#
#
#sub number_exons {
#  my ($self,$value) = @_;
#  if (defined($value)) {
#    $self->{'number'} = $value;
#  }
#  return $self->{'number'};
#}

=head

   Name : prev_exon 
   Arg  : Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended; 
   Func : points to previous exon (5'prim) 
   Return : Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended; 

=cut 

sub prev_exon {
  my ($self,$value) = @_;
  if (defined($value)) {
    $self->{'prev_exon'} = $value;
  }
  return $self->{'prev_exon'};
}


=head

   Name : next_exon  
   Arg  : Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended; 
   Func : points to next exon (3'prim) 
   Return : Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended; 

=cut 


sub next_exon {
  my ($self,$value) = @_;
  if (defined($value)) {
    $self->{'next_exon'} = $value;
  }
  return $self->{'next_exon'};
}

=head

   Name : is_terminal_exon 
   Arg  : none 
   Func : returns 1 if exon is terminal and 0 othewise 
   Return : int (0 || 1) 

=cut 

sub is_terminal_exon {
  my ($self) = @_;
  if ($self->next_exon && $self->prev_exon) {
    return 0 ; 
  } 
  return  1 ; 
}

=head

   Name : is_3prim_exon 
   Arg  : none
   Func : returns 1 if exon is at 3prim end of transcript 
   Return : int (0 || 1) 

=cut 

sub is_3prim_exon {
  my ($self) = @_;
  return 1 if ($self->prev_exon && !$self->next_exon) ; 
  return 0 ;   
}

=head

   Name : is_5prim_exon 
   Arg  : none 
   Func : returns 1 if exon is at 5prim end of transcript 
   Return : int (0 || 1) 

=cut 

sub is_5prim_exon {
  my ($self) = @_;
  return 1 if ($self->next_exon && !$self->prev_exon) ; 
  return 0 ;   
}

=head

   Name : cluster
   Arg  : Bio::EnsEMBL::Transcript
   Func : getter/setter for the Bio::EnsEMBL::Transcript the Exon belongs to 
   Return : Bio::EnsEMBL::Transcript

=cut 

sub cluster {
  my ($self,$value) = @_;
  if (defined($value)) {
    $self->{'cluster'} = $value;
  }
  return $self->{'cluster'};
}


=head visited ( $val ) 

  Function : marks if exon has been visited in recursion procedure or not 
  Arg : integer (true / false ) 
  Returnval : integer (true || false ) 
  Caller : recursion procedure in Condense_EST.pm

=cut
 
sub visited  {
  my ($self,$value) = @_;
  if (defined($value)) {
    $self->{'visited'} = $value;
  }
  return $self->{'visited'};
}


  


