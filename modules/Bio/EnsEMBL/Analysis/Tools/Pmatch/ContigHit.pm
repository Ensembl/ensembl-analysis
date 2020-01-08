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



# holds a pmatch contig hit - simply the name(identifier) of the contig
# and a list of start-end positions

package Bio::EnsEMBL::Analysis::Tools::Pmatch::ContigHit;
use warnings ;
use strict ;
use vars qw(@ISA);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
@ISA = qw();


=head2 new

 Title   : new
 Usage   :
 Function: constructor
 Example :
 Returns : 
 Args    : 


=cut

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($id) = rearrange(['ID'], @args);

  throw("No id") unless defined $id;
  $self->id($id);
 
  $self->{'_forward_pairs'} = [];
  $self->{'_reverse_pairs'} = [];

  return $self;

}

=head2 id

 Title   : id
 Usage   :
 Function: get/set for contig id
 Example :
 Returns : 
 Args    : 


=cut

sub id {
  my ($self,$id) = @_;
  if ($id) {
    $self->{'id'} = $id;
  }
  return $self->{'id'};
}

=head2 add_CoordPair

 Title   : add_CoordPair
 Usage   :
 Function: adds a CoordPair to the list making up this hit
 Example :
 Returns : 
 Args    : 


=cut

sub add_CoordPair {
  my ($self,$pair) = @_;
  throw('No coord pair') unless defined $pair;
  throw('$pair is not a Bio::EnsEMBL::Analysis::Tools::Pmatch::CoordPair') unless $pair->isa("Bio::EnsEMBL::Analysis::Tools::Pmatch::CoordPair");
  if($pair->strand == 1) {
    push(@{$self->{_forward_pairs}},$pair);
  }
  else {
    push(@{$self->{_reverse_pairs}},$pair);
  }
}

=head2 each_ForwardPair

 Title   : each_ForwardPair
 Usage   :
 Function: returns CoordPairs represeting hits between a prtein and the forward strand of the contig
 Example :
 Returns : 
 Args    : 


=cut

sub each_ForwardPair {
  my ($self) = @_;
  return $self->{_forward_pairs};
}

=head2 each_ReversePair

 Title   : each_Reverseair
 Usage   :
 Function: returns CoordPairs representing hits between a protein and the reverse strand of the contig
 Example :
 Returns : 
 Args    : 


=cut

sub each_ReversePair {
  my ($self) = @_;
  return $self->{_reverse_pairs};
}
1;
