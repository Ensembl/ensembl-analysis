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
# holds a pmatch protein hit - simply the name(idenitifier) of the protein
# and a list of ContigHit objects

package Bio::EnsEMBL::Analysis::Tools::Pmatch::ProteinHit;

use warnings ;
use strict ;

use Bio::EnsEMBL::Utils::Argument qw (rearrange);

=head2 new

 Title   : new
 Usage   :
 Function:
 Example :
 Returns : 
 Args    : constructor


=cut

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  my ($id) = rearrange(['ID'], @args);

  $self->throw("No id") unless defined $id;
  $self->id($id);
  
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

=head2 add_ContigHit

 Title   : add_ContigHit
 Usage   :
 Function: adds a ContigHit into $self->{_contig_hits}
 Example :
 Returns : 
 Args    :


=cut

sub add_ContigHit {
  my ($self,$hit) = @_;
  $self->throw('No contig hit') unless defined $hit;
  $self->throw('$contig is not a Bio::EnsEMBL::Analysis::Tools::Pmatch::ContigHit') unless $hit->isa("Bio::EnsEMBL::Analysis::Tools::Pmatch::ContigHit");
  $self->{_contig_hits}{$hit->id()} = $hit;
}

=head2 each_ContigHit

 Title   : each_ContigHit
 Usage   :
 Function: returns all entries in $self->{_contig_hits}
 Example :
 Returns : 
 Args    :


=cut

sub each_ContigHit {
  my ($self) = @_;
  return values %{$self->{_contig_hits}};
}


=head2 get_ContigHit

 Title   : get_ContigHit
 Usage   :
 Function: returns entries for a particular contig in $self->{_contig_hits}
 Example :
 Returns : 
 Args    :


=cut

sub get_ContigHit {
  my ($self,$contig) = @_;
  return ($self->{_contig_hits}{$contig}) if defined $contig;
  return undef;
}

1;
