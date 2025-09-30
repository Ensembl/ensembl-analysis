# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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
# holds a Coordinate Pair; name of the protein and contig that are hit (as a cross check) and 
# the start and end positions of the 

package Bio::EnsEMBL::Analysis::Tools::Pmatch::CoordPair;
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

  my ($query, $target, $qstart, $qend, $tstart, $tend, $percent, $strand) = 
    rearrange(['QUERY',
               'TARGET',
               'QSTART',
               'QEND',
               'TSTART',
               'TEND', 
               'PERCENT',
               'STRAND'],@args);
  
  throw("No query") unless defined $query;
  $self->query($query);

  throw("No target") unless defined $target;
  $self->target($target);

  throw("No query start") unless defined $qstart;
  $self->qstart($qstart);

  throw("No query end") unless defined $qend;
  $self->qend($qend);

  throw("No target start") unless defined $tstart;
  $self->tstart($tstart);

  throw("No target end") unless defined $tend;
  $self->tend($tend);

  throw("No percent") unless defined $percent;
  $self->percent($percent);

  throw("No strand") unless defined $strand;
  $self->strand($strand);

  return $self;

}

=head2 query

 Title   : query
 Usage   :
 Function: get/set for query (contig name)
 Example :
 Returns : 
 Args    : 


=cut

sub query {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'query'} = $arg;
  }
  return $self->{'query'};
}

=head2 target

 Title   : target
 Usage   :
 Function: get/set for target (protein name)
 Example :
 Returns : 
 Args    : 


=cut

sub target {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'target'} = $arg;
  }
  return $self->{'target'};
}

sub qstart {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'qstart'} = $arg;
  }
  return $self->{'qstart'};
}

sub qend {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'qend'} = $arg;
  }
  return $self->{'qend'};
}

sub tstart {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'tstart'} = $arg;
  }
  return $self->{'tstart'};
}

sub tend {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'tend'} = $arg;
  }
  return $self->{'tend'};
}

sub percent {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'percent'} = $arg;
  }
  return $self->{'percent'};
}

sub strand {
  my ($self,$arg) = @_;
  if ($arg) {
    $self->{'strand'} = $arg;
  }
  return $self->{'strand'};
}



1;
