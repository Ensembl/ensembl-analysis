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

package Bio::EnsEMBL::Analysis::Tools::Otter::DnaAlignFeature;

# EnsEMBL module for storing dna-dna pairwise alignments
#
# You may distribute this module under the same terms as perl itself
#

=head1 NAME

  Bio::EnsEMBL::Analysis::Tools::Otter::DnaAlignFeature - Ensembl specific dna-dna pairwise alignment feature

=head1 SYNOPSIS

  This module inherits from the standard Bio::EnsEMBL::DnaDnaAlignFeature
  module. This module is required when reading DNA align features out
  of the Otter databases when one would like to attach the 
  dna_align_feature_history object to a dna_align_feature.

  In this module:
  The 'new' method from Bio::EnsEMBL::DnaDnaAlignFeature has been called
  to fetch all standard dna_align_feature attributes, and then the
  dna_align_feature_history object is also attached.

=cut


use warnings ;
use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

@ISA = qw( Bio::EnsEMBL::DnaDnaAlignFeature );


=head2 new

  Arg [..]   : List of named arguments. (incl. -dna_align_feature_history) defined
               in this constructor, others defined in BaseFeaturePair and 
               SeqFeature superclasses.  
  Example    : $daf = new DnaDnaAlignFeature(-cigar_string => '3M3I12M');
  Description: Creates a new DnaDnaAlignFeature using either a cigarstring or
               a list of ungapped features.  
  Returntype : Bio::EnsEMBL::DnaDnaAlignFeature
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);

  my ($dna_align_feature_history) = rearrange([qw(DNA_ALIGN_FEATURE_HISTORY)], @_);

  # get/set the dna_align_feature_history
  if (!defined $dna_align_feature_history){
    throw("dna_align_feature_history not defined");
  }
  $self->dna_align_feature_history($dna_align_feature_history);

  return $self;
}

sub dna_align_feature_history{
  my $self = shift;

  if (@_) {
    my $dafh = shift;
    if (defined $dafh && (!ref $dafh || ! $dafh->isa("Bio::EnsEMBL::Analysis::Tools::Otter::DnaAlignFeatureHistory"))) {
      throw("dna_align_feature_history argument must be a Bio::EnsEMBL::Analysis::Tools::Otter::DnaAlignFeatureHistory");
    }
    $self->{'dna_align_feature_history'} = $dafh;
  }

  return $self->{'dna_align_feature_history'};
}


1;
