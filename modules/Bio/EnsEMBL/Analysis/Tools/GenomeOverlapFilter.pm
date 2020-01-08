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


package Bio::EnsEMBL::Analysis::Tools::GenomeOverlapFilter;

use strict;
use warnings;

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );



sub new{
  my ($class, @args) = @_;
  my $self = bless {},$class;

  if (scalar(@args)) {
    throw("GenomeOverlapFilter should have no args in new");
  }

  return $self;
}

#####################################
sub filter {
  my ($self, $these, $others) = @_;

  # interference is judged by overlap at genomic level
  # assumption is that @others is sorted by gene start

  my @filtered;

  foreach my $obj (@$these) {
    my ($left_bound, $genomic_overlap);

    for(my $i=0; $i < @$others && !$genomic_overlap; $i++) {
      my $o_obj = $others->[$i];

      next if $o_obj->strand != $obj->strand;

      if ($o_obj->end < $obj->start) {
        next;
      } elsif ($o_obj->start > $obj->end) {
        last;
      } else {
        $genomic_overlap = 1;
      }
    }

    if (not $genomic_overlap) {
      push @filtered, $obj;
    }
  }

  return \@filtered;
}

1;
