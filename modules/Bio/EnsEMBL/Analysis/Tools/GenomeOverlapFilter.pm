# Ensembl module for Bio::EnsEMBL::Analysis::Tools::GenomeOverlapFilter
#
# Copyright (c) 2004 Ensembl
#


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

