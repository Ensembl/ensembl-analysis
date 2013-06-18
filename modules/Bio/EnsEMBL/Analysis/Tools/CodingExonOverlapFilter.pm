# Ensembl module for Bio::EnsEMBL::Analysis::Tools::CodingExonOverlapFilter
#
# Copyright (c) 2004 Ensembl
#


# $Source: /tmp/ENSCOPY-ENSEMBL-ANALYSIS/modules/Bio/EnsEMBL/Analysis/Tools/CodingExonOverlapFilter.pm,v $
# $Revision: 1.1 $
package Bio::EnsEMBL::Analysis::Tools::CodingExonOverlapFilter;

use strict;
use warnings;

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );



sub new{
  my ($class, @args) = @_;
  my $self = bless {},$class;

  if (scalar(@args)) {
    throw("CodingExonOverlapFilter should have no args in new");
  }

  return $self;
}

#####################################
sub filter {
  my ($self, $these, $others) = @_;

  # interference is judged by overlap at exon level
  # assumption is that @others is sorted by gene start

  my @filtered;

  my $cur_idx = 0;

  foreach my $obj (@$these) {
    my (@genomic_overlap, $left_bound);


    for(my $i=$cur_idx; $i < @$others; $i++) {
      my $o_obj = $others->[$i];

      if ($o_obj->end >= $obj->start and not defined $left_bound) {
        $left_bound = $i;
      }

      if ($o_obj->end < $obj->start) {
        next;
      } elsif ($o_obj->start > $obj->end) {
        last;
      } else {
        push @genomic_overlap, $o_obj;
      }
    }

    $cur_idx = $left_bound if defined $left_bound;

    my $exon_overlap = 0;
    if (@genomic_overlap) {
      my @exons = @{$obj->get_all_translateable_Exons};
      OG: foreach my $o_obj (@genomic_overlap) {
        foreach my $oe (@{$o_obj->get_all_translateable_Exons}) {
          foreach my $e (@exons) {
            if ($oe->strand == $e->strand and
                $oe->end >= $e->start and
                $oe->start <= $e->end) {
              $exon_overlap = 1;
              last OG;
            }
          }
        }
      }
    }

    if (not $exon_overlap) {
      push @filtered, $obj;
    }
  }

  return \@filtered;
}
1;
