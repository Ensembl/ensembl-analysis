
=head1 NAME


=head1 SYNOPSIS

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Analysis::Tools::WGA2Genes::CoordUtils;

use strict;
use Exporter;

use vars qw (@ISA @EXPORT);

use Bio::EnsEMBL::Utils::Exception qw(verbose 
                                      throw 
                                      warning 
                                      stack_trace_dump);



@ISA = qw(Exporter);

@EXPORT = qw(check_consistent_coords
             distance_between_coords
             separate_coords
             merge_coords
            );



sub check_consistent_coords {
  my ($cj, $ci) = @_; 

  if ($ci->id ne $cj->id or
      $ci->strand != $cj->strand) {
    return 0;
  } 
  else {
    if ($ci->strand == 1) {
      # j should come before i in genomic coords
      if (not $cj->end < $ci->start) {
        return 0;
      } 
    } else {
      # j should come after i
      if (not $ci->end < $cj->start) {
        return 0;
      }
    }
  }
  
  return 1;
}


sub distance_between_coords {
  my ($cj, $ci) = @_; 

  if ($ci->strand == 1) {
    # j should come before i in genomic coords
    return $ci->start - $cj->end - 1;
  } else {
    return $cj->start - $ci->end - 1;
  }
}


sub merge_coords {
  my ($ci, $cj) = @_;

  # assumption: coords are consistent
  my ($start, $end);
  if ($cj->start > $ci->start) {
    ($start, $end) = ($ci->start, $cj->end); 
  } else {
    ($start, $end) = ($cj->start, $ci->end); 
  }

  return Bio::EnsEMBL::Mapper::Coordinate->new($ci->id,
                                               $start,
                                               $end,
                                               $ci->strand);
}



sub separate_coords {
  my ($ci, $cj, $target_slice) = @_;

  my ($ci_slice, $cj_slice, $between_slice);

  my $sa = $target_slice->adaptor; 

  if (defined $ci and defined $cj) {
    $ci_slice = $sa->fetch_by_region('toplevel',
                                     $target_slice->seq_region_name,
                                     $ci->start,
                                     $ci->end);
    
    $cj_slice = $sa->fetch_by_region('toplevel',
                                     $target_slice->seq_region_name,
                                     $cj->start,
                                     $cj->end);

    $between_slice = $sa->fetch_by_region('toplevel',
                                          $target_slice->seq_region_name,
                                          $ci->end + 1,
                                          $cj->start - 1);
  } elsif (defined $ci) {
    $ci_slice = $sa->fetch_by_region('toplevel',
                                     $target_slice->seq_region_name,
                                     $ci->start,
                                     $ci->end);

    $between_slice = $sa->fetch_by_region('toplevel',
                                          $target_slice->seq_region_name,
                                          $ci->end + 1,
                                          $target_slice->length);    

  } elsif (defined $cj) {
    $cj_slice = $sa->fetch_by_region('toplevel',
                                     $target_slice->seq_region_name,
                                     $cj->start,
                                     $cj->end);
    $between_slice = $sa->fetch_by_region('toplevel',
                                          $target_slice->seq_region_name,
                                          1,
                                          $cj->start - 1);

  }

  # following puts the Projection segments back into reference coords. 
  my (@seq_level_i, @seq_level_j, @between);

  if (defined $ci_slice) {
    @seq_level_i = map {
      {
        from_start => $_->from_start + $ci_slice->start - 1,
        from_end   => $_->from_end   + $ci_slice->start - 1,
        to_Slice   => $_->to_Slice,
      };
    } @{$ci_slice->project('seqlevel')};
  }
  if (defined $cj_slice) {
    @seq_level_j = map {
      {
        from_start => $_->from_start + $cj_slice->start - 1,
        from_end   => $_->from_end   + $cj_slice->start - 1,
        to_Slice   => $_->to_Slice,
      };
    } @{$cj_slice->project('seqlevel')};
  }
  @between = map {
    {
      from_start => $_->from_start + $between_slice->start - 1,
      from_end   => $_->from_end   + $between_slice->start - 1,
      to_Slice   => $_->to_Slice,
    };
  } @{$between_slice->project('seqlevel')};

  # if the whole between slice is accounted for by seq-level bits, 
  # then the ci and ck cannot be split either side of a gap
  my $between_seq_coverage = 0;
  foreach my $bit (@between) {
    $between_seq_coverage += $bit->{from_end} - $bit->{from_start} +1;
  }

  if (# no i coord given, so extend j
      not @seq_level_i or
      # no j coord given, so extend i      
      not @seq_level_j or
      # ci ends in a gap, so the two must be separable
      $seq_level_i[-1]->{from_end} < $ci->end or
      # cj starts with a gap so the two must be separable
      $seq_level_j[0]->{from_start} > $cj->start or
      # between region contained a sequence gap
      $between_seq_coverage < $between_slice->length) {
    
    # return a pair that extends to the first gap on each side
    my ($new_left, $new_right);

    if (defined $ci) {
      if (not @between or 
          $between[0]->{from_start} > $between_slice->start) {
        $new_left = Bio::EnsEMBL::Mapper::Coordinate->
            new($ci->id,
                $ci->start,
                $ci->end,
                $ci->strand);
      } else {
        $new_left = Bio::EnsEMBL::Mapper::Coordinate->
            new($ci->id,
                $ci->start,
                $between[0]->{from_end},
                $ci->strand);
      }
    }

    if (defined $cj) {
      if (not @between or $between[-1]->{from_end} < $between_slice->end) {
        $new_right = Bio::EnsEMBL::Mapper::Coordinate->
            new($cj->id,
                $cj->start,
                $cj->end,
                $cj->strand);
      } else {
        $new_right = Bio::EnsEMBL::Mapper::Coordinate->
            new($cj->id,
                $between[-1]->{from_start},
                $cj->end,
                $cj->strand);
      }
    }

    if (defined $new_left and defined $new_right) {
      return ($new_left, $new_right);
    } elsif (defined $new_left and not defined $new_right) {
      return ($new_left, undef);
    } elsif (defined $new_right and not defined $new_left) {
      return (undef, $new_right);
    }
  }

  return (undef, undef);
}



1;
