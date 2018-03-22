# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

=head1 NAME


=head1 SYNOPSIS

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Analysis::Tools::WGA2Genes::CoordUtils;

use warnings ;
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
             merge_coords
             extend_coord
             transcripts_to_coords
            );


######################################################

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

######################################################

sub distance_between_coords {
  my ($cj, $ci) = @_; 

  if ($ci->strand == 1) {
    # j should come before i in genomic coords
    return $ci->start - $cj->end - 1;
  } else {
    return $cj->start - $ci->end - 1;
  }
}

######################################################

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

######################################################

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

############################################################

sub extend_coord {
  my ($coord, $target_slice) = @_;

  # this method pushes $coord out to the left and right 
  # to the point of the next sequence-level gap
  
  throw("Cannot extend coord; wrong slice given")
      if $coord->id ne $target_slice->seq_region_name;

  my @proj = @{$target_slice->project('seqlevel')};

  my ($new_start, $new_end, $up_seg, $down_seg);
  foreach my $seg (@proj) {
    if ($seg->from_start <= $coord->start and
        $seg->from_end   >= $coord->start) {
      $new_start = $seg->from_start;
    } elsif ($seg->from_start <= $coord->start) {
      $up_seg = $seg;
    }
    if ($seg->from_start <= $coord->end and
        $seg->from_end   >= $coord->end) {
      $new_end = $seg->from_end;
    } elsif ($seg->from_end >= $coord->end) {
      $down_seg = $seg if not defined $down_seg;
    }
  }

  throw("Could not find sequence-level pieces for " .
        $coord->id . "/" .
        $coord->start . "-" . 
        $coord->end)
      if not defined $new_start or not defined $new_end;


  my $ex_coord = Bio::EnsEMBL::Mapper::Coordinate->new($coord->id,
                                                       $new_start,
                                                       $new_end,
                                                       $coord->strand);

  my ($up_coord, $down_coord);
  if (defined $up_seg) {
    $up_coord = Bio::EnsEMBL::Mapper::Coordinate
        ->new($coord->id,
              $up_seg->from_start,
              $up_seg->from_end,
              $coord->strand);
  }
  if (defined $down_seg) {
    $down_coord = Bio::EnsEMBL::Mapper::Coordinate
        ->new($coord->id,
              $down_seg->from_start,
              $down_seg->from_end,
              $coord->strand);
  }


  return ($ex_coord, $up_coord, $down_coord);
}

######################################################

sub transcripts_to_coords {
  my ($map, $from_cs_name, @trans) = @_;

  my (@all_good_bits, @results);

  foreach my $tr (@trans) {
    my @exons = sort {$a->start <=> $b->start} @{$tr->get_all_translateable_Exons};

    my @bits;

    for(my $i = 0; $i < @exons; $i++) {
      my $e = $exons[$i];

      my @c = $map->map_coordinates($e->slice->seq_region_name,
                                    $e->start,
                                    $e->end,
                                    1,
                                    $from_cs_name);
      my $current_pos = $e->start;
      foreach my $c (@c) {
        my $current_end = $current_pos + $c->length - 1;

        if ($c->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
          push @bits, Bio::EnsEMBL::Mapper::Coordinate->new($e->slice->seq_region_name,
                                                            $current_pos,
                                                            $current_end,
                                                            1);
        } else {
          push @bits, $c;
        }
            
        $current_pos += $c->length;
      }
    }

    # policy: only consider gaps as worth filling if they are flanked
    # by aligned sequence that contains at least one full codon

    my $total_bps = 0;
    my $first_good = -1;

    for(my $i=0; $i < @bits; $i++) {
      if ($bits[$i]->isa("Bio::EnsEMBL::Mapper::Coordinate")) {

        my $spare_start = (3 - ($total_bps % 3)) % 3;
        my $spare_end = ($bits[$i]->length - $spare_start) % 3;
        
        my $codons = ($bits[$i]->length - $spare_start - $spare_end) / 3;
        
        if ($codons >= 1) {
          $first_good = $i;
          last;
        }
      }

      $total_bps += $bits[$i]->length;
    }

    next if $first_good < 0;

    @bits = splice(@bits, $first_good, scalar(@bits) - $first_good);

    $total_bps = 0;
    $first_good = -1;

    for(my $i=scalar(@bits) - 1; $i >= 0; $i--) {
      if ($bits[$i]->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
        my $spare_end = (3 - ($total_bps % 3)) % 3;
        my $spare_start = ($bits[$i]->length - $spare_end) % 3;
        
        my $codons = ($bits[$i]->length - $spare_start - $spare_end) / 3;
        
        if ($codons >= 1) {
          $first_good = $i;
          last;
        }
      }

      $total_bps += $bits[$i]->length;
    }

    next if $first_good < 0;

    @bits = splice(@bits, 0, $first_good + 1);

    push @all_good_bits, @bits;
  }

  @all_good_bits = sort { $a->start <=> $b->start } @all_good_bits;

  foreach my $g (@all_good_bits) {
    if (not @results or
        $results[-1]->end < $g->start - 1) {
      push @results, Bio::EnsEMBL::Mapper::Coordinate->
          new($trans[0]->slice->seq_region_name,
              $g->start,
              $g->end,
              1);
    } else {
      if ($g->end > $results[-1]->end) {
        $results[-1]->end($g->end);
      }
    }
  }

  return \@results;
}



1;
