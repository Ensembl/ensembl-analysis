
=head1 NAME

Bio::EnsEMBL::Analysis::Tools::WGA2Genes::GeneScaffold 

=head1 SYNOPSIS

A Bio::EnsEMBL::Slice that is comprised
of pieces of different target sequences, inferred by an alignment
of the target to a finished, query genome sequence. This object
extends Slice with mappings to/from the query and target

=cut


package Bio::EnsEMBL::Analysis::Tools::WGA2Genes::GeneScaffold;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Utils::Argument qw(rearrange);

use Bio::EnsEMBL::Slice;

@ISA = qw(Bio::EnsEMBL::Slice);


sub new {
  my $caller = shift;
  my %given_args = @_;

  my $class = ref($caller) || $caller;

  $given_args{-coord_system} =  
      Bio::EnsEMBL::CoordSystem->new(-name => 'genescaffold',
                                     -rank => 1);

  my $self = $class->SUPER::new(%given_args);

  my ($qmapper,$tmapper) = rearrange([qw(QUERY_MAPPER
                                         TARGET_MAPPER)], @_);
  if (defined $qmapper) {
    $self->query_mapper($qmapper);
  }
  if (defined $tmapper) {
    $self->target_mapper($tmapper);
  }
  
  return $self;
}



###################################################################
# FUNCTION   : project_transcript
#
# Description:
#    Takes a transcript, and uses the mapping between 
#    query coords and gene scaffold coords to produces a transcript 
#    that is the result of "projecting" the original transcript, 
#    through alignment, onto the gene scaffold. 
###################################################################

sub project_transcript {
  my ($self, 
      $tran, 
      $max_readthrough) = @_;

  my ($tran_length, @all_coords, @new_exons);

  my @orig_exons = @{$tran->get_all_translateable_Exons};
  if ($tran->strand < 0) {
    @orig_exons = reverse @orig_exons;
  }
  map { $tran_length += $_->length } @orig_exons; 

  foreach my $orig_exon (@orig_exons) {    
    my @crds = $self->query_mapper->map_coordinates($orig_exon->slice->seq_region_name,
                                                    $orig_exon->start,
                                                    $orig_exon->end,
                                                    1,
                                                    "query");
    push @all_coords, @crds;
  }


  my $start_not_found = 0;
  my $end_not_found   = 0;

  # Replace coords at start and end that map down to gaps with gaps.
  # Although we have already trimmed back the gene scaffold for gap 
  # exons at the ends, individual transripts could begin/end anywhere
  # in the gene scaffold. 

  for(my $i=0; $i < @all_coords; $i++) {
    if ($all_coords[$i]->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
      my ($tcoord) = $self->target_mapper->map_coordinates($self->seq_region_name,
                                                           $all_coords[$i]->start,
                                                           $all_coords[$i]->end,
                                                           1,
                                                           $self->seq_region_name);
      if ($tcoord->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
        last;
      } else {
        $all_coords[$i] = 
            Bio::EnsEMBL::Mapper::Gap->new(1, $tcoord->length);
        $start_not_found = 1;
      }
    } else {
      $start_not_found = 1;
    }
  }
  for(my $i=scalar(@all_coords)-1; $i >= 0; $i--) {
    if ($all_coords[$i]->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
      my ($tcoord) = $self->target_mapper->map_coordinates($self->seq_region_name,
                                                           $all_coords[$i]->start,
                                                           $all_coords[$i]->end,
                                                           1,
                                                           $self->seq_region_name);
      if ($tcoord->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
        last;
      } else {
        $all_coords[$i] = 
            Bio::EnsEMBL::Mapper::Gap->new(1, $tcoord->length);
        $end_not_found = 1;
      }
    } else {
      $end_not_found = 1;
    }
  }

  my $need_another_pass;
  do {
    $need_another_pass = 0;

    my (@proc_coords, @gap_indices);
    # merge gaps
    foreach my $c (@all_coords) {
      if ($c->isa("Bio::EnsEMBL::Mapper::Gap")) {
        if (@proc_coords and 
            $proc_coords[-1]->isa("Bio::EnsEMBL::Mapper::Gap")) {
          $proc_coords[-1]->end( $proc_coords[-1]->end + $c->length );
        } else {
          push @proc_coords, Bio::EnsEMBL::Mapper::Gap->new(1, 
                                                            $c->length);
          push @gap_indices, scalar(@proc_coords) - 1;
        }
      } else {
        push @proc_coords, $c;
      }
    }

    GAP: foreach my $idx (@gap_indices) {
      my $gap = $proc_coords[$idx];
      my $frameshift = $gap->length % 3;

      if ($frameshift) {
        my $bases_to_remove = 3 - $frameshift;      

        # calculate "surplus" bases on incomplete codons to left and right
        my ($left_surplus, $right_surplus) = (0,0);
        for(my $j=$idx-1; $j >= 0; $j--) {
          $left_surplus += $proc_coords[$j]->length;
        }
        for(my $j=$idx+1; $j < @proc_coords; $j++) {
          $right_surplus += $proc_coords[$j]->length;
        }
        
        $left_surplus  = $left_surplus % 3;
        $right_surplus = $right_surplus % 3;

        if ($left_surplus) {
          # eat left
          $bases_to_remove = $left_surplus;
          
          my $left_coord = $proc_coords[$idx - 1];
          if ($left_coord->length > $bases_to_remove) {
            $gap->end($gap->end + $bases_to_remove);
            $left_coord->end( $left_coord->end - $bases_to_remove );
          } else {
            # we need to eat away the whole of this coord
            $proc_coords[$idx-1] = 
                Bio::EnsEMBL::Mapper::Gap->new(1,$left_coord->length);
          }
        }
        if ($right_surplus) {
          $bases_to_remove = $right_surplus;

          my $right_coord = $proc_coords[$idx + 1];
          if ($right_coord->length > $bases_to_remove) {
            $gap->end($gap->end + $bases_to_remove);
            $right_coord->start( $right_coord->start + $bases_to_remove);
          } else {
            # we need to eat away the whole of this coord
            $proc_coords[$idx+1] = 
                Bio::EnsEMBL::Mapper::Gap->new(1,$right_coord->length);
          }
        }
        
        $need_another_pass = 1;
        last GAP;
      }      
    }
    @all_coords = @proc_coords;    
  } while ($need_another_pass);

  my ($total_tran_bps, $real_seq_bps);
  foreach my $coord (@all_coords) {
    if ($coord->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
      push @new_exons, Bio::EnsEMBL::Exon->new(-start => $coord->start,
                                               -end   => $coord->end,
                                               -strand => $tran->strand,
                                               -slice => $self);

      $total_tran_bps += $coord->length;
      my ($tcoord) = $self->target_mapper->map_coordinates($self->seq_region_name,
                                                           $coord->start,
                                                           $coord->end,
                                                           1,
                                                           $self->seq_region_name);
      if ($tcoord->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
        $real_seq_bps += $coord->length;
      } 
    }
  }

  if (not @new_exons) {
    # the whole transcript mapped to gaps
    return 0;
  }

  #
  # sort exons into rank order 
  #
  if ($tran->strand < 0) {
    @new_exons = sort { $b->start <=> $a->start } @new_exons;
  } else {
    @new_exons = sort { $a->start <=> $b->start } @new_exons;
  }

  #
  # calculate phases, and add supporting features
  #
  my ($previous_exon);
  foreach my $exon (@new_exons) {

    if (defined $previous_exon) {
      $exon->phase($previous_exon->end_phase);
    } else {
      $exon->phase(0);
    }

    $exon->end_phase((($exon->end - $exon->start + 1) + $exon->phase)%3);

    # need to map back to the genomic coords to get the supporting feature
    # for this exon;
    my $extent_start = $exon->start;
    my $extent_end   = $exon->end;
    if ($exon->strand > 0) {
      $extent_start += 3 - $exon->phase if $exon->phase;
      $extent_end   -= $exon->end_phase if $exon->end_phase;
    } else {
      $extent_start += $exon->end_phase if $exon->end_phase;
      $extent_end   -=  3 - $exon->phase if $exon->phase;
    }

    if ($extent_end > $extent_start) {
      # if not, we've eaten away the whole exon, so there is no support

      my @gen_coords = $self->query_mapper->map_coordinates($self->seq_region_name,
                                                            $extent_start,
                                                            $extent_end,
                                                            1,
                                                            $self->seq_region_name);

      my @fps;
      my $cur_gs_start = $extent_start;
      foreach my $g_coord (@gen_coords) {
        my $cur_gs_end = $cur_gs_start + $g_coord->length - 1;
                
        if ($g_coord->isa("Bio::EnsEMBL::Mapper::Coordinate")) {

          my ($p_coord) = $tran->genomic2pep($g_coord->start, 
                                             $g_coord->end,
                                             $exon->strand);
          
          my $fp = Bio::EnsEMBL::FeaturePair->
              new(-seqname  => $self->seq_region_name,
                  -start    => $cur_gs_start,
                  -end      => $cur_gs_end,
                  -strand   => $exon->strand,
                  -score    => 100.0,
                  -hseqname => $tran->translation->stable_id,
                  -hstart   => $p_coord->start,
                  -hend     => $p_coord->end,
                  -hstrand => $p_coord->strand);
          push @fps, $fp;
        }
        
        $cur_gs_start += $g_coord->length;
      }
        
      if (@fps) {
        my $f = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@fps);
        $exon->add_supporting_features($f);
      }
    }

    $previous_exon = $exon;
  }

  #
  # merge abutting exons; deals with exon fusion events, and 
  # small, frame-preserving insertions in the target
  #
  my @merged_exons;
  foreach my $exon (@new_exons) {
    if (@merged_exons) {

      my $prev_exon = pop @merged_exons;
   
      my ($new_start, $new_end);

      if ($tran->strand < 0) {
        my $intron_len = $prev_exon->start - $exon->end - 1;
        if ($intron_len % 3 == 0 and 
            $intron_len <= $max_readthrough) { 
          $new_start = $exon->start;
          $new_end   = $prev_exon->end;
        }
      } else {
        my $intron_len = $exon->start - $prev_exon->end - 1;
        if ($intron_len % 3 == 0 and 
            $intron_len <= $max_readthrough) {
          $new_start = $prev_exon->start;
          $new_end   = $exon->end;
        }
      }

      if (defined $new_start and defined $new_end) {
        my $merged_exon = Bio::EnsEMBL::Exon->
            new(-start => $new_start,
                -end   => $new_end,
                -strand => $tran->strand,
                -phase => $prev_exon->phase,
                -end_phase => $exon->end_phase,
                -slice  => $exon->slice);
        
        my @ug_feats;
        if (@{$prev_exon->get_all_supporting_features}) {
          my ($sf) = @{$prev_exon->get_all_supporting_features};
          push @ug_feats, $sf->ungapped_features;
        }
        if (@{$exon->get_all_supporting_features}) {
          my ($sf) = @{$exon->get_all_supporting_features};
          push @ug_feats, $sf->ungapped_features;
        }
        if (@ug_feats) {
          my $new_sup_feat = Bio::EnsEMBL::DnaPepAlignFeature->
              new(-features => \@ug_feats);
          $merged_exon->add_supporting_features($new_sup_feat);
        }

        push @merged_exons, $merged_exon;
        next;
      } else {
        push @merged_exons, $prev_exon;
        push @merged_exons, $exon;
      }      
    } else {
      push @merged_exons, $exon;
    }
  }
  

  my $proj_tran = Bio::EnsEMBL::Transcript->new();

  my (@trans_fps);
  foreach my $exon (@merged_exons) {
    $proj_tran->add_Exon($exon);
    
    if (@{$exon->get_all_supporting_features}) {
      my ($sf) = @{$exon->get_all_supporting_features};
      my @e_fps = $sf->ungapped_features;
      push @trans_fps, @e_fps;
    }
  }

  #
  # do transcript-level supporting features/attributes
  #
  my $t_sf = Bio::EnsEMBL::DnaPepAlignFeature->
      new(-features => \@trans_fps);
  $t_sf->hcoverage( 100 * ($total_tran_bps / $tran_length) );
  $proj_tran->add_supporting_features($t_sf);

  #
  # set translation
  #
  my $translation = Bio::EnsEMBL::Translation->new();
  $translation->start_Exon($merged_exons[0]);
  $translation->start(1);
  $translation->end_Exon($merged_exons[-1]);
  $translation->end($merged_exons[-1]->end - $merged_exons[-1]->start + 1);

  $proj_tran->translation($translation);

  if (not defined $proj_tran->translate) {
    # this can happen if the transcript comprises a single stop codon only
    return 0;
  }

  #
  # finally, attributes
  #
  my @attributes;

  my $gap_attr = Bio::EnsEMBL::Attribute->
      new(-code => 'PropNonGap',
          -name => 'proportion non gap',
          -description => 'proportion non gap',
          -value => sprintf("%.1f", 
                            100 * ($real_seq_bps / $total_tran_bps)));
  push @attributes, $gap_attr;
  
  if ($start_not_found and $tran->strand > 0 or
      $end_not_found and $tran->strand < 0) {
    my $attr = Bio::EnsEMBL::Attribute->
        new(-code => 'StartNotFound',
            -name => 'start not found',
            -description => 'start not found',
            -value => 1);
    push @attributes, $attr;
  }
  if ($end_not_found and $tran->strand > 0 or
      $start_not_found and $tran->strand < 0) {
    my $attr = Bio::EnsEMBL::Attribute->
        new(-code => 'EndNotFound',
            -name => 'end not found',
            -description => 'end not found',
            -value => 1);
    push @attributes, $attr;
  }

  my $pep = $proj_tran->translate->seq;
  my $num_stops = $pep =~ tr/\*/\*/;

  my $stop_attr = Bio::EnsEMBL::Attribute->
      new(-code => 'NumStops',
          -name => 'number of stops',
          -desc => 'Number of stops before editing',
          -value => $num_stops);
  push @attributes, $stop_attr;

  # indentify gap exons
  my $gap_exons = 0;
  foreach my $e (@{$proj_tran->get_all_Exons}) {
    my ($coord) = $self->target_mapper->map_coordinates($self->seq_region_name,
                                                        $e->start,
                                                        $e->end,
                                                        1,
                                                        $self->seq_region_name);
    if ($coord->isa("Bio::EnsEMBL::Mapper::Gap")) {
      $gap_exons++;
    }
  }
  my $gap_exon_attr = Bio::EnsEMBL::Attribute->
      new(-code => 'GapExons',
          -name => 'gap exons',
          -description => 'number of gap exons',
          -value => $gap_exons);
  push @attributes, $gap_exon_attr;

  #if (defined $gene_id) {
  #  my $geneid_attr = Bio::EnsEMBL::Attribute->
  #      new(-code => 'SourceGene',
  #          -name => 'source gene',
  #          -description => 'human source gene',
  #          -value => $gene_id);
  #  push @attributes, $geneid_attr;
  #}
  my $tranid_attr = Bio::EnsEMBL::Attribute->
      new(-code => 'SourceTran',
          -name => 'source transcript',
          -description => 'source transcript',
          -value => $tran->stable_id);
  push @attributes, $tranid_attr;

  $proj_tran->add_Attributes(@attributes);

  return $proj_tran;
}


sub query_components {
  my ($self) = @_;

  return $self->query_mapper->map_coordinates($self->seq_region_name,
                                              1,
                                              $self->length,
                                              1,
                                              $self->seq_region_name);
}

sub target_components {
  my ($self) = @_;

  return $self->target_mapper->map_coordinates($self->seq_region_name,
                                               1,
                                               $self->length,
                                               1,
                                               $self->seq_region_name);
}


sub query_mapper {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_query_mapper} = $val;
  }

  return $self->{_query_mapper};
}


sub target_mapper {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_target_mapper} = $val;
  }
  return $self->{_target_mapper};

}



1;
