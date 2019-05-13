# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Tools::WGA2Genes::ChainUtils;

use warnings ;
use strict;
use Exporter;

use vars qw (@ISA @EXPORT);

use Bio::EnsEMBL::Utils::Exception qw(verbose 
                                      throw 
                                      warning); 


@ISA = qw(Exporter);

@EXPORT = qw(filter_irrelevant_chains
             filter_inconsistent_chains
             filter_block_overlap_chains             
             flatten_chains
             stringify_chains
            );



###################################################################
# FUNCTION: filter_irrelevant_chains
#
# Decription:
#  This function returns the chains that have a block overlap
#  with at least one feature in the given list (thus removing
#  "irrelevant" chains that do not interact with exons)
###################################################################
sub filter_irrelevant_chains {
  my ($chains, $feats) = @_;

  my @kept_chains;

  my @feats = sort { $a->start <=> $b->start } @$feats;

  foreach my $c (@$chains) {
    my $chain_overlaps_feat = 0;

    my $j=0;
    BLOCK: foreach my $b (@$c) {
      my $b_st = $b->reference_genomic_align->dnafrag_start;
      my $b_en = $b->reference_genomic_align->dnafrag_end;

      foreach my $this_f (@feats) {

        if ($this_f->end < $b_st) {
          next;
        } elsif ($this_f->start > $b_en) {
          last;
        } else {
          # block overlap
          $chain_overlaps_feat = 1;
          last BLOCK;
        }
      }
    }

    if ($chain_overlaps_feat) {
      push @kept_chains, $c;
    }
  }

  return \@kept_chains;
}



###################################################################
# FUNCTION: filter_block_overlap_chains
#
# Decription:
#
#  This function filters out chains from the first list that have
#  block overlap with any of the given blocks
###################################################################
sub filter_block_overlap_chains {
  my ($raw_chains, $given_blocks) = @_;

  my @chains_to_process;

  foreach my $chain (@$raw_chains) {
    # if this chain conatains a block that overlaps with any of the 
    # excluded regions on the target, ignore it.
    my $keep_chain = 1;

    CHAIN: foreach my $block (@$chain) {
      my ($tga) = @{$block->get_all_non_reference_genomic_aligns};
      foreach my $ex_block (@$given_blocks) {
        my ($ex_tga) = @{$ex_block->get_all_non_reference_genomic_aligns};

        if ($tga->dnafrag->name eq $ex_tga->dnafrag->name and
            $tga->dnafrag_end >= $ex_tga->dnafrag_start and
            $tga->dnafrag_start <= $ex_tga->dnafrag_end) {

          $keep_chain = 0;
          last CHAIN;
        }
      }
    }
    if ($keep_chain) {
      push @chains_to_process, $chain;
    }
  }

  return \@chains_to_process;
}


###################################################################
# FUNCTION: filter_inconsistent_chains
#
# Decription:
#
# Progressively projects chains on to the query sequence (in order
# of score), rejecting chains that "interfere" with the set of 
# chains retained so far. Interference is defined as an overlap
# in the query at the block level
###################################################################
sub filter_inconsistent_chains {
  my ($input_chains, $overlap_chain_thresh) = @_;

  my (@sorted_chains, @to_ignore_chains, @kept_chains);

  foreach my $bl_list (@$input_chains) {
    my $qga_l = $bl_list->[0]->reference_genomic_align;
    my $qga_r = $bl_list->[-1]->reference_genomic_align;
    my ($tga_l) = @{$bl_list->[0]->get_all_non_reference_genomic_aligns};
    my ($tga_r) = @{$bl_list->[-1]->get_all_non_reference_genomic_aligns};
    
    push @sorted_chains, {
      blocks => $bl_list,
      score  => $bl_list->[0]->score,
      qid    => $qga_l->dnafrag->name,
      qstart  => $qga_l->dnafrag_start,
      qend    => $qga_r->dnafrag_end,
      tid     => $tga_l->dnafrag->name,
      tstart  => $tga_l->dnafrag_strand > 0 ? $tga_l->dnafrag_start : $tga_r->dnafrag_start,
      tend    => $tga_l->dnafrag_strand > 0 ? $tga_r->dnafrag_end : $tga_l->dnafrag_end,
      tstrand => $tga_l->dnafrag_strand,
    };      

  }
  
  for(my $i=0; $i < @sorted_chains; $i++) {

    next if $to_ignore_chains[$i];

    my (@ov_chains, @block_ov_regs);

    my $c = $sorted_chains[$i];

    my @these_blocks = @{$c->{blocks}};

    foreach my $kc (@kept_chains) {
      if ($c->{qstart} <= $kc->{qend} and 
          $c->{qend} >= $kc->{qstart}) {

        push @ov_chains, $kc;
      }
    }

    for (my $i=0; $i < @these_blocks; $i++) {
      my $b = $these_blocks[$i];
      my $qga = $b->reference_genomic_align;

      $block_ov_regs[$i] = []; 

      KCHAIN: foreach my $kc (@ov_chains) {
        foreach my $kb (@{$kc->{blocks}}) {
          
          my $k_qga = $kb->reference_genomic_align;
          
          if ($qga->dnafrag_start <= $k_qga->dnafrag_end and
              $qga->dnafrag_end >= $k_qga->dnafrag_start) {

            my $ov_start = $k_qga->dnafrag_start;
            $ov_start = $qga->dnafrag_start 
                if $qga->dnafrag_start > $ov_start;
            my $ov_end = $k_qga->dnafrag_end;
            $ov_end = $qga->dnafrag_end 
                if $qga->dnafrag_end < $ov_end;

            push @{$block_ov_regs[$i]}, {
              start => $ov_start,
              end   => $ov_end,
            };
          }
        }
      }
    }

    my $rejected = 0;

    if (not grep { scalar(@{$_}) > 0 } @block_ov_regs) {
      # no overlapping blocks, so keep chain as is
      push @kept_chains, $c;
    } else {
      # if the amount of overlap only constitutes a small proportion
      # of the chain, keep the chain (with the overlapping blocks
      # removed)

      my $total_chain_len = 0;

      for(my $i=0; $i < @these_blocks; $i++) {
        my $b = $these_blocks[$i];
        my $ga = $b->reference_genomic_align;

        $total_chain_len += $ga->dnafrag_end - $ga->dnafrag_start + 1;
      }

      # flatten the overlapping regions
      for (my $i=0; $i < @block_ov_regs; $i++) {
        my @regs = @{$block_ov_regs[$i]};

        my @new_regs;
        foreach my $reg (sort { $a->{start} <=> $b->{start} } @regs) {
          if (not @new_regs or $reg->{start} > $new_regs[-1]->{end}) {
            push @new_regs, $reg;
          } else {
            if ($reg->{end} > $new_regs[-1]->{end}) {
              $new_regs[-1]->{end} = $reg->{end};
            }
          }
        }

        $block_ov_regs[$i] = \@new_regs;
      }

      # we're only going to allow simple truncation of the chain
      # to resolve the overlap, at one end, or both. If by
      # truncating the chain in this way we lose too much of
      # it, again we reject the whole thing. 
      
      # find the longest, contigous, non-overlapping stretch of chain
      my @contig_chain_parts;
      my $last_was_non_ov = 0;
      for(my $i=0; $i < @these_blocks; $i++) {
        
        my $block = $these_blocks[$i];
        my $qga = $block->reference_genomic_align;
        my $ovs = $block_ov_regs[$i];
        
        if (@$ovs) {
          if ($ovs->[0]->{start} > $qga->dnafrag_start) {
            if ($last_was_non_ov) {
              push @{$contig_chain_parts[-1]}, {
                index => $i,
                start => $qga->dnafrag_start,
                end   => $ovs->[0]->{start} - 1,
              };
            } else {
              push @contig_chain_parts, [{
                index => $i,
                start => $qga->dnafrag_start,
                end   => $ovs->[0]->{start} - 1,
              }];
            }
          }
          
          for(my $j=1; $j < @$ovs; $j++) {
            push @contig_chain_parts, [{
              index => $i,
              start => $ovs->[$j-1]->{end} + 1,
              end   => $ovs->[$j]->{start} - 1,
            }];
          }
          if ($ovs->[-1]->{end} < $qga->dnafrag_end) {
            push @contig_chain_parts, [{
              index => $i,
              start => $ovs->[-1]->{end} + 1,
              end   => $qga->dnafrag_end,
            }];
            $last_was_non_ov = 1;
          } else {
            $last_was_non_ov = 0;
          }
        } else {
          if ($last_was_non_ov) {
            push @{$contig_chain_parts[-1]}, {
              index => $i,
              start => $qga->dnafrag_start,
              end   => $qga->dnafrag_end,
            };
          } else {
            push @contig_chain_parts, [{
              index => $i,
              start => $qga->dnafrag_start,
              end   => $qga->dnafrag_end,
            }];
          }
          $last_was_non_ov = 1;
        }
      }
      
      # find longest 
      my ($longest, $longest_len);

      if (not @contig_chain_parts) {
        $longest_len = 0;
      } else {
        foreach my $reg_list (@contig_chain_parts) {
          my $len = 0;
          foreach my $reg (@$reg_list) {
            $len += $reg->{end} - $reg->{start} + 1;
          }
          if (not defined $longest or $len > $longest_len) {
            $longest = $reg_list;
            $longest_len = $len;
          }
        }
      }

      # we've identifed a sub-region of the chain to keep.
      # But is it worth keeping?
      if ($longest_len / $total_chain_len > $overlap_chain_thresh) {
        my @bs;
        foreach my $el (@$longest) {
          my $block = $these_blocks[$el->{index}];
          my $st = $el->{start};
          my $en = $el->{end};
          my $qga = $block->reference_genomic_align;

          if ($st != $qga->dnafrag_start or
              $en != $qga->dnafrag_end) {
            my $new_block = $block->restrict_between_reference_positions($st, $en);
            $new_block->score($block->score);
            push @bs, $new_block;
          } else {
            push @bs, $block;
          }
        };

        my $qga_l = $bs[0]->reference_genomic_align;
        my $qga_r = $bs[-1]->reference_genomic_align;
        my ($tga_l) = @{$bs[0]->get_all_non_reference_genomic_aligns};
        my ($tga_r) = @{$bs[-1]->get_all_non_reference_genomic_aligns};

        $c->{blocks} = \@bs;
        $c->{qstart} = $qga_l->dnafrag_start;
        $c->{qend} = $qga_r->dnafrag_end;
        $c->{tstart} = ($tga_l->dnafrag_strand > 0) 
            ? $tga_l->dnafrag_start 
            : $tga_r->dnafrag_start;
        $c->{tend} = ($tga_l->dnafrag_strand > 0) 
            ? $tga_r->dnafrag_end 
            : $tga_l->dnafrag_end;

        push @kept_chains, $c;
      } else {
        $rejected = 1;
      }
    }
    if ($rejected) {
      # also reject all lower scoring chains of the same target
      for(my $j=$i+1; $j < @sorted_chains; $j++) {
        if ($sorted_chains[$j]->{tid} eq $c->{tid}) {
          $to_ignore_chains[$j] = 1;
        }
      }
    } else {
      # also keep all lower scoring chains of the same target that 
      # are (a) consistent with this one in terms of coords and
      # (b) do not interfere with any kept chains so far
      for(my $j=$i+1; $j < @sorted_chains; $j++) {
        my $oc = $sorted_chains[$j];

        my $consistent = 0;
        if ($oc->{tid} eq $c->{tid} and
            $oc->{tstrand} eq $c->{tstrand}) {
          if ($oc->{qstart} > $c->{qend}) {
            if ($oc->{tstrand} > 0 and $oc->{tstart} > $c->{tend} or
                $oc->{tstrand} < 0 and $oc->{tend} < $c->{tstart}) {
              $consistent = 1;
            }
          } elsif ($oc->{qend} < $c->{qstart}) {
            if ($oc->{tstrand} > 0 and $c->{tstart} > $oc->{tend} or
                $oc->{tstrand} < 0 and $c->{tend} < $oc->{tstart}) {
              $consistent = 1;
            }
          }
        }
        if ($consistent) {
          # any query overlap with other retained chains?
          for(my $k=0; $k < @kept_chains; $k++) {
            if ($oc->{qstart} <= $kept_chains[$k]->{qend} and
                $oc->{qend} >= $kept_chains[$k]->{qstart}) {
              $consistent = 0;
              last;
            }
          }
          if ($consistent) {
            push @kept_chains, $oc;
            $to_ignore_chains[$j] = 1;
          }
        }          
      }
    }
  }
  
  return [map { $_->{blocks} } @kept_chains];
}


###################################################################
# FUNCTION: flatten_chains
#
# Description:
#    Takes a list of chains of GenomicAlignBlocks and computes 
#    and flattens it into a list of GenomicAlignBlocks
#    NOTE: The method assumes that the input chains are non-
#    overlapping in the reference sequence at the block level. If
#    the '$check' flag is supplied, this is checked and an exception
#    thrown if this is not the case
####################################################################

sub flatten_chains {
  my ($chains, $check) = @_;

  my @blocks;
  foreach my $chain (@$chains) {
    push @blocks, @$chain;
  }
  
  @blocks = sort {
    $a->reference_genomic_align->dnafrag_start <=>
        $b->reference_genomic_align->dnafrag_start;
  } @blocks;

  if ($check) {
    for(my $i=1; $i < @blocks; $i++) {
      my $this = $blocks[$i];
      my $prev = $blocks[$i-1];

      if ($this->reference_genomic_align->dnafrag_start < 
          $prev->reference_genomic_align->dnafrag_end) {
        throw("Error when calculating simple net; overlapping blocks");
      }
    }
  }

  return \@blocks;
}



sub stringify_chains {
  my ($chains) = @_;

  my $str = "";

  foreach my $c (@$chains) {
    $str .= sprintf "CHAIN %d\n", $c->[0]->score;
    foreach my $b (@$c) {
      my $q = $b->reference_genomic_align;
      my ($t) = @{$b->get_all_non_reference_genomic_aligns};
      $str .= sprintf(" %s/%d-%d, %s/%d-%d %d (%d)\n", 
                      $q->dnafrag->name, 
                      $q->dnafrag_start, 
                      $q->dnafrag_end, 
                      $t->dnafrag->name, 
                      $t->dnafrag_start, 
                      $t->dnafrag_end, 
                      $t->dnafrag_strand, 
                      $t->level_id);
    }
  }


  return $str;
}



1;
