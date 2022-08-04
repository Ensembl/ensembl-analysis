# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Tools::WGA2Genes::GeneScaffold 

=head1 SYNOPSIS

A Bio::EnsEMBL::Slice that is comprised
of pieces of different target sequences, inferred by an alignment
of the target to a finished, query genome sequence. This object
extends Slice with mappings to/from the query and target

Assumptions:

- that the given GenomicAlignBlocks are sorted with respect
  to the reference and non-overlapping on query and target

- that the given list of features is sorted with respect to
  the reference

=cut


package Bio::EnsEMBL::Analysis::Tools::WGA2Genes::GeneScaffold;

use warnings ;
use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);

use Bio::EnsEMBL::Slice;

use Bio::EnsEMBL::Analysis::Tools::WGA2Genes::CoordUtils;

@ISA = qw(Bio::EnsEMBL::Slice);

my $FROM_CS_NAME = 'chromosome';
my $TO_CS_NAME   = 'scaffold';
my $GENE_SCAFFOLD_CS_NAME  = 'genescaffold';

my $INTERPIECE_PADDING      = 100;
my $NEAR_CONTIG_END         = 15;
###############################################

sub new {
  my ($caller, %given_args) = @_;

  my $class = ref($caller) || $caller;

  my ($name,
      $genomic_align_blocks,
      $transcripts,
      $from_slice,
      $to_slices,
      $max_readthrough_dist,
      $extend_into_gaps,
      $add_gaps,
      $direct_target_slice,
      ) = rearrange([qw(
                        NAME
                        GENOMIC_ALIGN_BLOCKS
                        TRANSCRIPTS
                        FROM_SLICE
                        TO_SLICES
                        MAX_READTHROUGH_DIST
                        EXTEND_INTO_GAPS
                        ADD_GAPS
                        DIRECT_TARGET_COORDS
                        )], %given_args);

  my $coding_transcripts = _check_transcripts($transcripts);

  $name = "GeneScaffold" if not defined $name;
  $max_readthrough_dist = 15 if not defined $max_readthrough_dist;  

  if ($direct_target_slice) {
    if ($add_gaps or $extend_into_gaps) {
      warning("GeneScaffold: cannot use -direct_target_coords with either of " . 
              "-add_gaps or -extend_into_gaps; unsetting both");
      $add_gaps = 0;
      $extend_into_gaps = 0;
    }

    # check that the genomic align blocks will allow this
    $direct_target_slice = _check_direct_target_coordinates($genomic_align_blocks,
                                                           $to_slices);
  }

  my $aln_map = _make_alignment_mapper($genomic_align_blocks);

  my ($gs_seq, $from_mapper, $to_mapper) = 
      _construct_sequence($genomic_align_blocks,
                          $aln_map,
                          $coding_transcripts,
                          $from_slice,
                          $to_slices,
                          $max_readthrough_dist,
                          $extend_into_gaps,
                          $add_gaps,
                          $direct_target_slice);

  return undef if not defined $from_mapper;

#  print "1: ", $direct_target_slice,"\n2: ",$direct_target_slice->length,"\n3: ",length($gs_seq),"\n";

  my $gs_end = $direct_target_slice 
      ? $direct_target_slice->length
      : length($gs_seq);


  my $self = $class->SUPER::new(-coord_system => 
                                   Bio::EnsEMBL::CoordSystem->new(-name => $GENE_SCAFFOLD_CS_NAME,
                                                                  -rank => 1),
                                -seq_region_name => $name,
                                -seq => ( $gs_seq eq "" ? undef : $gs_seq),
                                -start => 1,
                                -end   => $gs_end,
                                );


  $self->direct_target_slice($direct_target_slice)
      if defined $direct_target_slice;  
  $self->max_readthrough_dist($max_readthrough_dist)
      if defined $max_readthrough_dist;
  $self->from_slice($from_slice);
  $self->to_slices($to_slices);
  $self->alignment_mapper($aln_map);
  $self->from_mapper($from_mapper);
  $self->to_mapper($to_mapper);

  return $self;
}


###################################################################
# FUNCTION   : place_transcript
#
# Description:
#    Takes a transcript, and uses the mapping between 
#    query coords and gene scaffold coords to produces a transcript 
#    that is the result of "projecting" the original transcript, 
#    through alignment, onto the gene scaffold. 
###################################################################

sub place_transcript {
  my ($self, 
      $tran,
      $add_attributes,
      $external_db_id) = @_;


  my (@all_coords, @new_exons);
 
  # Some human protein_coding genes contain both protein_coding and non-coding transcripts.  This "if" statement
  # is to catch those non-coding transcripts so they don't get placed in the projected genome.  Otherwise, the code
  # would die in the next few lines as we can't get translateable exons from non-coding transcripts.

  if(! $tran->translation){
    print "Transcript ".$tran->stable_id." doesn't have a translation. Not using it in projection.\n"; 
    return undef;             
  }
   
  my @orig_exons = @{$tran->get_all_translateable_Exons};

  @orig_exons = sort { $a->start <=> $b->start } @orig_exons;
  
  my $orig_exon_coords = transcripts_to_coords($self->alignment_mapper,
                                               $FROM_CS_NAME,
                                               $tran);

  
  my $restricted_range;
  if (@$orig_exon_coords) {

   # print "ORIGINAL EXON COORDS: First start: ",$orig_exon_coords->[0]->start, " last end: ",$orig_exon_coords->[-1]->end," length: ",$orig_exon_coords->[-1]->end-$orig_exon_coords->[0]->start+1 ,"\n";
    $restricted_range = Bio::EnsEMBL::Mapper::Coordinate->new($orig_exon_coords->[0]->id,
                                                              $orig_exon_coords->[0]->start,
                                                              $orig_exon_coords->[-1]->end,
                                                              $orig_exon_coords->[0]->strand);
  } else {

    #print "ORIGINAL EXON COORDS: First start: ",$orig_exons[0]->start, " last end: ",$orig_exons[-1]->end," length ",$orig_exons[-1]->end-$orig_exons[0]->start+1,"\n";
    $restricted_range = Bio::EnsEMBL::Mapper::Coordinate->new($orig_exons[0]->slice->seq_region_name,
                                                              $orig_exons[0]->start,
                                                              $orig_exons[-1]->end,
                                                              $orig_exons[0]->strand);
  }
    
  my %filled_coords;
  
  # look for regions that map down to filled gaps. If any of these
  # lie outside the restricted range (i.e. extent of the transcript
  # within which gaps were potentially filled), then they must 
  # be gaps that were filled for a different transcript. In that
  # case, discard these pieces

  foreach my $orig_exon (@orig_exons) {

    #print "Single exon original coords:  START: ", $orig_exon->start," END: ", $orig_exon->end," length ",$orig_exon->end-$orig_exon->start+1,"\n";
    my @crds = $self->from_mapper->map_coordinates($orig_exon->slice->seq_region_name,
                                                   $orig_exon->start,
                                                   $orig_exon->end,
                                                   1,
                                                   $FROM_CS_NAME);

    foreach my $c (@crds) {
      if ($c->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
      
        #print "Exon Mapper from coordinates: START:", $c->start, " END: ",$c->end," length ",$c->end-$c->start+1,"\n"; 
        my ($tc) = $self->to_mapper->map_coordinates($GENE_SCAFFOLD_CS_NAME,
                                                     $c->start,
                                                     $c->end,
                                                     1,
                                                     $GENE_SCAFFOLD_CS_NAME);

        if ($tc->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
          # this piece may partially extend into a real sequence gap;
          # break it up into its constituent pieces    

          #print "TC Mapper from coordinates: START:", $tc->start, " END: ",$tc->end," length ",$tc->end-$tc->start+1,"\n"; 
          my @alncrds = $self->alignment_mapper->map_coordinates($tc->id,
                                                                 $tc->start,
                                                                 $tc->end,
                                                                 $tc->strand,
                                                                 $TO_CS_NAME);
          my $current_start = $c->start;
          foreach my $ac (@alncrds) {
            #print "AC coords START: ", $ac->start, "  END: ".$ac->end," length ",$ac->end-$ac->start+1,"\n";

            my $current_end = $current_start + $ac->length - 1;
            my $newc = Bio::EnsEMBL::Mapper::Coordinate->new($c->id,
                                                             $current_start,
                                                             $current_end,
                                                             $c->strand);
            $current_start += $ac->length;

            if ($ac->isa("Bio::EnsEMBL::Mapper::Gap")) {
              #print "\nAC is a GAP\n\n";
              $filled_coords{$newc} = 1;
            }
            push @all_coords, $newc;
          }
        } else {
          $filled_coords{$c} = 1;
          push @all_coords, $c;
        }
      } else {
        push @all_coords, $c;
      }
    }
  }
  
  my @kept_coords;
  foreach my $c (@all_coords) {
    if ($c->isa("Bio::EnsEMBL::Mapper::Coordinate") and 
                exists $filled_coords{$c}) {
      # this is a filled coord, so MUST correspond to exactly one piece of the
      # original query sequence
      #print "KEPT coord START: ",$c->start,"  END: ",$c->end," length ",$c->end-$c->start+1,"\n";

      my ($oc) = $self->from_mapper->map_coordinates($GENE_SCAFFOLD_CS_NAME,
                                                     $c->start,
                                                     $c->end,
                                                     1,
                                                     $GENE_SCAFFOLD_CS_NAME);
      if ($oc->start >= $restricted_range->start and
          $oc->end <= $restricted_range->end) {
        push @kept_coords, $c;
      } else {
       
        push @kept_coords, Bio::EnsEMBL::Mapper::Gap->new(1, $c->length);
      }
    } else {
      #print "KEPT coord START: ",$c->start,"  END: ",$c->end," length ",$c->end-$c->start+1,"\n";
      push @kept_coords, $c;
    }
  }
  @all_coords = @kept_coords;
  

  #
  # now massage the gaps to account for frameshifts
  #
  my $need_another_pass;
  my $loop_counter = 0;
  do {

    $need_another_pass = 0;

    # loop_counter is just a hack for a stuff that is getting stuck
    # This is an optimisation issue that affects a very low number of transcripts, but it does cause the pipeline to
    # before stuck. This is a poor solution but again time dictates that it will go in for the foreseeable future
    # On the plus side it seems to work fine as a hacky solution
    # 5000 is just some arbitrary number
    if($loop_counter >= 5000) {
      warning("Projection of transcript has become stuck, skipping");
      return undef;
    }
    $loop_counter++;

    my (@proc_coords, @gap_indices);
    # merge gaps
    foreach my $c (@all_coords) {
      if ($c->isa("Bio::EnsEMBL::Mapper::Gap")) {
         #print "GAP PRESENT START: ",$c->start,"  END: ",$c->end," length ",$c->end-$c->start+1,"\n";


        if (@proc_coords and $proc_coords[-1]->isa("Bio::EnsEMBL::Mapper::Gap")) {
          $proc_coords[-1]->end( $proc_coords[-1]->end + $c->length );
          #print "Extending existing GAP\n";
        } else {
          #print "ADDing new GAP with length ",$c->length ,"\n";
          push @proc_coords, Bio::EnsEMBL::Mapper::Gap->new(1, $c->length);
          push @gap_indices, scalar(@proc_coords) - 1;
        }
      } else {
        push @proc_coords, $c;
      }
    }
    
    GAP: foreach my $idx (@gap_indices) {
      #print "I have to handle a GAP\n";
      my $gap = $proc_coords[$idx];
      my $frameshift = $gap->length % 3;
     
      if ($frameshift) {
        #print "!!! Have a frameshift at " . $gap->start . "\n";
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
  
  my $start_not_found = $all_coords[0]->isa("Bio::EnsEMBL::Mapper::Gap") ? 1 : 0;
  my $end_not_found = $all_coords[-1]->isa("Bio::EnsEMBL::Mapper::Gap") ? 1 : 0;



  #
  # The remaining non-gap pieces are the exons
  #
  @all_coords = grep { $_->isa("Bio::EnsEMBL::Mapper::Coordinate") } @all_coords;
  return 0 if not @all_coords;

  my $sl =  $self->direct_target_slice
       ? $self->direct_target_slice
       : $self;

  foreach my $coord (@all_coords) {
    push @new_exons, Bio::EnsEMBL::Exon->new(-start => $coord->start,
                                             -end   => $coord->end,
                                             -strand => $tran->strand,
                                             -slice => $sl);
  }

  if ($tran->strand < 0) {
    @new_exons = sort { $b->start <=> $a->start } @new_exons;
  } else {
    @new_exons = sort { $a->start <=> $b->start } @new_exons;
  }
  
  #
  # calculate phases, and add supporting features
  #
  my $source_pep = $tran->translate->seq;
  my $source_cds_len = 0; 
  map { $source_cds_len  += $_->length } @orig_exons;
  if (length($source_pep) + 1 == ($source_cds_len / 3)) {
    # stop was removed
    $source_pep .= "*";
  }

  my $transcript_aligned_aas = 0;
  my $transcript_identical_aas = 0;
  my $incomplete_codon_bps = 0;

  my ($previous_exon);

  foreach my $exon (@new_exons) {
		#	if ($exon->start == 16141929) {
		#		  print "Modifying exon->start\n";
		#			$exon->start($exon->start+2);
    #  }
    if (defined $previous_exon) {
      $exon->phase($previous_exon->end_phase);
    } else {
      $exon->phase(0);
    }
  #  $exon->end_phase((($exon->end - $exon->start + 1) + $exon->phase)%3);
	  #print " Exon length = " . $exon->length . " phase = " . $exon->phase . "\n";
    $exon->end_phase((($exon->length) + $exon->phase)%3);
	  #print " Setting Exon end_phase to  " . $exon->end_phase . "\n";

    my $exon_aligned_aas = 0;
    my $exon_identical_aas = 0;

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
    #print "Exon phase : ",$exon->phase,"\n";
    #print "Exon end phase : ",$exon->end_phase,"\n";
    #print "Previous exon end phase : ",$previous_exon->end_phase,"\n" if ($previous_exon);
    if ($extent_end > $extent_start) {
			#print "extent start = " . $extent_start . " extent_end = " . $extent_end ."\n";
      # if not, we've eaten away the whole exon, so there is no support
      $incomplete_codon_bps += ($extent_start - $exon->start); 
      $incomplete_codon_bps += ($exon->end - $extent_end); 

      my @gen_coords = $self->from_mapper->map_coordinates($GENE_SCAFFOLD_CS_NAME,
                                                           $extent_start,
                                                           $extent_end,
                                                           1,
                                                           $GENE_SCAFFOLD_CS_NAME);

      #print " Have " . scalar(@gen_coords) . " gen_coords\n";
      my @fps;
      my $cur_gs_start = $extent_start;
      #print "Exon coords = " . $exon->start . " " . $exon->end ."\n";
      foreach my $g_coord (@gen_coords) {
        #print " g_coord = " . $g_coord->start . " " . $g_coord->end ."\n";
        #print "EXTEND LENGTH: ",($extent_end - $extent_start+1),"\n";
        #print "G LENGTH: ",($g_coord->end- $g_coord->start+1),"\t length: ",$g_coord->length,"\n";

        my $cur_gs_end = $cur_gs_start + $g_coord->length - 1;
        
        if ($g_coord->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
          
          my ($p_coord) = $tran->genomic2pep($g_coord->start, 
                                             $g_coord->end,
                                             $tran->strand);

          

          my $p_substr = uc(substr($source_pep, $p_coord->start - 1, $p_coord->length));
          my $gs_reg = $sl->subseq($cur_gs_start, $cur_gs_end, $exon->strand);
          #print "SEQ: ",length($gs_reg),"\n";

          my $gs_p_substr = uc(Bio::PrimarySeq->new(-seq => $gs_reg)->translate->seq);

            #print $p_substr, "\n AND: ",$gs_p_substr,"\n";

          if (length($p_substr) != length ($gs_p_substr)) {
            #print "Length from source_pep = ", length($p_substr), "\n AND: from genomic translation ",length($gs_p_substr),"\n";
            warning("Pep segments differ; cannot calculate percentage identity");

          } else {
            $exon_aligned_aas += length($p_substr);

            my @p1 = split(//, $p_substr);
            my @p2 = split(//, $gs_p_substr);
            for(my $i=0; $i < @p1; $i++) {
              if ($p1[$i] eq $p2[$i]) {
                $exon_identical_aas++;
              }
            }
          }

          my $fp;
          my $target_feature_length =  $cur_gs_end - $cur_gs_start + 1;
          my $query_feature_length = $p_coord->end - $p_coord->start + 1;
          # There is an issue that crops up for a some features, where the query end is off
          # by 1 aa. I don't have time to fully check this now, for the moment I will just warn and correct the
          # sf here
          if($query_feature_length * 3 != $target_feature_length) {
            warning("Something has gone wrong with the feature creation, feature will fail to store so throwing.\n".
                    "Will shorted by 1 aa to try and compensate\n".
                    "Genomic feature length: ".$target_feature_length."\n".
                    "Query feature length: ".$query_feature_length."\n".
                    "Expected length ratio: 3 1\n".
                    "Query length time 3: ".($query_feature_length * 3)."\n".
                    "Exon start: ".$exon->start."\n".
                    "Exon end: ".$exon->end."\n".
                    "Exon strand: ".$exon->strand."\n".
                    "Gene segment start: ".$cur_gs_start."\n".
                    "Gene segment end: ".$cur_gs_end."\n".
                    "Hit start: ".$p_coord->start."\n".
                    "Hit end: ".$p_coord->end
                   );

              $fp = Bio::EnsEMBL::FeaturePair->
              new(-seqname  => $self->seq_region_name,
                  -start    => $cur_gs_start,
                  -end      => $cur_gs_end,
                  -strand   => $exon->strand,
                  -score    => 100.0,
                  -hseqname => $tran->translation->stable_id,
                  -hstart   => $p_coord->start,
                  -external_db_id => $external_db_id,
                  -hend     => ($p_coord->end - 1),
                  -hstrand  => $p_coord->strand,
                  -slice    => $sl);
          } else {
            $fp = Bio::EnsEMBL::FeaturePair->
              new(-seqname  => $self->seq_region_name,
                  -start    => $cur_gs_start,
                  -end      => $cur_gs_end,
                  -strand   => $exon->strand,
                  -score    => 100.0,
                  -hseqname => $tran->translation->stable_id,
                  -hstart   => $p_coord->start,
                  -external_db_id => $external_db_id, 
                  -hend     => $p_coord->end,
                  -hstrand  => $p_coord->strand,
                  -slice    => $sl);
         }
          push @fps, $fp;
        }
        
        $cur_gs_start += $g_coord->length;
      }
      
      if (@fps) {
        my $f = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@fps, -align_type => 'ensembl');
        if($exon_aligned_aas == 0) {
          $f->percent_id(0);
        } else {
          $f->percent_id(100 * ($exon_identical_aas / $exon_aligned_aas));
        }
        $exon->add_supporting_features($f);
      }
    } else {
      #print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Short exon of " . $exon->length . " bases\n";
      $incomplete_codon_bps += $exon->length;
    }
    
    $transcript_aligned_aas += $exon_aligned_aas;
    $transcript_identical_aas += $exon_identical_aas;
    
    $previous_exon = $exon;
  }
  # optimistically, we are going to assume that split codons
  # are identical and covered. This of course may not always
  # be the case, but it's too fiddly (for now) to do it properly. 
  $transcript_aligned_aas += int($incomplete_codon_bps / 3);
  $transcript_identical_aas += int($incomplete_codon_bps / 3);

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
            $intron_len <= $self->max_readthrough_dist) { 
          $new_start = $exon->start;
          $new_end   = $prev_exon->end;
        }
      } else {
        my $intron_len = $exon->start - $prev_exon->end - 1;
        if ($intron_len % 3 == 0 and 
            $intron_len <= $self->max_readthrough_dist) {
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
        my $len_l = 0;
        my $pid_l = 0;
        if (@{$prev_exon->get_all_supporting_features}) {
          my ($sf) = @{$prev_exon->get_all_supporting_features};
          $pid_l = $sf->percent_id;
          foreach my $ug ($sf->ungapped_features) {
            push @ug_feats, $ug;
            $len_l += $ug->length;
          }
        }
        my $len_r = 0;
        my $pid_r = 0;
        if (@{$exon->get_all_supporting_features}) {
          my ($sf) = @{$exon->get_all_supporting_features};
          $pid_r = $sf->percent_id;
          foreach my $ug ($sf->ungapped_features) {
            push @ug_feats, $ug;
            $len_r += $ug->length;
          }
        }
        if (@ug_feats) {
          map { $_->percent_id(100) } @ug_feats;

          my $new_sup_feat = Bio::EnsEMBL::DnaPepAlignFeature->
              new(-features => \@ug_feats, -align_type => 'ensembl');
          $new_sup_feat->percent_id((($len_l * $pid_l) + ($len_r * $pid_r)) 
                                    / ($len_l + $len_r));
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
    
  my $proj_tran = Bio::EnsEMBL::Transcript->new(-analysis => $tran->analysis);
  
  map { $proj_tran->add_Exon($_) } @merged_exons;
  

  #
  # do transcript-level supporting features/attributes
  #
  
  my (@trans_fps);
  foreach my $exon (@merged_exons) {
    
    if (@{$exon->get_all_supporting_features}) {
      my ($sf) = @{$exon->get_all_supporting_features};
      my @e_fps = $sf->ungapped_features;
      # need to reset the pids here, otherwise the API complains
      map { $_->percent_id(100) } @e_fps;
      push @trans_fps, @e_fps;
    }
  }

  if (@trans_fps) {
    my $t_sf = Bio::EnsEMBL::DnaPepAlignFeature->
        new(-features => \@trans_fps, -align_type => 'ensembl');
    # use score to hold coverage, so that it is stored somewhere when written to db
    $t_sf->score(100 * ($transcript_aligned_aas / length($source_pep)));
    $t_sf->hcoverage( 100 * ($transcript_aligned_aas / length($source_pep)));
    $t_sf->percent_id(100 * ($transcript_identical_aas / $transcript_aligned_aas));
    $proj_tran->add_supporting_features($t_sf);

    # now copy this hcoverage value to all the exon supporting features
    foreach my $exon (@{$proj_tran->get_all_Exons}) {
      foreach my $exon_sf (@{$exon->get_all_supporting_features}) {
        $exon_sf->hcoverage($t_sf->hcoverage); 
      }
    }

  } else {
    # there are no complete codons in any exon!
    # this is clearly rubbish, so bail out
    return 0;
  }
  

  #
  # set translation
  #
  my $translation = Bio::EnsEMBL::Translation->new();
  $translation->start_Exon($merged_exons[0]);
  $translation->start(1);
  $translation->end_Exon($merged_exons[-1]);
  $translation->end($merged_exons[-1]->end - $merged_exons[-1]->start + 1);
  
  $proj_tran->translation($translation);
  
  my $pep = $proj_tran->translate;
  
  if (not defined $pep) {
    # this can happen if the transcript comprises a single stop codon only
    return 0;
  }

  my $prop_non_gap = 100 - (100 * (($pep->seq =~ tr/X/X/) / $pep->length));
  my $num_stops = $pep->seq =~ tr/\*/\*/;
  my $num_gaps  = $pep->seq =~ tr/X/X/;

  #
  # finally, attributes
  #
  my @attributes;
  
  my $perc_id_attr = Bio::EnsEMBL::Attribute->
      new(-code => 'HitSimilarity',
          -name => 'hit similarity',
          -description => 'percentage id to parent transcripts',
          -value => sprintf("%.1f",
                            100 * ($transcript_identical_aas / $transcript_aligned_aas)));
  push @attributes, $perc_id_attr;

  my $cov_attr = Bio::EnsEMBL::Attribute->
      new(-code => 'HitCoverage',
          -name => 'hit coverage',
          -description => 'coverage of parent transcripts',
          -value => sprintf("%.1f",
                            100 * ($transcript_aligned_aas / length($source_pep))));
  push @attributes, $cov_attr;
  
  my $gap_attr = Bio::EnsEMBL::Attribute->
      new(-code => 'PropNonGap',
          -name => 'proportion non gap',
          -description => 'proportion non gap',
          -value => sprintf("%.1f", 
                            $prop_non_gap));
  push @attributes, $gap_attr;
  
  my $prop_non_gap_of_coverage_attr = Bio::EnsEMBL::Attribute->
    new(-code => 'NonGapHCov',
        -name => 'proportion non gap of hit coverage',
        -description => 'proportion non gap of hit coverage',
        -value => sprintf("%.1f",
                          100 * (($transcript_aligned_aas - $num_gaps) / length($source_pep))));
  push @attributes, $prop_non_gap_of_coverage_attr;

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
    
  my $stop_attr = Bio::EnsEMBL::Attribute->
      new(-code => 'NumStops',
          -name => 'number of stops',
          -desc => 'Number of stops before editing',
          -value => $num_stops);
  push @attributes, $stop_attr;
  
  # indentify gap exons
  my $gap_exons = 0;
  foreach my $e (@{$proj_tran->get_all_Exons}) {
    my ($coord) = $self->to_mapper->map_coordinates($GENE_SCAFFOLD_CS_NAME,
                                                    $e->start,
                                                    $e->end,
                                                    1,
                                                    $GENE_SCAFFOLD_CS_NAME);
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
  
  my $tranid_attr = Bio::EnsEMBL::Attribute->
      new(-code => 'SourceTran',
          -name => 'source transcript',
          -description => 'source transcript',
          -value => $tran->stable_id);
  push @attributes, $tranid_attr;

  if ($add_attributes) {
    $proj_tran->add_Attributes(@attributes);
  }
  if ($self->direct_target_slice and
      $self->direct_target_slice->strand < 0) {
    $proj_tran = $proj_tran->transfer($self->direct_target_slice->invert);
  }
  return $proj_tran;
}


sub stringify_alignment {
  my ($self) = @_;

  my @coords = $self->alignment_mapper->map_coordinates($self->from_slice->seq_region_name,
                                                        $self->from_slice->start,
                                                        $self->from_slice->end,
                                                        1,
                                                        $FROM_CS_NAME);
  my $string = "";
  my $position = $self->from_slice->start;
  foreach my $c (@coords) {
    my $ref_start = $position;
    my $ref_end = $position + $c->length - 1;
    $position += $c->length;
    
    $string .= sprintf("%s %d %d ", $self->from_slice->seq_region_name, $ref_start, $ref_end);

    if ($c->isa("Bio::EnsEMBL::Mapper::Gap")) {
      $string .= "[GAP]\n";
    } else {
      $string .= "[%s %d %d %s]\n", $c->id, $c->start, $c->end, $c->strand;
    }
  }

  return $string;
}

sub project_up {
  my ($self) = @_;

  my $tlsl = $self->from_slice;

  my @comps = $self->from_mapper->map_coordinates($GENE_SCAFFOLD_CS_NAME,
                                                  1,
                                                  $self->length,
                                                  1,
                                                  $GENE_SCAFFOLD_CS_NAME);
                                                  
  my @segments;
  
  my $current_pos = 1;
  foreach my $c (@comps) {
    my $start = $current_pos;
    my $end   = $current_pos + $c->length - 1;
    $current_pos += $c->length;
    
    if ($c->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
      
      my $slice = $tlsl->adaptor->fetch_by_region($tlsl->coord_system->name,
                                                  $tlsl->seq_region_name,
                                                  $c->start,
                                                  $c->end,
                                                  1);
      
      push @segments, bless([$start, $end, $slice],
                            "Bio::EnsEMBL::ProjectionSegment");
    }
  }

  return @segments;
}

sub project_down {
  my ($self) = @_;

  my @comps = $self->to_mapper->map_coordinates($GENE_SCAFFOLD_CS_NAME,
                                                1,
                                                $self->length,
                                                1,
                                                $GENE_SCAFFOLD_CS_NAME);
                                                  
  my @segments;
  
  my $current_pos = 1;
  foreach my $c (@comps) {
    my $start = $current_pos;
    my $end   = $current_pos + $c->length - 1;
    $current_pos += $c->length;
    
    if ($c->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
      my $tlsl = $self->to_slices->{$c->id};
      
      my $slice = $tlsl->adaptor->fetch_by_region($tlsl->coord_system->name,
                                                  $tlsl->seq_region_name,
                                                  $c->start,
                                                  $c->end,
                                                  $c->strand);
      
      push @segments, bless([$start, $end, $slice],
                            "Bio::EnsEMBL::ProjectionSegment");
    }
  }

  return @segments;
}


###############################################
# Internal helper methods
###############################################

#################################################

sub _check_transcripts {
  my ($trans) = @_;

  my @coding_transcripts;
  my $transcript_is_good = 1;

  if (scalar(@$trans) == 0) {
    throw("Attempt to create GeneScaffold with empty transcript list");
  }

  foreach my $t (@$trans) {
    if (not $t->translation) {
     # throw("Attempt to create GeneScaffold with non-coding Transcript (".$t->stable_id.")");
warn("Attempt to create GeneScaffold with non-coding Transcript (".$t->stable_id.")");
       $transcript_is_good = 0;
    }
    if (length($t->translateable_seq) % 3 != 0) {
    #  throw("Attempt to create GeneScaffold with non-mod-3 coding length Transcript (".$t->stable_id.")");
warn("Attempt to create GeneScaffold with non-mod-3 coding length Transcript (".$t->stable_id.")");
       $transcript_is_good = 0;
    }
    if($transcript_is_good == 1){
      push (@coding_transcripts, $t);
    }
  }
  return \@coding_transcripts;
}


#################################################

sub _construct_sequence {
  my ($gen_al_blocks,
      $map,
      $transcripts,
      $from_slice,
      $to_slices,
      $max_readthrough_dist,
      $extend_into_gaps,
      $add_gaps,
      $direct_target_slice) = @_;

  # Basic gene-scaffold structure is taken directly from the given block list
  my @block_coord_pairs;
  foreach my $bl (@{$gen_al_blocks}) {
    my $qy_al = $bl->reference_genomic_align;
    my ($tg_al) = @{$bl->get_all_non_reference_genomic_aligns};

    my $from_coord = Bio::EnsEMBL::Mapper::Coordinate->new($qy_al->dnafrag->name,
                                                           $qy_al->dnafrag_start,
                                                           $qy_al->dnafrag_end,
                                                           1);
    my $to_coord   = Bio::EnsEMBL::Mapper::Coordinate->new($tg_al->dnafrag->name,
                                                           $tg_al->dnafrag_start,
                                                           $tg_al->dnafrag_end,
                                                           $tg_al->dnafrag_strand);
    my $pair = Bio::EnsEMBL::Mapper::Pair->new($from_coord, 
                                               $to_coord);
    push @block_coord_pairs, $pair;
  }
  @block_coord_pairs = sort { $a->from->start <=> $b->from->start } @block_coord_pairs;

  # we now proceed to amend the structure with inserted gaps

  # step 1: flatten the exons from the transcripts into a non-overlapping
  # list, ommitting terminal exons that map completely to gaps
  my $tran_coords = transcripts_to_coords($map,
                                          $FROM_CS_NAME,
                                          @$transcripts);
  

  # step 2: infer list of exon regions that map to gaps. We are only 
  # interested, at this stage, in regions outside the blocks, because
  # these are the ones with potential to be "filled"
  my @fillable_gaps;

  foreach my $tc (sort { $a->start <=> $b->start } @$tran_coords) {
    my $current_pos = $tc->start;
    
    foreach my $c ($map->map_coordinates($tc->id,
                                         $tc->start,
                                         $tc->end,
                                         1,
                                         $FROM_CS_NAME)) {
      my $from_coord = Bio::EnsEMBL::Mapper::Coordinate->new($tc->id,
                                                             $current_pos,
                                                             $current_pos + $c->length - 1,
                                                             1);
      
      
      $current_pos += $c->length;
      
      if ($c->isa("Bio::EnsEMBL::Mapper::Gap")) {
        # only consider gaps that lies outside blocks
        my $overlaps_block = 0;
        foreach my $bl (@block_coord_pairs) {
          if ($bl->from->start <= $from_coord->end and
              $bl->from->end   >= $from_coord->start) {
            $overlaps_block = 1;
            last;
          }
        }
        if (not $overlaps_block) {
          push @fillable_gaps, Bio::EnsEMBL::Mapper::Pair->new($from_coord,
                                                               $c);
        }
      }
    }
  }

  # Non-fillable gaps:
  # 1. Gaps before the first or after the last block
  # 2. If one of the flanking coords is on the same CDS region
  #    as the gap, and the end of the aligned region does
  #    not align to a sequence-level gap
  # 3. If the 2 flanking coords are consistent and not
  #    separated by a sequence-level gap
  # 

  my @all_coord_pairs = (@block_coord_pairs, @fillable_gaps);
  @all_coord_pairs = sort { $a->from->start <=> $b->from->start } @all_coord_pairs;

  my @extend_pairs;

  if ($extend_into_gaps) {
    my (%pairs_to_remove, @replacements);
    
    for(my $i=0; $i < @all_coord_pairs; $i++) {
      my $this_pair = $all_coord_pairs[$i];
      
      if ($this_pair->to->isa("Bio::EnsEMBL::Mapper::Gap")) {
        # if it's gap that can be filled, leave it. Otherwise, remove it
        my ($left_non_gap, $right_non_gap);
        for(my $j=$i-1; $j>=0; $j--) {
          if ($all_coord_pairs[$j]->to->
              isa("Bio::EnsEMBL::Mapper::Coordinate")) {
            $left_non_gap = $all_coord_pairs[$j];
            last;
          }
        }
        for(my $j=$i+1; $j < @all_coord_pairs; $j++) {
          if ($all_coord_pairs[$j]->to->
              isa("Bio::EnsEMBL::Mapper::Coordinate")) {
            $right_non_gap = $all_coord_pairs[$j];
            last;
          }
        }
        
        
        my ($ex_left, $ex_left_up, $ex_left_down) = 
            extend_coord($left_non_gap->to,
                         $to_slices->{$left_non_gap->to->id});
        my ($ex_right, $ex_right_up, $ex_right_down) = 
            extend_coord($right_non_gap->to,
                         $to_slices->{$right_non_gap->to->id});

        # flanking coords are inconsistent,
        # which means that they come from different chains.
        # By chain filtering then, they must either come
        # from different target sequences, or be separable in the
        # same target sequence by a sequence-level gap. Either
        # way, we can nominally "fill" the gap. 
        #
        # However, if the gap coincides with the end of a block,
        # and furthermore if the block end conincides with the
        # end of a sequence-level piece in the target, it is
        # more appropriate to extend the exon into the existing
        # gap; in that case, we replace the gap with a fake
        # piece of alignment
        
        if (not check_consistent_coords($left_non_gap->to,
                                        $right_non_gap->to) or
            $ex_left->start > $ex_right->end or
            $ex_left->end   < $ex_right->start) {
          
          if ($left_non_gap->from->end == $this_pair->from->start - 1 and
              $right_non_gap->from->start == $this_pair->from->end + 1) {
            
            my $remove_coord = 1;
            my (@replace_coord);
            
            if ($left_non_gap->to->strand > 0 and 
                $ex_left->end - $left_non_gap->to->end <= $NEAR_CONTIG_END) {
              $remove_coord = 0;
              
              if (defined $ex_left_down and
                  $this_pair->to->length <= $ex_left_down->start - $ex_left->end - 1) {            
                push @replace_coord, Bio::EnsEMBL::Mapper::Coordinate
                    ->new($ex_left->id,
                          $ex_left->end + 1,
                          $ex_left->end + $this_pair->to->length,
                          $ex_left->strand);                  
              }
            } elsif ($left_non_gap->to->strand < 0 and 
                     $ex_left->start - $left_non_gap->to->start <= $NEAR_CONTIG_END) {
              $remove_coord = 0;
              
              if (defined $ex_left_up and 
                  $this_pair->to->length <= $ex_left->start - $ex_left_up->end - 1) {
                push @replace_coord, Bio::EnsEMBL::Mapper::Coordinate
                    ->new($ex_left->id,
                          $ex_left->start - $this_pair->to->length,
                          $ex_left->start - 1,
                          $ex_left->strand);
              }            
            } 
            
            if ($right_non_gap->to->strand > 0 and
                $right_non_gap->to->start - $ex_right->start <= $NEAR_CONTIG_END) {
              $remove_coord = 0;
              
              if (defined $ex_right_up and
                  $this_pair->to->length <= $ex_right->start - $ex_right_up->end - 1) {
                push @replace_coord, Bio::EnsEMBL::Mapper::Coordinate
                    ->new($ex_right->id,
                          $ex_right->start - $this_pair->to->length,
                          $ex_right->start - 1,
                          $ex_right->strand);
              }
            } elsif ($right_non_gap->to->strand < 0 and
                     $ex_right->end - $right_non_gap->to->end <= $NEAR_CONTIG_END) {
              $remove_coord = 0;
              
              if (defined $ex_right_down and
                  $this_pair->to->length <= $ex_right_down->start - $ex_right->end - 1) {
                push @replace_coord, Bio::EnsEMBL::Mapper::Coordinate
                    ->new($ex_right->id,
                          $ex_right->end + 1,
                          $ex_right->end + $this_pair->to->length,
                          $ex_right->strand);
              }
            } 
            
            if ($remove_coord) {
              # gap does not align with the end of a contig; junk it
              $pairs_to_remove{$this_pair} = 1;
            } elsif (@replace_coord) {
              # arbitrarily chose the first one
              push @replacements, [$this_pair,
                                   $replace_coord[0]];
            }
            
          } elsif ($left_non_gap->from->end == $this_pair->from->start - 1) {
            if ($left_non_gap->to->strand > 0 and 
                $ex_left->end - $left_non_gap->to->end <= $NEAR_CONTIG_END) {
              
              if (defined $ex_left_down and
                  $this_pair->to->length <= $ex_left_down->start - $ex_left->end - 1) {
                push @replacements, [$this_pair,
                                     Bio::EnsEMBL::Mapper::Coordinate
                                     ->new($ex_left->id,
                                           $ex_left->end + 1,
                                           $ex_left->end + $this_pair->to->length,
                                           $ex_left->strand)];
              }
            } elsif ($left_non_gap->to->strand < 0 and 
                     $ex_left->start - $left_non_gap->to->start <= $NEAR_CONTIG_END) {
              
              if (defined $ex_left_up and 
                  $this_pair->to->length <= $ex_left->start - $ex_left_up->end - 1) {
                push @replacements, [$this_pair,
                                     Bio::EnsEMBL::Mapper::Coordinate
                                     ->new($ex_left->id,
                                           $ex_left->start - $this_pair->to->length,
                                           $ex_left->start - 1,
                                           $ex_left->strand)];
              }
            } else {
              # gap does not align with the end of a contig; junk it
              $pairs_to_remove{$this_pair} = 1;
            }
          } elsif ($right_non_gap->from->start == $this_pair->from->end + 1) {
            if ($right_non_gap->to->strand > 0 and 
                $right_non_gap->to->start - $ex_right->start <= $NEAR_CONTIG_END) {
              
              if (defined $ex_right_up and
                  $this_pair->to->length <= $ex_right->start - $ex_right_up->end - 1) {
                push @replacements, [$this_pair,
                                     Bio::EnsEMBL::Mapper::Coordinate
                                     ->new($ex_right->id,
                                           $ex_right->start - $this_pair->to->length,
                                           $ex_right->start - 1,
                                           $ex_right->strand)];
              }
            } elsif ($right_non_gap->to->strand < 0 and 
                     $ex_right->end - $right_non_gap->to->end <= $NEAR_CONTIG_END) {
              
              if (defined $ex_right_down and 
                  $this_pair->to->length <= $ex_right_down->start - $ex_right->end - 1) {
                push @replacements, [$this_pair, 
                                     Bio::EnsEMBL::Mapper::Coordinate
                                     ->new($ex_right->id,
                                           $ex_right->end + 1,
                                           $ex_right->end + $this_pair->to->length,
                                           $ex_right->strand)];
              }
            } else {
              # gap does not align with the end of a contig; junk it
              $pairs_to_remove{$this_pair} = 1;
            }
          }
          # else this gap is an isolate. It can be kept iff the coords are
          # on different chains, or on the same chain but on different
          # contigs; we've already determined that the components can
          # be separated, so fine          
          
        }
        else {
          $pairs_to_remove{$this_pair} = 1;
        }
      }
    }
    
    @all_coord_pairs = grep { not exists $pairs_to_remove{$_} } @all_coord_pairs;
    
    @replacements = sort { 
      $a->[1]->id cmp $b->[1]->id or
          $a->[1]->start <=> $b->[1]->start;    
    } @replacements;
    
    while(@replacements) {
      my $el = shift @replacements;
      my ($pair, $rep) = @$el;
      
      my $fill = 1;
      
      if (@replacements) {
        my $nel = shift @replacements;
        my ($npair, $nrep) = @$nel;
        
        if ($rep->id eq $nrep->id and
            $rep->start <= $nrep->end and
            $rep->end   >= $nrep->start) {
          # we have over-filled a gap. Remove these fills
          $fill = 0;
        } else {
          unshift @replacements, $nel;
        }
      }
      if ($fill) {
        $pair->to($rep);      
        push @extend_pairs, $pair;
      }
    }
  }

  if (not $add_gaps) {
    @all_coord_pairs = grep { not $_->to->isa("Bio::EnsEMBL::Mapper::Gap") } @all_coord_pairs;
  }

  # merge adjacent targets   
  #  we want to be able to account for small, frame-preserving 
  #  insertions in the target sequence with respect to the query. 
  #  To give the later, gene-projection code the opportunity to 
  #  "read through" these insertions, we have to merge togther 
  #  adjacent, consistent targets that are within this "maximum 
  # read-through" distance

  my @merged_pairs;

  for(my $i=0; $i<@all_coord_pairs; $i++) {
    my $this_pair = $all_coord_pairs[$i];

    if ($this_pair->to->isa("Bio::EnsEMBL::Mapper::Coordinate") and
        @merged_pairs and
        $merged_pairs[-1]->to->isa("Bio::EnsEMBL::Mapper::Coordinate") and 
        check_consistent_coords($merged_pairs[-1]->to,
                                $this_pair->to)) {

      my $dist = distance_between_coords($merged_pairs[-1]->to,
                                         $this_pair->to);
      
      if ($dist <= $max_readthrough_dist) {
        
        my $last_pair = pop @merged_pairs;

        my $new_from = merge_coords($last_pair->from,
                                    $this_pair->from);
        my $new_to = merge_coords($last_pair->to,
                                  $this_pair->to);
        
        # check that the new merged coord will not result in an overlap
        my $overlap = 0;
        foreach my $tg (@merged_pairs) {
          if ($tg->to->isa("Bio::EnsEMBL::Mapper::Coordinate") and
              $tg->to->id eq $new_to->id and
              $tg->to->start < $new_to->end and
              $tg->to->end   > $new_to->start) {
            $overlap = 1;
            last;
          }
        }
        if (not $overlap) {
          for (my $j=$i+1; $j < @all_coord_pairs; $j++) {
            my $tg = $all_coord_pairs[$j];
            if ($tg->to->isa("Bio::EnsEMBL::Mapper::Coordinate") and
                $tg->to->id eq $new_to->id and 
                $tg->to->start < $new_to->end and
                $tg->to->end   > $new_to->start) {
              $overlap = 1;
              last;
            }
          }
        }
        
        if ($overlap) {
          push @merged_pairs, $last_pair, $this_pair;
        } else {
          push @merged_pairs, Bio::EnsEMBL::Mapper::Pair->new($new_from, 
                                                                $new_to);
        }
      } else {
        push @merged_pairs, $this_pair;
      }
    }
    else {
      push @merged_pairs, $this_pair;
    }
  }
    
  #########################################################

  my $t_map = Bio::EnsEMBL::Mapper->new($TO_CS_NAME,
                                        $GENE_SCAFFOLD_CS_NAME); 
  my $q_map = Bio::EnsEMBL::Mapper->new($FROM_CS_NAME,
                                        $GENE_SCAFFOLD_CS_NAME);
  
  my ($seq, $last_end_pos) = ("", 0);
  for(my $i=0; $i < @merged_pairs; $i++) {
    my $pair = $merged_pairs[$i];

    if ($direct_target_slice) {
      if ($pair->to->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
        my ($gs_start, $gs_end);
        if ($direct_target_slice->strand < 0) {
          $gs_start = $direct_target_slice->length - $pair->to->end + 1;
          $gs_end   = $direct_target_slice->length - $pair->to->start + 1;
        } else {
          $gs_start = $pair->to->start;
          $gs_end = $pair->to->end;
        }
        
        $t_map->add_map_coordinates($pair->to->id,
                                    $pair->to->start,
                                    $pair->to->end,
                                    $pair->to->strand,
                                    $GENE_SCAFFOLD_CS_NAME,
                                    $gs_start,
                                    $gs_end);

      }
    } else{
      if ($pair->to->isa("Bio::EnsEMBL::Mapper::Coordinate")) {            
        # the sequence itself        
        my $slice = $to_slices->{$pair->to->id};
        my $this_seq = $slice->subseq($pair->to->start, $pair->to->end);
        if ($pair->to->strand < 0) {
          reverse_comp(\$this_seq);
        }      
        $seq .= $this_seq;
        
        $t_map->add_map_coordinates($pair->to->id,
                                    $pair->to->start,
                                    $pair->to->end,
                                    $pair->to->strand,
                                    $GENE_SCAFFOLD_CS_NAME,
                                    $last_end_pos + 1,
                                    $last_end_pos + $pair->to->length);
      } else {
        # the sequence itself
        $seq .= ('n' x $pair->from->length);
        
        # and the map. This is a target gap we have "filled", so no position 
        # in target, but a position in query
        
        $q_map->add_map_coordinates($from_slice->seq_region_name,
                                    $pair->from->start,
                                    $pair->from->end,
                                    1,
                                    $GENE_SCAFFOLD_CS_NAME,
                                    $last_end_pos + 1,
                                    $last_end_pos + $pair->from->length);
      }

      # add padding between the pieces
      if ($i < @merged_pairs - 1) {
        $last_end_pos += 
            $pair->to->length + $INTERPIECE_PADDING;
        $seq .= ('n' x $INTERPIECE_PADDING);
      }
    }
  }

  #
  # add the gaps we have extended into to the query map
  #
  foreach my $pair (@extend_pairs) {
    my ($coord) = $t_map->map_coordinates($pair->to->id,
                                          $pair->to->start,
                                          $pair->to->end,
                                          $pair->to->strand,
                                          $TO_CS_NAME);
    $q_map->add_map_coordinates($from_slice->seq_region_name,
                                $pair->from->start,
                                $pair->from->end,
                                1,
                                $GENE_SCAFFOLD_CS_NAME,
                                $coord->start,
                                $coord->end);
  }
  
  #
  # finally add of the original alignment pieces to the query map
  #
  my @aln_coords = $map->map_coordinates($from_slice->seq_region_name,
                                         $from_slice->start,
                                         $from_slice->end,
                                         1,
                                         $FROM_CS_NAME);
  my $current_start = $from_slice->start;
  foreach my $c (@aln_coords) {
    my $current_end = $current_start + $c->length - 1;

    if ($c->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
      # get the gene_scaffold position from the target map
      my ($coord) = $t_map->map_coordinates($c->id,
                                            $c->start,
                                            $c->end,
                                            $c->strand,
                                            $TO_CS_NAME);

      $q_map->add_map_coordinates($from_slice->seq_region_name,
                                  $current_start,
                                  $current_end,
                                  1,
                                  $GENE_SCAFFOLD_CS_NAME,
                                  $coord->start,
                                  $coord->end);
    }
    $current_start += $c->length;
  }

  return ($seq, $q_map, $t_map);
}

#################################################

sub _make_alignment_mapper {
  my ($gen_al_blocks) = @_;

  my $mapper = Bio::EnsEMBL::Mapper->new($FROM_CS_NAME,
                                         $TO_CS_NAME);

  foreach my $bl (@$gen_al_blocks) {
    foreach my $ugbl (@{$bl->get_all_ungapped_GenomicAlignBlocks}) {      
      my ($from_bl) = $ugbl->reference_genomic_align;
      my ($to_bl)   = @{$ugbl->get_all_non_reference_genomic_aligns};

      $mapper->add_map_coordinates($from_bl->dnafrag->name,
                                   $from_bl->dnafrag_start,
                                   $from_bl->dnafrag_end,
                                   $from_bl->dnafrag_strand * $to_bl->dnafrag_strand,
                                   $to_bl->dnafrag->name,
                                   $to_bl->dnafrag_start,
                                   $to_bl->dnafrag_end);
    }
  }

  return $mapper;
}


##############################################

sub _check_direct_target_coordinates {
  my ($blocks, $slices) = @_;

  my ($tname, $tstrand);

  my @blocks = @$blocks;
  for(my $i=0; $i < @blocks; $i++) {
    my ($right_to) = @{$blocks[$i]->get_all_non_reference_genomic_aligns};
    if ($i > 0) { 
      my ($left_to) = @{$blocks[$i-1]->get_all_non_reference_genomic_aligns};

      if ($left_to->dnafrag->name ne $right_to->dnafrag->name or
          $left_to->dnafrag_strand != $right_to->dnafrag_strand or
          ($left_to->dnafrag_strand > 0 and $left_to->dnafrag_end >= $right_to->dnafrag_start) or
          ($left_to->dnafrag_strand < 0 and $left_to->dnafrag_start <= $right_to->dnafrag_end)) {
        
        throw("Cannot use direct target coordinates with inconsistent set of blocks");
      }
    }

    $tname = $right_to->dnafrag->name if not defined $tname;
    $tstrand = $right_to->dnafrag_strand if not defined $tstrand;
  }

  my $target_slice;
  if ($tstrand < 0) {
    $target_slice = $slices->{$tname}->invert;
  } else {
    $target_slice = $slices->{$tname};
  }

  return $target_slice;
}


##############################################
# Get/Sets
##############################################

sub alignment_mapper {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_alignment_mapper} = $val;
  }

  return $self->{_alignment_mapper};
}


sub direct_target_slice {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_direct_target} = $val;
  }

  return $self->{_direct_target};
}

sub max_readthrough_dist {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_max_readthrough} = $val;
  }

  return $self->{_max_readthrough};
}


sub from_slice {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_from_slice} = $val;
  }

  return $self->{_from_slice};
}


sub from_mapper {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_from_mapper} = $val;
  }

  return $self->{_from_mapper};
}


sub to_slices {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_to_slices} = $val;
  }

  return $self->{_to_slices};
}


sub to_mapper {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_to_mapper} = $val;
  }
  return $self->{_to_mapper};

}



1;
