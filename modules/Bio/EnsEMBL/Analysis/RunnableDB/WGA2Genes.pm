# Cared for by Ensembl
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::WGA2Genes

=head1 SYNOPSIS

  my $db      = Bio::EnsEMBL::DBAdaptor->new($locator);
  my $genscan = Bio::EnsEMBL::Analysis::RunnableDB::WGA2Genes->new (
                                                    -db         => $pipelinedb,
                                                    -input_id   => $input_id
                                                    -analysis   => $analysis );
  $genscan->fetch_input();
  $genscan->run();
  $genscan->write_output(); #writes to DB


=head1 DESCRIPTION


=cut
package Bio::EnsEMBL::Analysis::RunnableDB::WGA2Genes;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Config::General;

use Bio::EnsEMBL::Analysis::Config::WGA2Genes;

use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use Bio::EnsEMBL::Pipeline::Runnable::Genewise;
use Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript;
use Bio::EnsEMBL::Analysis::Runnable::AlignmentNets;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;


@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB);

############################################################
sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  $self->read_and_check_config($WGA2GENES_CONFIG_BY_LOGIC);

  return $self;
}




=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
  my( $self) = @_; 
  
  my $input_id = $self->input_id;  
  throw("No input id") unless defined($input_id);

  my $q_dbh = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%{$self->QUERY_CORE_DB});
  my $t_dbh = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%{$self->TARGET_CORE_DB});
  my $compara_dbh = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(%{$self->COMPARA_DB});
  
  my $query_species = $q_dbh->get_MetaContainerAdaptor->get_Species->binomial;
  my $target_species = $t_dbh->get_MetaContainerAdaptor->get_Species->binomial;
  
  my $gdb_adap = $compara_dbh->get_GenomeDBAdaptor;
  my $q_gdb = $gdb_adap->fetch_by_name_assembly($query_species);
  my $t_gdb = $gdb_adap->fetch_by_name_assembly($target_species);

  ################################################
  # check that the default assembly for the query and target agrees with that
  # for the method_link_species_set GenomeDBs
  #################################################

  my ($q_assembly_version, $t_assembly_version);
  eval {
    $q_assembly_version = $q_dbh->get_CoordSystemAdaptor->fetch_by_name
        ('toplevel',
         $q_gdb->assembly);

    $t_assembly_version = $t_dbh->get_CoordSystemAdaptor->fetch_by_name
        ('toplevel',
         $t_gdb->assembly);
  };
  $@ and do {
    throw("Had trouble fetching coord systems for ". 
          $q_gdb->assembly . " and " . 
          $t_gdb->assembly . " from core dbs: $@");
  };

  #####
  # fetch the gene; need to work in the coordinate space of the top-level slice to 
  # be consistent with compara
  #####
  
  my ($gene);
  eval {
    $gene = $q_dbh->get_GeneAdaptor->fetch_by_stable_id($input_id);
  };
  if ($@ or not defined $gene) {
    throw("Could not find gene '$input_id' in query database");
  }
  my $tl_q_slice = $q_dbh->get_SliceAdaptor->fetch_by_region('toplevel', 
                                                                    $gene->seq_region_name);
  $gene = $gene->transfer($tl_q_slice);

  ################################################################
  # get the compara data: MethodLinkSpeciesSet, reference DnaFrag, 
  # and all GenomicAlignBlocks
  ################################################################
  my $mlss = $compara_dbh->get_MethodLinkSpeciesSetAdaptor
      ->fetch_by_method_link_type_GenomeDBs($self->INPUT_METHOD_LINK_TYPE,
                                            [$q_gdb, $t_gdb]);
  throw("No MethodLinkSpeciesSet for :\n" .
        $self->INPUT_METHOD_LINK_TYPE . "\n" . 
        $query_species . "\n" . 
        $target_species)
      if not $mlss;

  my $dnafrag = $compara_dbh->
      get_DnaFragAdaptor->fetch_by_GenomeDB_and_name($q_gdb,$tl_q_slice->seq_region_name);

  my @gen_al_blocks = @{$compara_dbh->get_GenomicAlignBlockAdaptor->
                            fetch_all_by_MethodLinkSpeciesSet_DnaFrag($mlss,
                                                                      $dnafrag,
                                                                      $gene->start,
                                                                      $gene->end)};

  my (%chains, @chains);
  foreach my $block (@gen_al_blocks) {
    my $qga = $block->reference_genomic_align;    
    my ($tga) = @{$block->get_all_non_reference_genomic_aligns};
    
    # fetch the target slice for later reference
    if (not exists $self->target_slices->{$tga->dnafrag->name}) {

      $self->target_slices->{$tga->dnafrag->name} = 
          $t_dbh->get_SliceAdaptor->fetch_by_region('toplevel',
                                                    $tga->dnafrag->name);
    }

    my $chain_id =  $qga->genomic_align_group_id_by_type("chain");

    push @{$chains{$chain_id}}, $block;
  }
  foreach my $chain_id (keys %chains) {
    push @chains, [
                   sort {
                     $a->reference_genomic_align->dnafrag_start <=> 
                         $b->reference_genomic_align->dnafrag_start;
                   } @{$chains{$chain_id}}
                   ];
  }

  $self->query_top_level_slice($tl_q_slice);
  $self->gene($gene);      
  $self->genomic_align_block_chains(\@chains);
}


=head2 run

    Title   :   run
    Usage   :   $self->run
    Returns :   nothing
    Args    :   none

=cut


sub run {
  my ($self) = @_;

  my @blocks_used;
  my @alignment_chains = @{$self->genomic_align_block_chains};

  my @cds_feats = @{$self->get_all_transcript_cds_features(@{$self->get_all_Transcripts($self->gene)})};


  printf("%s \#SUMMARY (%s/%d-%d %d)\n",
         $self->gene->stable_id,
         $self->gene->seq_region_name,
         $self->gene->start,
         $self->gene->end, 
         $self->gene->strand);

  for(my $iteration=0; ;$iteration++) {
    my $net_blocks = $self->make_alignment_net(\@alignment_chains, \@blocks_used);
    my (@projected_cds_feats);
    
    my $blocks_used_before = scalar(@blocks_used);

    foreach my $exon (@cds_feats) {
      my ($proj_ex_els, $blocks_used) = $self->map_feature_to_target($exon, $net_blocks);
      
      push @projected_cds_feats, $proj_ex_els;
      push @blocks_used, @$blocks_used;
    }  

    my $additional_blocks_used = scalar(@blocks_used) - $blocks_used_before;

    if ($additional_blocks_used == 0) {
      # no alignment coverage, so stop;
      last;
    }

    my ($gene_scaffold, 
        $query_to_genescaf_map,
        $target_to_genescaf_map) = 
            $self->gene_scaffold_from_projection($self->gene->stable_id . "-" . $iteration, 
                                                 \@projected_cds_feats);

    if (not defined $gene_scaffold) {
      # inconsistent assembly; be conservative and dont process further down for now
      last;
    }
    
    my @realigned_trans;
    foreach my $tran (@{$self->get_all_Transcripts($self->gene)}) {
      my @trans_frags;
      #@trans_frags = $self->align_protein($gene_scaffold,
      #                                    $self->get_translation($tran),
      #                                    $tran->strand);  
      my $proj_trans = $self->make_projected_transcript($tran,
                                                        $gene_scaffold,
                                                        $query_to_genescaf_map);
                                                        

      if (defined $proj_trans) {
        push @realigned_trans,  {
          original => $tran,
          #projected => $self->process_transcripts($gene_scaffold, [$proj_trans]),
          projected => $proj_trans,
          realigned => $self->process_transcripts($gene_scaffold, \@trans_frags),
        };
      }
    }
    
    $self->write_report($iteration,
                        $self->gene,
                        $gene_scaffold,
                        $target_to_genescaf_map,
                        \@projected_cds_feats,
                        \@realigned_trans);

    # record which alignment blocks were used in this iteration for the next
    # for now though, just do one iteration
    #last;
  }
}

sub write_output {
  my ($self) = @_;

  # do nothing

  return;
}


######################################
# internal methods
#####################################

sub make_alignment_net_old {
  my ($self, $raw_chains, $excluded_regions) = @_;

  my @blocks;
  foreach my $chain (@{$raw_chains}) {
    foreach my $block (@{$chain}) {
      push @blocks, $block;
    }
  }

  @blocks = sort { 
    $a->reference_genomic_align->dnafrag_start <=> $b->reference_genomic_align->dnafrag_start;
  } @blocks;

  return \@blocks;
}



=head2 make_alignment_net

 Title   : make_alignment_net
 Description:
    Takes a list of chains of GenomicAlignBlocks and computes a query-oriented
    "Net" so that no bp on the query sequence is covered by more than one
    target sequence. 

=cut

sub make_alignment_net {
  my ($self, $raw_chains, $excluded_regions) = @_; 

  my @chains_to_process;

  foreach my $chain (@$raw_chains) {
    # if this chain conatains a block that overlaps with any of the 
    # excluded regions on the target, ignore it.
    my $keep_chain = 1;

    CHAIN: foreach my $block (@$chain) {
      my ($tga) = @{$block->get_all_non_reference_genomic_aligns};
      foreach my $ex_block (@$excluded_regions) {
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

  # make hashes for lengths of query and target
  my (%q_len, %t_lens);
  $q_len{$self->query_top_level_slice->seq_region_name} = 
      $self->query_top_level_slice->length;

  foreach my $k (keys %{$self->target_slices}) {
    my $ts = $self->target_slices->{$k};
    $t_lens{$ts->seq_region_name} = $ts->length;
  }

  my $run = Bio::EnsEMBL::Analysis::Runnable::AlignmentNets->
      new(-analysis       => Bio::EnsEMBL::Analysis->new(),
          -chains         => \@chains_to_process,
          -query_lengths  => \%q_len,
          -target_lengths => \%t_lens,
          -chainnet       => '/usr/local/ensembl/bin/chainNet');
  $run->run;

  my @blocks;
  foreach my $chain (@{$run->output}) {
    foreach my $block (@$chain) {
      push @blocks, $block;
    }
  }

  @blocks = sort { 
    $a->reference_genomic_align->dnafrag_start <=> $b->reference_genomic_align->dnafrag_start;
  } @blocks;
  
  return \@blocks;
}



=head2 align_protein

 Title   : align_protein
 Description:
    Aligns the given protein to the given genomic slice 
    (using either Exonerate or Genewise depending on option)
                                                        
=cut

sub align_protein {
  my ($self, $genomic, $protein, $strand) = @_;

  my $runnable;

  my $analysis = Bio::EnsEMBL::Analysis->new();

  if ($self->USE_GENEWISE) {
    $runnable = Bio::EnsEMBL::Pipeline::Runnable::Genewise->
      new(
          -analysis  => $analysis,
          -query     => $genomic,
          -protein   => $protein,
          -gap       => 12,
          -extension => 2,
          -memory    => 100000,
          -subs      => 0.0000001,
          -matrix    => 'BLOSUM62.bla',
          -program   => 'genewise-2.2.3-optlin623',
          -options   => '-quiet ',
          -reverse   => $strand < 0 ? 1 : 0,
          -verbose   => 1,
          );

  } else {
    my $EXONERATE = "/usr/local/ensembl/bin/exonerate-0.9.0";
    my $EXONERATE_OPTIONS = "--model protein2genome --softmasktarget TRUE ";
    #my $EXONERATE_OPTIONS = "--model protein2genome --exhaustive TRUE --softmasktarget TRUE ";
    
    $runnable = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript->
        new(-analysis     => $analysis,
            -target_seqs    => [$genomic],
            -query_seqs     => [$protein],
            -query_type     => 'protein',
            -coverage_by_aligned => 0,
            -exonerate      => $EXONERATE,
            -options        => $EXONERATE_OPTIONS,
            -verbose        => 1,
            );
  }    
  
  $runnable->run;
  return @{$runnable->output};
}




=head2 map_feature_to_target

 Title   : map_feature_to_target
 Decription:
    Takes an feature and a list of GenomicAlignBlocks and returns a list of 
    {query => [], target = []}, representing the gapped alignment of the query 
    (in the region of the feature) to the target. The alignment is represented 
    a list of Bio::EnsEMBL::Mapper::Coordinates and Bio::EnsEMBL::Mapper::Gap
    for each of the query and target

=cut

sub map_feature_to_target {
  my ($self, $feat, $gen_al_blocks) = @_;

  my @projected_feat_segments;
  my @overlapping_blocks;

  # binary search the block list
  #my ($left, $j) = (0, scalar(@$gen_al_blocks)-1);
  #while ($j - $left > 0) {
  #  my $mid = int(($left + $j) / 2);
  #  if ($gen_al_blocks->[$mid]->reference_genomic_align->dnafrag_end < $feat->start) {
  #    $left = $mid + 1;
  #  } else {
  #    $j = $mid;
  #  }
  #}
  my $j=0;

  for(my $i=$j; $i < @$gen_al_blocks; $i++) {
    my $block = $gen_al_blocks->[$i];

    my $qga = $block->reference_genomic_align;
    my ($tga) = @{$block->get_all_non_reference_genomic_aligns};

    next if $qga->dnafrag_end < $feat->start;
    last if $qga->dnafrag_start > $feat->end;

    my $q_feat_start = $qga->dnafrag_start < $feat->start ? $feat->start : $qga->dnafrag_start;
    my $q_feat_end   = $qga->dnafrag_end   > $feat->end   ? $feat->end   : $qga->dnafrag_end;
    
    push @overlapping_blocks, $block;

    my $current_pos = $q_feat_start;
    foreach my $coord ($qga->get_Mapper->map_coordinates("sequence",
                                                         $q_feat_start,
                                                         $q_feat_end,
                                                         $qga->dnafrag_strand,
                                                         "sequence")) {
      
      if ($coord->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
        my $current_pos2 = $current_pos;
        
        foreach my $coord2 ($tga->get_Mapper->map_coordinates("alignment",
                                                              $coord->start,
                                                              $coord->end,
                                                              1,
                                                              "alignment")) {
          
          if ($coord2->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
            push @projected_feat_segments, { 
              query  => Bio::EnsEMBL::Mapper::Coordinate->new($feat->slice->seq_region_name,
                                                              $current_pos2,
                                                              $current_pos2 + ($coord2->end - $coord2->start),
                                                              1),
              target => Bio::EnsEMBL::Mapper::Coordinate->new($tga->dnafrag->name,
                                                              $coord2->start,
                                                              $coord2->end,
                                                              $qga->dnafrag_strand * $tga->dnafrag_strand),
            };
          }
          
          $current_pos2 += ($coord2->end - $coord2->start + 1);
        }          
      }
      
      $current_pos += ($coord->end - $coord->start + 1);
    }
  }
  
  #########
  # fill in Gap segments between the coordinates so that each original feat
  # is completely accounted for by scaffold segments and gap segments
  #########
  
  if (not @projected_feat_segments) {
    
    push @projected_feat_segments, {
      query  => Bio::EnsEMBL::Mapper::Coordinate->new($feat->slice->seq_region_name,
                                                      $feat->start,
                                                      $feat->end,
                                                      $feat->slice->strand),
      target => Bio::EnsEMBL::Mapper::Gap->new($feat->start, $feat->end),
    };
  } 
  else {
    # start
    if ($feat->start < $projected_feat_segments[0]->{query}->start) {
      my $st = $feat->start;
      my $en = $projected_feat_segments[0]->{query}->start - 1;

      unshift @projected_feat_segments, {
        query  => Bio::EnsEMBL::Mapper::Coordinate->new($feat->slice->seq_region_name,
                                                        $st,
                                                        $en,
                                                        $feat->slice->strand),
        target => Bio::EnsEMBL::Mapper::Gap->new($st, $en),
      };
    }
    
    # end
    if ($feat->end > $projected_feat_segments[-1]->{query}->end) {
      my $st = $projected_feat_segments[-1]->{query}->end + 1;
      my $en = $feat->end;
      
      push @projected_feat_segments, {
        query  => Bio::EnsEMBL::Mapper::Coordinate->new($feat->slice->seq_region_name,
                                                        $st,
                                                        $en,
                                                        $feat->slice->strand),
        target => Bio::EnsEMBL::Mapper::Gap->new($st, $en),
      };
    }
    
    # middle
    my @new_bits;
    for(my $i=1; $i < @projected_feat_segments; $i++) {
      my ($l, $r) = ($projected_feat_segments[$i-1], $projected_feat_segments[$i]);
      if ($r->{query}->start > $l->{query}->end + 1) {
        my $st = $l->{query}->end + 1;
        my $en = $r->{query}->start - 1;
        
        push @new_bits, {
          query  => Bio::EnsEMBL::Mapper::Coordinate->new($feat->slice->seq_region_name,
                                                          $st,
                                                          $en,
                                                          $feat->slice->strand),
          target => Bio::EnsEMBL::Mapper::Gap->new($st, $en),
        };
      }
    }
    
    push @projected_feat_segments, @new_bits;
  }

  @projected_feat_segments = 
      sort {$a->{query}->start <=> $b->{query}->start} @projected_feat_segments;
  
  return (\@projected_feat_segments, \@overlapping_blocks);
}




=head2 gene_scaffold_from_projection

 Title   : gene_scaffold_from_projection
 Description:
    Takes a list of exon projection elements (representing the gapped
    alignment between the query and target in the region of the exons)
    and returns 3 things:
    1. A "gene scaffold" which is a mini-assembly of the target sequences
       necessary to present the whole gene on a single sequence
    2. Query mapper: a map between the coords of the gene scaffold and the
       original query sequence
    3. Target mapper: a map between the coords of the gene scaffold and the
       comprised target sequences.                                                        
=cut

sub gene_scaffold_from_projection {
  my ($self, $header, $projected_exon_elements) = @_;

  my @projected_exon_elements = @$projected_exon_elements;

  my (@target_coords, @new_targets);

  # step 0: flatten list
  my $ex_num = $self->gene->strand > 0 ? 1 : scalar(@projected_exon_elements);
  foreach my $ex (@projected_exon_elements) {
    foreach my $coord_pair (@$ex) {
      my $coord = $coord_pair->{target};

      push @target_coords, {
        exon  => $ex_num,
        exons => { $ex_num => 1 },
        coord => $coord,
      };
    }
    $ex_num += $self->gene->strand;
  }

  ############ DEBUG
  #if (0) {
    foreach my $fragp (@target_coords) {
      my $frag = $fragp->{coord};
      if ($frag->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
        printf "$header BEFORE SEG: %s %d %d %s ", $frag->id, $frag->start, $frag->end, $frag->strand;
      } else {
        printf "$header BEFORE GAP: %d %d (%d) ", $frag->start, $frag->end, $frag->end - $frag->start + 1; 
      }
      print "[";
      foreach my $k (sort {$a <=> $b} keys %{$fragp->{exons}}) {
        print " $k ";
      }
      print "]\n";
    }
  #}
  

  # step 1 remove gaps which cannot be filled
  @new_targets = ();
  for(my $i=0; $i < @target_coords; $i++) {
    my $this_coord = $target_coords[$i];

    if ($this_coord->{coord}->isa("Bio::EnsEMBL::Mapper::Gap")) {
      # if it's gap that can be filled, leave it. Otherwise, remove it
      my ($left_coord, $right_coord);
      for(my $j=$i-1; $j >=0; $j--) {
        if ($target_coords[$j]->{coord}->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
          $left_coord = $target_coords[$j];
          last;
        }
      }
      for(my $j=$i+1; $j < @target_coords; $j++) {
        if ($target_coords[$j]->{coord}->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
          $right_coord = $target_coords[$j];
          last;
        }
      }

      # under what circumstances can this gap NOT be filled?
      # 1. If one of the flanking coords is on the same exon
      #    as the gap, and the end of the aligned region does
      #    not align to a gap
      # 2. If the 2 flanking coords are consistent and not
      #    separated by a sequence-level gap
      my $keep_gap = 1;

      if (defined $left_coord and 
          defined $right_coord and
          $self->check_consistent_coords($left_coord->{coord}, 
                                         $right_coord->{coord})) {
        
        my ($left, $right) = ($left_coord->{coord}, 
                              $right_coord->{coord});
        if ($left->strand < 0) {
          ($left, $right) = ($right, $left);
        }
        my ($new_left, $new_right) = 
            $self->separate_coords($left, $right,
                                   $self->target_slices->{$left_coord->{coord}->id});
        
        if (not defined $new_left and not defined $new_right) {
          $keep_gap = 0;
        }
      } 


      if ($keep_gap and 
          defined $left_coord and 
          $left_coord->{exon} == $this_coord->{exon}) {
        
        my ($new_left, $dummy);
        if ($left_coord->{coord}->strand < 0) {
          ($dummy, $new_left) = 
              $self->separate_coords(undef,
                                     $left_coord->{coord}, 
                                     $self->target_slices->{$left_coord->{coord}->id});

          if (abs($new_left->start - $left_coord->{coord}->start) > 10) {
            # we're not near the end of a contig. This gap cannot be filled;
            $keep_gap = 0;
          }
        } else {
          ($new_left, $dummy) = 
              $self->separate_coords($left_coord->{coord}, 
                                     undef,
                                     $self->target_slices->{$left_coord->{coord}->id});

          if (abs($new_left->end - $left_coord->{coord}->end) > 10) {
            # we're not near the end of a contig. This gap cannot be filled;
            $keep_gap = 0;
          }
        }
      }

      if ($keep_gap and 
          defined $right_coord and
          $right_coord->{exon} == $this_coord->{exon}) {

        my ($dummy, $new_right);
        if ($right_coord->{coord}->{strand} < 0) {
          ($new_right, $dummy) = 
              $self->separate_coords($right_coord->{coord}, 
                                     undef,
                                     $self->target_slices->{$right_coord->{coord}->id});

          if (abs($right_coord->{coord}->end - $new_right->end) > 10) {
            # we're not near the end of a contig. This gap cannot be filled;
            $keep_gap = 0;
          }
        } else {
          ($dummy, $new_right) = 
              $self->separate_coords(undef,
                                     $right_coord->{coord}, 
                                     $self->target_slices->{$right_coord->{coord}->id});

          if (abs($right_coord->{coord}->start - $new_right->start) > 10) {
            # we're not near the end of a contig. This gap cannot be filled;
            $keep_gap = 0;
          }
        }
      }
      
      if ($keep_gap) {
        push @new_targets, $target_coords[$i];
      }

    } else {
      push @new_targets, $target_coords[$i];
    }
  }
  @target_coords = @new_targets;


  ############ DEBUG
  if (0) {
    print "\nAFTER REMOVING UNPLUGGABE GAPS\n";
    foreach my $fragp (@target_coords) {
      my $frag = $fragp->{coord};
      if ($frag->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
        printf "SEG: %s %d %d %s ", $frag->id, $frag->start, $frag->end, $frag->strand;
      } else {
        printf "GAP: %d %d (%d) ", $frag->start, $frag->end, $frag->end - $frag->start + 1; 
      }
      print "[";
      foreach my $k (sort {$a <=> $b} keys %{$fragp->{exons}}) {
        print " $k";
      }
      print "]\n";
    }
  }    
  ##########################

  # step 2: merge adjacent targets   
  @new_targets = ();
  foreach my $target_coord (@target_coords) {

    if ($target_coord->{coord}->isa("Bio::EnsEMBL::Mapper::Coordinate") and
        @new_targets and
        $new_targets[-1]->{coord}->isa("Bio::EnsEMBL::Mapper::Coordinate") and
        $self->check_consistent_coords($new_targets[-1]->{coord}, 
                                       $target_coord->{coord})) {
      
      my $last_coord = pop @new_targets;
      my $new_coord = $self->merge_coords($last_coord->{coord},
                                          $target_coord->{coord});
      push @new_targets,  {
        coord => $new_coord,
        exons => {},
      };
      map { $new_targets[-1]->{exons}->{$_} = 1 } (keys %{$last_coord->{exons}}, 
                                                   keys %{$target_coord->{exons}});
      
    } else {
      push @new_targets, $target_coord;
    }
  }    
  @target_coords = @new_targets;


  ############ DEBUG
  if (0) {
    print "\nAFTER MERGING ADJACENT THINGS\n";
    foreach my $fragp (@target_coords) {
      my $frag = $fragp->{coord};
      if ($frag->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
        printf "SEG: %s %d %d %s ", $frag->id, $frag->start, $frag->end, $frag->strand;
      } else {
        printf "GAP: %d %d (%d) ", $frag->start, $frag->end, $frag->end - $frag->start + 1; 
      }
      print "[";
      foreach my $k (sort {$a <=> $b} keys %{$fragp->{exons}}) {
        print " $k";
      }
      print "]\n";
    }
  }
  ##########################

  
  
  # step 3: extend the pieces as far as possible
  @new_targets = ();
  for(my $i=0; $i < @target_coords; $i++) {
    my $this_coord = $target_coords[$i];

    if ($this_coord->{coord}->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
      # scan left
      my $left_consistent;
      for(my $j = $i-1; $j >= 0; $j--) {
        my $that_coord = $target_coords[$j];
        if ($that_coord->{coord}->isa("Bio::EnsEMBL::Mapper::Coordinate") and
            $self->check_consistent_coords($that_coord->{coord}, 
                                           $this_coord->{coord})) {
          $left_consistent = $that_coord;
          last;
        }
      }
      # scan right
      my $right_consistent;
      for(my $j = $i+1; $j < @target_coords; $j++) {
        my $that_coord = $target_coords[$j];
        if ($that_coord->{coord}->isa("Bio::EnsEMBL::Mapper::Coordinate") and
            $self->check_consistent_coords($this_coord->{coord}, 
                                           $that_coord->{coord})) {
          $right_consistent = $that_coord;
          last;
        }
      }
      
      my $extended_left = 1;
      my $extended_right = $self->target_slices->{$this_coord->{coord}->{id}}->length;
      if (defined $left_consistent) {
        if ($this_coord->{coord}->strand < 0) {
          # left_coord will be downtream of this_coord in target
          my ($new_right, $new_left) = 
              $self->separate_coords($this_coord->{coord},
                                     $left_consistent->{coord},
                                     $self->target_slices->{$this_coord->{coord}->id});
        
          if (defined $new_right) {
            $extended_right = $new_right->end;
          }
        } else {
          # left_coord will be uptream of this_coord in target
          my ($new_left, $new_right) = 
              $self->separate_coords($left_consistent->{coord},
                                     $this_coord->{coord},
                                     $self->target_slices->{$this_coord->{coord}->id});
        
          if (defined $new_right) {
            $extended_left = $new_right->start;
          }
        }
      }
      if (defined $right_consistent) {
        if ($this_coord->{coord}->strand < 0) {
          # right_coord will be uptream of this_coord in target
          my ($new_right, $new_left) = 
              $self->separate_coords($right_consistent->{coord},
                                     $this_coord->{coord},
                                     $self->target_slices->{$this_coord->{coord}->id});
          if (defined $new_left) {
            $extended_left = $new_left->start;
          }
        } else {
          # right_coord will be downtream of this_coord in target
          my ($new_left, $new_right) = 
              $self->separate_coords($this_coord->{coord},
                                     $right_consistent->{coord},
                                     $self->target_slices->{$this_coord->{coord}->id});
          if (defined $new_left) {
            $extended_right = $new_left->end;
          }          
        }      
      }

      push @new_targets, {
        coord => Bio::EnsEMBL::Mapper::Coordinate->new($this_coord->{coord}->id,
                                                       $extended_left, 
                                                       $extended_right,
                                                       $this_coord->{coord}->strand),
        exons => $this_coord->{exons},
      };
    }
    else {
      push @new_targets, $this_coord;
    }
  }
  @target_coords = @new_targets;

  ############ DEBUG
  #if (0) {
    foreach my $fragp (@target_coords) {
      my $frag = $fragp->{coord};
      if ($frag->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
        printf "$header AFTER SEG: %s %d %d %s ", $frag->id, $frag->start, $frag->end, $frag->strand;
      } else {
        printf "$header AFTER GAP: %d %d (%d) ", $frag->start, $frag->end, $frag->end - $frag->start + 1; 
      }
      print "[";
      foreach my $k (sort {$a <=> $b} keys %{$fragp->{exons}}) {
        print " $k";
      }
      print "]\n";
    }
  #}
  ##########################

  # check for inconsistencies; if the overlapping scaffold fragments remain, 
  # we have an inconsistent mini-assembly (usually caused by a scaffold that 
  # could not be sensibly split
  FRAGS: for(my $i=0; $i<@target_coords; $i++) {
    my $coord_i = $target_coords[$i]->{coord}; 
    next if not $coord_i->isa("Bio::EnsEMBL::Mapper::Coordinate");
    for(my $j=$i-1; $j>=0; $j--) {
      my $coord_j = $target_coords[$j]->{coord}; 
      next if not $coord_j->isa("Bio::EnsEMBL::Mapper::Coordinate");
      if ($coord_i->id eq $coord_j->id and
          ($coord_i->strand ne $coord_j->strand or            
           ($coord_i->end >= $coord_j->start and
            $coord_i->start <= $coord_j->end))) {
        print "$header AFTER\tINCONSISTENT\n";
        #last FRAGS;
        return ();
      }           
    }
  }
  print "\n";
  #############


  # Reseparate according to exons. This will not result in one element 
  # per exon because some assembly frgaments (hopefully many) will have more
  # than one exon
  my (@fragments, %exons_seen);
  foreach my $coord (@target_coords) {
    my $common = 0;
    foreach my $enum (keys %{$coord->{exons}}) {
      if (exists($exons_seen{$enum})) {
        $common = 1;
      }
      $exons_seen{$enum} = 1;
    }
    if ($common) {
      push @{$fragments[-1]}, $coord;
    } else {
      push @fragments, [$coord];
    }
  }


  ###################################
  # build the gene scaffold, and mappings between the new gene
  # scaffold and each of query and target coords

  my $target_map = Bio::EnsEMBL::Mapper->new('target',
                                             'genescaffold'); 
  my $query_map = Bio::EnsEMBL::Mapper->new('query',
                                            'genescaffold');
  my $gene_scaffold_id = "the_gene_scaffold";
  
  my ($seq, $last_end_pos) = ("", 0);
  for(my $i=0; $i < @fragments; $i++) {
    my $frag = $fragments[$i];

    foreach my $coord_exons (@$frag) {
      my $coord = $coord_exons->{coord};

      if ($coord->isa("Bio::EnsEMBL::Mapper::Coordinate")) {            
        # the sequence itself        
        my $slice = $self->target_slices->{$coord->id};
        my $this_seq = $slice->get_repeatmasked_seq(['RepeatMask'],1)->seq;        
        $this_seq = substr($this_seq,
                           $coord->start  - 1,
                           $coord->length);        
        if ($coord->strand < 0) {
          reverse_comp(\$this_seq);
        }      
        $seq .= $this_seq;
        
        # and the map
        $target_map->add_map_coordinates($coord->id,
                                         $coord->start,
                                         $coord->end,
                                         $coord->strand,
                                         $gene_scaffold_id,
                                         $last_end_pos + 1,
                                         $last_end_pos + $coord->length);
      } else {
        # seq sequence itself
        $seq .= ('n' x $coord->length);

        # and the map. This is a target gap we have "filled", so no position 
        # in target, but a position in query
        $query_map->add_map_coordinates($self->gene->slice->seq_region_name,
                                        $coord->start,
                                        $coord->end,
                                        1,
                                        $gene_scaffold_id,
                                        $last_end_pos + 1,
                                        $last_end_pos + $coord->length);
      }
      
      $last_end_pos += $coord->length;
    }

    if ($i < @fragments - 1) {
      $seq .= ('n' x $self->GENE_SCAFFOLD_COMPONENT_PADDING);
      $last_end_pos += $self->GENE_SCAFFOLD_COMPONENT_PADDING;
    }
  }

  # now add of the original exon pieces to the query map
  foreach my $el (@projected_exon_elements) {
    foreach my $pair (@$el) {
      my ($q, $t) = ($pair->{query}, $pair->{target});

      if ($t->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
        # get the gene_scaffold position from the target map
        my ($coord) = $target_map->map_coordinates($t->id,
                                                   $t->start,
                                                   $t->end,
                                                   $t->strand,
                                                   'target');
        
        $query_map->add_map_coordinates($self->gene->slice->seq_region_name,
                                        $q->start,
                                        $q->end,
                                        1,
                                        $gene_scaffold_id,
                                        $coord->start,
                                        $coord->end);
      }
    }
  }

  # finally, make the gene scaffold itself
  my $gene_scaffold = Bio::EnsEMBL::Slice
      ->new(-coord_system => Bio::EnsEMBL::CoordSystem->new(-name => 'genescaffold',
                                                            -rank => 1),
            -seq_region_name => $gene_scaffold_id,
            -seq  => $seq,
            -start => 1,
            -end   => length($seq));

  return ($gene_scaffold, $query_map, $target_map);
}




=head2 make_projected_transcript

 Title   : make_projected_transcript
 Description:
    Takes a transcript, a gene scaffold, and a mapping between query coords and
    gene scaffold coords. Produces the transcript that is the result of "projecting"
    the original transcript, through alignment, onto the gene scaffold. 
=cut

sub make_projected_transcript {
  my ($self, $transcript, $gene_scaf, $mapper) = @_;

  my (@all_coords, @new_exons);

  my @orig_exons = @{$transcript->get_all_translateable_Exons};
  if ($transcript->strand < 0) {
    @orig_exons = reverse @orig_exons;
  }

  foreach my $orig_exon (@orig_exons) {    
    my @these_coords = $mapper->map_coordinates($orig_exon->slice->seq_region_name,
                                                $orig_exon->start,
                                                $orig_exon->end,
                                                1,
                                                "query");
    push @all_coords, @these_coords;
  }


  my $need_another_pass;
  do {
    $need_another_pass = 0;

    my (@processed_coords, @gap_indices);
    # merge gaps
    foreach my $c (@all_coords) {
      if ($c->isa("Bio::EnsEMBL::Mapper::Gap")) {
        if (@processed_coords and $processed_coords[-1]->isa("Bio::EnsEMBL::Mapper::Gap")) {
          $processed_coords[-1]->end( $processed_coords[-1]->end + $c->length );
        } else {
          push @processed_coords, Bio::EnsEMBL::Mapper::Gap->new(1, $c->length);
          push @gap_indices, scalar(@processed_coords) - 1;
        }
      } else {
        push @processed_coords, $c;
      }
    }

    GAP: foreach my $idx (@gap_indices) {
      my $gap = $processed_coords[$idx];
      my $frameshift = $gap->length % 3;

      # if gap is at the start, look right, otherwise look left
      if ($frameshift) {
        my $bases_to_remove = 3 - $frameshift;      

        # calculate "surplus" bases on incomplete codes to left and right
        my ($left_surplus, $right_surplus) = (0,0);
        for(my $j=$idx-1; $j >= 0; $j--) {
          $left_surplus += $processed_coords[$j]->length;
        }
        for(my $j=$idx+1; $j < @processed_coords; $j++) {
          $right_surplus += $processed_coords[$j]->length;
        }
        
        $left_surplus  = $left_surplus % 3;
        $right_surplus = $right_surplus % 3;

        if ($left_surplus) {
          # eat left
          $bases_to_remove = $left_surplus;
          
          my $left_coord = $processed_coords[$idx - 1];
          if ($left_coord->length > $bases_to_remove) {
            $gap->end($gap->end + $bases_to_remove);
            $left_coord->end( $left_coord->end - $bases_to_remove );
          } else {
            # we need to eat away the whole of this coord
            $processed_coords[$idx-1] = Bio::EnsEMBL::Mapper::Gap->new(1,$left_coord->length);
          }
          $need_another_pass = 1;
          last GAP;
        } else {
          # eat right; $bases_to_remove must all be in $right_surplus
          my $right_coord = $processed_coords[$idx + 1];
          if ($right_coord->length > $bases_to_remove) {
            $gap->end($gap->end + $bases_to_remove);
            $right_coord->start( $right_coord->start + $bases_to_remove);
          } else {
            # we need to eat away the whole of this coord
            $processed_coords[$idx+1] = Bio::EnsEMBL::Mapper::Gap->new(1,$right_coord->length);
          }
        }

        $need_another_pass = 1;
        last GAP;
      }      
    }
    @all_coords = @processed_coords;    
  } while ($need_another_pass);


  foreach my $coord (@all_coords) {
    if ($coord->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
      push @new_exons, Bio::EnsEMBL::Exon->new(-start     => $coord->start,
                                               -end       => $coord->end,
                                               -strand    => $transcript->strand,
                                               -slice     => $gene_scaf);
    }
  }

  if (not @new_exons) {
    return undef;
  }

  # sort exons into rank order 
  if ($transcript->strand < 0) {
    @new_exons = sort { $b->start <=> $a->start } @new_exons;
  } else {
    @new_exons = sort { $a->start <=> $b->start } @new_exons;
  }

  # calculate phases and set translation
  my ($previous_exon);
  foreach my $exon (@new_exons) {

    if (defined $previous_exon) {
      $exon->phase($previous_exon->end_phase);
    } else {
      $exon->phase(0);
    }

    $exon->end_phase( (($exon->end - $exon->start + 1) + $exon->phase) % 3 );



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

      my @gen_coords = $mapper->map_coordinates($gene_scaf->seq_region_name,
                                                $extent_start,
                                                $extent_end,
                                                1,
                                                "genescaffold");       

      my @fps;
      my $cur_gs_start = $extent_start;
      foreach my $g_coord (@gen_coords) {
        my $cur_gs_end = $cur_gs_start + $g_coord->length - 1;
                
        if ($g_coord->isa("Bio::EnsEMBL::Mapper::Coordinate")) {

          my ($p_coord) = $transcript->genomic2pep($g_coord->start, 
                                                   $g_coord->end,
                                                   $exon->strand);
          
          my $fp = Bio::EnsEMBL::FeaturePair->new(-seqname  => $gene_scaf->seq_region_name,
                                                  -start    => $cur_gs_start,
                                                  -end      => $cur_gs_end,
                                                  -strand   => $exon->strand,
                                                  -hseqname => $transcript->stable_id,
                                                  -hstart   => $p_coord->start,
                                                  -hend     => $p_coord->end,
                                                  -hstrand => $p_coord->strand);
          push @fps, $fp;
        }
        
        $cur_gs_start += $g_coord->length;
      }
        
      if (@fps) {
        $exon->add_supporting_features(Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@fps));
      }
    }

    $previous_exon = $exon;
  }


  # merge abutting exons; deals with exon fusion events, and 
  # small, frame-preserving insertions in the target
  my @merged_exons;
  foreach my $exon (@new_exons) {
    if (@merged_exons) {

      my $prev_exon = pop @merged_exons;
   
      my ($new_start, $new_end);

      if ($transcript->strand < 0) {
        my $intron_len = $prev_exon->start - $exon->end - 1;
        if ($intron_len % 3 == 0 and $intron_len <= 15) { 
          $new_start = $exon->start;
          $new_end   = $prev_exon->end;
        }
      } else {
        my $intron_len = $exon->start - $prev_exon->end - 1;
        if ($intron_len % 3 == 0 and $intron_len <= 15) {
          $new_start = $prev_exon->start;
          $new_end   = $exon->end;
        }
      }

      if (defined $new_start and defined $new_end) {
        my $merged_exon = Bio::EnsEMBL::Exon->new(-start => $new_start,
                                                  -end   => $new_end,
                                                  -strand => $transcript->strand,
                                                  -phase => $prev_exon->phase,
                                                  -end_phase => $exon->end_phase,
                                                  -slice  => $exon->slice);
        
        my @ug_feats;
        if (@{$prev_exon->get_all_supporting_features}) {
          push @ug_feats, $prev_exon->get_all_supporting_features->[0]->ungapped_features;
        }
        if (@{$exon->get_all_supporting_features}) {
          push @ug_feats, $exon->get_all_supporting_features->[0]->ungapped_features;
        }
        if (@ug_feats) {
          my $new_sup_feat = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@ug_feats);
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

  my ($hcoverage, @trans_fps);
  foreach my $exon (@merged_exons) {
    $proj_tran->add_Exon($exon);
    
    if (@{$exon->get_all_supporting_features}) {
      my @e_fps = $exon->get_all_supporting_features->[0]->ungapped_features;
      foreach my $fp (@e_fps) {
        $hcoverage += $fp->hend - $fp->hstart + 1;
        push @trans_fps, $fp;
      }
    }
  }
  
  my $t_sf = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@trans_fps);
  $t_sf->hcoverage(sprintf("%.1f", 100 * ($hcoverage / $transcript->translate->length)));
  $proj_tran->add_supporting_features($t_sf);

  my $translation = Bio::EnsEMBL::Translation->new();
  $translation->start_Exon($merged_exons[0]);
  $translation->start(1);
  $translation->end_Exon($merged_exons[-1]);
  $translation->end($merged_exons[-1]->end - $merged_exons[-1]->start + 1);

  $proj_tran->translation($translation);


  return $proj_tran;
}




=head2 process_transcripts

 Title   : process_transcripts
 Description:
    Filters the transcripts into a self-consistent set and modifies them 
    to skip over any in-frame stops
=cut

sub process_transcripts{
  my ($self,$gene_scaffold, $trans) = @_;

  foreach my $t (@$trans) {
    foreach my $e (@{$t->get_all_Exons}) {
      $e->slice($gene_scaffold);
    }
  }
  
  my @filtered_trans = @{$self->filter_transcripts($trans)};

  my @processed_transcripts;

  foreach my $tran (@filtered_trans) {
    my @exons = @{$tran->get_all_Exons};

    my $pep = $tran->translate->seq;

    my $had_stops = 0;
    while($pep =~ /\*/g) {
      my $position = pos($pep);
      my @coords = $tran->pep2genomic($position, $position);

      if (@coords > 1) {
        # the codon is split by an intron. Messy. Leave these for now
        next;
      } 

      my ($stop) = @coords;
 
      # locate the exon that this stop lies in
      my @new_exons;
      foreach my $exon (@exons) {

        if ($stop->start > $exon->start and $stop->end < $exon->end) {
          # this stop lies completely within an exon. We split the exon
          # into two
          my $exon_left = Bio::EnsEMBL::Exon->new(-slice     => $exon->slice,
                                                  -start     => $exon->start,
                                                  -end       => $stop->start - 1,
                                                  -strand    => $exon->strand,
                                                  -phase     => $exon->phase,
                                                  -end_phase => 0);
          my $exon_right = Bio::EnsEMBL::Exon->new(-slice     => $exon->slice,
                                                   -start     => $stop->end + 1,
                                                   -end       => $exon->end,
                                                   -strand    => $exon->strand,
                                                   -phase     => 0,
                                                   -end_phase => $exon->end_phase);
          
          # need to split the supporting features between the two
          my @sfs = @{$exon->get_all_supporting_features};
          my (@ug_left, @ug_right);
          foreach my $f (@sfs) {
            foreach my $ug ($f->ungapped_features) {

              if ($ug->start >= $exon_left->start and $ug->end <= $exon_left->end) {
                #completely within the left-side of the split
                push @ug_left, $ug;
              } elsif ($ug->start >= $exon_right->start and $ug->end <= $exon_right->end) {
                #completely within the right-side of the split
                push @ug_right, $ug;
              } else {
                # this ug must span the split
                my $fp_left = Bio::EnsEMBL::FeaturePair->new();
                $fp_left->seqname   ($ug->seqname);
                $fp_left->strand    ($ug->strand);
                $fp_left->hseqname  ($ug->hseqname);
                $fp_left->score     ($ug->score);
                $fp_left->percent_id($ug->percent_id);
                $fp_left->start     ($ug->start);
                $fp_left->end       ($stop->start - 1);
                
                my $fp_right = Bio::EnsEMBL::FeaturePair->new();
                $fp_right->seqname   ($ug->seqname);
                $fp_right->strand    ($ug->strand);
                $fp_right->hseqname  ($ug->hseqname);
                $fp_right->score     ($ug->score);
                $fp_right->percent_id($ug->percent_id);
                $fp_right->start     ($stop->end + 1);
                $fp_right->end       ($ug->end);
                
                if ($exon->strand > 1) {
                  $fp_left->hstart  ($ug->hstart);
                  $fp_left->hend    ($fp_left->hstart + (($fp_left->end - $fp_left->start + 1)/3) - 1);

                  $fp_right->hend   ($ug->hend);
                  $fp_right->hstart ($ug->hend - (($fp_right->end - $fp_right->start + 1)/3) + 1);
                } else {
                  $fp_left->hend    ($ug->hend);
                  $fp_left->hstart  ($ug->hend - (($fp_left->end - $fp_left->start + 1)/3) + 1);

                  $fp_right->hstart ($ug->hstart);
                  $fp_right->hend    ($fp_right->hstart + (($fp_right->end - $fp_right->start + 1)/3) - 1);
                }

                if ($fp_left->end >= $fp_left->start) { 
                  push @ug_left, $fp_left;
                }
                if ($fp_right->end >= $fp_right->start) {
                  push @ug_right, $fp_right;
                }
              }
            }
          }
          
          $exon_left->add_supporting_features(Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@ug_left));
          $exon_right->add_supporting_features(Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@ug_right));

          if ($exon->strand < 0) {
            if ($exon_right->end >= $exon_right->start) {
              push @new_exons, $exon_right;
            }
            if ($exon_left->end >= $exon_left->start) {
              push @new_exons, $exon_left;
            }
          } else {
            if ($exon_left->end >= $exon_left->start) {
              push @new_exons, $exon_left;
            }
            if ($exon_right->end >= $exon_right->start) {
              push @new_exons, $exon_right;
            } 
          }
        } else {
          # this exon is unaffected by this stop
          push @new_exons, $exon;
        }
      }

      @exons = @new_exons;
      $had_stops = 1;
    }

    if ($had_stops) {
      $tran->flush_Exons;
      foreach my $exon (@exons) {
        $tran->add_Exon($exon);
      }
      my $translation = Bio::EnsEMBL::Translation->new();
      $translation->start_Exon($exons[0]);
      $translation->end_Exon($exons[-1]);
      $translation->start(1);
      $translation->end($exons[-1]->end - $exons[-1]->start + 1);
      $tran->translation($translation);
    }

    push @processed_transcripts, $tran;
  }

  return \@processed_transcripts;
}




=head2 filter_transcripts

 Title   : filter_transcripts
    Takes an initial "raw" set of transcript fragments and filters them to form a
    mutually consistent set (i.e. where all members could feasibly be part
    part of the same transcript).    
    NOTE: not necessary for Genewise, which places each a.a. of the target
    protein only once
=cut

sub filter_transcripts {
  my ($self, $transcripts) = @_;

  my @sorted_trans;

  foreach my $tran (@$transcripts) {
    # transcript will only have one supporting feature for use cases of this filter
    my ($sf) = @{$tran->get_all_supporting_features};

    push @sorted_trans, [$tran, $sf];
  }

  @sorted_trans = sort { $b->[1]->score <=> $a->[1]->score } @sorted_trans;
    
  my (@filtered_trans);

  TRANSCRIPT:
  foreach my $tran_arr (@sorted_trans) {
    my ($tran, $sf) = @$tran_arr;
    my $keep_tran = 1;

    foreach my $ftran_arr (@filtered_trans) {
      my ($ftran, $fsf) = @$ftran_arr;

      # must be same strand
      $keep_tran = 0 if $sf->strand != $fsf->strand;

      
      # must not overlap in query or target

      $keep_tran = 0 if $sf->hstart <= $fsf->hend and $sf->hend >= $fsf->hstart;
      $keep_tran = 0 if $sf->start <= $fsf->end and $sf->end >= $fsf->start;

      if ($sf->strand < 0) {
        if ($sf->start > $fsf->end) {
          $keep_tran = 0 if $fsf->hstart < $sf->hend;
        } else {
          $keep_tran = 0 if $sf->hstart < $fsf->hend;
        }
      } else {
        if ($sf->start > $fsf->end) {
          $keep_tran = 0 if $sf->hstart < $fsf->hend;
        } else {
          $keep_tran = 0 if $fsf->hstart < $sf->hend;
        }
      }
    }

    push @filtered_trans, $tran_arr if $keep_tran;
  }

  return [map { $_->[0] } @filtered_trans];
}


#################
# checks
#################


sub check_consistent_coords {
  my ($self, $cj, $ci) = @_; 

  my $consistent = 1;

  if ($ci->id ne $cj->id or
      $ci->strand != $cj->strand) {
    $consistent = 0;
  } 
  else {
    if ($ci->strand == 1) {
      # j should come before i in genomic coords
      if (not $cj->end < $ci->start) {
        $consistent = 0;
      }
    } else {
      # j should come after i
      if (not $ci->end < $cj->start) {
        $consistent = 0;
      }
    }
  }

  return $consistent;
}


sub separate_coords {
  my ($self, $ci, $cj, $target_slice) = @_;

  my ($ci_slice, $cj_slice, $between_slice);

  if (defined $ci and defined $cj) {
    $ci_slice = $target_slice->adaptor->fetch_by_region('toplevel',
                                                        $target_slice->seq_region_name,
                                                        $ci->start,
                                                        $ci->end);

    $cj_slice = $target_slice->adaptor->fetch_by_region('toplevel',
                                                        $target_slice->seq_region_name,
                                                        $cj->start,
                                                        $cj->end);
    $between_slice = $target_slice->adaptor->fetch_by_region('toplevel',
                                                             $target_slice->seq_region_name,
                                                             $ci->end + 1,
                                                             $cj->start - 1);
  } elsif (defined $ci) {
    $ci_slice = $target_slice->adaptor->fetch_by_region('toplevel',
                                                        $target_slice->seq_region_name,
                                                        $ci->start,
                                                        $ci->end);
    $between_slice = $target_slice->adaptor->fetch_by_region('toplevel',
                                                             $target_slice->seq_region_name,
                                                             $ci->end + 1,
                                                             $target_slice->length);    

  } elsif (defined $cj) {
    $cj_slice = $target_slice->adaptor->fetch_by_region('toplevel',
                                                        $target_slice->seq_region_name,
                                                        $cj->start,
                                                        $cj->end);
    $between_slice = $target_slice->adaptor->fetch_by_region('toplevel',
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
      if (not @between or $between[0]->{from_start} > $between_slice->start) {
        $new_left = Bio::EnsEMBL::Mapper::Coordinate->new($ci->id,
                                                          $ci->start,
                                                          $ci->end,
                                                          $ci->strand);
      } else {
        $new_left = Bio::EnsEMBL::Mapper::Coordinate->new($ci->id,
                                                          $ci->start,
                                                          $between[0]->{from_end},
                                                          $ci->strand);
      }
    }

    if (defined $cj) {
      if (not @between or $between[-1]->{from_end} < $between_slice->end) {
        $new_right = Bio::EnsEMBL::Mapper::Coordinate->new($cj->id,
                                                           $cj->start,
                                                           $cj->end,
                                                           $cj->strand);
      } else {
        $new_right = Bio::EnsEMBL::Mapper::Coordinate->new($cj->id,
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


sub merge_coords {
  my ($self, $ci, $cj) = @_;

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



###########################################################
# write report; for development purposes
###########################################################

sub write_report {
  my ($self,
      $iteration, 
      $gene,
      $gene_scaffold,
      $gs_mapper,
      $projected_cds_regs,
      $new_trans) = @_;

  my $header = $gene->stable_id . "-" . $iteration;

  ####################################################
  #
  # first, summarise some properties of the input data
  # 
  #####################################################
   
  ###############
  # Initial candidate min-assembly
  ###############
  my @regions = $gs_mapper->map_coordinates($gene_scaffold->seq_region_name,
                                            1,
                                            $gene_scaffold->length,
                                            1,
                                            "genescaffold");

  if (0) {
    print "\n$header \#DETAILED_PRIOR_ASSEMBLY\n";
    my $last_pos = 0;
    foreach my $coord (@regions) {
      if ($coord->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
        printf("$header   %s %s : %s %d %d %s\n", 
               $last_pos + 1,
               $last_pos + $coord->length,
               $coord->id, 
               $coord->start, 
               $coord->end, 
               $coord->strand);
      } else {
        printf "$header   GAP/%d-%d\n", $coord->start, $coord->end;
      }
      $last_pos += $coord->length;
    }
    print "\n";
  }
    

  ##############
  # summarise where the CDS regions in the query lie on the target(s)
  ##############
  my ($total_cds_len, $total_coverage) = (0,0);

  if (0) {
    print "$header \#PROJECTED CDS regions\n";
    

    foreach my $ex_els (@$projected_cds_regs) {
      $total_cds_len += $ex_els->[-1]->{query}->end - $ex_els->[0]->{query}->start + 1;
      
      printf("$header  %s/%d-%d\t",
             $gene->seq_region_name,
             $ex_els->[0]->{query}->start,
             $ex_els->[-1]->{query}->end);
      
      foreach my $pair (@$ex_els) {
        if ($pair->{target}->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
          printf(" : %s\t%d\t%d\t%d",
                 $pair->{target}->id,
                 $pair->{target}->start,
                 $pair->{target}->end,
                 $pair->{target}->strand);
          $total_coverage += $pair->{query}->end - $pair->{query}->start + 1;
        } else {
          printf(" : %s\t%d\t%d", "GAP : ", $pair->{query}->start, $pair->{query}->end);
        }
      }
      print "\n";
    }
    printf("$header  Total coverage of source CDS regs by alignment = %.2f\n", ($total_coverage / $total_cds_len) * 100);
  }

  printf("$header \#PER-TRANSCRIPT RESULTS\n");
  foreach my $tran_obj (@$new_trans) {
    my $original_t = $tran_obj->{original};
    my $projected = $tran_obj->{projected};
    my @realigned = @{$tran_obj->{realigned}};

    my ($t_sf) = @{$projected->get_all_supporting_features};

    ####################################################
    # Summarise the result of projected the transcript
    #####################################################
   
    print "$header  PROJECTION:\n";
    printf("$header  Transcript (%s %d-%d) : %s/%d-%d %d\n",
          $t_sf->hseqname,
          $t_sf->hstart, 
          $t_sf->hend,
          $gene_scaffold->seq_region_name,
          $projected->start,
          $projected->end,
          $projected->strand);


    my @exons = sort {$a->start <=> $b->start} @{$projected->get_all_Exons};

    for(my $j=0; $j < @exons; $j++) {
      my $exon = $exons[$j];

      my $ex_num = $projected->strand < 0 ? @exons - $j : $j+1;

      my ($e_sf) = @{$exon->get_all_supporting_features};

      my @map_segs = $gs_mapper->map_coordinates($exon->slice->seq_region_name,
                                                 $exon->start,
                                                 $exon->end,
                                                 1,
                                                 "genescaffold");
        
      foreach my $map_seg (@map_segs) {
        if ($map_seg->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
          printf("$header  EXON $ex_num (pep %d-%d, gs %d %d) : %s\t%d\t%d\t%d\n",
                 defined $e_sf ? $e_sf->hstart : 0,
                 defined $e_sf ? $e_sf->hend   : 0,
                 $exon->start,
                 $exon->end,
                 $map_seg->id, 
                 $map_seg->start,
                 $map_seg->end,
                 $map_seg->strand);

        } else {
          printf("$header  EXON $ex_num (pep %d-%d, gs %d %d) : GAP/%s\n",
                 defined $e_sf ? $e_sf->hstart : 0,
                 defined $e_sf ? $e_sf->hend   : 0,
                 $exon->start,
                 $exon->end,
                 $map_seg->length);
        }
      }  
    }

    my $translation = $projected->translate;
    my $gaps = $translation->seq =~ tr/X/X/;
    my $stops = $translation->seq =~ tr/\*/\*/;

    printf "$header Total coverage of parent peptide = %.2f\n", $t_sf->hcoverage;
    printf "$header Total number of stops = %d\n", $stops;
    if ($gaps < $translation->length) {
      printf "$header Total proportion of non-gap      = %.2f\n", 100 - 100 * ($gaps / $translation->length);

      my $seqio = Bio::SeqIO->new(-format => 'fasta',
                                  -fh => \*STDOUT);
      my $seq = Bio::PrimarySeq->new();
      $seq->seq($translation->seq);
      $seq->id($original_t->stable_id);
      $seq->desc($gene->stable_id);

      $seqio->write_seq($seq);
      $seqio->close;
    } else {
      print "Total proportion of non-gap = 0.00\n";
    }


    ###################################################
    #
    # Now summarise the result of realigning the peptide
    #
    #####################################################

    if (0) {
      printf "$header Realignment:\n";
      printf("$header  Total aligned transcripts = %d\n",  scalar(@realigned));
      
      my @sorted_ts = sort {$a->start <=> $b->start} @realigned;
      
      $total_coverage = 0;
      for (my $i=0; $i < @sorted_ts; $i++) {
        my $tran = $sorted_ts[$i];
        
        ($t_sf) = @{$tran->get_all_supporting_features};
        
        $total_coverage += $t_sf->hcoverage;
        
        printf("$header  Trancript $i (%s %d-%d) : %s/%d-%d %d\n",
               $t_sf->hseqname,
               $t_sf->hstart, 
               $t_sf->hend,
               $gene_scaffold->seq_region_name,
               $tran->start,
               $tran->end,
               $tran->strand);
        
        @exons = sort {$a->start <=> $b->start} @{$tran->get_all_Exons};
        
        my @implicated_scaffolds;
      
        for (my $j = 0; $j < @exons; $j++) {
          my $exon = $exons[$j];
          
          my ($e_sf) = @{$exon->get_all_supporting_features};
          
          my @map_segs = $gs_mapper->map_coordinates($gene_scaffold->seq_region_name,
                                                     $exon->start,
                                                     $exon->end,
                                                     1,
                                                     "genescaffold");
          
          foreach my $map_seg (@map_segs) {
            if ($map_seg->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
              printf("$header  EXON $j (pep %d-%d, gs %d %d) : %s\t%d\t%d\t%d\n",
                     defined $e_sf ? $e_sf->hstart : 0,
                     defined $e_sf ? $e_sf->hend   : 0,
                     $exon->start,
                     $exon->end,
                     $map_seg->id, 
                     $map_seg->start,
                     $map_seg->end,
                     $map_seg->strand);
              
            } else {
              printf("$header  EXON $j (pep %d-%d) : GAP/%s\n",
                     defined $e_sf ? $e_sf->hstart : 0,
                     defined $e_sf ? $e_sf->hend   : 0,
                     $map_seg->length);
            }
          }      
        }
      }    
      
      printf("$header Total coverage of parent peptide of %s by all transcripts: %.2f\n", 
             $original_t->stable_id,
             $total_coverage);
    }
  }
}



###########################
# utility methods
###########################

sub get_all_transcript_cds_features {
  my ($self, @transcripts) = @_;

  my (@orig_cds_feats, @merged_cds_feats);

  foreach my $tran (@transcripts) {
    foreach my $exon (@{$tran->get_all_translateable_Exons}) {
      push @orig_cds_feats, Bio::EnsEMBL::Feature->
          new(-start  => $exon->start,
              -end    => $exon->end,
              -slice  => $exon->slice);
    }
  }

  foreach my $feat (sort {$a->start <=> $b->start} @orig_cds_feats) {
    if (not @merged_cds_feats or $feat->start > $merged_cds_feats[-1]->end) {
      push @merged_cds_feats, $feat;
    } else {
      if ($feat->end > $merged_cds_feats[-1]->end) {
        $merged_cds_feats[-1]->end($feat->end);
      }
    }
  }

  return \@merged_cds_feats;  
}



#####################################
# get/sets
#####################################


sub get_all_Transcripts {
  my ($self, $gene, $single) = @_;

  # if the singel flag is given returns the longest
  # coding transcript. Otherwise returns all coding
  # transcripts. 

  if (not exists $self->{_transcripts}) {
    $self->{_transcripts} = {};
  }
  
  if (not exists $self->{_transcripts}->{$gene}) {

    my ($longest, $longest_len, @transcripts);
    
    foreach my $t (@{$gene->get_all_Transcripts}) {
      my $translation = $t->translate;
      
      if (defined $translation) {
        if (not defined $longest or
            $translation->length > $longest_len) {
          
          $longest = $t;
          $longest_len = $translation->length;
        }

        push @transcripts, $t;
      }
    }
    
    if ($single) {
      @transcripts = ($longest);
    }
    $self->{_transcripts}->{$gene} = \@transcripts;
  } 

  return $self->{_transcripts}->{$gene};

}


sub get_translation {
  my ($self, $trans) = @_;

  if (not exists $self->{_translations}) {
    $self->{_translations} = {};
  }

  if (not exists $self->{_translations}->{$trans}) {

    my $pep_seq = $trans->translation->seq;
    # convert selenosysteine gene to a aa that 
    # can be gracefully handled by alignment programs!
    $pep_seq =~ s/U/C/g;
    
    my $pep = Bio::PrimarySeq->new(-id => $trans->stable_id,
                                   -seq => $pep_seq);
    
    $self->{_translations}->{$trans} = $pep;
  }

  return $self->{_translations}->{$trans};

}


sub projected_cds {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_projected_cds} = $val;
  }
  
  return $self->{_projected_cds};
}

sub gene {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_gene} = $val;
  }

  return $self->{_gene};
}


sub genomic_align_block_chains {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_gen_al_chains} = $val;
  }

  return $self->{_gen_al_chains};
}


sub query_top_level_slice {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_query_tl_slice} = $val;
  }

  return $self->{_query_tl_slice};
}


sub target_slices {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_target_slices} = $val;
  }
  
  if (not exists $self->{_target_slices}) {
    $self->{_target_slices} = {};
  }

  return $self->{_target_slices};
}


sub mapper {
  my ($self, $val) = @_;

    if (defined $val) {
    $self->{_mapper} = $val;
  }

  return $self->{_mapper};

}



####################################
# config variable holders
####################################

sub read_and_check_config {
  my ($self, $hash) = @_;

  $self->SUPER::read_and_check_config($hash);
 
  my $logic = $self->analysis->logic_name;

  foreach my $var (qw(INPUT_METHOD_LINK_TYPE
                      QUERY_CORE_DB
                      TARGET_CORE_DB
                      COMPARA_DB)) {

    throw("You must define $var in config for logic '$logic'" . 
          " or in the DEFAULT entry")
        if not $self->$var;
  }
}

sub USE_GENEWISE {
  my ($self, $val) = @_;
 
  if (defined $val) {
    $self->{_use_genewise} = $val;
  }

  return $self->{_use_genewise};
}


sub GENE_SCAFFOLD_COMPONENT_PADDING {
  my ($self, $val) = @_;
  
  if (defined $val) {
    $self->{_inter_scaffold_padding} = $val;
  }

  return $self->{_inter_scaffold_padding};
}



sub INPUT_METHOD_LINK_TYPE {
  my ($self, $type) = @_;

  if (defined $type) {
    $self->{_input_method_link_type} = $type;
  }

  return $self->{_input_method_link_type};
}


sub COMPARA_DB {
  my ($self, $db) = @_;

  if (defined $db) {
    $self->{_compara_db} = $db;
  }

  return $self->{_compara_db};
}



sub QUERY_CORE_DB {
  my ($self, $db) = @_;

  if (defined $db) {
    $self->{_query_core_db} = $db;
  }

  return $self->{_query_core_db};
}


sub TARGET_CORE_DB {
  my ($self, $db) = @_;

  if (defined $db) {
    $self->{_target_core_db} = $db;
  }

  return $self->{_target_core_db};
}


1;
