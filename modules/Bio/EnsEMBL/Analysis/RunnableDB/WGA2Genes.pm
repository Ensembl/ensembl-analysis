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
  my $genscan = Bio::EnsEMBL::Analysis::RunnableDB::WGA2Genes
    ->new (-db         => $pipelinedb,
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

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use Bio::EnsEMBL::Analysis::Tools::WGA2Genes::GeneScaffold;
#use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;

use Bio::EnsEMBL::Analysis::Tools::WGA2Genes::CoordUtils;
use Bio::EnsEMBL::Analysis::Tools::WGA2Genes::ChainUtils;

use Bio::EnsEMBL::Analysis::Tools::Logger;

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

  my $q_dbh = Bio::EnsEMBL::DBSQL::DBAdaptor->
      new(%{$self->QUERY_CORE_DB});
  my $t_dbh = Bio::EnsEMBL::DBSQL::DBAdaptor->
      new(%{$self->TARGET_CORE_DB});
  my $compara_dbh = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->
      new(%{$self->COMPARA_DB});
  
  my $query_species = 
      $q_dbh->get_MetaContainerAdaptor->get_Species->binomial;
  my $target_species = 
      $t_dbh->get_MetaContainerAdaptor->get_Species->binomial;
  
  my $gdb_adap = $compara_dbh->get_GenomeDBAdaptor;
  my $q_gdb = $gdb_adap->fetch_by_name_assembly($query_species);
  my $t_gdb = $gdb_adap->fetch_by_name_assembly($target_species);

  ################################################
  # check that the default assembly for the query and target agrees 
  # with that for the method_link_species_set GenomeDBs
  #################################################

  my ($q_assembly_version, $t_assembly_version);
  eval {
    $q_assembly_version = $q_dbh->get_CoordSystemAdaptor->
        fetch_by_name('toplevel',
                      $q_gdb->assembly);
    
    $t_assembly_version = $t_dbh->get_CoordSystemAdaptor->
        fetch_by_name('toplevel',
                      $t_gdb->assembly);
  };
  $@ and do {
    throw("Had trouble fetching coord systems for ". 
          $q_gdb->assembly . " and " . 
          $t_gdb->assembly . " from core dbs: $@");
  };

  #####
  # fetch the gene; need to work in the coordinate space of the 
  # top-level slice to be consistent with compara
  #####
  my $sa = $q_dbh->get_SliceAdaptor;

  if ($input_id =~ /:/) {
    # assume slice name

    my $slice = $sa->fetch_by_name($input_id);
    $self->query_slice($sa->fetch_by_region('toplevel',
                                            $slice->seq_region_name));


    my @genes;
    foreach my $g (@{$slice->get_all_Genes}) {
      next if $g->biotype ne 'protein_coding';
      $g = $g->transfer($self->query_slice);

      foreach my $t (@{$self->get_all_Transcripts($g)}) {
        next if $t->coding_region_start < $slice->start;
        next if $t->coding_region_end   > $slice->end;
      }

      push @genes, $g;
    }

    $self->genes(\@genes);

  } else {
    # assume iid is a gene stable id
    my ($gene);
    eval {
      $gene = $q_dbh->get_GeneAdaptor->fetch_by_stable_id($input_id);
    };
    if ($@ or not defined $gene) {
      throw("Could not find gene '$input_id' in query database");
    }
    $self->query_slice($sa->fetch_by_region('toplevel', 
                                            $gene->seq_region_name));
    $gene = $gene->transfer($self->query_slice);

    $self->genes([$gene]);
  }

  if ($self->REJECT_BAD_QUERY_TRANSCRIPTS) {
    foreach my $g (@{$self->genes}) {
      $self->filter_bad_Transcripts($g);
    }
  }

  my ($reg_start, $reg_end);
  foreach my $g (@{$self->genes}) {
    foreach my $t (@{$self->get_all_Transcripts($g)}) {
      my @e = @{$t->get_all_Exons};

      $reg_start = $t->start 
          if not defined $reg_start or $t->start < $reg_start;
      $reg_end   = $t->end   
          if not defined $reg_end   or $t->end   > $reg_end;
    }
  }

  if ($self->PSEUDOGENE_CHAIN_FILTER) {
    my @repeat_features;
    my $res_slice = $sa->fetch_by_region('toplevel',
                                         $self->query_slice->seq_region_name,
                                         $reg_start,
                                         $reg_end);
    foreach my $rf (@{$res_slice->get_all_RepeatFeatures('RepeatMask')}) {
      $rf = $rf->transfer($self->query_slice);
      
      next if $rf->start < $self->query_slice->start;
      next if $rf->end > $self->query_slice->end;
      
      push @repeat_features, $rf;
    }  
    $self->query_repeats(\@repeat_features);
  }

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

  my $dnafrag = $compara_dbh->get_DnaFragAdaptor->
      fetch_by_GenomeDB_and_name($q_gdb,
                                 $self->query_slice->seq_region_name);

  my $gaba = $compara_dbh->get_GenomicAlignBlockAdaptor;

  my $gen_al_blocks = 
      $gaba->fetch_all_by_MethodLinkSpeciesSet_DnaFrag($mlss,
                                                       $dnafrag,
                                                       $reg_start,
                                                       $reg_end);

  my (%chains, @chains);
  foreach my $block (@$gen_al_blocks) {
    my $qga = $block->reference_genomic_align;    
    my ($tga) = @{$block->get_all_non_reference_genomic_aligns};

    # fetch the target slice for later reference
    if (not exists $self->target_slices->{$tga->dnafrag->name}) {

      $self->target_slices->{$tga->dnafrag->name} = 
          $t_dbh->get_SliceAdaptor->fetch_by_region('toplevel',
                                                    $tga->dnafrag->name);
    }

    my $chain_id =  $qga->genomic_align_group_id_by_type("chain");

    if ($block->reference_genomic_align->dnafrag_strand < 0) {
      $block->reverse_complement;
    }

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
  @chains = sort { $b->[0]->score <=> $a->[0]->score } @chains;  
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

  my (@results);

  my $alignment_chains = $self->genomic_align_block_chains;

  logger_info("ORIGINAL_CHAINS\n" . stringify_chains($alignment_chains));

  if ($self->PSEUDOGENE_CHAIN_FILTER) {
    $alignment_chains = 
        filter_pseudogene_chains($alignment_chains,
                                 $self->query_repeats,
                                 $self->PCF_MAX_FUSION_INTERVAL,
                                 $self->PCF_MAX_REPEAT_IN_INTERVAL
                                 );
    
    logger_info( "NON_PG_CHAINS\n" . stringify_chains($alignment_chains));
  }

  # segment the genes into overlapping groups
  my @gene_sets = @{$self->get_individual_gene_sets};

  foreach my $geneset (@gene_sets) {
    my (%blocks_used_so_far);
    
    ITER:for(my $iteration=0; ;$iteration++) {

      logger_info("** Iteration $iteration\n");

      my $max_editable_stops = ($iteration == 0)
          ? $self->MAX_EDITABLE_STOPS_PRIMARY
          : $self->MAX_EDITABLE_STOPS_NON_PRIMARY;

      my $filtered_chains = filter_block_overlap_chains($alignment_chains,
                                                        [values %blocks_used_so_far]);

      logger_info("NON_USED_CHAINS\n" . stringify_chains($filtered_chains));
      
      my $gs_name = $geneset->[0]->{stable_id} . "-" . $iteration;

      my @these_genes;
      foreach my $g (@$geneset) {
        push @these_genes, {
          stable_id => $g->stable_id, 
          transcripts => $self->get_all_Transcripts($g),
          proj_transcripts => [],
          good_sources => [],
        }
      }

      # Repeat transcript contruction until we had no rejected transcripts
      # This is because if a transcript is rejected, the underlying 
      # GeneScaffold components are possible no longer all necessary
      for(;;) {
        my @these_trans = map { @{$_->{transcripts}} } @these_genes;

        logger_info("Working with transcripts " . join(" ", map { $_->stable_id} @these_trans));

        my @cds_feats = @{$self->get_all_transcript_cds_features(@these_trans)};

        $filtered_chains = 
            filter_irrelevant_chains($filtered_chains,
                                     \@cds_feats);        
        logger_info("RELEVANT_CHAINS\n" . stringify_chains($filtered_chains));
        
        $filtered_chains = filter_inconsistent_chains($filtered_chains, 
                                                      $self->OVERLAP_CHAIN_FILTER);
        logger_info("CONSISTENT_CHAINS\n" . stringify_chains($filtered_chains));
        
        if ($self->NO_CONTIG_SPLITS) { 
          $filtered_chains = $self->remove_contig_split_chains($filtered_chains);
          logger_info("NON_CONTIG_SPLIT_CHAINS\n" . stringify_chains($filtered_chains));
        }

        my $net_blocks = flatten_chains($filtered_chains, 1);

        my ($projected_cds_feats, $blocks_used_this_iter) = 
            $self->map_features_to_target(\@cds_feats, $net_blocks);

        # if there is no coverage of the CDS regions, finish
        last ITER if not keys %$blocks_used_this_iter;
        
        my $gene_scaffold = 
                $self->gene_scaffold_from_projection($projected_cds_feats,
                                                     $gs_name);
                                                     
      
        last ITER if not defined $gene_scaffold;

        my $result = {
          gene_scaffold   => $gene_scaffold,
          genes           => [],
        };
        
        my $had_to_reject_one = 0;
        foreach my $gene (@these_genes) {
          
          foreach my $tran (@{$gene->{transcripts}}) {
            my $proj_trans = 
                $gene_scaffold->project_transcript($tran,
                                                   $self->MAX_EXON_READTHROUGH_DIST);
            $proj_trans = 
                $self->process_transcript($proj_trans, 
                                          $max_editable_stops,
                                          $self->MIN_COVERAGE,
                                          $self->MIN_NON_GAP,
                                          $tran->stable_id);
            
            if ($proj_trans) {
              push @{$gene->{good_sources}}, $tran;
              push @{$gene->{proj_transcripts}}, $proj_trans;
            } else {
              $had_to_reject_one = 1;
            }
          }
        }

        # check that no source projections were rejected; if so, go again
        if ($had_to_reject_one) {
          my @kept_genes;
          foreach my $gene (@these_genes) {
            $gene->{transcripts} = $gene->{good_sources};
            $gene->{good_sources} = [];
            $gene->{proj_transcripts} = [];
            if (@{$gene->{transcripts}}) {
              push @kept_genes, $gene;
            }
          }
          @these_genes = @kept_genes;
          next;
        }

        foreach my $gene (@these_genes) {
          my @transcripts =  @{$self->make_nr_transcript_set($gene->{proj_transcripts},
                                                            $gene_scaffold)};
            
          my $new_gene_name = 
              $gene_scaffold->seq_region_name . "." . 
              $gene->{stable_id};
          
          push @{$result->{genes}}, { 
            name        => $new_gene_name,
            transcripts => \@transcripts,
          };
        }
        
        push @results, $result;

        foreach my $block_id (keys %$blocks_used_this_iter) {
          $blocks_used_so_far{$block_id} = $blocks_used_this_iter->{$block_id};
        }
        last;
      }
    }
  }

  $self->output(\@results);
}


=head2 write_output
    Title   :   write_output
    Usage   :   $self->write_output
    Returns :   nothing
    Args    :   none
=cut

sub write_output {
  my ($self) = @_;

  print  "#\n";
  printf("# WGA2Genes output for %s\n", $self->input_id);
  print "#\n";

  foreach my $obj (@{$self->output}) {
    my $gs = $obj->{gene_scaffold};

    my @genes = @{$obj->{genes}};

    $self->write_agp($gs);
    foreach my $g (@genes) {
      $self->write_gene($gs, 
                        $g->{name}, 
                        @{$g->{transcripts}});      
    }
  }

  return;
}


######################################
# internal methods
#####################################

###################################################################
# FUNCTION: get_individual_gene_sets
#
# Decription:
#  segments the reference gene set into sets that are 
#  treated independently by the method. This is only done
#  if the SEGMENT_BY_GENE_OVERLAP is set
#  non-overlapping on the query sequence (using our own internal 
#  get_all_Transcripts method that removes transcripts with 
#  outlier introns
###################################################################
sub get_individual_gene_sets {
  my ($self) = @_;

  my (@gene_objs, $max_end);
  foreach my $g (@{$self->genes}) {
    my ($g_start, $g_end);
    foreach my $t (@{$self->get_all_Transcripts($g)}) {
      foreach my $e (@{$t->get_all_translateable_Exons}) {
        if (not defined $g_start or $e->start < $g_start) {
          $g_start = $e->start;
        }
        if (not defined $g_end or $e->end > $g_end) {
          $g_end = $e->end; 
        }
      }
    }

    if (defined $g_start and defined $g_end) {
      push @gene_objs, { 
        gene  => $g,
        start => $g_start,
        end   => $g_end,
      };
    }
  }
  @gene_objs = sort { $a->{start} <=> $b->{start} } @gene_objs;

  if ($self->SEGMENT_SLICE_BY_GENE_OVERLAP) {
    my @gene_lists;

    foreach my $g_obj (@gene_objs) {      
      if (defined $max_end and $g_obj->{start} < $max_end) {
        push @{$gene_lists[-1]}, $g_obj->{gene};
        $max_end = $g_obj->{end} if $g_obj->{end} > $max_end;
      } else {
        push @gene_lists, [$g_obj->{gene}];
        $max_end = $g_obj->{end};
      }
    }

    return \@gene_lists;
  } else {
    my @genes = map { $_->{gene} } @gene_objs;

    return [\@genes];
  }
}


###################################################################
# FUNCTION: map_features_to_target
#
# Decription:
#    Takes a feature list and a GenomicAlignBlock list and returns 
#    a list of elements (one for each feature) comprising 
#    {query => [], target = []}, segments, representing the gapped 
#    alignment of the query (in the region of the feature) to the 
#    target. The alignment is represented a list of 
#    Bio::EnsEMBL::Mapper::Coordinates and Bio::EnsEMBL::Mapper::Gap
#    for each of the query and target
#
###################################################################

sub map_features_to_target {
  my ($self, $feats, $gen_al_blocks) = @_;

  my (@proj_feats, %overlapping_blocks);

  foreach my $feat (@$feats) {

    my @projected_feat_segments;

    # binary search the block list
    #my ($left, $j) = (0, scalar(@$gen_al_blocks)-1);
    #while ($j - $left > 0) {
    #  my $mid = int(($left + $j) / 2);
    #  if ($gen_al_blocks->[$mid]->reference_genomic_align->dnafrag_end < 
    #      $feat->start) {
    #    $left = $mid + 1;
    #  } else {
    #    $j = $mid;
    #  }
    #}
    my $j=0;
    
    my (@overlapping_blocks);

    for(my $i=$j; $i < @$gen_al_blocks; $i++) {
      my $block = $gen_al_blocks->[$i];
      
      my $qga = $block->reference_genomic_align;
      my ($tga) = @{$block->get_all_non_reference_genomic_aligns};
      
      next if $qga->dnafrag_end < $feat->start;
      last if $qga->dnafrag_start > $feat->end;
      
      push @overlapping_blocks, $block;
    }

    foreach my $block (@overlapping_blocks) {
      my $qga = $block->reference_genomic_align;
      my ($tga) = @{$block->get_all_non_reference_genomic_aligns};

      my ($q_feat_start, $q_feat_end);
      if ($qga->dnafrag_start < $feat->start) { 
        $q_feat_start = $feat->start; 
      } else {
        $q_feat_start = $qga->dnafrag_start;
      }
      if ($qga->dnafrag_end > $feat->end) {
        $q_feat_end = $feat->end;
      } else {
        $q_feat_end = $qga->dnafrag_end;
      }
      
      $overlapping_blocks{$block} = $block;
      
      my $current_pos = $q_feat_start;
      my @coords = $qga->get_Mapper->
          map_coordinates("sequence",
                          $q_feat_start,
                          $q_feat_end,
                          $qga->dnafrag_strand,
                          "sequence");
      
      foreach my $coord (@coords) {
        
        if ($coord->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
          my $current_pos2 = $current_pos;
          my @coords2 = $tga->get_Mapper->
              map_coordinates("alignment",
                              $coord->start,
                              $coord->end,
                              1,
                              "alignment");
          
          foreach my $coord2 (@coords2) {
            
            if ($coord2->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
              push @projected_feat_segments, { 
                query  => Bio::EnsEMBL::Mapper::Coordinate->
                    new($feat->slice->seq_region_name,
                        $current_pos2,
                        $current_pos2 + ($coord2->end - $coord2->start),
                        1),
                target => Bio::EnsEMBL::Mapper::Coordinate->
                    new($tga->dnafrag->name,
                        $coord2->start,
                        $coord2->end,
                        $qga->dnafrag_strand * $tga->dnafrag_strand),
                level  => $tga->level_id,
              };
            }
            
            $current_pos2 += ($coord2->end - $coord2->start + 1);
          }          
        }
        
        $current_pos += ($coord->end - $coord->start + 1);
      }
    }
  
    #########
    # fill in Gap segments between the coordinates so that each original 
    # feat is completely accounted for by scaffold segments and gap 
    # segments
    #########
    
    if (not @projected_feat_segments) {
      
      push @projected_feat_segments, {
        query  => Bio::EnsEMBL::Mapper::Coordinate->
            new($feat->slice->seq_region_name,
                $feat->start,
                $feat->end,
                1),
            target => Bio::EnsEMBL::Mapper::Gap->new($feat->start, 
                                                     $feat->end),
          };
    } 
    else {
      # start
      if ($feat->start < $projected_feat_segments[0]->{query}->start) {
        my $st = $feat->start;
        my $en = $projected_feat_segments[0]->{query}->start - 1;
        
        unshift @projected_feat_segments, {
          query  => Bio::EnsEMBL::Mapper::Coordinate->
              new($feat->slice->seq_region_name,
                  $st,
                  $en,
                  1),
              target => Bio::EnsEMBL::Mapper::Gap->new($st, $en),
            };
      }
      
      # end
      if ($feat->end > $projected_feat_segments[-1]->{query}->end) {
        my $st = $projected_feat_segments[-1]->{query}->end + 1;
        my $en = $feat->end;
        
        push @projected_feat_segments, {
          query  => Bio::EnsEMBL::Mapper::Coordinate->
              new($feat->slice->seq_region_name,
                  $st,
                  $en,
                1),
              target => Bio::EnsEMBL::Mapper::Gap->new($st, $en),
            };
      }
      
      # middle
      my @new_bits;
      for(my $i=1; $i < @projected_feat_segments; $i++) {
        my ($l, $r) = ($projected_feat_segments[$i-1], 
                       $projected_feat_segments[$i]);
        
        if ($r->{query}->start > $l->{query}->end + 1) {
          my $st = $l->{query}->end + 1;
          my $en = $r->{query}->start - 1;
          
          push @new_bits, {
            query  => Bio::EnsEMBL::Mapper::Coordinate->
                new($feat->slice->seq_region_name,
                    $st,
                    $en,
                  1),
                target => Bio::EnsEMBL::Mapper::Gap->new($st, $en),
              };
        }
      }
      
      push @projected_feat_segments, @new_bits;
    }
    
    @projected_feat_segments = sort {
      $a->{query}->start <=> $b->{query}->start;
    } @projected_feat_segments;
   
    push @proj_feats, \@projected_feat_segments;
  }

  return (\@proj_feats, \%overlapping_blocks);
}


###################################################################
# FUNCTION:  gene_scaffold_from_projection
#
# Description:
#   Takes a list of exon projection elements (representing the gapped
#   alignment between the query and target in the region of the exons)
#   and returns a GeneScaffold object which is a mini-assembly of the 
#   target sequences necessary to present the whole gene on a single 
#   sequence
###################################################################

sub gene_scaffold_from_projection {
  my ($self, $projected_cds_elements, $gene_scaffold_id) = @_;

  my (@projected_cds_elements, @targets, @new_targets);

  $gene_scaffold_id = "NoName" if not defined $gene_scaffold_id;

  # find max length for each component scaffolds used in the projection
  #
  #my %max_lengths;
  #foreach my $cds (@$projected_cds_elements) {
  #  foreach my $coord_pair (@$cds) {
  #    my $tcoord = $coord_pair->{target};
  #    if ($tcoord->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
  #      if (not exists $max_lengths{$tcoord->id} or
  #          $max_lengths{$tcoord->id} < $tcoord->length) {
  #        $max_lengths{$tcoord->id} = $tcoord->length;
  #      }
  #    }
  #  }
  #}

  # flatten list and add exon number to each projected component
  my $cds_id = 0;
  foreach my $cds (@$projected_cds_elements) {
    push @projected_cds_elements, [];
    foreach my $coord_pair (@$cds) {
      push @targets, {
        cds_id  => $cds_id,
        coord   => $coord_pair->{target},
      };
      push @{$projected_cds_elements[-1]}, $coord_pair;
    }
    $cds_id++;
  }

  # remove gaps which cannot be filled. Non-fillable gaps:
  # 1. If one of the flanking coords is on the same CDS region
  #    as the gap, and the end of the aligned region does
  #    not align to a sequence-level gap
  # 2. If the 2 flanking coords are consistent and not
  #    separated by a sequence-level gap
  @new_targets = ();
  for(my $i=0; $i < @targets; $i++) {
    my $this_target = $targets[$i];

    if ($this_target->{coord}->isa("Bio::EnsEMBL::Mapper::Gap")) {
      # if it's gap that can be filled, leave it. Otherwise, remove it
      my ($left_target, $right_target);
      for(my $j=$i-1; $j>=0; $j--) {
        if ($targets[$j]->{coord}->
            isa("Bio::EnsEMBL::Mapper::Coordinate")) {
          $left_target = $targets[$j];
          last;
        }
      }
      for(my $j=$i+1; $j < @targets; $j++) {
        if ($targets[$j]->{coord}->
            isa("Bio::EnsEMBL::Mapper::Coordinate")) {
          $right_target = $targets[$j];
          last;
        }
      }

      my $keep_gap = 1;

      if (defined $left_target and 
          defined $right_target and
          check_consistent_coords($left_target->{coord}, 
                                  $right_target->{coord})) {
        
        my ($lc, $rc) = ($left_target->{coord}, 
                         $right_target->{coord});
        if ($lc->strand < 0) {
          ($lc, $rc) = ($rc, $lc);
        }
        my ($new_left, $new_right) = 
            separate_coords($lc, 
                            $rc,
                            $self->target_slices->{$lc->id});
        
        if (not defined $new_left and not defined $new_right) {
          $keep_gap = 0;
        }
      } 


      if ($keep_gap and 
          defined $left_target and 
          $left_target->{cds_id} == $this_target->{cds_id}) {
        
        my ($new_left, $dummy);
        my $lc = $left_target->{coord};
        if ($lc->strand < 0) {
          ($dummy, $new_left) = 
              separate_coords(undef,
                              $lc,
                              $self->target_slices->{$lc->id});
          
          if (abs($new_left->start - $lc->start) > 10) {
            # we're not near the end of a contig. 
            # This gap cannot be filled;
            $keep_gap = 0;
          }
        } else {
          ($new_left, $dummy) = 
              separate_coords($lc,
                              undef,
                              $self->target_slices->{$lc->id});
          
          if (abs($new_left->end - $lc->end) > 10) {
            # we're not near the end of a contig. 
            # This gap cannot be filled;
            $keep_gap = 0;
          }
        }
      }

      if ($keep_gap and 
          defined $right_target and
          $right_target->{cds_id} == $this_target->{cds_id}) {

        my ($dummy, $new_right);
        my $rc = $right_target->{coord};
        if ($rc->{strand} < 0) {
          ($new_right, $dummy) = 
              separate_coords($rc,
                              undef,
                              $self->target_slices->{$rc->id});
          
          if (abs($rc->end - $new_right->end) > 10) {
            # we're not near the end of a contig. 
            # This gap cannot be filled;
            $keep_gap = 0;
          }
        } else {
          ($dummy, $new_right) = 
              separate_coords(undef,
                              $rc,
                              $self->target_slices->{$rc->id});
          
          if (abs($rc->start - $new_right->start) > 10) {
            # we're not near the end of a contig. 
            # This gap cannot be filled;
            $keep_gap = 0;
          }
        }
      }
      
      if ($keep_gap) {
        push @new_targets, $targets[$i];
      }

    } else {
      push @new_targets, $targets[$i];
    }
  }
  @targets = @new_targets;

  # prune away gaps at the start and end; we won't "fill" these
  while(@targets and 
        $targets[0]->{coord}->isa("Bio::EnsEMBL::Mapper::Gap")) {
    shift @targets;
  }
  while(@targets and 
        $targets[-1]->{coord}->isa("Bio::EnsEMBL::Mapper::Gap")) {
    pop @targets;
  }

  return () if not @targets;

  # merge adjacent targets   
  #  we want to be able to account for small, frame-preserving 
  #  insertions in the target sequence with respect to the query. 
  #  To give the later, gene-projection code the opportunity to 
  #  "read through" these insertions, we have to merge togther 
  #  adjacent, consistent targets that are within this "maximum 
  # read-through" distance
  @new_targets = ();
  for(my $i=0; $i<@targets; $i++) {
    my $target = $targets[$i];

    if ($target->{coord}->isa("Bio::EnsEMBL::Mapper::Coordinate") and
        @new_targets and
        $new_targets[-1]->{coord}->
           isa("Bio::EnsEMBL::Mapper::Coordinate") and 
        check_consistent_coords($new_targets[-1]->{coord}, 
                                $target->{coord})) {
      
      my $dist = distance_between_coords($new_targets[-1]->{coord}, 
                                         $target->{coord});
      
      if ($dist <= $self->MAX_EXON_READTHROUGH_DIST) {
        
        my $last_target = pop @new_targets;
        my $new_coord = merge_coords($last_target->{coord},
                                     $target->{coord});
        
        # check that the new merged coord will not result in an overlap
        my $overlap = 0;
        foreach my $tg (@new_targets) {
          my $tg_c = $tg->{coord};
          if ($tg_c->isa("Bio::EnsEMBL::Mapper::Coordinate") and
              $tg_c->id eq $new_coord->id and
              $tg_c->start < $new_coord->end and
              $tg_c->end   > $new_coord->start) {
            $overlap = 1;
            last;
          }
        }
        if (not $overlap) {
          for (my $j=$i+1; $j < @targets; $j++) {
            my $tg_c = $targets[$j]->{coord};
            if ($tg_c->isa("Bio::EnsEMBL::Mapper::Coordinate") and
                $tg_c->id eq $new_coord->id and 
                $tg_c->start < $new_coord->end and
                $tg_c->end   > $new_coord->start) {
              $overlap = 1;
              last;
            }
          }
        }
        
        if ($overlap) {
          push @new_targets, $last_target, $target;
        } else {
          push @new_targets,  {
            coord => $new_coord,
          };
        }
      } else {
        push @new_targets, $target;
      }
    }
    else {
      push @new_targets, $target;
    }
  }
    
  @targets = @new_targets;

  my @components = map { $_->{coord} } @targets;
  my @ref_comp_coord_pairs;
  foreach my $cds_el (@projected_cds_elements) {
    push @ref_comp_coord_pairs, map { [$_->{query}, $_->{target}] } @$cds_el;
  }

  my $gene_scaffold = Bio::EnsEMBL::Analysis::Tools::WGA2Genes::GeneScaffold
      ->new(-name => $gene_scaffold_id,
            -component_slices => $self->target_slices,
            -component_coords => \@components,
            -reference_slice => $self->query_slice,
            -alignment_coords => \@ref_comp_coord_pairs);
            
  return $gene_scaffold;
}



#################################################################
# FUNCTION  : process_transcript
#
# Description:
#    Subjects the given transcript to a few tests, returning 
#    the transcript if they succeed, undef if not. If the
#    transcript contains less than $max_stops stops, these
#    are "spliced out"; otherwise the transcripts is rejected
#################################################################

sub process_transcript {
  my ($self, 
      $tran, 
      $max_stops,
      $min_coverage,
      $min_non_gap,
      $source_id) = @_;
  
  return 0 if not $tran;

  my (@processed_transcripts, $num_stops, $prop_non_gap);

  my @exons = @{$tran->get_all_Exons};
  my ($tsf) = @{$tran->get_all_supporting_features};
  my $pep = $tran->translate->seq;

  ##################
  # reject transcripts that have:
  #   zero length
  #   less than minimum coverage of parent peptide
  #   higher that maximum proportion of gap residues
  ##################

  if (length($pep) == 0) {
    logger_info("Rejecting proj of $source_id because was just a stop codon");
    return 0;
  }

  if ($tsf->hcoverage < $min_coverage) {
    logger_info("Rejecting proj of $source_id because coverage (" . 
                $tsf->hcoverage . 
                ") is below threshold");
    return 0;
  }

  foreach my $attr (@{$tran->get_all_Attributes}) {
    if ($attr->code eq "PropNonGap") {
      $prop_non_gap = $attr->value;
    } elsif ($attr->code eq "NumStops") {
      $num_stops = $attr->value;
    }
  }

  if ($prop_non_gap < $min_non_gap) {
    logger_info("Rejecting proj of $source_id due to too many Ns ($prop_non_gap)");
    return 0;
  }

  if ($num_stops > $max_stops) {
    logger_info("Rejecting proj of $source_id due to too many stops ($num_stops)");
    return 0;
  } elsif ($num_stops == 0) {
    return $tran;
  } 

  ##################
  # number of stops is non-zero but acceptable. Need to 
  # operate on the transcript to jump over the stops
  ##################

  while($pep =~ /\*/g) {
    my $position = pos($pep);

    my @coords = $tran->pep2genomic($position, $position);

    if (@coords > 1) {
      # the codon is split by an intron. Messy. Leave these for now
      logger_info("Rejecting proj of $source_id due to stop interrupted by intron");
      return 0;
    } 
    my ($stop) = @coords;

    # locate the exon that this stop lies in
    my @new_exons;
    foreach my $exon (@exons) {
      if ($stop->start >= $exon->start and $stop->end <= $exon->end) {
        # this stop lies completely within an exon. We split the exon
        # into two
        my $exon_left = Bio::EnsEMBL::Exon->
            new(-slice     => $exon->slice,
                -start     => $exon->start,
                -end       => $stop->start - 1,
                -strand    => $exon->strand,
                -phase     => $exon->strand < 0 ? 0 : $exon->phase,
                -end_phase => $exon->strand < 0 ? $exon->end_phase  :0);
        my $exon_right = Bio::EnsEMBL::Exon->
            new(-slice     => $exon->slice,
                -start     => $stop->end + 1,
                -end       => $exon->end,
                -strand    => $exon->strand,
                -phase     => $exon->strand < 0 ? $exon->phase : 0,
                -end_phase => $exon->strand < 0 ? 0 : $exon->end_phase);
        # need to split the supporting features between the two
        my @sfs = @{$exon->get_all_supporting_features};
        my (@ug_left, @ug_right);
        foreach my $f (@sfs) {
          foreach my $ug ($f->ungapped_features) {
            if ($ug->start >= $exon_left->start and 
                $ug->end <= $exon_left->end) {
              #completely within the left-side of the split
              push @ug_left, $ug;
            } elsif ($ug->start >= $exon_right->start and 
                     $ug->end <= $exon_right->end) {
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
              
              if ($exon->strand > 0) {
                $fp_left->hstart($ug->hstart);
                $fp_left->hend($fp_left->hstart +
                               ($fp_left->length / 3) - 
                               1);
                
                $fp_right->hend ($ug->hend);
                $fp_right->hstart($ug->hend - 
                                  ($fp_right->length / 3) + 
                                  1);
              } else {
                $fp_left->hend ($ug->hend);
                $fp_left->hstart($ug->hend - 
                                 ($fp_left->length / 3) + 
                                 1);
                
                $fp_right->hstart($ug->hstart);
                $fp_right->hend($fp_right->hstart +
                                ($fp_right->length / 3) - 
                                1);
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

        if (@ug_left) {
          my $f = Bio::EnsEMBL::DnaPepAlignFeature->
              new(-features => \@ug_left);
          $exon_left->add_supporting_features($f);
        }
        if (@ug_right) {
          my $f = Bio::EnsEMBL::DnaPepAlignFeature->
              new(-features => \@ug_right);
          $exon_right->add_supporting_features($f);
        }
        
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
  }
  
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

  return $tran;
}




###################################################################
# FUNCTION   : make_nr_transcript_set
#
# Description : 
#    Takes an initial "raw", ordered set of transcripts and proceeds 
#    through the list, rejecting redundant transcripts. The aim
#    is to remove transcript that are redundant when considering
#    gap exons (for exmaple, 2 transcripts that differ in the location 
#    of the 4th exon, but where that exon is a gap exon in one of 
#    those cases, should be considered the same).
###################################################################

sub make_nr_transcript_set {
  my ($self, $transcripts, $gene_scaf) = @_;

  my (@t_objs, @kept_t_objs);

  TRAN:foreach my $tran (@$transcripts) {
    # reject if transcript has an identical structure to one already seen

    my (@non_gap_exons, @non_gap_introns);
    my $non_gap_length = 0;
    my $gap_length = 0;

    my $last_exon;
    foreach my $e (sort {$a->start <=> $b->start} @{$tran->get_all_Exons}) {      
      my ($coord) = $gene_scaf->component_mapper->map_coordinates($gene_scaf->seq_region_name,
                                                                  $e->start,
                                                                  $e->end,
                                                                  1,
                                                                  $gene_scaf->seq_region_name);

      if ($coord->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
        push @non_gap_exons, $e;
        $non_gap_length += $e->length;
        
        if (defined $last_exon) {
          my $intron = Bio::EnsEMBL::Feature->
              new(-start => $last_exon->end + 1,
                  -end   => $e->start - 1);
          push @non_gap_introns, $intron;
        }
        $last_exon = $e;

      } else {
        $gap_length += $e->length;
        $last_exon = undef;
      }
    }

    push @t_objs, {
      transcript => $tran,
      non_gap_length => $non_gap_length, 
      gap_length => $gap_length,
      non_gap_exons => \@non_gap_exons,
      non_gap_introns => \@non_gap_introns,
    };
  }
  
  @t_objs = sort { 
    $b->{non_gap_length} <=> $a->{non_gap_length} or
    $a->{gap_length} <=> $b->{gap_length};
  } @t_objs;

  foreach my $t_obj (@t_objs) {

    my @exons = @{$t_obj->{non_gap_exons}};
    my @introns = @{$t_obj->{non_gap_introns}};

    my $found_subsuming_transcript = 0;

    KTRANS: foreach my $k_t_obj (@kept_t_objs) {    
      my @k_exons = @{$k_t_obj->{non_gap_exons}};
      my @k_introns = @{$k_t_obj->{non_gap_introns}};

      foreach my $e (@exons) {
        my $found_match = 0;

        foreach my $ke (@k_exons) {
          if ($e->start == $ke->start and
              $e->end   == $ke->end) {
            $found_match = 1;
            last;
          }
        }

        if (not $found_match) {
          next KTRANS;
        }
      }
      # if we get here, the transcript had no unique exons;
      # now check the introns
      foreach my $i (@introns) {
        my $found_match = 0;
        foreach my $ki (@k_introns) {
          if ($i->start == $ki->start and
              $i->end == $ki->end) {
            $found_match = 1;
            last;
          }
        }

        if (not $found_match) {
          next KTRANS;
        }
      }
      
      # if we get here, the transcript had no unique introns
      # either. 
      $found_subsuming_transcript = 1;
      last;
    }

    if (not $found_subsuming_transcript) {
      push @kept_t_objs, $t_obj;
    }
  }
  
  return [map {$_->{transcript}} @kept_t_objs];
}


#################################################################
# FUNCTION   : get_all_transcript_cds_features
#
# Description : 
#    Returns a list of Bio::EnsEMBL::Features representing the 
#    projection of the coding regions of the given transcripts 
#    onto the query 
#################################################################

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
    if (not @merged_cds_feats or 
        $feat->start > $merged_cds_feats[-1]->end) {
      push @merged_cds_feats, $feat;
    } else {
      if ($feat->end > $merged_cds_feats[-1]->end) {
        $merged_cds_feats[-1]->end($feat->end);
      }
    }
  }

  return \@merged_cds_feats;  
}



#################################################################
# FUNCTION   : get_all_Transcripts
#
# Description : 
#   Returns the set of transcripts for the given gene. We store 
#   these explicitly here for each gene so that individual ones
#   can be removed if necessary
#################################################################

sub get_all_Transcripts {
  my ($self, $gene) = @_;

  # Return all coding transcripts for the gene
  # pruning "problematic" transcripts (defintion
  # of problematic: causes the gene to overlap
  # with another gene, and removal of the trans
  # leaves at least one remaining non-problematic
  # transcript)

  if (not exists $self->{_transcripts}) {
    $self->{_transcripts} = {};
  }
  
  if (not exists $self->{_transcripts}->{$gene}) {
    $self->{_transcripts}->{$gene} = [];

    my @trans = @{$gene->get_all_Transcripts};

    # even if we dont want to reject potentially dodgy transcripts,
    # we have to remove the ones that are not mod-3 in length because
    # thy screw up downstream code
    my @kept_trans;

    foreach my $t (@trans) {
      next if not defined $t->translate;
      
      my @exons = sort {
        $a->start <=> $b->start;
      } @{$t->get_all_translateable_Exons};
      
      my $total_cds_len = 0;
      map { $total_cds_len += $_->length } @exons;

      if ($total_cds_len % 3 == 0) {
        push @kept_trans, $t;
      }
    }      

    $self->{_transcripts}->{$gene} = \@kept_trans;
  }

  return $self->{_transcripts}->{$gene};
}


#################################################################
# FUNCTION   : filter_bad_Transcripts
#
# Description : 
#   Attempts to remove "dodgy" source transcripts from the input gene
#   (i.e. those can cause overlap with another gene, where the
#   other transcripts do not). 
#################################################################


sub filter_bad_Transcripts {
  my ($self, $gene) = @_;
    
  my (@transcripts);
      
  my @orig_trans = @{$self->get_all_Transcripts($gene)};

  foreach my $t (@orig_trans) {
    # record large introns in the transcript
    my @introns;
    
    my @exons = sort {
      $a->start <=> $b->start;
    } @{$t->get_all_translateable_Exons};
    
    for(my $i=1; $i<@exons; $i++) {
      my $intron_start = $exons[$i-1]->end + 1;
      my $intron_end   = $exons[$i]->start - 1;
      
      push @introns, { 
        start => $intron_start,
        end   => $intron_end,
      };
    }
    
    push @transcripts, { 
      tran   => $t,
      introns => \@introns,
    };
  }
  
  my (@kept_transcripts, @maybe_keep_transcripts);
  foreach my $t (@transcripts) {
    if (not @{$t->{introns}}) {
      push @kept_transcripts, $t;
    } else {
      # also keep if the large introns do not overlap other genes
      my $intron_contains_gene = 0;
      OTHERGENE: foreach my $og (@{$self->genes}) {
        next if $og->stable_id eq $gene->stable_id;
        
        foreach my $li (@{$t->{introns}}) {
          if ($og->start <= $li->{end} and 
              $og->end   >=  $li->{start}) {
            $intron_contains_gene = 1;
            last OTHERGENE;
          }
        }
      }
      if (not $intron_contains_gene) {
        push @kept_transcripts, $t;
      } else {
        # this transcript has intron-containing genes, but we might
        # have to keep it if all such transcripts have this property
        push @maybe_keep_transcripts, $t;
      }
    }
  }
  
  if (not @kept_transcripts) {
    @kept_transcripts = @maybe_keep_transcripts;
  }
  
  $self->{_transcripts}->{$gene} = [map {$_->{tran}} @kept_transcripts];

}


###################################################################
# FUNCTION: remove_contig_split_chains
#
# Decription:
#  This function removes chains that (a) interlock with higher
#  scoring chains, and (b) this interlocking cannot be explained
#  by contig mis-ordering within the involved scaffolds
###################################################################
sub remove_contig_split_chains {
  my ($self, $input_chains) = @_;

  # plan:
  # for each chain
  #  add to list of chains retained so far
  #  sort chains by target id
  #  merge consistent chains with the same target id
  #  separate chains that are interrupted by some other block
  #  check that all separated chains of the same target id can be 
  #   split at the contig level

  my @kept_chains;

  foreach my $c (@$input_chains) {
    my (%coords_by_tid, @all_coords, @merged_coords);

    my $consistent = 1;

    my @these_chains = (@kept_chains, $c);

    my $blocks = flatten_chains(\@these_chains, 1);

    foreach my $b (@$blocks) {
      my ($tga) = @{$b->get_all_non_reference_genomic_aligns};
      my $tgac = Bio::EnsEMBL::Mapper::Coordinate->new
          ($tga->dnafrag->name,
           $tga->dnafrag_start,
           $tga->dnafrag_end,
           $tga->dnafrag_strand);
      
      push @{$coords_by_tid{$tgac->id}}, $tgac;
      push @all_coords, $tgac;
    }

    foreach my $tgac (@all_coords) {
      my $merged = 0;

      if (@merged_coords and 
          $merged_coords[-1]->id eq $tgac->id and
          $merged_coords[-1]->strand eq $tgac->strand) {

        if ($tgac->strand > 0) {
          if ($tgac->start > $merged_coords[-1]->end) {
            my $can_merge = 1;
            foreach my $o_c (@{$coords_by_tid{$tgac->id}}) {
              if ($o_c->start > $merged_coords[-1]->end and
                  $o_c->end   < $tgac->start) {
                $can_merge = 0;
                last;
              }
            }
            if ($can_merge) {
              $merged_coords[-1]->end($tgac->end);
              $merged = 1;
            }
          }
        } else {
          if ($tgac->end < $merged_coords[-1]->start) {
            my $can_merge = 1;
            foreach my $o_c (@{$coords_by_tid{$tgac->id}}) {
              if ($o_c->start > $tgac->end and
                  $o_c->end < $merged_coords[-1]->start) {
                $can_merge = 0;
                last;
              }
            }
            if ($can_merge) {
              $merged_coords[-1]->start($tgac->start);
              $merged = 1;
            }
          }
        }
      }

      if (not $merged) {
        push @merged_coords, $tgac;
      }
    }

    %coords_by_tid = ();
    foreach my $c (@merged_coords) {
      push @{$coords_by_tid{$c->id}}, $c;
    }

    TID: foreach my $id (keys %coords_by_tid) {
      my @chains = @{$coords_by_tid{$id}};
      @chains = sort { $a->start <=> $b->start } @chains; 

      # each chain should be separable from the last at contig
      # level
      for(my $i=1; $i < @chains; $i++) {
        my $prev = $chains[$i-1];
        my $this = $chains[$i];

        my ($left, $right) = 
            separate_coords($prev, 
                            $this,
                            $self->target_slices->{$id});
        
        if (not defined $left or not defined $right) {
          $consistent = 0;
          last TID;
        }
      }
    }

    if ($consistent) {
      push @kept_chains, $c;
    }
  }
    
  return \@kept_chains;
}



###########################################################
# write methods
###########################################################

sub write_agp {
  my ($self, 
      $gscaf) = @_;

  my $fh = \*STDOUT;

  my $prefix = "##-AGP";

  my @tpieces = $gscaf->component_coords;
  my @qpieces = $gscaf->reference_coords;
  
  print $fh "$prefix \#\n";
  printf($fh "$prefix \#\# AGP for gene scaffold %s source-region=%s/%d-%d\n", 
         $gscaf->seq_region_name,
         $qpieces[0]->id, 
         $qpieces[0]->start <= $qpieces[-1]->start ? $qpieces[0]->start : $qpieces[-1]->start,
         $qpieces[-1]->end >= $qpieces[0]->end ? $qpieces[-1]->end : $qpieces[0]->end);

  print $fh "$prefix \#\n";

  my $last_end = 0;
  for(my $i=0; $i < @tpieces; $i++) {
    my $piece = $tpieces[$i];

    printf($fh "$prefix %s\t%d\t%d\t",
           $gscaf->seq_region_name,
           $last_end + 1,
           $last_end + $piece->length,
           $i + 1);
        
    if ($piece->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
      printf($fh "%s\t%s\t%d\t%d\t%s\n", 
             "W",
             $piece->id,
             $piece->start,
             $piece->end,
             $piece->strand < 0 ? "-" : "+");
    } else {
      printf($fh "%s\t%d\n", 
             "N",
             $piece->length);
    }

    $last_end += $piece->length;
  }
}


sub write_gene {
  my ($self, 
      $gene_scaf,
      $gene_name,
      @trans) = @_;

  my $fh = \*STDOUT;

  my $prefix = "##-GENES ";
  
  my $seq_id = $gene_scaf->seq_region_name;
  my $gene_id = $gene_name;

  printf $fh "$prefix \# Gene report for $seq_id\n";
  
  foreach my $tran (@trans) {
    my ($sf) = @{$tran->get_all_supporting_features};
    my $tran_id = $gene_id . "_" . $sf->hseqname; 

    printf($fh "$prefix \##\-ATTRIBUTE transcript=$tran_id code=%s value=%.2f\n", 
           "HitCoverage",
           $sf->hcoverage);
    
    foreach my $attr (@{$tran->get_all_Attributes}) {
      printf($fh "$prefix \##\-ATTRIBUTE transcript=$tran_id code=%s value=%s\n",
             $attr->code, 
             $attr->value);
    }

    my @exons = @{$tran->get_all_Exons};
    
    for (my $i=0; $i < @exons; $i++) {
      my $exon = $exons[$i];
      my $exon_id = $tran_id . "_" . $i;
      
      printf($fh "%s %s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s=%s; %s=%s; %s=%s\n", 
             $prefix,
             $seq_id,
             "WGA2Genes",
             "Exon",
             $exon->start,
             $exon->end,
             $exon->strand,
             $exon->phase,
             $exon->end_phase,
             "exon",
             $exon_id,
             "transcript",
             $tran_id,
             "gene",
             $gene_id);
      foreach my $sup_feat (@{$exon->get_all_supporting_features}) {
        foreach my $ug ($sup_feat->ungapped_features) {
          printf($fh "%s %s\t%s\t%s\t%d\t%d\t%d\t%s=%s; %s=%s; %s=%s; %s=%s\n", 
                 $prefix,
                 $seq_id,
                 "WGA2Genes",
                 "Supporting",
                 $ug->start,
                 $ug->end,
                 $ug->strand,
                 "exon",
                 $exon_id,
                 "hname",
                 $ug->hseqname,
                 "hstart",
                 $ug->hstart,
                 "hend",
                 $ug->hend);
        }
      }      
    }
  }
}


###########################
# gets/sets
###########################

sub genomic_align_block_chains {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_gen_al_chains} = $val;
  }

  return $self->{_gen_al_chains};
}

sub genes {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_genes} = $val;
  }

  return $self->{_genes};
}

sub query_repeats {
  my ($self, $val) = @_;

  if (defined $val) {
    # make sure given list is non-overlapping
    my @nr;
    foreach my $rf (sort {$a->start <=> $b->start} @$val) {
      if (not @nr or $nr[-1]->end < $rf->start) {
        push @nr, $rf;
      } else {
        if ($rf->end > $nr[-1]->end) {
          $nr[-1]->end($rf->end);
        }
      }
    }

    $self->{_query_repeat_feats} = \@nr;
  }
  
  return $self->{_query_repeat_feats};
}

sub query_slice {
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

#
# core options
#

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


#
# Dictates whether we process slice as a whole, 
#
sub SEGMENT_SLICE_BY_GENE_OVERLAP {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_segment_slice_by_gene_overlap} = $val;
  }

  return $self->{_segment_slice_by_gene_overlap};
}


#
# Initial query gene-set filtering
#

sub REJECT_BAD_QUERY_TRANSCRIPTS {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_reject_bad_query_transcripts} = $val;
  }

  return $self->{_reject_bad_query_transcripts};
}


#
# chain filtering
#

sub PSEUDOGENE_CHAIN_FILTER {
  my ($self, $val) = @_;
  
  if (defined $val) {
    $self->{_pseudogene_chain_filter} = $val;
  }
}

sub PCF_MAX_FUSION_INTERVAL {
  my ($self, $val) = @_;
  
  if (defined $val) {
    $self->{_pcf_max_fusion_interval} = $val;
  }
  return $self->{_pcf_max_fusion_interval};
}

sub PCF_MAX_REPEAT_IN_INTERVAL {
  my ($self, $val) = @_;
  
  if (defined $val) {
    $self->{_pcf_max_repeat_in_interval} = $val;
  }
  return $self->{_pcf_max_repeat_in_interval};
}

sub OVERLAP_CHAIN_FILTER {
  my ($self, $val) = @_;
  
  if (defined $val) {
    $self->{_overlap_chain_filter} = $val;
  }
  if (not exists($self->{_overlap_chain_filter})) {
    $self->{_overlap_chain_filter} = 0;
  }

  return $self->{_overlap_chain_filter};
}


#
# transcript editing and filtering
#

sub MAX_EXON_READTHROUGH_DIST {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_max_exon_readthrough_dist} = $val; 
  }

  return $self->{_max_exon_readthrough_dist};
}

sub MAX_EDITABLE_STOPS_PRIMARY {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_max_editable_stops} = $val; 
  }

  return $self->{_max_editable_stops};
}

sub MAX_EDITABLE_STOPS_NON_PRIMARY {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_max_editable_stops_non_primary} = $val; 
  }

  return $self->{_max_editable_stops_non_primary};
}

sub MIN_COVERAGE {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_min_coverage} = $val;
  }

  return $self->{_min_coverage};
}

sub MIN_NON_GAP {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_min_non_gap} = $val;
  }

  return $self->{_min_non_gap};
}

#
# Gene scaffold options
#
sub NO_CONTIG_SPLITS {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_no_contig_splits} = $val;
  }

  return $self->{_no_contig_splits};
}

sub GENE_SCAFFOLD_INTER_PIECE_PADDING {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_interpiece_padding} = $val;
  }

  return $self->{_interpiece_padding};
}



1;
