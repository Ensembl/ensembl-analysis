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

use Bio::EnsEMBL::Utils::Sequence qw(reverse_comp);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript;
use Bio::EnsEMBL::Analysis::Runnable::AlignmentNets;
use Bio::EnsEMBL::Pipeline::Tools::TranscriptUtils;

use IO::String;

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

  #$self->print_chains($alignment_chains, "ORIGINAL CHAINS");

  if ($self->PSEUDOGENE_CHAIN_FILTER) {
    $alignment_chains = 
        $self->remove_pseudogene_chains($alignment_chains,
                                        $self->PCF_MAX_FUSION_INTERVAL,
                                        $self->PCF_MAX_REPEAT_IN_INTERVAL
                                        );
    
    #$self->print_chains($alignment_chains, "AFTER PG REMOVAL");
  }

  # segment the genes into overlapping groups
  my @gene_sets = @{$self->get_non_overlapping_gene_sets};

  foreach my $geneset (@gene_sets) {
    my (%blocks_used_so_far);
    
    ITER:for(my $iteration=0; ;$iteration++) {

      my $max_editable_stops = ($iteration == 0)
          ? $self->MAX_EDITABLE_STOPS_PRIMARY
          : $self->MAX_EDITABLE_STOPS_NON_PRIMARY;

      my $filtered_chains = $self->remove_used_chains($alignment_chains,
                                                      \%blocks_used_so_far);
      #$self->print_chains($filtered_chains, "CHAINS LEFT AFTER USED REMOVED");
      
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

        my @cds_feats = @{$self->get_all_transcript_cds_features(@these_trans)};

        $filtered_chains = 
            $self->remove_irrelevant_chains($filtered_chains,
                                            \@cds_feats);
        
        #$self->print_chains($alignment_chains, "RELEVANT CHAINS");
        
        $filtered_chains = $self->remove_interfering_chains($filtered_chains);
        
        # $self->print_chains($filtered_chains, "CONSISTENT CHAINS");
        
        if ($self->NO_CONTIG_SPLITS) { 
          $filtered_chains = $self->remove_contig_split_chains($filtered_chains);
          # $self->print_chains($filtered_chains, "NO CONTIG SPLIT CHAINS");
        }

        my $net_blocks = $self->flatten_chains($filtered_chains, 1);

        my ($projected_cds_feats, $blocks_used_this_iter) = 
            $self->map_features_to_target(\@cds_feats, $net_blocks);

        # if there is no coverage of the CDS regions, finish
        last ITER if not keys %$blocks_used_this_iter;
        
        my ($gene_scaffold, 
            $qy_to_gs_map,
            $tg_to_gs_map) = 
                $self->gene_scaffold_from_projection($gs_name,
                                                     $projected_cds_feats);
      
        last ITER if not defined $gene_scaffold;
        
        my $result = {
          gene_scaffold   => $gene_scaffold,
          mapper          => $tg_to_gs_map,
          genes           => [],
        };
        
        my $had_to_reject_one = 0;
        foreach my $gene (@these_genes) {
          my (@proj_trans, @good_sources);
          
          foreach my $tran (@{$gene->{transcripts}}) {
            my $proj_trans = 
                $self->make_projected_transcript($tran,
                                                 $gene_scaffold,
                                                 $qy_to_gs_map,
                                                 $tg_to_gs_map,
                                                 $gene->{stable_id});
            $proj_trans = 
                $self->process_transcript($proj_trans, $max_editable_stops);
            
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
                                                            $gene_scaffold,
                                                            $tg_to_gs_map)};
            
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
    my $map = $obj->{mapper};
    my @genes = @{$obj->{genes}};

    $self->write_agp(\*STDOUT, $gs, $map);
    foreach my $g (@genes) {
      $self->write_gene(\*STDOUT, 
                        $gs, 
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
# FUNCTION: get_non_overlapping_get_sets
#
# Decription:
#  segments the reference gene set into sets that are 
#  non-overlapping on the query sequence (using our own internal 
#  get_all_Transcripts method that removes transcripts with 
#  outlier introns
###################################################################
sub get_non_overlapping_gene_sets {
  my ($self) = @_;

  my (@gene_objs, @gene_lists, $max_end);
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
#   and returns 3 things:
#   1. A "gene scaffold" which is a mini-assembly of the target sequences
#      necessary to present the whole gene on a single sequence
#   2. Query mapper: a map between the coords of the gene scaffold and 
#      the original query sequence
#   3. Target mapper: a map between the coords of the gene scaffold and 
#      the comprised target sequences.
###################################################################

sub gene_scaffold_from_projection {
  my ($self, $gene_scaffold_id, $projected_cds_elements) = @_;

  my (@projected_cds_elements, @targets, @new_targets);

  # find max length for each component scaffolds used in the projection
  #
  my %max_lengths;
  foreach my $cds (@$projected_cds_elements) {
    foreach my $coord_pair (@$cds) {
      my $tcoord = $coord_pair->{target};
      if ($tcoord->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
        if (not exists $max_lengths{$tcoord->id} or
            $max_lengths{$tcoord->id} < $tcoord->length) {
          $max_lengths{$tcoord->id} = $tcoord->length;
        }
      }
    }
  }

  # flatten list and replace scaffolds that only occur in small 
  # components with gaps
  my $cds_id = 0;
  foreach my $cds (@$projected_cds_elements) {
    push @projected_cds_elements, [];
    foreach my $coord_pair (@$cds) {
      #my $tcoord = $coord_pair->{target};
      #my $qcoord = $coord_pair->{query};

      #if ($tcoord->isa("Bio::EnsEMBL::Mapper::Coordinate") and
      #    $max_lengths{$tcoord->id} < $self->MIN_COMPONENT_SIZE) {
      #  $coord_pair->{target} = 
      #      Bio::EnsEMBL::Mapper::Gap->new($coord_pair->{target}->start,
      #                                     $coord_pair->{target}->end);
      #}
      #push @targets, {
      #  cds_id  => $cds_id,
      #  coord   => $coord_pair->{target},
      #};

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
          $self->check_consistent_coords($left_target->{coord}, 
                                         $right_target->{coord})) {
        
        my ($lc, $rc) = ($left_target->{coord}, 
                         $right_target->{coord});
        if ($lc->strand < 0) {
          ($lc, $rc) = ($rc, $lc);
        }
        my ($new_left, $new_right) = 
            $self->separate_coords($lc, 
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
              $self->separate_coords(undef,
                                     $lc,
                                     $self->target_slices->{$lc->id});

          if (abs($new_left->start - $lc->start) > 10) {
            # we're not near the end of a contig. 
            # This gap cannot be filled;
            $keep_gap = 0;
          }
        } else {
          ($new_left, $dummy) = 
              $self->separate_coords($lc,
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
              $self->separate_coords($rc,
                                     undef,
                                     $self->target_slices->{$rc->id});

          if (abs($rc->end - $new_right->end) > 10) {
            # we're not near the end of a contig. 
            # This gap cannot be filled;
            $keep_gap = 0;
          }
        } else {
          ($dummy, $new_right) = 
              $self->separate_coords(undef,
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
        $self->check_consistent_coords($new_targets[-1]->{coord}, 
                                       $target->{coord})) {
      
      my $dist = $self->distance_between_coords($new_targets[-1]->{coord}, 
                                                $target->{coord});
      
      if ($dist <= $self->MAX_EXON_READTHROUGH_DIST) {
        
        my $last_target = pop @new_targets;
        my $new_coord = $self->merge_coords($last_target->{coord},
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

  ###################################
  # build the gene scaffold, and mappings between the new gene
  # scaffold and each of query and target coords

  my $t_map = Bio::EnsEMBL::Mapper->new('target',
                                        $gene_scaffold_id); 
  my $q_map = Bio::EnsEMBL::Mapper->new('query',
                                        $gene_scaffold_id);
  
  my ($seq, $last_end_pos) = ("", 0);
  for(my $i=0; $i < @targets; $i++) {
    my $target = $targets[$i];
    my $coord = $target->{coord};

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
      $t_map->add_map_coordinates($coord->id,
                                  $coord->start,
                                  $coord->end,
                                  $coord->strand,
                                  $gene_scaffold_id,
                                  $last_end_pos + 1,
                                  $last_end_pos + $coord->length);
    } else {
      # the sequence itself
      $seq .= ('n' x $coord->length);
      
      # and the map. This is a target gap we have "filled", so no position 
      # in target, but a position in query
      
      $q_map->add_map_coordinates($self->query_slice->seq_region_name,
                                  $coord->start,
                                  $coord->end,
                                  1,
                                  $gene_scaffold_id,
                                  $last_end_pos + 1,
                                  $last_end_pos + $coord->length);
    }

    # add padding between the pieces
    if ($i < @targets - 1) {
      $last_end_pos += 
          $coord->length + $self->GENE_SCAFFOLD_INTER_PIECE_PADDING;
      $seq .= ('n' x $self->GENE_SCAFFOLD_INTER_PIECE_PADDING);
    }
  }
  

  # now add of the original exon pieces to the query map
  foreach my $el (@projected_cds_elements) {
    foreach my $pair (@$el) {
      my ($q, $t) = ($pair->{query}, $pair->{target});

      if ($t->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
        # get the gene_scaffold position from the target map
        my ($coord) = $t_map->map_coordinates($t->id,
                                              $t->start,
                                              $t->end,
                                              $t->strand,
                                              'target');

        $q_map->add_map_coordinates($self->query_slice->seq_region_name,
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
      ->new(-coord_system => 
              Bio::EnsEMBL::CoordSystem->new(-name => 'genescaffold',
                                             -rank => 1),
            -seq_region_name => $gene_scaffold_id,
            -seq  => $seq,
            -start => 1,
            -end   => length($seq));

  return ($gene_scaffold, $q_map, $t_map);
}


###################################################################
# FUNCTION   : make_projected_transcript
#
# Description:
#    Takes a transcript, a gene scaffold, and a mapping between 
#    query coords and gene scaffold coords. Produces the transcript 
#    that is the result of "projecting" the original transcript, 
#    through alignment, onto the gene scaffold. 
###################################################################

sub make_projected_transcript {
  my ($self, $tran, $gene_scaf, $qmap, $tmap, $gene_id) = @_;

  my ($tran_length, @all_coords, @new_exons);

  my @orig_exons = @{$tran->get_all_translateable_Exons};
  if ($tran->strand < 0) {
    @orig_exons = reverse @orig_exons;
  }
  map { $tran_length += $_->length } @orig_exons; 

  foreach my $orig_exon (@orig_exons) {    
    my @crds = $qmap->map_coordinates($orig_exon->slice->seq_region_name,
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
      my ($tcoord) = $tmap->map_coordinates($gene_scaf->seq_region_name,
                                            $all_coords[$i]->start,
                                            $all_coords[$i]->end,
                                            1,
                                            $gene_scaf->seq_region_name);
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
      my ($tcoord) = $tmap->map_coordinates($gene_scaf->seq_region_name,
                                            $all_coords[$i]->start,
                                            $all_coords[$i]->end,
                                            1,
                                            $gene_scaf->seq_region_name);
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
                                               -slice => $gene_scaf);

      $total_tran_bps += $coord->length;
      my ($tcoord) = $tmap->map_coordinates($gene_scaf->seq_region_name,
                                            $coord->start,
                                            $coord->end,
                                            1,
                                            $gene_scaf->seq_region_name);
      if ($tcoord->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
        $real_seq_bps += $coord->length;
      } 
    }
  }

  if (not @new_exons) {
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

      my @gen_coords = $qmap->map_coordinates($gene_scaf->seq_region_name,
                                              $extent_start,
                                              $extent_end,
                                              1,
                                              $gene_scaf->seq_region_name);

      my @fps;
      my $cur_gs_start = $extent_start;
      foreach my $g_coord (@gen_coords) {
        my $cur_gs_end = $cur_gs_start + $g_coord->length - 1;
                
        if ($g_coord->isa("Bio::EnsEMBL::Mapper::Coordinate")) {

          my ($p_coord) = $tran->genomic2pep($g_coord->start, 
                                             $g_coord->end,
                                             $exon->strand);
          
          my $fp = Bio::EnsEMBL::FeaturePair->
              new(-seqname  => $gene_scaf->seq_region_name,
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
            $intron_len <= $self->MAX_EXON_READTHROUGH_DIST) { 
          $new_start = $exon->start;
          $new_end   = $prev_exon->end;
        }
      } else {
        my $intron_len = $exon->start - $prev_exon->end - 1;
        if ($intron_len % 3 == 0 and 
            $intron_len <= $self->MAX_EXON_READTHROUGH_DIST) {
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
    my ($coord) = $tmap->map_coordinates($gene_scaf->seq_region_name,
                                         $e->start,
                                         $e->end,
                                         1,
                                         $gene_scaf->seq_region_name);
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

  if (defined $gene_id) {
    my $geneid_attr = Bio::EnsEMBL::Attribute->
        new(-code => 'SourceGene',
            -name => 'source gene',
            -description => 'human source gene',
            -value => $gene_id);
    push @attributes, $geneid_attr;
  }
  my $tranid_attr = Bio::EnsEMBL::Attribute->
      new(-code => 'SourceTran',
          -name => 'source transcript',
          -description => 'source transcript',
          -value => $tran->stable_id);
  push @attributes, $tranid_attr;

  $proj_tran->add_Attributes(@attributes);

  return $proj_tran;
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
  my ($self, $tran, $max_stops) = @_;
  
  return 0 if not $tran;

  my (@processed_transcripts, $num_stops);

  my @exons = @{$tran->get_all_Exons};
  my ($tsf) = @{$tran->get_all_supporting_features};
  my $pep = $tran->translate->seq;

  ##################
  # reject transcripts that have:
  #   zero length
  #   less than minimum coverage of parent peptide
  #   higher that maximum proportion of gap residues
  ##################

  return 0 if length($pep) == 0;
  return 0 if $tsf->hcoverage < $self->MIN_COVERAGE;

  foreach my $attr (@{$tran->get_all_Attributes}) {
    if ($attr->code eq "PropNonGap") {
      return 0 if $attr->value < $self->MIN_NON_GAP;
    } elsif ($attr->code eq "NumStops") {
      $num_stops = $attr->value;
    }
  }

  if ($num_stops == 0) {
    return $tran;
  } elsif ($num_stops > $max_stops) {
    return 0;
  }

  ##################
  # surgery; split transcript across stops
  ##################

  while($pep =~ /\*/g) {
    my $position = pos($pep);

    my @coords = $tran->pep2genomic($position, $position);

    if (@coords > 1) {
      # the codon is split by an intron. Messy. Leave these for now
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
  my ($self, $transcripts, $gene_scaf, $map) = @_;

  my (@t_objs, @kept_t_objs);

  TRAN:foreach my $tran (@$transcripts) {
    # reject if transcript has an identical structure to one already seen

    my (@non_gap_exons, @non_gap_introns);
    my $non_gap_length = 0;
    my $gap_length = 0;

    my $last_exon;
    foreach my $e (sort {$a->start <=> $b->start} @{$tran->get_all_Exons}) {      
      my ($coord) = $map->map_coordinates($gene_scaf->seq_region_name,
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
# CHAIN / NET processing functions
###################################################################


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
  my ($self, $chains, $check) = @_;

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



###################################################################
# FUNCTION: remove_irrelevant_chains
#
# Decription:
#  This function returns the chains that have a block overlap
#  with at least one feature in the given list (thus removing
#  "irrelevant" chains that do not interact with exons)
###################################################################
sub remove_irrelevant_chains {
  my ($self, $chains, $feats) = @_;

  my @kept_chains;

  foreach my $c (@$chains) {
    my $chain_overlaps_feat = 0;

    my $j=0;
    BLOCK: foreach my $b (@$c) {
      my $b_st = $b->reference_genomic_align->dnafrag_start;
      my $b_en = $b->reference_genomic_align->dnafrag_end;

      for(; $j < @$feats; $j++) {
        my $this_f = $feats->[$j];

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
# FUNCTION: remove_used_chains
#
# Decription:
#
#  This function filters out chains used so far
###################################################################
sub remove_used_chains {
  my ($self, $raw_chains, $used_so_far) = @_;

  my @chains_to_process;

  foreach my $chain (@$raw_chains) {
    # if this chain conatains a block that overlaps with any of the 
    # excluded regions on the target, ignore it.
    my $keep_chain = 1;

    CHAIN: foreach my $block (@$chain) {
      my ($tga) = @{$block->get_all_non_reference_genomic_aligns};
      foreach my $ex_block (values %$used_so_far) {
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
# FUNCTION: remove_interfering_chains
#
# Decription:
#
# Progressively projects chains on to the query sequence (in order
# of score), rejecting chains that "interfere" with the set of 
# chains retained so far. Interference is defined as an overlap
# in the query at the block level
###################################################################
sub remove_interfering_chains {
  my ($self, $input_chains) = @_;

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

      my $total_chain_cov = 0;

      for(my $i=0; $i < @these_blocks; $i++) {
        my $b = $these_blocks[$i];
        my $ga = $b->reference_genomic_align;

        $total_chain_cov += $ga->dnafrag_end - $ga->dnafrag_start + 1;
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
      if ($longest_len / $total_chain_cov > $self->OVERLAP_CHAIN_FILTER) {
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

    my $blocks = $self->flatten_chains(\@these_chains, 1);

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
            $self->separate_coords($prev, 
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


###################################################################
# FUNCTION: remove_pseudogene_chains
#
# Decription:
#  This function removes chains that look like processed pseudogenes
#  in the target. These can be spotted by having a large extent in
#  the query sequence but a relatively small extent in the target
#  sequence (which gaps between target blocks < 10bp in all cases).
#  A conservative strategy is adopted whereby every block interval
#  in the query has to be small or repeat-filled, and separations
#  in the target are not considered if they are repeat filled
###################################################################
sub remove_pseudogene_chains {
  my ($self, $chains, $min_sep, $max_repeat) = @_;

  my @kept_chains;


  foreach my $c (@$chains) {
    my ($last_b);

    my $real_separations = 0;
    my $total_fusions = 0;
    my $t_name;

    foreach my $b (@$c) {

      if (defined $last_b) {
        my $last_q = $last_b->reference_genomic_align;
        my ($last_t) = @{$last_b->get_all_non_reference_genomic_aligns};

        my $this_q = $b->reference_genomic_align;
        my ($this_t) = @{$b->get_all_non_reference_genomic_aligns};

        $t_name = $this_t->dnafrag->name;

        my $q_interval_start =  $last_q->dnafrag_end + 1;
        my $q_interval_end   =  $this_q->dnafrag_start - 1;
        my $q_interval_length = $q_interval_end - $q_interval_start + 1;

        my ($t_interval_start, $t_interval_end);
        if ($last_t->dnafrag_strand < 0) {
          $t_interval_start = $this_t->dnafrag_end   + 1;
          $t_interval_end   = $last_t->dnafrag_start - 1;
        } else {
          $t_interval_start = $last_t->dnafrag_end   + 1;
          $t_interval_end   = $this_t->dnafrag_start - 1;
        }
        my $t_interval_length = $t_interval_end - $t_interval_start + 1;

        # is this a "real" block separation? We define a real separation
        # as one in which the gap on both sides is large, or if the gap
        # on one side is large and that gap is not occupied completely
        # by repeat

        if ($q_interval_length >= $min_sep and
            $t_interval_length >= $min_sep) {
          $real_separations++;
        } elsif ($q_interval_length >= $min_sep and 
            $t_interval_length < $min_sep) {

          # calculate how much of the interval contains repeats
          my $overlap = 0;
          foreach my $rf (@{$self->query_repeats}) {
            if ($rf->end >= $q_interval_start and 
                $rf->start <= $q_interval_end) {
              # overlap
              my $ov_start = $rf->start; 
              $ov_start = $q_interval_start if $q_interval_start > $ov_start;

              my $ov_end = $rf->end; 
              $ov_end = $q_interval_end if $q_interval_end < $ov_end;
 
              $overlap += $ov_end - $ov_start + 1;
            }
          }

          if ($overlap / $q_interval_length < $max_repeat) {
            $real_separations++;
            $total_fusions++;
          } else {
            # seems to be an artefactual separation caused by repeat
            # insertion since divergence. Ignore this block pair
          }
        }
      }
      $last_b = $b;
    }

    # if there are real separations, and they are ALL fusions in the 
    # target, then we assert that this is a pseudogene chain
    if ($real_separations == 0 or $total_fusions < $real_separations) {
      push @kept_chains, $c;
    }

  }

  return \@kept_chains;
}


#############################
# coordinate utility methods
#############################


sub check_consistent_coords {
  my ($self, $cj, $ci) = @_; 

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
  my ($self, $cj, $ci) = @_; 

  if ($ci->strand == 1) {
    # j should come before i in genomic coords
    return $ci->start - $cj->end - 1;
  } else {
    return $cj->start - $ci->end - 1;
  }
}


sub separate_coords {
  my ($self, $ci, $cj, $target_slice) = @_;

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
# write methods
###########################################################

sub write_agp {
  my ($self, 
      $fh, 
      $gene_scaf, 
      $map) = @_;

  my $prefix = "##-AGP";
  
  print $fh "$prefix \#\n";
  printf($fh "$prefix \# AGP for gene scaffold %s\n", 
         $gene_scaf->seq_region_name);
  print $fh "$prefix \#\n";

  my @pieces = $map->map_coordinates($gene_scaf->seq_region_name,
                                     1,
                                     $gene_scaf->length,
                                     1,
                                     $gene_scaf->seq_region_name);
  my $last_end = 0;
  for(my $i=0; $i < @pieces; $i++) {
    my $piece = $pieces[$i];

    printf($fh "$prefix %s\t%d\t%d\t",
           $gene_scaf->seq_region_name,
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
      $fh,
      $gene_scaf,
      $gene_name,
      @trans) = @_;


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



sub print_chains {
  my ($self, $chains, $text) = @_;

  print "\n$text\n";
  foreach my $c (@$chains) {
    printf "CHAIN %d\n", $c->[0]->score;
    foreach my $b (@$c) {
      my $q = $b->reference_genomic_align;
      my ($t) = @{$b->get_all_non_reference_genomic_aligns};
      printf " %s/%d-%d, %s/%d-%d %d (%d)\n", $q->dnafrag->name, $q->dnafrag_start, $q->dnafrag_end, $t->dnafrag->name, $t->dnafrag_start, $t->dnafrag_end, $t->dnafrag_strand, $t->level_id;
    }
  }

}



1;
