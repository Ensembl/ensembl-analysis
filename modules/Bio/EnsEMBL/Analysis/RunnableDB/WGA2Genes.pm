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

use IO::String;

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB);


my $GENE_SCAFFOLD_INTER_PIECE_PADDING = 100;
my $MAX_EXON_READTHROUGH_DIST = 15;

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

  my (@blocks_used, @results);
  my @alignment_chains = @{$self->genomic_align_block_chains};

  my @cds_feats = @{$self->get_all_transcript_cds_features(@{$self->get_all_Transcripts($self->gene)})};

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
        $target_to_genescaf_map,
        $warnings) = 
            $self->gene_scaffold_from_projection($self->gene->stable_id . "-" . $iteration, 
                                                 \@projected_cds_feats);

    my $result = {
      gene_scaffold         => $gene_scaffold,
      mapper                => $target_to_genescaf_map,
      warnings              => $warnings,
      transcripts           => [],
    };

    $self->write_agp(\*STDOUT, $gene_scaffold, $target_to_genescaf_map, $warnings);

    my @transcripts;
    foreach my $tran (@{$self->get_all_Transcripts($self->gene)}) {

      my $proj_trans = $self->make_projected_transcript($tran,
                                                        $gene_scaffold,
                                                        $query_to_genescaf_map,
                                                        $target_to_genescaf_map);
      if (defined $proj_trans) {
        # for the first, top level net, we split the transcript across stops. 
        # For subsequent Nets, we reject if the transcript contains stops
        $proj_trans = $self->process_transcript($gene_scaffold, 
                                                $proj_trans, 
                                                0);
        if (defined $proj_trans) {
          push @transcripts, $proj_trans;
        }
      }
    }
    
    if (@transcripts) {
      push @{$result->{transcripts}}, @{$self->make_nr_transcript_set(\@transcripts, 
                                                                      $gene_scaffold,
                                                                      $target_to_genescaf_map)};
    }
    
    push @results, $result;
  }

  $self->output(\@results);
}


sub write_output {
  my ($self) = @_;

  my $gene = $self->gene;
  print  "#\n";
  printf("# WGA2Genes output for gene %s (%s/%d-%d %s)\n", 
         $gene->stable_id, 
         $gene->slice->seq_region_name,
         $gene->start,
         $gene->end,
         $gene->strand);
  print "#\n";

  foreach my $obj (@{$self->output}) {
    my $gs = $obj->{gene_scaffold};
    my $map = $obj->{mapper};
    my $warns = $obj->{warnings};
    my @transcripts = @{$obj->{transcripts}};

    # determine whether there is any output for this gene scaffold
    if (@transcripts > 0) {

      $self->write_agp(\*STDOUT, $gs, $map, $warns);
      $self->write_gene(\*STDOUT, $gs, @transcripts);
    }
  }

  # do nothing

  return;
}


######################################
# internal methods
#####################################



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

  my @ret_trans;
  foreach my $t (@{$runnable->output}) {
    foreach my $e (@{$t->get_all_Exons}) {
      $e->slice($genomic);
    }
    push @ret_trans, $t;
  }

  @ret_trans = @{$self->make_consistent_transcript_set(\@ret_trans)};

  return \@ret_trans;
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
  my ($self, $gene_scaffold_id, $projected_cds_elements) = @_;

  my @projected_cds_elements = @$projected_cds_elements;

  my (@targets, @new_targets);

  # step 0: flatten list
  my $cds_num = $self->gene->strand > 0 ? 1 : scalar(@projected_cds_elements);
  foreach my $cds (@projected_cds_elements) {
    foreach my $coord_pair (@$cds) {
      my $coord = $coord_pair->{target};

      push @targets, {
        cds_id  => $cds_num,
        coord   => $coord,
      };
    }
    $cds_num += $self->gene->strand;
  }


  # step 1 remove gaps which cannot be filled
  @new_targets = ();
  for(my $i=0; $i < @targets; $i++) {
    my $this_target = $targets[$i];

    if ($this_target->{coord}->isa("Bio::EnsEMBL::Mapper::Gap")) {
      # if it's gap that can be filled, leave it. Otherwise, remove it
      my ($left_target, $right_target);
      for(my $j=$i-1; $j>=0; $j--) {
        if ($targets[$j]->{coord}->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
          $left_target = $targets[$j];
          last;
        }
      }
      for(my $j=$i+1; $j < @targets; $j++) {
        if ($targets[$j]->{coord}->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
          $right_target = $targets[$j];
          last;
        }
      }

      # under what circumstances can this gap NOT be filled?
      # 1. If one of the flanking coords is on the same CDS region
      #    as the gap, and the end of the aligned region does
      #    not align to a gap
      # 2. If the 2 flanking coords are consistent and not
      #    separated by a sequence-level gap
      my $keep_gap = 1;

      if (defined $left_target and 
          defined $right_target and
          $self->check_consistent_coords($left_target->{coord}, 
                                         $right_target->{coord})) {
        
        my ($left_coord, $right_coord) = ($left_target->{coord}, 
                                          $right_target->{coord});
        if ($left_coord->strand < 0) {
          ($left_coord, $right_coord) = ($right_coord, $left_coord);
        }
        my ($new_left, $new_right) = 
            $self->separate_coords($left_coord, $right_coord,
                                   $self->target_slices->{$left_coord->id});
        
        if (not defined $new_left and not defined $new_right) {
          $keep_gap = 0;
        }
      } 


      if ($keep_gap and 
          defined $left_target and 
          $left_target->{cds_id} == $this_target->{cds_id}) {
        
        my ($new_left, $dummy);
        if ($left_target->{coord}->strand < 0) {
          ($dummy, $new_left) = 
              $self->separate_coords(undef,
                                     $left_target->{coord}, 
                                     $self->target_slices->{$left_target->{coord}->id});

          if (abs($new_left->start - $left_target->{coord}->start) > 10) {
            # we're not near the end of a contig. This gap cannot be filled;
            $keep_gap = 0;
          }
        } else {
          ($new_left, $dummy) = 
              $self->separate_coords($left_target->{coord}, 
                                     undef,
                                     $self->target_slices->{$left_target->{coord}->id});

          if (abs($new_left->end - $left_target->{coord}->end) > 10) {
            # we're not near the end of a contig. This gap cannot be filled;
            $keep_gap = 0;
          }
        }
      }

      if ($keep_gap and 
          defined $right_target and
          $right_target->{cds_id} == $this_target->{cds_id}) {

        my ($dummy, $new_right);
        if ($right_target->{coord}->{strand} < 0) {
          ($new_right, $dummy) = 
              $self->separate_coords($right_target->{coord}, 
                                     undef,
                                     $self->target_slices->{$right_target->{coord}->id});

          if (abs($right_target->{coord}->end - $new_right->end) > 10) {
            # we're not near the end of a contig. This gap cannot be filled;
            $keep_gap = 0;
          }
        } else {
          ($dummy, $new_right) = 
              $self->separate_coords(undef,
                                     $right_target->{coord}, 
                                     $self->target_slices->{$right_target->{coord}->id});

          if (abs($right_target->{coord}->start - $new_right->start) > 10) {
            # we're not near the end of a contig. This gap cannot be filled;
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

  # step 2: prune away gaps at the start and end; we won't "fill" these
  while(@targets and $targets[0]->{coord}->isa("Bio::EnsEMBL::Mapper::Gap")) {
    shift @targets;
  }
  while(@targets and $targets[-1]->{coord}->isa("Bio::EnsEMBL::Mapper::Gap")) {
    pop @targets;
  }

  # step 3: merge adjacent targets   
  #  we want to be able to account for small, frame-preserving insertions
  #  in the target sequence with respect to the query. To give the later,
  #  gene-projection code the opportunity to "read through" these insertions,
  #  we have to merge togther adjacent, consistent targets that are within
  #  this "maximum read-through" distance
  @new_targets = ();
  for(my $i=0; $i<@targets; $i++) {
    my $target = $targets[$i];

    if ($target->{coord}->isa("Bio::EnsEMBL::Mapper::Coordinate") and
        @new_targets and
        $new_targets[-1]->{coord}->isa("Bio::EnsEMBL::Mapper::Coordinate") and 
        $self->check_consistent_coords($new_targets[-1]->{coord}, 
                                       $target->{coord})) {
      
      my $dist = $self->distance_between_coords($new_targets[-1]->{coord}, 
                                                $target->{coord});
      
      if ($dist <= $MAX_EXON_READTHROUGH_DIST) {
        
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
                $tg_c->id eq $new_coord->id and$tg_c->start < $new_coord->end and
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



  my $warn_strings = {
    1 => "inconsistent strand",
    2 => "inconsistent coordinate order",
    3 => "split across single sequence-level contig",
  };
  my (%warnings, @warnings);

  for(my $i=0; $i < @targets; $i++) {
    my $this_coord = $targets[$i]->{coord};
    
    next if not $this_coord->isa("Bio::EnsEMBL::Mapper::Coordinate");
    # scan left
    my $seen_something_else = 0;
    for(my $j=$i-1; $j>=0; $j--) {
      my $left_coord = $targets[$j]->{coord};

      if ($left_coord->isa("Bio::EnsEMBL::Mapper::Gap") or
          $left_coord->id ne $this_coord->id) {
        $seen_something_else = 1;
        next;
      }
          
      # assert: both coords, with same id
      if ($left_coord->strand != $this_coord->strand) {
        $warnings{$left_coord->id}->{'1'} = 1;
      } elsif (not $self->check_consistent_coords($left_coord, $this_coord)) {
        $warnings{$left_coord->id}->{'2'} = 1;
      } elsif ($seen_something_else) {
        if ($left_coord->strand < 0) {
          # $left_coord will be downstream of $this_coord in target
          ($left_coord, $this_coord) = ($this_coord, $left_coord);
        }
        my ($new_left, $new_this) = 
            $self->separate_coords($left_coord, $this_coord,
                                   $self->target_slices->{$left_coord->id});
        if (not defined $new_left and not defined $new_this) {
          $warnings{$left_coord->id}->{'3'} = 1;
        }
      }
      
      last;
    }
  }
  foreach my $id (keys %warnings) {
    foreach my $code (keys %{$warnings{$id}}) {
          
      push @warnings, sprintf("%s : %s", $id, $warn_strings->{$code}); 
    }
  }

  ###################################
  # build the gene scaffold, and mappings between the new gene
  # scaffold and each of query and target coords

  my $target_map = Bio::EnsEMBL::Mapper->new('target',
                                             $gene_scaffold_id); 
  my $query_map = Bio::EnsEMBL::Mapper->new('query',
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
      $target_map->add_map_coordinates($coord->id,
                                       $coord->start,
                                       $coord->end,
                                       $coord->strand,
                                       $gene_scaffold_id,
                                       $last_end_pos + 1,
                                       $last_end_pos + $coord->length);
    } else {
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

    # add padding between the pieces
    if ($i < @targets - 1) {
      $last_end_pos += $coord->length + $GENE_SCAFFOLD_INTER_PIECE_PADDING;
      $seq .= ('n' x $GENE_SCAFFOLD_INTER_PIECE_PADDING);
    }
  }
  

  # now add of the original exon pieces to the query map
  foreach my $el (@projected_cds_elements) {
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

  return ($gene_scaffold, $query_map, $target_map, \@warnings);
}




=head2 make_projected_transcript

 Title   : make_projected_transcript
 Description:
    Takes a transcript, a gene scaffold, and a mapping between query coords and
    gene scaffold coords. Produces the transcript that is the result of "projecting"
    the original transcript, through alignment, onto the gene scaffold. 
=cut

sub make_projected_transcript {
  my ($self, $transcript, $gene_scaf, $qmapper, $tmapper) = @_;

  my (@all_coords, @new_exons);

  my @orig_exons = @{$transcript->get_all_translateable_Exons};
  if ($transcript->strand < 0) {
    @orig_exons = reverse @orig_exons;
  }

  foreach my $orig_exon (@orig_exons) {    
    my @these_coords = $qmapper->map_coordinates($orig_exon->slice->seq_region_name,
                                                 $orig_exon->start,
                                                 $orig_exon->end,
                                                 1,
                                                 "query");
    push @all_coords, @these_coords;
  }

  #
  # replace coords at start and end that map down to gaps with gaps
  # 
  for(my $i=0; $i < @all_coords; $i++) {
    if ($all_coords[$i]->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
      my ($tcoord) = $tmapper->map_coordinates($gene_scaf->seq_region_name,
                                               $all_coords[$i]->start,
                                               $all_coords[$i]->end,
                                               1,
                                               $gene_scaf->seq_region_name);
      if ($tcoord->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
        last;
      } else {
        $all_coords[$i] = Bio::EnsEMBL::Mapper::Gap->new(1, $tcoord->length);
      }
    }
  }
  for(my $i=scalar(@all_coords)-1; $i >= 0; $i--) {
    if ($all_coords[$i]->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
      my ($tcoord) = $tmapper->map_coordinates($gene_scaf->seq_region_name,
                                               $all_coords[$i]->start,
                                               $all_coords[$i]->end,
                                               1,
                                               $gene_scaf->seq_region_name);
      if ($tcoord->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
        last;
      } else {
        $all_coords[$i] = Bio::EnsEMBL::Mapper::Gap->new(1, $tcoord->length);
      }
    }
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

        # calculate "surplus" bases on incomplete codons to left and right
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

  my ($total_transcript_bps, $real_seq_bps);
  foreach my $coord (@all_coords) {
    if ($coord->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
      push @new_exons, Bio::EnsEMBL::Exon->new(-start     => $coord->start,
                                               -end       => $coord->end,
                                               -strand    => $transcript->strand,
                                               -slice     => $gene_scaf);

      $total_transcript_bps += $coord->length;
      my ($tcoord) = $tmapper->map_coordinates($gene_scaf->seq_region_name,
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
    return undef;
  }

  #
  # sort exons into rank order 
  #
  if ($transcript->strand < 0) {
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

      my @gen_coords = $qmapper->map_coordinates($gene_scaf->seq_region_name,
                                                 $extent_start,
                                                 $extent_end,
                                                 1,
                                                 $gene_scaf->seq_region_name);

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
                                                  -hseqname => $transcript->translation->stable_id,
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

  #
  # merge abutting exons; deals with exon fusion events, and 
  # small, frame-preserving insertions in the target
  #
  my @merged_exons;
  foreach my $exon (@new_exons) {
    if (@merged_exons) {

      my $prev_exon = pop @merged_exons;
   
      my ($new_start, $new_end);

      if ($transcript->strand < 0) {
        my $intron_len = $prev_exon->start - $exon->end - 1;
        if ($intron_len % 3 == 0 and $intron_len <= $MAX_EXON_READTHROUGH_DIST) { 
          $new_start = $exon->start;
          $new_end   = $prev_exon->end;
        }
      } else {
        my $intron_len = $exon->start - $prev_exon->end - 1;
        if ($intron_len % 3 == 0 and $intron_len <= $MAX_EXON_READTHROUGH_DIST) {
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
  
  #
  # do transcript-level supporting features/attributes
  #
  my $t_sf = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@trans_fps);
  $t_sf->hcoverage(100 * ($hcoverage / $transcript->translate->length));
  # misuse score field as the proportion of the transcript that maps to non-gap sequence
  $t_sf->score( 100 * ($real_seq_bps / $total_transcript_bps) ); 
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


  return $proj_tran;
}




=head2 process_transcript

 Title   : process_transcript
 Description:

=cut

sub process_transcript {
  my ($self,$gene_scaffold, $tran, $splice_out_stops) = @_;
  
  my @processed_transcripts;

  my @exons = @{$tran->get_all_Exons};
  my ($tsf) = @{$tran->get_all_supporting_features};

  ##################
  # reject transcripts that have:
  #   less than minimum coverage of parent peptide
  #   higher that maximum proportion of gap residues
  ##################
  return undef if $tsf->hcoverage < $self->MIN_COVERAGE;
  return undef if $tsf->score < $self->MIN_NON_GAP;

  my $pep = $tran->translate->seq;

  if ($pep !~ /\*/) {
    return $tran;
  } elsif (not $splice_out_stops) {
    return undef;
  }

  ##################
  # surgery; split transcript across stops
  ##################

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
  
  return $tran;
}




=head2 make_nr_transcript_set

 Title   : make_nr_transcript_set
    Takes an initial "raw", ordered set of transcripts and proceeds through
    the list, rejecting transcripts that have no "unique" introns with
    respect to the previous one. For this, "gap" exons are ignored.

=cut

sub make_nr_transcript_set {
  my ($self, $transcripts, $gene_scaf, $map) = @_;

  my (@all_introns, @kept_transcripts);

  foreach my $tran (@$transcripts) {
    my (@exons, @introns);
    foreach my $exon (sort {$a->start <=> $b->start} @{$tran->get_all_Exons}) {
      my ($coord) = $map->map_coordinates($gene_scaf->seq_region_name,
                                          $exon->start,
                                          $exon->end,
                                          1,
                                          $gene_scaf->seq_region_name);

      if ($coord->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
        push @exons, $exon;
      }
    }
    
    # get introns
    for(my$i=1; $i < @exons; $i++) {
      push @introns, Bio::EnsEMBL::Feature->new(-start => $exons[$i-1]->end + 1,
                                                -end   => $exons[$i]->start - 1);
    }

    # if none of the introns are unique, reject. 
    my $at_least_one_unique = 0;
    foreach my $this_intron (@introns) {
      my $unique_intron = 1;
      foreach my $other_intron (@all_introns) {
        if ($this_intron->start == $other_intron->start and
            $this_intron->end   == $other_intron->end) {
          $unique_intron = 0;
          last;
        }
      }
      if ($unique_intron) {
        $at_least_one_unique = 1;
        last;
      }
    }
    if ($at_least_one_unique) {
      push @all_introns, @introns;
      push @kept_transcripts, $tran;
    }
  }

  return \@kept_transcripts;
}


#################
# checks
#################


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
# write methods
###########################################################

sub write_agp {
  my ($self, 
      $fh, 
      $gene_scaf, 
      $map,
      $warns) = @_;

  my $prefix = "##-AGP";
  
  print $fh "$prefix \#\n";
  printf $fh "$prefix \# AGP for gene scaffold %s\n", $gene_scaf->seq_region_name;
  foreach my $warning (@$warns) {
    print $fh "$prefix \# WARNING : $warning\n";
  }
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
      @trans) = @_;


  my $prefix = "##-GENES ";
  
  my $seq_id = $gene_scaf->seq_region_name;
  my $gene_id = $seq_id;
  
  my $fasta_string;
  my $stringio = IO::String->new($fasta_string);
  my $seqio = Bio::SeqIO->new(-format => 'fasta',
                              -fh => $stringio);

  printf $fh "$prefix \# Gene report for $seq_id\n";
  
  foreach my $tran (@trans) {
    my ($sf) = @{$tran->get_all_supporting_features};
    my $tran_id = $gene_id . "_" . $sf->hseqname; 
 
    print $fh "$prefix \# Transcript $tran_id:\n";
    printf($fh "$prefix \# Coverage = %.2f\n", $sf->hcoverage);
    printf($fh "$prefix \# Proportion non-gap = %.2f\n", $sf->score);
   
    my $pep = Bio::PrimarySeq->new(-id => $tran_id,
                                   -seq => $tran->translate->seq);
    $seqio->write_seq($pep);
            
    my @exons = @{$tran->get_all_Exons};
    
    for (my $i=0; $i < @exons; $i++) {
      my $exon = $exons[$i];
      my $exon_id = $tran_id . "_" . $i;
      
      printf($fh "$prefix %s\t%s\t%d\t%d\t%d\t%d\t%d\texon=%s; transcript=%s; gene=%s\n", 
             $seq_id,
             "Exon",
             $exon->start,
             $exon->end,
             $exon->strand,
             $exon->phase,
             $exon->end_phase,
             $exon_id,
             $tran_id,
             $gene_id);
      foreach my $sup_feat (@{$exon->get_all_supporting_features}) {
        foreach my $ug ($sup_feat->ungapped_features) {
          printf($fh "$prefix %s\t%s\t%d\t%d\t%d\texon=%s; hname=%s; hstart=%s; hend=%s\n", 
                 $seq_id,
                 "Supporting",
                 $ug->start,
                 $ug->end,
                 $ug->strand,
                 $exon_id,
                 $ug->hseqname,
                 $ug->hstart,
                 $ug->hend);
        }
      }      
    }
  }

  $seqio->close;
  foreach my $line (split /\s+/, $fasta_string) {
    print $fh "##-FASTA $line\n";
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

  # if the single flag is given returns the longest
  # coding transcript. Otherwise returns all coding
  # transcripts. 

  if (not exists $self->{_transcripts}) {
    $self->{_transcripts} = {};
  }
  
  if (not exists $self->{_transcripts}->{$gene}) {

    my (@transcripts);
    
    foreach my $t (@{$gene->get_all_Transcripts}) {
      my $translation = $t->translate;
      
      if (defined $translation) {
        push @transcripts, { 
          tran   => $t,
          length => $translation->length,
        };
      }
    }
    
    @transcripts = sort { $b->{length} <=> $a->{length} } @transcripts;

    $self->{_transcripts}->{$gene} = [map {$_->{tran}} @transcripts];
  } 
  
  if ($single) {
    [$self->{_transcripts}->{$gene}->[0]];
  } else {
    return $self->{_transcripts}->{$gene};
  }
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



1;
