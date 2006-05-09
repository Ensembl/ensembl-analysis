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
  $genscan->write_output(); #writes to stdout


=head1 DESCRIPTION


=cut
package Bio::EnsEMBL::Analysis::RunnableDB::WGA2Genes;

use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Analysis::RunnableDB;

use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Analysis::Config::WGA2Genes;

use Bio::EnsEMBL::Analysis::Tools::Logger;
use Bio::EnsEMBL::Analysis::Tools::WGA2Genes::GeneScaffold;
use Bio::EnsEMBL::Analysis::Tools::WGA2Genes::CoordUtils;
use Bio::EnsEMBL::Analysis::Tools::WGA2Genes::ChainUtils;



@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB);


############################################################
=head2 new
    Title   :   fetch_input
=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  $self->read_and_check_config($WGA2GENES_CONFIG_BY_LOGIC);

  return $self;
}

############################################################
=head2 fetch_input
    Title   :   fetch_input
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

  ########
  # check that the default assembly for the query and target agrees 
  # with that for the method_link_species_set GenomeDBs
  ########

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

  ##############
  # populate a list of genes and transcripts to ignore
  ##############
  my $kill_list = {};
  if ($self->KILL_LIST) {
    open(KILL, $self->KILL_LIST) or 
        throw("Could not open kill list " . $self->KILL_LIST . " for reading");
    while(<KILL>) {
      /^\#/ and next;
      /^(\S+)/ and $kill_list->{$1} = 1;
    }
    close(KILL);
  }

  ########
  # fetch the genes; need to work in the coordinate space of the 
  # top-level slice to be consistent with compara
  ########
  my $sa = $q_dbh->get_SliceAdaptor;

  my (@gene_records, $reg_start, $reg_end);

  if ($input_id =~ /:/) {
    # assume slice name

    my $slice = $sa->fetch_by_name($input_id);
    $reg_start = $slice->start;
    $reg_end   = $slice->end;
    
    $self->query_slice($sa->fetch_by_region('toplevel',
                                            $slice->seq_region_name));

    foreach my $g (@{$slice->get_all_Genes}) {
      next if $g->biotype ne 'protein_coding';
      next if exists $kill_list->{$g->stable_id};

      $g = $g->transfer($self->query_slice);

      my @good_trans;

      foreach my $t (@{$g->get_all_Transcripts}) {
        next if $t->coding_region_start < $slice->start;
        next if $t->coding_region_end   > $slice->end;
        next if exists $kill_list->{$t->stable_id};
       
        push @good_trans, $t;
      }

      if (@good_trans) {
        if ($self->LONGEST_SOURCE_TRANSCRIPT) {
          @good_trans = sort { 
            length($b->translateable_seq) <=> length($a->translateable_seq);
          } @good_trans;
          @good_trans = ($good_trans[0]);
        }
        push @gene_records, Bio::EnsEMBL::Analysis::RunnableDB::WGA2Genes::GeneRecord
            ->new(-gene => $g,
                  -source_transcripts => \@good_trans);
      }
    }

  } else {
    # assume iid is a gene stable id; ignore the kill-list
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

    $reg_start = $gene->start;
    $reg_end   = $gene->end;

    my @good_trans;
    foreach my $t (@{$gene->get_all_Transcripts}) {
      if (not exists($kill_list->{$t->stable_id})) {
        push @good_trans, $t;
      }
    }

    if (@good_trans) {
      @good_trans = sort { 
        length($b->translateable_seq) <=> length($a->translateable_seq);
      } @good_trans;
      @good_trans = ($good_trans[0]);

      push @gene_records, Bio::EnsEMBL::Analysis::RunnableDB::WGA2Genes::GeneRecord
          ->new(-gene => $gene,
                -source_transcripts => \@good_trans);
    }
  }

  $self->gene_records( \@gene_records );
  

  #########
  # get the compara data: MethodLinkSpeciesSet, reference DnaFrag, 
  # and all GenomicAlignBlocks
  #########
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


############################################################
=head2 run
    Title   :   run
=cut

sub run {
  my ($self) = @_;

  my (@results, @original_genes, @chains_generecs_pairs);

  logger_info("INITIAL GENES: " . 
              join(" ", map {$_->gene->stable_id} @{$self->gene_records}) . "\n");

  logger_info("ORIGINAL_CHAINS\n" . 
              stringify_chains($self->genomic_align_block_chains));

  @chains_generecs_pairs = 
      $self->segregate_chains_and_generecords($self->genomic_align_block_chains,
                                              $self->gene_records);
      
  foreach my $chains_generecs (@chains_generecs_pairs) {
    my ($alignment_chains, $generecs) = @$chains_generecs;
    
    my $pair_id = join(":", 
                       map { $_->gene->stable_id } @$generecs);
    
    logger_info("CHAINS_FOR_GENESET $pair_id\n" . 
                stringify_chains($alignment_chains));
    
    my %blocks_from_used_chains;
    
    ITER:for(my $iteration=0; $iteration < 1 ;$iteration++) {

      logger_info("** Iteration $iteration $pair_id\n");
    
      my $filtered_chains = filter_block_overlap_chains($alignment_chains,
                                                        [values %blocks_from_used_chains]);

      logger_info("NON_USED_CHAINS\n" . stringify_chains($filtered_chains));
      
      my @these_generecs = map { $_->clone } @$generecs;

      # Repeat transcript contruction until we had no rejected transcripts.
      # This is because if a transcript is rejected, some of the underlying 
      # GeneScaffold components may be unnecessary
      for(;;) {
        logger_info("Working with genes " . join(" ", map { $_->gene->stable_id} @these_generecs));

        my @cds_feats = map { @{$_->transcript_cds_feats} } @these_generecs;

        $filtered_chains =             
            filter_irrelevant_chains($filtered_chains,
                                     \@cds_feats);

        logger_info("RELEVANT_CHAINS\n" . stringify_chains($filtered_chains));
        
        $filtered_chains = filter_inconsistent_chains($filtered_chains, 
                                                      $self->OVERLAP_CHAIN_FILTER);
        logger_info("CONSISTENT_CHAINS\n" . stringify_chains($filtered_chains));
        
        $filtered_chains = $self->remove_contig_split_chains($filtered_chains);
        logger_info("NON_CONTIG_SPLIT_CHAINS\n" . stringify_chains($filtered_chains));
        
        my @these_pairs = $self->segregate_chains_and_generecords($filtered_chains,
                                                                  \@these_generecs);
      
        my @these_results;
        
        foreach my $pair (@these_pairs) {
          my ($subset_chains, $subset_generecs) = @$pair;
          
          my $net_blocks = flatten_chains($subset_chains, 1);
          
          my $gs_name = $subset_generecs->[0]->gene->stable_id . "-" . $iteration;
          
          my $gene_scaffold = $self->make_gene_scaffold_and_project_genes($net_blocks,
                                                                          $subset_generecs,
                                                                          $gs_name);
          
          push @these_results, [$gene_scaffold, @$subset_generecs];
        }
        
        # check that all transcripts were successfully projected. If not,
        # the gene scaffolds may contain unnecessary pieces, in which case
        # we remove the unsuccesful transcript and repeat. 
        
        my @kept_generecs;
        my $need_to_repeat = 0;

        foreach my $res (@these_results) {
          my ($gs, @res_genes) = @$res;
          
          foreach my $res_gene (@res_genes) {
            if (scalar(@{$res_gene->source_transcripts}) == 
                scalar(@{$res_gene->projected_transcripts})) {
              # all okay
            } else {
              $res_gene->source_transcripts($res_gene->good_source_transcripts);
              $res_gene->projected_transcripts([]);
              $res_gene->good_source_transcripts([]);
              
              if (@{$res_gene->source_transcripts}) {
                push @kept_generecs, $res_gene;
              }
              
              $need_to_repeat = 1;
            }
          }
        }
        
        if ($need_to_repeat) {
          if (@kept_generecs) {
            @these_generecs = @kept_generecs;
            next;
          } else {
            last;
          }
        } else { 
          if (@these_results) {
            push @results, @these_results; 
            
            foreach my $pair (@these_pairs) {
              my ($these_chains) = @$pair;
              foreach my $chain (@$these_chains) {
                foreach my $block (@$chain) {
                  $blocks_from_used_chains{$block} = $block;
                }
              }
            }
          }
          
          last;
        }
      }
    }
  }

  $self->output(\@results);
}


############################################################
=head2 write_output
    Title   :   write_output
=cut

sub write_output {
  my ($self) = @_;

  print  "#\n";
  printf("# WGA2Genes output for %s\n", $self->input_id);
  printf("#  Genes are:");
  foreach my $grec (@{$self->gene_records}) {
    print " " . $grec->gene->stable_id;
  }
  print "\n#\n";

  foreach my $obj (@{$self->output}) {
    my ($gs, @res_genes) = @$obj;

    $self->write_agp($gs);
    foreach my $g (@res_genes) {
      $self->write_gene($gs, 
                        $g);
    }
  }

  return;
}


######################################
# internal methods
#####################################



#################################################################
# FUNCTION  : make_gene_scaffold_and_project_genes
#
# Description:
#  Obvious
#################################################################
sub make_gene_scaffold_and_project_genes {
  my ($self, $blocks, $res_genes, $name) = @_;
  
  my $gene_scaffold = Bio::EnsEMBL::Analysis::Tools::WGA2Genes::GeneScaffold->
      new(-name => $name,
          -genomic_align_blocks => $blocks,
          -from_slice => $self->query_slice,
          -to_slices  => $self->target_slices,
          -transcripts   => [map {@{$_->source_transcripts}} @$res_genes],
          -max_readthrough_dist => $self->MAX_EXON_READTHROUGH_DIST,
          -add_gaps => 1,
          );
  
  foreach my $res_gene (@$res_genes) {
    $res_gene->name($name . "." . $res_gene->gene->stable_id);
    
    foreach my $tran (@{$res_gene->source_transcripts}) {
      my $proj_trans = 
          $gene_scaffold->place_transcript($tran);
      
      $proj_trans = 
          $self->process_transcript($proj_trans, 
                                    $self->MAX_EDITABLE_STOPS_PRIMARY,
                                    $self->MIN_COVERAGE,
                                    $self->MIN_NON_GAP,
                                    $tran->stable_id);
      
      if ($proj_trans) {
        push @{$res_gene->good_source_transcripts}, $tran;
        push @{$res_gene->projected_transcripts}, $proj_trans;
      }
    }
  }

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
# FUNCTION: remove_contig_split_chains
#
# Decription:
#  This function removes chains that (a) interlock with higher
#  scoring chains, and (b) this interlocking cannot be explained
#  by contig mis-ordering within the involved scaffolds
###################################################################
sub remove_contig_split_chains {
  my ($self, $input_chains) = @_;

  my @kept_chains;

  foreach my $c (@$input_chains) {
    my (%coords_by_tid, @all_coords, @merged_coords);

    my $consistent = 1;

    my @these_chains = (@kept_chains, $c);

    my $blocks = flatten_chains(\@these_chains, 1);

    foreach my $bl (@$blocks) {
      my ($tga) = @{$bl->get_all_non_reference_genomic_aligns};
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
      my @chain = @{$coords_by_tid{$id}};
      @chain = sort { $a->start <=> $b->start } @chain; 

      # each chain should be separable from the last at contig
      # level
      for(my $i=1; $i < @chain; $i++) {
        my $prev = $chain[$i-1];
        my $this = $chain[$i];

        my ($ex_prev) = extend_coord($prev,
                                   $self->target_slices->{$id});

        my ($ex_this) = extend_coord($this,
                                   $self->target_slices->{$id});

        if ($ex_prev->end >= $ex_this->start) {
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
# FUNCTION: segregate_chains_and_generecords
#
# Decription:
#   partition the chains and genes into (chain_list, gene_list)
#   pairs, where each chain_list contains chains that touch
#   one or more genes in gene_list, and each gene
###################################################################
sub segregate_chains_and_generecords {
  my ($self, $chains, $generecs) = @_;

  my (%generecs_by_id, %chains_by_id);

  foreach my $grec (@$generecs) {
    $generecs_by_id{$grec} = $grec;
  }

  my (%genes_per_chain);

  foreach my $c (@$chains) {
    $chains_by_id{$c} = $c;

    #print "CONSIDERING CHAIN $c\n";
    
    my @ref_gas = map { $_->reference_genomic_align } @$c;
    @ref_gas = sort { $a->dnafrag_start <=> $b->dnafrag_start } @ref_gas;
    
    my %overlapping;
    
    foreach my $grec (@$generecs) {
      my $overlaps = 0;

    BLOCK: 
      foreach my $bl (@ref_gas) {
        #printf " looking at block %s %d %d\n", $bl->dnafrag->name, $bl->dnafrag_start, $bl->dnafrag_end;

      FEAT:
        foreach my $f (@{$grec->transcript_cds_feats}) {
          #printf "  looking at cds feature %d %d\n", $f->start, $f->end;
          if ($f->start > $bl->dnafrag_end) {
            # no exons in this transcript will overlap the block
            next BLOCK;
          } elsif ($f->end < $bl->dnafrag_start) {
            next FEAT;
          } else {
            $overlaps = 1;
            last BLOCK;
          }
        }
      }
      
      if ($overlaps) {
        $overlapping{$grec} = 1;
      }
    }
    
    foreach my $grecid (keys %overlapping) {
      #print "FOUND A GENE OVERLAP\n";
      $genes_per_chain{$c}->{$grecid} = 1;
    }
  }
    
  my (@clusters, @name_clusters, @res_clusters);
    
  foreach my $cref (keys %genes_per_chain) {
    my $clust = {
      chains => { $cref => 1 },
      genes  => {},
    };
    
    foreach my $grecid (keys %{$genes_per_chain{$cref}}) {
      $clust->{genes}->{$grecid} = 1;
    }

    push @clusters, $clust;
  }
    
  while(@clusters) {
    my $this_clust = shift @clusters;

    my @merged;
    
    do {
      @merged = ();

      for (my $i=0; $i < @clusters; $i++) {
        # possibly merge with this clust
        # if merge, remove from @clusters and re-search
        my $other_clust = $clusters[$i];

        my $in_common = 0;
        foreach my $this_g (keys %{$this_clust->{genes}}) {
          if (exists $other_clust->{genes}->{$this_g}) {
            $in_common = 1;
            last;
          }
        }
        
        if ($in_common) {
          foreach my $cref (keys %{$other_clust->{chains}}) {
            $this_clust->{chains}->{$cref} = 1;
          }
          foreach my $gref (keys %{$other_clust->{genes}}) {
            $this_clust->{genes}->{$gref} = 1;
          }

          push @merged, $i;
        }
      }

      if (@merged) {
        foreach my $i (reverse @merged) {
          splice (@clusters, $i, 1);
        }
      }
    } while @merged;

    push @name_clusters, $this_clust;
  }

  # finally: form the clusters
  foreach my $c (@name_clusters) {
    my (@chains, @genes);

    foreach my $cref (keys %{$c->{chains}}) {
      push @chains, $chains_by_id{$cref};
    }
    foreach my $gref (keys %{$c->{genes}}) {
      push @genes, $generecs_by_id{$gref};
    }

    @chains = sort { $b->[0]->score <=> $a->[0]->score } @chains;
    @genes  = sort { $a->gene->start <=> $b->gene->start } @genes;

    push @res_clusters, [\@chains, \@genes];
  }

  return @res_clusters;
}      


###########################################################
# write methods
###########################################################

sub write_agp {
  my ($self, 
      $gscaf) = @_;

  my $fh = \*STDOUT;

  my $prefix = "##-AGP";

  my @tsegs = $gscaf->project_down;
  my @qsegs = $gscaf->project_up;
  
  print $fh "$prefix \#\n";
  printf($fh "$prefix \#\# AGP for gene scaffold %s source-region=%s/%d-%d\n", 
         $gscaf->seq_region_name,
         $qsegs[0]->to_Slice->seq_region_name,
         $qsegs[0]->to_Slice->start,
         $qsegs[-1]->to_Slice->end);

  print $fh "$prefix \#\n";

  my $piece_count = 1;

  for(my $i=0; $i < @tsegs; $i++) {
    my $seg = $tsegs[$i];

    printf($fh "$prefix %s\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%s\n",
           $gscaf->seq_region_name,
           $seg->from_start,
           $seg->from_end,
           $piece_count++,
           "W",
           $seg->to_Slice->seq_region_name,
           $seg->to_Slice->start,
           $seg->to_Slice->end,
           $seg->to_Slice->strand < 0 ? "-" : "+"
           );

    if ($i < @tsegs - 1) {
      printf($fh "$prefix %s\t%d\t%d\t%s\t%s\t%d\n",
             $gscaf->seq_region_name,
             $seg->from_end + 1,
             $tsegs[$i+1]->from_start - 1,
             $piece_count++,
             "N",
             $tsegs[$i+1]->from_start - $seg->from_end - 1,
             );
    }
  }
}


sub write_gene {
  my ($self, 
      $gene_scaf,
      $g) = @_;

  my $fh = \*STDOUT;

  my $prefix = "##-GENES ";
  
  my $seq_id = $gene_scaf->seq_region_name;
  my $gene_id = $g->name;

  printf $fh "$prefix \# Gene report for $seq_id\n";
  
  foreach my $tran (@{$g->projected_transcripts}) {
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


sub gene_records {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_gene_recs} = $val;
  }

  return $self->{_gene_recs};
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
# chain filtering
#


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


sub LONGEST_SOURCE_TRANSCRIPT {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_longest_source} = $val;
  }

  return $self->{_longest_source};
}


sub KILL_LIST {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_kill_list} = $val;
  }

  return $self->{_kill_list};
}


###########################################
# Local class Result to encapsulate results
############################################

package Bio::EnsEMBL::Analysis::RunnableDB::WGA2Genes::GeneRecord;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );

sub new {
  my ($class, @args) = @_;

  my $self = bless {}, $class;

  my ($name, $gene, $source_trans) = rearrange
    (['NAME', 'GENE', 'SOURCE_TRANSCRIPTS'], @args);

  $self->name($name) if defined $name;
  $self->gene($gene) if defined $gene;
  $self->source_transcripts($source_trans) if defined $source_trans;

  $self->good_source_transcripts([]);
  $self->projected_transcripts([]);

  return $self;
}


sub clone {
  my ($self) = @_;

  my $new_one = Bio::EnsEMBL::Analysis::RunnableDB::WGA2Genes::GeneRecord
      ->new(-source_transcripts => $self->source_transcripts,
            -gene               => $self->gene);

  return $new_one;
}



#################################################################
# FUNCTION   : get_all_transcript_cds_features
#
# Description : 
#    Returns a list of Bio::EnsEMBL::Features representing the 
#    projection of the coding regions of the given transcripts 
#    onto the query 
#################################################################

sub _get_all_transcript_cds_features {
  my ($self) = @_;

  my (@orig_cds_feats, @merged_cds_feats);

  foreach my $tran (@{$self->source_transcripts}) {
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

  $self->transcript_cds_feats(\@merged_cds_feats);
}

##################
# get/sets
#################

sub gene {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_gene} = $val;
  }

  return $self->{_gene};
}


sub name {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_name} = $val;
  }

  return $self->{_name};
}


sub source_transcripts {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_source_transcripts} = $val;

    $self->_get_all_transcript_cds_features;
  }

  return $self->{_source_transcripts};
}


sub transcript_cds_feats {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_transcript_cds} = $val;
  }

  return $self->{_transcript_cds};
}



sub good_source_transcripts {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_good_source_transcripts} = $val;
  }

  return $self->{_good_source_transcripts};
}


sub projected_transcripts {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_projected_transcripts} = $val;
  }

  return $self->{_projected_transcripts};
}







1;
