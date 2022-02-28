=head1 LICENSE

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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::WGA2Genes - 

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


=head1 METHODS

=cut
package Bio::EnsEMBL::Analysis::RunnableDB::WGA2Genes;

require Exporter;
use vars qw(@ISA @EXPORT);
use strict;

use warnings;
use Data::Dumper;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;

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

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils 
    qw(replace_stops_with_introns);

use IO::String;
use Bio::SeqIO;

use DBI qw(:sql_types);
use DBD::mysql;
@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB 
          Exporter Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);
@EXPORT = (@{$DBI::EXPORT_TAGS{'sql_types'}});


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


  my $t_dbh = $self->get_dbadaptor($self->TARGET_CORE_DB) ;
  my $q_dbh = $self->get_dbadaptor($self->QUERY_CORE_DB, '', 1) ;
  my $compara_dbh = $self->get_dbadaptor($self->COMPARA_DB, 'compara') ;
  
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
          # load seq up front. Otherwise get warnings in the sort
          map { $_->translateable_seq } @good_trans;
          
          @good_trans = sort { 
            CORE::length($b->translateable_seq) <=> CORE::length($a->translateable_seq);
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
    my ($gene, $tran);
    eval {
      $gene = $q_dbh->get_GeneAdaptor->fetch_by_stable_id($input_id);
    };
    eval {
      $tran = $q_dbh->get_TranscriptAdaptor->fetch_by_stable_id($input_id);
    };
    if (not defined $gene and not defined $tran) {
      throw("Could not find gene or transcript with '$input_id' in query database");
    }

    if (not defined $gene) {
      $gene = Bio::EnsEMBL::Gene->new;
      $gene->stable_id("None");
      $gene->add_Transcript($tran);
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
      if ($self->LONGEST_SOURCE_TRANSCRIPT) {
        # load seq up front. Otherwise get warnings in the sort
        map { $_->translateable_seq } @good_trans;

        @good_trans = sort { 
          CORE::length($b->translateable_seq) <=> CORE::length($a->translateable_seq);
        } @good_trans;
        @good_trans = ($good_trans[0]);
      }

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

    if ($block->reference_genomic_align->dnafrag_strand < 0) {
      $block->reverse_complement;
    }

    push @{$chains{$block->group_id}}, $block;
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

  my $external_db_id = $self->get_external_db_id();

  logger_info("INITIAL GENES: " . join(" ", map {$_->gene->stable_id}
  @{$self->gene_records}) . "\n");

  logger_info("ORIGINAL_CHAINS\n" . 
              stringify_chains($self->genomic_align_block_chains));

  @chains_generecs_pairs = 
      $self->segregate_chains_and_generecords($self->genomic_align_block_chains,
                                              $self->gene_records,
                                              $self->MIN_CHAIN_ANNOTATION_OVERLAP);
      
  foreach my $chains_generecs (@chains_generecs_pairs) {
    my ($alignment_chains, $generecs) = @$chains_generecs;
    
    my $pair_id = join(":", 
                       map { $_->gene->stable_id } @$generecs);
    
    logger_info("CHAINS_FOR_GENESET $pair_id\n" . 
                stringify_chains($alignment_chains));
    
    my %blocks_from_used_chains;
    
    my @these_generecs = map { $_->clone } @$generecs;
    
    # Repeat transcript contruction until we had no rejected transcripts.
    # This is because if a transcript is rejected, some of the underlying 
    # GeneScaffold components may be unnecessary
    for(;;) {
      logger_info("Working with genes " . join(" ", map { $_->gene->stable_id} @these_generecs));
      
      my @cds_feats = map { @{$_->transcript_cds_feats} } @these_generecs;
      
      my $filtered_chains =             
          filter_irrelevant_chains($alignment_chains,
                                   \@cds_feats);
      
      logger_info("RELEVANT_CHAINS\n" . stringify_chains($filtered_chains));
      
      $filtered_chains = filter_inconsistent_chains($filtered_chains, 
                                                    $self->OVERLAP_CHAIN_FILTER);
      logger_info("CONSISTENT_CHAINS\n" . stringify_chains($filtered_chains));
      
      $filtered_chains = $self->remove_contig_split_chains($filtered_chains);
      logger_info("NON_CONTIG_SPLIT_CHAINS\n" . stringify_chains($filtered_chains));

      $filtered_chains = $self->trim_exon_split_chains($filtered_chains, 
                                                       \@cds_feats);
      logger_info("TRIMMED_CHAINS\n" . stringify_chains($filtered_chains));
      
      my @these_pairs = $self->segregate_chains_and_generecords($filtered_chains,
                                                                \@these_generecs,
                                                                $self->MIN_CHAIN_ANNOTATION_OVERLAP);
      
      my @these_results;
      
      foreach my $pair (@these_pairs) {
        my ($subset_chains, $subset_generecs) = @$pair;
        
        my $net_blocks = flatten_chains($subset_chains, 1);

        my $gs_name = $subset_generecs->[0]->gene->stable_id . "-0";
        my $gene_scaffold = $self->make_gene_scaffold_and_project_genes($net_blocks,
                                                                        $subset_generecs,
                                                                        $gs_name,
                                                                        $self->MAX_EDITABLE_STOPS_PRIMARY,
                                                                        $self->MIN_COVERAGE,
                                                                        $self->MIN_NON_GAP,
                                                                        $self->FULLY_GAP_EXONS,
                                                                        $self->PARTIAL_GAP_EXONS,
                                                                        $external_db_id);
        
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
            push @kept_generecs, $res_gene;
          } else {
            if (@{$res_gene->good_source_transcripts}) {
              push @kept_generecs, $res_gene;
            }
            $need_to_repeat = 1;
          }
        }
      }
      
      if ($need_to_repeat) {
        if (@kept_generecs) {
          @these_generecs = @kept_generecs;
          foreach my $grec (@these_generecs) {
            $grec->source_transcripts($grec->good_source_transcripts);
            $grec->projected_transcripts([]);
            $grec->good_source_transcripts([]);
          }
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
    
    my $filtered_chains = filter_block_overlap_chains($alignment_chains,
                                                      [values %blocks_from_used_chains]);
    
    logger_info("ITERATING over remaining single chains\n");

    foreach my $chain (@$filtered_chains) {
      my ($single_pair) = $self->segregate_chains_and_generecords([$chain],
                                                                  $generecs,
                                                                  $self->MIN_CHAIN_ANNOTATION_OVERLAP);

      if (defined $single_pair) {
        logger_info("CHAIN_WITH_GENES\n" . stringify_chains([$chain]));

        my @grs = map { $_->clone } @{$single_pair->[1]};;

        my ($tg_al) = @{$chain->[0]->get_all_non_reference_genomic_aligns};

        my $gs_name = $tg_al->dnafrag->name . "-" . $tg_al->dnafrag_start;

        
        my $gene_scaffold = $self->make_gene_scaffold_and_project_genes($chain,
                                                                        \@grs,
                                                                        $gs_name,
                                                                        $self->MAX_EDITABLE_STOPS_NON_PRIMARY,
                                                                        $self->MIN_COVERAGE,
                                                                        $self->MIN_NON_GAP,
                                                                        $self->FULLY_GAP_EXONS,
                                                                        $self->PARTIAL_GAP_EXONS,
                                                                        $external_db_id);
        my @kept;
        foreach my $gr (@grs) {
          if (scalar(@{$gr->projected_transcripts})) {
            push @kept, $gr;
          }
        }
        if (@kept) {
          push @results, [$gene_scaffold, @kept];
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

  my $outfile = $self->create_output_dir;

  # This is just to make sure that the output files written properly
  # and contain all the information that they should. This is done by
  # counting the number of lines in each file and comparing them to
  # what is actually printed out.

  open(FH, ">$outfile") or throw("Failed to open $outfile!\n$!");

  my $fh = \*FH;
  my $line_numbers = 0;

  print ($fh "#\n");
  $line_numbers++;
  printf($fh "# WGA2Genes output for %s\n", $self->input_id);
  $line_numbers++;
  printf($fh "#  Genes are:");
  $line_numbers++;
  foreach my $grec (@{$self->gene_records}) {
    print ($fh " " . $grec->gene->stable_id);
  }
  print ($fh "\n#\n");
  $line_numbers++;
  #print "OUTFILE ".$outfile." used\n";

  foreach my $obj (@{$self->output}) {
    my ($gs, @res_genes) = @$obj;

    $line_numbers += $self->write_agp($gs, $fh);
    foreach my $g (@res_genes) {
      $line_numbers += $self->write_gene($gs, 
                        $g,
                        $fh);
    }
  }

  if (!-e $outfile) { throw("The $outfile was not created.\n"); }

  #print $line_numbers, "<-- counted number of lines\n";

  my $cmd = "wc -l $outfile";

  # Run the command and capture the output
  my $cmd_string = `$cmd`;
  my @wc = (split(/ /, $cmd_string));
  #print $wc[0], "<-- what should have been written\n";
  if ($wc[0] !~ $line_numbers) {
    throw("The number of lines in the output files differ to what should have
    been saved. You might need to rerun the pipeline\n");
  }

  close($fh);
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
  my ($self, 
      $blocks, 
      $res_genes, 
      $name,
      $max_stops,
      $min_coverage,
      $min_non_gap,
      $add_gaps,
      $extend_into_gaps,
      $external_db_id) = @_;

  my $gene_scaffold = Bio::EnsEMBL::Analysis::Tools::WGA2Genes::GeneScaffold->
      new(-name => $name,
          -genomic_align_blocks => $blocks,
          -from_slice => $self->query_slice,
          -to_slices  => $self->target_slices,
          -transcripts   => [map {@{$_->source_transcripts}} @$res_genes],
          -max_readthrough_dist => $self->MAX_EXON_READTHROUGH_DIST,
          -add_gaps => $add_gaps,
          -extend_into_gaps => $extend_into_gaps,
          );
  
  foreach my $res_gene (@$res_genes) {
    $res_gene->name($name . "." . $res_gene->gene->stable_id);
    
    foreach my $tran (@{$res_gene->source_transcripts}) {
      my $proj_trans = $gene_scaffold->place_transcript($tran, 1, $external_db_id); 

     if ( $tran->analysis) { 
        $proj_trans->analysis($tran->analysis);     
        for my $e( @{$proj_trans->get_all_Exons } ) {  
            for my $sf ( @{ $e->get_all_supporting_features } ) {  
                 $sf->analysis($tran->analysis); 
             } 
        }
     } 
      $proj_trans = 
          $self->process_transcript($proj_trans, 
                                    $max_stops,
                                    $min_coverage,
                                    $min_non_gap,
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

  if (CORE::length($pep) == 0) {
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
  $tran = replace_stops_with_introns($tran);

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
# FUNCTION: trim_exon_split_chains
#
# Decription:
#  When a single exon overlaps blocks in more than one chain, the
#  minority chains are trimmed back.
#  This prevents the situation of a single-exon gene in the 
#  reference mapping down to a multi-scaffold GeneScaffold, and
#  generally gives cleaner results.
###################################################################
sub trim_exon_split_chains {
  my ($self, $chains, $feats) = @_;

  # make nr list of coding regions
  my (@nr_cds);

  foreach my $cds (sort { $a->start <=> $b->start } @$feats) {
    if (not @nr_cds or $cds->start > $nr_cds[-1]->end + 1) {
      push @nr_cds, $cds;
    } elsif ($cds->end > $nr_cds[-1]->end) {
      $nr_cds[-1]->end($cds->end);
    }
  }

  foreach my $cds (@nr_cds) {
    my (@ov_chains, @other_chains);

    CH: foreach my $ch (@$chains) {
      my $ov_bps = 0;
      foreach my $bl (@$ch) {
        my $ga = $bl->reference_genomic_align;
        if ($ga->dnafrag_start <= $cds->end and
            $ga->dnafrag_end >= $cds->start) {
          my ($st, $en) = ($cds->start, $cds->end);
          $st = $ga->dnafrag_start if $ga->dnafrag_start > $st;
          $en = $ga->dnafrag_end if $ga->dnafrag_end < $en;

          $ov_bps += ($en - $st + 1);          
        }
      }
      if ($ov_bps) {
        push @ov_chains, {
          chain => $ch,
          ov_bps => $ov_bps,
        };
      } else {
        push @other_chains, $ch;
      }
    }
    
    if (scalar(@ov_chains) > 1) {
      # more than one chain implicated. Find the "dominant" chain
      # for this exon and trim back the blocks for the other ones
      @ov_chains = sort { $b->{ov_bps} <=> $a->{ov_bps} } @ov_chains;
      my @new_chains = (@other_chains, $ov_chains[0]->{chain});

      for (my $i=1; $i < @ov_chains; $i++) {
        my @retained_blocks;
        foreach my $bl (@{$ov_chains[$i]->{chain}}) {
          my $ga = $bl->reference_genomic_align;
          if ($ga->dnafrag_end < $cds->start or $ga->dnafrag_start > $cds->end) {
            push @retained_blocks, $bl;
          } else {
            if ($ga->dnafrag_start < $cds->start) {
              my $cut_bl = $bl->restrict_between_reference_positions($ga->dnafrag_start, 
                                                                     $cds->start - 1);
              $cut_bl->score($bl->score);
              push @retained_blocks, $cut_bl; 

            }
            if ($ga->dnafrag_end > $cds->end) {
              my $cut_bl = $bl->restrict_between_reference_positions($cds->end + 1, 
                                                                     $ga->dnafrag_end);
              $cut_bl->score($bl->score);
              push @retained_blocks, $cut_bl;
            }
          }
        }
        if (@retained_blocks) { 
          push @new_chains, \@retained_blocks;
        }
      }
      $chains = \@new_chains;
    }
  }

  return $chains;
}



###################################################################
# FUNCTION: segregate_chains_and_generecords
#
# Decription:
#   partition the chains and genes into (chain_list, gene_list)
#   pairs, where each chain_list contains chains that touch
#   one or more genes in gene_list
###################################################################
sub segregate_chains_and_generecords {
  my ($self, $chains, $generecs, $min_overlap) = @_;

  my (%generecs_by_id, %chains_by_id);

  foreach my $grec (@$generecs) {
    $generecs_by_id{$grec} = $grec;
  }

  my (%genes_per_chain);

  foreach my $c (@$chains) {
    $chains_by_id{$c} = $c;

    my @ref_gas = map { $_->reference_genomic_align } @$c;
    @ref_gas = sort { $a->dnafrag_start <=> $b->dnafrag_start } @ref_gas;
    
    my %overlapping;
    
    foreach my $grec (@$generecs) {
      my $overlap_bps = 0;

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
            my $ov_start = $f->start;
            my $ov_end   = $f->end;

            $ov_start = $bl->dnafrag_start if $bl->dnafrag_start > $ov_start;
            $ov_end  = $bl->dnafrag_end if $bl->dnafrag_end < $ov_end;

            $overlap_bps += ($ov_end - $ov_start + 1);
          }
        }
      }

      if ($overlap_bps >= $min_overlap) {
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
      $gscaf,
      $fh) = @_;

  $fh = \*STDOUT if(!$fh);

  my $prefix = "##-AGP";

  my $line_count = 0;
  my @tsegs = $gscaf->project_down;
  my @qsegs = $gscaf->project_up;

  print $fh "$prefix \#\n";
  $line_count++;

  printf($fh "$prefix \#\# AGP for %s source-region=%s/%d-%d (iid=%s)\n", 
         $gscaf->seq_region_name,
         $qsegs[0]->to_Slice->seq_region_name,
         $qsegs[0]->to_Slice->start,
         $qsegs[-1]->to_Slice->end,
         $self->input_id);
  $line_count++;

  print $fh "$prefix \#\n";
  $line_count++;

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
    $line_count++;

    if ($i < @tsegs - 1) {
      printf($fh "$prefix %s\t%d\t%d\t%s\t%s\t%d\n",
             $gscaf->seq_region_name,
             $seg->from_end + 1,
             $tsegs[$i+1]->from_start - 1,
             $piece_count++,
             "N",
             $tsegs[$i+1]->from_start - $seg->from_end - 1,
             );
      $line_count++;
    }
  }
  return $line_count;
}


sub write_gene {
  my ($self,
      $gene_scaf,
      $g,
      $fh ) = @_;

  $fh = \*STDOUT if(!$fh);

  my $prefix = "##-GENES ";

  my $line_count = 0;

  my $seq_id = $gene_scaf->seq_region_name;
  my $gene_id = $g->name;

  printf $fh "$prefix \# Gene report for $seq_id\n";
  $line_count++;

  foreach my $tran (@{$g->projected_transcripts}) {
    my ($sf) = @{$tran->get_all_supporting_features};
    my $tran_id = $gene_id . "_" . $sf->hseqname; 

    foreach my $attr (@{$tran->get_all_Attributes}) {
      printf($fh "$prefix \##\-ATTRIBUTE transcript=$tran_id code=%s value=%s\n",
             $attr->code, 
             $attr->value);
      $line_count++;
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
      $line_count++;

      foreach my $sup_feat (@{$exon->get_all_supporting_features}) {
        foreach my $ug ($sup_feat->ungapped_features) {
          printf($fh "%s %s\t%s\t%s\t%d\t%d\t%d\t%s=%s; %s=%s; %s=%s; %s=%s; %s=%s; %s=%s\n", 
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
                 $ug->hend,
                 "external_db_id",
                 $ug->external_db_id,
                 "hcoverage",
                 $ug->hcoverage);
          $line_count++;
        }
      }
    }

    my $pep = $tran->translate;
    $pep->id($tran_id);

    my $fasta_string;
    my $stringio = IO::String->new($fasta_string);
    my $seqio = Bio::SeqIO->new(-format => 'fasta',
                                -fh => $stringio);
    $seqio->write_seq($pep);
    logger_info("PROJECTED-PEPTIDE:\n$fasta_string\n"); 
    if ($pep->seq =~ /^X/ or $pep->seq =~ /X$/) {
      logger_info("Transcript $tran_id has Xs at the ends");
    }
  }
  return $line_count;
}

sub create_output_dir {
  my ($self, $output_dir_number) = @_;
  $output_dir_number = 10 if(!$output_dir_number);
  my $num = int(rand($output_dir_number));
  my $output_dir = $self->OUTPUT_DIR;
  if (! -e $output_dir) {
    warning("Your output directory '" . $output_dir . "' does not exist");
    warning("- it will be created now\n");
    eval{
      system("mkdir -p " . $output_dir);
    };
    if($@){
      throw("Failed to make output directory " . $output_dir . "$@" );
    }
  }
  my $dir = $output_dir."/".$num;
  if(! -e $dir){
    my $command = "mkdir $dir";
    eval{
      system($command);
    };
    if($@){
      throw("Failed to make $dir $@");
    }
  }
  my $filename = $dir. "/" . $self->input_id."_".$self->analysis->logic_name."_";
  $filename .= int(rand(1000));
  $filename .= ".out";
  return $filename;
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
                      COMPARA_DB
                      TARGET_SPECIES_EXTERNAL_DBNAME)) {

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

sub TARGET_SPECIES_EXTERNAL_DBNAME {
  my ($self, $target_species_external_dbname) = @_;

  if (defined $target_species_external_dbname) {
    $self->{_target_species_external_dbname} = $target_species_external_dbname;
  }

  return $self->{_target_species_external_dbname};
}

#
# chain filtering
#

sub MIN_CHAIN_ANNOTATION_OVERLAP {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_min_chain_annotation_overlap} = $val;
  }

  if (not exists($self->{_min_chain_annotation_overlap})) {
    $self->{_min_chain_annotation_overlap} = 1;
  }

  return $self->{_min_chain_annotation_overlap};
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

sub FULLY_GAP_EXONS {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_fully_gap_exons} = $val; 
  }

  return $self->{_fully_gap_exons};
}


sub PARTIAL_GAP_EXONS {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_partial_gap_exons} = $val; 
  }

  return $self->{_partial_gap_exons};
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

sub OUTPUT_DIR {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_output_dir} = $val;
  }

  return $self->{_output_dir};
}

sub KILL_LIST {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_kill_list} = $val;
  }

  return $self->{_kill_list};
}

sub get_external_db_id {
  my ($self) = @_;

  my $query_db = $self->get_dbadaptor($self->QUERY_CORE_DB, '', 1) ;
  my $external_dbname = $self->TARGET_SPECIES_EXTERNAL_DBNAME;

  my $sth = $query_db->prepare(
            "SELECT external_db_id ".
            "FROM external_db ".
            "WHERE db_name = ?"
                          );
  $sth->bind_param(1, $external_dbname, SQL_VARCHAR);
  $sth->execute();

  my ($external_db_id) = $sth->fetchrow();
  $sth->finish();

  return $external_db_id;
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
