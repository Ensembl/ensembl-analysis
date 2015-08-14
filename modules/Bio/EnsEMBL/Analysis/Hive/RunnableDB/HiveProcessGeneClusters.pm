#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProcessGeneClusters;

use strict;
use warnings;
use feature 'say';
use Data::Dumper;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

# This module is built with the view of in some ways being the opposite to our current
# set of modules which include TranscriptConsensus, LayerAnnotation and GeneBuilder.
# These modules have always been run with the view of reducing a large set of transcripts
# to one or very few well supported transcripts. This is not biologically realistic. This
# module will instead attempt to find as many well supported transcripts as possible

sub fetch_input {
  my $self = shift;

  # This will take a set of input dbs and a set of logic_name and allowed biotypes
  # It will then retrieve all the transcripts that fit this set. If no allowed biotypes
  # hash exists for the set it will retrieve all the transcripts for that logic_name from
  # every db it exists in (though I would hope for logic_names to not be spread across
  # several dbs, as the set associated with each logic_name should be non-redundant).
  # Also I guess if no logic names are provided then just take all from the input dbs
  # Once fetch input is done there should be a geneset associated with some of the major
  # groupings

  my $analysis = Bio::EnsEMBL::Analysis->new(
                                              -logic_name => $self->param('logic_name'),
                                              -module => $self->param('module'),
                                            );
  $self->analysis($analysis);

  my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
  $self->hrdb_set_con($dna_dba,'dna_db');
  my $in_dba = $self->hrdb_get_dba($self->param('cluster_db'));
  $in_dba->dnadb($dna_dba);
  $self->hrdb_set_con($in_dba,'cluster_db');
  my $out_dba = $self->hrdb_get_dba($self->param('processed_cluster_db'));
  $out_dba->dnadb($dna_dba);
  $self->hrdb_set_con($out_dba,'processed_cluster_db');

  my $input_id_array = $self->param('iid');
  my $gene_adaptor = $in_dba->get_GeneAdaptor();
  my $unprocessed_genes = [];

  for my $db_id (@{$input_id_array}) {
    my $unprocessed_gene = $gene_adaptor->fetch_by_dbID($db_id);
    push(@{$unprocessed_genes},$unprocessed_gene);
  }

  $self->unprocessed_genes($unprocessed_genes);

  return 1;
}

sub run {
  my $self = shift;

  my $output_genes = $self->process_rough_genes();
  $self->output_genes($output_genes);

  return 1;
}

sub write_output {
  my $self = shift;

  my $adaptor = $self->hrdb_get_con('processed_cluster_db')->get_GeneAdaptor;
  my @output = @{$self->output_genes};

  say "Writing genes to output db";
  foreach my $gene (@output){
    empty_Gene($gene);
    $adaptor->store($gene);
  }
  say "...finished writing genes to output db";

  return 1;
}

sub process_rough_genes {
  my ($self) = @_;

  my $output_genes = [];
  my $unprocessed_genes = $self->unprocessed_genes();
  foreach my $unprocessed_gene (@{$unprocessed_genes}) {
    my $exon_pairs = $self->generate_exon_pairs($unprocessed_gene);
    my $candidate_transcripts = $self->score_transcript_support($unprocessed_gene,$exon_pairs);
    my $final_transcripts = $self->remove_redundant_transcripts($candidate_transcripts);

    unless(scalar(@{$final_transcripts})) {
      next;
    }

    my $output_gene = Bio::EnsEMBL::Gene->new();
    $output_gene->analysis($self->analysis);
    $output_gene->biotype($self->analysis->logic_name);
    foreach my $final_transcript (@{$final_transcripts}) {
      $output_gene->add_Transcript($final_transcript);
    }
    push(@{$output_genes},$output_gene);
  }

  return($output_genes);
}

sub generate_exon_pairs {
  my ($self,$gene) = @_;

  my $exon_pairs = {};

  my $transcripts = $gene->get_all_Transcripts();
  foreach my $transcript (@{$transcripts}) {
    my $exons = $transcript->get_all_Exons();
    if(scalar(@{$exons}) == 1) {
      # single exon transcript
    } else {
      my $i=0;
      for($i=0; $i<(scalar(@{$exons})-1); $i++) {
        my $exon_left = ${$exons}[$i];
        my $exon_right = ${$exons}[$i+1];
        my $els = $exon_left->seq_region_start;
        my $ele = $exon_left->seq_region_end;
        my $ers = $exon_right->seq_region_start;
        my $ere = $exon_right->seq_region_end;
        my $coord_string = $els.":".$ele.":".$ers.":".$ere;
        $exon_pairs->{$coord_string}->{$transcript->dbID} = 1;
      }
    } # end else
  }
  return($exon_pairs);
}

sub score_transcript_support {
  my ($self,$gene,$exon_pairs) = @_;

  my $output_transcripts = [];
  my $support_cut_off = 4;
  my $transcripts = $gene->get_all_Transcripts();
  foreach my $transcript (@{$transcripts}) {
    my $keep_transcript = 1;
    my $exons = $transcript->get_all_Exons();
    if(scalar(@{$exons}) == 1) {
      # single exon transcript
    } else {
      foreach my $exon_pair_key (keys(%$exon_pairs)) {
        my $exon_pair = $exon_pairs->{$exon_pair_key};
        if($exon_pair->{$transcript->dbID}) {
          my $support_count = scalar(keys(%$exon_pair));
          unless($support_count >= $support_cut_off) {
            $keep_transcript = 0;
          }
        }
      }
      if($keep_transcript) {
        say "Keeping transcript: ".$transcript->stable_id();
        push(@{$output_transcripts},$transcript);
      }
    } # end else
  }
  return($output_transcripts);
}

sub remove_redundant_transcripts {
  my ($self,$transcripts) = @_;

  my $final_transcripts = [];
  my $transcript_redundancy = {};
  my $i=0;
  for($i=0; $i<scalar@{$transcripts}; $i++) {
    my $t1 = ${$transcripts}[$i];
    my $t1_exons = $t1->get_all_Exons;
    my $j = $i+1;
    for($j=$i+1; $j<scalar@{$transcripts}; $j++) {
      my $t2 = ${$transcripts}[$j];
      my $t2_exons = $t2->get_all_Exons;
      if(scalar(@{$t1_exons}) >= scalar(@{$t2_exons})) {
        my $is_redundant = $self->exon_subset($t1,$t2,$t1_exons,$t2_exons);
        if($is_redundant) {
          $transcript_redundancy->{$t2->dbID()} = 1;
        }
      } else {
         my $is_redundant = $self->exon_subset($t2,$t1,$t2_exons,$t1_exons);
         if($is_redundant) {
          $transcript_redundancy->{$t1->dbID()} = 1;
        }
      }
    }
  }

  foreach my $transcript (@{$transcripts}) {
    unless($transcript_redundancy->{$transcript->dbID}) {
      push(@{$final_transcripts},$transcript)
    }
  }

  return $final_transcripts;
}

sub exon_subset {
  my ($self,$transcript_a,$transcript_b,$exons_a,$exons_b) = @_;

  my $is_subset = 0;
  my $start_exon_b = ${$exons_b}[0];
  my $exon_match_count = 0;
  my $i=0;
  for($i=0; $i<scalar(@{$exons_a}); $i++) {
    if(${$exons_a}[$i]->seq_region_start == $start_exon_b->seq_region_start &&
       ${$exons_a}[$i]->seq_region_end == $start_exon_b->seq_region_end) {
      $exon_match_count++;
      my $j=1;
      for($j=1; $j<scalar(@{$exons_b}); $j++) {
        if(${$exons_a}[$i+$j]->seq_region_start == ${$exons_b}[$j]->seq_region_start &&
          ${$exons_a}[$i+$j]->seq_region_end == ${$exons_b}[$j]->seq_region_end) {
          $exon_match_count++;
        }
      }
    }
  }

  if($exon_match_count == scalar(@{$exons_b})) {
    say "Model ".$transcript_b->stable_id()." is redundant to model ".$transcript_a->stable_id();
    $is_subset = 1;
  }
  return $is_subset;
}
sub unprocessed_genes {
  my ($self,$val) = @_;

  if($val) {
    $self->param('_unprocessed_genes',$val);
  }

  return($self->param('_unprocessed_genes'));
}

sub output_genes {
  my ($self,$val) = @_;

  if($val) {
    $self->param('_output_genes',$val);
  }

  return($self->param('_output_genes'));
}

1;
