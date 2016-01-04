#!/usr/bin/env perl

# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

  my $input_id = $self->param('iid');
  my $unprocessed_genes = [];

  my $input_id_type = $self->param('iid_type');

  # If the the input is a slice, then fetch everything on the slice, but remove any gene that crosses
  # the 3' boundary, as this will be picked up on another slice
  if($input_id_type eq 'slice') {
     my $slice = $self->fetch_sequence($input_id,$dna_dba);
     $self->query($slice);
     my $gene_adaptor = $in_dba->get_GeneAdaptor();
     my $unfiltered_genes = $gene_adaptor->fetch_all_by_Slice($slice);
     $unprocessed_genes = $self->remove_3_prime_boundry_genes($unfiltered_genes);
  }

  # If the input is an array of gene ids then loop through them and fetch them from the input db
  elsif($input_id_type eq 'gene_id') {
    my $gene_adaptor = $in_dba->get_GeneAdaptor();

    for my $db_id (@{$input_id}) {
      my $unprocessed_gene = $gene_adaptor->fetch_by_dbID($db_id);
      push(@{$unprocessed_genes},$unprocessed_gene);
    }
  }

  $self->unprocessed_genes($unprocessed_genes);

#  my $transcript_weight_threshold = $self->param('transcript_weight_threshold');
#  $self->transcript_weight_threshold($transcript_weight_threshold);

  my $logic_name_weights = $self->param('logic_name_weights');
  $self->logic_name_weights($logic_name_weights);

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

  # This sub takes in a set of cluster genes (basically where each gene is set of transcripts that have been clustered
  # together based on exon overlap. Transcripts may have come from many input dbs/sources)
  # For each gene: 1) calculate the set of unique exon pairs (coordinates of consecutive exons pairs)
  #                2) score the transcripts against these pairs, transcripts were all the exon pairs pass the support
  #                   cut-off are inlcuded in the candidate transcript set
  #                3) removed redundant transcripts from the candidate set. For all transcripts that pass the support
  #                   cut-off values calculated when any transcripts are redundant and remove them. A transcript is
  #                   redundant if the structure is completely contained in another transcript. It must match exactly
  #                   in terms of being a contiguous subset of exons.
  my $output_genes = [];
  my $unprocessed_genes = $self->unprocessed_genes();
  foreach my $unprocessed_gene (@{$unprocessed_genes}) {
    my $exon_pairs = $self->generate_exon_pairs($unprocessed_gene);
    my $candidate_transcripts = $self->score_transcript_support($unprocessed_gene,$exon_pairs);
    my $final_transcripts = $self->remove_redundant_transcripts($candidate_transcripts);
    unless(scalar(@{$final_transcripts})) {
      $final_transcripts = $self->find_best_remaining_transcript($unprocessed_gene,$exon_pairs);
    }

    unless(scalar(@{$final_transcripts})) {
      next;
    }

    # Need to put a call to a sub here that reclusters the transcript set. If we get multiple clusters
    # then we need multiple genes. The problems is that they are all one gene and I would need to split
    # into multiple genes to do to. Maybe a fast check for split genes and then a call to break and
    # recluster if there is a break
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

  # This is a sub to loop through every transcript in our initial transcript cluster (the 'gene') and
  # generate a set of exon pairs. Exon pairs are a bit like intron supporting evidence, we want to
  # use them as evidence of support for sections of different transcripts. The idea is that you take
  # the start and end coords of each exon in the pair and define that as a single structure. A bit
  # like a stricter version of an intron. Then as you loop over the transcripts in the cluster you
  # count the number of times each exon pair exists. Maybe it's in every single transcript, maybe it's
  # unique to one. After it's done you have a list of pairs of exons for the gene and a list of times
  # each pair was observed. Later on the code will score the transcripts against this list and any
  # transcript that does now meet the observation cut-off for its associated logic_name/biotype for
  # every exon pair it contains will be droppped. The count itself is done by assigning all transcript
  # dbIDs associated with the pair to that pair as a set of keys. Later the count of these keys it
  # used to determine the support for a particular exon pair across the transcript cluster

  my $exon_pairs = {};

  my $transcripts = $gene->get_all_Transcripts();
  foreach my $transcript (@{$transcripts}) {
    my $exons = $transcript->get_all_translateable_Exons();
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

  # This sub takes a clustered set of transcripts (as a gene) and a set of observed exon pair coordinates
  # and uses this information to score each transcript in the gene against a weight matrix. The weight
  # matrix deserves a little explanation. It should be set in the config file and can have various levels
  # of detail. Here is a simple example:
  #
  # 'logic_name_weights' => { 'rnaseq_blast' => 4 },
  #
  # The above has a single logic_name, 'rnaseq_blast', which represents a set of models associated with this
  # logic_name. There is a single associated value '4', which is the transcript_weight_threshold. This weight
  # needs to be supported across all the exons pairs in these transcripts for a particular transcript to be
  # kept. This means that every exon pair needs to be observed at least 4 times across all the transcripts
  # in this gene. As the 'gene' represents a cluster of transcripts that might come from many input sources
  # the support need not come from other 'rnaseq_blast' transcripts. So for example if you have models made
  # by genblast or genewise that also support that exact exon pair structure, then they count towards the
  # score.
  #
  # Now for a slightly more complicated example. Say you have a set of models associated with a logic_name,
  # that you want to treat differently to the others, i.e. they are separated by biotype. So the example
  # would be the 'rnaseq_blast' logic_name having a set of models that are high quality that you want to
  # definitely include. The biotype could be 'rnaseq_80_100', representing RNA-seq models that have been verified
  # with a protein alignment of >= 80 percent coverage and identity. This could be set in the config as:
  # 'logic_name_weights' => { 'rnaseq_blast' => { 'rnaseq_80_100' => 1,
  #                                               'default_biotype_weight' => 4 } },
  #
  # In the above the 'rnaseq_80_100' biotype is given a transcript_weight_threshold of '1' ('0' would also work)
  # Note also that there is 'default_biotype_weight', which is set to 4. This should be set anytime you want
  # to specify a particular biotype (or set of biotypes) like the example above. Any biotypes you have not directly
  # assigned a weight to will get this transcript_weight_threshold. This is handy instead of having to specify
  # a weight for every single biotype. You can specifiy specific ones and have everything else use 'default_biotype_weight'
  #
  # The matrix can then be as complicated as you like in the config. Here is an example:
  # 'logic_name_weights' => { 'swiss_prot_pe12_genewise' => 1,
  #                           'rnaseq_blast'             => { 'rnaseq_80_100' => 1,
  #                                                           'default_biotype_weight' => 4 }',
  #                           'genblast'                 => { 'human_pe12' => 2,
  #                                                           'primates_pe12' => 3,
  #                                                           'mammals_pe12' => 5,
  #                                                           'vert_pe12' => 7,
  #                                                           'primates_pe345' => 10,
  #                                                         },
  #                          },
  #
  # The above has 'swiss_prot_pe12_genewise', with a threshold of '1', so keep everything from that logic_name,
  # 'rnaseq_blast' says take everything with biotype 'rnaseq_80_100', and for everything else set the threshold
  # to '4'. For 'genblast' all biotypes are defined with different cut-offs based on increasing distance from
  # the target species (we'll say baboon for that example). We trust the human protein set most, so it has a
  # modest cut-off of 2 with 'primates_pe345' being the opposite end of the spectrum, where models are only
  # chosen if all exon pairs are found 10 or more times across all the transcripts in this transcript cluster
  # Note that no 'default_biotype_weight' was set in the example for 'genblast' as all the biotypes were full
  # defined. This is okay, but in general it is safer to always define a default once you've defined a score
  # for any amount of biotypes for a particular logic name. The reason being that if you forget about a
  # particular biotype it will automatically get a cut-off of 1 (at present anyway) and thus get kept

  my $output_transcripts = [];
  my $logic_name_weights = $self->logic_name_weights();

  my $transcripts = $gene->get_all_Transcripts();
  foreach my $transcript (@{$transcripts}) {
    my $logic_name = $transcript->analysis->logic_name();
    my $biotype = $transcript->biotype();
    my $transcript_weight_threshold = 1;
    # Check if the logic_name points to a hashref
    if(ref($logic_name_weights->{$logic_name}) eq 'HASH') {
      # If we find a weight for the biotype of the transcript, then assign it
      if(exists($logic_name_weights->{$logic_name}->{$biotype})) {
        $transcript_weight_threshold = $logic_name_weights->{$logic_name}->{$biotype};
      }
      # Else if we don't find a weight for that biotype then look for a default weight key
      # If we don't find it then the biotype will just a default of 1
      elsif (exists($logic_name_weights->{$logic_name}->{'default_biotype_weight'})) {
        $transcript_weight_threshold = $logic_name_weights->{$logic_name}->{'default_biotype_weight'};
      }
    }
    # Else the logic name (hopefully) points to a weight, so set the weight cut off to that
    elsif($logic_name_weights->{$logic_name}) {
      $transcript_weight_threshold = $logic_name_weights->{$logic_name};
    }

    # Just in case someone sets something weird in the config, like pairing a biotype to undef in the hash...
    # This will cause the transcript to be included in the final set
    unless($transcript_weight_threshold) {
      $transcript_weight_threshold = 1;
    }

    my $keep_transcript = 1;
    my $exons = $transcript->get_all_Exons();
    if(scalar(@{$exons}) == 1) {
      # single exon transcript
    } else {
      # Loop through all exon pairs in the exon pair hash for the gene
      foreach my $exon_pair_key (keys(%$exon_pairs)) {
        my $exon_pair = $exon_pairs->{$exon_pair_key};
        # If this exon pair was present in the transcript then count the number to times it was observed in
        # total by counting the keys for it (which is the set of transcript dbIDs it was found for)
        if($exon_pair->{$transcript->dbID}) {
          my $support_count = scalar(keys(%$exon_pair));
          # If the any exon pair fails to pass the threshold then drop the transcript
          unless($support_count >= $transcript_weight_threshold) {
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
    my $t1_exons = $t1->get_all_translateable_Exons;
    my $j = $i+1;
    for($j=$i+1; $j<scalar@{$transcripts}; $j++) {
      my $t2 = ${$transcripts}[$j];
      my $t2_exons = $t2->get_all_translateable_Exons;

      # This section is too simple at the moment. I want a way of keeping models with UTR when the the models
      # are identical except for UTR. This will probably work, but doesn't seem like a great solution. The
      # other issue is figuring out heuristically when to skip the comparison. There's probably some way to
      # do it. Less important as this is relatively fast anyway
      if(scalar(@{$t1_exons}) > scalar(@{$t2_exons})) {
        my $is_redundant = $self->exon_subset($t1,$t2,$t1_exons,$t2_exons);
        if($is_redundant) {
          $transcript_redundancy->{$t2->dbID()} = 1;
        }
      }
      elsif(scalar(@{$t1_exons}) == scalar(@{$t2_exons})) {
        my $is_redundant = 0;
        if($t1->length >= $t2->length) {
          $is_redundant = $self->exon_subset($t1,$t2,$t1_exons,$t2_exons);
          if($is_redundant) {
            $transcript_redundancy->{$t2->dbID()} = 1;
          }
        } else {
          $is_redundant = $self->exon_subset($t2,$t1,$t2_exons,$t1_exons);
          if($is_redundant) {
            $transcript_redundancy->{$t1->dbID()} = 1;
          }
        }
      }
      else {
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
  say "Transcript A: ".$transcript_a->analysis->logic_name().", ".$transcript_a->dbID();
  say "Transcript B: ".$transcript_b->analysis->logic_name().", ".$transcript_b->dbID();

  my $i=0;
  for($i=0; $i<scalar(@{$exons_a}); $i++) {
 #   say "i index: ".$i;
    if(${$exons_a}[$i]->seq_region_start == $start_exon_b->seq_region_start &&
       ${$exons_a}[$i]->seq_region_end == $start_exon_b->seq_region_end) {
      $exon_match_count++;
#       say "In j, i index: ".$i;
      my $j=1;
      for($j=1; $j<scalar(@{$exons_b}) && ($i+$j)<scalar(@{$exons_a}); $j++) {
#        say "j index: ".$j;
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

sub remove_3_prime_boundry_genes {
  my ($self,$unfiltered_genes) = @_;

  my $filtered_genes = [];
  my $slice = $self->query();
  my $slice_3_prime_end = $slice->end;

  foreach my $gene (@{$unfiltered_genes}) {
    unless($gene->seq_region_end > $slice_3_prime_end) {
      push(@{$filtered_genes},$gene);
    } else {
      say "Skipping gene with dbID ".$gene->dbID." as it will be picked up on the adjoining another slice";
    }
  }
  return($filtered_genes);
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

#sub transcript_weight_threshold {
#  my ($self,$val) = @_;
#
#  if($val) {
#    $self->param('_transcript_weight_threshold',$val);
#  }
#
#  return($self->param('_transcript_weight_threshold'));
#}

sub logic_name_weights {
  my ($self,$val) = @_;

  if($val) {
    $self->param('_logic_name_weights',$val);
  }

  return($self->param('_logic_name_weights'));
}

sub find_best_remaining_transcript {
  my ($self,$gene,$exon_pairs) = @_;

  my $best_transcripts = [];
  my $max_support_score = 0;
  my $min_allowed_score = 2;

  # This will go through all the transcripts and pick the one with the highest level of support
  # I'm not sure how that would work. I want to avoid biasing it towards two exon genes, but
  # also don't want to pick long models just because they're long. There might be the possibility
  # of averaging across exon pairs and then adding on some constant based on exon pair count. The
  # constant could be small, like 0.5 * num_exon_pairs. But maybe by adding it like that I could
  # select longer models over shorter models, without always picking the longer model. I think if
  # a longer model has say 5 exons pairs and 4 have one support and the other has 5 then average
  # is 9/5 = 1.8. Now 1.8 + (5 * 0.5) = 4.3, whereas if you had a model that had 1 pair and 5
  # observations it would be  5 + (1 * 0.5) = 5.5, I think that's fair. What I want are longer
  # medium to well supported models to be chosen over shorter medium to well supported models
  # I might make the constant smaller than 0.5 though
  my $transcripts = $gene->get_all_Transcripts();
  foreach my $transcript (@{$transcripts}) {
    my $support_score = $self->calculate_support_average($transcript,$exon_pairs);
    my $biotype = $transcript->biotype();
    $biotype .= "_br";
    $transcript->biotype($biotype);
    if($support_score == $max_support_score && $support_score >= $min_allowed_score) {
      push(@{$best_transcripts},$transcript);
    } elsif($support_score > $max_support_score && $support_score >= $min_allowed_score) {
      $best_transcripts = [];
      push(@{$best_transcripts},$transcript);
      $max_support_score = $support_score;
    }
  }

  return $best_transcripts;
}

sub calculate_support_average {
  my ($self,$transcript,$exon_pairs) = @_;

  my $total_support_amount = 0;
  my $total_exon_pairs = 0;
  my $exon_pair_count_bonus = 0.2;
  my $final_support = 0;
  my $exons = $transcript->get_all_translateable_Exons();

  if(scalar(@{$exons}) == 1) {
     # single exon transcript
  } else {
    # Loop through all exon pairs in the exon pair hash for the gene
    foreach my $exon_pair_key (keys(%$exon_pairs)) {
      my $exon_pair = $exon_pairs->{$exon_pair_key};
      if($exon_pair->{$transcript->dbID}) {
        $total_exon_pairs++;
        $total_support_amount += scalar(keys(%$exon_pair));
      }
    }

    my $average_support = $total_support_amount / $total_exon_pairs;
    $final_support = $average_support + ($total_exon_pairs * $exon_pair_count_bonus);
  } # end else

  return($final_support);
}

1;
