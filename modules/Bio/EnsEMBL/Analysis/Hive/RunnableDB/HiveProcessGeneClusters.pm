#!/usr/bin/env perl

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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
  my $out_dba = $self->hrdb_get_dba($self->param('cluster_output_db'));
  $out_dba->dnadb($dna_dba);
  $self->hrdb_set_con($out_dba,'cluster_output_db');

  my $input_id = $self->param('iid');
  my $unprocessed_genes = [];

  my $input_id_type = $self->param('iid_type');

  # If the the input is a slice, then fetch everything on the slice, but remove any gene that crosses
  # the 3' boundary, as this will be picked up on another slice
  if($input_id_type eq 'slice') {
    my $in_dba = $self->hrdb_get_dba($self->param('cluster_db'));
    $in_dba->dnadb($dna_dba);
    $self->hrdb_set_con($in_dba,'cluster_db');
    my $slice = $self->fetch_sequence($input_id,$dna_dba);
    $self->query($slice);
    my $gene_adaptor = $in_dba->get_GeneAdaptor();
    my $unfiltered_genes = $gene_adaptor->fetch_all_by_Slice($slice);
    $unprocessed_genes = $self->remove_3_prime_boundry_genes($unfiltered_genes);
  }

  # If the input is an array of gene ids then loop through them and fetch them from the input db
  elsif($input_id_type eq 'gene_id') {
    my $in_dba = $self->hrdb_get_dba($self->param('cluster_db'));
    $in_dba->dnadb($dna_dba);
    $self->hrdb_set_con($in_dba,'cluster_db');
    my $gene_adaptor = $in_dba->get_GeneAdaptor();

    for my $db_id (@{$input_id}) {
      my $unprocessed_gene = $gene_adaptor->fetch_by_dbID($db_id);
      push(@{$unprocessed_genes},$unprocessed_gene);
    }
  }

  # If the input is an array of gene ids then loop through them and fetch them from the input db
  elsif($input_id_type eq 'cluster') {
     my $unprocessed_gene = $self->load_cluster($input_id);
     push(@{$unprocessed_genes},$unprocessed_gene);
  }

  $self->unprocessed_genes($unprocessed_genes);

  my $logic_name_weights = $self->param('logic_name_weights');
  $self->logic_name_weights($logic_name_weights);

  my $single_exon_support_penalty = 1;
  if($self->param('single_exon_support_penalty')) {
    $single_exon_support_penalty = $self->param('single_exon_support_penalty');
  }
  $self->single_exon_support_penalty($single_exon_support_penalty);

  return 1;
}

sub run {
  my $self = shift;

  my $processed_genes = $self->process_rough_genes();
  my $final_genes = $self->recluster_genes($processed_genes);
  $self->output_genes($final_genes);

  return 1;
}

sub write_output {
  my $self = shift;

  my $adaptor = $self->hrdb_get_con('cluster_output_db')->get_GeneAdaptor;
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
  #                2) score the transcripts against these pairs, transcripts where all the exon pairs pass the support
  #                   cut-off are inlcuded in the candidate transcript set
  #                3) remove redundant transcripts from the candidate set. For all transcripts that pass the support
  #                   cut-off values calculate when any transcripts are redundant and remove them. A transcript is
  #                   redundant if the structure is completely contained in another transcript. It must match exactly
  #                   in terms of being a contiguous subset of exons.
  my $output_genes = [];
  my $unprocessed_genes = $self->unprocessed_genes();

  foreach my $unprocessed_gene (@{$unprocessed_genes}) {
    my $input_transcript_count = scalar(@{$unprocessed_gene->get_all_Transcripts});
    my $exon_pairs = $self->generate_exon_pairs($unprocessed_gene);
    my $candidate_transcripts = $self->score_transcript_support($unprocessed_gene,$exon_pairs);
    my $candidate_transcript_count = scalar(@{$candidate_transcripts});
#    my $near_final_transcripts = $self->remove_redundant_transcripts($candidate_transcripts);
    my $final_transcripts = $self->compare_transcripts($candidate_transcripts);     #$near_final_transcripts);
    unless(scalar(@{$final_transcripts})) {
      $final_transcripts = $self->find_best_remaining_transcript($unprocessed_gene,$exon_pairs);
    }

    unless(scalar(@{$final_transcripts})) {
      next;
    }

    my $final_transcript_count = scalar(@{$final_transcripts});
    say "Inital transcript count for gene: ".$input_transcript_count;
    say "Candidate transcript count for gene: ".$candidate_transcript_count;
    say "Final transcript count for gene: ".$final_transcript_count;

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
      my $exon = shift(@{$exons});
      my $exon_start = $exon->seq_region_start;
      my $exon_end = $exon->seq_region_end;
      # For a single exon model we will just build a pair with the exon coords repeated
      # Later, in the scoring phase we can add a penalty to the cut off for these
      my $coord_string = $exon_start.":".$exon_end.":".$exon_start.":".$exon_end;
      $exon_pairs->{$coord_string}->{$transcript->{'internal_transcript_id'}} = 1;
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
        $exon_pairs->{$coord_string}->{$transcript->{'internal_transcript_id'}} = 1;
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
  # Note that no 'default_biotype_weight' was set in the example for 'genblast' as all the biotypes were fully
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
    my $single_exon_support_penalty = $self->single_exon_support_penalty();

    # Check if the logic_name points to a hashref
    if(ref($logic_name_weights->{$logic_name}) eq 'HASH') {
      # If we find a weight for the biotype of the transcript, then assign it
      if(exists($logic_name_weights->{$logic_name}->{$biotype})) {
        $transcript_weight_threshold = $logic_name_weights->{$logic_name}->{$biotype};
        $transcript->{'transcript_weight_threshold'} = $transcript_weight_threshold;
      }
      # Else if we don't find a weight for that biotype then look for a default weight key
      # If we don't find it then the biotype will just a default of 1
      elsif (exists($logic_name_weights->{$logic_name}->{'default_biotype_weight'})) {
        $transcript_weight_threshold = $logic_name_weights->{$logic_name}->{'default_biotype_weight'};
        $transcript->{'transcript_weight_threshold'} = $transcript_weight_threshold;
      }
    }
    # Else the logic name (hopefully) points to a weight, so set the weight cut off to that
    elsif($logic_name_weights->{$logic_name}) {
      $transcript_weight_threshold = $logic_name_weights->{$logic_name};
      $transcript->{'transcript_weight_threshold'} = $transcript_weight_threshold;
    }

    # Just in case someone sets something weird in the config, like pairing a biotype to undef in the hash...
    # This will cause the transcript to be included in the final set
    unless($transcript_weight_threshold) {
      $transcript_weight_threshold = 1;
      $transcript->{'transcript_weight_threshold'} = $transcript_weight_threshold;
    }

    my $keep_transcript = 1;
    my $exons = $transcript->get_all_Exons();

    # Loop through all exon pairs in the exon pair hash for the gene
    foreach my $exon_pair_key (keys(%$exon_pairs)) {
      my $exon_pair = $exon_pairs->{$exon_pair_key};
      # If this exon pair was present in the transcript then count the number to times it was observed in
      # total by counting the keys for it (which is the set of transcript dbIDs it was found for)
      if($exon_pair->{$transcript->{'internal_transcript_id'}}) {
        my $support_count = scalar(keys(%$exon_pair));
        # If we have a single exon transcript we add a penalty for the support threshold
        if(scalar(@{$exons}) == 1) {
          unless($support_count >= ($transcript_weight_threshold + $single_exon_support_penalty)) {
            $keep_transcript = 0;
          }
        }
        else {
          # If the any exon pair fails to pass the threshold then drop the transcript
          unless($support_count >= $transcript_weight_threshold) {
            $keep_transcript = 0;
          }
        }
      }
    }

    if($keep_transcript) {
      say "Keeping transcript";
      push(@{$output_transcripts},$transcript);
    }
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
          $transcript_redundancy->{$t2->{'internal_transcript_id'}} = 1;
        }
      }
      elsif(scalar(@{$t1_exons}) == scalar(@{$t2_exons})) {
        my $is_redundant = 0;
        if($t1->length >= $t2->length) {
          $is_redundant = $self->exon_subset($t1,$t2,$t1_exons,$t2_exons);
          if($is_redundant) {
            $transcript_redundancy->{$t2->{'internal_transcript_id'}} = 1;
          }
        } else {
          $is_redundant = $self->exon_subset($t2,$t1,$t2_exons,$t1_exons);
          if($is_redundant) {
            $transcript_redundancy->{$t1->{'internal_transcript_id'}} = 1;
          }
        }
      }
      else {
         my $is_redundant = $self->exon_subset($t2,$t1,$t2_exons,$t1_exons);
         if($is_redundant) {
          $transcript_redundancy->{$t1->{'internal_transcript_id'}} = 1;
        }
      }
    }
  }

  foreach my $transcript (@{$transcripts}) {
    unless($transcript_redundancy->{$transcript->{'internal_transcript_id'}}) {
      push(@{$final_transcripts},$transcript)
    }
  }

  return $final_transcripts;
}



sub compare_transcripts {
  my ($self,$transcripts) = @_;

  my $final_transcripts = [];
  my $transcript_redundancy = {};
  my $i=0;
  for($i=0; $i<scalar@{$transcripts}; $i++) {
    my $t1 = ${$transcripts}[$i];
    my $j = $i+1;
    for($j=$i+1; $j<scalar@{$transcripts}; $j++) {
      my $t2 = ${$transcripts}[$j];

      say "COMPARING: ".$t1->{'internal_transcript_id'}." TO ".$t2->{'internal_transcript_id'};
      my $unique_exons_t1 = $self->find_unique_exons($t1,$t2);
      my $unique_exons_t2 = $self->find_unique_exons($t2,$t1);

      # In this case all exons overlap, so select one of the two and mark the other as redundant
      unless($unique_exons_t1 || $unique_exons_t2) {
        my $redundant_internal_transcript_id = $self->choose_best_transcript($t1,$t2);
        $transcript_redundancy->{$redundant_internal_transcript_id} = 1;
      }
    }
  }

  foreach my $transcript (@{$transcripts}) {
    unless($transcript_redundancy->{$transcript->{'internal_transcript_id'}}) {
      push(@{$final_transcripts},$transcript)
    }
  }

  return $final_transcripts;
}


sub choose_best_transcript {
  my ($self,$transcript_a,$transcript_b) = @_;

  my $redundant_internal_transcript_id;

  my $transcript_a_cov = $transcript_a->{'cov'};
  my $transcript_b_cov = $transcript_b->{'cov'};

  my $transcript_a_pid = $transcript_a->{'pid'};
  my $transcript_b_pid = $transcript_a->{'pid'};

  my $transcript_a_support = $transcript_a_cov + $transcript_a_pid;
  my $transcript_b_support = $transcript_b_cov + $transcript_b_pid;

  my $transcript_a_weight_threshold = $transcript_a->{'transcript_weight_threshold'};
  my $transcript_b_weight_threshold = $transcript_b->{'transcript_weight_threshold'};

  # First pick the one with the lowest weight as this should correspond to the most reliable source
  # Then if the weights are the same pick the one with the best combined coverage and pid score
  if($transcript_a_weight_threshold < $transcript_b_weight_threshold) {
    $redundant_internal_transcript_id = $transcript_b->{'internal_transcript_id'};
    say "Removing transcript ".$transcript_b->{'internal_transcript_id'}."\nCov: ".$transcript_b_cov.
    "\nPid: ".$transcript_b_pid."\n\n in favour of transcript ".$transcript_a->{'internal_transcript_id'}.
    "\nCov: ".$transcript_a_cov."\nPid: ".$transcript_a_pid."\n";
  } elsif($transcript_a_weight_threshold > $transcript_b_weight_threshold) {
    $redundant_internal_transcript_id = $transcript_a->{'internal_transcript_id'};
    say "Removing transcript ".$transcript_a->{'internal_transcript_id'}."\nCov: ".$transcript_a_cov.
    "\nPid: ".$transcript_a_pid."\n\n in favour of transcript ".$transcript_b->{'internal_transcript_id'}.
    "\nCov: ".$transcript_b_cov."\nPid: ".$transcript_b_pid."\n";
  } elsif($transcript_a_weight_threshold == $transcript_b_weight_threshold) {
    if($transcript_a_support >= $transcript_b_support) {
      $redundant_internal_transcript_id = $transcript_b->{'internal_transcript_id'};
      say "Removing transcript ".$transcript_b->{'internal_transcript_id'}."\nCov: ".$transcript_b_cov.
          "\nPid: ".$transcript_b_pid."\n\n in favour of transcript ".$transcript_a->{'internal_transcript_id'}.
          "\nCov: ".$transcript_a_cov."\nPid: ".$transcript_a_pid."\n";
    } else {
      $redundant_internal_transcript_id = $transcript_a->{'internal_transcript_id'};
      say "Removing transcript ".$transcript_a->{'internal_transcript_id'}."\nCov: ".$transcript_a_cov.
          "\nPid: ".$transcript_a_pid."\n\n in favour of transcript ".$transcript_b->{'internal_transcript_id'}.
          "\nCov: ".$transcript_b_cov."\nPid: ".$transcript_b_pid."\n";
    }
  }

  return($redundant_internal_transcript_id);
}


sub features_overlap {
  my ($self,$feature_a,$feature_b) = @_;

  if (($feature_a->seq_region_start() <= $feature_b->seq_region_end()) &&
      ($feature_a->seq_region_end() >= $feature_b->seq_region_start())) {
    return 1;
  }

  return 0;
}

sub find_unique_exons {
  my ($self,$transcript_a,$transcript_b) = @_;

  foreach my $exon_a (@{$transcript_a->get_all_translateable_Exons()}) {
    my $overlap = 0;
    say "E_A: ".$exon_a->start().", ".$exon_a->end();
    foreach my $exon_b (@{$transcript_b->get_all_translateable_Exons()}) {
      say "E_B: ".$exon_b->start().", ".$exon_b->end();
      $overlap = $self->features_overlap($exon_a,$exon_b);
      if($overlap) {
        say "OVERLAP FOUND!!!!";
        say "OVERLAP INNER: ".$overlap;
        last;
      }
    }

    say "OVERLAP OUTER: ".$overlap;
    unless($overlap) {
      say "NO OVERLAP FOUND, UNIQUE EXON PRESENT";
      return(1);
    }
  }

  return 0;
}


sub exon_subset {
  my ($self,$transcript_a,$transcript_b,$exons_a,$exons_b) = @_;

  my $is_subset = 0;
  my $start_exon_b = ${$exons_b}[0];
  my $exon_match_count = 0;
  say "Transcript A: ".$transcript_a->analysis->logic_name().", ".$transcript_a->{'internal_transcript_id'};
  say "Transcript B: ".$transcript_b->analysis->logic_name().", ".$transcript_b->{'internal_transcript_id'};

  my $i=0;
  for($i=0; $i<scalar(@{$exons_a}); $i++) {
    if(${$exons_a}[$i]->seq_region_start == $start_exon_b->seq_region_start &&
       ${$exons_a}[$i]->seq_region_end == $start_exon_b->seq_region_end) {
      $exon_match_count++;

      my $j=1;
      for($j=1; $j<scalar(@{$exons_b}) && ($i+$j)<scalar(@{$exons_a}); $j++) {
        if(${$exons_a}[$i+$j]->seq_region_start == ${$exons_b}[$j]->seq_region_start &&
          ${$exons_a}[$i+$j]->seq_region_end == ${$exons_b}[$j]->seq_region_end) {
          $exon_match_count++;
        }
      }
    }
  }

  if($exon_match_count == scalar(@{$exons_b})) {
    say "Model ".$transcript_b->{'internal_transcript_id'}." is redundant to model ".$transcript_a->{'internal_transcript_id'};
    $is_subset = 1;
  }
  return $is_subset;
}


sub recluster_genes {
  my ($self,$processed_genes) = @_;

  my $single_transcript_genes = [];
  my $output_genes = [];
  foreach my $processed_gene (@{$processed_genes}) {
    my $transcripts = $processed_gene->get_all_Transcripts();

    # If the gene has a single transcript then there is no need to recluster
    if(scalar(@{$transcripts}) == 1) {
      push(@{$single_transcript_genes},$processed_gene);
      next;
    }

    # If there is more than one transcript, then make a new gene for each
    elsif(scalar(@{$transcripts}) > 1) {
      foreach my $transcript (@{$transcripts}) {
        my $single_transcript_gene = Bio::EnsEMBL::Gene->new();
        $single_transcript_gene->add_Transcript($transcript);
        push(@{$single_transcript_genes},$single_transcript_gene);
      }
    }

    else {
      $self->throw("Got a gene without any transcripts");
    }
  }


  my $biotypes_hash = $self->get_all_biotypes($single_transcript_genes);
  my $biotypes_array = [keys(%$biotypes_hash)];
  my $types_hash;
  $types_hash->{genes} = $biotypes_array;

  say "Reclustering gene models based on final transcript set";
  my ($clusters, $unclustered) = cluster_Genes($single_transcript_genes,$types_hash);
  say "...finished reclustering genes models";
  say "Found clustered sets: ".scalar(@{$clusters});
  say "Found unclustered sets: ".scalar(@{$unclustered});

  say "Processing gene clusters into single gene models...";
  $output_genes = $self->process_clusters($clusters,$unclustered);
  say "...finished processing gene clusters into single gene models";
  say "Made output genes: ".scalar(@{$output_genes});

  $self->output_genes($output_genes);

}


sub process_clusters {

  my ($self,$clustered,$unclustered) = @_;

  my $output_genes = [];
  foreach my $single_cluster (@{$unclustered}) {
    my $cluster_genes = $single_cluster->get_Genes();
    foreach my $single_gene (@{$cluster_genes}) {
      my $output_gene = Bio::EnsEMBL::Gene->new();
      $output_gene->slice($self->query());
      $output_gene->analysis($self->analysis);
      $output_gene->biotype($self->analysis->logic_name);
      my $transcripts = $single_gene->get_all_Transcripts();
      my $single_transcript = shift(@{$transcripts});
      $output_gene->add_Transcript($single_transcript);
      push(@{$output_genes},$output_gene);
    }
  }

  foreach my $single_cluster (@{$clustered}) {
    my $combined_gene = Bio::EnsEMBL::Gene->new();
    $combined_gene->slice($self->query());
    $combined_gene->analysis($self->analysis);
    $combined_gene->biotype($self->analysis->logic_name);

    my $cluster_genes = $single_cluster->get_Genes();

    foreach my $single_gene (@{$cluster_genes}) {
      my $transcripts = $single_gene->get_all_Transcripts();
      my $single_transcript = shift(@{$transcripts});
      $combined_gene->add_Transcript($single_transcript);
    }
    push(@{$output_genes},$combined_gene);
  }

  return($output_genes);

}

sub load_cluster {
  my ($self,$cluster) = @_;

  my $gene = new Bio::EnsEMBL::Gene();
  my $dna_dba = $self->hrdb_get_con('dna_db');
  my $input_dbs = $self->param('input_gene_dbs');

  # Use the var below to give each transcript an internal id, since they are coming from multiple input dbs
  my $internal_transcript_id = 0;
  foreach my $adaptor_name (keys(%{$cluster})) {
    unless($input_dbs->{$adaptor_name}) {
      $self->throw("You are using a cluster type input id, but one of the the adaptor names in the input id did not match a corresponding db hash".
                   "All adaptor names must have a correspondingly named db hash passed in. Offending adaptor name:\n".$adaptor_name);
    }
    my $dba = $self->hrdb_get_dba($input_dbs->{$adaptor_name});
    unless($dba) {
      $self->throw("You are using a cluster type input id, but one of the the adaptor names in the input id did not match a corresponding db hash".
                   "All adaptor names must have a correspondingly named db hash passed in. Offending adaptor name:\n".$adaptor_name);
    }
    $dba->dnadb($dna_dba);
    $self->hrdb_set_con($dba,$adaptor_name);
    my $transcript_adaptor = $dba->get_TranscriptAdaptor();
    my $transcript_id_array = $cluster->{$adaptor_name};
    foreach my $transcript_id (@{$transcript_id_array}) {
      my $transcript = $transcript_adaptor->fetch_by_dbID($transcript_id);
      my $transcript_sf = $transcript->get_all_supporting_features();
      my $transcript_cov = ${$transcript_sf}[0]->hcoverage();
      my $transcript_pid = ${$transcript_sf}[0]->percent_id();

      # These are some internal values that get used elsewhere for convenience. The internal id is needed to avoid id conflicts
      $transcript->{'cov'} = $transcript_cov;
      $transcript->{'pid'} = $transcript_pid;
      $transcript->{'internal_transcript_id'} = $internal_transcript_id;
      $internal_transcript_id++;
      $gene->add_Transcript($transcript);
    }
  }

  say "Created a new gene, transcript count: ".scalar(@{$gene->get_all_Transcripts()});
  return $gene;
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
      say "Skipping gene as it will be picked up on the adjoining another slice";
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


sub logic_name_weights {
  my ($self,$val) = @_;

  if($val) {
    $self->param('_logic_name_weights',$val);
  }

  return($self->param('_logic_name_weights'));
}


sub single_exon_support_penalty {
  my ($self,$val) = @_;

  if($val) {
    $self->param('_single_exon_support_penalty',$val);
  }

  return($self->param('_single_exon_support_penalty'));
}


sub find_best_remaining_transcript_old {
  my ($self,$gene,$exon_pairs) = @_;

  my $best_transcripts = [];
  my $max_support_score = 0;
  my $min_allowed_score = $self->param('min_backup_score');
  unless(defined($min_allowed_score)) {
    $min_allowed_score = 10;
  }

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


sub find_best_remaining_transcript {
  my ($self,$gene,$exon_pairs) = @_;

  my $lowest_weight_transcripts = [];
  my $max_support_score = 0;
  my $lowest_weight_threshold = 999;
  my $min_allowed_score = $self->param('min_backup_score');
  my $max_allowed_weight = $self->param('max_backup_weight');
  unless(defined($min_allowed_score)) {
    $min_allowed_score = 10;
  }
  unless(defined($max_allowed_weight)) {
    $max_allowed_weight = 10;
  }

  my $transcripts = $gene->get_all_Transcripts();
  foreach my $transcript (@{$transcripts}) {
    my $support_score = $self->calculate_support_average($transcript,$exon_pairs);
    unless($support_score >= $min_allowed_score) {
      next;
    }

    my $transcript_weight_threshold = $transcript->{'transcript_weight_threshold'};
    if(($transcript_weight_threshold <= $lowest_weight_threshold) && ($transcript_weight_threshold <= $max_allowed_weight)) {
      $lowest_weight_threshold = $transcript_weight_threshold;
      push(@{$lowest_weight_transcripts},$transcript);
    }
  }

  my $best_length = 0;
  my $best_combined_cov_and_pid = 0;
  my $final_transcript;
  foreach my $transcript (@{$lowest_weight_transcripts}) {
    my $cov = $transcript->{'cov'};
    my $pid = $transcript->{'pid'};
    my $combined_cov_and_pid = $cov + $pid;
    if($combined_cov_and_pid > $best_combined_cov_and_pid) {
      $best_combined_cov_and_pid = $combined_cov_and_pid;
      $best_length = length($transcript->translate());
      $final_transcript = $transcript;
    } elsif($combined_cov_and_pid == $best_combined_cov_and_pid) {
      if(length($transcript->translate() > $best_length)) {
        $best_length = length($transcript->translate());
        $final_transcript = $transcript;
      }
    }
  }

  if($final_transcript) {
    $final_transcript->biotype($final_transcript->biotype."_br");
    return [$final_transcript];
  }

  return [];
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
      if($exon_pair->{$transcript->{'internal_transcript_id'}}) {
        $total_exon_pairs++;
        $total_support_amount += scalar(keys(%$exon_pair));
      }
    }

    my $average_support = $total_support_amount / $total_exon_pairs;
    $final_support = $average_support + ($total_exon_pairs * $exon_pair_count_bonus);
  } # end else

  return($final_support);
}


sub get_all_biotypes {
  my ($self,$master_genes_array) = @_;

  my $master_biotypes_hash = {};

  foreach my $gene (@{$master_genes_array}) {
    unless($master_biotypes_hash->{$gene->biotype}) {
      $master_biotypes_hash->{$gene->biotype} = 1;
    }
  }
  return($master_biotypes_hash);
}


1;
