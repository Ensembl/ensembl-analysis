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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveUTRAddition;

use strict;
use warnings;
use feature 'say';
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;


  my $analysis = Bio::EnsEMBL::Analysis->new(
                                              -logic_name => $self->param('logic_name'),
                                              -module => $self->param('module'),
                                            );
  $self->analysis($analysis);

  my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
  $self->hrdb_set_con($dna_dba,'dna_db');

#  my $no_utr_dba = $self->hrdb_get_dba($self->param('no_utr_db'));
#  $no_utr_dba->dnadb($dna_dba);
#  $self->hrdb_set_con($no_utr_dba,'no_utr_input_db');

#  my $utr_source_dba = $self->hrdb_get_dba($self->param('utr_source_db'));
#  $utr_source_dba->dnadb($dna_dba);
#  $self->hrdb_set_con($utr_source_dba,'utr_source_db');

  my $utr_output_dba = $self->hrdb_get_dba($self->param('utr_output_db'));
  $utr_output_dba->dnadb($dna_dba);
  $self->hrdb_set_con($utr_output_dba,'utr_output_db');

  my $input_id = $self->param('iid');
  my $input_id_type = $self->param('iid_type');

  # If the input is an array of gene ids then loop through them and fetch them from the input db
  if($input_id_type eq 'slice') {
     $self->load_cluster($input_id);
  } else {
    $self->throw("You have selected an input id type using iid_type that is not supported by the module. Offending type: ".$input_id_type);
  }

  return 1;
}

sub run {
  my $self = shift;

  my $utr_transcripts = $self->utr_transcripts();
  my $no_utr_transcripts = $self->no_utr_transcripts();

  unless(scalar(@{$utr_transcripts})) {
    say "Found no overlapping transcripts from the utr transcript db";
    return 1;
  }

  unless(scalar(@{$no_utr_transcripts})) {
    say "Found no overlapping transcripts from the no utr transcript db";
    return 1;
  }

  say "Attempting to add utr to non utr transcripts";
#  $self->add_utr($utr_transcripts,$no_utr_transcripts);
  my $rough_gene = $self->transcript_merge($utr_transcripts,$no_utr_transcripts);
  my $final_genes = $self->recluster_genes([$rough_gene]);

  $self->output_genes($final_genes);

  return 1;
}

sub write_output {
  my $self = shift;

  my $adaptor = $self->hrdb_get_con('utr_output_db')->get_GeneAdaptor;
  my @output = @{$self->output_genes};

  say "Writing genes to output db";
  foreach my $gene (@output){
    empty_Gene($gene);
    $adaptor->store($gene);
  }
  say "...finished writing genes to output db";

  return 1;
}


sub load_cluster {
  my ($self,$cluster) = @_;

  my $gene = new Bio::EnsEMBL::Gene();
  my $dna_dba = $self->hrdb_get_con('dna_db');
  my $input_dbs = $self->param('input_gene_dbs');

  my $utr_transcripts = [];
  my $no_utr_transcripts = [];

  # Use the var below to give each transcript an internal id, since they are coming from multiple input dbs
  my $internal_transcript_id = 0;
  foreach my $adaptor_name (keys(%{$cluster})) {
    unless($input_dbs->{$adaptor_name}) {
      $self->throw("You are using a cluster type input id, but one of the the adaptor names in the input id did not match a corresponding db hash\n".
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

      $transcript->{'internal_transcript_id'} = $internal_transcript_id;
      $internal_transcript_id++;

      if($adaptor_name eq 'no_utr_db') {
        push(@{$no_utr_transcripts},$transcript);
      } else {
        push(@{$utr_transcripts},$transcript);
      }
    }
  }


    say "utr transcript count: ".scalar(@{$utr_transcripts});
    $self->utr_transcripts($utr_transcripts);

    say "no utr transcript count: ".scalar(@{$no_utr_transcripts});
    $self->no_utr_transcripts($no_utr_transcripts);

}

sub add_utr {
  my ($self,$utr_transcripts,$no_utr_transcripts) = @_;

  foreach my $no_utr_transcript (@{$no_utr_transcripts}) {
    foreach my $utr_transcript (@{$utr_transcripts}) {
      my $five_prime_utrs = $utr_transcript->get_all_five_prime_UTRs();
      my $three_prime_utrs = $utr_transcript->get_all_three_prime_UTRs();
      unless(scalar(@{$five_prime_utrs}) || scalar(@{$three_prime_utrs})) {
        say "Found no UTR on either 5 prime or 3 prime ends of utr transcript, so skipping future comparison";
        next;
      }

      my $no_utr_transcript_cds = $no_utr_transcript->get_all_CDS();
      my $utr_transcript_cds = $utr_transcript->get_all_CDS();
      if($self->check_cds_introns($utr_transcript_cds,$no_utr_transcript_cds)) {
 #       if(scalar(@{$five_prime_utrs})) {
 #         $self->add_five_prime_utr($five_prime_utrs,$utr_transcript_cds,$no_utr_transcript_cds);
 #       }
 #        if(scalar(@{$three_prime_utrs})) {
 #          $self->add_three_prime_utr($three_prime_utrs,$utr_transcript_cds,$no_utr_transcript_cds);
 #       }

        $self->merge($utr_transcript,$no_utr_transcript);
      }
    }
  }


}


sub transcript_merge {
  my ($self,$utr_transcripts,$no_utr_transcripts) = @_;


  my $output_gene = Bio::EnsEMBL::Gene->new();
  $output_gene->analysis($self->analysis);
  $output_gene->biotype($self->analysis->logic_name);

  my $final_transcripts = [];
  foreach my $no_utr_transcript (@{$no_utr_transcripts}) {
    foreach my $utr_transcript (@{$utr_transcripts}) {
      my $no_utr_transcript_cds = $no_utr_transcript->get_all_CDS();
      my $utr_transcript_cds = $utr_transcript->get_all_CDS();
      if($self->check_cds_introns($utr_transcript_cds,$no_utr_transcript_cds)) {
        my $final_transcript = $self->merge_supporting_features($utr_transcript,$no_utr_transcript);
        $final_transcript->biotype($no_utr_transcript->biotype());
        push(@{$final_transcripts},$final_transcript);
        last;
      } else {
        push(@{$final_transcripts},$no_utr_transcript);
      }
    }
  }

  foreach my $final_transcript (@{$final_transcripts}) {
    $output_gene->add_Transcript($final_transcript);
  }
  return($output_gene);
}


sub add_three_prime_utr {
  my ($self,$utrs,$cds_a,$cds_b) = @_;

  my $final_exons = $cds_b;
  my $terminal_exon_a = $$cds_a[scalar(@{$cds_a})-1];
  my $terminal_exon_b = $$cds_b[scalar(@{$cds_b})-1];

  my $first_utr = $$utrs[0];
  say "FIRST UTR START: ".$first_utr->start();
  say "TERMINAL EXON A END: ".$terminal_exon_a->end();
  say "TERMINAL EXON B END: ".$terminal_exon_b->end();
}


sub add_five_prime_utr {
  my ($self,$utrs,$cds_a,$cds_b) = @_;

  my $final_exons = $cds_b;
  my $terminal_exon_a = $$cds_a[0];
  my $terminal_exon_b = $$cds_b[0];

  my $first_utr = $$utrs[scalar(@{$utrs})-1];
  say "LAST UTR END: ".$first_utr->end();
  say "TERMINAL EXON A START: ".$terminal_exon_a->start();
  say "TERMINAL EXON B START: ".$terminal_exon_b->start();
}


sub check_cds_introns {
  my ($self,$cds_a,$cds_b) = @_;

  my $exon_count_a = scalar(@{$cds_a});
  my $exon_count_b = scalar(@{$cds_b});

  unless($exon_count_a && $exon_count_b) {
    return 0;
  }

  unless($exon_count_a == $exon_count_b) {
    return 0;
  }

  if($exon_count_a == 1) {
    my $exon_a = shift(@{$cds_a});
    my $exon_b = shift(@{$cds_b});
    if($self->features_overlap($exon_a,$exon_b)) {
      return 1;
    } else {
      return 0;
    }
  } else {
    my $intron_coords_a = $self->intron_coords($cds_a);
    my $intron_coords_b = $self->intron_coords($cds_b);

    if($intron_coords_a eq $intron_coords_b) {
      return 1;
    } else {
      return 0;
    }
  }

}

sub features_overlap {
  my ( $featureA, $featureB ) = @_;

  if (
     ( $featureA->seq_region_start() <= $featureB->seq_region_end() ) &&
     ( $featureA->seq_region_end() >= $featureB->seq_region_start() ) )
    {
    return 1;
  }

  return 0;
}

sub features_overlap_new {
  my ($self,$feature_a,$feature_b) = @_;

  if(($feature_a->seq_region_start() <= $feature_b->seq_region_end() ) && ($feature_a->seq_region_end() >= $feature_b->seq_region_start())) {
    return 1;
  }

  return 0;
}

sub intron_coords {
  my ($self,$exons) = @_;

  my $intron_coords = "";
  for(my $i=0; $i<(scalar@{$exons})-1; $i++) {
    my $exon_a = $$exons[$i];
    my $exon_b = $$exons[$i+1];
    $intron_coords .= $exon_a->end()."..".$exon_b->start().":";
  }

  unless($intron_coords) {
    $self->throw("Issue with finding the intron coordinates!");
  }

  return($intron_coords);
}


sub merge_supporting_features {
  my ($self,$transcript_a,$transcript_b) = @_;

  # Copy all transcript related features from $source_transcript into
  # $target_transcript:
  #     supporting features
  #     intron supporting evidence
  {
    my @supporting_features = ();

    # only transfer the supporting features that overlap.
    # as the merge is based on intron match, there can be cases where
    # a longer Secondary transcript evidence would have been transferred
    # beyond the exon boundaries of the Primary target transcript
    foreach my $sf (@{$transcript_b->get_all_supporting_features()}) {
      if (features_overlap($sf,$transcript_a)) {
        push(@supporting_features,$sf);
      }
    }

    $transcript_a->add_supporting_features(@supporting_features);

    say "Merge: transferred transcript supporting features";

    my @intron_support = @{ $transcript_b->get_all_IntronSupportingEvidence() };

    foreach my $intron_support (@intron_support) {
      $transcript_a->add_IntronSupportingEvidence($intron_support);
    }

    say "Merge: transferred transcript supporting features";

  }

  # Transfer all exon related features from $source_transcript into
  # $target_transcript:
  #     supporting features
  {
    my @supporting_features;

    foreach my $source_exon (@{$transcript_b->get_all_Exons()} ) {
      my @exon_sf = ();
      # only transfer the supporting features that overlap.
      # as the merge is based on intron match, there can be cases where
      # a longer Secondary transcript evidence would have been transferred
      # beyond the exon boundaries of the Primary target transcript
      foreach my $sf (@{ $source_exon->get_all_supporting_features() }) {
        if (features_overlap($sf,$transcript_a)) {
          push(@exon_sf,$sf);
        }
      }
      push(@supporting_features,[@exon_sf]);
    }

    my $exon_index    = 0;
    my $feature_count = 0;
    foreach my $target_exon (@{$transcript_a->get_all_Exons()}) {
      $feature_count += scalar( @{$supporting_features[$exon_index]});
      $target_exon->add_supporting_features(@{$supporting_features[$exon_index++]});
#      add_logic_name_suffix( $target_exon, 'merged' );
    }

    say "Merge: transferred exon supporting evidence";

  }

  $transcript_a->source($transcript_a->source().'_'.$transcript_b->source());
  return($transcript_a);
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


sub utr_transcripts {
  my ($self,$val) = @_;

  if($val) {
    $self->param('_utr_transcripts',$val);
  }

  return($self->param('_utr_transcripts'));
}

sub no_utr_transcripts {
  my ($self,$val) = @_;

  if($val) {
    $self->param('_no_utr_transcripts',$val);
  }

  return($self->param('_no_utr_transcripts'));
}


sub output_genes {
  my ($self,$val) = @_;

  if($val) {
    $self->param('_output_genes',$val);
  }

  return($self->param('_output_genes'));
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


1;
