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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClusterSourceGenes;

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

  my $master_genes_hash = {};

  my $input_gene_dbs =  $self->param('input_gene_dbs');
  my $allowed_input_sets = $self->param('allowed_input_sets');

  my $dba = $self->hrdb_get_dba($self->param('dna_db'));
  $self->hrdb_set_con($dba,'dna_db');

  my $out_dba = $self->hrdb_get_dba($self->param('output_db'));
  $self->hrdb_set_con($out_dba,'output_db');

  my $input_id = $self->param('iid');
  my $slice = $self->fetch_sequence($input_id,$dba);
  $self->query($slice);

  if($input_gene_dbs) {
    my $input_genes = $self->fetch_source_genes($input_gene_dbs,$allowed_input_sets);
    $master_genes_hash->{'input_genes'} = $input_genes;
  }

  $self->master_genes_hash($master_genes_hash);

  # Now need to think how this would be represented in config
  # input_gene_dbs => ['rnaseq_blast_db','genblast_db','exonerate_db'],
  # rnaseq_allowed_sets => {'baboon_liver_models' => {'protein_coding_80_100' => 1},
  #                         'baboon_lung_models' => {'protein_coding_80_100' => 1,
  #                         'protein_coding_50_80' => 1}}
  #
  # e.g unless($rnaseq_allowed_biotypes->{$transcript->logic_name}->{$transcript->biotype})
  #
  # I've made the above felixed, but for autonmation the min expected input will just be
  # the source dbs and it will be assumed that everything in them should be taken


  return 1;
}

sub run {
  my $self = shift;

  my $master_genes_array = $self->get_all_genes();

  my $master_biotypes_hash = $self->get_all_biotypes($master_genes_array);
  my $master_biotypes_array = [keys(%$master_biotypes_hash)];
  my $types_hash;
  $types_hash->{genes} = $master_biotypes_array;

  say "Clustering genes from input_dbs...";
  my ($clusters, $unclustered) = cluster_Genes($master_genes_array,$types_hash);
  say "...finished clustering genes from input_dbs";
  say "Found clustered sets: ".scalar(@{$clusters});
  say "Found unclustered sets: ".scalar(@{$unclustered});

  say "Processing gene clusters into single gene models...";
  my $output_genes = $self->process_clusters($clusters,$unclustered);
  say "...finished processing gene clusters into single gene models";
  say "Made output genes: ".scalar(@{$output_genes});

  $self->output_genes($output_genes);

  return 1;
}

sub write_output {
  my $self = shift;

  my $adaptor = $self->hrdb_get_con('output_db')->get_GeneAdaptor;
  my @output = @{$self->output_genes};

  say "Writing genes to output db";
  foreach my $gene (@output){
    empty_Gene($gene);
    $adaptor->store($gene);
  }
  say "...finished writing genes to output db";

  return 1;
}

sub fetch_source_genes {
  my ($self,$gene_source_dbs,$gene_allowed_sets) = @_;

  my $final_filtered_gene_set = [];
  my $slice = $self->query;

  foreach my $db_con_hash (@{$gene_source_dbs}) {
    say "Getting con for: ".$db_con_hash->{'-dbname'};
    my $db_adaptor = $self->hrdb_get_dba($db_con_hash);
    my $gene_adaptor = $db_adaptor->get_GeneAdaptor();
    my $all_genes = $gene_adaptor->fetch_all_by_Slice($slice);
    my $filtered_genes = [];
    if($gene_allowed_sets) {
      foreach my $gene (@{$all_genes}) {
        my $logic_name = $gene->analysis->logic_name;
        if($gene_allowed_sets->{$logic_name}) {
          if(ref($gene_allowed_sets->{$logic_name}) eq 'HASH') {
            my $biotype = $gene->biotype;
            if($gene_allowed_sets->{$logic_name}->{$biotype}) {
              push(@{$filtered_genes},$gene);
            }
          } elsif($gene_allowed_sets->{$logic_name}) {
            push(@{$filtered_genes},$gene);
          }
        }
      }
      my $set_gene_count = scalar(@{$filtered_genes});
      say "Set genes fetched: ".$set_gene_count;
      push(@{$final_filtered_gene_set},@{$filtered_genes});
    } else {
      my $set_gene_count = scalar(@{$all_genes});
      say "Set genes fetched: ".$set_gene_count;
      push(@{$final_filtered_gene_set},@{$all_genes});
    }
  }

  my $final_gene_count = scalar(@{$final_filtered_gene_set});
  say "Total genes fetched: ".$final_gene_count;
  return $final_filtered_gene_set;
}


sub get_all_genes {
  my ($self) = @_;

  my $master_genes_hash = $self->master_genes_hash();
  my $master_genes_array = [];
  foreach my $genes_set (keys(%$master_genes_hash)) {
    push(@{$master_genes_array},@{$master_genes_hash->{$genes_set}});
  }

  return($master_genes_array);
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

sub master_genes_hash {
  my ($self,$val) = @_;

  if($val) {
    $self->param('_master_gene_hash',$val);
  }

  return($self->param('_master_gene_hash'));
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

sub output_genes {
  my ($self,$val) = @_;

  if($val) {
    $self->param('_output_genes',$val);
  }

  return($self->param('_output_genes'));
}

1;
