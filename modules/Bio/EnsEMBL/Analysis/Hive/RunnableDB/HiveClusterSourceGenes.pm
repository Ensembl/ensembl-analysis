#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);
use Bio::EnsEMBL::DBSQL::SliceAdaptor;

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

  my $input_gene_dbs =  $self->param('input_gene_dbs');
  my $allowed_input_sets = $self->param('allowed_input_sets');

  my $dba = $self->hrdb_get_dba($self->param('dna_db'));
  $self->hrdb_set_con($dba,'dna_db');

  my $input_id = $self->param('iid');

  unless($self->param('input_gene_dbs')) {
    $self->throw("You must specify an arrayref of input dbs via the 'input_gene_dbs' parameter");
  }

  my $slice;
  if($self->param('iid_type') eq 'slice') {
    $slice = $self->fetch_sequence($input_id,$dba);
    $self->query($slice);
  } else {
    $self->throw("You must specify an input_id type in the config using the 'iid_type' parameter");
  }

  # Set the slice at this point
  $self->query($slice);

  say "Fetching input genes...";
  my $final_input_genes = $self->filter_input_genes($input_gene_dbs,$allowed_input_sets);
  say "finished fetching input genes";
  $self->input_genes($final_input_genes);
  return 1;
}

sub run {
  my $self = shift;

  my $input_genes = $self->input_genes();
  my $biotypes_hash = $self->get_all_biotypes($input_genes);
  my $biotypes_array = [keys(%$biotypes_hash)];

  my $types_hash;
  $types_hash->{genes} = $biotypes_array;

  say "Clustering genes from input_dbs...";
  my ($clusters, $unclustered) = cluster_Genes($input_genes,$types_hash);
  say "...finished clustering genes from input_dbs";
  say "Found clustered sets: ".scalar(@{$clusters});
  say "Found unclustered sets: ".scalar(@{$unclustered});

  say "Processing gene clusters into single gene models...";
  my $output_ids = $self->output_ids($clusters,$unclustered);
  say "...finished processing gene clusters into single gene models";

  $self->output($output_ids);

  return 1;
}


sub write_output {
  my $self = shift;

  my $output_ids = $self->output();

  unless(scalar(@{$output_ids})) {
    $self->warning("No input ids generated for this analysis!");
  }

  foreach my $output_id (@{$output_ids}) {

    # Yet to implement this, will allow for slices that don't have any of the specified features to be skipped
    if($self->param('check_slice_for_features')) {
      unless($self->check_slice_for_features()) {
        next;
      }
    }

    my $output_hash = {};
    $output_hash->{'iid'} = $output_id;
    $self->dataflow_output_id($output_hash,2);
  }

  return 1;

}

sub filter_input_genes {
  my ($self,$gene_source_dbs,$allowed_transcript_sets) = @_;

  my $final_input_genes = [];
  my $slice = $self->query;

  foreach my $adaptor_name (keys(%{$gene_source_dbs})) {
    say "  Found adaptor: ".$adaptor_name;
    my $db_con_hash = $gene_source_dbs->{$adaptor_name};
    my $db_adaptor = $self->hrdb_get_dba($db_con_hash);
    my $transcript_adaptor = $db_adaptor->get_TranscriptAdaptor();
    my $all_transcripts = $transcript_adaptor->fetch_all_by_Slice($slice);
    my $input_genes = [];

    if($allowed_transcript_sets) {
        foreach my $transcript (@{$all_transcripts}) {
          my $logic_name = $transcript->analysis->logic_name;
          if($allowed_transcript_sets->{$logic_name}) {
            if(ref($allowed_transcript_sets->{$logic_name}) eq 'HASH') {
              my $biotype = $transcript->biotype;
              if($allowed_transcript_sets->{$logic_name}->{$biotype}) {
                my $single_transcript_gene = new Bio::EnsEMBL::Gene;
                $transcript->{'adaptor_name'} = $adaptor_name;
                $single_transcript_gene->{'adaptor_name'} = $adaptor_name;
                $single_transcript_gene->add_Transcript($transcript);
                push(@{$final_input_genes},$single_transcript_gene);
              }
            } elsif($allowed_transcript_sets->{$logic_name}) {
              my $single_transcript_gene = new Bio::EnsEMBL::Gene;
              $transcript->{'adaptor_name'} = $adaptor_name;
              $single_transcript_gene->{'adaptor_name'} = $adaptor_name;
              $single_transcript_gene->add_Transcript($transcript);
              push(@{$final_input_genes},$single_transcript_gene);
            }
          }
        }

       # Should maybe use scalar($input_genes) == 0 here to sanity check that something was retrieved from the db
    } else {
      foreach my $transcript (@{$all_transcripts}) {
        my $single_transcript_gene = new Bio::EnsEMBL::Gene;
        $transcript->{'adaptor_name'} = $adaptor_name;
        $single_transcript_gene->{'adaptor_name'} = $adaptor_name;
        $single_transcript_gene->add_Transcript($transcript);
        push(@{$final_input_genes},$single_transcript_gene);
      }

      # Should maybe use scalar($input_genes) == 0 here to sanity check that something was retrieved from the db
    }

  }

  return $final_input_genes;
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


sub output_ids {

  my ($self,$clustered,$unclustered) = @_;

  my $all_clusters = [@{$clustered},@{$unclustered}];

  my $output_ids = [];
  foreach my $single_cluster (@{$all_clusters}) {
    my $cluster_genes = $single_cluster->get_Genes();
    my $cluster_output_id = {};
    foreach my $single_gene (@{$cluster_genes}) {
      my $transcripts = $single_gene->get_all_Transcripts();
      my $single_transcript = shift(@{$transcripts});

      my $adaptor_name = $single_transcript->{'adaptor_name'};
      unless($cluster_output_id->{$adaptor_name}) {
        $cluster_output_id->{$adaptor_name} = [];
      }
      push(@{$cluster_output_id->{$adaptor_name}},$single_transcript->dbID());
    }
    push(@{$output_ids},$cluster_output_id);
  }

  return($output_ids);

}


sub input_genes {
  my ($self,$val) = @_;

  if($val) {
    $self->param('_input_genes',$val);
  }

  return($self->param('_input_genes'));
}


1;
