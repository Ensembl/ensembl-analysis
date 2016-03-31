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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveFindIntergenicRegions;

use strict;
use warnings;
use feature 'say';
use Data::Dumper;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);
use Bio::EnsEMBL::DBSQL::SliceAdaptor;


use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub fetch_input {
  my $self = shift;

  my $input_gene_dbs =  $self->param('input_gene_dbs');

  my $dba = $self->hrdb_get_dba($self->param('dna_db'));
  $self->hrdb_set_con($dba,'dna_db');

  my $input_id = $self->param('iid');

  unless($self->param('input_gene_dbs')) {
    $self->throw("You must specify an arrayref of input dbs via the 'input_gene_dbs' parameter");
  }

  my $slice;
  my $cluster_strand;
  if($self->param('iid_type') eq 'slice') {
    $slice = $self->fetch_sequence($input_id,$dba);
  } else {
    $self->throw("You must specify an input_id type in the config using the 'iid_type' parameter");
  }

  # Set the slice at this point
  $self->query($slice);

  my $input_gene_coords = $self->calculate_input_gene_coords($input_gene_dbs);
  $self->input_gene_coords($input_gene_coords);

  return 1;
}

sub run {
  my $self = shift;

  my $input_gene_coords = $self->input_gene_coords();

  say "Processing input genes to find intergenic regions...";
  my $output_ids = $self->split_slice_on_intergenic($input_gene_coords);
  say "...finished finding intergenic regions";

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


sub calculate_input_gene_coords {
  my ($self,$gene_source_dbs) = @_;

  my @sorted_start_end_coords = ();
  my @unsorted_start_end_coords = ();
  my $slice = $self->query;
  my $skip_duplicates = {};
  foreach my $adaptor_name (keys(%{$gene_source_dbs})) {
    my $db_con_hash = $gene_source_dbs->{$adaptor_name};
    my $db_adaptor = $self->hrdb_get_dba($db_con_hash);
    my $gene_adaptor = $db_adaptor->get_GeneAdaptor();
    my $all_genes = $gene_adaptor->fetch_all_by_Slice($slice);

    foreach my $gene (@{$all_genes}) {
      my $start = $gene->start;
      my $end = $gene->end;
      if($skip_duplicates->{$start.":".$end}) {
        next;
      } else {
        push(@unsorted_start_end_coords,{'s' => $start, 'e' => $end});
        $skip_duplicates->{$start.":".$end} = 1;
      }
    }
  }

  @sorted_start_end_coords = sort {
                                    $a->{'s'} <=> $b->{'s'} ||
                                    $a->{'e'} <=> $b->{'e'}
                                  } @unsorted_start_end_coords;

  return \@sorted_start_end_coords;
}


sub split_slice_on_intergenic {
  my ($self,$input_gene_coords) = @_;

  my $slice = $self->query;
  my @slice_eles = split(':',$slice->name());

  my $output_ids = [];
  my $in_overlap = 0;
  my $overlap_start = 0;
  my $overlap_end = 0;
  for(my $i=0; $i<scalar(@{$input_gene_coords}); $i++) {
    my $right_start = ${$input_gene_coords}[$i]->{'s'};
    my $right_end = ${$input_gene_coords}[$i]->{'e'};
    if($i==0) {
      $overlap_start = $right_start;
      $overlap_end = $right_end;
      next;
    }

    my $left_start = ${$input_gene_coords}[$i-1]->{'s'};
    my $left_end = ${$input_gene_coords}[$i-1]->{'e'};

    # If there is no overlap, print out the previous region
    if($right_start > $overlap_end) {
      say "Overlap: ".$overlap_start."..".$overlap_end;
      $slice_eles[3] = $overlap_start;
      $slice_eles[4] = $overlap_end;
      push(@{$output_ids},join(':',@slice_eles));
      $overlap_start = $right_start;
      $overlap_end = $right_end;
    } else {
      if($right_end > $overlap_end) {
        $overlap_end = $right_end;
      }
    }
  }

  $slice_eles[3] = $overlap_start;
  $slice_eles[4] = $overlap_end;
  push(@{$output_ids},join(':',@slice_eles));

  return($output_ids);
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


sub input_gene_coords {
  my ($self,$val) = @_;

  if($val) {
    $self->param('_input_gene_coords',$val);
  }

  return($self->param('_input_gene_coords'));
}

1;
