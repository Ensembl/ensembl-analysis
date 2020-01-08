# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

#!/usr/bin/env perl

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveTranscriptExtractor;

use warnings;
use strict;
use feature 'say';

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub fetch_input {
  my $self = shift;


  my $input_dba = $self->hrdb_get_dba($self->param('input_db'));
  $self->hrdb_set_con($input_dba,'input_db');


  my $output_dba = $self->hrdb_get_dba($self->param('output_db'));
  $self->hrdb_set_con($output_dba,'output_db');

  $self->fetch_input_genes($self->input_id);
  say "Starting gene count: ".scalar(@{$self->input_genes});

  $self->allowed_biotypes($self->param('allowed_biotypes'));
  return 1;
}


sub run {
  my $self = shift;

  if($self->input_genes && scalar(@{$self->input_genes}) > 0) {
    my $output_genes = $self->build_single_transcript_genes();
    $self->output_genes($output_genes);
  } else {
    $self->warning("No genes found for the following input_id:\n".$self->input_id);
  }
  return 1;
}


sub write_output {
  my $self = shift;

  # Only output for slices that have valid transcripts
  $self->input_job->autoflow(0);

  my $output_gene_adaptor = $self->hrdb_get_con('output_db')->get_GeneAdaptor();

  if($self->output_genes && scalar(@{$self->output_genes}) > 0) {
    foreach my $gene (@{$self->output_genes}) {

      empty_Gene($gene);
      $output_gene_adaptor->store($gene);

      my $output_hash = {};
      $output_hash->{'iid'} = $gene->stable_id();
      $self->dataflow_output_id($output_hash,1);
    }

  }

  say "Output gene count: ".scalar(@{$self->output_genes});
  return 1;
}


sub allowed_biotypes {
  my ($self,$biotypes) = @_;

  if($biotypes) {
    $self->param('_allowed_biotypes',$biotypes);
  }

  return($self->param('_allowed_biotypes'));

}


sub fetch_input_genes {
  my ($self,$slice_name) = @_;

  my $sa = $self->hrdb_get_con('input_db')->get_SliceAdaptor();
  my $slice = $sa->fetch_by_name($slice_name);

  $self->input_genes($slice->get_all_Genes());
}


sub input_genes {
  my ($self,$input_genes) = @_;

  if($input_genes) {
    $self->param('_input_genes',$input_genes);
  }

  return($self->param('_input_genes'));

}


sub build_single_transcript_genes {
  my $self = shift;

  my $output_genes = [];
  foreach my $gene (@{$self->input_genes}) {
    foreach my $transcript (@{$gene->get_all_Transcripts()}) {

      unless(scalar(@{$transcript->get_all_Attributes('gencode_basic')})) {
        next;
      }

      unless($transcript->stable_id()) {
        $self->warning("Transcript must have a stable id for the next step, skipping transcript");
        next;
      }

      unless($self->allowed_biotypes->{$transcript->biotype}) {
        next;
      }

      my $new_gene = Bio::EnsEMBL::Gene->new(
              -STABLE_ID => $transcript->stable_id(),
              -START     => $transcript->start(),
              -END       => $transcript->end(),
              -STRAND    => $transcript->strand(),
              -SLICE     => $transcript->slice(),
              -ANALYSIS  => $transcript->analysis(),
      );

      $new_gene->add_Transcript($transcript);
      push(@{$output_genes},$new_gene);
    }
  }
  return $output_genes;
}


sub output_genes {
  my ($self,$output_genes) = @_;

  unless($self->param_is_defined('_output_genes')) {
    $self->param('_output_genes',[]);
  }

  if($output_genes) {
    $self->param('_output_genes',$output_genes);
  }

  return $self->param('_output_genes');
}

1;
