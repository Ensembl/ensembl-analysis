#!/usr/bin/env perl

# Copyright [2017-2024] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCollapseIGTR;

use strict;
use warnings;
use feature 'say';

use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene attach_Analysis_to_Gene);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');



sub fetch_input {
  my $self = shift;

  $self->setup_fasta_db;
  my $dba = $self->get_database_by_name('target_db');
  $self->hrdb_set_con($dba,'target_db');

}

sub run {
  my $self = shift;
  $self->process_genes();
}

sub write_output {
  my $self = shift;
  my $output_genes = $self->output();
  my $analysis = new Bio::EnsEMBL::Analysis(-logic_name => 'ig_tr_collapse');

  my $ga = $self->hrdb_get_con('target_db')->get_GeneAdaptor();
  foreach my $gene (@{$output_genes}) {
    empty_Gene($gene);
    attach_Analysis_to_Gene($gene,$analysis);
    $ga->store($gene);
  }
}


sub process_genes {
  my ($self) = @_;

  my $genes = [];
  my $dba = $self->hrdb_get_con('target_db');
  my $ga = $dba->get_GeneAdaptor();
  my $logic_names = $self->param_required('logic_names_to_cluster');
  foreach my $logic_name (@{$logic_names}) {
    my $current_genes = $ga->fetch_all_by_logic_name($logic_name);
    push(@{$genes},@{$current_genes});
  }

  my $biotypes_hash = $self->get_all_biotypes($genes);
  my $biotypes_array = [keys(%$biotypes_hash)];
  my $types_hash;
  $types_hash->{genes} = $biotypes_array;
  my ($clustered, $unclustered) = cluster_Genes($genes,$types_hash);
  my $output_genes = $self->process_clusters($clustered,$unclustered);
  $self->output($output_genes);
}

sub process_clusters {
  my ($self,$clustered,$unclustered) = @_;

  my $output_genes = [];
  my $all_clusters = [@{$clustered},@{$unclustered}];
  say "Total clusters: ".scalar(@{$clustered});
  say "Total unclustered: ".scalar(@{$unclustered});

  foreach my $cluster (@{$all_clusters}) {
    my $cluster_genes = $cluster->get_Genes();
    my $best_combined_cov_pid = 0;
    my $best_gene;
    foreach my $single_gene (@{$cluster_genes}) {
      say "Cluster gene 1";
      my $transcript = pop(@{$single_gene->get_all_Transcripts});
      say "Transcript: ".$transcript->dbID();
      my $tsf = pop(@{$transcript->get_all_supporting_features});
      unless($tsf) {
        next;
      }

      my $combined_cov_pid = $tsf->hcoverage + $tsf->percent_id;
      if($combined_cov_pid > $best_combined_cov_pid) {
        $best_gene = $single_gene;
        $best_combined_cov_pid = $combined_cov_pid;
      }
    }
    push(@{$output_genes},$best_gene);
  }

  return($output_genes);
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
