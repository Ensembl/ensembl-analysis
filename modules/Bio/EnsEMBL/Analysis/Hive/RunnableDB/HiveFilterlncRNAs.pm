#!/usr/bin/env perl

# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveFilterlncRNAs;

use strict;
use warnings;
use feature 'say';
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene attach_Slice_to_Gene);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(calculate_exon_phases);
use Bio::EnsEMBL::Analysis::Runnable::GeneBuilder;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');
use Data::Dumper;

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    input_biotype  => 'pre_lncRNA',
    output_biotype => 'lncRNA',
    # lncRNAs can be kept based on passing either min sequence or min exon cutoffs
    min_sequence_size => 200,
    min_exons => 3,
    min_sequence_if_less_min_exons => 1000,

    # Genebuilder options for collapse
    max_transcripts_per_cluster => 5,
    min_short_intron_len => 7,
    max_short_intron_len => 15,
    blessed_biotypes => {},
    coding_only => $self->param('coding_only'),
    skip_readthrough_check => 1,
  }
}

sub fetch_input {
  my $self = shift;

  $self->setup_fasta_db;
  my $test_case = 0;
  my $analysis = Bio::EnsEMBL::Analysis->new(
                                              -logic_name => $self->param('logic_name'),
                                              -module => $self->param('module'),
                                            );
  $self->analysis($analysis);

  my $input_dba = $self->get_database_by_name('input_gene_db');
  $self->hrdb_set_con($input_dba,'input_gene_db');
}

sub run {
  my ($self) = @_;

  my $input_dba = $self->hrdb_get_con('input_gene_db');
  my $slice_adaptor = $input_dba->get_SliceAdaptor();
  my $slices = $slice_adaptor->fetch_all('toplevel');
  my $final_genes = [];
  my $genes_to_remove = [];

  foreach my $slice (@$slices) {
    say "Slice: ".$slice->name;
    my $genes = $slice->get_all_Genes_by_type($self->param_required('input_biotype'));
    unless($genes) {
      next;
    }
    $self->query($slice);
    push(@$genes_to_remove,@$genes);
    say "Initial lncRNA genes: ".scalar(@$genes);
    say "Filtering intial genes on exons, sequence length and protein coding overlap";
    my $filtered_genes = $self->filter_lncrnas($genes);
    say "Genes remaining after filtering: ".scalar(@$filtered_genes);
    say "Collapsing remaining genes:";
    my $collapsed_genes = $self->collapse_remaining_lncrnas($filtered_genes);
    say "Genes remaining after collapse: ".scalar(@$collapsed_genes);
    my $lncrna_genes = $self->remove_translations($collapsed_genes);
    push(@$final_genes,@$lncrna_genes);
#     push(@$final_genes,@$collapsed_genes);
  }
  $self->output($final_genes);
  $self->original_genes_to_remove($genes_to_remove);
}

sub write_output {
  my ($self) = @_;

  my $input_dba = $self->hrdb_get_con('input_gene_db');
  my $gene_adaptor = $input_dba->get_GeneAdaptor();
  my $output_genes = $self->output;
  foreach my $output_gene (@$output_genes) {
    $output_gene->analysis($self->analysis());
    $output_gene->biotype($self->param_required('output_biotype'));
    my $transcripts = $output_gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts) {
      $transcript->biotype($self->param('output_biotype'));
    }
    empty_Gene($output_gene);
    $gene_adaptor->store($output_gene);
  }

  my $genes_to_remove = $gene_adaptor->fetch_all_by_biotype($self->param('input_biotype'));
  foreach my $gene_to_remove (@$genes_to_remove) {
     $gene_adaptor->remove($gene_to_remove);
  }
}

sub load_genes {
  my ($self) = @_;

  my $dba = $self->hrdb_get_con('input_gene_db');
  my $gene_adaptor = $dba->get_GeneAdaptor();
  my $genes = $gene_adaptor->fetch_all_by_biotype($self->param('input_biotype'));
  $self->genes_to_qc($genes);
}


sub original_genes_to_remove {
  my ($self,$genes) = @_;
  if($genes) {
    $self->param('_genes_to_remove',$genes);
  }
  return($self->param('_genes_to_remove'));
}


sub genes_to_qc {
  my ($self,$genes) = @_;
  if($genes) {
    $self->param('_genes_to_qc',$genes);
  }

  return($self->param('_genes_to_qc'));
}

sub filter_lncrnas {
  my ($self,$genes) = @_;

  my $filtered_genes = [];
  say "Have ".scalar(@{$genes})." genes to process";
  my $dba = $self->hrdb_get_con('input_gene_db');
  my $gene_adaptor = $dba->get_GeneAdaptor();
  foreach my $gene (@{$genes}) {
    my $remove = 0;
    # First look for any overlapping same stranded protein coding genes, if there is any overlap then we remove the gene
    my $overlapping_genes = $gene->get_overlapping_Genes(1);
    foreach my $overlapping_gene (@$overlapping_genes) {
      if($overlapping_gene->biotype eq 'protein_coding') {
        say "Removing gene due to stranded overlap";
        $remove = 1;
        last;
      }
    }

    # This is skipped if we have found overlap
    unless($remove) {
      # In here we will loop through the transcripts, if any of the transcripts pass the checks then we will just keep the gene
      $remove = 1;
      my $transcripts = $gene->get_all_Transcripts();
      foreach my $transcript (@{$transcripts}) {
        unless($self->remove_gene($transcript)) {
          $remove = 0;
        }
      }
    }

    unless($remove) {
      push(@$filtered_genes,$gene);
    }
  } # End foreach my $gene (@{$genes})

  return($filtered_genes);
}


sub collapse_remaining_lncrnas {
  my ($self,$genes)  = @_;

  my $collapsed_genes = [];
  my $types_hash;
  $types_hash->{'remaining_lncrnas'} = [$self->param('input_biotype')];


  say "Clustering remaining lncRNA genes genes...";
  my ($clusters, $unclustered) = cluster_Genes($genes,$types_hash);
  say "...finished clustering genes";

  say "Found ".scalar(@$clusters)." clusters";
  say "Found ".scalar(@$unclustered)." unclustered genes, will not change";

  foreach my $unclustered (@$unclustered) {
    my $unclustered_genes = $unclustered->get_Genes_by_Set('remaining_lncrnas');
    foreach my $unclustered_gene (@$unclustered_genes) {
      push(@$collapsed_genes,$unclustered_gene);
    }
  }

  foreach my $cluster (@$clusters) {
    my $lncrna_cluster_genes = $cluster->get_Genes_by_Set('remaining_lncrnas');
    if(scalar(@$lncrna_cluster_genes)) {
      say "Inputting ".scalar(@$lncrna_cluster_genes)." single transcript genes to new genebuilder run";
      say "Running genebuilder...";
      my $runnable = Bio::EnsEMBL::Analysis::Runnable::GeneBuilder->new(
                      -query => $self->query,
                      -analysis => $self->analysis,
                      -genes => $lncrna_cluster_genes,
                      -output_biotype => $self->param('output_biotype'),
                      -max_transcripts_per_cluster => $self->param('max_transcripts_per_cluster'),
                      -min_short_intron_len => $self->param('min_short_intron_len'),
                      -max_short_intron_len => $self->param('max_short_intron_len'),
                      -blessed_biotypes => $self->param('blessed_biotypes'),
                      -coding_only => $self->param('coding_only'),
                      -skip_readthrough_check => $self->param('skip_readthrough_check'),
                     );
      $runnable->run;
      say "Created ".scalar(@{$runnable->output})." output genes";
      push(@$collapsed_genes,@{$runnable->output});
    } else {
      say "Not lncRNA clusters to process on cluster genes";
    }
  }

  return($collapsed_genes);
}


sub remove_gene {
  my ($self,$transcript) = @_;

  my $min_sequence_if_less_min_exons = $self->param_required('min_sequence_if_less_min_exons');
  my $min_exon_count = $self->param_required('min_exons');
  my $exons = $transcript->get_all_Exons;

  # If the transcript is below the smallest allowed size for a lncRNA (currently 200bp) then mark as bad
  if($transcript->length < $self->param_required('min_sequence_size')) {
    return(1);
  } elsif(scalar(@$exons) >= $min_exon_count) {
    return(0);
  } elsif($transcript->length >= $min_sequence_if_less_min_exons) {
    return(0);
  } else {
    return(1);
  }
}


sub remove_translations {
  my ($self,$genes) = @_;

  my $input_dba = $self->hrdb_get_con('input_gene_db');
  my $translation_adaptor = $input_dba->get_TranslationAdaptor;

  foreach my $gene (@$genes) {
    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts) {
      if($transcript->translation) {
#         $translation_adaptor->remove($transcript->translation);
         $transcript->translation(undef);
         $transcript->adaptor(undef);
      }

      my $exons = $transcript->get_all_Exons();
      foreach my $exon (@{$exons}) {
        $exon->phase(-1);
        $exon->end_phase(-1);
      }
    }
  }
  return($genes);
}


sub calculate_coverage {
  my ($self, $transcript_a, $transcript_b) = @_;

  my $coverage = 0;
  my $overlap = 0;
  foreach my $exon_a (@{$transcript_a->get_all_Exons()}) {
    foreach my $exon_b (@{$transcript_b->get_all_Exons()}) {
      $overlap += overlap_length($exon_a, $exon_b);
    }
  }

  $coverage = ($overlap / $transcript_b->length) * 100;

  return $coverage;
}


sub overlap_length {
  my ($feature_a, $feature_b) = @_;

  if (!features_overlap($feature_a, $feature_b)) {
    return 0;
  }

  my $min_end = $feature_a->seq_region_end();
  if ($feature_b->seq_region_end() < $min_end) {
    $min_end = $feature_b->seq_region_end();
  }

  my $max_start = $feature_a->seq_region_start();
  if ($feature_b->seq_region_start() > $max_start) {
    $max_start = $feature_b->seq_region_start();
  }

  return $min_end - $max_start + 1;
}


sub features_overlap {
  my ($feature_a, $feature_b) = @_;

  if(($feature_a->seq_region_start() <= $feature_b->seq_region_end() ) && ($feature_a->seq_region_end() >= $feature_b->seq_region_start())) {
    return 1;
  }

  return 0;
}

1;
