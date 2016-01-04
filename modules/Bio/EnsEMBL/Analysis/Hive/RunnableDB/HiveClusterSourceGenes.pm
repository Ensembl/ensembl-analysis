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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClusterSourceGenes;

use strict;
use warnings;
use feature 'say';
use Data::Dumper;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);
use Bio::EnsEMBL::DBSQL::SliceAdaptor;


use vars qw (@EXPORT);

use Exporter 'import';

@EXPORT = qw(
fetch_source_genes
get_all_genes
get_all_biotypes
process_clusters
master_genes_hash
);

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
  $out_dba->dnadb($dba);
  $self->hrdb_set_con($out_dba,'output_db');

  my $input_id = $self->param('iid');

  unless($self->param('input_gene_dbs')) {
    $self->throw("You must specify an arrayref of input dbs via the 'input_gene_dbs' parameter");
  }

  my $slice;
  my $cluster_strand;
  if($self->param('iid_type') eq 'slice') {
    $slice = $self->fetch_sequence($input_id,$dba);
  } elsif($self->param('iid_type') eq 'cluster_slice') {
    unless(split(':',$input_id) == 7) {
      $self->throw("You have selected to use the 'cluster_slice' input id type, but it was not possible to parse ".
                   "the input_id correctly. Cluster slice input ids are similar to slice ids, but with an additional ".
                   "column at the end to define the strand of the genes the cluster represents, e.g.:\n".
                   "\"chromosome:PapAnu2.0:18:808824:882888:1:-1\"\nOffending input id:\n".$input_id);
    }

    unless($input_id =~ /^(.+)\:(\-?1)$/) {
      $self->throw("You have selected to use the 'cluster_slice'!!!!!!!!!!!!, but it was not possible to parse ".
                   "the input_id correctly. Cluster slice input ids are similar to slice ids, but with an additional ".
                      "column at the end to define the strand of the genes the cluster represents, e.g.:\n".
                   "\"chromosome:PapAnu2.0:18:808824:882888:1:-1\"\nOffending input id:\n".$input_id);
    }

    my $slice_name = $1;
    $cluster_strand = $2;
    $slice = $self->fetch_sequence($slice_name,$dba);
  } else {
    $self->throw("You must specify an input_id type in the config using the 'iid_type' parameter");
  }

  # Set the slice at this point
  $self->query($slice);


  my $input_genes = $self->fetch_source_genes($input_gene_dbs,$allowed_input_sets);
  if($self->param('iid_type') eq 'cluster_slice') {
    $input_genes = $self->filter_by_strand($input_genes,$cluster_strand);
  }

  $master_genes_hash->{'input_genes'} = $input_genes;
  $self->master_genes_hash($master_genes_hash);

  # Now need to think how this would be represented in config
  # input_gene_dbs => ['rnaseq_blast_db','genblast_db','exonerate_db'],
  # allowed_gene_sets => {'baboon_liver_models' => {'protein_coding_80_100' => 1},
  #                       'baboon_lung_models' => {'protein_coding_80_100' => 1,
  #                                                'protein_coding_50_80' => 1}}
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

  my $test_transcripts = $$master_genes_array[0]->get_all_Transcripts;
  my $test_transcript = $$test_transcripts[0];
      say "TTRANSCRIPT START: ".$test_transcript->start;
      say "TTRANSCRIPT END: ".$test_transcript->end;
      say "TTRANSCRIPT STRAND: ".$test_transcript->strand;
      say "TTRANSCRIPT SLICE: ".$test_transcript->slice->name;
      my $exons = $test_transcript->get_all_Exons;
      say "TT EXONS: ".scalar(@{$exons});

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

#  say "FM2 DUMPER OUTGENES: ".Dumper($output_genes);
  $self->output_genes($output_genes);

  return 1;
}

sub write_output {
  my $self = shift;

  my $adaptor = $self->hrdb_get_con('output_db')->get_GeneAdaptor;
  my @output = @{$self->output_genes};

  say "Writing genes to output db";
  foreach my $gene (@output){
    $gene->slice($self->query);
#    say "GENE START: ".$gene->start;
#    say "GENE END: ".$gene->end;
 #   say "GENE STRAND: ".$gene->strand;
#    say "GENE SLICE: ".$gene->slice->name;
  #  $gene->slice($self->query->seq_region_Slice);

#    my $transcript = ${$gene->get_all_Transcripts}[0];
#    say "TRANSCRIPT START: ".$transcript->seq_region_start;
#    say "TRANSCRIPT END: ".$transcript->seq_region_end;
#    say "TRANSCRIPT STRAND: ".$transcript->seq_region_strand;
#    say "TRANSCRIPT SLICE: ".$transcript->slice->name;
#    my $exons = $transcript->get_all_Exons;

#    foreach my $exon (@{$exons}) {
 #     say "EXON START: ".$exon->seq_region_start;
 #     say "EXON END: ".$exon->seq_region_end;
#      say "EXON STRAND: ".$exon->seq_region_strand;
#      say "EXON SLICE: ".$exon->slice->name;
#    }

    empty_Gene($gene);
    $adaptor->store($gene);
  }
  say "...finished writing genes to output db";

  return 1;
}

sub fetch_source_genes {
  my ($self,$gene_source_dbs,$allowed_gene_sets) = @_;

  my $final_filtered_gene_set = [];
  my $slice = $self->query;

  foreach my $db_con_hash (@{$gene_source_dbs}) {
    say "Getting con for: ".$db_con_hash->{'-dbname'};
    my $db_adaptor = $self->hrdb_get_dba($db_con_hash);
    my $gene_adaptor = $db_adaptor->get_GeneAdaptor();
    my $all_genes = $gene_adaptor->fetch_all_by_Slice($slice);
    my $filtered_genes = [];
    if($allowed_gene_sets) {
      foreach my $gene (@{$all_genes}) {
        my $logic_name = $gene->analysis->logic_name;
        if($allowed_gene_sets->{$logic_name}) {
          if(ref($allowed_gene_sets->{$logic_name}) eq 'HASH') {
            my $biotype = $gene->biotype;
            if($allowed_gene_sets->{$logic_name}->{$biotype}) {
              push(@{$filtered_genes},$gene);
            }
          } elsif($allowed_gene_sets->{$logic_name}) {
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


sub filter_by_strand {
  my ($self,$unfiltered_genes,$strand) = @_;

  my $filtered_genes = [];
  my $slice = $self->query();
  my $cluster_start = $slice->start;
  my $cluster_end = $slice->end;
  foreach my $gene (@{$unfiltered_genes}) {
    # If the gene is not completely within the cluster boundary or is on the wrong slice then skip it
    if($gene->seq_region_start < $cluster_start || $gene->seq_region_end > $cluster_end || $gene->strand != $strand) {
      next;
    }
    push(@{$filtered_genes},$gene);
  }
  return($filtered_genes);
}


sub process_clusters {

  my ($self,$clustered,$unclustered) = @_;


  my $all_clusters = [@{$clustered},@{$unclustered}];

  my $output_genes = [];
#  foreach my $single_cluster (@{$unclustered}) {
#    my $cluster_genes = $single_cluster->get_Genes();
#    foreach my $single_gene (@{$cluster_genes}) {
#      my $output_gene = Bio::EnsEMBL::Gene->new();
#      $output_gene->adaptor($self->hrdb_get_con('dna_db'));
#      $output_gene->analysis($self->analysis);
#      $output_gene->biotype($self->analysis->logic_name);
#      my $transcripts = $single_gene->get_all_Transcripts();
#      my $single_transcript = shift(@{$transcripts});
#      my $exons = $single_transcript->get_all_Exons();

#      say "TRANSCRIPT START: ".$single_transcript->start;
#      say "TRANSCRIPT END: ".$single_transcript->end;
#      say "TRANSCRIPT STRAND: ".$single_transcript->strand;
#      say "TRANSCRIPT SLICE: ".$single_transcript->slice->name;
#      say "FM2 DUMPER TSLICE: ".Dumper($single_transcript->slice);
#      say "EXON COUNT: ".scalar(@{$exons});
#      my $slice = $self->query;
#      $single_transcript->flush_Exons();
#      foreach my $exon (@{$exons}) {
#        $exon->slice($slice);
#        $single_transcript->add_Exon($exon);
#      }

#      $output_gene->add_Transcript($single_transcript);
#      push(@{$output_genes},$output_gene);
#    }
#  }

  foreach my $single_cluster (@{$all_clusters}) {
    my $combined_gene = Bio::EnsEMBL::Gene->new();
    $combined_gene->analysis($self->analysis);
    $combined_gene->biotype($self->analysis->logic_name);

    my $cluster_genes = $single_cluster->get_Genes();

    foreach my $single_gene (@{$cluster_genes}) {
      my $transcripts = $single_gene->get_all_Transcripts();
      foreach my $transcript (@{$transcripts}) {
        my $exons = $transcript->get_all_Exons();
        my $slice = $self->query;
        $transcript->flush_Exons();
        foreach my $exon (@{$exons}) {
          $exon->slice($slice);
          $transcript->add_Exon($exon);
        }
        $combined_gene->add_Transcript($transcript);
      }
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
