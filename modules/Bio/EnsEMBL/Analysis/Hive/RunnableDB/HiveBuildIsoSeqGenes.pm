#!/usr/bin/env perl

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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBuildIsoSeqGenes;

use strict;
use warnings;
use feature 'say';
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);
use Bio::EnsEMBL::DBSQL::SliceAdaptor;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;

  $self->create_analysis;

  my $input_gene_dbs =  $self->param('input_gene_dbs');
  my $allowed_input_sets = $self->param('allowed_input_sets');

  my $dba = $self->hrdb_get_dba($self->param('dna_db'));
  $self->hrdb_set_con($dba,'dna_db');


  my $output_dba = $self->hrdb_get_dba($self->param('target_db'));
  $self->hrdb_set_con($output_dba,'target_db');

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

  say "Reclustering unclustered set on stranded genomic overlap to remove problematic genes";

  my $clustered_genes = $self->extract_genes($clusters);
  my $unclustered_genes = $self->extract_genes($unclustered);


  my ($second_clusters, $second_unclustered) = cluster_Genes($unclustered_genes,$types_hash,0,1);
  say "...finished reclustering genes";

  say "Found second clustered sets: ".scalar(@{$second_clusters});
  say "Found second unclustered sets: ".scalar(@{$second_unclustered});

  my $second_unclustered_genes = $self->extract_genes($second_unclustered);
  foreach my $gene (@$second_unclustered_genes) {
    $gene->biotype('third_unclustered');
  }

  $types_hash->{third_unclustered} = ['third_unclustered'];
  say "Doing final round on remaining unclustered set to see if they fall within any of the original clusters";
  my ($third_clusters, $third_unclustered) = cluster_Genes([@$clustered_genes,@$second_unclustered_genes],$types_hash,0,1);

  my $final_unclustered_genes = $self->extract_genes($third_unclustered,'third_unclustered');
  say "Found final unclustered sets: ".scalar(@{$final_unclustered_genes});

  my $final_clustered_genes = $self->build_genes($clusters);
  my $output_genes = [@$final_clustered_genes,@$final_unclustered_genes];
  $self->output($output_genes);
  return 1;
}


sub write_output {
  my $self = shift;

  my $adaptor = $self->hrdb_get_con('target_db')->get_GeneAdaptor;

  foreach my $gene (@{$self->output}) {
    $gene->analysis($self->analysis);
    $gene->biotype('isoseq');
    empty_Gene($gene);
    $adaptor->store($gene);
  }

  return 1;
}

sub filter_input_genes {
  my ($self,$gene_source_dbs,$allowed_transcript_sets) = @_;

  my $min_transcript_length = 500;
  my $max_intron_size = 75000;
  my $min_exons = 4;
  my $biotype = 'isoseq';

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
            unless($transcript->length() >= $min_transcript_length) {
              next;
  	    }
            unless($self->check_introns($transcript,$max_intron_size)) {
              next;
            }
            unless(scalar(@$transcript->get_all_Exons()) >= $min_exons) {
              next;
	    }

            if(ref($allowed_transcript_sets->{$logic_name}) eq 'HASH') {
              my $biotype = $transcript->biotype;
              if($allowed_transcript_sets->{$logic_name}->{$biotype}) {
                my $single_transcript_gene = new Bio::EnsEMBL::Gene;
                $transcript->{'adaptor_name'} = $adaptor_name;
                $single_transcript_gene->{'adaptor_name'} = $adaptor_name;
                $transcript->biotype($biotype);
                $single_transcript_gene->add_Transcript($transcript);
                push(@{$final_input_genes},$single_transcript_gene);
              }
            } elsif($allowed_transcript_sets->{$logic_name}) {
              my $single_transcript_gene = new Bio::EnsEMBL::Gene;
              $transcript->{'adaptor_name'} = $adaptor_name;
              $single_transcript_gene->{'adaptor_name'} = $adaptor_name;
              $transcript->biotype($biotype);
              $single_transcript_gene->add_Transcript($transcript);
              push(@{$final_input_genes},$single_transcript_gene);
            }
          }
        }

       # Should maybe use scalar($input_genes) == 0 here to sanity check that something was retrieved from the db
    } else {
      foreach my $transcript (@{$all_transcripts}) {
        unless($transcript->length() >= $min_transcript_length) {
          next;
        }
        unless($self->check_introns($transcript,$max_intron_size)) {
          next;
        }
        unless(scalar(@{$transcript->get_all_Exons()}) >= $min_exons) {
          next;
	}

        my $single_transcript_gene = new Bio::EnsEMBL::Gene;
        $transcript->{'adaptor_name'} = $adaptor_name;
        $single_transcript_gene->{'adaptor_name'} = $adaptor_name;
        $transcript->biotype($biotype);
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



sub check_introns {
  my ($self,$transcript,$max_intron_size) = @_;

  my $introns = $transcript->get_all_Introns();
  foreach my $intron (@$introns) {
    if($intron->length() > $max_intron_size) {
      return(0);
    }
  }

  return(1);
}


sub extract_genes {
  my ($self,$clusters,$type) = @_;
  my $genes = [];
  foreach my $cluster (@$clusters) {
    if($type) {
      push(@$genes,@{$cluster->get_Genes_of_Type($type)});
    } else {
      push(@$genes,@{$cluster->get_Genes});
    }
  }
  return($genes);
}

sub build_genes {
  my ($self,$clusters) = @_;

  my $final_genes = [];
  foreach my $cluster (@$clusters) {
    my $genes = $cluster->get_Genes();
    my $gene_string_hash = {};
    my $final_gene = new Bio::EnsEMBL::Gene;
    foreach my $gene (@$genes) {
      my $gene_string = $self->gene_string($gene);
      if($gene_string_hash->{$gene_string}) {
        say "Skipping gene/transcript because it was already seen";
       next;
      } else {
        $gene_string_hash->{$gene_string} = 1;
        my $transcripts = $gene->get_all_Transcripts();
        foreach my $transcript (@$transcripts) {
          $final_gene->add_Transcript($transcript);
        }
      }
    } # end foreach my $gene
    push(@$final_genes,$final_gene);
  } # end foreach my $cluster

  return($final_genes);
}

sub gene_string {
  my ($self,$gene) = @_;

  my $gene_string = $gene->biotype.":".$gene->start.":".$gene->end.":".$gene->strand.":".$gene->seq_region_name;
  my $exons = $gene->get_all_Exons();
  $gene_string .= $self->generate_exon_string($exons);

  return($gene_string);
}

sub generate_exon_string {
  my ($self,$exon_array) = @_;

  my $exon_string = "";
  foreach my $exon (@{$exon_array}) {
    my $start = $exon->start();
    my $end = $exon->end();
    $exon_string .= $start."..".$end.":";
  }

  return($exon_string);
}

1;

