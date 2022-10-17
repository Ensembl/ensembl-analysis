=head1 LICENSE

Copyright [2021] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::Minimap2Remap

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::CollapseParalogues;

use warnings;
use strict;
use feature 'say';
use Data::Dumper;
use File::Spec::Functions qw(tmpdir catfile);

use Bio::EnsEMBL::Analysis::Runnable::GeneBuilder;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name align_nucleotide_seqs);
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);
use Bio::EnsEMBL::Analysis::Runnable::Minimap2Remap;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene clone_Gene);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my($self) = @_;
  $self->create_analysis;

  my $target_dba = $self->hrdb_get_dba($self->param('target_db'));

  my $dna_dba;
  if($self->param('use_genome_flatfile')) {
    unless($self->param_required('genome_file') && -e $self->param('genome_file')) {
      $self->throw("You selected to use a flatfile to fetch the genome seq, but did not find the flatfile. Path provided:\n".$self->param('genome_file'));
    }
    setup_fasta(
                 -FASTA => $self->param_required('genome_file'),
               );
  } else {
    $dna_dba = $self->hrdb_get_dba($self->param('target_dna_db'));
    $target_dba->dnadb($dna_dba);
  }

  $self->hrdb_set_con($target_dba,'target_db');


  my $all_new_genes = [];
  my $slice_adaptor = $target_dba->get_SliceAdaptor();
  my $slices = $slice_adaptor->fetch_all('toplevel');
  foreach my $slice (@$slices) {
    my $initial_genes = $slice->get_all_Genes();
    my $new_genes = [];
    foreach my $input_gene (@$initial_genes) {
      my $input_gene_description = $input_gene->description();
      if ($input_gene_description and $input_gene_description =~ /potential_paralogue/) {
        if (scalar(@{$input_gene->get_all_Transcripts()})) {
          push(@$new_genes,$input_gene);
          push(@$all_new_genes,$input_gene);
        }
      }
    }
    $self->make_runnables($new_genes,$slice);
  }

  $self->param('input_genes',$all_new_genes);
}


sub run {
  my ($self) = @_;

#  foreach my $runnable (@{$self->runnable()}) {
#    $runnable->run();
#    my $output_genes = $runnable->output();
#    my @parent_gene_ids = @{$runnable->{'parent_gene_ids'}};

    # For some reason there seems to be an issue with the same id getting added multiple times
    # and this is happening for every new paralogue
    # Should debug at some point, but for now it's fine to just make unique here
#    my $id_seen = {};
#    my $id_string = "";
#    foreach my $id (@parent_gene_ids) {
#      unless($id_seen->{$id}) {
#        $id_string .= $id.":";
#        $id_seen->{$id} = 1;
#      }
#    }
#    $id_string =~ s/\:$//;
#    foreach my $output_gene (@$output_genes) {
#      my $gene_description = $output_gene->description().";parent_gene=".$id_string.";mapping_type=potential_paralogue";
#      $output_gene->description($gene_description);
#      # add string of source genes stable ids as gene attribute
#      my $parent_attribute = Bio::EnsEMBL::Attribute->new(-CODE => 'proj_parent_g',-VALUE => $id_string);
#      $output_gene->add_Attributes($parent_attribute);
#      $self->output([$output_gene]);
#    }
#  }
}


sub write_output {
  my ($self) = @_;
  my $output_dba = $self->hrdb_get_con('target_db');
  my $output_gene_adaptor = $output_dba->get_GeneAdaptor();
  my $output_genes = $self->output();

  my $current_genes = $output_gene_adaptor->fetch_all();
  foreach my $output_gene (@$output_genes) {
    # Just cloning for safety there as the next bit will remove the existing genes from the db
    $output_gene->analysis($self->analysis());
    my $cloned_output_gene = clone_Gene($output_gene);
    empty_Gene($cloned_output_gene);
    $output_gene_adaptor->store($cloned_output_gene);
  }

  # Remove the initial set of paralogues prior to collapsing
  my $original_new_genes = $self->param('input_genes');
  foreach my $current_gene (@$original_new_genes) {
    $output_gene_adaptor->remove($current_gene);
  }
}


sub make_runnables {
  my ($self,$genes,$slice) = @_;

  my $biotypes_hash = $self->generate_biotypes_hash($genes);

  say "Processing ".scalar(@$genes)." input genes for slice: ".$slice->seq_region_name();
  say "Clustering new genes with  on parent slice...";
  my ($clusters, $unclustered) = cluster_Genes($genes,$biotypes_hash);

  foreach my $cluster (@$clusters) {
    my $clustered_genes = $cluster->get_Genes();
    my $output_gene = $self->simple_layer_genes($clustered_genes);
    $self->output([$output_gene]);

    # NOTE!: Disabiling this for the moment as there seem to be genes getting through that have multiple parents with incompatible biotypes
    # Instead we'll just take the parent biotype group with a basic heirarchy and then take the transcript with the most exons for that
    # group

#    my $layered_clustered_genes = $self->layer_cluster_genes($clustered_genes);
#    my $parent_gene_ids = [];

    # There is a slight chance that this cluster has genes from multiple parents, so track this
#    foreach my $gene (@$clustered_genes) {
#      my $description = $gene->description();
#      $description =~ /;parent_gene=([^;]+);/;
#      my $parent_stable_id = $1;
#      push(@$parent_gene_ids,$parent_stable_id)
#    }

#    my $output_biotype = ${$clustered_genes}[0]->biotype();
#    my $runnable = Bio::EnsEMBL::Analysis::Runnable::GeneBuilder->new(
#                     -query => $slice,
#                     -analysis => $self->analysis(),
#                     -genes => $clustered_genes,
#                     -output_biotype => $output_biotype,
#                     -max_transcripts_per_cluster => 1000,
#                     -min_short_intron_len => 0,
#                     -max_short_intron_len => 15,
#                     -blessed_biotypes => {},
#                     -skip_readthrough_check => 1,
#                     -max_exon_length => 50000,
#                     -coding_only => 1,
#                   );
#    $runnable->{'parent_gene_ids'} = $parent_gene_ids;
#    $self->runnable($runnable);
  }

  # For unclustered stuff there's no need to build a runnable, so just put directly on the output
  foreach my $unclustered (@$unclustered) {
    my $unclustered_genes = $unclustered->get_Genes();
    foreach my $unclustered_gene (@$unclustered_genes) {
      $self->output([$unclustered_gene]);
    }
  }
}


sub simple_layer_genes {
  my ($self,$genes) = @_;

  # This is just a simple filtering to make sure we don't make any weird complex clusters like sncRNA transcripts with protein coding ones
  # It does mean that some potentially interesting edge cases might be missed, but this is fine for now
  my $biotype_group_priorities = { 'coding'     => 5,
                                   'lnoncoding' => 4,
                                   'pseudogene' => 3,
                                   'snoncoding' => 2,
                                   'mnoncoding' => 1,
                                   'undefined'  => 0,
                                 };

  my $selected_gene;
  my $selected_transcript;
  my $max_priority = 0;
  foreach my $gene (@$genes) {
    my $transcript = ${$gene->get_all_Transcripts()}[0];
    my $biotype_group = $gene->get_Biotype->biotype_group();
    my $priority = $biotype_group_priorities->{$biotype_group};

    if(!($selected_gene)) {
      $selected_gene = $gene;
      $selected_transcript = $transcript;
      $max_priority = $priority;
    } elsif($priority > $max_priority) {
      $selected_gene = $gene;
      $selected_transcript = $transcript;
      $max_priority = $priority;
    } elsif($priority == $max_priority and $transcript->length() > $selected_transcript->length()) {
      $selected_gene = $gene;
      $selected_transcript = $transcript;
      $max_priority = $priority;
    }
  }

  if($selected_gene->get_Biotype->biotype_group() eq 'coding') {
    $selected_gene->biotype('protein_coding');
    $selected_transcript->biotype('protein_coding');
  } elsif($selected_gene->get_Biotype->biotype_group() eq 'lnoncoding') {
    $selected_gene->biotype('lncRNA');
    $selected_transcript->biotype('lncRNA');
  } else {
    $selected_transcript->biotype($selected_gene->biotype());
  }

  unless($selected_gene) {
    $self->throw("No gene selected from cluster, something has gone wrong");
  }
  return($selected_gene);
}


#sub layer_cluster_genes {
#  my ($self,$genes) = @_;

  # This is just a simple filtering to make sure we don't make any weird complex clusters like sncRNA transcripts with protein coding ones
  # It does mean that some potentially interesting edge cases might be missed, but this is fine for now
#  my $biotype_group_priorities = { 'coding'     => 5,
#                                   'lnoncoding' => 4,
#                                   'pseudogene' => 3,
#                                   'snoncoding' => 2,
#                                   'mnoncoding' => 1,
#                                   'undefined'  => 0,
#                                 };

#  my $filtered_genes = [];
#  my $max_priority = 0;
#  foreach my $gene (@$genes) {
#    my $biotype_group = $gene->get_Biotype->biotype_group();
#    my $priority = $biotype_group_priorities->{$biotype_group};
#    if($priority > $max_priority) {
#      $max_priority = $priority;
#      $filtered_genes = [];
#      push(@$filtered_genes,$gene);
#    } elsif($priority == $max_priority) {
#      push(@$filtered_genes,$gene);
#    }
#  }

#  return($filtered_genes)
#}


sub generate_biotypes_hash {
  my ($self,$genes) = @_;

  my $unique_biotypes;
  my $biotypes_array = [];
  my $biotypes_hash = {};

  foreach my $gene (@$genes) {
    unless($unique_biotypes->{$gene->biotype()}) {
      push(@$biotypes_array,$gene->biotype());
      $unique_biotypes->{$gene->biotype()} = 1;
    }
  }

  $biotypes_hash->{'slice_genes'} = $biotypes_array;
  return($biotypes_hash);
}


1;
