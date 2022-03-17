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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::FinaliseRemappedGeneset

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::FinaliseRemappedGeneset;

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

  # Define the source and target gene dbs
  my $source_gene_dba = $self->hrdb_get_dba($self->param('source_gene_db'));
  $self->hrdb_set_con($source_gene_dba,'source_gene_db');


  # Get all the genes and do some potential finalisations
  # 1) Remove genes with transcripts with low percent ids that overlap another gene. These are likely just paralogues misaligned (like Y genes on X)
  # 2) Examine genes that have duplicate ids, where both are the primary mapping. See if there is a way to remove one if it looks bad/holds little information
  # 3) Need to check biotypes of the genes, particularly in the case of the duplicates. It may be that the dup has an incorrect biotype since it's not connected to the main locus
  # 4) Examine CDS sequences for stop codons, perhaps employing the stop codon removal code for now
  # 5) Add gene symbols based on the symbol of the parent if present
  # 6) Add parent gene/transcript attribs like in the mouse strains
  # 7) Could consider something about pseudogenes/polymorphic pseudogenes now being coding, but maybe out of scope


  my $gene_adaptor = $target_dba->get_GeneAdaptor();
  my $all_initial_genes = $gene_adaptor->fetch_all();

  # First we'll parse the descriptions of all the gene/transcripts and filter out the potential paralogues, as these have already been qc'd at this point
  my $primary_genes = $self->parse_descriptions($all_initial_genes);

  say "Have a total of ".scalar(@$primary_genes)." initial genes to process";

  # We'll need the genes organised by slice for later
  my $genes_by_slice_hash = $self->sort_genes_by_slice($primary_genes);

  $self->param('primary_genes',$primary_genes);
  $self->param('genes_by_slice',$genes_by_slice_hash);
}


sub run {
  my ($self) = @_;

  my $initial_genes = $self->param('primary_genes');
  my $first_filtering = $self->filter_overlapping($self->param('genes_by_slice'));

  # Next we want to build a list of the duplicate ids to process them separately
  my $duplicate_genes_hash = $self->generate_duplicate_genes_hash($first_filtering);
  my $duplicate_filtering = $self->filter_duplicates($first_filtering,$duplicate_genes_hash);

  my $genes_to_remove = $self->genes_to_remove($initial_genes,$duplicate_filtering);
  $self->output($genes_to_remove);
}


sub write_output {
  my ($self) = @_;

  my $output_dba = $self->hrdb_get_con('target_db');
  my $output_gene_adaptor = $output_dba->get_GeneAdaptor();
  my $output_slice_adaptor = $output_dba->get_SliceAdaptor();
  my $genes_to_remove = $self->output();

  # Remove the initial set of paralogues prior to collapsing
  foreach my $gene (@$genes_to_remove) {
    $output_gene_adaptor->remove($gene);
  }

  $self->add_gene_symbols($output_gene_adaptor,$output_slice_adaptor);
}


sub parse_descriptions {
  my ($self,$genes) = @_;

  # This will parse the gene/transcript descriptions and store the info separately on the gene/transcripts under temp keys
  # This will make deciding what to keep/remove more straightforward later
  # It will also calculate an average cov/perc id for the gene from all the transcripts
  # The returned array will only have primary transcripts as the the potential paralogues have essentially been qc'd by now
  my $avg_cov_threshold = 80;
  my $avg_perc_id_threshold = 95;
  my $primary_genes = [];
  foreach my $gene (@{$genes}) {
    my $gene_description = $gene->description();
    $gene_description =~ /;parent_gene=([^\,]+);mapping_type=(.+)$/;
    my $gene_stable_id = $1;
    my $gene_type = $2;
    $gene->{'parent_stable_id'} = $gene_stable_id; # Not really needed, since the stable id should be assigned, but doing it as it may be useful
    $gene->{'type'} = $gene_type;

    my $total_cov = 0;
    my $total_perc_id = 0;
    my $transcripts = $gene->get_all_Transcripts();
    if (scalar(@$transcripts) == 0) {
      next;
    }
    foreach my $transcript (@$transcripts) {
      # Parse the cov and percent id out, also track the total cov/perc id to get average for gene
      my $transcript_description = $transcript->description();

      $transcript_description =~ /;parent_transcript=([^\,]+);mapping_coverage=([^\,]+);mapping_identity=(.+)$/;
      my $transcript_stable_id = $1;
      my $cov = $2;
      my $perc_id = $3;
      $transcript->{'parent_stable_id'} = $transcript_stable_id;
      $transcript->{'cov'} = $cov;
      $transcript->{'perc_id'} = $perc_id;
      $total_cov += $cov;
      $total_perc_id += $perc_id;
    }

    my $avg_cov = $total_cov/scalar(@$transcripts);
    my $avg_perc_id = $total_perc_id/scalar(@$transcripts);
    $gene->{'avg_cov'} = $avg_cov;
    $gene->{'avg_perc_id'} = $avg_perc_id;

    if($avg_cov >= $avg_cov_threshold and $avg_perc_id >= $avg_perc_id_threshold) {
      $gene->{'quality_call'} = 'good';
    } else {
      $gene->{'quality_call'} = 'bad';
    }

    if($gene_type eq 'primary_mapping') {
      push(@$primary_genes,$gene);
    }
  }

  return($primary_genes);
}


sub sort_genes_by_slice {
  my ($self,$genes) = @_;

  my $genes_by_slice_hash = {};
  foreach my $gene (@$genes) {
    unless($genes_by_slice_hash->{$gene->seq_region_name()}) {
      $genes_by_slice_hash->{$gene->seq_region_name()} = [];
    }
    push(@{$genes_by_slice_hash->{$gene->seq_region_name()}},$gene);
  }

  return($genes_by_slice_hash);
}


sub generate_duplicate_genes_hash {
  my ($self,$genes) = @_;

  say "Searching for duplicate parent ids in the set of primary genes";
  my $initial_genes_hash = {};
  my $duplicate_genes_hash = {};

  foreach my $gene (@$genes) {
    my $parent_stable_id = $gene->{'parent_stable_id'};

    unless($initial_genes_hash->{$parent_stable_id}) {
      $initial_genes_hash->{$parent_stable_id} = [];
    }
    push(@{$initial_genes_hash->{$parent_stable_id}},$gene);
  }

  foreach my $id (keys(%{$initial_genes_hash})) {
    if(scalar(@{$initial_genes_hash->{$id}}) > 1) {
      $duplicate_genes_hash->{$id} = $initial_genes_hash->{$id};
    }
  }

  say "Found a total of ".scalar(keys(%{$duplicate_genes_hash}))." sets of duplicates";
  return($duplicate_genes_hash);
}


sub filter_overlapping {
  my ($self,$genes_by_slice_hash) = @_;

  # This is a tricky one, we want to get overlapping genes on the current gene, and if the there are good genes of certain biotype groups, then remove the bad one
  my $filtered_genes = [];
  foreach my $slice_name (keys(%{$genes_by_slice_hash})) {
    my $genes = $genes_by_slice_hash->{$slice_name};
    my $biotypes_hash = $self->generate_biotypes_hash($genes);
    say "Processing ".scalar(@$genes)." input genes for slice: ".$slice_name;

    my ($clusters, $unclustered) = cluster_Genes($genes,$biotypes_hash);

    foreach my $cluster (@$clusters) {
      my $clustered_genes = $cluster->get_Genes();
      my $good_biotype_groups = {};

      # First we'll track the good gene's biotype groups. Basically if there's a protein coding/long non-coding/pseudogene then if there's an overlapping
      # bad gene we should just remove the bad one. Otherwise the good genes are small ones, probably not worth excluding the bad gene for
      foreach my $gene (@$clustered_genes) {
        if($gene->{'quality_call'} eq 'good') {
          $good_biotype_groups->{$gene->get_Biotype->biotype_group()} = 1;
        }
      }

      foreach my $gene (@$clustered_genes) {
        if($gene->{'quality_call'} eq 'good') {
          push(@$filtered_genes,$gene);
        } else {
          # Skip if there's a good one of the following in the cluster
          if($good_biotype_groups->{'coding'} or $good_biotype_groups->{'lnoncoding'} or $good_biotype_groups->{'pseudogene'}) {
            next;
          }
        } # End else
      } # End foreach my $gene (@$clustered_genes)
    }

    # Just push singletons onto the filtered array
    foreach my $singleton (@$unclustered) {
      my $unclustered_genes = $singleton->get_Genes();
      push(@$filtered_genes,@$unclustered_genes);
    }
  } # End foreach my $slice_name

  say "Have a total of ".scalar(@$filtered_genes)." genes after removing weak overlapping genes";
  return($filtered_genes);
}



sub filter_duplicates {
  my ($self,$genes,$duplicate_genes_hash) = @_;

  my $filtered_genes = [];
  my $source_gene_db = $self->hrdb_get_con('source_gene_db');
  my $source_gene_adaptor = $source_gene_db->get_GeneAdaptor();

  # Here we want to figure out if there's a straightforward way to pick between loci, particularly in a split gene
  # There are a few simple options. One would be to select any good loci and remove and bad loci. The next bit is
  # What to do if there's multiple good loci. In that case we could remove any less important loci. So if we had
  # one with protein coding transcripts and the other with non-coding transcripts, then pick the protein coding one

  my $processed_duplicate_hash = {};
  foreach my $gene (@$genes) {
    my $gene_parent_id = $gene->{'parent_stable_id'};
    if($processed_duplicate_hash->{$gene_parent_id}) {
      next;
    }

    if($duplicate_genes_hash->{$gene_parent_id}) {
      $processed_duplicate_hash->{$gene_parent_id} = 1;
      my $gene_set = $duplicate_genes_hash->{$gene_parent_id};

      my $good_genes = [];
      my $bad_genes = [];
      foreach my $gene (@$gene_set) {
        if($gene->{'quality_call'} eq 'good') {
          push(@$good_genes,$gene);
        } else {
          push(@$bad_genes,$gene);
        }
      }

      # For some reason asking to fetch by versioned stable id didn't work using the method. So unversioning it to get the gene
      my $canonical_gene;
      my $unversioned_id = $gene_parent_id;
      $unversioned_id =~ s/\.\d+$//;
      my $parent_gene = $source_gene_adaptor->fetch_by_stable_id($unversioned_id);
      unless($parent_gene) {
        $self->throw("Couldn't fetch the parent source gene from the source db using versioned stable id. Stable id used: ".$unversioned_id);
      }

      my $canonical_transcript_id = $parent_gene->canonical_transcript->stable_id_version();
      foreach my $duplicate_gene (@$gene_set) {
        my $transcripts = $duplicate_gene->get_all_Transcripts();
        foreach my $transcript (@$transcripts) {
          if($transcript->{'parent_stable_id'} eq $canonical_transcript_id) {
            $canonical_gene = $duplicate_gene;
          }

          if($transcript->translation()) {
            $duplicate_gene->{'coding'} = 1;
          }
        } # foreach my $transcript
      } # foreach my $duplicate_gene

      # This actually doesn't work, as sometimes the canonical is in two genes
      # If the canonical is good, then just go with that gene
 #     if($canonical_gene and $canonical_gene->{'quality_call'} eq 'good') {
 #       push(@$filtered_genes,$gene);
 #     } else {

      my $good_coding = 0;
      foreach my $gene (@$good_genes) {
        if($gene->{'coding'}) {
           push(@$filtered_genes,$gene);
           $good_coding = 1;
        }
      } # End foreach my $gene (@$good_genes)

      # At this point we have no coding good genes, so take any good non-coding and bad coding
      # Then failing that, there's basically nothing we can say, so just keep them all
      unless($good_coding) {
        my $final_selections = [];
        foreach my $gene (@$good_genes) {
          push(@$final_selections,$gene);
        }

        foreach my $gene (@$bad_genes) {
          if($gene->{'coding'}) {
            push(@$final_selections,$gene);
          }
        } # foreach my $gene (@$bad_genes)

        if(scalar(@$final_selections)) {
          push(@$filtered_genes,@$final_selections);
        } else {
          # Just give up at this point
          push(@$filtered_genes,@$gene_set);
        }
      } # End unless($good_coding)
    } else {
      push(@$filtered_genes,$gene);
    }
  }

  say "Have a total of ".scalar(@$filtered_genes)." genes after processing duplicates";
  return($filtered_genes);
}


sub genes_to_remove {
  my ($self,$initial_genes,$filtered_genes) = @_;

  my $genes_to_remove = [];
  # This is just to figure out what has been removed along the way and flag those up for deletion
  my $filtered_gene_ids = {};
  foreach my $gene (@$filtered_genes) {
    $filtered_gene_ids->{$gene->dbID()} = 1;
  }

  foreach my $gene (@$initial_genes) {
    unless($filtered_gene_ids->{$gene->dbID()}) {
      push(@$genes_to_remove,$gene);
    }
  }

  say "Have a total of ".scalar(@$genes_to_remove)." genes to delete";
  return($genes_to_remove);
}


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


sub add_gene_symbols {
  my ($self,$target_gene_adaptor,$target_slice_adaptor) = @_;

  say "Setting gene symbols using source genes in the target gene set";
  my $source_gene_db = $self->hrdb_get_con('source_gene_db');
  my $target_gene_db = $self->hrdb_get_con('target_db');
  my $source_gene_adaptor = $source_gene_db->get_GeneAdaptor();
  my $target_slices = $target_slice_adaptor->fetch_all('toplevel');
  foreach my $slice (@$target_slices) {
    my $genes = $slice->get_all_Genes();
    foreach my $gene (@$genes) {
      my $gene_description = $gene->description();
      $gene_description =~ /;parent_gene=([^\,]+);mapping_type=(.+)$/;
      my $gene_stable_id = $1;
      my $gene_type = $2;
      my $unversioned_parent_stable_id = $gene_stable_id;
      $unversioned_parent_stable_id =~ s/\.(\d+)$//;
      my $source_gene = $source_gene_adaptor->fetch_by_stable_id($unversioned_parent_stable_id);
      unless($source_gene) {
        say "Couldn't find parent gene with stable id: ".$unversioned_parent_stable_id;
        next;
      }

      my $xref = $source_gene->display_xref();
      if($xref) {
        $gene->description($gene->description().";parent_gene_display_id=".$source_gene->display_id());
        $target_gene_adaptor->update($gene);
        my $dbea = $target_gene_db->get_DBEntryAdaptor();
        $dbea->store($xref,$gene->dbID(),'Gene',1); # 1 to ignore the external db version
      }
    } # End foreach my $gene
  } # End foreach my $slice
}


1;
