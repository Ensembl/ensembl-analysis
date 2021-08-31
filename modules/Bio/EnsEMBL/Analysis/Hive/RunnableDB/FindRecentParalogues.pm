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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::FindRecentParalogues;

use warnings;
use strict;
use feature 'say';
use Data::Dumper;
use File::Spec::Functions qw(tmpdir catfile);

use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name align_nucleotide_seqs);
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);
use Bio::EnsEMBL::Analysis::Runnable::Minimap2Remap;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);
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

  my $genome_index = $self->param_required('genome_index');
  my $genome_file = $self->param_required('genome_file');
  unless(-e $genome_file) {
    $self->throw("Could not find the genome file. Path used:\n".$genome_file);
  }

  my $input_genes = $self->fetch_input_genes_by_id($self->param('iid'),$target_dba);

  say "Processing ".scalar(@$input_genes)." genes for mapping";
  $self->param('input_genes',$input_genes);

#  my $parent_genes_hash = {};
#  foreach my $gene (@$input_genes) {
#    $parent_genes_hash->{$gene->dbID()} = $gene;
#  }
#  $self->param('parent_genes_hash',$parent_genes_hash);

  my $program = $self->param('minimap2_path');
  my $paftools =  $self->param('paftools_path');

  unless($program) {
    $program = "minimap2";
  }

  unless($paftools) {
    $paftools = "paftools.js";
  }

  my $source_transcript_hash = {};
  foreach my $gene (@$input_genes) {
    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts) {
      $source_transcript_hash->{$transcript->dbID()} = $transcript;
    }
  }

  $self->param('source_transcripts',$source_transcript_hash);

  my $input_file = $self->create_input_file($input_genes);
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Minimap2->new(
       -analysis          => $self->analysis(),
       -program           => $program,
       -paftools_path     => $paftools,
       -genome_index      => $genome_index,
       -input_file        => $input_file,
       -database_adaptor  => $target_dba,
       -bestn             => 10,
       -delete_input_file => 1, # NB!! only set this when creating ranged files, not when using the original input file
  );

  $self->runnable($runnable);

}


sub run {
  my ($self) = @_;

  foreach my $runnable (@{$self->runnable()}) {
    $runnable->run();
    my $mapped_genes = $runnable->output();
    my $parent_genes_hash = $self->param('parent_genes_hash');
    my $filtered_new_genes = $self->filter_new_genes($mapped_genes);
    $self->output($filtered_new_genes);
  }
}



sub write_output {
  my ($self) = @_;
  my $output_dba = $self->hrdb_get_con('target_db');
  my $output_gene_adaptor = $output_dba->get_GeneAdaptor;
  my $output_genes = $self->output();
  foreach my $output_gene (@$output_genes) {
    say "Final gene; ".$output_gene->stable_id();
    $output_gene->biotype($output_gene->biotype());
    empty_Gene($output_gene);
    $output_gene_adaptor->store($output_gene);
  }

  return 1;
}


sub filter_new_genes {
  my ($self,$genes) = @_;

  my $source_transcripts_hash = $self->param('source_transcripts');

  my $target_db = $self->hrdb_get_con('target_db');
  my $genes_by_slice = {};

  say "Number of initial genes from mapping: ".scalar(@$genes);

  my $secondary_filtered_genes = $self->select_secondary_mappings($genes,$source_transcripts_hash);
  say "Number of genes after removing the primary mapping: ".scalar(@$secondary_filtered_genes);

  my $genes_after_cutoffs = $self->filter_by_cutoffs($secondary_filtered_genes,$source_transcripts_hash);
  say "Number of genes after initial cut-offs: ".scalar(@$genes_after_cutoffs);

  my $final_filtered_genes = $self->cluster_source_genes($genes_after_cutoffs);
  say "Number of genes after clustering against main gene set: ".scalar(@$final_filtered_genes);

  return($final_filtered_genes);
}


sub select_secondary_mappings {
  my ($self,$genes,$source_transcripts_hash) = @_;

  my $secondary_genes = [];
  my $mapped_transcript_ids_hash = {};
  # This is a bit of a mess because of things being genes at this point, but basically we want to go through everything
  # at this point and only select true secondary mappings. So to do this there are a few possibilities, but the simplest
  # one is to just look to see if the gene overlaps the source transcript. That way we should get rid of all the primary
  # mappings, and also the vast amount of the genes to process
  foreach my $gene (@$genes) {
    my $transcript = ${$gene->get_all_Transcripts()}[0];
    my $source_transcript = $source_transcripts_hash->{$transcript->stable_id()};
    unless($source_transcript) {
      $self->throw("Issue finding source transcript for transcript with dbID: ".$transcript->stable_id());
    }

    unless($self->features_overlap($transcript,$source_transcript)) {
      push(@$secondary_genes,$gene);
    }
  }

  return($secondary_genes);
}


sub features_overlap {
  my ($self,$feature_a, $feature_b) = @_;

  if(($feature_a->seq_region_start() <= $feature_b->seq_region_end() ) && ($feature_a->seq_region_end() >= $feature_b->seq_region_start())) {
    return 1;
  }

  return 0;
}


sub cluster_source_genes {
  my ($self,$new_genes) = @_;

  my $genes_by_slice = {};
  my $filtered_genes = [];

  foreach my $new_gene (@$new_genes) {
    my $slice = $new_gene->slice();

    unless($genes_by_slice->{$slice->seq_region_name()}) {
      $genes_by_slice->{$slice->seq_region_name()} = $slice->get_all_Genes();
    }

    $new_gene->biotype("new_".$new_gene->biotype);

    my $intial_slice_genes = $genes_by_slice->{$slice->seq_region_name()};
    my $slice_genes = [];

    foreach my $slice_gene (@$intial_slice_genes) {
      unless($self->features_overlap($new_gene,$slice_gene)) {
        next;
      }

      # This is mostly just there if testing to ignore test transcripts
      unless($slice_gene->biotype() =~ /^new/) {
        push(@$slice_genes,$slice_gene);
      }
    }

    my $biotypes_hash = $self->generate_biotypes_hash($slice_genes);
    $biotypes_hash->{'new_gene'} = [$new_gene->biotype()];

    my $all_genes = [$new_gene,@$slice_genes];
    say "Clustering new gene with ".scalar(@$slice_genes)." genes on parent slice...";
    my ($clusters, $unclustered) = cluster_Genes($all_genes,$biotypes_hash);
    say "...finished clustering genes";

    say "Found ".scalar(@$clusters)." clusters";
    say "Found ".scalar(@$unclustered)." unclustered genes";

    # For the moment we'll just keep things that have no exon overlap with an existing gene
    # We probably want to put something smarter in, for example it could be that a new protein coding gene overlaps with a pseudogene
    # In that case we might want to keep the new gene under the assumption that it works in the target genome
    # Similarly, may want to process clusters to find clusters where the overlap is between a new protein coding or lncRNA gene and an existing
    # small non-coding one. In this case, we'd probably want to keep the new gene even if there's an overlap
    foreach my $unclustered (@$unclustered) {
      my $unclustered_genes = $unclustered->get_Genes_by_Set('new_gene');
      foreach my $unclustered_gene (@$unclustered_genes) {
        push(@$filtered_genes,$unclustered_gene);
      }
    }
  }

  return($filtered_genes);
}


sub filter_by_cutoffs {
  my ($self,$genes,$source_transcripts_hash) = @_;

  my $filtered_genes = [];
  my $coverage_cutoff_groups = {'coding' => 95,
                                'pseudogene' => 80,
                                'snoncoding' => 95,
                                'mnoncoding' => 95,
                                'lnoncoding' => 90,
                                'undefined'  => 95};

  my $percent_identity_groups = {'coding' => 95,
                                 'pseudogene' => 90,
                                 'snoncoding' => 95,
                                 'mnoncoding' => 95,
                                 'lnoncoding' => 90,
                                 'undefined'  => 95};

  foreach my $gene (@$genes) {
    my $transcript = ${$gene->get_all_Transcripts()}[0];
    my $source_transcript = $source_transcripts_hash->{$transcript->stable_id()};
    unless($source_transcript) {
      $self->throw("Couldn't fetch the source transcript for mapped transcript with dbID: ".$transcript->stable_id());
    }

    my $source_gene = $source_transcript->get_Gene();
    $gene->biotype($source_gene->biotype());
    $gene->version($source_gene->version());
    $gene->stable_id($source_gene->stable_id());
    $transcript->biotype($source_transcript->biotype());
    my $gene_description = "Parent: ".$source_gene->stable_id().".".$source_gene->version().", Type: Potential paralogue";
    $gene->description($gene_description);

    my $source_transcript_seq;
    my $transcript_seq;
    if($source_transcript->translation()) {
      $source_transcript_seq = $source_transcript->translateable_seq();
      $transcript_seq = $transcript->translateable_seq();

      unless($transcript_seq) {
        compute_translation($transcript);
        $transcript_seq = $transcript->translateable_seq();
      }
    } else {
      $source_transcript_seq = $source_transcript->seq->seq();
      $transcript_seq = $transcript->seq->seq();
    }

    my ($coverage,$percent_id,$aligned_source_seq,$aligned_target_seq) = align_nucleotide_seqs($source_transcript_seq,$transcript_seq);
    $transcript->{'cov'} = $coverage;
    $transcript->{'perc_id'} = $percent_id;
    $transcript->{'source_stable_id'} = $source_transcript->stable_id();
    $transcript->{'parent_gene_stable_id'} = $source_transcript->get_Gene->stable_id().".".$source_transcript->get_Gene->version();
    $transcript->{'source_biotype_group'} = $source_transcript->get_Biotype->biotype_group();
    $transcript->{'source_length'} = $source_transcript->length();
    my $biotype_group = $transcript->{'source_biotype_group'};
    my $coverage_cutoff = $coverage_cutoff_groups->{$biotype_group};
    my $perc_id_cutoff = $percent_identity_groups->{$biotype_group};
    my $transcript_description = "Parent: ".$source_transcript->stable_id().".".$source_transcript->version().", Coverage: ".$coverage.", Perc id: ".$percent_id;
    $transcript->description($transcript_description);

    if($transcript->{'cov'} >= $coverage_cutoff and $transcript->{'perc_id'} >= $perc_id_cutoff) {
      say "Transcript ".$transcript->{'source_stable_id'}." (".$transcript->stable_id().", ".$biotype_group.") passed check: ".$transcript->{'cov'}." cov, ".
          $transcript->{'perc_id'}." perc_id, ".$transcript->{'cds_length_diff'}." length diff";
      push(@$filtered_genes,$gene);
    }
  }

  return($filtered_genes);
}


sub runnable_failed {
  my ($self,$runnable_failed) = @_;
  unless ($self->param_is_defined('_runnable_failed')) {
    $self->param('_runnable_failed',[]);
  }
  if ($runnable_failed) {
    push (@{$self->param('_runnable_failed')},$runnable_failed);
  }
  return ($self->param('_runnable_failed'));
}


sub write_input_file {
  my ($self,$genomic_reads) = @_;

  my $output_file = $self->create_filename();
  open(OUT,">".$output_file);
  foreach my $genomic_read (@$genomic_reads) {
    say OUT $genomic_read;
  }
  close OUT;

  return($output_file);
}


sub create_filename{
  my ($self, $stem, $ext, $dir, $no_clean) = @_;
  return create_file_name($stem, $ext, $dir, $no_clean);
}


sub create_input_file {
  my ($self,$genes,$source_transcript_hash) = @_;

  # In this we'll just take the canonical of the input genes and write then to file
  # Output the db id of the transcript in the header to help figure out later what the overlap is
  my $output_file = $self->create_filename();

  open(OUT,">".$output_file);
  foreach my $gene (@$genes) {
    my $transcript = $gene->canonical_transcript();
    unless($transcript) {
      $transcript = $self->set_canonical($gene);
    }
    $source_transcript_hash->{$transcript->dbID()} = $transcript;
    my $seq = $transcript->seq->seq();
    say OUT ">".$transcript->dbID();
    say OUT $seq;
  }
  close OUT;

  return($output_file);
}


sub set_canonical {
  my ($self,$gene) = @_;

  my $longest_translation;
  my $longest_transcript;

  my $transcripts = $gene->get_all_Transcripts();
  foreach my $transcript (@$transcripts) {
    if($transcript->translation) {
      if($longest_translation) {
        if($transcript->translation->length() > $longest_translation->translation->length()) {
          $longest_translation = $transcript;
        } elsif(($transcript->translation->length() > $longest_translation->translation->length()) and ($transcript->length() > $longest_translation->length())) {
         $longest_translation = $transcript;
        }
      } else {
        $longest_translation = $transcript;
      }
    } else {
      if($longest_transcript) {
        if($transcript->length() > $longest_transcript->length()) {
          $longest_transcript = $transcript;
        }
      } else{
        $longest_transcript = $transcript;
      }
    } # End else
  } # End foreach

  if($longest_translation) {
    $gene->canonical_transcript($longest_translation);
  } else {
    $gene->canonical_transcript($longest_transcript);
  }
}

sub fetch_input_genes_by_id {
  my ($self,$gene_ids,$source_gene_dba) = @_;

  # TEST!!!!
#  $gene_ids = [34154,7703];
#  $gene_ids = [35728];
  # TEST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  my $input_genes = [];
  my $source_gene_adaptor = $source_gene_dba->get_GeneAdaptor();
  foreach my $gene_id (@$gene_ids) {
    my $gene = $source_gene_adaptor->fetch_by_dbID($gene_id);

    unless($gene) {
      $self->throw("Couldn't fetch gene from source db with dbID: ".$gene_id);
    }

    if($gene->seq_region_name() eq 'MT') {
      next;
    }

    my $is_readthrough = 0;
    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts) {
      if($is_readthrough) {
         last;
      }

      my $attributes = $transcript->get_all_Attributes();
      foreach my $attribute (@{$attributes}) {
        if($attribute->value eq 'readthrough') {
          $is_readthrough = 1;
          last;
        }
      } # foreach my $attribute
    } # End foreach my $transcript

    unless($is_readthrough) {
      push(@$input_genes,$gene);
    }
  }

  return($input_genes);
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


1;
