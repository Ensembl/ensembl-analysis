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
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(compute_translation);
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

  my $program = $self->param('minimap2_path');
  my $paftools =  $self->param('paftools_path');

  unless($program) {
    $program = "minimap2";
  }

  unless($paftools) {
    $paftools = "paftools.js";
  }

  my $source_transcript_hash = {};

  $self->param('source_transcripts',$source_transcript_hash);

  # Potentially should add in something to calculate the max intron size, but since these are already somewhat predictive
  # it also might just be better to leave it at default
  my $input_file = $self->create_input_file($input_genes,$source_transcript_hash);
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Minimap2->new(
       -analysis          => $self->analysis(),
       -program           => $program,
       -paftools_path     => $paftools,
       -genome_index      => $genome_index,
       -input_file        => $input_file,
       -database_adaptor  => $target_dba,
       -bestn             => 10,
       -skip_compute_translation => 1,
       -skip_introns_check       => 1,
       -delete_input_file => 1, # NB!! only set this when creating ranged files, not when using the original input file
  );

  $self->runnable($runnable);

}


sub run {
  my ($self) = @_;

  foreach my $runnable (@{$self->runnable()}) {
    $runnable->run();
    my $mapped_genes = $runnable->output();
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
    say "Final gene: ".$output_gene->stable_id();
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

    my $intial_slice_genes = $genes_by_slice->{$slice->seq_region_name()};
    my $slice_genes = [];

    foreach my $slice_gene (@$intial_slice_genes) {
      unless($self->features_overlap($new_gene,$slice_gene)) {
        next;
      }

      # This is mostly just there if testing to ignore test transcripts
      unless($slice_gene->description =~ /potential_paralogue/) {
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

    my $source_transcript_seq;
    my $transcript_seq;
    if(defined($source_transcript->translation())) {
      $source_transcript_seq = $source_transcript->translateable_seq();
      $transcript_seq = $transcript->translateable_seq();

      unless($transcript_seq) {
        compute_translation($transcript);
        $transcript_seq = $transcript->translateable_seq();
      }

      unless($transcript_seq) {
        next;
      }
      
      if ($transcript->coding_region_start() and
          $transcript->coding_region_end() and
          $transcript->coding_region_end()-$transcript->coding_region_start()+1 < 3) {
        # By convention, the coding_region_end is always higher than the
        # value returned by the coding_region_start method.
        say "Transcript CDS is too short (< 3 bp). Transcript dbID: ".$transcript->dbID();
        next;
      }

      if ($transcript->biotype() eq 'protein_coding' and !($transcript->translation())) {
        say "Transcript biotype is protein_coding but it does not have any translation. Transcript dbID: ".$transcript->dbID();
        next;
      }

    } else {
      $source_transcript_seq = $source_transcript->seq->seq();
      $transcript_seq = $transcript->seq->seq();
      if($transcript->translation()) {
        $self->throw("Found a translation on mapped transcript even though parent has none. Mapped transcript dbID: ".$transcript->dbID());
      }
    }

    my ($coverage,$percent_id,$aligned_source_seq,$aligned_target_seq) = align_nucleotide_seqs($source_transcript_seq,$transcript_seq);
    $transcript->{'cov'} = $coverage;
    $transcript->{'perc_id'} = $percent_id;
    $transcript->{'source_stable_id'} = $source_transcript->stable_id();
    $transcript->{'parent_gene_stable_id'} = $source_transcript->{'parent_gene_stable_id'};
    $transcript->{'source_biotype_group'} = $source_transcript->get_Biotype->biotype_group();
    $transcript->{'source_length'} = $source_transcript->length();
    
    my ($parent_t_stable_id_att) = @{$source_transcript->get_all_Attributes('proj_parent_t')};
    if (!$parent_t_stable_id_att) {
      $self->throw("Issue getting the proj_parent_t attribute for transcript with dbID ".$transcript->dbID());
    }
    my $parent_t_stable_id = $parent_t_stable_id_att->value();
    if (!$parent_t_stable_id) {
      $self->throw("Issue getting the parent transcript stable id from transcript attribute for transcript with dbID ".$transcript->dbID());
    }
    
    my $biotype_group = $transcript->{'source_biotype_group'};
    my $coverage_cutoff = $coverage_cutoff_groups->{$biotype_group};
    my $perc_id_cutoff = $percent_identity_groups->{$biotype_group};
    my $transcript_description = ";parent_transcript=".$parent_t_stable_id.";mapping_coverage=".$coverage.";mapping_identity=".$percent_id;
    if ($source_transcript->translation()) {
      my ($cds_coverage,$cds_percent_id,$aligned_source_seq,$aligned_target_seq) = align_nucleotide_seqs($source_transcript->translateable_seq(),$transcript->translateable_seq());
      my $cds_description = ";cds_coverage=".$cds_coverage.";cds_identity=".$cds_percent_id;
      my $aligned_source_seq_copy = $aligned_source_seq;
      my $aligned_target_seq_copy = $aligned_target_seq;
      $aligned_source_seq_copy =~ s/\-\-\-//g;
      $aligned_target_seq_copy =~ s/\-\-\-//g;

      if ($aligned_source_seq_copy =~ /\-/ or $aligned_target_seq_copy =~ /\-/) {
        $cds_description .= ";cds_gap=1";
#        my $transcript_attrib = Bio::EnsEMBL::Attribute->new(-CODE => 'proj_parent_t',
#                                                             -VALUE => ">source_cds_align\n".$aligned_source_seq."\n>target_cds_align\n".$aligned_target_seq."\n");
#        $transcript->add_Attributes($transcript_attrib);
#        my $translation_attrib = Bio::EnsEMBL::Attribute->new(-CODE => 'proj_parent_t',
#                                                              -VALUE => ">source_translation\n".$source_transcript->translation->seq().
#                                                                        "\n>target_translation\n".$transcript->translation->seq()."\n");
#        $transcript->translation->add_Attributes($translation_attrib);
      } else {
        $cds_description .= ";cds_gap=0";
      }
      $transcript->{'cds_description'} = $cds_description;
    }

    if ($transcript->{'cds_description'}) {
      $transcript_description .= $transcript->{'cds_description'};
    }

    $transcript_description .= ";annotation_method=minimap_paralogue";
    $transcript->description($transcript_description);

    # add source transcript stable id as transcript attribute
    my $parent_attribute = Bio::EnsEMBL::Attribute->new(-CODE => 'proj_parent_t',-VALUE => $source_transcript->stable_id().".".$source_transcript->version());
    $transcript->add_Attributes($parent_attribute);

    #my $gene_description = ";parent_gene=".$source_transcript->{'parent_gene_stable_id'}.";mapping_type=potential_paralogue";
    #$gene->description($gene_description);
    if ($gene->description() !~ /;parent_gene=([^;]+);/) {
      $gene->description($source_gene->description().";parent_gene=".$source_transcript->{'parent_gene_stable_id'});
    }
    $gene->description($gene->description().";mapping_type=potential_paralogue");

    # add source gene stable id as gene attribute
    $parent_attribute = Bio::EnsEMBL::Attribute->new(-CODE => 'proj_parent_g',-VALUE => $source_transcript->{'parent_gene_stable_id'});
    $gene->add_Attributes($parent_attribute);

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
    my $gene_description = $gene->description();
    if ($gene_description) {
      my ($type) = $gene_description =~ /;mapping_type=(.+)$/;
      if ($type and $type eq 'potential_paralogue') {
        # Just in case there's accidental re-runs or testing
        next;
      }
    }
    my ($parent_stable_id_att) = @{$gene->get_all_Attributes('proj_parent_g')};
    if (!$parent_stable_id_att) {
      $self->throw("Issue getting the proj_parent_g attribute for gene with dbID ".$gene->dbID().". Description: ".$gene_description);
    }
    my $parent_stable_id = $parent_stable_id_att->value();
    if (!$parent_stable_id) {
      $self->throw("Issue getting the parent stable id from gene attribute for gene with dbID ".$gene->dbID().". Description: ".$gene_description);
    }

    # Just going to use the canonical
    my $transcript = $gene->canonical_transcript();
    $transcript->{'parent_gene_stable_id'} = $parent_stable_id;
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
    if($transcript->translation()) {
      if($longest_translation) {
        if(length($transcript->translateable_seq) > length($longest_translation->translateable_seq())) {
          $longest_translation = $transcript;
        } elsif((length($transcript->translateable_seq()) > length($longest_translation->translateable_seq())) and ($transcript->length() > $longest_translation->length())) {
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

  my $coverage_cutoff = 98;
  my $perc_id_cutoff = 95;

  my $input_genes = [];
  my $source_gene_adaptor = $source_gene_dba->get_GeneAdaptor();
  foreach my $gene_id (@$gene_ids) {
    my $gene = $source_gene_adaptor->fetch_by_dbID($gene_id);
    unless($gene) {
      next;
      $self->throw("Couldn't fetch gene from source db with dbID: ".$gene_id);
    }

   $self->set_canonical($gene);
   my $transcript = $gene->canonical_transcript();
    if ($transcript) {
      my $transcript_description = $transcript->description();
      $transcript_description =~ /;mapping_coverage=([0-9.]+);.*mapping_identity=([0-9.]+)/;
      my $coverage = $1;
      my $perc_id = $2;
      unless(defined($coverage) and defined($perc_id)) {
        $self->throw("Issue parsing coverage and percent id from transcript description. Transcript description:\n".$transcript_description);
      }

      # Skip over canonical transcripts that don't have near perfect mappings
      unless($coverage >= $coverage_cutoff and $perc_id >= $perc_id_cutoff) {
        next;
      }
    } else {
      next;
    }

    push(@$input_genes,$gene);
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
