=head1 LICENSE

Copyright [2020] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::Minimap2Remap;

use warnings;
use strict;
use feature 'say';
use Data::Dumper;
use File::Spec::Functions qw(tmpdir catfile);

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name);
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);
use Bio::EnsEMBL::Analysis::Runnable::Minimap2Remap;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw (align_proteins);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my($self) = @_;
  $self->create_analysis;

  $self->param('region_padding',500);

  # Define the dna dbs
  my $source_dna_dba = $self->hrdb_get_dba($self->param('source_dna_db'));
  my $target_dna_dba = $self->hrdb_get_dba($self->param('target_dna_db'));
  $self->hrdb_set_con($source_dna_dba,'source_dna_db');
  $self->hrdb_set_con($target_dna_dba,'target_dna_db');

  # Define the source and target gene dbs
  my $source_gene_dba = $self->hrdb_get_dba($self->param('source_gene_db'));
  my $target_gene_dba = $self->hrdb_get_dba($self->param('target_gene_db'));
  $source_gene_dba->dnadb($source_dna_dba);
  $target_gene_dba->dnadb($target_dna_dba);
  $self->hrdb_set_con($source_gene_dba,'source_gene_db');
  $self->hrdb_set_con($target_gene_dba,'target_gene_db');

  my $input_ids = [];
  my $slice_adaptor = $source_gene_dba->get_SliceAdaptor();
  my $test_slice = $slice_adaptor->fetch_by_region('toplevel','1');
  my $genes = $test_slice->get_all_Genes();
  foreach my $gene (@$genes) {
 #   unless($gene->biotype eq 'protein_coding') {
 #     next;
 #   }
    push(@$input_ids,$gene->dbID);
  }


#  $input_ids = [1330904,1331879,1332102,1333521,1335796];

  my $sequence_adaptor = $source_gene_dba->get_SequenceAdaptor();

  my $parent_gene_id_hash = {};
  my $genomic_reads = [];
  say "Processing ".scalar(@$input_ids)." genes";
  foreach my $input_id (@$input_ids) {
    my $gene = $source_gene_dba->get_GeneAdaptor->fetch_by_dbID($input_id);
    my $gene_id = $gene->dbID();
    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts) {
      my $transcript_id = $transcript->dbID();
      my $biotype = $transcript->get_Biotype();
      my $biotype_group = $biotype->biotype_group();
      $parent_gene_id_hash->{$transcript_id}->{'gene_id'} = $gene_id;
      $parent_gene_id_hash->{$transcript_id}->{'biotype_group'} = $biotype_group;
    }

    my $slice = $gene->slice();
    my $stable_id = $gene->stable_id.".".$gene->version;


    my $region_start = $gene->seq_region_start - $self->param('region_padding');
    if($region_start <= 0) {
      $region_start = 1;
    }

    my $region_end = $gene->seq_region_end + $self->param('region_padding');
    if($region_end > $slice->length()) {
      $region_end = $slice->length();
    }

    my $strand = $gene->strand;

    my $genomic_seq = ${ $sequence_adaptor->fetch_by_Slice_start_end_strand($slice, $region_start, $region_end, $strand) };

    my $fasta_record = ">".$gene_id."\n".$genomic_seq;
    push(@$genomic_reads, $fasta_record);
  } #close foreach input_id

  my $input_file = $self->write_input_file($genomic_reads);
  $self->param('input_file',$input_file);

  my $analysis = $self->create_analysis;
  $analysis->logic_name("minimap2remap");
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Minimap2Remap->new(
       -analysis          => $analysis,
       -program           => $self->param('minimap2_path'),
       -paftools_path     => $self->param('paftools_path'),
       -genome_index      => $self->param_required('genome_index'),
       -input_file        => $input_file,
       -source_adaptor  => $source_gene_dba,
       -target_adaptor  => $target_gene_dba,
       -parent_gene_ids => $parent_gene_id_hash,
  );
  $self->runnable($runnable);
}


sub write_output {
  my ($self) = @_;
  my $output_dba = $self->hrdb_get_con('target_gene_db');
  my $output_gene_adaptor = $output_dba->get_GeneAdaptor;
  my $output_genes = $self->output();
  foreach my $output_gene (@$output_genes) {
    $output_gene_adaptor->store($output_gene);
  }

  return 1;
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


sub process_transcript {
  my ($self,$transcript,$compara_dba,$mlss,$source_genome_db,$source_transcript_dba) = @_;

  my $max_cluster_gap_length = $self->max_cluster_gap_length($transcript);
  say "Max align gap: ".$max_cluster_gap_length;
  my $all_target_genomic_aligns = [];
  my $exons = $transcript->get_all_Exons;
  my $exon_region_padding = $self->param('exon_region_padding');
  foreach my $exon (@{$exons}) {
##    my $exon_padded_start = $exon->seq_region_start - $exon_region_padding;
##   if($exon_padded_start < 0) {
##      $exon_padded_start = 0;
##    }

##    my $exon_padded_end = $exon->seq_region_end + $exon_region_padding;
##    if($exon_padded_end > $transcript->slice->length) {
##      $exon_padded_end = $transcript->slice->length;
##    }

##    my $slice_adaptor = $source_transcript_dba->get_SliceAdaptor();
##    my $exon_slice = $slice_adaptor->fetch_by_region($exon->slice->coord_system_name, $exon->slice->seq_region_name, $exon_padded_start, $exon_padded_end);


    my $exon_region_padding = $self->param('exon_region_padding');
    my $exon_padded_start = $exon->start - $exon_region_padding;
    my $exon_padded_end = $exon->end + $exon_region_padding;
    if($exon_padded_end > $exon->slice->length) {
      $exon_padded_end = $exon->slice->length;
    }
    if($exon_padded_start < 1) {
      $exon_padded_start = 1;
    }

    my $dna_fragments = $compara_dba->get_DnaFragAdaptor->fetch_by_GenomeDB_and_name($source_genome_db,$exon->slice->seq_region_name);
    my $genomic_align_block_adaptor = $compara_dba->get_GenomicAlignBlockAdaptor;
    my $genomic_align_blocks = $genomic_align_block_adaptor->fetch_all_by_MethodLinkSpeciesSet_DnaFrag($mlss,$dna_fragments,$exon_padded_start,$exon_padded_end);
##    my $genomic_align_blocks = $genomic_align_block_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice($mlss, $exon_slice);

    foreach my $genomic_align_block (@$genomic_align_blocks) {
##      my $restricted_genomic_align_block = $genomic_align_block->restrict_between_reference_positions($exon_padded_start, $exon_padded_end);
      push(@{$all_target_genomic_aligns},@{$genomic_align_block->get_all_non_reference_genomic_aligns});
##      push(@{$all_target_genomic_aligns},@{$restricted_genomic_align_block->get_all_non_reference_genomic_aligns});
    }
  }

  my $unique_target_genomic_aligns = $self->unique_genomic_aligns($all_target_genomic_aligns);
  my @sorted_target_genomic_aligns = sort { $a->get_Slice->seq_region_name cmp $b->get_Slice->seq_region_name ||
                                            $a->get_Slice->start <=> $b->get_Slice->start
                                          } @{$unique_target_genomic_aligns};

  unless(scalar(@sorted_target_genomic_aligns)) {
    say "No align blocks so skipping transcript";
    next;
  }

  my $transcript_slices = $self->make_cluster_slices(\@sorted_target_genomic_aligns,$max_cluster_gap_length);
  unless($transcript_slices) {
    $self->throw("The sorted align blocks were not converted into transcript slices, something went wrong");
  }

  foreach my $transcript_slice (@{$transcript_slices}) {
    say "Created transcript slice: ".$transcript_slice->name."\n";
  }
  return $transcript_slices;
}


sub make_runnables {
  my ($self, $transcript_seq, $transcript_slices, $input_id, $annotation_features) = @_;
  my %parameters = %{$self->parameters_hash};
  if (not exists($parameters{-options}) and defined $self->OPTIONS) {
    $parameters{-options} = $self->OPTIONS;
  }
  if (not exists($parameters{-coverage_by_aligned}) and defined $self->COVERAGE_BY_ALIGNED) {
    $parameters{-coverage_by_aligned} = $self->COVERAGE_BY_ALIGNED;
  }

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript->new(
              -program  => $self->param('exonerate_path'),
              -analysis => $self->analysis,
              -query_type     => $self->QUERYTYPE,
              -calculate_coverage_and_pid => $self->param('calculate_coverage_and_pid'),
              -annotation_features => $annotation_features,
              %parameters,
              );
  $runnable->target_seqs($transcript_slices);
  $runnable->query_seqs([$transcript_seq]);
  $runnable->{'_transcript_id'} = $input_id;
  $self->runnable($runnable);
}


sub max_cluster_gap_length {
  my ($self,$transcript) = @_;

  # This sub will loop through the introns and decide on how long the max allowed value should be
  # This will be used to decide when to break up clusters on the same seq_region
  # The max value multiplier will change based on how long the biggest intron is
  # Have 50 as a min default value to handle small breaks in single exon genes
  my $longest_intron = 50;
  my $introns = $transcript->get_all_Introns();
  unless(scalar(@{$introns}) > 0) {
    return($longest_intron);
  }

  foreach my $intron (@{$introns}) {
    if($intron->length > $longest_intron) {
      $longest_intron = $intron->length;
    }
  }

  # This is only an intial guess, using orthologs would allow more accurate estimations for these
  # values in future. It would make sense to have different values for different clades
  if($longest_intron <= 10000) {
    return(int($longest_intron * 2));
  } elsif($longest_intron <= 50000) {
    return(int($longest_intron * 1.5));
  } else {
    return(int($longest_intron * 1.25))
  }
}


sub make_cluster_slices {
  my ($self,$genomic_aligns,$max_cluster_gap_length) = @_;

  my $slice_adaptor = $self->hrdb_get_con('target_dna_db')->get_SliceAdaptor;
  my $cluster_slices = [];
  my $previous_genomic_align = shift(@{$genomic_aligns});
  my $cluster_start = $previous_genomic_align->get_Slice->start();
  my $cluster_end = $previous_genomic_align->get_Slice->end();
  unless($previous_genomic_align) {
    return;
  }

  foreach my $current_genomic_align (@{$genomic_aligns}) {
    if($previous_genomic_align->get_Slice->seq_region_name ne $current_genomic_align->get_Slice->seq_region_name) {
      my $slice = $slice_adaptor->fetch_by_region($previous_genomic_align->get_Slice->coord_system_name,
                                                  $previous_genomic_align->get_Slice->seq_region_name, $cluster_start, $cluster_end);
      push(@{$cluster_slices},$slice);
      $cluster_start = $current_genomic_align->get_Slice->start();
      $cluster_end = $current_genomic_align->get_Slice->end();
      $previous_genomic_align = $current_genomic_align;
    } elsif($current_genomic_align->get_Slice->start - $previous_genomic_align->get_Slice->end > $max_cluster_gap_length) {
      my $slice = $slice_adaptor->fetch_by_region($previous_genomic_align->get_Slice->coord_system_name,
                                                  $previous_genomic_align->get_Slice->seq_region_name, $cluster_start, $cluster_end);
      push(@{$cluster_slices},$slice);
      $cluster_start = $current_genomic_align->get_Slice->start();
      $cluster_end = $current_genomic_align->get_Slice->end();
      $previous_genomic_align = $current_genomic_align;
    } else {
      $cluster_end = $current_genomic_align->get_Slice->end();
      $previous_genomic_align = $current_genomic_align;
    }
  } # end foreach my $current_genomic_align

  # push the final slice
  my $slice = $slice_adaptor->fetch_by_region($previous_genomic_align->get_Slice->coord_system_name,
                                              $previous_genomic_align->get_Slice->seq_region_name, $cluster_start, $cluster_end);
  push(@{$cluster_slices},$slice);
  return($cluster_slices);
}


sub unique_genomic_aligns {
  my ($self,$genomic_aligns) = @_;

  my $seen_id_hash = {};
  my $unique_genomic_aligns = [];
  foreach my $genomic_align (@{$genomic_aligns}) {
    unless($seen_id_hash->{$genomic_align->dbID}) {
      push(@{$unique_genomic_aligns},$genomic_align);
      $seen_id_hash->{$genomic_align->dbID} = 1;
    }
  }
  return($unique_genomic_aligns);
}



sub select_best_transcripts {
  my ($self,$preliminary_transcripts) = @_;

  my $selected_transcripts = [];
  my $best_score = 0;
  foreach my $preliminary_transcript (@{$preliminary_transcripts}) {
    my @tsfs = @{$preliminary_transcript->get_all_supporting_features};
    my $cov = $tsfs[0]->hcoverage;
    my $pid = $tsfs[0]->percent_id;
    my $combined_score = $cov + $pid;
    if($combined_score > $best_score) {
      $best_score = $combined_score;
    }
    $preliminary_transcript->{'combined_score'} = $combined_score;
  }

  unless($best_score) {
    $self->throw('Issue with calculating the best score from projection result set. This should not happen.');
  }

  # At this point we know the highest value in terms of combined percent identity and coverage. Now
  # we want all transcripts that fall within 5 percent of this value
  foreach my $preliminary_transcript (@{$preliminary_transcripts}) {
    my $combined_score = $preliminary_transcript->{'combined_score'};
    if($combined_score/$best_score >= 0.95) {
      push(@{$selected_transcripts},$preliminary_transcript);
    }
  }

  unless(scalar(@{$selected_transcripts})) {
    $self->throw("No transcripts selected in select_best_transcripts, something went wrong");
  }

  return($selected_transcripts);

}


sub filter_transcript {
  my ($self,$transcript) = @_;

  my $supporting_features = $transcript->get_all_supporting_features;
  my $supporting_feature = ${$supporting_features}[0];
  my $transcript_coverage = $supporting_feature->hcoverage;
  my $transcript_identity = $supporting_feature->percent_id;

  unless($transcript_identity >= $self->param_required('exonerate_percent_id') && $transcript_coverage >= $self->param_required('exonerate_coverage')) {
    print("Transcript failed coverage (".$self->param_required('exonerate_coverage').") and/or percent id (".$self->param_required('exonerate_percent_id').") filter, will not store. Transcript coverage and percent id are: ".$transcript_coverage." ".$transcript_identity."\n");
    return(1);
  }

  if($transcript->translation) {
    if($transcript->translation->seq =~ /\*/) {
      say "Transcript translation seq has stop, will not store";
      return(1);
    }

    my $original_transcript = $self->hrdb_get_con('source_transcript_db')->get_TranscriptAdaptor->fetch_by_dbID($transcript->{'_old_transcript_id'});
    my ($translation_coverage,$translation_identity) = align_proteins($original_transcript->translation->seq, $transcript->translation->seq);
    unless($translation_identity >= $self->param_required('exonerate_percent_id') && $translation_coverage >= $self->param_required('exonerate_coverage')) {
      print("Translation failed coverage (".$self->param_required('exonerate_coverage').") and/or percent id (".$self->param_required('exonerate_percent_id').") filter, will not store. Translation coverage and percent id are: ".$translation_coverage." ".$translation_identity."\n");
      return(1);
    }
  } # end if($transcript->translation)
  return(0);
}


sub create_annotation_features {
  my ($self,$transcript) = @_;

  my $cds_start  = $transcript->cdna_coding_start;
  my $cds_end = $transcript->cdna_coding_end;
  my $stable_id  = $transcript->stable_id.".".$transcript->version;

  my $annotation_feature = Bio::EnsEMBL::Feature->new(-seqname => $stable_id,
                                                      -strand  => 1,
                                                      -start   => $cds_start,
                                                      -end     => $cds_end);

 my $annotation_features->{$stable_id} = $annotation_feature;
 return($annotation_features);
}


sub retrieve_biotype {
  my ($self, $transcript) = @_;

  my $supporting_features = $transcript->get_all_supporting_features;
  my $supporting_feature = ${$supporting_features}[0];
  my $stable_id = $supporting_feature->hseqname;

  my $biotype = $self->param('_transcript_biotype')->{$stable_id};
  unless($biotype) {
    $self->throw("Failed to retieve biotype for output transcript. Should have been set in fetch_input");
  }
  return($biotype);
}


sub IIDREGEXP {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('IIDREGEXP',$value);
  }

  if ($self->param_is_defined('IIDREGEXP')) {
    return $self->param('IIDREGEXP');
  } else {
    return undef;
  }
}

sub OPTIONS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('OPTIONS',$value);
  }

  if ($self->param_is_defined('OPTIONS')) {
    return $self->param('OPTIONS');
  } else {
    return undef;
  }
}

sub COVERAGE_BY_ALIGNED {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('COVERAGE_BY_ALIGNED',$value);
  }

  if ($self->param_is_defined('COVERAGE_BY_ALIGNED')) {
    return $self->param('COVERAGE_BY_ALIGNED');
  } else {
    return undef;
  }
}

sub QUERYTYPE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('QUERYTYPE',$value);
  }

  if ($self->param_is_defined('QUERYTYPE')) {
    return $self->param('QUERYTYPE');
  } else {
    return undef;
  }
}

sub filter {
  my ($self, $val) = @_;
  if ($val) {
    $self->param('_transcript_filter',$val);
  }
  elsif ($self->param_is_defined('FILTER')) {
    $self->require_module($self->param('FILTER')->{OBJECT});
    $self->param('_transcript_filter', $self->param('FILTER')->{OBJECT}->new(%{$self->param('FILTER')->{FILTER_PARAMS}}));
  }
  if ($self->param_is_defined('_transcript_filter')) {
    return $self->param('_transcript_filter');
  }
  else {
    return;
  }
}

sub KILL_TYPE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->param('KILL_TYPE',$value);
  }

  if ($self->param_is_defined('KILL_TYPE')) {
    return $self->param('KILL_TYPE');
  } else {
    return undef;
  }
}

1;
