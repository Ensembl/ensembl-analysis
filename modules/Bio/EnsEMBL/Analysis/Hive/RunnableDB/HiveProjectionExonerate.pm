=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProjectionExonerate

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProjectionExonerate;

use warnings;
use strict;
use feature 'say';
use Data::Dumper;
use File::Spec::Functions qw(tmpdir catfile);

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(attach_Analysis_to_Gene attach_Slice_to_Gene);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(attach_Slice_to_Transcript empty_Transcript);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw (align_proteins);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my($self) = @_;
  $self->create_analysis;
  $self->analysis->logic_name($self->param('logic_name')) if ($self->param('logic_name'));
  my $input_ids = $self->param('iid');
  $self->param('exon_region_padding',100);
  $self->param('_branch_to_flow_on_fail',-3);
  $self->param('_branch_to_flow_on_filter',-4);
  # Define the dna dbs
  my $source_dna_dba = $self->hrdb_get_dba($self->param('source_dna_db'));
  my $target_dna_dba = $self->hrdb_get_dba($self->param('target_dna_db'));
  $self->hrdb_set_con($source_dna_dba,'source_dna_db');
  $self->hrdb_set_con($target_dna_dba,'target_dna_db');

  # Define the source transcript and target transcript dbs
  my $source_transcript_dba = $self->hrdb_get_dba($self->param('source_db'));
  my $target_transcript_dba = $self->hrdb_get_dba($self->param('target_db'));
  $target_transcript_dba->dnadb($target_dna_dba);
  $self->hrdb_set_con($source_transcript_dba,'source_transcript_db');
  $self->hrdb_set_con($target_transcript_dba,'target_transcript_db');

  # Define the compara db
  my $compara_dba = $self->hrdb_get_dba($self->param('compara_db'),undef,'Compara');
  $self->hrdb_set_con($compara_dba,'compara_db');

  # Get the genome db adpator
  my $genome_dba = $compara_dba->get_GenomeDBAdaptor;

  # Retrieve the production names for the query and target species
  my $source_species = $source_transcript_dba->get_MetaContainerAdaptor->get_production_name();
  my $target_species = $target_transcript_dba->get_MetaContainerAdaptor->get_production_name();

  my $source_genome_db = $genome_dba->fetch_by_core_DBAdaptor($source_transcript_dba);
  my $target_genome_db = $genome_dba->fetch_by_core_DBAdaptor($target_transcript_dba);

  ########
  # check that the default assembly for the query and target agrees
  # with that for the method_link_species_set GenomeDBs
  ########

  my $source_assembly = $source_genome_db->assembly;
  my $target_assembly = $target_genome_db->assembly;

  my ($source_assembly_version, $target_assembly_version);
  eval {
    $source_assembly_version = $source_transcript_dba->get_CoordSystemAdaptor->fetch_by_name('toplevel',$source_genome_db->assembly);
    $target_assembly_version = $target_transcript_dba->get_CoordSystemAdaptor->fetch_by_name('toplevel',$target_genome_db->assembly);
  };
  if($@) {
    $self->throw("Had trouble fetching coord systems for ".$source_genome_db->assembly . " and " .$target_genome_db->assembly . " from core dbs:\n".$@);
  }

  my $transcript_align_slices = {};
  my $mlss = $compara_dba->get_MethodLinkSpeciesSetAdaptor->fetch_by_method_link_type_GenomeDBs($self->param('method_link_type'),
                                                                                                [$source_genome_db,
                                                                                                $target_genome_db]);
  unless($mlss) {
    $self->throw("No MethodLinkSpeciesSet for :\n" .$self->param('method_link_type') . "\n" .$source_species . "\n" .$target_species);
  }

  $self->param('_transcript_biotype',{});
  foreach my $input_id (@$input_ids) {

    my $transcript = $source_transcript_dba->get_TranscriptAdaptor->fetch_by_dbID($input_id);
    my $transcript_seq;
    my $biotype = $transcript->biotype;
    my $stable_id = $transcript->stable_id.".".$transcript->version;
    my $annotation_features;
    $self->param('_transcript_biotype')->{$stable_id} = $biotype;

    say "Processing source transcript: ".$transcript->stable_id;
    if($self->QUERYTYPE eq 'protein') {
      unless($transcript->translation) {
        $self->warning("You have protein exonerate selected, but transcript does not have a translation, skipping");
        next;
      }
      $transcript_seq = $transcript->translation->seq;
    } else {
      if($self->param('generate_annotation_file')) {
        $annotation_features = $self->create_annotation_features($transcript);
      }
      $transcript_seq = $transcript->seq->seq;
    }
    my $transcript_slices = $self->process_transcript($transcript,$compara_dba,$mlss,$source_genome_db,$source_transcript_dba);
    my $transcript_header = $transcript->stable_id.'.'.$transcript->version;
    my $transcript_seq_object = Bio::Seq->new(-display_id => $transcript_header, -seq => $transcript_seq);
    $self->make_runnables($transcript_seq_object, $transcript_slices, $input_id, $annotation_features);
  } #close foreach input_id

}


sub run {
  my ($self) = @_;
  $self->runnable_failed(0);
  foreach my $runnable (@{$self->runnable}) {
    eval {
      $runnable->run;
    };

    if($@) {
      my $except = $@;
      $self->runnable_failed($runnable->{'_transcript_id'});
      $self->warning("Issue with running exonerate, will dataflow input id on fail branch. Exception:\n".$except);
      #$self->param('_branch_to_flow_on_fail', -3);
    } else {
      # Store original transcript id for realignment later. Should implement a cleaner solution at some point
      foreach my $output (@{$runnable->output}) {
        $output->{'_old_transcript_id'} = $runnable->{'_transcript_id'};
      }

      # If the transcript has bee projected to multiple places then select the best one in terms of combined coverage and
      # percent identity but also any that fall within 5 percent of this value
      my $preliminary_transcripts = $runnable->output;
      my $selected_transcripts;
      unless(scalar(@{$preliminary_transcripts})) {
        next;
      }

      if(scalar(@{$preliminary_transcripts}) == 1) {
        $selected_transcripts = $preliminary_transcripts;
      } else {
        $selected_transcripts = $self->select_best_transcripts($preliminary_transcripts);
      }
      $self->output($selected_transcripts);
    }
  }

  return 1;
}


sub write_output{
  my ($self) = @_;

  my $adaptor = $self->hrdb_get_con('target_transcript_db')->get_GeneAdaptor;
  my $slice_adaptor = $self->hrdb_get_con('target_transcript_db')->get_SliceAdaptor;

  my @output = @{$self->output};
  my $analysis = $self->analysis;

  foreach my $transcript (@output){

    my $slice_id = $transcript->start_Exon->seqname;
    my $slice = $slice_adaptor->fetch_by_name($slice_id);
    my $biotype = $self->retrieve_biotype($transcript);
    $transcript->biotype($biotype);
    attach_Slice_to_Transcript($transcript, $slice);

    if ($self->filter_transcript($transcript)) {
      # consider it as filtered out to be passed to the next analysis
      # via the _branch_to_flow_on_filter branch
      $self->filtered_out($transcript->{'_old_transcript_id'});
      $self->warning("The transcript has been filtered out. Dataflow input id on filtered out branch. Transcript ".$transcript->{'_old_transcript_id'});
      next;
    }

    my $gene = Bio::EnsEMBL::Gene->new();
    $gene->biotype($biotype);
    $gene->add_Transcript($transcript);
    $gene->slice($slice);
    attach_Analysis_to_Gene($gene, $analysis);
    $adaptor->store($gene);
  }
  my $output_hash = {};
  my $failure_branch_code = $self->param('_branch_to_flow_on_fail');
  my $failed_transcript_ids = $self->runnable_failed;
  if (scalar @$failed_transcript_ids ) {
    $output_hash->{'iid'} = $failed_transcript_ids;
    $self->dataflow_output_id($output_hash, $failure_branch_code);
  }

  my $filtered_out_hash = {};
  my $filtered_out_branch_code = $self->param('_branch_to_flow_on_filter');
  my $filtered_out_transcript_ids = $self->filtered_out();
  if (scalar @$filtered_out_transcript_ids ) {
    $filtered_out_hash->{'iid'} = $filtered_out_transcript_ids;
    $self->dataflow_output_id($filtered_out_hash,$filtered_out_branch_code);
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

sub filtered_out {
  my ($self,$filtered_out) = @_;
  if (not $self->param_is_defined('_filtered_out')) {
    $self->param('_filtered_out',[]);
  }
  if ($filtered_out) {
    push (@{$self->param('_filtered_out')},$filtered_out);
  }
  return ($self->param('_filtered_out'));
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
