=head1 LICENSE

 Copyright [2019-2020] EMBL-European Bioinformatics Institute

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

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::LowDivergenceProjection

=head1 SYNOPSIS

 my $runnable =
    Bio::EnsEMBL::Analysis::Runnable::LowDivergenceProjection->new();

 $runnable->run;
 my @results = $runnable->output;

=head1 DESCRIPTION


=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::Runnable::LowDivergenceProjection;

use warnings;
use strict;
use feature 'say';
use Data::Dumper;
use POSIX;
use List::Util qw(min max);

use Bio::EnsEMBL::Analysis::Runnable::Minimap2;
use Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(attach_Slice_to_Gene empty_Gene);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(compute_best_translation);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(set_alignment_supporting_features attach_Slice_to_Transcript attach_Analysis_to_Transcript calculate_exon_phases replace_stops_with_introns);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name align_nucleotide_seqs map_cds_location align_proteins);
use File::Spec;
use Bio::DB::HTS::Faidx;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use parent ('Bio::EnsEMBL::Analysis::Runnable::Minimap2Remap');


=head2 new

 Arg [DECOMPRESS]           : String as a command like 'gzip -c -'
 Arg [EXPECTED_ATTRIBUTES]  : String specify the attribute expected for the output, see STAR manual
 Description                : Creates a  object to align reads to a genome using STAR
 Returntype                 :
 Exceptions                 : Throws if WORKDIR does not exist
                              Throws if the genome has not been indexed

=cut

sub new {
  my ( $class, @args ) = @_;

  my $self = $class->SUPER::new(@args);
  my ($genome_index, $input_file, $paftools_path, $source_adaptor, $target_adaptor, $delete_input_file, $parent_genes, $parent_gene_ids, $gene_genomic_seqs_hash, $no_projection, $coverage_cutoff, $perc_id_cutoff, $extended_length_variation_cutoff, $anchor_coverage_cutoff, $anchor_perc_id_cutoff) = rearrange([qw (GENOME_INDEX INPUT_FILE PAFTOOLS_PATH SOURCE_ADAPTOR TARGET_ADAPTOR DELETE_INPUT_FILE PARENT_GENES PARENT_GENE_IDS GENE_GENOMIC_SEQS_HASH NO_PROJECTION COVERAGE_CUTOFF PERC_ID_CUTOFF EXTENDED_LENGTH_VARIATION_CUTOFF ANCHOR_COVERAGE_CUTOFF ANCHOR_PERC_ID_CUTOFF)],@args);
  $self->genome_index($genome_index);
  $self->input_file($input_file);
  $self->paftools_path($paftools_path);
  $self->source_adaptor($source_adaptor);
  $self->target_adaptor($target_adaptor);
  $self->delete_input_file($delete_input_file);
  $self->genes_to_process($parent_genes);
  $self->parent_gene_ids($parent_gene_ids);
  $self->gene_genomic_seqs_hash($gene_genomic_seqs_hash);
  $self->no_projection($no_projection);
  $self->coverage_cutoff($coverage_cutoff);
  $self->perc_id_cutoff($perc_id_cutoff);
  $self->extended_length_variation_cutoff($extended_length_variation_cutoff);
  $self->anchor_coverage_cutoff($anchor_coverage_cutoff);
  $self->anchor_perc_id_cutoff($anchor_perc_id_cutoff);
  return $self;
}

=head2 run

 Arg [1]    : None
 Description: Run Star to align reads to an indexed genome. The resulting output file will be stored in $self->output
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;

  my $leftover_genes = [];
  my $paf_file = $self->create_filename(undef,'paf');
  $self->files_to_delete($paf_file);

  my $genome_index  = $self->genome_index;
  my $input_file    = $self->input_file;
  my $options = $self->options;

  # TEST!!!!!!!!!!!!
#  my $source_adaptor = $self->source_adaptor();
#  my $source_gene_adaptor = $source_adaptor->get_GeneAdaptor();
#  my $source_genes = [$source_gene_adaptor->fetch_by_stable_id('ENSG00000175170')];

  # Process the batches first to map co-localised genes first
  my $source_genes = $self->genes_to_process();

  my $all_source_transcripts_id_hash = {};
  foreach my $source_gene (@$source_genes) {
    my $source_transcripts = $source_gene->get_all_Transcripts();
    foreach my $source_transcript (@$source_transcripts) {
      $all_source_transcripts_id_hash->{$source_transcript->dbID()} = $source_transcript;
    }
  }

  $self->process_gene_batches($source_genes,$genome_index,$all_source_transcripts_id_hash,$self->coverage_cutoff(),$self->perc_id_cutoff(),$self->extended_length_variation_cutoff(),$self->anchor_coverage_cutoff(),$self->anchor_perc_id_cutoff());
} # End run


sub process_gene_batches {
  my ($self,$source_genes,$genome_index,$all_source_transcripts_id_hash,$coverage_cutoff,$perc_id_cutoff,$extended_length_variation_cutoff,$anchor_coverage_cutoff,$anchor_perc_id_cutoff) = @_;

  my $source_adaptor = $self->source_adaptor();
  my $source_sequence_adaptor = $source_adaptor->get_SequenceAdaptor();
  my $target_adaptor = $self->target_adaptor();

  my $batched_input_genes = $self->batch_input_genes($source_genes);
  $self->print_batch_details($batched_input_genes);
  $self->calculate_anchors($batched_input_genes,$source_sequence_adaptor);
  $self->map_anchors($batched_input_genes,$genome_index,$anchor_coverage_cutoff,$anchor_perc_id_cutoff);

  # At this point we should have everything we need to project the batch
  my $target_genes_by_id = $self->project_batch_genes($batched_input_genes,$source_sequence_adaptor);
  my $target_genes = [];
  foreach my $id (keys(%$target_genes_by_id)) {
    my $genes = $target_genes_by_id->{$id};
    push(@$target_genes,@$genes);
  }

  my $source_genes_by_stable_id = $self->genes_by_stable_id($source_genes);
  my $target_genes_by_stable_id = $self->genes_by_stable_id($target_genes);


  my $problematic_genes = [];
  foreach my $gene (@$target_genes) {
    say "Checking if ".$gene->stable_id()." has problematic transcripts";
    $self->check_for_problematic_transcripts($gene,$source_genes_by_stable_id,$source_adaptor,$coverage_cutoff,$perc_id_cutoff,$extended_length_variation_cutoff);
    if($gene->{'is_problematic'}) {
      say "Gene has problematic transcripts";
      push(@$problematic_genes,$gene);
    } else {
      $self->output([$gene]);
    }
  }

  if(scalar(@$problematic_genes)) {
    say "Found ".scalar(@$problematic_genes)." problematic genes, will try to recover problematic transcripts";
  } else {
    say "No problematic genes found after projection";
  }

  my $output_genes = [];
  foreach my $gene (@$problematic_genes) {
    my $target_genomic_name = $gene->slice->seq_region_name();
    my $target_slice_adaptor = $target_adaptor->get_SliceAdaptor;
    my $target_parent_slice = $target_slice_adaptor->fetch_by_region('toplevel',$target_genomic_name);#$target_gene->slice->seq_region_Slice();

    say "Recovering transcripts for ".$gene->stable_id();
    my $updated_transcripts = $self->recover_transcripts($gene,$source_genes_by_stable_id,$source_adaptor,$target_parent_slice,$target_adaptor,$target_slice_adaptor,$coverage_cutoff,$perc_id_cutoff);
    if($updated_transcripts) {
      say "Gene has transcript set modified, will build new gene";
      my $source_gene = ${$source_genes_by_stable_id->{$gene->stable_id()}}[0];
      unless($source_gene) {
        $self->throw("Could not fund source gene for target gene with stable id: ".$gene->stable_id());
      }
      my $new_gene = $self->build_new_gene($source_gene,$updated_transcripts,$all_source_transcripts_id_hash,$target_parent_slice);
      if ($new_gene and $new_gene->get_all_Transcripts()) {
        push(@$output_genes,$new_gene);
      } else {
        say "The new gene has not been built because it does not contain any transcript."
      }
    } else {
      say "Transcript set was not updated, so will keep gene as-is";
      push(@$output_genes,$gene);
    }
  }

  $self->output($output_genes);
}


sub build_new_gene {
  my ($self,$source_gene,$target_transcripts,$all_source_transcripts_id_hash,$target_parent_slice) = @_;
  my $new_gene = Bio::EnsEMBL::Gene->new(-analysis => $self->analysis);
  $new_gene->{'parent_gene_id'} = $source_gene->dbID();
  $new_gene->{'parent_gene_versioned_stable_id'} = $source_gene->stable_id().".".$source_gene->version;
  $new_gene->stable_id($source_gene->stable_id);
  $new_gene->version($source_gene->version);
  $new_gene->biotype($source_gene->biotype);
  $new_gene->description($source_gene->description());
  # We take the strand of the first transcript to replicate what add_Transcript/recalculate_coordinates would do
  $new_gene->strand($target_transcripts->[0]->strand);

  # add source gene stable id as gene attribute
  my $parent_attribute = Bio::EnsEMBL::Attribute->new(-CODE => 'proj_parent_g',-VALUE => $source_gene->stable_id_version);
  $new_gene->add_Attributes($parent_attribute);
  $new_gene->{'is_new'} = 1;

  foreach my $target_transcript (@$target_transcripts) {
    my $source_transcript = $all_source_transcripts_id_hash->{$target_transcript->stable_id()};
    unless($source_transcript) {
      $self->throw("Couldn't retrieve source transcript for target transcript with stable id (dbID of source transcript): ".$target_transcript->stable_id());
    }

    $target_transcript->biotype($source_transcript->biotype);
    $target_transcript->{'parent_transcript_versioned_stable_id'} = $source_transcript->stable_id().".".$source_transcript->version();

    my $projected_transcript = $target_transcript->transfer($target_parent_slice);
    $projected_transcript->{annotation_method} = $target_transcript->{annotation_method};
    $projected_transcript->{cov} = $target_transcript->{cov};
    $projected_transcript->{perc_id} = $target_transcript->{perc_id};
    $projected_transcript->{parent_transcript_versioned_stable_id} = $target_transcript->{parent_transcript_versioned_stable_id};

    $self->set_transcript_description($target_transcript,$source_transcript);

    if ($projected_transcript->strand() == $new_gene->strand()) {
      $new_gene->add_Transcript($projected_transcript);
    } else {
      $self->warning('The target transcript (stable ID = '.$projected_transcript->stable_id().') strand and the new gene '.$new_gene->stable_id().' strand do not match.');
    }
  }
  return($new_gene);
}


sub check_for_problematic_transcripts {
  my ($self,$target_gene,$source_genes_by_stable_id,$source_adaptor,$coverage_cutoff,$perc_id_cutoff,$extended_length_variation_cutoff) = @_;

  my $source_transcript_adaptor = $source_adaptor->get_TranscriptAdaptor();
  my $problematic_genes = [];

  # The input is source genes, so only process target genes linked to the input source genes via a stable id
  my $source_gene = ${$source_genes_by_stable_id->{$target_gene->stable_id()}}[0];
  unless($source_gene) {
    return;
  }

  $target_gene->{'problematic_transcripts'} = {};
  my $source_transcripts = $source_gene->get_all_Transcripts();
  my $source_transcripts_id_hash = {};

  # This is used later to keep track of any missing transcripts
  foreach my $source_transcript (@$source_transcripts) {
    $source_transcripts_id_hash->{$source_transcript->stable_id()} = $source_transcript;
  }

  my $target_transcripts = $target_gene->get_all_Transcripts();
  foreach my $target_transcript (@$target_transcripts) {
    my $description = $target_transcript->description();
    $description =~ /parent_transcript=((ENS|GeneID_|LOC|XR_|XM_).+)\.\d*;mapping_coverage=(\d+\.\d+);mapping_identity=(\d+\.\d+)/;
    my $stable_id = $1;
    my $coverage = $3;
    my $perc_id = $4;
    $target_transcript->{'cov'} = $coverage;
    $target_transcript->{'perc_id'} = $perc_id;
    unless($stable_id and defined($coverage) and defined($perc_id)) {
      $self->throw("Issue finding info for parent stable id and coverage/identity from transcript description. Description: ".$description);
    }

    # This overwrites the transcript ref, so if the value is 1 we won't consider it missing

    my $source_transcript = $source_transcripts_id_hash->{$stable_id};
    $source_transcripts_id_hash->{$stable_id} = 1;
    my $length_diff = $target_transcript->length() - $source_transcript->length();

    my $problematic_transcript = 0;
    my $problem_string = "Problematic transcript ".$stable_id.", issue:";
    if($length_diff >= 0 and ($length_diff/$source_transcript->length() > $extended_length_variation_cutoff)) {
      $problem_string .= " length";
      $problematic_transcript = 1;
    }

    if($coverage < $coverage_cutoff) {
      $problem_string .= " coverage";
      $problematic_transcript = 1;
    }

    if($perc_id < $perc_id_cutoff) {
      $problem_string .= " perc_id";
      $problematic_transcript = 1;
    }

    if($problematic_transcript) {
      say $problem_string;
      $target_gene->{'is_problematic'} = 1;
      $target_gene->{'problematic_transcripts'}->{$stable_id} = [$source_transcript,$target_transcript];
    }
  } # foreach my $target_transcript

  # This adds in any missing transcripts as problematic transcripts
  foreach my $source_transcript_stable_id (keys(%$source_transcripts_id_hash)) {
    if($source_transcripts_id_hash->{$source_transcript_stable_id} == 1) {
      next;
    }
    say "Problematic transcript ".$source_transcript_stable_id.", issue: missing";
    $target_gene->{'problematic_transcripts'}->{$source_transcript_stable_id} = [$source_transcripts_id_hash->{$source_transcript_stable_id}];
  }
}


sub recover_transcripts {
  my ($self,$target_gene,$source_genes_by_stable_id,$source_adaptor,$target_parent_slice,$target_adaptor,$target_slice_adaptor,$coverage_cutoff,$perc_id_cutoff) = @_;


  #my $coverage_cutoff = 98;
  #my $perc_id_cutoff = 99;
  my $region_scaling = 0.1;
  my $min_buffer = 500;
  my $source_gene_adaptor = $source_adaptor->get_TranscriptAdaptor();
  my $problematic_genes = [];

  # We now have a list of missing and problematic transcripts. We want to look for these, but in the region
  # the gene seems to lie in. The issue is that the transcripts that are present may not represent the true
  # boundaries of the gene, even if there are fully mapped transcript they might only represent part of the
  # are covered. To deal with this, we'll pad the region the gene covers with whatever the difference in length
  # is either side. This is not perfect, but should cover most of the issues
  my $stable_id = $target_gene->stable_id();
  my $source_gene = ${$source_genes_by_stable_id->{$target_gene->stable_id()}}[0];
  my $source_transcripts = $source_gene->get_all_Transcripts();

  my $source_gene_length = $source_gene->length();
  my $target_gene_length = $target_gene->length();
  my $length_diff = $source_gene_length - $target_gene_length;
  my $region_buffer = 0;
  if($length_diff >= 0) {
    $region_buffer = $length_diff + ($length_diff * $region_scaling) + $min_buffer;
  } else {
    $region_buffer = $min_buffer;
  }

  $region_buffer = ceil($region_buffer);
  my $target_region_left_boundary = $target_gene->seq_region_start() - $region_buffer;
  my $target_region_right_boundary = $target_gene->seq_region_end() + $region_buffer;
  if($target_region_left_boundary <= 0) {
    $target_region_left_boundary = 1;
  }

  if($target_region_right_boundary > $target_gene->slice->length()) {
    $target_region_right_boundary = $target_gene->slice->length();
  }

  my $target_genomic_name = $target_parent_slice->seq_region_name();
  my $target_genomic_start = $target_region_left_boundary;
  my $target_genomic_end = $target_region_right_boundary;

  say "Mapping gene ".$source_gene->stable_id()." to ".$target_genomic_start.":".$target_genomic_end.":".$target_gene->strand.":".$target_genomic_name;

  # Note that we are going to use the forward strand regardless of what strand the PAF hit is on here because minimap2/exonerate assume the region is on the forward strand
  my $target_strand = $target_gene->strand();
  my $target_sequence_adaptor = $target_adaptor->get_SequenceAdaptor;
  my $target_genomic_seq = ${ $target_sequence_adaptor->fetch_by_Slice_start_end_strand($target_parent_slice, $target_genomic_start, $target_genomic_end, 1) };
  my $target_region_slice = $target_slice_adaptor->fetch_by_region('toplevel', $target_genomic_name, $target_genomic_start, $target_genomic_end, 1);
  say "Target slice identified to search for transcripts in after first pass: ".$target_region_slice->name();
  my $target_genomic_fasta = ">".$target_genomic_name."\n".$target_genomic_seq;
  my $target_genome_file = $self->write_input_file([$target_genomic_fasta]);

  say "Projecting gene: ".$source_gene->stable_id();

  my $best_transcripts_by_id = {};
  my $source_transcript_id_hash = {}; # A hash to organise transcripts by dbID (note that the dbID is somewhat confusingly saved in the stable id field here, safer than using the realstable id)
  say "Source transcript list:";
  my $problematic_source_transcripts = [];
  foreach my $source_transcript (@$source_transcripts) {
    my $source_transcript_id = $source_transcript->dbID();
    $source_transcript_id_hash->{$source_transcript_id} = $source_transcript;

    unless($target_gene->{'problematic_transcripts'}->{$source_transcript->stable_id()}) {
      next;
    }

    say "  ".$source_transcript->stable_id();
    push(@$problematic_source_transcripts,$source_transcript);
  }

  my $max_intron_size = $self->calculate_max_intron_size($problematic_source_transcripts);
  $max_intron_size = ceil($max_intron_size);

  say "FERGAL DEBUG TARGETGS: ".$target_genomic_start;
  my $minimap_transcripts_by_id = $self->map_gene_minimap($problematic_source_transcripts,$target_genomic_start,$target_region_slice,$target_strand,
                                                          $target_genome_file,$source_transcript_id_hash,$max_intron_size,$target_adaptor,$target_slice_adaptor,
                                                          $best_transcripts_by_id);
  $self->print_transcript_stats($minimap_transcripts_by_id,'minimap local');
  $self->update_best_transcripts($best_transcripts_by_id,$minimap_transcripts_by_id);
  my $exonerate_transcripts_by_id = $self->map_gene_exonerate($problematic_source_transcripts,$target_genomic_start,$target_region_slice,$target_strand,
                                                              $target_genome_file,$source_transcript_id_hash,$max_intron_size,$target_adaptor,
                                                              $target_slice_adaptor,$best_transcripts_by_id);
  $self->print_transcript_stats($exonerate_transcripts_by_id,'exonerate local');
  $self->update_best_transcripts($best_transcripts_by_id,$exonerate_transcripts_by_id);

  $self->set_cds_sequences($best_transcripts_by_id,$source_transcript_id_hash);
  $self->qc_cds_sequences($best_transcripts_by_id,$source_transcript_id_hash);
  $self->set_transcript_descriptions($best_transcripts_by_id,$source_transcript_id_hash);

  # Now need to loop through the best recovered transcripts and overwrite the original ones (if present)
  my $final_transcripts = [];
  my $target_transcripts = $target_gene->get_all_Transcripts();
  foreach my $target_transcript (@$target_transcripts) {
    # I cannot remember why the target transcripts had dbID in the stable id field at this point, but this is something that should be refactored
    my $target_transcript_stable_id = $source_transcript_id_hash->{$target_transcript->stable_id()}->stable_id();
    unless($target_transcript_stable_id) {
      $self->throw("Failed to retrieve a source stable id transcript with dbID ".$target_transcript->stable_id());
    }
    unless($target_gene->{'problematic_transcripts'}->{$target_transcript_stable_id}) {
      push(@$final_transcripts,$target_transcript);
    }
  }

  # Problem: best_transcripts_by_id uses source_transcript->dbID for id, problematic_transcripts uses stable id

  my $gene_modified = 0;
  foreach my $source_transcript (@$source_transcripts) {
    if($target_gene->{'problematic_transcripts'}->{$source_transcript->stable_id()}) {
      my $problematic_transcripts = $target_gene->{'problematic_transcripts'}->{$source_transcript->stable_id()};
      # If there's only one, then this was a missing transcript, so just take whatever is in the best transcript hash
      if(scalar(@$problematic_transcripts) == 1) {
        my $best_transcript = $best_transcripts_by_id->{${$problematic_transcripts}[0]->dbID()};
        if($best_transcript) {
          $best_transcript->{'is_new'} = 1;
          $best_transcript->{'parent_transcript_versioned_stable_id'} = $source_transcript->stable_id().".".$source_transcript->version();
          say "  Adding missing transcript ".$best_transcript->stable_id();
          push(@$final_transcripts,$best_transcript);
          $gene_modified = 1;
        }
      } else {
        # In this case there's an existing projected transcript
        my $best_transcript = $best_transcripts_by_id->{${$problematic_transcripts}[0]->dbID()};
        unless($best_transcript){
          say "  Warning: did not find a best transcript for ".${$problematic_transcripts}[0]->stable_id();
        }
        $best_transcript->{'is_new'} = 1;
        $best_transcript->{'parent_transcript_versioned_stable_id'} = $source_transcript->stable_id().".".$source_transcript->version();
        my $projected_transcript = ${$problematic_transcripts}[1];
        my $best_coverage = $best_transcript->{'cov'} || 0;
        my $best_perc_id = $best_transcript->{'perc_id'} || 0;
        my $best_total = $best_coverage + $best_perc_id;
        my $projected_transcript_coverage = $projected_transcript->{'cov'};
        my $projected_transcript_perc_id = $projected_transcript->{'perc_id'};
        my $projected_total = $projected_transcript_coverage + $projected_transcript_perc_id;
        if($best_total >= $projected_total or ($best_perc_id >= $perc_id_cutoff and $best_coverage >= $coverage_cutoff)) {
          say "  Replacing projected transcript ".$best_transcript->stable_id()." with best mapped transcript";
          push(@$final_transcripts,$best_transcript);
          $gene_modified = 1;
        } else {
          say "  Keeping projected transcript ".$projected_transcript->stable_id()." over best mapped transcript";
          push(@$final_transcripts,$projected_transcript);
        }
      }
    }
  } # End foreach my $source_transcript


  if(scalar(@$final_transcripts) < scalar(@$target_transcripts)) {
    $self->throw("After recoverey have ".scalar(@$final_transcripts)." transcripts, which is less than the initial amount of ".scalar(@$target_transcripts));
  }

  say "After recoverey have ".scalar(@$final_transcripts)." transcripts, initially had ".scalar(@$target_transcripts);

  unless($gene_modified) {
    return;
  } else {
    return($final_transcripts);
  }
}

sub batch_input_genes {
  my ($self,$genes) = @_;

  my $max_batch_span = 100000;
  my $source_flank = 1000;
  my $batched_input_genes = {};
  my $batch_id = 1;

  my $current_batch_start = 0;
  my $current_batch_end = 0;
  my $current_batch_slice_name;
  my $current_batch_slice;
  my $current_batch_genes = [];
  for(my $i=0; $i<scalar(@$genes); $i++) {
    my $gene = ${$genes}[$i];

    unless($current_batch_slice) {
      $current_batch_slice_name = $gene->slice->seq_region_name();
      $current_batch_slice = $gene->slice->seq_region_Slice();
      $current_batch_start = $gene->seq_region_start();
      $current_batch_end = $gene->seq_region_end();
    }

    if((($gene->seq_region_end - $current_batch_start + 1) > $max_batch_span and $gene->seq_region_end() > $current_batch_end) or $gene->slice->seq_region_name ne $current_batch_slice_name) {
      $batched_input_genes->{$batch_id}->{'genes'} = $current_batch_genes;
      $batched_input_genes->{$batch_id}->{'start'} = $current_batch_start;
      $batched_input_genes->{$batch_id}->{'end'} = $current_batch_end;
      $batched_input_genes->{$batch_id}->{'slice_name'} = $current_batch_slice_name;
      $batched_input_genes->{$batch_id}->{'slice'} = $current_batch_slice;
      $batch_id++;
      $current_batch_slice_name = $gene->slice->seq_region_name();
      $current_batch_slice = $gene->slice->seq_region_Slice();
      $current_batch_start = $gene->seq_region_start();
      $current_batch_end = $gene->seq_region_end();
      $current_batch_genes = [$gene];
      next;
    }

    if($gene->seq_region_end() > $current_batch_end) {
      $current_batch_end = $gene->seq_region_end();
    }

    push(@$current_batch_genes,$gene);
  }

  if(scalar(@$current_batch_genes)) {
    $batched_input_genes->{$batch_id}->{'genes'} = $current_batch_genes;
    $batched_input_genes->{$batch_id}->{'start'} = $current_batch_start;
    $batched_input_genes->{$batch_id}->{'end'} = $current_batch_end;
    $batched_input_genes->{$batch_id}->{'slice_name'} = $current_batch_slice_name;
    $batched_input_genes->{$batch_id}->{'slice'} = $current_batch_slice;
  }

  # Add flanking region to each batch
  foreach my $id (keys(%$batched_input_genes)) {
    my $batch = $batched_input_genes->{$id};
    my $slice = $batch->{'slice'};
    $batch->{'start'} = $batch->{'start'} - $source_flank;
    $batch->{'end'} = $batch->{'end'} + $source_flank;
    if($batch->{'start'} < 1) {
      $batch->{'start'} = 1;
    }

    if($batch->{'end'} > $slice->length()) {
      $batch->{'end'} = $slice->length();
    }

    $batch->{'length'} = $batch->{'end'} - $batch->{'start'} + 1;
    $batch->{'midpoint'} = $batch->{'start'} + ceil($batch->{'length'}/2);
  }

  return($batched_input_genes);
}


sub print_batch_details {
  my ($self,$batched_input_genes) = @_;

  foreach my $id (keys(%$batched_input_genes)) {
    my $batch = $batched_input_genes->{$id};
    say "Batch info: ID: ".$id.", Slice: ".$batch->{'slice_name'}.", Start: ".$batch->{'start'}.", End: ".$batch->{'end'}.", Midpoint: ".$batch->{'midpoint'}.", Length: ".$batch->{'length'};
    my $genes = $batch->{'genes'};
    my @gene_db_ids;
    foreach my $gene (@$genes) {
      say "  Gene: ID: ".$gene->stable_id().", Slice: ".$gene->seq_region_name().", Start: ".$gene->seq_region_start().", End: ".$gene->seq_region_end().", Strand: ".$gene->strand();
      push(@gene_db_ids,$gene->dbID());
    }
    my $gene_id_string = "[".join(',',@gene_db_ids)."]";
    say "  Gene dbIDs: ".$gene_id_string;
  }
}


sub calculate_anchors {
  my ($self,$batched_input_genes,$sequence_adaptor) = @_;

  # This offset will be added to either side of an anchor point, so effectively the anchor will be double this size
  my $batch_offset = 5000;

  foreach my $id (keys(%$batched_input_genes)) {
    say "Finding Anchors for batch ID: ".$id;
    my $batch = $batched_input_genes->{$id};
    my $batch_start = $batch->{'start'};
    my $batch_midpoint = $batch->{'midpoint'};
    my $batch_end = $batch->{'end'};
    my $batch_length = $batch_end - $batch_start + 1;
    my $batch_slice_name = $batch->{'slice_name'};
    my $batch_slice = $batch->{'slice'};

    my $anchor_left_start = $batch_start - $batch_offset;
    my $anchor_left_end = $batch_start + $batch_offset;

    if($anchor_left_start < 1) {
      $anchor_left_start = 1;
    }

    if($anchor_left_end > $batch_slice->length()) {
      $anchor_left_end = $batch_slice->length();
    }

    my $anchor_left_seq = ${ $sequence_adaptor->fetch_by_Slice_start_end_strand($batch_slice,$anchor_left_start,$anchor_left_end,1) };

    my $anchor_middle_start = $batch_midpoint - $batch_offset;
    my $anchor_middle_end = $batch_midpoint + $batch_offset;

    if($anchor_middle_start < 1) {
      $anchor_middle_start = 1;
    }

    if($anchor_middle_end > $batch_slice->length()) {
      $anchor_middle_end = $batch_slice->length();
    }

    my $anchor_middle_seq = ${ $sequence_adaptor->fetch_by_Slice_start_end_strand($batch_slice,$anchor_middle_start,$anchor_middle_end,1) };
    my $anchor_right_start = $batch_end - $batch_offset;
    my $anchor_right_end = $batch_end + $batch_offset;

    if($anchor_right_start < 1) {
      $anchor_right_start = 1;
    }

    if($anchor_right_end > $batch_slice->length()) {
      $anchor_right_end = $batch_slice->length();
    }

    my $anchor_right_seq = ${ $sequence_adaptor->fetch_by_Slice_start_end_strand($batch_slice,$anchor_right_start,$anchor_right_end,1) };

    $batch->{'anchor_seqs'} = [[$anchor_left_start,$anchor_left_end,$anchor_left_seq],
                               [$anchor_middle_start,$anchor_middle_end,$anchor_middle_seq],
                               [$anchor_right_start,$anchor_right_end,$anchor_right_seq]];
  } # foreach my $id
}


sub map_anchors {
  my ($self,$batched_input_genes,$genome_index,$anchor_coverage_cutoff,$anchor_perc_id_cutoff) = @_;

  say "Creating fasta records";
  my $fasta_records = [];
  foreach my $id (keys(%$batched_input_genes)) {
    my $batch = $batched_input_genes->{$id};
    my $anchors = $batch->{'anchor_seqs'};
    my $fasta_record = ">batch_".$id."_al\n".${$anchors}[0][2]."\n>batch_".$id."_am\n".${$anchors}[1][2]."\n>batch_".$id."_ar\n".${$anchors}[2][2]."\n";
#    say "DEBUG:\n".$fasta_record;
    push(@$fasta_records,$fasta_record);
  }

  my $input_file = $self->write_input_file($fasta_records);
  my $paf_file = $input_file.".paf";
  my $minimap2_command = $self->program()." --cs --secondary=yes -x map-ont -N 10 ".$genome_index." ".$input_file." > ".$paf_file;
  system($minimap2_command);

  say "DEBUG: ".$minimap2_command;

  unless(-e $paf_file) {
    $self->throw("Could not find paf file from running minimap on anchor seqs");
  }

  open(IN,$paf_file);
  while(<IN>) {
    my $line = $_;
    print "PAF: ".$line;
    unless($line =~ /^batch\_(\d+)\_(.+)$/) {
      next;
    }
    my $batch_id = $1;
    my $result_line = $2;
    unless($batched_input_genes->{$batch_id}->{'paf_results'}) {
      $batched_input_genes->{$batch_id}->{'paf_results'} = [$result_line];
    } else {
      push(@{$batched_input_genes->{$batch_id}->{'paf_results'}},$result_line);
    }
  }
  close IN;

  foreach my $id (keys(%$batched_input_genes)) {
    say "Mapping Anchors for batch ID: ".$id;
    my $batch = $batched_input_genes->{$id};
    $self->calculate_target_regions($batch,$anchor_coverage_cutoff,$anchor_perc_id_cutoff);
  }
}


sub calculate_target_regions {
  my ($self,$batch,$anchor_coverage_cutoff,$anchor_perc_id_cutoff) = @_;

  my $anchors = $batch->{'anchor_seqs'};
  my $paf_results = $batch->{'paf_results'};
  unless($paf_results and scalar(@$paf_results ) > 1) {
    return;
  }

  # First calculate the source distances between the chosen anchor mid
  my $source_anchor1_midpoint = ${$anchors}[0][1] + ceil((${$anchors}[0][1]-${$anchors}[0][0])/2);
  my $source_anchor2_midpoint = ${$anchors}[1][1] + ceil((${$anchors}[1][1]-${$anchors}[1][0])/2);
  my $source_anchor3_midpoint = ${$anchors}[2][1] + ceil((${$anchors}[2][1]-${$anchors}[2][0])/2);


  say "  Source anchor 1 midpoint: ".$source_anchor1_midpoint;
  say "  Source anchor 2 midpoint: ".$source_anchor2_midpoint;
  say "  Source anchor 3 midpoint: ".$source_anchor3_midpoint;

  my $anchors_left_res = [];
  my $anchors_middle_res = [];
  my $anchors_right_res = [];

  foreach my $line (@$paf_results) {
    my @eles = split("\t",$line);
    $eles[2]++;
    $eles[3]++;
    $eles[7]++;
    $eles[8]++;

    my $anchor_name = $eles[0];
    my $source_length = $eles[1];
    my $source_hit_start = $eles[2];
    my $source_hit_end = $eles[3];
    my $target_strand = $eles[4];
    my $target_genomic_name = $eles[5];
    my $target_genomic_length = $eles[6];
    my $target_genomic_start = $eles[7];
    my $target_genomic_end = $eles[8];
    my $alignment_identities = $eles[9];
    my $alignment_length = $eles[10];
    my $source_hit_midpoint = ceil($source_length/2);
    my $source_hit_length = $source_hit_end - $source_hit_start + 1;
    my $failed_coverage = 0;
    my $failed_identity = 0;

    # Skip hits that are below 99 percent identity, should probably reconsider this in terms of more divergent stuff like mouse strains
    unless($alignment_identities/$source_hit_length >= $anchor_perc_id_cutoff) {
      say "Hit fails the identity cutoff: (".$target_genomic_start."/".$target_genomic_end.")";
      $failed_identity = 1;
    }

    unless($source_hit_length >= ($source_length * $anchor_coverage_cutoff)) {
      say "Hit fails the coverage cutoff: (".$target_genomic_start."/".$target_genomic_end."), Hit start: ".$source_hit_start.", Hit end: ".$source_hit_end;
      $failed_coverage = 1;
    }

    my $target_midpoint_equivalent = $self->calculate_target_midpoint($source_hit_start,$source_hit_end,$source_hit_midpoint,$target_genomic_start,$target_genomic_end,$target_strand,$target_genomic_length);
    if($anchor_name eq 'al') {
      push(@{$anchors_left_res},[$target_midpoint_equivalent,\@eles,$failed_coverage,$failed_identity]);
    } elsif($anchor_name eq 'am') {
      push(@{$anchors_middle_res},[$target_midpoint_equivalent,\@eles,$failed_coverage,$failed_identity]);
    } elsif($anchor_name eq 'ar') {
      push(@{$anchors_right_res},[$target_midpoint_equivalent,\@eles,$failed_coverage,$failed_identity]);
    } else {
      $self->throw('Found an unexpected anchor name: '.$anchor_name);
    }
    say "Calculated target genomic midpoint for: Anchor type: ".$anchor_name.", Target region: ".$target_genomic_name.", Target start: ".$target_genomic_start.", Target end: ".$target_genomic_end.", Target strand: ".$target_strand.", Target midpoint: ".$target_midpoint_equivalent;
  }

  # This really only works if all three are found and even then it's completely biased towards the left anchor, for the moment just going to go with this, but
  # it means a separate method if we want to take advantage of partial matches
  my $anchor_sets = [];
  if(scalar(@{$anchors_left_res}) and scalar(@{$anchors_middle_res}) and scalar(@{$anchors_right_res})) {
    $anchor_sets = $self->get_best_anchor_sets($source_anchor1_midpoint,$source_anchor2_midpoint,$source_anchor3_midpoint,$anchors_left_res,$anchors_middle_res,$anchors_right_res);
  }

  foreach my $anchor_set (@$anchor_sets) {
    my @a1_eles = @{${$anchor_set}[0][1]};
    my @a2_eles = @{${$anchor_set}[1][1]};
    my @a3_eles = @{${$anchor_set}[2][1]};
    my $dist = ${$anchor_set}[3];
    my $issue_count = ${$anchor_set}[4];
    say "Chosen anchor: (".$a1_eles[7].":".$a1_eles[8].") (".$a2_eles[7].":".$a2_eles[8].") (".$a3_eles[7].":".$a3_eles[8]."). Dist: ".$dist.", Issue count: ".$issue_count;
  }

  if(scalar(@$anchor_sets)) {
    $self->collapse_anchor_sets($anchor_sets,$batch);
    return;
  }

  say "No suitable anchor set found from initial results! Looking at pairwise anchors";
  $anchor_sets = $self->get_pairwise_sets($source_anchor1_midpoint,$source_anchor2_midpoint,$source_anchor3_midpoint,$anchors_left_res,$anchors_middle_res,$anchors_right_res);

  foreach my $anchor_set (@$anchor_sets) {
    my @a1_eles = @{${$anchor_set}[0][1]};
    my @a2_eles = @{${$anchor_set}[1][1]};
    my $dist = ${$anchor_set}[2];
    say "Chosen anchor: (".$a1_eles[7].":".$a1_eles[8].") (".$a2_eles[7].":".$a2_eles[8].") Dist: ".$dist;
  }

  if(scalar(@$anchor_sets)) {
    $self->collapse_anchor_sets($anchor_sets,$batch);
    return;
  }

  say "No anchor set found!!!";
  return;
}


sub collapse_anchor_sets {
  my ($self,$anchor_sets,$batch) = @_;

  my $clustered_indices = {};
  my $current_boundary_sets = [];
  my $merged_anchor_sets = [];
  my $anchor_info = {};
  say "Collapsing anochor sets";
  my $anchor_set_id = 0;
  foreach my $anchor_set (@$anchor_sets) {
    my $anchor_set1 = ${$anchor_set}[0];

    # Should probably standarise the anchor sets prior to this, but at the moment if there are all three anchors
    # then we should pick the first and third, otherwise pick the first and second
    my $anchor_set2;
    if(scalar(@$anchor_set) > 3) {
      $anchor_set2 = ${$anchor_set}[2];
    } else {
      $anchor_set2 = ${$anchor_set}[1];
    }

    my @a1_eles = @{${$anchor_set1}[1]};
    my @a2_eles = @{${$anchor_set2}[1]};
    my $left_boundary = min($a1_eles[7],$a2_eles[7]);
    my $right_boundary = max($a1_eles[8],$a2_eles[8]);
    # These two have to be the same in terms of a1/a2 at this point
    my $strand = $a1_eles[4];
    my $region = $a1_eles[5];
    say "  Anchor set boundaries: (".$left_boundary.":".$right_boundary.") ".$region." ".$strand;
    $anchor_info->{$region}->{$strand}->{$anchor_set_id}->{'start'} = $left_boundary;
    $anchor_info->{$region}->{$strand}->{$anchor_set_id}->{'end'} = $right_boundary;
    $anchor_set_id++;
  }

  # At this point every anchor is in the anchor_info hash, divided by region, strand and then id
  # So then going through a particular region/strand should allow the sorting of the ids, which will
  # then allow them to be clustered into final anchor based on overlap. There are various potential issues
  # with doing this, particularly making giant regions, but it's worth trying this without specialised cutoffs
  foreach my $region (keys(%$anchor_info)) {
    my $region_hash = $anchor_info->{$region};
    foreach my $strand (keys(%$region_hash)) {
      my $strand_hash = $region_hash->{$strand};
      my $anchors = [];
      foreach my $id (keys(%$strand_hash)) {
        push(@$anchors,$strand_hash->{$id});
      }
      my @sorted_anchors = sort { $a->{'start'} <=> $b->{'start'} } @{$anchors};
      for(my $i=0; $i<scalar(@sorted_anchors)-1; $i++) {
        my $a1 = $sorted_anchors[$i];
        unless($a1) {
          next;
        }
        for(my $j=$i+1; $j<scalar(@sorted_anchors); $j++) {
          my $a2 = $sorted_anchors[$j];
          unless($a2) {
            next;
          }
          if($self->coords_overlap($a1->{'start'},$a1->{'end'},$a2->{'start'},$a2->{'end'})) {
            say "Found overlapping anchors, will merge";
            say "Original boundaries: (".$a1->{'start'}."/".$a1->{'end'}.") (".$a2->{'start'}."/".$a2->{'end'}.")";
            my $region_start = min($a1->{'start'},$a2->{'start'});
            my $region_end = max($a1->{'end'},$a2->{'end'});
            $a1->{'start'} = $region_start;
            $a1->{'end'} = $region_end;
            $sorted_anchors[$j] = 0;
            say "  Merged anchor boundaries: (".$a1->{'start'}."/".$a1->{'end'}.")";
          }
        }
      }
      foreach my $merged_anchor (@sorted_anchors) {
        if($merged_anchor) {
          say "  Final anchor boundaries: (".$merged_anchor->{'start'}.":".$merged_anchor->{'end'}.") ".$region." ".$strand;
          unless($batch->{'target_regions'}) {
            $batch->{'target_regions'} = [[$region,$strand,$merged_anchor->{'start'},$merged_anchor->{'end'}]];
          } else {
            push(@{$batch->{'target_regions'}},[$region,$strand,$merged_anchor->{'start'},$merged_anchor->{'end'}]);
          }
        }
      }
    }
  }
}


sub get_best_anchor_sets {
  my ($self,$source_anchor1_midpoint,$source_anchor2_midpoint,$source_anchor3_midpoint,$anchors_left_res,$anchors_middle_res,$anchors_right_res) = @_;
  my $best_anchors = [];

  my $distance_limit = 250000;
  my $source_dist1 = abs($source_anchor2_midpoint - $source_anchor1_midpoint);
  my $source_dist2 = abs($source_anchor3_midpoint - $source_anchor2_midpoint);
  my $source_dist3 = abs($source_anchor3_midpoint - $source_anchor1_midpoint);

  # Could consider changing this from skipping things with issues to an all paths calc and then if there are paths with no issues, pick only those
  # Otherwise report the one with the least issues
  # What are issues:
  # - fail identity
  # - fail coverage
  # - not colinear
  # - large dist calculated (or dist does not make sense, e.g. one anchor is on a different region)
  # - strand issues
  #
  # It seems like these numbers should be calculated and then a sanity check should be done on any potential candidates
  # For example, we would expect sets to have at least two anchors on the same region and the same strand and we would
  # remove the problematic anchor
  # This also might miss the cases where there's single issue with the true anchor, like in the case of a large insertion
  # messing up the distance
  my $existing_anchor_coords = {};
  my $anchor_sets = [];
  for(my $e=0; $e<scalar(@$anchors_left_res); $e++) {
    my $issue_tracker = 0;
    my $a1 = ${$anchors_left_res}[$e];
    my $a1_mid = ${$a1}[0];
    my @a1_eles = @{${$a1}[1]};
    my $a1_strand = $a1_eles[4];
    my $a1_region = $a1_eles[5];
    my $a1_start = $a1_eles[7];
    my $a1_end = $a1_eles[8];
    my $a1_failed_coverage = ${$a1}[2];
    my $a1_failed_identity = ${$a1}[3];
    if($a1_failed_coverage) {
      $issue_tracker++;
    }

    if($a1_failed_identity) {
      $issue_tracker++;
    }

    for(my $f=0; $f<scalar(@$anchors_middle_res); $f++) {
      my $a2 = ${$anchors_middle_res}[$f];
      my $a2_mid = ${$a2}[0];
      my @a2_eles = @{${$a2}[1]};
      my $a2_strand = $a2_eles[4];
      my $a2_region = $a2_eles[5];
      my $a2_start = $a2_eles[7];
      my $a2_end = $a2_eles[8];
      my $a2_failed_coverage = ${$a2}[2];
      my $a2_failed_identity = ${$a2}[3];

      if($a2_failed_coverage) {
        $issue_tracker++;
      }

      if($a2_failed_identity) {
        $issue_tracker++;
      }

      if($a2_strand ne $a1_strand) {
        $issue_tracker++;
      }

      if($a2_region ne $a1_region) {
        $issue_tracker++;
      }

      my $dist1 = abs($a1_mid - $a2_mid);
      for(my $g=0; $g<scalar(@$anchors_right_res); $g++) {
        my $a3 = ${$anchors_right_res}[$g];
        my $a3_mid = ${$a3}[0];
        my @a3_eles = @{${$a3}[1]};
        my $a3_strand = $a3_eles[4];
        my $a3_region = $a3_eles[5];
        my $a3_start = $a3_eles[7];
        my $a3_end = $a3_eles[8];
        my $a3_failed_coverage = ${$a3}[2];
        my $a3_failed_identity = ${$a3}[3];

        if($a3_failed_coverage) {
          $issue_tracker++;
        }

        if($a3_failed_identity) {
          $issue_tracker++;
        }

        if($a3_strand ne $a2_strand) {
          $issue_tracker++;
        }

        if($a3_region ne $a2_region) {
          $issue_tracker++;
        }

        my $dist2 = abs($a2_mid - $a3_mid);
        my $dist3 = abs($a1_mid - $a3_mid);
        my $final_dist = abs($source_dist1 - $dist1) + abs($source_dist2 - $dist2) + abs($source_dist3 - $dist3);
        my $co_linear = 0;

        # There is no perfect way of doing this, but the idea here is if a path isn't on the same region and it's not
        # co-linear then the path should be skipped. There's an argument for adding strand to this, but small local inversions
        # could happen and if things are still co-linear then it's still work considering the path
        # We can be stricter here and consider this the gold standard set, whereas a pairwise comparison after this could
        # pick up on good partial anchor sets
        unless($a1_region eq $a2_region and $a2_region eq $a3_region) {
          next;
        }


        if($final_dist > $distance_limit) {
          next;
        }

        if($a1_strand eq '+' and $a2_mid >= $a1_mid and $a3_mid >= $a2_mid) {
          $co_linear = 1;
        } elsif($a1_strand eq '-' and $a2_mid <= $a1_mid and $a3_mid <= $a2_mid) {
          $co_linear = 1;
        } else {
          next;
        }


        say "  Current anchors: (".$a1_eles[7].":".$a1_eles[8].") (".$a2_eles[7].":".$a2_eles[8].") (".$a3_eles[7].":".$a3_eles[8]."). Dist: ".$final_dist.
            ", Co-linear: ".$co_linear.", Issue count: ".$issue_tracker;

        my $anchor_set_start = min(($a1_start,$a2_start,$a3_start));
        my $anchor_set_end = max(($a1_end,$a2_end,$a3_end));
        my $combined_coord_string = $anchor_set_start.":".$anchor_set_end;
        unless($existing_anchor_coords->{$combined_coord_string}) {
          $existing_anchor_coords->{$combined_coord_string} = [$a1,$a2,$a3,$final_dist,$issue_tracker];
        } else {
          my $existing_anchor = $existing_anchor_coords->{$combined_coord_string};
          my $existing_dist = ${$existing_anchor}[3];
          my $existing_issue_count = ${$existing_anchor}[4];
          if($issue_tracker < $existing_issue_count) {
            $existing_anchor_coords->{$combined_coord_string} = [$a1,$a2,$a3,$final_dist,$issue_tracker];
          } if($issue_tracker == $existing_issue_count and $final_dist < $existing_dist) {
            $existing_anchor_coords->{$combined_coord_string} = [$a1,$a2,$a3,$final_dist,$issue_tracker];
          }
        }
        $issue_tracker = 0;
      } # end for(my $g
    } # end for(my $f
  } # end for(my $e

  my $issue_limit = 2;
  my $current_min_issue_count;
  foreach my $anchor_coord (keys(%$existing_anchor_coords)) {
    my $anchor_set = $existing_anchor_coords->{$anchor_coord};
    my $issue_count = ${$anchor_set}[4];
    unless(defined($current_min_issue_count)) {
      $current_min_issue_count = $issue_count;
    }

    if($issue_count > $issue_limit or $issue_count > $current_min_issue_count) {
      next;
    }

    if($issue_count < $current_min_issue_count) {
      $anchor_sets = [$anchor_set];
      $current_min_issue_count = $issue_count;
    } elsif($issue_count == $current_min_issue_count) {
      push(@$anchor_sets,$anchor_set);
    }
  }
  return($anchor_sets);
}


sub get_pairwise_sets {
  my ($self,$source_anchor1_midpoint,$source_anchor2_midpoint,$source_anchor3_midpoint,$anchors_left_res,$anchors_middle_res,$anchors_right_res) = @_;

  my $source_dist1 = abs($source_anchor2_midpoint - $source_anchor1_midpoint);
  my $source_dist2 = abs($source_anchor3_midpoint - $source_anchor2_midpoint);
  my $source_dist3 = abs($source_anchor3_midpoint - $source_anchor1_midpoint);

  # This is a bit like the get_best_anchor sets, but it's stricter since it's only looking at two anchors
  my $pairwise_anchors = [];
  $self->filter_pairwise_anchors($anchors_left_res,$anchors_middle_res,$source_dist1,$pairwise_anchors);
  $self->filter_pairwise_anchors($anchors_middle_res,$anchors_right_res,$source_dist2,$pairwise_anchors);
  $self->filter_pairwise_anchors($anchors_left_res,$anchors_right_res,$source_dist3,$pairwise_anchors);

  return($pairwise_anchors);
}


sub filter_pairwise_anchors {
  my ($self,$anchors1,$anchors2,$source_dist,$pairwise_anchors) = @_;
  my $existing_anchor_coords = {};
  for(my $e=0; $e<scalar(@$anchors1); $e++) {
    my $a1 = ${$anchors1}[$e];
    my $a1_mid = ${$a1}[0];
    my @a1_eles = @{${$a1}[1]};
    my $a1_strand = $a1_eles[4];
    my $a1_region = $a1_eles[5];
    my $a1_start = $a1_eles[7];
    my $a1_end = $a1_eles[8];
    my $a1_failed_coverage = ${$a1}[2];
    my $a1_failed_identity = ${$a1}[3];
    if($a1_failed_coverage or $a1_failed_identity) {
      next;
    }

    for(my $f=0; $f<scalar(@$anchors2); $f++) {
      my $a2 = ${$anchors2}[$f];
      my $a2_mid = ${$a2}[0];
      my @a2_eles = @{${$a2}[1]};
      my $a2_strand = $a2_eles[4];
      my $a2_region = $a2_eles[5];
      my $a2_start = $a2_eles[7];
      my $a2_end = $a2_eles[8];
      my $a2_failed_coverage = ${$a2}[2];
      my $a2_failed_identity = ${$a2}[3];

      if($a2_failed_coverage or $a2_failed_identity) {
        next;
      }

      if($a2_strand ne $a1_strand or $a2_region ne $a1_region) {
        next;
      }

      my $dist = abs($a1_mid - $a2_mid);
      unless($dist <= $source_dist * 1.2) {
        next;
      }

      my $anchor_set_start = min(($a1_start,$a2_start));
      my $anchor_set_end = max(($a1_end,$a2_end));
      my $combined_coord_string = $anchor_set_start.":".$anchor_set_end;
      unless($existing_anchor_coords->{$combined_coord_string}) {
        $existing_anchor_coords->{$combined_coord_string} = [$a1,$a2,$dist];
      } else {
        my $existing_anchor = $existing_anchor_coords->{$combined_coord_string};
        my $existing_dist = ${$existing_anchor}[2];
        if($dist < $existing_dist) {
          $existing_anchor_coords->{$combined_coord_string} = [$a1,$a2,$dist];
        }
      }

    } # end for(my $f
  } # end for(my $e

  foreach my $anchor_coord (keys(%$existing_anchor_coords)) {
    my $anchor_set = $existing_anchor_coords->{$anchor_coord};
    push(@$pairwise_anchors,$anchor_set);
  }
}


sub calculate_target_midpoint {
  my ($self,$source_hit_start,$source_hit_end,$source_hit_midpoint,$target_genomic_start,$target_genomic_end,$target_strand,$target_genomic_length) = @_;

  # In this case the midpoint is within the hit boundaries, so no adjustment is needed. Note that the midpoint is only quite approximate as there
  # could be changes in length between the sounce and target
  my $target_midpoint;
  if($source_hit_midpoint >= $source_hit_start and $source_hit_midpoint <= $source_hit_end) {
    if($target_strand eq '+') {
      my $source_hit_offset = $source_hit_midpoint - $source_hit_start;
      # Could put in a target offset here in terms of the difference in the length of the source and target sets
      $target_midpoint = $target_genomic_start + $source_hit_offset;
      return($target_midpoint);
    } else {
      my $source_hit_offset = $source_hit_end - $source_hit_midpoint;
      $target_midpoint = $target_genomic_start + $source_hit_offset;
    }
  } else{
    # This means the source midpoint is outside of the hit, so we want to guess where it is, though the guess is going to be wrong since the implication is
    # that the reason the hit is incomplete is because there's some sort of gap
    if($target_strand eq '+') {
      if($source_hit_midpoint < $source_hit_start) {
        my $source_hit_offset = $source_hit_start - $source_hit_midpoint;
        $target_midpoint = $target_genomic_start - $source_hit_offset;
        if($target_midpoint < 1) {
          $target_midpoint = 1;
        }
        return($target_midpoint);
      } else {
        my $source_hit_offset = $source_hit_midpoint - $source_hit_end;
        $target_midpoint = $target_genomic_end + $source_hit_offset;
        if($target_midpoint > $target_genomic_length) {
          $target_midpoint = $target_genomic_length;
        }
        return($target_midpoint);
      }

    } else {
      if($source_hit_midpoint < $source_hit_start) {
        my $source_hit_offset = $source_hit_start - $source_hit_midpoint;
        $target_midpoint = $target_genomic_end + $source_hit_offset;
        if($target_midpoint > $target_genomic_length) {
          $target_midpoint = $target_genomic_length;
        }
        return($target_midpoint);
      } else {
        my $source_hit_offset = $source_hit_midpoint - $source_hit_end;
        $target_midpoint = $target_genomic_start - $source_hit_offset;
        if($target_midpoint < 1) {
          $target_midpoint = 1;
        }
        return($target_midpoint);
      }
    }
  } # End outer else
}





sub project_batch_genes {
  my ($self,$batched_input_genes,$source_sequence_adaptor) = @_;

  my $target_genes_by_id = {};
  my $target_adaptor = $self->target_adaptor();
  my $target_sequence_adaptor = $target_adaptor->get_SequenceAdaptor();
  my $target_slice_adaptor = $target_adaptor->get_SliceAdaptor();
  foreach my $id (keys(%$batched_input_genes)) {
    say "Projecting batch ID: ".$id;

    my $batch = $batched_input_genes->{$id};
    my $target_regions = $batch->{'target_regions'};
    unless($target_regions and scalar(@$target_regions)) {
      say "No target regions found for batch ".$id.", skipping";
      next;
    }

    my $source_genomic_start = $batch->{'start'};
    my $source_genomic_end = $batch->{'end'};
    my $source_parent_slice = $batch->{'slice'};
    my $source_region_name = $batch->{'slice_name'};
    my $source_region_strand = 1;
    say "  Source region: (".$source_genomic_start."/".$source_genomic_end.") ".$source_region_name." ".$source_region_strand;
    my $source_genomic_sequence = ${$source_sequence_adaptor->fetch_by_Slice_start_end_strand($source_parent_slice,$source_genomic_start,$source_genomic_end,$source_region_strand)};
    foreach my $target_region (@$target_regions) {
      my $target_region_name = ${$target_region}[0];
      my $target_region_strand = ${$target_region}[1];
      my $target_genomic_start = ${$target_region}[2];
      my $target_genomic_end = ${$target_region}[3];
      if($target_region_strand eq '+') {
        $target_region_strand = 1;
      } else {
        $target_region_strand = -1;
      }
      ${$target_region}[1] = $target_region_strand;

      say "  Target region: (".$target_genomic_start."/".$target_genomic_end.") ".$target_region_name." ".$target_region_strand;
      # At this point we want to get the source and target sequence, the target sequence should be stranded on the region

      my $target_parent_slice = $target_slice_adaptor->fetch_by_region('toplevel',$target_region_name);
      my $target_genomic_sequence = ${$target_sequence_adaptor->fetch_by_Slice_start_end_strand($target_parent_slice, $target_genomic_start, $target_genomic_end, $target_region_strand)};

      my $target_region_slice = $target_slice_adaptor->fetch_by_region('toplevel',$target_region_name,$target_genomic_start,$target_genomic_end,1);

      my $coverage = 0;
      my $percent_id = 0;
      my $aligned_source_seq;
      my $aligned_target_seq;

      eval {
        ($coverage,$percent_id,$aligned_source_seq,$aligned_target_seq) = align_nucleotide_seqs($source_genomic_sequence,$target_genomic_sequence);
      };
      if ($@) {
        $self->warning("Issue with running MAFFT on target region");
        next;
      } else {
        say "Aligned source and target regions with MAFFT: Coverage: ".$coverage.", Percent id: ".$percent_id;
      }

      ${$target_region}[4] = $aligned_source_seq;
      ${$target_region}[5] = $aligned_target_seq;
      ${$target_region}[6] = $target_region_slice;
    }
    $self->build_batch_genes($batch,$target_genes_by_id);
  } # foreach my $id
  return($target_genes_by_id);
}


sub build_batch_genes {
  my ($self,$batch,$target_genes_by_id) = @_;

  my $batch_source_genes = $batch->{'genes'};
  unless(scalar(@$batch_source_genes)) {
    $self->throw("No source genes found for the batch");
  }

  my $source_genomic_start = $batch->{'start'};
  my $source_genomic_end = $batch->{'end'};
  my $source_parent_slice = $batch->{'slice'};
  my $source_region_name = $batch->{'slice_name'};
  my $source_region_strand = 1;

  my $target_regions = $batch->{'target_regions'};
  foreach my $target_region (@$target_regions) {
    my $target_transcripts_by_id = {};
    my $target_region_name = ${$target_region}[0];
    my $target_region_strand = ${$target_region}[1];
    my $target_genomic_start = ${$target_region}[2];
    my $target_genomic_end= ${$target_region}[3];
    my $aligned_source_seq = ${$target_region}[4];
    my $aligned_target_seq = ${$target_region}[5];
    my $target_region_slice = ${$target_region}[6];
    foreach my $source_gene (@$batch_source_genes) {
      my $source_transcripts = $source_gene->get_all_Transcripts();
      my $projected_exons_by_id = {};
      my $exons = $source_gene->get_all_Exons();
      foreach my $exon (@$exons) {
        my $exon_genomic_start = $exon->seq_region_start();
        my $exon_genomic_end = $exon->seq_region_end();
        say "Source region start/end: ".$source_genomic_start."/".$source_genomic_end;
        say "Source exon region start/end: ".$exon_genomic_start."/".$exon_genomic_end;
        say "Source exon strand: ".$exon->strand();
        say "Target strand: ".$target_region_strand;
        my $projected_exon = $self->project_batch_feature(undef,$exon,$source_genomic_start,$exon_genomic_start,$exon_genomic_end,$aligned_source_seq,$aligned_target_seq,
                                                          $target_region_slice,$target_region_strand);

        if($projected_exon) {
          if ($exon->seq->seq eq $projected_exon->seq->seq) {
            $projected_exon->{'cov'} = 100.00;
            $projected_exon->{'perc_id'} = 100.00;
          }
          else {
            my ($proj_coverage,$proj_percent_id,$aligned_source_seq,$aligned_target_seq) = align_nucleotide_seqs($exon->seq->seq(),$projected_exon->seq->seq());
            $projected_exon->{'cov'} = $proj_coverage;
            $projected_exon->{'perc_id'} = $proj_percent_id;
          }
          $projected_exon->{'source_stable_id'} = $exon->stable_id();
          $projected_exon->{'source_length'} = $exon->length();
#        say "Projected exon (".$projected_exon->{'source_stable_id'}."): Coverage: ".$projected_exon->{'cov'}.", Percent id: ".$projected_exon->{'perc_id'};
#        say "Alignment:\n".$aligned_source_seq."\n".$aligned_target_seq;
          $projected_exons_by_id->{$projected_exon->{'source_stable_id'}} = $projected_exon;
        } else {
          say "Failed to project exon (".$exon->stable_id.")";
        }
      } # end foreach my $exon


      my $target_gene = Bio::EnsEMBL::Gene->new(-analysis => $self->analysis);
    print join(' ', __LINE__, $target_gene), "\n";
      $target_gene->{'parent_gene_id'} = $source_gene->dbID();
      $target_gene->{'parent_gene_versioned_stable_id'} = $source_gene->stable_id().".".$source_gene->version();
      $target_gene->stable_id($source_gene->stable_id);
      $target_gene->version($source_gene->version);
      $target_gene->biotype($source_gene->biotype);
      $target_gene->description($source_gene->description());
      #my $gene_description = "Parent: ".$target_gene->stable_id().".".$target_gene->version().", Type: Primary mapping";
      #$target_gene->description($gene_description);

      # add source gene stable id as gene attribute
      my $parent_attribute = Bio::EnsEMBL::Attribute->new(-CODE => 'proj_parent_g',-VALUE => $source_gene->stable_id_version);
      $target_gene->add_Attributes($parent_attribute);

      my $complete_projections = 0;
      foreach my $source_transcript (@$source_transcripts) {
        my $projected_transcript = $self->reconstruct_transcript($source_transcript,$projected_exons_by_id);
        if($projected_transcript) {
          $projected_transcript->{'parent_transcript_versioned_stable_id'} = $source_transcript->stable_id().".".$source_transcript->version();
          $projected_transcript->{'annotation_method'} = 'alignment_projection';
          if($source_transcript->translation()) {
            $self->project_cds($projected_transcript,$source_transcript);
            say join(' ', __LINE__, $projected_transcript);
            $projected_transcript = $self->qc_cds_sequence($projected_transcript,$source_transcript);
            say join(' ', __LINE__, $projected_transcript);
          }
          $self->set_transcript_description($projected_transcript,$source_transcript);
          
          if ($projected_transcript->cdna_coding_start() and
              $projected_transcript->cdna_coding_end() and
              $projected_transcript->cdna_coding_end()-$projected_transcript->cdna_coding_start()+1 < 3) {
            # By convention, the coding_region_end is always higher than the
            # value returned by the coding_region_start method.
            say "Projected transcript CDS is too short (< 3 bp). Parent transcript stable id: ".$projected_transcript->{'parent_transcript_versioned_stable_id'};
          } elsif ($projected_transcript->biotype() eq 'protein_coding' and !($projected_transcript->translation())) {
            say "Projected transcript biotype is protein_coding but it does not have any translation. Parent transcript stable id: ".$projected_transcript->{'parent_transcript_versioned_stable_id'};
          } else {
            $target_gene->add_Transcript($projected_transcript);
            if($projected_transcript->{'cov'} >= 99 and $projected_transcript->{'perc_id'} >= 99) {
              $complete_projections++;
            }
          }
        }
      } # end foreach my $source_transcript

      unless($target_gene) {
        say "Unable to project source gene with stable id: ".$source_gene->stable_id();
        next;
      }
      if(scalar(@$source_transcripts) == $complete_projections) {
        say "Complete projection for gene with stable_id ".$source_gene->stable_id();
      }
      if($target_gene->get_all_Transcripts) {
        if($target_genes_by_id->{$source_gene->dbID()}) {
          push(@{$target_genes_by_id->{$source_gene->dbID()}},$target_gene);
        } else {
          $target_genes_by_id->{$source_gene->dbID()} = [$target_gene];
        }
      } else {
        $self->warning("Target gene has no transcripts, so will not keep. Source gene stable id: ".$source_gene->stable_id());
      }
    } # end foreach my $source_gene
  } # end foreach my $target_region
}


sub project_batch_feature {
  my ($self,$transcript,$exon,$source_region_start,$feature_start,$feature_end,$aligned_source_seq,$aligned_target_seq,$target_region_slice,$target_strand) = @_;

  # This is just a copy of project_feature to test on batches, as the batch projection is fundamentally different to regular projection
  # The regunlar projection was done with the source region being stranded on the gene strand, with the batches all the source regions
  # are on the forward strand, and the genes themselves can vary across strands. This is then further complicated by the target region
  # also potentially being on the opposite strand
  my $source_seq = $aligned_source_seq;
  $source_seq =~ s/\-//g;
  my $target_seq = $aligned_target_seq;
  $target_seq =~ s/\-//g;

  my @source_align_array = split('',$aligned_source_seq);
  my @target_align_array = split('',$aligned_target_seq);
  my $source_seq_pos = 0;
  my $source_align_seq_pos = 0;
  my $target_seq_pos = 0;
  my $target_align_seq_pos = 0;

  # Seen to have some issue in terms of getting the real start/end for some features. Was using +1, but seems off
  my $source_feature_seq_start =  $feature_start - $source_region_start + 1;
  my $source_feature_seq_end = $feature_end - $source_region_start + 1;
  say "Source feature seq start: ".$source_feature_seq_start;
  say "Source feature seq end: ".$source_feature_seq_end;
  say "Source region start: ".$source_region_start;
  say "Feature start: ".$feature_start;
  say "Feature end: ".$feature_end;

  my $source_seq_start_index;
  my $source_seq_end_index;
  my $target_seq_start_index;
  my $target_seq_end_index;
  my $in_feature_alignment = 0;
  my $out_of_feature_alignment = 0;
  my $align_feature_start;
  my $align_feature_end;
  for(my $i=0; $i < scalar(@source_align_array); $i++) {
    my $source_value = $source_align_array[$i];
    my $target_value = $target_align_array[$i];
    if($source_value ne '-') {
      $source_seq_pos++;
    }

    if($target_value ne '-') {
      $target_seq_pos++;
    }

    if($source_seq_pos == $source_feature_seq_start) {
      $source_seq_start_index = $source_seq_pos;
      $align_feature_start = $i;
      say "Index of the feature start in source seq: ".$source_seq_start_index;
      $in_feature_alignment = 1;
    }

    if($in_feature_alignment and $target_value ne '-' and !(defined($target_seq_start_index))) {
      $target_seq_start_index = $target_seq_pos;
    }

    if($source_seq_pos == $source_feature_seq_end) {
      $source_seq_end_index = $source_seq_pos;
      $target_seq_end_index = $target_seq_pos;
      $align_feature_end = $i;
      say "Index of the feature end in source seq: ".$source_seq_end_index;
      $in_feature_alignment = 0;
    }
  }

  unless(defined($target_seq_start_index) and defined($target_seq_end_index)) {
    $self->warning("Issue with recovering start/end of the exon feature in target alignment sequence, not building exon");
    return;
  }

  my $source_feature_length = $source_seq_end_index - $source_seq_start_index + 1;
  my $target_feature_length = $target_seq_end_index - $target_seq_start_index + 1;
  if ($source_feature_length*1.5 < $target_feature_length) {
    $self->warning("Target exon is at least 1.5 bigger than source: $source_feature_length < $target_feature_length");
    return;
  }
  my $recovered_source_feature_seq = substr($source_seq,$source_seq_start_index-1,$source_feature_length);
  if($exon->strand != 1) {
    $recovered_source_feature_seq = $self->revcomp($recovered_source_feature_seq);
  }


  my $align_feature_length = $align_feature_end - $align_feature_start + 1;
  my $source_feature_align_seq = substr($aligned_source_seq,$align_feature_start,$align_feature_length);
  my $target_feature_align_seq = substr($aligned_target_seq,$align_feature_start,$align_feature_length);
  say "Source feature in alignment:\n".$source_feature_align_seq;
  say "Target feature in alignment:\n".$target_feature_align_seq;

  my $recovered_target_feature_seq = substr($target_seq,$target_seq_start_index-1,$target_feature_length);
  # There's something wrong with this, but I'm not sure what the right condition is. The feature seq is correct in the end, but
  # I think this is wrong if the target strand is -1
  if($target_strand * $exon->strand != 1) {
    $recovered_target_feature_seq = $self->revcomp($recovered_target_feature_seq);
  }

  say "Source seq start index/end index/length: ".$source_seq_start_index.":".$source_seq_end_index.":".$source_feature_length;
  say "Target seq start index/end index/length: ".$target_seq_start_index.":".$target_seq_end_index.":".$target_feature_length;
  say "Recovered feature source seq:";
  say $recovered_source_feature_seq;
  say "Recovered feature target seq:";
  say $recovered_target_feature_seq;

  my $projected_exon = $self->build_projected_exon($transcript,$exon,$target_seq_start_index,$target_seq_end_index,$target_region_slice,$target_strand);
  return($projected_exon);
}


sub coverage_cutoff {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_coverage_cutoff} = $val;
  }

  return $self->{_coverage_cutoff};
}


sub perc_id_cutoff {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_perc_id_cutoff} = $val;
  }

  return $self->{_perc_id_cutoff};
}


sub extended_length_variation_cutoff {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_extended_length_variation_cutoff} = $val;
  }

  return $self->{_extended_length_variation_cutoff};
}


sub anchor_coverage_cutoff {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_anchor_coverage_cutoff} = $val;
  }

  return $self->{_anchor_coverage_cutoff};
}


sub anchor_perc_id_cutoff {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_anchor_perc_id_cutoff} = $val;
  }

  return $self->{_anchor_perc_id_cutoff};
}


1;
