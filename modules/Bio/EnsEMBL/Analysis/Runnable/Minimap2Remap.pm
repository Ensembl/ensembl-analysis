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

Bio::EnsEMBL::Analysis::Runnable::Star

=head1 SYNOPSIS

  my $runnable =
    Bio::EnsEMBL::Analysis::Runnable::Star->new();

 $runnable->run;
 my @results = $runnable->output;

=head1 DESCRIPTION

This module uses Star to align fastq to a genomic sequence. Star is a splice aware
aligner. It creates output files with the reads overlapping splice sites and the reads
aligning on the exons. Some reads are aligned multiple times in the genome.

=head1 METHODS

=cut


package Bio::EnsEMBL::Analysis::Runnable::Minimap2Remap;

use warnings;
use strict;
use feature 'say';
use Data::Dumper;
use POSIX;

use Bio::EnsEMBL::Analysis::Runnable::Minimap2;
use Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(compute_best_translation);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(set_alignment_supporting_features attach_Slice_to_Transcript attach_Analysis_to_Transcript calculate_exon_phases replace_stops_with_introns);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name align_nucleotide_seqs map_cds_location align_proteins execute_with_wait);
use File::Spec;
use Bio::DB::HTS::Faidx;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use parent ('Bio::EnsEMBL::Analysis::Runnable');


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
  my ($genome_index, $input_file, $paftools_path, $source_adaptor, $target_adaptor, $delete_input_file, $parent_genes, $parent_gene_ids, $gene_synteny_hash, $gene_genomic_seqs_hash, $no_projection) = rearrange([qw (GENOME_INDEX INPUT_FILE PAFTOOLS_PATH SOURCE_ADAPTOR TARGET_ADAPTOR DELETE_INPUT_FILE PARENT_GENES PARENT_GENE_IDS GENE_SYNTENY_HASH GENE_GENOMIC_SEQS_HASH NO_PROJECTION)],@args);
  $self->genome_index($genome_index);
  $self->input_file($input_file);
  $self->paftools_path($paftools_path);
  $self->source_adaptor($source_adaptor);
  $self->target_adaptor($target_adaptor);
  $self->delete_input_file($delete_input_file);
  $self->genes_to_process($parent_genes);
  $self->parent_gene_ids($parent_gene_ids);
  $self->gene_synteny_hash($gene_synteny_hash);
  $self->gene_genomic_seqs_hash($gene_genomic_seqs_hash);
  $self->no_projection($no_projection);
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

  # run minimap2
  my $minimap2_command = $self->program." --cs --secondary=no -x map-ont ".$genome_index." ".$input_file." > ".$paf_file;
  $self->warning("Command:\n".$minimap2_command."\n");
  execute_with_wait($minimap2_command);

  my $paf_results = [];
  open(IN, $paf_file) or $self->throw("Could not open $paf_file for reading");
  while(<IN>) {
    chomp($_);
    push(@$paf_results,$_);
  }
  close(IN) or $self->throw("Could not close $paf_file");

  my $high_confidence = 0;
  my $total_results = 0;
  my $processed_gene_ids = {};
  my $paf_results_hash = {};
  foreach my $paf_result (@$paf_results) {
    say "PAF results for first pass:\n".$paf_result;
    my @result_cols = split("\t",$paf_result);
    my $gene_id = $result_cols[0];
    if($paf_results_hash->{$gene_id}) {
      push(@{$paf_results_hash->{$gene_id}},\@result_cols)
    } else {
      $paf_results_hash->{$gene_id} = [\@result_cols];
    }
  }

  my $genes_to_process = $self->genes_to_process();
  foreach my $gene (@$genes_to_process) {
    say "Processing source gene: ".$gene->stable_id();
    $high_confidence += $self->process_results($gene,$paf_results_hash);
    $total_results++;
  }

  say "TOTAL RESULTS: ".$total_results;
  say "HIGH CONFIDENCE: ".$high_confidence;
} # End run


sub process_results {
  my ($self,$source_gene,$paf_results) = @_;

  my $high_confidence = 0;

  my $gene_seq = $source_gene->seq();

  my $source_transcripts = $source_gene->get_all_Transcripts();
  my $source_transcript_id_hash = {}; # A hash to organise transcripts by dbID (note that the dbID is somewhat confusingly saved in the stable id field here, safer than using the realstable id)

  say "Source transcript list:";
  foreach my $source_transcript (@$source_transcripts) {
    say "  ".$source_transcript->stable_id();
    my $source_transcript_id = $source_transcript->dbID();
    $source_transcript_id_hash->{$source_transcript_id} = $source_transcript;
  }

  my $max_intron_size = $self->calculate_max_intron_size($source_transcripts);
  $max_intron_size = ceil($max_intron_size);
  say "Max intron size of transcripts in the cluster: ".$max_intron_size;

  my $gene_paf_results = $paf_results->{$source_gene->dbID()};
  my $good_transcripts = []; # Transcripts that pass the cut-off thresholds
  my $bad_source_transcripts = []; # The source transcripts for mapped Transcripts that don't pass the threshold
  my $bad_transcripts = []; # When both the minimap and exonerate mappings fail the cut-offs, this will store the version with the highest combined coverage and percent id

  my $target_adaptor = $self->target_adaptor();

  my $final_genes = [];

  # This is a paf result from the alignment of a genomic read representing a gene against the target genome
  # Only one such mapped read per gene will be passed into this based on the top hit in the paf results
  # This method will update the mapped area if it is shorter than the genomic span of the original gene
  # Once the updated mapped area is calculated:
  # 1) Run minimap2 on the transcripts from the source gene onto the mapped region
  #    Foreach source transcript:
  #      - Run mapping
  #      - If coding, calculate coverage and identity of the CDS
  #      - If coding, check whether CDS is longer than expected
  #      - If non-coding, calculated coverage and identity of full transcript seq
  # 2) Once transcripts have been mapped, assess quality of mapping
  #    Foreach minimap2 transcript:
  #      - If the coverage or percent identity of the seq is < 95, or if there's a CDS sequence and it is > 5 percent longer
  #        -> Run exonerate the soruce transcript on the region
  #        -> If the cov/perc id of the resulting transcript passes the above criteria mark as good
  #      - Else mark the minimap2 transcript as good
  # 3) For the transcripts not marked as good, examine the minimap2 and exonerate versions:
  #    Foreach bad transcript:
  #      - Check the other transcripts, if available
  #      - If there is another transcript from the gene marked as good
  #        -> Examine the span of all good transcripts in the source gene
  #        -> If the bad transcript is within the boundary of that span in the source gene then there is an issue, but just add or filter out
  #        -> If the bad transcript extends outside the boundary, calculate the offset and try and align to the offset equivalent in the target, then add or filter out
  #      - If there are no other transcripts, then the gene is bad and this needs to be dealt with separately
  # 4) Deal with bad genes (gene with no good transcripts)
  #    Foreach bad gene:
  #      - Take the source canonical transcript and align to the whole genome using minimap2 and exonerate
  #      - For the best of the two alignment, check if they are better than the original combined cov/identity (if there is an existing model)
  #      - If better or no existing then:
  #        -> Calculate offsets based on the source canonical transcript for the expected gene boundaries in the target
  #        -> Using the calculated boundaries, get the region in target and align remaining transcripts
  #      - If worse, then do no more
  #      - Once you have the final set of transcripts decide whether to keep the gene or not

  my $good_transcripts_hash = {};
  my $bad_transcripts_hash = {};
  my $best_transcripts_by_id = {};
  if($gene_paf_results) {
    my $chained_paf_result = $self->calculate_region_boundary($gene_paf_results);

    my $target_genomic_start = ${$chained_paf_result}[0];
    my $target_genomic_end = ${$chained_paf_result}[1];
    my $target_strand = ${$chained_paf_result}[2];
    my $target_genomic_name = ${$chained_paf_result}[3];

    if($target_strand eq '+') {
      $target_strand = 1;
    } else {
      $target_strand = -1;
    }

    say "First pass genomic start/end: ".$target_genomic_start."/".$target_genomic_end;


    say "First pass adjusted genomic start/end: ".$target_genomic_start."/".$target_genomic_end;

    my $target_slice_adaptor = $target_adaptor->get_SliceAdaptor();
    my $target_parent_slice = $target_slice_adaptor->fetch_by_region('toplevel',$target_genomic_name);

    unless($target_parent_slice) {
      $self->throw("Could not fetch the slice in the target assembly. Slice name: ".$target_genomic_name);
    }

    # Note that we are going to use the forward strand regardless of what strand the PAF hit is on here because minimap2/exonerate assume the region is on the forward strand
    my $target_sequence_adaptor = $target_adaptor->get_SequenceAdaptor;
    my $target_genomic_seq = ${ $target_sequence_adaptor->fetch_by_Slice_start_end_strand($target_parent_slice, $target_genomic_start, $target_genomic_end, 1) };
    my $target_region_slice = $target_slice_adaptor->fetch_by_region('toplevel', $target_genomic_name, $target_genomic_start, $target_genomic_end, 1);
    say "Target slice identified to search for transcripts in after first pass: ".$target_region_slice->name();
    my $target_genomic_fasta = ">".$target_genomic_name."\n".$target_genomic_seq;
    my $target_genome_file = $self->write_input_file([$target_genomic_fasta]);

    say "Projecting gene: ".$source_gene->stable_id();

    my $coverage_threshold = 98;
    my $perc_id_threshold = 98;
    my $projected_transcripts_by_id = $self->project_gene_coords($source_gene,$source_transcripts,$target_genomic_start,$target_region_slice,$target_strand);
    $self->print_transcript_stats($projected_transcripts_by_id,'projection');
    $self->update_best_transcripts($best_transcripts_by_id,$projected_transcripts_by_id);
    my $minimap_transcripts_by_id = $self->map_gene_minimap($source_gene,$source_transcripts,$target_genomic_start,$target_region_slice,$target_strand,
                                                            $target_genome_file,$source_transcript_id_hash,$max_intron_size,$target_adaptor,$target_slice_adaptor,
                                                            $best_transcripts_by_id);
    $self->print_transcript_stats($minimap_transcripts_by_id,'minimap local');

    $self->update_best_transcripts($best_transcripts_by_id,$minimap_transcripts_by_id);
    my $exonerate_transcripts_by_id = $self->map_gene_exonerate($source_transcripts,$target_genomic_start,$target_region_slice,$target_strand,
                                                                $target_genome_file,$source_transcript_id_hash,$max_intron_size,$target_adaptor,
                                                                $target_slice_adaptor,$best_transcripts_by_id);
    $self->print_transcript_stats($exonerate_transcripts_by_id,'exonerate local');

    $self->update_best_transcripts($best_transcripts_by_id,$exonerate_transcripts_by_id);
  } else {
    say "Did not get a paf result from the first pass, will resort to global mapping of the transcripts";
  }


  # Now split into good and bad transcripts
# $self->split_initial_transcripts($best_transcripts_by_id,$good_transcripts,$bad_transcripts);

  # At this point we have the selected set of good/bad transcripts. There may also be missing transcripts at the moment, i.e.
  # those that did no get aligned by either method. We now want to take the original sequences for the bad and missing transcripts
  # and run a genome-wide alignment via minimap. We then have two scenarios to work out. The first is if there were any good transcripts
  # to begin with. If there were then were going to make the assumption (correctly or incorrectly), that the good transcripts represent
  # where the gene is likely to be. We need to decide once we do the global mapping what to do with any mapped transcripts. The first thing
  # is to check whether a global transcript actually scores better than the original mapping (if there was one). If not, then we just add
  # the selected bad transcript to the gene. If it is better then we have to work out what to do. We basically need to cluster the transcripts
  # and figure out if they're all in the same region.

  say "Checking for missing transcripts";
  my $transcripts_for_global_mapping = $self->filter_transcripts_to_map($source_transcripts,$best_transcripts_by_id);
  say "Found ".scalar(@$transcripts_for_global_mapping)." missing transcripts";

  # Now run a global mapping of the missing/bad transcripts
  my $transcripts_to_map_fasta_seqs = [];
  foreach my $transcript (@$transcripts_for_global_mapping) {
    my $transcript_sequence = $transcript->seq->seq();
    my $fasta_record = ">".$transcript->dbID()."\n".$transcript_sequence;
    push(@$transcripts_to_map_fasta_seqs,$fasta_record);
  }

  say "Preparing to map bad and missing transcripts globally to the genome";
  say "Number of transcripts to map globally: ".scalar(@$transcripts_for_global_mapping);
  my $target_global_genome_index = $self->genome_index();
  my $global_mapped_transcripts = $self->generate_minimap_transcripts($transcripts_to_map_fasta_seqs,$target_global_genome_index,$target_adaptor,1,$max_intron_size);
  say "Number of globally mapped transcripts: ".scalar(@$global_mapped_transcripts);

  my $global_transcripts_by_id = {};
  foreach my $transcript (@$global_mapped_transcripts) {
    my $source_transcript = $source_transcript_id_hash->{$transcript->stable_id};
    my ($mapped_coverage,$mapped_percent_id,$aligned_source_seq,$aligned_target_seq) = align_nucleotide_seqs($source_transcript->seq->seq(),$transcript->seq->seq());
    $transcript->{'cov'} = $mapped_coverage;
    $transcript->{'perc_id'} = $mapped_percent_id;
    $transcript->{'aligned_source_seq'} = $aligned_source_seq;
    $transcript->{'aligned_target_seq'} = $aligned_target_seq;

    $transcript->{'annotation_method'} = 'minimap_global';
    my $db_id = $transcript->stable_id();
    $global_transcripts_by_id->{$db_id} = $transcript;
  }

  # Want to use global mapping if the the current best model is below the coverage thresholds and the global model is better
  # The coverage thresholds are lower in this instance as we would give more weight to any of the other approaches since they use the alignment
  $self->print_transcript_stats($global_transcripts_by_id,'minimap global');
  $self->update_best_transcripts($best_transcripts_by_id,$global_transcripts_by_id);

  $self->set_cds_sequences($best_transcripts_by_id,$source_transcript_id_hash);
  $self->qc_cds_sequences($best_transcripts_by_id,$source_transcript_id_hash);
  $self->fix_cds_issues($best_transcripts_by_id);
  $self->set_transcript_descriptions($best_transcripts_by_id,$source_transcript_id_hash);

  # You want ot cluster all transcripts at this point. Clusters of the original good transcripts determine where we believe the locus is
  # If you have such a cluster then you keep all good/bad transcripts in this cluster (taking care to remove any duplicates since some bad
  # transcripts will have global versions too)
  # If good cluster:
  #   Add bad transcripts
  #   Add good/bad global transcripts
  #   If the global is good, remove original bad version if present
  #   If the global is bad and there's a bad original, pick the best of the pair
  # If global good cluster:
  #   Add good/bad transcripts
  #   If the global is good, remove original bad version if present
  #   If the global is bad and there's a bad original, pick the best of the pair
  # If bad cluster and no good/global good cluster:
  #   Add clustered bad genes
  # If bad cluster and good/global good:
  #   Skip

  my $all_transcripts = $self->label_transcript_status($best_transcripts_by_id);
  my $biotypes_hash = $self->generate_biotypes_hash($all_transcripts);

  # Need to create single transcript genes for clustering purposes
  say "Building single transcript genes ahead of clustering";
  my $single_transcript_genes = $self->generate_single_transcript_genes($all_transcripts);

  my $genes_by_seq_region = {};
  foreach my $gene (@$single_transcript_genes) {
    my $seq_region_name = $gene->seq_region_name();
    unless($genes_by_seq_region->{$seq_region_name}) {
      $genes_by_seq_region->{$seq_region_name} = [];
    }
    push(@{$genes_by_seq_region->{$seq_region_name}},$gene);
  }

  say "Clustering genes";
  my $found_good_cluster = 0;
  my $all_clusters = [];
  foreach my $seq_region_name (keys(%{$genes_by_seq_region})) {
    my $genes = $genes_by_seq_region->{$seq_region_name};
    my ($clusters, $unclustered) = cluster_Genes($genes,$biotypes_hash);
    # There's an issue here in terms of genes
    say "Found ".scalar(@$clusters)." transcript clusters";
    say "Found ".scalar(@$unclustered)." unclustered transcripts";
    push(@$all_clusters,@$clusters);
    push(@$all_clusters,@$unclustered);
  }

  say "Found ".scalar(@$all_clusters)." initial clusters";

  # The all_clusters array now has every cluster associated with mapped transcripts from the current source gene
  # Some of these clusters may have multiple mappings for a source transcript, or even multiple copies across clusters
  # So there needs to be a way to take the best copies, ideally only having one mapped transcript per source transcript
  foreach my $cluster (@$all_clusters) {
    $self->check_cluster_status($cluster);
    if($cluster->{'status'} eq 'good') {
      $found_good_cluster = 1;
    }
  }

  # To get the final list of clusters, you take all good clusters. If there are no good clusters you take the bad clusters
  my $final_clusters = [];
  foreach my $cluster (@$all_clusters) {
    if($cluster->{'status'} eq 'good') {
      push(@$final_clusters,$cluster);
    } elsif($cluster->{'status'} eq 'bad' and !$found_good_cluster) {
      push(@$final_clusters,$cluster);
    }
  }

  say "Found ".scalar(@$final_clusters)." final clusters";

  say "Building genes from final clusters";
  # Now remove duplicates and form genes
  my $parent_gene_ids = $self->parent_gene_ids();

  foreach my $cluster (@$final_clusters) {
    my $gene = $self->create_gene_from_cluster($cluster,$parent_gene_ids,$source_transcript_id_hash);
    say "Created gene: ".$gene->stable_id()." ".$gene->seq_region_name().":".$gene->seq_region_start.":".$gene->seq_region_end.":".$gene->strand();
    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts) {
      say "  Transcript: ".$transcript->stable_id()." ".$transcript->seq_region_start.":".$transcript->seq_region_end.":".$transcript->strand();
      my $updated_description = $transcript->description().";annotation_method=".$transcript->{'annotation_method'};
      $transcript->description($updated_description);
    }

    push(@$final_genes,$gene);
  }

  say "Build ".scalar(@$final_genes)." final genes";

  $self->output($final_genes);
  return($high_confidence);
}


sub write_input_file {
  my ($self,$fasta_records) = @_;

  my $output_file = $self->create_filename();
  open(OUT,">".$output_file);
  foreach my $fasta_record (@$fasta_records) {
    say OUT $fasta_record;
  }
  close OUT;

  return($output_file);
}


sub set_transcript_descriptions {
  my ($self,$transcripts_by_id,$source_transcript_id_hash) = @_;

  foreach my $id (keys(%$transcripts_by_id)) {
    my $transcript = $transcripts_by_id->{$id};
    my $source_transcript = $source_transcript_id_hash->{$transcript->stable_id()};

    # add source transcript stable id as transcript attribute
    my $parent_attribute = Bio::EnsEMBL::Attribute->new(-CODE => 'proj_parent_t',-VALUE => $source_transcript->stable_id().".".$source_transcript->version());
    $transcript->add_Attributes($parent_attribute);

    my $description_string = ";parent_transcript=".$source_transcript->stable_id().".".$source_transcript->version().";mapping_coverage=".$transcript->{'cov'}.";mapping_identity=".$transcript->{'perc_id'};
    if($transcript->{'cds_description'}) {
      $description_string .= $transcript->{'cds_description'};
    }
    $transcript->description($description_string);
  }
}


sub qc_cds_sequences {
  my ($self,$transcripts_by_id,$source_transcript_id_hash) = @_;

  foreach my $id (keys(%$transcripts_by_id)) {
    my $transcript = $transcripts_by_id->{$id};
    my $source_transcript = $source_transcript_id_hash->{$transcript->stable_id()};

    if($source_transcript->translation()) {
      my ($cds_coverage,$cds_percent_id,$aligned_source_seq,$aligned_target_seq) = align_nucleotide_seqs($source_transcript->translateable_seq(),$transcript->translateable_seq());
      my $cds_description = ";cds_coverage=".$cds_coverage.";cds_identity=".$cds_percent_id;
      my $aligned_source_seq_copy = $aligned_source_seq;
      my $aligned_target_seq_copy = $aligned_target_seq;
      $aligned_source_seq_copy =~ s/\-\-\-//g;
      $aligned_target_seq_copy =~ s/\-\-\-//g;

      if($aligned_source_seq_copy =~ /\-/ or $aligned_target_seq_copy =~ /\-/) {
        $cds_description .= ";cds_gap=1";
        my $transcript_attrib = Bio::EnsEMBL::Attribute->new(-CODE => 'hidden_remark',
                                                             -VALUE => ">source_cds_align\n".$aligned_source_seq."\n>target_cds_align\n".$aligned_target_seq."\n");
        $transcript->add_Attributes($transcript_attrib);
        my $translation_attrib = Bio::EnsEMBL::Attribute->new(-CODE => 'hidden_remark',
                                                              -VALUE => ">source_translation\n".$source_transcript->translation->seq().
                                                                        "\n>target_translation\n".$transcript->translation->seq()."\n");
        $transcript->translation->add_Attributes($translation_attrib);
      } else {
        $cds_description .= ";cds_gap=0";
      }
      $transcript->{'cds_description'} = $cds_description;
    }
  }
}


sub set_cds_sequences {
  my ($self,$transcripts_by_id,$source_transcript_id_hash) = @_;

  foreach my $id (keys(%$transcripts_by_id)) {
    my $transcript = $transcripts_by_id->{$id};
    my $source_transcript = $source_transcript_id_hash->{$transcript->stable_id()};
    if($source_transcript->translation()) {
      $self->project_cds($transcript,$source_transcript);
    }
  }
}


sub project_cds {
  my ($self,$transcript,$source_transcript) = @_;

  # Here we want to utilise the alignment between the parent and target sequences to map the parent CDS start/end
  # Each transcript has the alignment stored on the object regardless of the reconstruction method, so this can be
  # a generic process. Note that the alignment stored is done on the transcipt seqs relative to whatever strand they
  # are on, so the two seqs will be in 5' to 3' orientation
  my $aligned_source_seq = $transcript->{'aligned_source_seq'};
  my $aligned_target_seq = $transcript->{'aligned_target_seq'};
  unless($aligned_source_seq and $aligned_target_seq) {
    $self->throw("Issue fetching alignment for transcript with stable_id ".$source_transcript->stable_id().", expected target transcript to have the alignment stored on it");
  }

  my $source_seq = $aligned_source_seq;
  $source_seq =~ s/\-//g;
  my $target_seq = $aligned_target_seq;
  $target_seq =~ s/\-//g;

  my @source_align_array = split('',$aligned_source_seq);
  my @target_align_array = split('',$aligned_target_seq);
  my $source_cds_pos = 0;
  my $source_align_seq_pos = 0;
  my $target_cds_pos = 0;
  my $target_align_seq_pos = 0;

  # Seen to have some issue in terms of getting the real start/end for some features. Was using +1, but seems off
  my $source_cds_start = $source_transcript->cdna_coding_start();
  my $source_cds_end = $source_transcript->cdna_coding_end();
  say "Source CDS start: ".$source_cds_start;
  say "Source CDS end: ".$source_cds_end;

  my $source_cds_start_index;
  my $source_cds_end_index;
  my $target_cds_start_index;
  my $target_cds_end_index;
  my $in_feature_alignment = 0;
  my $out_of_feature_alignment = 0;
  for(my $i=0; $i < scalar(@source_align_array); $i++) {
    my $source_value = $source_align_array[$i];
    my $target_value = $target_align_array[$i];
    if($source_value ne '-') {
      $source_cds_pos++;
    }

    if($target_value ne '-') {
      $target_cds_pos++;
    }

    if($source_cds_pos == $source_cds_start) {
      $source_cds_start_index = $source_cds_pos - 1;
      say "Index of the feature start in source seq: ".$source_cds_start_index;
      $in_feature_alignment = 1;
    }

    if($in_feature_alignment and $target_value ne '-' and !(defined($target_cds_start_index))) {
      $target_cds_start_index = $target_cds_pos - 1;
    }

    if($source_cds_pos == $source_cds_end) {
      $source_cds_end_index = $source_cds_pos - 1;
      $target_cds_end_index = $target_cds_pos - 1;
      say "Index of the feature end in source seq: ".$source_cds_end_index;
      $in_feature_alignment = 0;
    }
  }

  my $source_feature_length = $source_cds_end_index - $source_cds_start_index + 1;
  my $target_feature_length = $target_cds_end_index - $target_cds_start_index + 1;
  my $recovered_source_feature_seq = substr($source_seq,$source_cds_start_index,$source_feature_length);
  my $recovered_target_feature_seq = substr($target_seq,$target_cds_start_index,$target_feature_length);
  say "For transcript ".$source_transcript->stable_id()." original and mapped CDS sequences:\n".$recovered_source_feature_seq."\n".$recovered_target_feature_seq;

  my $cds_start_exon;
  my $cds_end_exon;
  my $cds_start_offset = 0;
  my $cds_end_offset = 0;
  my $cumulative_length = 0;
  my $exons = $transcript->get_all_Exons();
  foreach my $exon (@$exons) {
    say "Checking exons with the following start/end for CDS coords: ".$exon->start()."/".$exon->end();
    if(($target_cds_start_index >= $cumulative_length) and ($target_cds_start_index <= $cumulative_length + $exon->length())) {
      $cds_start_exon = $exon;
      if($transcript->strand == 1) {
        $cds_start_offset = $target_cds_start_index - $cumulative_length + 1;
      } else {
        say "FERGAL CDS START DEBUG: ".$exon->length." - (".($exon->length - $target_cds_start_index).") + 1";
        $cds_start_offset = $exon->length - ($exon->length - $target_cds_start_index) - $cumulative_length + 1;
      }
    }

    if(($target_cds_end_index >= $cumulative_length) and ($target_cds_end_index <= $cumulative_length + $exon->length())) {
      $cds_end_exon = $exon;
      if($transcript->strand == 1) {
        $cds_end_offset = $target_cds_end_index - $cumulative_length + 1;
      } else {
        say "FERGAL CDS END DEBUG: ".$exon->length." - (".($exon->length - $target_cds_end_index).") + 1";
        $cds_end_offset = $exon->length - ($exon->length - $target_cds_end_index) - $cumulative_length + 1;
      }
    }
    $cumulative_length += $exon->length();
  }

  unless($cds_start_exon and $cds_end_exon) {
    # Note this is only for testing, if this happens we should just put a warning and do a compute_translation
    # as there will probably be plenty of cases where a bad mapping means the exons are missing
    $self->warning("For ".$source_transcript->stable_id()." couldn't find the equivalent CDS start/end exon in the alignment of the transcripts. ".
                   "Will compute translation instead");
    compute_best_translation($transcript);
    return;
  }

  say "Orig CDS start/end exon lengths: ".$source_transcript->translation->start_Exon->length."/".$source_transcript->translation->end_Exon->length;
  say "Orig CDS start/end exon phases: ".$source_transcript->translation->start_Exon->phase()."/".$source_transcript->translation->end_Exon->end_phase();
  say "Orig CDS start exon frame: ".$source_transcript->translation->start_Exon->frame();
  say "Orig CDS start/end exon offsets: ".$source_transcript->translation->start()."/".$source_transcript->translation->end();
  say "CDS start/end index: ".$target_cds_start_index."/".$target_cds_end_index;
  say "CDS start/end exon lengths: ".$cds_start_exon->length()."/".$cds_end_exon->length();
  say "CDS start/end exon offsets: ".$cds_start_offset."/".$cds_end_offset;
  my $translation = Bio::EnsEMBL::Translation->new();
  $translation->start_Exon($cds_start_exon);
  $translation->start($cds_start_offset);
  $translation->end_Exon($cds_end_exon);
  $translation->end($cds_end_offset);
  $transcript->translation($translation);

  # Set the phases
  calculate_exon_phases($transcript,0);
  # Set the phase of the start exon to match the start exon of the source. We probably want some checks on this in general to ensure there's a
  # valid alignment at the start of the CDS. If the start is not conserved the best option would probably be to clip the start in the target to
  # the nearest aligned codon
  $transcript->translation->start_Exon->phase($source_transcript->translation->start_Exon->phase());

  say "For transcript ".$source_transcript->stable_id()." translateable seq:\n".$transcript->translateable_seq();
  say "For transcript ".$source_transcript->stable_id()." translation seq:\n".$transcript->translation->seq();
  say "Source transcript translation seq:\n".$source_transcript->translation->seq();
}


sub fix_cds_issues {
  my ($self,$transcripts_by_id) = @_;

  # First we should consider any transcript that has a non-perfect cds alignment


  # After the CDS has been updated in terms of potential issues with frameshifts or gaps, then examine for inframe stops and remove them
  # as needed. Polymorphic pseudogenes are expected to have inframe stops, so ignore these
  my $max_stops = 999;
  foreach my $id (keys(%$transcripts_by_id)) {
    my $transcript = $transcripts_by_id->{$id};
    if ($transcript->translation and $transcript->translation->seq() =~ /\*/ and $transcript->biotype ne 'polymorphic_pseudogene') {
      $transcript = replace_stops_with_introns($transcript,$max_stops);
      if ($transcript and $transcript->translate()->seq !~ /\*/) {
          #;#push(@$minimap_transcripts,$transcript_after_replaced_stops);
      } elsif ($transcript and !$transcript->translate()) {
        print STDERR "minimap transcript (seq_region_start,seq_region_end,seq_region_strand,seq_region_name) (".$transcript->seq_region_start().
                     ",".$transcript->seq_region_end().",".$transcript->seq_region_strand().",".$transcript->seq_region_name().") does not translate after replacing a maximum of $max_stops stops. Discarded.\n";
        $transcript->translation(undef);
        $transcript->biotype("processed_transcript");
      } elsif ($transcript) {
        print STDERR "minimap transcript (seq_region_start,seq_region_end,seq_region_strand,seq_region_name) (".$transcript->seq_region_start().
                     ",".$transcript->seq_region_end().",".$transcript->seq_region_strand().",".$transcript->seq_region_name().") has more than the maximum of $max_stops stops or it has stops after replace_stops_with_introns. Discarded.\n";
        $transcript->translation(undef);
        $transcript->biotype("processed_transcript");
      } else {
        print STDERR "No transcript defined in 'fix_cds_issues()' after 'replace_stops_with_introns()'.\n";
      }
    }

    if ($transcript and $transcript->translation()) {
      $self->check_and_fix_translation_boundaries($transcript);
    }
  }
}


sub calculate_region_boundary {
  my ($self,$paf_results) = @_;

  my $top_paf_result = shift(@$paf_results);
  my $top_paf_strand = ${$top_paf_result}[4];
  my $top_paf_source_genomic_start = ${$top_paf_result}[2];
  my $top_paf_source_genomic_end = ${$top_paf_result}[3];
  my $top_paf_target_genomic_start = ${$top_paf_result}[7];
  my $top_paf_target_genomic_end = ${$top_paf_result}[8];
  my $top_paf_target_genomic_name = ${$top_paf_result}[5];
  my $source_genomic_length = ${$top_paf_result}[1];
  my $target_genomic_length = ${$top_paf_result}[6];

  say "Top paf source start/end: ".$top_paf_source_genomic_start."/".$top_paf_source_genomic_end;
  say "Top paf target start/end: ".$top_paf_target_genomic_start."/".$top_paf_target_genomic_end;
  # This will control the variability of the gap between two hits on the target relative to the
  # gap in coverage on the source sequence

  my $cluster_source_genomic_start = $top_paf_source_genomic_start;
  my $cluster_source_genomic_end = $top_paf_source_genomic_end;
  my $cluster_target_genomic_start = $top_paf_target_genomic_start;
  my $cluster_target_genomic_end = $top_paf_target_genomic_end;

  my $target_flanking = 500;
  my $cluster_scaling_ratio = 1.5;

  my $source_missing_coverage = $source_genomic_length - ($top_paf_source_genomic_end - $top_paf_source_genomic_start) + 1;
  my $source_missing_left = $top_paf_source_genomic_start;
  my $source_missing_right = $source_genomic_length - $top_paf_source_genomic_end;
  my $source_scaled_left = ceil($cluster_scaling_ratio * $source_missing_left) + $target_flanking;
  my $source_scaled_right = ceil($cluster_scaling_ratio * $source_missing_right) + $target_flanking;

  if($top_paf_strand eq '-') {
    my $tmp = $source_scaled_left;
    $source_scaled_left = $source_scaled_right;
    $source_scaled_right = $tmp;
  }

  say "Say source missing coverage/left/right: ".$source_missing_coverage."/".$source_missing_left."/".$source_missing_right;

  say "Source cluster coverage start-end: ".$cluster_source_genomic_start."-".$cluster_source_genomic_end;
  say "Unadjusted cluster details: Start: ".$cluster_target_genomic_start.", End: ".$cluster_target_genomic_end.", Strand: ".$top_paf_strand.", Name: ".$top_paf_target_genomic_name;

  $cluster_target_genomic_start -= $source_scaled_left;
  $cluster_target_genomic_end += $source_scaled_right;

  if($cluster_target_genomic_start <= 0) {
    $cluster_target_genomic_start = 1;
  }

  if($cluster_target_genomic_end > $target_genomic_length) {
    $cluster_target_genomic_end = $target_genomic_length;
  }

  say "Final cluster details: Start: ".$cluster_target_genomic_start.", End: ".$cluster_target_genomic_end.", Strand: ".$top_paf_strand.", Name: ".$top_paf_target_genomic_name;
  return([$cluster_target_genomic_start,$cluster_target_genomic_end,$top_paf_strand,$top_paf_target_genomic_name]);
}



sub chain_paf_results {
  my ($self,$paf_results) = @_;

  my $top_paf_result = shift(@$paf_results);
  my $top_paf_strand = ${$top_paf_result}[4];
  my $top_paf_source_genomic_start = ${$top_paf_result}[2];
  my $top_paf_source_genomic_end = ${$top_paf_result}[3];
  my $top_paf_target_genomic_start = ${$top_paf_result}[7];
  my $top_paf_target_genomic_end = ${$top_paf_result}[8];
  my $top_paf_target_genomic_name = ${$top_paf_result}[5];
  my $source_genomic_length = ${$top_paf_result}[1];

  say "Top paf source start/end: ".$top_paf_source_genomic_start."/".$top_paf_source_genomic_end;
  say "Top paf target start/end: ".$top_paf_target_genomic_start."/".$top_paf_target_genomic_end;
  # This will control the variability of the gap between two hits on the target relative to the
  # gap in coverage on the source sequence

  my $cluster_source_genomic_start = $top_paf_source_genomic_start;
  my $cluster_source_genomic_end = $top_paf_source_genomic_end;
  my $cluster_target_genomic_start = $top_paf_target_genomic_start;
  my $cluster_target_genomic_end = $top_paf_target_genomic_end;

  my $min_allowed_gap = 500;
  my $small_cluster_gap_size = 5000;
  my $small_cluster_gap_ratio = 2.5;
  my $large_cluster_gap_ratio = 1.2;
  foreach my $paf_result (@$paf_results) {
    my $hit_strand = ${$paf_result}[4];
    my $hit_genomic_name = ${$paf_result}[5];
    my $hit_source_genomic_start = ${$paf_result}[2];
    my $hit_source_genomic_end = ${$paf_result}[3];
    my $hit_target_genomic_start = ${$paf_result}[7];
    my $hit_target_genomic_end = ${$paf_result}[8];
    my $hit_target_genomic_name = ${$paf_result}[5];

    unless($hit_strand eq $top_paf_strand and $hit_genomic_name eq $top_paf_target_genomic_name) {
      next;
    }

    # Check if there's a feature overlap with the current cluster boundaries. If so skip
    if($self->coords_overlap($hit_source_genomic_start,$hit_source_genomic_end,$top_paf_source_genomic_start,$top_paf_source_genomic_end) or
       $self->coords_overlap($hit_target_genomic_start,$hit_target_genomic_end,$top_paf_target_genomic_start,$top_paf_target_genomic_end)) {
      next;
    }

    say "Hit source start/end: ".$hit_source_genomic_start."/".$hit_source_genomic_end;
    say "Hit target start/end: ".$hit_target_genomic_start."/".$hit_target_genomic_end;

    # At this point find the gap between the source region and the top hit source region, then do the same with the target regions
    # and determine if the source gap is similar to the target gap and if it is then adjust the cluster boundaries on the target
    my $source_gap = $self->calculate_coord_distance($top_paf_source_genomic_start,$top_paf_source_genomic_end,$hit_source_genomic_start,$hit_source_genomic_end);
    my $scaled_source_gap = $source_gap;
    if($source_gap <= $small_cluster_gap_size) {
      $scaled_source_gap = ($source_gap * $small_cluster_gap_ratio) + $min_allowed_gap;
    } else {
      $scaled_source_gap = ($source_gap * $large_cluster_gap_ratio) + $min_allowed_gap;;
    }

    my $target_gap = $self->calculate_coord_distance($top_paf_target_genomic_start,$top_paf_target_genomic_end,$hit_target_genomic_start,$hit_target_genomic_end);

    say "Source gap: ".$source_gap;
    say "Scaled source gap: ".$scaled_source_gap;
    say "Target gap: ".$target_gap;

    if($target_gap <= $scaled_source_gap) {
      if($hit_target_genomic_start < $cluster_target_genomic_start) {
        $cluster_target_genomic_start = $hit_target_genomic_start;
      }

      if($hit_target_genomic_end > $cluster_target_genomic_end) {
        $cluster_target_genomic_end = $hit_target_genomic_end;
      }

      if($hit_source_genomic_start < $cluster_source_genomic_start) {
        $cluster_source_genomic_start = $hit_source_genomic_start;
      }

      if($hit_source_genomic_end > $cluster_source_genomic_end) {
        $cluster_source_genomic_end = $hit_source_genomic_end;
      }
    }
  }

  say "Source cluster coverage start-end: ".$cluster_source_genomic_start."-".$cluster_source_genomic_end;
  say "Unadjuster cluster details: Start: ".$cluster_target_genomic_start.", End: ".$cluster_target_genomic_end.", Strand: ".$top_paf_strand.", Name: ".$top_paf_target_genomic_name;
  my $adjust_left = $cluster_target_genomic_start - $cluster_source_genomic_start;
  my $adjust_right = $cluster_target_genomic_end + ($source_genomic_length - $cluster_source_genomic_end);
  say "Final cluster details: Start: ".$cluster_target_genomic_start.", End: ".$cluster_target_genomic_end.", Strand: ".$top_paf_strand.", Name: ".$top_paf_target_genomic_name;
  return([$cluster_target_genomic_start,$cluster_target_genomic_end,$top_paf_strand,$top_paf_target_genomic_name]);
}


sub coords_overlap {
  my ($self,$s1,$e1,$s2,$e2) = @_;

  if (($s1 <= $e2) and ($e1 >= $s2)) {
    return 1;
  }
  return 0;
}


sub calculate_coord_distance {
  my ($self,$s1,$e1,$s2,$e2) = @_;

  my $dist = 0;
  if($s1 < $s2) {
    $dist = $s2-$e1;
  } else {
    $dist = $s1-$e2;
  }

  if($dist < 0) {
    $dist = 0;
  }
  return $dist;
}



sub check_for_small_gene {
  my ($self,$transcripts) = @_;

  # For the moment it's just going to check for genes with single exon short transcripts
  my $small_length = 100;
  my $is_small = 1;
  foreach my $transcript (@$transcripts) {
    if($transcript->length > $small_length) {
      $is_small = 0;
    }

    my $exons = $transcript->get_all_Exons();
    if(scalar(@$exons) > 1) {
      $is_small = 0;
    }
  }
  return($is_small);
}


sub print_transcript_stats {
  my ($self,$transcripts_by_id,$tag) = @_;

  say "Transcript stats for ".$tag;
  foreach my $id (keys(%$transcripts_by_id)) {
    my $transcript = $transcripts_by_id->{$id};
    say "  ID: ".$transcript->stable_id().", Cov: ".$transcript->{'cov'}.", Perc id: ".$transcript->{'perc_id'};
  }
}


sub map_gene_minimap {
  my ($self,$gene,$source_transcripts,$target_genomic_start,$target_region_slice,$target_strand,
      $target_genome_file,$source_transcript_id_hash,$max_intron_size,$target_adaptor,$target_slice_adaptor,$best_transcripts_by_id) = @_;

  my $source_transcript_fasta_seqs = [];
  my $source_transcripts_to_map = $self->filter_transcripts_to_map($source_transcripts,$best_transcripts_by_id);

  say "Processing a total of ".scalar(@$source_transcripts_to_map)." source transcripts via local minimap";

  foreach my $source_transcript (@$source_transcripts_to_map) {
    say "Writing ".$source_transcript->stable_id()." to file for mapping";
    my $source_transcript_sequence = $source_transcript->seq->seq();
    my $fasta_record = ">".$source_transcript->dbID()."\n".$source_transcript_sequence;
    push(@$source_transcript_fasta_seqs,$fasta_record);
  }

  say "Generating initial set of minimap2 mappings";
  my $minimap_transcripts = $self->generate_minimap_transcripts($source_transcript_fasta_seqs,$target_genome_index,$target_adaptor,$target_genomic_start,$max_intron_size);

  my $transcripts_by_id = {};
  foreach my $transcript (@$minimap_transcripts) {
    $transcripts_by_id->{$transcript->stable_id} = $transcript;
    my $source_transcript = $source_transcript_id_hash->{$transcript->stable_id};

    my ($mapped_coverage,$mapped_percent_id,$aligned_source_seq,$aligned_target_seq) = align_nucleotide_seqs($source_transcript->seq->seq(),$transcript->seq->seq());
    $transcript->{'cov'} = $mapped_coverage;
    $transcript->{'perc_id'} = $mapped_percent_id;
    $transcript->{'aligned_source_seq'} = $aligned_source_seq;
    $transcript->{'aligned_target_seq'} = $aligned_target_seq;
    say "Mapped transcript (".$source_transcript->stable_id()."): Coverage: ".$transcript->{'cov'}.", Percent id: ".$transcript->{'perc_id'};

#    my $description_string = "Parent: ".$source_transcript->stable_id().".".$source_transcript->version().", Coverage: ".$transcript->{'cov'}.", Perc id: ".$transcript->{'perc_id'};
#    $transcript->description($description_string);
  }
  return($transcripts_by_id);

#  my $bad_minimap_transcripts = $self->check_mapping_quality($output_minimap_transcripts,$source_transcript_id_hash,$good_transcripts);
#  say "Number of good transcripts after minimap2: ".scalar(@$good_transcripts);
#  say "Number of bad transcripts after minimap2: ".scalar(@$bad_minimap_transcripts);

#  say "Running exonerate on bad transcripts";
#  my $output_exonerate_transcripts = $self->generate_exonerate_transcripts($bad_minimap_transcripts,$source_transcript_id_hash,$target_region_slice,$target_slice_adaptor,$max_intron_size);
#  my $bad_exonerate_transcripts = $self->check_mapping_quality($output_exonerate_transcripts,$source_transcript_id_hash,$good_transcripts);
#  say "Number of good transcripts after exonerate: ".scalar(@$good_transcripts);
#  say "Number of bad transcripts after exonerate: ".scalar(@$bad_exonerate_transcripts);

  # Store these in a hash for making selecting the best based on dbID easier
#  foreach my $transcript (@$good_transcripts) {
#    $good_transcripts_hash->{$transcript->stable_id()} = $transcript;
#  }

#  my $bad_minimap_transcripts_hash = {};
#  foreach my $transcript (@$bad_minimap_transcripts) {
#    $bad_minimap_transcripts_hash->{$transcript->stable_id()} = $transcript;
#  }

#  my $bad_exonerate_transcripts_hash = {};
#  foreach my $transcript (@$bad_exonerate_transcripts) {
#    $bad_exonerate_transcripts_hash->{$transcript->stable_id()} = $transcript;
#  }

  # This is bit of a mess as I only realised that many of the hardest edge cases can be recovered with exonerate, but might totally fail on minimap2
  # Thus there needed to be a step added for when the minimap transcript is missing to begin wit, because the way the rest of the code is designed
  # a transcript might just not be found via minimap and then all that will happen is that it gets globally aligned via minimap later and likely not
  # found again. This just ensures that a local exonerate is run on the region for totally missing transcripts, not just ones that fail the minimap
  # cut-offs. To refactor this at some point it would be better to rely less on arrays of the mapped transcripts and instead go through the id hashes
#  say "Checking for missing transcripts to perform exonerate on region";
#  my $initial_missing_source_transcripts = $self->list_missing_transcripts($good_transcripts_hash,$bad_exonerate_transcripts_hash,$source_transcripts);
#  say "Found ".scalar(@$initial_missing_source_transcripts)." initial missing transcripts";
#  my $missing_exonerate_transcripts = [];
#  foreach my $missing_source_transcript (@$initial_missing_source_transcripts) {
#    say "Initial missing source transcript: ".$missing_source_transcript->stable_id();
#    my $exonerate_transcripts = $self->run_exonerate($missing_source_transcript,$target_region_slice,$target_slice_adaptor,$max_intron_size);
#    push(@$missing_exonerate_transcripts,@$exonerate_transcripts);
#  }

#  say "After mapping initial missing transcripts with exonerate in region found ".scalar(@$missing_exonerate_transcripts)." exonerate transcripts";

#   my $bad_initial_missing_transcripts = $self->check_mapping_quality($missing_exonerate_transcripts,$source_transcript_id_hash,$good_transcripts);
#   foreach my $transcript (@$bad_initial_missing_transcripts) {
#     $bad_exonerate_transcripts_hash->{$transcript->stable_id()} = $transcript;
#   }

   # Select the best transcript foreach pair across both hashes, and any transcripts that are unique to either hash
   # Then put the corresponding source transcripts into an array to do a global mapping later
#   say "Selecting best transcripts out of bad minimap2/exonerate transcripts";
#   $best_bad_transcripts = $self->select_best_transcripts($bad_minimap_transcripts_hash,$bad_exonerate_transcripts_hash);
#   foreach my $transcript (@$best_bad_transcripts) {
     # If a transcript was bad from minimap and good after exonerate, then we want to skip over it
#     if($good_transcripts_hash->{$transcript->stable_id()}) {
#       next;
#     }

#     $best_bad_transcripts_hash->{$transcript->stable_id()} = $transcript;
#     push(@$bad_source_transcripts,$source_transcript_id_hash->{$transcript->stable_id()});
#   }

#   say "Number of bad transcripts after selection: ".scalar(@$bad_source_transcripts);
}


sub map_gene_exonerate {
  my ($self,$source_transcripts,$target_genomic_start,$target_region_slice,$target_strand,
      $target_genome_file,$source_transcript_id_hash,$max_intron_size,$target_adaptor,$target_slice_adaptor, $best_transcripts_by_id) = @_;

  my $source_transcripts_to_map = $self->filter_transcripts_to_map($source_transcripts,$best_transcripts_by_id);

  say "Processing a total of ".scalar(@$source_transcripts_to_map)." source transcripts using Exonerate";

  my $exonerate_transcripts = $self->generate_exonerate_transcripts($source_transcripts_to_map,$source_transcript_id_hash,$target_region_slice,$target_slice_adaptor,$max_intron_size);

  my $transcripts_by_id = {};
  foreach my $transcript (@$exonerate_transcripts) {
    $transcripts_by_id->{$transcript->stable_id} = $transcript;
  }

  return($transcripts_by_id);
}


sub update_best_transcripts {
  my ($self,$best_transcripts_by_id,$new_transcripts_by_id,$coverage_threshold,$perc_id_threshold) = @_;

#  say "FERGAL DUMPER: ".Dumper($new_transcripts_by_id);
#  say "FERGAL REF: ".ref($new_transcripts_by_id);
#  foreach my $test_id (keys(%$new_transcripts_by_id)) {
#    say "TEST ID:".$test_id;
#  }

  my @ids = keys(%{$new_transcripts_by_id});
  foreach my $id (@ids) {
    unless($best_transcripts_by_id->{$id}) {
      $best_transcripts_by_id->{$id} = $new_transcripts_by_id->{$id};
    } else {
      my $best_coverage = $best_transcripts_by_id->{$id}->{'cov'};
      my $best_perc_id = $best_transcripts_by_id->{$id}->{'perc_id'};
      my $best_total = $best_coverage + $best_perc_id;
      my $transcript_coverage = $new_transcripts_by_id->{$id}->{'cov'};
      my $transcript_perc_id = $new_transcripts_by_id->{$id}->{'perc_id'};
      my $transcript_total = $transcript_coverage + $transcript_perc_id;

      if(($coverage_threshold and $perc_id_threshold) and ($transcript_coverage < $coverage_threshold or $transcript_perc_id < $perc_id_threshold)) {
        return;
      }

      if($transcript_total > $best_total) {
        # If a coverage and perc_id are specified, then we will not replace the current best model if it passes the thresholds
        # The use case is to not replace a model placed via the 2-pass approach with a globally aligned model unless the 2-pass
        # model fails the threshold. This is unlikely, but in the cases of close paralogues, a global incorrectly placed model
        # could score better than a correctly placed 2-pass model and this is slightly gatekeeping that
        if(($coverage_threshold and $perc_id_threshold) and ($best_coverage >= $coverage_threshold and $best_perc_id >= $perc_id_threshold)) {
          return;
        }
        $best_transcripts_by_id->{$id} = $new_transcripts_by_id->{$id};
      }
    }
  }
}


sub filter_transcripts_to_map {
  my ($self,$source_transcripts,$best_transcripts_by_id) = @_;

  my $coverage_cutoff = 98;
  my $identity_cutoff = 98;

  my $filtered_transcripts = [];
  foreach my $source_transcript (@$source_transcripts) {
    my $source_transcript_id = $source_transcript->dbID();
    unless($best_transcripts_by_id->{$source_transcript_id}) {
      push(@$filtered_transcripts,$source_transcript);
      next;
    }
    my $best_coverage = $best_transcripts_by_id->{$source_transcript_id}->{'cov'};
    my $best_perc_id = $best_transcripts_by_id->{$source_transcript_id}->{'perc_id'};
    unless($best_coverage >= $coverage_cutoff and $best_perc_id >= $identity_cutoff) {
      push(@$filtered_transcripts,$source_transcript);
    }
  }
  return($filtered_transcripts);
}


sub project_gene_coords {
  my ($self,$gene,$source_transcripts,$target_genomic_start,$target_region_slice,$target_strand) = @_;

  my $transcripts_by_id = {};

  if($self->no_projection()) {
    return($transcripts_by_id);
  }

  # As this uses a MAFFT alignment, there are some limits in terms of getting that to run, for the handful of really long genes
  # just skip and rely on the minimap approach
  # Note, skipping as I'm not sure big genes are really much of an issue for MAFFT, it's likely other factors that determine if
  # something gets stuck and projection is likely to be better in the case of big genes
  my $gene_length_threshold = 250000;
#  if($gene->length >= $gene_length_threshold) {
#    say "The source gene is longer than the gene length threshold, will align via minimap in the target region instead of projecting";
#    return($transcripts_by_id);
#  }

  my $gene_genomic_seqs_hash = $self->gene_genomic_seqs_hash();
  my $source_genomic_seq_info = $gene_genomic_seqs_hash->{$gene->dbID()};

  unless($source_genomic_seq_info) {
    $self->throw("Could not find the genomic seq info for source gene with dbID: ".$gene->dbID());
  }

  my $source_region_start = ${$source_genomic_seq_info}[0];
  my $source_region_end = ${$source_genomic_seq_info}[1];
  my $source_genome_seq = ${$source_genomic_seq_info}[2];

  my $exons = $gene->get_all_Exons();
  my $target_genome_seq = $target_region_slice->seq();

  # Put back onto forward strand for simplicity
  if($gene->strand != 1) {
    $source_genome_seq = $self->revcomp($source_genome_seq);
  }

  if($target_strand != $gene->strand) {
    $target_genome_seq = $self->revcomp($target_genome_seq);
  }

  my $projected_exons_by_id = {};
  # Sometimes issues with the MySQL server disconnecting on long alignments
  $self->source_adaptor->dbc->disconnect_when_inactive(1);
  $self->target_adaptor->dbc->disconnect_when_inactive(1);

  my $coverage = 0;
  my $percent_id = 0;
  my $aligned_source_seq;
  my $aligned_target_seq;
  eval {
    ($coverage,$percent_id,$aligned_source_seq,$aligned_target_seq) = align_nucleotide_seqs($source_genome_seq,$target_genome_seq);
  };
  if ($@) {
    $self->warning("Issue with running MAFFT on target region");
  } else {
    say "Aligned source and target regions with MAFFT: Coverage: ".$coverage.", Percent id: ".$percent_id;
  }

  # Check if there's an issue with the alignment, MAFFT can sometimes get regions with duplications completely wrong
  if($coverage < 90 or $percent_id < 95 and $gene->length < $gene_length_threshold) {
    say "MAFFT coverage/id failed threshold, will re-align with MUSCLE and take the alignment with the highest combined coverage/id";

    my $muscle_coverage = 0;
    my $muscle_percent_id = 0;
    my $muscle_aligned_source_seq;
    my $muscle_aligned_target_seq;
#    eval {
#      ($muscle_coverage,$muscle_percent_id,$muscle_aligned_source_seq,$muscle_aligned_target_seq) = align_nucleotide_seqs($source_genome_seq,$target_genome_seq,'muscle');
#    }; if ($@) {
#      $self->warning("Issue with running MUSCLE on target region");
#    } else {
#      say "Aligned source and target regions with MUSCLE: Coverage: ".$muscle_coverage.", Percent id: ".$muscle_percent_id;
#      if (($muscle_coverage + $muscle_percent_id) > ($coverage + $percent_id)) {
#         $aligned_source_seq = $muscle_aligned_source_seq;
#         $aligned_target_seq = $muscle_aligned_target_seq;
#      }
#    }
  }

  unless($coverage and $percent_id) {
    say "No coverage/percent id calculated therefore not proceeding with projection";
    return($transcripts_by_id);
  }

  foreach my $exon (@$exons) {
    my $exon_region_start = $exon->seq_region_start();
    my $exon_region_end = $exon->seq_region_end();
    say "Source region start/end: ".$source_region_start."/".$source_region_end;
    say "Source exon region start/end: ".$exon_region_start."/".$exon_region_end;
    say "Source exon strand: ".$exon->strand();
    say "Target strand: ".$target_strand;

    my $projected_exon = $self->project_feature(undef,$exon,$source_region_start,$exon_region_start,$exon_region_end,$aligned_source_seq,$aligned_target_seq,$target_region_slice,$target_strand);

    if($projected_exon) {
      my ($proj_coverage,$proj_percent_id,$aligned_source_seq,$aligned_target_seq) = align_nucleotide_seqs($exon->seq->seq(),$projected_exon->seq->seq());
      $projected_exon->{'cov'} = $proj_coverage;
      $projected_exon->{'perc_id'} = $proj_percent_id;
      $projected_exon->{'source_stable_id'} = $exon->stable_id();
      $projected_exon->{'source_length'} = $exon->length();
#      say "Projected exon (".$projected_exon->{'source_stable_id'}."): Coverage: ".$projected_exon->{'cov'}.", Percent id: ".$projected_exon->{'perc_id'};
#      say "Alignment:\n".$aligned_source_seq."\n".$aligned_target_seq;
      $projected_exons_by_id->{$projected_exon->{'source_stable_id'}} = $projected_exon;
    } else {
      say "Failed to project exon (".$projected_exon->{'source_stable_id'}.")";
    }
  }

  foreach my $source_transcript (@$source_transcripts) {
#    unless($source_transcript->stable_id() eq 'ENST00000556119') {
#      next;
#    }
    my $projected_transcript = $self->reconstruct_transcript($source_transcript,$projected_exons_by_id);
    if($projected_transcript) {
      $projected_transcript->{'annotation_method'} = 'alignment_projection';
      $transcripts_by_id->{$projected_transcript->stable_id()} = $projected_transcript;
    }
  }
  return($transcripts_by_id);
}


sub reconstruct_transcript {
  my ($self,$source_transcript,$projected_exons_by_id) = @_;

  my $source_exons = $source_transcript->get_all_Exons();

  my $projected_exons = [];
  foreach my $source_exon (@$source_exons) {
    my $projected_exon = $projected_exons_by_id->{$source_exon->stable_id};
    if($projected_exon) {
      push(@$projected_exons,$projected_exon);
    }
  }

  unless(scalar(@$projected_exons)) {
    return;
  }

  my $projected_transcript = Bio::EnsEMBL::Transcript->new(-EXONS => $projected_exons);
  $projected_transcript->slice(${$projected_exons}[0]->slice());
  $projected_transcript->analysis($self->analysis());
  $projected_transcript->biotype($source_transcript->biotype());
  $projected_transcript->stable_id($source_transcript->dbID());
  $projected_transcript->version($source_transcript->version());
  $projected_transcript->{'source_stable_id'} = $source_transcript->stable_id();
  $projected_transcript->{'source_biotype_group'} = $source_transcript->get_Biotype->biotype_group();
  $projected_transcript->{'source_length'} = $source_transcript->length();

  my ($proj_coverage,$proj_percent_id,$aligned_source_seq,$aligned_target_seq) = align_nucleotide_seqs($source_transcript->seq->seq(),$projected_transcript->seq->seq());
  $projected_transcript->{'cov'} = $proj_coverage;
  $projected_transcript->{'perc_id'} = $proj_percent_id;
  $projected_transcript->{'aligned_source_seq'} = $aligned_source_seq;
  $projected_transcript->{'aligned_target_seq'} = $aligned_target_seq;
  say "Projected transcript (".$projected_transcript->{'source_stable_id'}."): Coverage: ".$projected_transcript->{'cov'}.", Percent id: ".$projected_transcript->{'perc_id'};
  say "Alignment:\n".$aligned_source_seq."\n".$aligned_target_seq;

#  my $description_string = "Parent: ".$source_transcript->stable_id().".".$source_transcript->version().", Coverage: ".$proj_coverage.", Perc id: ".$proj_percent_id;
#  $projected_transcript->description($description_string);

  return($projected_transcript);
}


sub project_feature {
  my ($self,$transcript,$exon,$source_region_start,$feature_start,$feature_end,$aligned_source_seq,$aligned_target_seq,$target_region_slice,$target_strand) = @_;

  # At this point we have the aligned seqs and many coords
  # The source region start/feature start/feature end are used to calculate the offset of the feature start/end in the aligned source seq
  # Once they are located, then the equivalent positions (if they exist) can be calculated in the target seq
  # After they've been located, we can then calculate the coords of the feature in the target and build it
  # There are various issues here. If the positions are at gaps in the alignment, then what to do? Most straightfoward thing might be to
  # scan and include all aligned bases between the source start/end and fill in with equivalent number of missing bases
  # Maybe look for the closest non-gap position to the coord of the start end that's within the genomic region of the source feature and then
  # calculate the base offset from there and include the equivalent number of missing bases from the target seq
  # Actually the above is probably just overcomplicated, just take all the bases aligned to positions within the feature in the source
  # then when you calculate coverage it will help sort out bad ones versus good ones
  # Also need to do something in terms of if there's a cds, since this is likely to cover IG genes with weird cds features
  # Could decide to drop the UTR in these cases (if there is UTR, since that's very unlikely), then that would make it more
  # straightforward in terms of just making the whole feature coding
#  say "Aligned seqs:";
#  say $aligned_source_seq;
#  say $aligned_target_seq;

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
  my $source_feature_seq_start =  $feature_start - $source_region_start;
  my $source_feature_seq_end = $feature_end - $source_region_start;
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
      say "Index of the feature start in source seq: ".$source_seq_start_index;
      $in_feature_alignment = 1;
    }

    if($in_feature_alignment and $target_value ne '-' and !(defined($target_seq_start_index))) {
      $target_seq_start_index = $target_seq_pos;
    }

    if($source_seq_pos == $source_feature_seq_end) {
      $source_seq_end_index = $source_seq_pos;
      $target_seq_end_index = $target_seq_pos;
      say "Index of the feature end in source seq: ".$source_seq_end_index;
      $in_feature_alignment = 0;
    }
  }

  my $source_feature_length = $source_seq_end_index - $source_seq_start_index + 1;
  my $target_feature_length = $target_seq_end_index - $target_seq_start_index + 1;
  my $recovered_source_feature_seq = substr($source_seq,$source_seq_start_index,$source_feature_length);
  if($exon->strand != 1) {
    $recovered_source_feature_seq = $self->revcomp($recovered_source_feature_seq);
  }

  unless(defined($target_seq_start_index) and defined($target_seq_end_index)) {
    $self->warning("Issue with recovering start/end of the exon feature in target alignment sequence, not building exon");
    return;
  }

  my $recovered_target_feature_seq = substr($target_seq,$target_seq_start_index,$target_feature_length);
  if($target_strand != 1) {
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


sub build_projected_exon {
  my ($self,$transcript,$exon,$seq_start_index,$seq_end_index,$region_slice,$target_strand) = @_;

  say "Region slice: ".$region_slice->name();
  my $region_start = $region_slice->seq_region_start();
  my $region_end = $region_slice->seq_region_end();
  say "Region start: ".$region_start;
  say "Region end: ".$region_end;

  my $parent_slice = $region_slice->seq_region_Slice();
#  say "Parent slice: ".$parent_slice->name();

  my $exon_start;
  my $exon_end;

  say "Start/End index: ".$seq_start_index."/".$seq_end_index;
  if($target_strand == $exon->strand) {
    $exon_start = $region_start + $seq_start_index;
    $exon_end = $region_start + $seq_end_index;
  } else {
    # In this case we need to reverse the coords
    $exon_end = $region_end - $seq_start_index;
    $exon_start = $exon_end - ($seq_end_index - $seq_start_index);
  }


  my $phase = -1;
  my $end_phase = -1;

  say "Projected exon start/end slice coords: ".$exon_start."/".$exon_end;
  say "Parent slice stard/end genomic coords: ".$parent_slice->seq_region_start."/".$parent_slice->seq_region_end;
  my $projected_exon = Bio::EnsEMBL::Exon->new(-start     => $exon_start,
                                               -end       => $exon_end,
                                               -strand    => $target_strand,
                                               -phase     => $phase,
                                               -end_phase => $end_phase,
                                               -analysis  => $self->analysis,
                                               -slice     => $parent_slice);

  if($transcript and $exon->is_coding($transcript)) {
    $projected_exon->phase($exon->phase());
    $projected_exon->end_phase($exon->end_phase());
  }

  if($exon_start > $exon_end) {
    $self->throw("Created an exon where the start > than the end, this shouldn't be possible: ".$parent_slice->name." ".$exon->start."..".$exon->end." ".$target_strand);
  }

#  say "New exon slice: ".$projected_exon->slice->name();
#  say "New exon start/end/strand: ".$projected_exon->start().":".$projected_exon->end().":".$projected_exon->strand();
  say "New exon seq:";
  say $projected_exon->seq->seq();
  say "Original exon seq:";
  say $exon->seq->seq();

  return($projected_exon);
}


sub set_complete_transcript_cds {
  my ($self,$transcript) = @_;

  # This will take a transcript and assume it is a completely coding sequence, so will just apply a translation across the whole thing
  # This is mostly to just replicate the cds for the small, coding, single exon genes that are projected at the moment
  my $exons = $transcript->get_all_Exons();
  my $start_exon = ${$exons}[0];
  my $end_exon = ${$exons}[scalar(@$exons)-1];
  my $translation = Bio::EnsEMBL::Translation->new();
  $translation->start_Exon($start_exon);
  $translation->start(1);
  $translation->end_Exon($end_exon);
  $translation->end($end_exon->length());
  $transcript->translation($translation);

#  say "New translation info:";
#  say "  Translation start: ".$transcript->translation->start();
#  say "  Translation end: ".$transcript->translation->end();
#  say "  Start exon phase: ".$transcript->start_Exon->phase();
#  say "  Start exon end phase: ".$transcript->start_Exon->end_phase();
#  say "  Start exon frame: ".$transcript->start_Exon->frame();

}


sub generate_minimap_transcripts {
  my ($self,$source_transcript_fasta_seqs,$target_genome_index,$target_adaptor,$target_genomic_start,$max_intron_size) = @_;

  # This will take in the target genomic region (via the index)
  my $coverage_cutoff = 80;
  my $perc_id_cutoff = 80;
  my $max_stops = 999;
  my $minimap_transcripts = [];
  my $analysis = $self->analysis();
  my $program = $self->program();
  my $paftools = $self->paftools_path();
  my $source_input_file = $self->write_input_file($source_transcript_fasta_seqs);
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Minimap2->new(
       -analysis                 => $analysis,
       -program                  => $program,
       -paftools_path            => $paftools,
       -genome_index             => $target_genome_index,
       -input_file               => $source_input_file,
       -database_adaptor         => $target_adaptor,
       -skip_introns_check       => 1,
       -add_offset               => $target_genomic_start - 1,
       -skip_compute_translation => 1,
       -max_intron_size          => $max_intron_size,
       -perc_id                  => $perc_id_cutoff,
       -coverage                 => $coverage_cutoff,
  );

  $runnable->run();

  my $output_genes = $runnable->output();
  foreach my $gene (@$output_genes) {
    say "FERGAL DEBUG MINIMAP LOCAL OUTPUT!!!!!!!";
    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts) {
#      my $transcript_after_replaced_stops = $transcript;
#      if ($transcript->translate() and $transcript->translate()->seq()) {
#        $transcript_after_replaced_stops = replace_stops_with_introns($transcript,$max_stops);
#        if ($transcript_after_replaced_stops and $transcript_after_replaced_stops->translate()->seq !~ /\*/) {
#          ;#push(@$minimap_transcripts,$transcript_after_replaced_stops);
#        } elsif (!$transcript_after_replaced_stops->translate()) {
#          print STDERR "minimap transcript (seq_region_start,seq_region_end,seq_region_strand,seq_region_name) (".$transcript->seq_region_start().",".$transcript->seq_region_end().",".$transcript->seq_region->strand().",".$transcript->seq_region_name().") does not translate after replacing a maximum of $max_stops stops. Discarded.\n";
#          $transcript_after_replaced_stops->translation(undef);
#          $transcript_after_replaced_stops->biotype("processed_transcript");
#        } else {
#          print STDERR "minimap transcript (seq_region_start,seq_region_end,seq_region_strand,seq_region_name) (".$transcript->seq_region_start().",".$transcript->seq_region_end().",".$transcript->seq_region->strand().",".$transcript->seq_region_name().") has more than the maximum of $max_stops stops or it has stops after replace_stops_with_introns. Discarded.\n";
#          $transcript_after_replaced_stops->translation(undef);
#          $transcript_after_replaced_stops->biotype("processed_transcript");
#        }
#      }
#      if ($transcript_after_replaced_stops->translation()) {
#        $self->check_and_fix_translation_boundaries($transcript_after_replaced_stops);
#      }
      $transcript->{'annotation_method'} = 'minimap_local';
      push(@$minimap_transcripts,$transcript);
    }
  }

  return($minimap_transcripts);
}


sub generate_exonerate_transcripts {
  my ($self,$transcripts_to_map,$source_transcript_id_hash,$target_region_slice,$target_slice_adaptor,$max_intron_size) = @_;

  my $output_transcripts = [];

  foreach my $transcript_to_map (@$transcripts_to_map) {
    my $source_transcript = $source_transcript_id_hash->{$transcript_to_map->dbID()};
    unless($source_transcript) {
      $self->throw("Couldn't find the dbID of the transcript in the source transcript hash when attempting to run exonerate");
    }

    my $exonerate_transcripts = $self->run_exonerate($source_transcript,$target_region_slice,$target_slice_adaptor,$max_intron_size);
    say "Got ".scalar(@$exonerate_transcripts)." from Exonerate";
    foreach my $transcript (@$exonerate_transcripts) {

      my ($mapped_coverage,$mapped_percent_id,$aligned_source_seq,$aligned_target_seq) = align_nucleotide_seqs($source_transcript->seq->seq(),$transcript->seq->seq());
      $transcript->{'cov'} = $mapped_coverage;
      $transcript->{'perc_id'} = $mapped_percent_id;
      $transcript->{'aligned_source_seq'} = $aligned_source_seq;
      $transcript->{'aligned_target_seq'} = $aligned_target_seq;
      $transcript->{'annotation_method'} = 'exonerate_local';
      say "Mapped transcript (".$source_transcript->stable_id()."): Coverage: ".$transcript->{'cov'}.", Percent id: ".$transcript->{'perc_id'};
#      my $description_string = "Parent: ".$source_transcript->stable_id().".".$source_transcript->version().", Coverage: ".$transcript->{'cov'}.", Perc id: ".$transcript->{'perc_id'};
#      $transcript->description($description_string);
      push(@$output_transcripts,$transcript);
    }
  }
  return($output_transcripts);
}


sub run_exonerate {
  my ($self,$source_transcript,$target_region_slice,$target_slice_adaptor,$max_intron_size) = @_;

  my $output_transcripts = [];
  my $exonerate_length_cutoff = 15000;
  my $max_stops = 999;
  if($source_transcript->length() > $exonerate_length_cutoff) {
    say "Not running exonerate on trascript ".$source_transcript->stable_id()." as length (".$source_transcript->length().") is greater than length cut-off (".$exonerate_length_cutoff.")";
    return($output_transcripts);
  }

  say "Running Exonerate on ".$source_transcript->stable_id();
  my $annotation_features;
  my $non_coding_transcript = 0;
  if($source_transcript->translation()) {
    $annotation_features = $self->create_annotation_features($source_transcript);
  } else {
    $non_coding_transcript = 1;
  }

  my $source_transcript_header = $source_transcript->stable_id.'.'.$source_transcript->version;
  my $source_transcript_seq_object = Bio::Seq->new(-display_id => $source_transcript_header, -seq => $source_transcript->seq->seq());


  my %parameters = ();
  $parameters{-options} = "--model cdna2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE --score 500 ".
                          "--saturatethreshold 100 --dnawordlen 15 --codonwordlen 15 --dnahspthreshold 60 --bestn 1 --maxintron ".$max_intron_size;
  $parameters{-coverage_by_aligned} = 1;

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript->new(
                    -program  => '/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/opt/exonerate22/bin/exonerate',
                    -analysis => $self->analysis(),
                    -query_type     => "dna",
                    -calculate_coverage_and_pid => 1,
                    -annotation_features => $annotation_features,
                    -non_coding_transcript => $non_coding_transcript,
                     %parameters,
                 );

  $runnable->target_seqs([$target_region_slice]);
  $runnable->query_seqs([$source_transcript_seq_object]);
  $runnable->run();
  if (scalar(@{$runnable->output()})) {
    my $initial_output_transcript = ${$runnable->output()}[0];
    my $output_transcript = $self->update_exonerate_transcript_coords($initial_output_transcript,$target_slice_adaptor);
    if ($source_transcript->translation() and $output_transcript->translation()) {
      $self->check_exonerate_translation($source_transcript,$output_transcript);
    }

#    my $transcript_after_replaced_stops = $output_transcript;
#    if ($output_transcript->translate() and $output_transcript->translate()->seq()) {
#      $transcript_after_replaced_stops = replace_stops_with_introns($output_transcript,$max_stops);
#      if ($transcript_after_replaced_stops and $transcript_after_replaced_stops->translate()->seq() !~ /\*/) {
#        ;#push(@$output_transcripts,$transcript_after_replaced_stops);
#      } elsif ($transcript_after_replaced_stops and !($transcript_after_replaced_stops->translate())) {
#        print STDERR "exonerate transcript (seq_region_start,seq_region_end,seq_region_strand,seq_region_name) (".$output_transcript->seq_region_start().",".$output_transcript->seq_region_end().",".$output_transcript->seq_region_strand().",".$output_transcript->seq_region_name().") does not translate after repl#acing a maximum of $max_stops stops. Removing translation and setting biotype to processed_transcript.\n";
#        $transcript_after_replaced_stops->translation(undef);
#        $transcript_after_replaced_stops->biotype("processed_transcript");
#      } elsif ($transcript_after_replaced_stops) {
#        print STDERR "exonerate transcript (seq_region_start,seq_region_end,seq_region_strand,seq_region_name) (".$output_transcript->seq_region_start().",".$output_transcript->seq_region_end().",".$output_transcript->seq_region_strand().",".$output_transcript->seq_region_name().") has more than the maximum of $max_stops stops or it has stops after replace_stops_with_introns. Removing translation and setting biotype to processed_transcript.\n";
#        $transcript_after_replaced_stops->translation(undef);
#        $transcript_after_replaced_stops->biotype("processed_transcript");
#      } else {
#        print STDERR "exonerate transcript (seq_region_start,seq_region_end,seq_region_strand,seq_region_name) (".$output_transcript->seq_region_start().",".$output_transcript->seq_region_end().",".$output_transcript->seq_region_strand().",".$output_transcript->seq_region_name().") replace_stops_with_introns failed. Pushing transcript as it was before the attempt to replace stops.\n";
#        $transcript_after_replaced_stops = $output_transcript;
#      }
#    }
#    if ($transcript_after_replaced_stops->translation()) {
#      $self->check_and_fix_translation_boundaries($transcript_after_replaced_stops);
#    }
    $output_transcript->stable_id($source_transcript->dbID());
    push(@$output_transcripts,$output_transcript);
  }

  return($output_transcripts);
}


sub check_exonerate_translation {
  my ($self,$source_transcript,$output_transcript) = @_;

  my $source_se = $source_transcript->translation->start_Exon();
  my $source_ee = $source_transcript->translation->end_Exon();

  my $output_se = $output_transcript->translation->start_Exon();
  my $output_ee = $output_transcript->translation->end_Exon();

  # This will sort an incomplete start codon out. Exonerate will always make a CDS that's a multiple
  # of three even if the annotation features is not a multiple of three. So this code will cut the
  # equivalent amount off the end of the CDS. Note that at the moment there's no code to deal with
  # edge cases like if this moved the cds into another exon (start would be < 1)
  if($source_se->phase()) {
    my $translation = $output_transcript->translation();
    $output_se->phase($source_se->phase());
    my $end_offset = $translation->end() - $source_se->phase();
    $translation->end($end_offset);
    $output_transcript->translation($translation);
  }

  # Need some code to sort out the situation where the cds end is wrong because exonerate won't make a
  # complete cds. This can happen if the end codon of the original cds is not complete. Exonerate will
  # just truncate it to the closest codon
  my $translation_offset = length($output_transcript->translateable_seq()) % 3;
  if($translation_offset) {
    my $translation = $output_transcript->translation();
    my $end_offset = $translation->end() - $translation_offset;
    $translation->end($end_offset);
    $output_transcript->translation($translation)
  }
}

sub update_exonerate_transcript_coords {
  my ($self,$transcript,$target_slice_adaptor) = @_;

  $transcript->flush_supporting_features();
  my $slice_id = $transcript->start_Exon->seqname;

  $slice_id =~ /[^\:]+\:[^\:]+\:([^\:]+)\:([^\:]+)\:([^\:]+)\:[^\:]+$/;
  my $region_name = $1;
  my $start_offset = $2 - 1;
  my $slice = $target_slice_adaptor->fetch_by_region('toplevel',$region_name);

  my $exons = $transcript->get_all_Exons();
  foreach my $exon (@$exons) {
    $exon->start($exon->start() + $start_offset);
    $exon->end($exon->end() + $start_offset);
  }

  $transcript->start($transcript->start() + $start_offset);
  $transcript->end($transcript->end() + $start_offset);
  attach_Slice_to_Transcript($transcript,$slice);
  attach_Analysis_to_Transcript($transcript,$self->analysis());

  return($transcript);
}


sub check_mapping_quality {
  my ($self,$target_transcripts,$source_transcript_id_hash,$good_transcripts) = @_;
  my $bad_transcripts = [];
  my $coverage_cutoff_groups = {'coding' => 95,
                                'pseudogene' => 80,
                                'snoncoding' => 90,
                                'mnoncoding' => 80,
                                'lnoncoding' => 80,
                                'undefined'  => 50};

  my $percent_identity_groups = {'coding' => 95,
                                 'pseudogene' => 80,
                                 'snoncoding' => 80,
                                 'mnoncoding' => 80,
                                 'lnoncoding' => 80,
                                 'undefined'  => 50};

  my $cds_length_diff_cutoff = 0.05;
  my $genomic_span_diff_cutoff = 0.80;
  my $exonerate_length_cutoff = 15000;

  my $processed_transcripts = {};

  foreach my $transcript (@$target_transcripts) {
    say "Checking mapping quality for mapped transcript with original dbID: ".$transcript->stable_id();
    my $source_transcript = $source_transcript_id_hash->{$transcript->stable_id()};
    unless($source_transcript) {
      $self->throw("Issue with fetching source transcript for target transcript with dbID: ".$transcript->stable_id());
    }

    my $source_genomic_span = $source_transcript->seq_region_end() - $source_transcript->seq_region_start() + 1;
    my $target_genomic_span = $transcript->seq_region_end() - $transcript->seq_region_start() + 1;
    my $source_transcript_seq;
    my $transcript_seq;
    my $cds_length_diff = 0;

    # Set the description now on the minor chance the transcript doesn't have a cds that can be calculated
    my $transcript_description = ";parent_transcript=".$source_transcript->stable_id().".".$source_transcript->version();
    if ($source_transcript->display_xref()) {
      if ($source_transcript->display_xref()->display_id()) {
        $transcript_description .= ";parent_transcript_display_xref=".$source_transcript->display_xref()->display_id();
      }
    }

    $transcript->description($transcript_description);

    # add source transcript stable id as transcript attribute
    my $parent_attribute = Bio::EnsEMBL::Attribute->new(-CODE => 'proj_parent_t',-VALUE => $source_transcript->stable_id().".".$source_transcript->version());
    $transcript->add_Attributes($parent_attribute);

    if($source_transcript->translation()) {
      $source_transcript_seq = $source_transcript->translateable_seq();
      $transcript_seq = $transcript->translateable_seq();

      unless($transcript_seq) {
        compute_best_translation($transcript);
        $transcript_seq = $transcript->translateable_seq();
      }

      unless($transcript_seq) {
        $self->warning("Couldn't find an ORF in transcript mapped from ".$source_transcript->stable_id().". Source transcript has an ORF");
#        push(@$bad_transcripts,$transcript);
        next;
      }

      # This is for cds seqs only, so this is only run if there's a translation. Meaning that cds_length_diff is 0 for all other transcripts, and thus
      # they automatically pass the length check. There isn't really a risk of a non-coding transcript being longer than expected, only shorter than
      # expected (and that's handled by coverage), whereas an ORF could be longer due to selecting an incorrect upstream methionine
      if(length($transcript_seq) > length($source_transcript_seq)) {
        $cds_length_diff = 1 - (length($source_transcript_seq)/length($transcript_seq));
      }
    } else {
      $source_transcript_seq = $source_transcript->seq->seq();
      $transcript_seq = $transcript->seq->seq();
    }

    my ($coverage,$percent_id,$aligned_source_seq,$aligned_target_seq) = align_nucleotide_seqs($source_transcript_seq,$transcript_seq);
    $transcript->{'cov'} = $coverage;
    $transcript->{'perc_id'} = $percent_id;
    $transcript->{'cds_length_diff'} = $cds_length_diff;
    $transcript->{'source_stable_id'} = $source_transcript->stable_id();
    $transcript->{'source_biotype_group'} = $source_transcript->get_Biotype->biotype_group();
    $transcript->{'source_length'} = $source_transcript->length();
    say "FERGAL SOURCE SPAN: ".$source_genomic_span;
    say "FERGAL TARGET SPAN: ".$target_genomic_span;
    my $transcript_genomic_span_diff = $target_genomic_span/$source_genomic_span;
    $transcript_genomic_span_diff = sprintf("%.2f", $transcript_genomic_span_diff);
    $transcript->{'transcript_genomic_span_diff'} = $transcript_genomic_span_diff;

    $transcript_description .= ";mapping_coverage=".$coverage.";mapping_identity=".$percent_id;
    $transcript->description($transcript_description);

    # I added this in because even when minimap is explicitly told not to output secondary alignments, it very occasionally does
    # The example case was ENST00000624628, which is a large lncRNA (24kb) and for some reason it seems to split into one large
    # alignment and one small one and outputs both
    if($processed_transcripts->{$transcript->stable_id()}) {
      if(($transcript->{'cov'} + $transcript->{'perc_id'}) > ($processed_transcripts->{$transcript->stable_id()}->{'cov'} + $processed_transcripts->{$transcript->stable_id()}->{'perc_id'})) {
        say "Found two transcripts for the same dbID. Selecting transcript with highest combined coverage and identity: ".$transcript->{'cov'}." cov, ".$transcript->{'perc_id'}." perc_id";
        $processed_transcripts->{$transcript->stable_id()} = $transcript;
      }
    } else {
      $processed_transcripts->{$transcript->stable_id()} = $transcript;
    }
  }

  foreach my $transcript_id (keys(%{$processed_transcripts})) {
    my $transcript = $processed_transcripts->{$transcript_id};
    my $biotype_group = $transcript->{'source_biotype_group'};
    my $coverage_cutoff = $coverage_cutoff_groups->{$biotype_group};
    my $perc_id_cutoff = $percent_identity_groups->{$biotype_group};

    unless($coverage_cutoff and $perc_id_cutoff) {
      $self->throw("Issue fetching coverage and percent id cutoffs for the biotype group of the parent transcript. Biotype group: ".$biotype_group);
    }

    if($transcript->{'cov'} >= $coverage_cutoff and $transcript->{'perc_id'} >= $perc_id_cutoff and
       $transcript->{'cds_length_diff'} <= $cds_length_diff_cutoff and $transcript->{'transcript_genomic_span_diff'} >= $genomic_span_diff_cutoff) {
      say "Transcript ".$transcript->{'source_stable_id'}." (".$transcript->stable_id().", ".$biotype_group.") passed check: ".$transcript->{'cov'}." cov, ".
          $transcript->{'perc_id'}." perc_id, ".$transcript->{'cds_length_diff'}." length diff, ".$transcript->{'transcript_genomic_span_diff'}." genomic span diff";
      push(@$good_transcripts,$transcript);
    } else {
      say "Transcript ".$transcript->{'source_stable_id'}." (".$transcript->stable_id().", ".$biotype_group.") failed check: ".$transcript->{'cov'}." cov, ".
          $transcript->{'perc_id'}." perc_id, ".$transcript->{'cds_length_diff'}." length diff, ".$transcript->{'transcript_genomic_span_diff'}." genomic span diff";
      push(@$bad_transcripts,$transcript);
    }
  }

  # Add any missing transcripts to the bad pile also


  return($bad_transcripts);
}


sub select_best_transcripts {
  my ($self,$bad_minimap_transcripts_hash,$bad_exonerate_transcripts_hash) = @_;

  # This will take two hashes of transcripts and check the keys, which are dbIDs to see which
  # hash has the highest combined coverage and perc identity in cases where there's a dbID across both
  # There probably shouldn't be dbIDs unique to each, but to make the code more robust it will handle
  # those cases too

  my $processed_ids_hash = {};
  my $selected_transcripts = [];

  foreach my $db_id (keys(%{$bad_minimap_transcripts_hash})) {
    my $minimap_transcript = $bad_minimap_transcripts_hash->{$db_id};
    my $exonerate_transcript = $bad_exonerate_transcripts_hash->{$db_id};

    # If there's no exonerate transcript just move on
    unless($exonerate_transcript) {
      push(@$selected_transcripts,$minimap_transcript);
      next;
    }

    # Just take whatever has the best combined perc id/cov for now
    if(($minimap_transcript->{'cov'} + $minimap_transcript->{'perc_id'}) > ($exonerate_transcript->{'cov'} + $exonerate_transcript->{'perc_id'})) {
      push(@$selected_transcripts,$minimap_transcript);
    } else {
      push(@$selected_transcripts,$exonerate_transcript);
    }

    $processed_ids_hash->{$db_id} = 1;
  }

  # Find and deal with any unique dbIDs in the exonerate hash
  foreach my $db_id (keys(%{$bad_exonerate_transcripts_hash})) {
    if($processed_ids_hash->{$db_id}) {
      next;
    }

    # In this scenario the only possiblity is that the transcript isn't in the minimap hash, so include it in the selected transcripts array
    my $exonerate_transcript = $bad_exonerate_transcripts_hash->{$db_id};
    push(@$selected_transcripts,$exonerate_transcript);
  }

  return($selected_transcripts);
}


sub list_missing_transcripts {
  my ($self,$transcripts_by_id,$source_transcripts) = @_;

  # This just gets the dbIDs from both hashes and then looks at the dbIDs of the source transcripts to see any ones that are not present
  my $missing_transcripts = [];

  my $db_id_hash = {};
  foreach my $db_id (keys(%{$transcripts_by_id})) {
    $db_id_hash->{$db_id} = 1;
  }

  foreach my $source_transcript (@$source_transcripts) {
    my $db_id = $source_transcript->dbID();
    unless($db_id_hash->{$db_id}) {
      say "Add source transcript ".$source_transcript->stable_id()." to missing transcripts list";
      push(@$missing_transcripts,$source_transcript);
    }
  }
  return($missing_transcripts);
}


sub label_transcript_status {
  my ($self,$transcripts_by_id,$status) = @_;

  my $coverage_threshold = 95;
  my $percent_id_threshold = 95;
  my $labelled_transcripts = [];
  # This labels a set of transcripts to help with deciding what to do when examining clusters
  foreach my $id (keys(%{$transcripts_by_id})) {
    my $transcript = $transcripts_by_id->{$id};
    my $transcript_coverage = $transcript->{'cov'};
    my $transcript_perc_id = $transcript->{'perc_id'};
    if ($transcript_coverage >= $coverage_threshold and $transcript_perc_id >= $percent_id_threshold) {
      $transcript->{'status'} = 'good';
    } else {
      $transcript->{'status'} = 'bad';
    }
    push(@$labelled_transcripts,$transcript);
  }

  return($labelled_transcripts);
}


sub generate_biotypes_hash {
  my ($self,$transcripts) = @_;

  my $unique_biotypes;
  my $biotypes_array = [];
  my $biotypes_hash = {};
  # This labels a set of transcripts to help with deciding what to do when examining clusters
  foreach my $transcript (@$transcripts) {
    unless($unique_biotypes->{$transcript->biotype()}) {
      push(@$biotypes_array,$transcript->biotype());
      $unique_biotypes->{$transcript->biotype()} = 1;
    }
  }

  $biotypes_hash->{'genes'} = $biotypes_array;
  return($biotypes_hash);
}


sub generate_single_transcript_genes {
  my ($self,$transcripts) = @_;

  my $single_transcript_genes = [];
  foreach my $transcript (@$transcripts) {
    my $gene = Bio::EnsEMBL::Gene->new(-analysis => $self->analysis);
    $gene->add_Transcript($transcript);
    $gene->biotype($transcript->biotype());
    $gene->stable_id($transcript->stable_id());
    $gene->{'status'} = $transcript->{'status'};
    $gene->{'cov'} = $transcript->{'cov'};
    $gene->{'perc_id'} = $transcript->{'perc_id'};
    $gene->slice($transcript->slice());
    push(@$single_transcript_genes,$gene);
  }

  return($single_transcript_genes);
}


sub check_cluster_status {
  my ($self,$cluster) = @_;

  my $genes = $cluster->get_Genes();
  $cluster->{'status'} = 'bad';
  foreach my $gene (@$genes) {
    if($gene->{'status'} eq 'good') {
      $cluster->{'status'} = 'good';
    }
  }
}


sub create_gene_from_cluster {
  my ($self,$cluster,$parent_gene_ids,$source_transcript_id_hash) = @_;

  my $cluster_genes = $cluster->get_Genes();
  my $processed_transcript_ids = {};
  my $selected_transcripts = {};

  # Loop through the transcripts, if there's any transcript that occurs twice, pick the one with the best combined coverage and percent id
  foreach my $cluster_gene (@$cluster_genes) {
    my $transcripts = $cluster_gene->get_all_Transcripts();
    my $transcript = ${$transcripts}[0];
    if($selected_transcripts->{$transcript->stable_id()}) {
      my $current_selected_transcript = $selected_transcripts->{$transcript->stable_id()};
      if(($transcript->{'cov'} + $transcript->{'perc_id'}) > ($current_selected_transcript->{'cov'} + $current_selected_transcript->{'perc_id'})) {
        $selected_transcripts->{$transcript->stable_id()} = $transcript;
      }
    } else {
      $selected_transcripts->{$transcript->stable_id()} = $transcript;
    }
  } # End foreach my $cluster_gene

  my $final_transcripts = [];
  foreach my $transcript_id (keys(%{$selected_transcripts})) {
    my $transcript = $selected_transcripts->{$transcript_id};
    push(@$final_transcripts,$transcript);
  }

  my $gene = Bio::EnsEMBL::Gene->new(-analysis => $self->analysis);
  my $transcript_id = ${$final_transcripts}[0]->stable_id();
  my $parent_gene_id = $parent_gene_ids->{$transcript_id}->{'gene_id'};
  my $parent_gene_stable_id = $parent_gene_ids->{$transcript_id}->{'gene_stable_id'};
  my $parent_gene_version = $parent_gene_ids->{$transcript_id}->{'gene_version'};
  my $parent_gene_biotype = $parent_gene_ids->{$transcript_id}->{'gene_biotype'};
  my $parent_gene_description = $parent_gene_ids->{$transcript_id}->{'gene_description'};

  foreach my $transcript (@$final_transcripts) {
    my $source_transcript = $source_transcript_id_hash->{$transcript->stable_id()};
    unless($source_transcript) {
      $self->throw("Issue with finding source transcript. The following dbID was not found in the source transcript id hash: ".$transcript->stable_id());
    }

    my $parent_transcript_stable_id = $parent_gene_ids->{$transcript->stable_id()}->{'transcript_stable_id'};
    my $parent_transcript_version = $parent_gene_ids->{$transcript->stable_id()}->{'transcript_version'};
    my $parent_gene_id = $parent_gene_ids->{$transcript->stable_id()}->{'gene_id'};
    my $is_canonical = $parent_gene_ids->{$transcript->stable_id()}->{'is_canonical'};
    my $source = $parent_gene_ids->{$transcript->stable_id()}->{'source'};

    $transcript->biotype($source_transcript->biotype());
    $transcript->stable_id($parent_transcript_stable_id);
    $transcript->version($parent_transcript_version);

    if($is_canonical) {
      $transcript->is_canonical(1);
    }

    $transcript->source($source);
    $gene->add_Transcript($transcript);
  }

  $gene->{'parent_gene_id'} = $parent_gene_id;
  $gene->stable_id($parent_gene_stable_id);
  $gene->version($parent_gene_version);
  $gene->biotype($parent_gene_biotype);
  #my $gene_description = ";parent_gene=".$parent_gene_stable_id.".".$parent_gene_version.";mapping_type=primary_mapping";
  #$gene->description($gene_description);
  $gene->description($parent_gene_description);

  # add source gene stable id as gene attribute
  my $parent_attribute = Bio::EnsEMBL::Attribute->new(-CODE => 'proj_parent_g',-VALUE => $parent_gene_stable_id.".".$parent_gene_version);
  $gene->add_Attributes($parent_attribute);

  return($gene);
}


sub access_transcripts_for_recovery {
  my ($self,$bad_transcripts,$source_transcript_id_hash,$good_transcripts,$target_parent_slice,$target_slice_adaptor,$target_sequence_adaptor,$target_genomic_name) = @_;

  # This will look at the bad transcripts and decide if there are any good transcripts from the same gene
  # If there are it will determine whether the bad transcript is contained within the boundaries covered
  # by the good transcript in the origninal source gene. If it is then the implication is that the region
  # if likely correct and there's just some fundamental issue with the transcript, so nothing can be done
  # If the transcript only partially overlapped with the good transcripts in the source gene, there's the
  # possibility that the target region was truncated and thus the transcript could only partially align to
  # it. If this is the case, then we want to attempt to create a new extended target region and try again

  my ($good_transcripts_start,$good_transcripts_end) = $self->get_transcript_boundaries($good_transcripts,$source_transcript_id_hash);
  my ($bad_transcripts_start,$bad_transcripts_end) = $self->get_transcript_boundaries($bad_transcripts,$source_transcript_id_hash);

  # Check if the bad transcript boundaries lie outside the good transcript boundaries in the source gene
  my $start_offset = 0;
  my $end_offset = 0;
  my $offset_buffer = 1000;
  if($bad_transcripts_start < $good_transcripts_start) {
    $start_offset = $good_transcripts_start - $bad_transcripts_start;
  }

  if($bad_transcripts_end > $good_transcripts_end) {
    $end_offset = $bad_transcripts_end - $good_transcripts_end;
  }

  my $target_transcripts_start;
  my $target_transcripts_end;
  foreach my $good_transcript (@$good_transcripts) {
    if(!$target_transcripts_start or $good_transcript->seq_region_start() < $target_transcripts_start) {
      $target_transcripts_start = $good_transcript->seq_region_start();
    }

    if(!$target_transcripts_end or $good_transcript->seq_region_end() > $target_transcripts_end) {
      $target_transcripts_end = $good_transcript->seq_region_end();
    }
  } # End foreach my $good_transcript

  $target_transcripts_start -= $start_offset - $offset_buffer;
  $target_transcripts_end += $end_offset + $offset_buffer;
  if($target_transcripts_start < 1) {
    $target_transcripts_start = 1;
  }


  my $target_strand = 1;
  my $target_genomic_seq = ${ $target_sequence_adaptor->fetch_by_Slice_start_end_strand($target_parent_slice, $target_transcripts_start, $target_transcripts_end, $target_strand) };
  my $target_region_slice = $target_slice_adaptor->fetch_by_region('toplevel', $target_genomic_name, $target_transcripts_start, $target_transcripts_end, $target_strand);
  my $target_genomic_fasta = ">".$target_genomic_name."\n".$target_genomic_seq;
  my $target_genome_file = $self->write_input_file([$target_genomic_fasta]);
  my $target_genome_index = $target_genome_file.".mmi";
  my $target_index_command = $self->program()." -d ".$target_genome_index." ".$target_genome_file;
  my $index_result = system($target_index_command);
  if($index_result) {
    $self->throw('The minimap2 index command returned a non-zero exit code. Commandline used:\n'.$target_index_command);
  }
}


sub calculate_max_intron_size {
  my ($self,$transcripts) = @_;

  my $scaling_factor = 1.5;
  my $max_intron_size = 100000;
  foreach my $transcript (@$transcripts) {
    my $introns = $transcript->get_all_Introns();
    foreach my $intron (@$introns) {
      my $scaled_intron_length = $scaling_factor * $intron->length();
      if($scaled_intron_length > $max_intron_size) {
        $max_intron_size = $scaled_intron_length;
      }
    }
  }

  return($max_intron_size);
}


sub get_transcript_boundaries {
  my ($self,$transcripts,$source_transcript_id_hash) = @_;

  my $transcripts_start;
  my $transcripts_end;

  foreach my $transcript (@$transcripts) {
    my $source_transcript = $source_transcript_id_hash->{$transcript->stable_id()};
    if(!$transcripts_start or $source_transcript->seq_region_start() < $transcripts_start) {
      $transcripts_start = $source_transcript->seq_region_start();
    }

    if(!$transcripts_end or $source_transcript->seq_region_end() > $transcripts_end) {
      $transcripts_end = $source_transcript->seq_region_end();
    }
  } # End foreach my $good_transcript

  return($transcripts_start,$transcripts_end);
}


sub create_annotation_features {
  my ($self,$transcript) = @_;

  my $cds_start  = $transcript->cdna_coding_start;
  my $cds_end    = $transcript->cdna_coding_end;
  my $stable_id  = $transcript->stable_id.".".$transcript->version;


  my $start_phase = $transcript->translation->start_Exon->phase();
  my $end_phase = $transcript->translation->end_Exon->end_phase();

  my $annotation_feature = Bio::EnsEMBL::Feature->new(-seqname => $stable_id,
                                                      -strand  => 1,
                                                      -start   => $cds_start,
                                                      -end     => $cds_end);

 my $annotation_features->{$stable_id} = $annotation_feature;
 return($annotation_features);
}


sub create_filename{
  my ($self, $stem, $ext, $dir, $no_clean) = @_;
  return create_file_name($stem, $ext, $dir, $no_clean);
}



sub genome_index {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_genome_index} = $val;
  }

  return $self->{_genome_index};
}


sub genes_to_process {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_genes_to_process} = $val;
  }

  return $self->{_genes_to_process};
}


sub input_file {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_input_file} = $val;
  }

  return $self->{_input_file};
}


sub paftools_path {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_paftools_path} = $val;
  }

  return $self->{_paftools_path};
}


sub source_adaptor {
  my ($self, $val) = @_;

  if (defined $val) {
    throw(ref($val).' is not a Bio::EnsEMBL::DBSQL::DBAdaptor')
      unless ($val->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'));
    $self->{_source_adaptor} = $val;
  }

  return $self->{_source_adaptor};
}


sub target_adaptor {
  my ($self, $val) = @_;

  if (defined $val) {
    throw(ref($val).' is not a Bio::EnsEMBL::DBSQL::DBAdaptor')
      unless ($val->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'));
    $self->{_target_adaptor} = $val;
  }

  return $self->{_target_adaptor};
}


sub delete_input_file {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_delete_input_file} = $val;
  }

  return $self->{_delete_input_file};
}


sub parent_gene_ids {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_parent_gene_ids} = $val;
  }

  return $self->{_parent_gene_ids};
}


sub gene_synteny_hash {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_gene_synteny_hash} = $val;
  }

  return $self->{_gene_synteny_hash};
}


sub gene_genomic_seqs_hash {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_gene_genomic_seqs_hash} = $val;
  }

  return $self->{_gene_genomic_seqs_hash};
}


sub no_projection {
  my ($self, $val) = @_;

  if ($val) {
    $self->{_no_projection} = $val;
  }

  return $self->{_no_projection};
}


=head2 exon_rank
  Arg [1]    : Bio::EnsEMBL::Exon to get the rank from
  Arg [2]    : Bio::EnsEMBL::Transcript where the exon in Arg [1] is found
  Description: It searches for the Arg [1] exon in the Arg [2] transcript exon array by using the exon coordinates.
               Note it is assumed that a transcript cannot have overlapping exons nor exons having the same coordinates.
  Returntype : int
               It returns the index (1 <= index <= number_of_exons) of the Arg [1] exon in the Arg [2] transcript exon array.
  Exceptions : It throws if the Arg [1] exon cannot be found in the Arg [2] transcript.
=cut

sub exon_rank {
  my ($self,$exon,$transcript) = @_;

  my $rank = 0;
  foreach my $t_exon (@{$transcript->get_all_Exons()}) {
    $rank++;
    if ($t_exon->seq_region_start() == $exon->seq_region_start() and
        $t_exon->seq_region_end() == $exon->seq_region_end() and
        $t_exon->seq_region_strand() == $exon->seq_region_strand()) {
      last;
    }
  }

  if ($rank) {
    return $rank;
  } else {
    $self->throw("Exon(seq_region_start,seq_region_end,seq_region_strand) ".$exon->seq_region_start()." ".$exon->seq_region_end()." ".$exon->seq_region_strand()." does not belong to transcript(seq_region_start,seq_region_end,seq_region_strand) ".$transcript->seq_region_start()." ".$transcript->seq_region_end()." ".$transcript->seq_region_strand());
  }
}

=head2 check_and_fix_translation_boundaries
  Arg [1]    : Bio::EnsEMBL::Transcript whose translation needs to be fixed
  Description: Fix the end of the translation if it's set to be beyond the end of the end exon
               or if it's set to be before the start of the end exon.
  Returntype : N/A
  Exceptions : N/A
=cut

sub check_and_fix_translation_boundaries {
  my ($self,$transcript) = @_;

  my $translation = $transcript->translation();
  my $end_exon = $translation->end_Exon(); # 'end_exon_id' exon in 'translation' table
  my $end_exon_length = $end_exon->length();

  if ($translation->end() > $end_exon_length) {
    print STDERR "Fixing the end of the translation (seq_start,seq_end,start_Exon->seq_region_start,end_Exon->seq_region_end) - (".$translation->start().",".$translation->end().",".$translation->start_Exon()->seq_region_start().",".$translation->end_Exon()->seq_region_end()." from ".$translation->end()." to ".$end_exon_length." because it is beyond the end of the end exon. Setting it to the maximum length of the end exon.\n";
    $translation->end($end_exon_length);
    $transcript->translation($translation);
  }

  $translation = $transcript->translation();
  my $end_exon_rank = $self->exon_rank($end_exon,$transcript);
  if ($translation->end() <= 0) {
    print STDERR "Transcript(seq_region_start,seq_region_end) ".$transcript->seq_region_start().",".$transcript->seq_region_end()." has translation end <= 0, number of exons = ".scalar(@{$transcript->get_all_Exons()}."\n"),
    print STDERR "Fixing the end of the translation (seq_start,seq_end,start_Exon->seq_region_start,end_Exon->seq_region_end) - (".$translation->start().",".$translation->end().",".$translation->start_Exon()->seq_region_start().",".$translation->end_Exon()->seq_region_end()." because it is set to 0. Set it to the end of the previous exon. End exon rank = ".$end_exon_rank."\n";
    foreach my $exon (@{$transcript->get_all_Exons()}) {
      my $exon_rank = $self->exon_rank($exon,$transcript);
      print STDERR "Translation end <= 0, foreach my exon (seq_region_start,seq_region_end,seq_region_strand,rank): ".$exon->seq_region_start().",".$exon->seq_region_end().",".$exon->seq_region_strand().",".$exon_rank."\n";
      if ($exon_rank == $end_exon_rank-1) {
        print STDERR "end_exon_rank-1 found.\n";
        $translation->end_Exon($exon);
        $translation->end($exon->length());
        print STDERR "End of the previous exon is ".$exon->length()."\n";
        $transcript->translation($translation);
        last;
      }
    }
  }
}

sub revcomp {
  my ($self,$seq) = @_;
  my $revcomp = reverse $seq;
  $revcomp =~ tr/ATGCatgc/TACGtacg/;
  return($revcomp);
}

1;
