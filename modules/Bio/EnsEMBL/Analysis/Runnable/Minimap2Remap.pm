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
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(compute_translation);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(set_alignment_supporting_features attach_Slice_to_Transcript attach_Analysis_to_Transcript calculate_exon_phases replace_stops_with_introns);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name align_nucleotide_seqs map_cds_location align_proteins);
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
  my ($genome_index, $input_file, $paftools_path, $source_adaptor, $target_adaptor, $delete_input_file, $parent_genes, $parent_gene_ids, $gene_synteny_hash, $gene_genomic_seqs_hash) = rearrange([qw (GENOME_INDEX INPUT_FILE PAFTOOLS_PATH SOURCE_ADAPTOR TARGET_ADAPTOR DELETE_INPUT_FILE PARENT_GENES PARENT_GENE_IDS GENE_SYNTENY_HASH GENE_GENOMIC_SEQS_HASH)],@args);
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
  if(system($minimap2_command)) {
    $self->throw("Error running minimap2\nError code: $?\nCommand line used:\n".$minimap2_command);
  }

  my $paf_results = [];
  open(IN,$paf_file);
  while(<IN>) {
    chomp($_);
    push(@$paf_results,$_);
  }
  close IN;

  my $high_confidence = 0;
  my $total_results = 0;
  my $processed_gene_ids = {};
  my $paf_results_hash = {};
  foreach my $paf_result (@$paf_results) {
    say "PAF results for first pass:\n".$paf_result;
    my @result_cols = split("\t",$paf_result);
    my $gene_id = $result_cols[0];
    if($paf_results_hash->{$gene_id}) {
      next;
    }
    $paf_results_hash->{$gene_id} = \@result_cols;
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

  my $paf_result = $paf_results->{$source_gene->dbID()};
  my $good_transcripts = []; # Transcripts that pass the cut-off thresholds
  my $bad_source_transcripts = []; # The source transcripts for mapped Transcripts that don't pass the threshold
  my $best_bad_transcripts = []; # When both the minimap and exonerate mappings fail the cut-offs, this will store the version with the highest combined coverage and percent id

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
  my $best_bad_transcripts_hash = {};
  if($paf_result) {
    my $source_gene_id = ${$paf_result}[0];
    my $source_genomic_length = ${$paf_result}[1];
    my $source_genomic_start = ${$paf_result}[2];
    my $source_genomic_end = ${$paf_result}[3];
    my $target_strand = ${$paf_result}[4];
    if($target_strand eq '+') {
       $target_strand = 1;
    } else {
      $target_strand = -1;
    }

    my $target_genomic_name = ${$paf_result}[5];
    my $target_genomic_length = ${$paf_result}[6];
    my $target_genomic_start = ${$paf_result}[7];
    my $target_genomic_end = ${$paf_result}[8];
    my $matching_bases = ${$paf_result}[9];
    my $total_bases = ${$paf_result}[10];
    my $mapping_quality = ${$paf_result}[11];

    my $mapping_identity = ($matching_bases/$total_bases) * 100;
    my $mapping_coverage = ($target_genomic_end - $target_genomic_start + 1)/$source_genomic_length;
    if($mapping_identity >= 80 and $mapping_coverage >= 0.8) {
      $high_confidence++;
    }

    my $adjust_left = $source_genomic_start;
    my $adjust_right = $source_genomic_length - $source_genomic_end;
    $target_genomic_start -= $adjust_left;
    $target_genomic_end += $adjust_right;
    if($target_genomic_start <= 0) {
      $target_genomic_start = 1;
    }

    if($target_genomic_end > $target_genomic_length) {
      $target_genomic_end = $target_genomic_length;
    }

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

    # At this point make a decision as to how to proceed with the gene. If it's small, single exon gene we need to treat it differently to a longer one
    if($self->check_for_small_gene($source_transcripts)) {
      # Process the short gene through direct projection
      $self->process_small_gene($source_gene,$source_transcripts,$target_genomic_start,$target_region_slice,$target_strand,
                                $good_transcripts,$best_bad_transcripts,$good_transcripts_hash,$best_bad_transcripts_hash);
    } else {
      my $target_genome_index = $target_genome_file.".mmi";
      my $target_index_command = $self->program()." -d ".$target_genome_index." ".$target_genome_file;
      my $index_result = system($target_index_command);
      if($index_result) {
        $self->throw('The minimap2 index command returned a non-zero exit code. Commandline used:\n'.$target_index_command);
      }

      say "Processing a total of ".scalar(@$source_transcripts)." source transcripts";

      my $source_transcript_fasta_seqs = [];

      foreach my $source_transcript (@$source_transcripts) {
        say "Writing ".$source_transcript->stable_id()." to file for mapping";
        my $source_transcript_sequence = $source_transcript->seq->seq();
        my $fasta_record = ">".$source_transcript->dbID()."\n".$source_transcript_sequence;
        push(@$source_transcript_fasta_seqs,$fasta_record);
      }

      say "Generating initial set of minimap2 mappings";
      my $output_minimap_transcripts = $self->generate_minimap_transcripts($source_transcript_fasta_seqs,$target_genome_index,$target_adaptor,$target_genomic_start,$max_intron_size);
      my $bad_minimap_transcripts = $self->check_mapping_quality($output_minimap_transcripts,$source_transcript_id_hash,$good_transcripts);
      say "Number of good transcripts after minimap2: ".scalar(@$good_transcripts);
      say "Number of bad transcripts after minimap2: ".scalar(@$bad_minimap_transcripts);

      say "Running exonerate on bad transcripts";
      my $output_exonerate_transcripts = $self->generate_exonerate_transcripts($bad_minimap_transcripts,$source_transcript_id_hash,$target_region_slice,$target_slice_adaptor,$max_intron_size);
      my $bad_exonerate_transcripts = $self->check_mapping_quality($output_exonerate_transcripts,$source_transcript_id_hash,$good_transcripts);
      say "Number of good transcripts after exonerate: ".scalar(@$good_transcripts);
      say "Number of bad transcripts after exonerate: ".scalar(@$bad_exonerate_transcripts);

      # Store these in a hash for making selecting the best based on dbID easier
      foreach my $transcript (@$good_transcripts) {
        $good_transcripts_hash->{$transcript->stable_id()} = $transcript;
      }

      my $bad_minimap_transcripts_hash = {};
      foreach my $transcript (@$bad_minimap_transcripts) {
        $bad_minimap_transcripts_hash->{$transcript->stable_id()} = $transcript;
      }

      my $bad_exonerate_transcripts_hash = {};
      foreach my $transcript (@$bad_exonerate_transcripts) {
        $bad_exonerate_transcripts_hash->{$transcript->stable_id()} = $transcript;
      }

      # This is bit of a mess as I only realised that many of the hardest edge cases can be recovered with exonerate, but might totally fail on minimap2
      # Thus there needed to be a step added for when the minimap transcript is missing to begin wit, because the way the rest of the code is designed
      # a transcript might just not be found via minimap and then all that will happen is that it gets globally aligned via minimap later and likely not
      # found again. This just ensures that a local exonerate is run on the region for totally missing transcripts, not just ones that fail the minimap
      # cut-offs. To refactor this at some point it would be better to rely less on arrays of the mapped transcripts and instead go through the id hashes
      say "Checking for missing transcripts to perform exonerate on region";
      my $initial_missing_source_transcripts = $self->list_missing_transcripts($good_transcripts_hash,$bad_exonerate_transcripts_hash,$source_transcripts);
      say "Found ".scalar(@$initial_missing_source_transcripts)." initial missing transcripts";
      my $missing_exonerate_transcripts = [];
      foreach my $missing_source_transcript (@$initial_missing_source_transcripts) {
        say "Initial missing source transcript: ".$missing_source_transcript->stable_id();
        my $exonerate_transcripts = $self->run_exonerate($missing_source_transcript,$target_region_slice,$target_slice_adaptor,$max_intron_size);
        push(@$missing_exonerate_transcripts,@$exonerate_transcripts);
      }

      say "After mapping initial missing transcripts with exonerate in region found ".scalar(@$missing_exonerate_transcripts)." exonerate transcripts";

      my $bad_initial_missing_transcripts = $self->check_mapping_quality($missing_exonerate_transcripts,$source_transcript_id_hash,$good_transcripts);
      foreach my $transcript (@$bad_initial_missing_transcripts) {
        $bad_exonerate_transcripts_hash->{$transcript->stable_id()} = $transcript;
      }

      # Select the best transcript foreach pair across both hashes, and any transcripts that are unique to either hash
      # Then put the corresponding source transcripts into an array to do a global mapping later
      say "Selecting best transcripts out of bad minimap2/exonerate transcripts";
      $best_bad_transcripts = $self->select_best_transcripts($bad_minimap_transcripts_hash,$bad_exonerate_transcripts_hash);
      foreach my $transcript (@$best_bad_transcripts) {
        # If a transcript was bad from minimap and good after exonerate, then we want to skip over it
        if($good_transcripts_hash->{$transcript->stable_id()}) {
          next;
        }

        $best_bad_transcripts_hash->{$transcript->stable_id()} = $transcript;
        push(@$bad_source_transcripts,$source_transcript_id_hash->{$transcript->stable_id()});
      }

      say "Number of bad transcripts after selection: ".scalar(@$bad_source_transcripts);
    } # End else
  } else {
    say "Did not get a paf result from the first pass, will resort to global mapping of the transcripts";
  }

  # At this point we have the selected set of good/bad transcripts. There may also be missing transcripts at the moment, i.e.
  # those that did no get aligned by either method. We now want to take the original sequences for the bad and missing transcripts
  # and run a genome-wide alignment via minimap. We then have two scenarios to work out. The first is if there were any good transcripts
  # to begin with. If there were then were going to make the assumption (correctly or incorrectly), that the good transcripts represent
  # where the gene is likely to be. We need to decide once we do the global mapping what to do with any mapped transcripts. The first thing
  # is to check whether a global transcript actually scores better than the original mapping (if there was one). If not, then we just add
  # the selected bad transcript to the gene. If it is better then we have to work out what to do. We basically need to cluster the transcripts
  # and figure out if they're all in the same region.

  say "Checking for missing transcripts";
  my $missing_source_transcripts = $self->list_missing_transcripts($good_transcripts_hash,$best_bad_transcripts_hash,$source_transcripts);
  say "Found ".scalar(@$missing_source_transcripts)." missing transcripts";

  # Now run a global mapping of the missing/bad transcripts
  my $transcripts_to_map = [@$missing_source_transcripts,@$bad_source_transcripts];
  my $transcripts_to_map_fasta_seqs = [];
  foreach my $transcript (@$transcripts_to_map) {
    my $transcript_sequence = $transcript->seq->seq();
    my $fasta_record = ">".$transcript->dbID()."\n".$transcript_sequence;
    push(@$transcripts_to_map_fasta_seqs,$fasta_record);
  }

  say "Preparing to map bad and missing transcripts globally to the genome";
  say "Number of transcripts to map globally: ".scalar(@$transcripts_to_map);
  my $target_global_genome_index = $self->genome_index();
  my $global_mapped_transcripts = $self->generate_minimap_transcripts($transcripts_to_map_fasta_seqs,$target_global_genome_index,$target_adaptor,1,$max_intron_size);
  say "Number of globally mapped transcripts: ".scalar(@$global_mapped_transcripts);

  my $global_mapped_transcripts_hash = {};
  foreach my $transcript (@$global_mapped_transcripts) {
    my $db_id = $transcript->stable_id();
    $global_mapped_transcripts_hash->{$db_id} = $transcript;
  }

  my $good_global_transcripts = [];
  my $bad_global_transcripts = $self->check_mapping_quality($global_mapped_transcripts,$source_transcript_id_hash,$good_global_transcripts);
  say "Number of good globally mapped transcripts: ".scalar(@$good_global_transcripts);
  say "Number of bad globally mapped transcripts: ".scalar(@$bad_global_transcripts);

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

  $self->label_transcript_status($good_transcripts,'good');
  $self->label_transcript_status($good_global_transcripts,'good');
  $self->label_transcript_status($best_bad_transcripts,'bad');
  $self->label_transcript_status($bad_global_transcripts,'bad');

  my $all_transcripts = [@$good_transcripts,@$good_global_transcripts,@$best_bad_transcripts,@$bad_global_transcripts];
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


sub process_small_gene {
  my ($self,$gene,$source_transcripts,$target_genomic_start,$target_region_slice,$target_strand,$good_transcripts,$best_bad_transcripts,$good_transcripts_hash,$best_bad_transcripts_hash) = @_;

  my $coverage_cutoff = 95;
  my $perc_id_cutoff = 95;

  say "Processing small gene: ".$gene->stable_id().", Biotype: ".$gene->biotype().", Length: ".$gene->length();
  # This will use the source gene seq originally calculated in the RunnableDB and the coords for it
  # Then want to align it to the target region and try to work out the coords of the start end of the feature in the alignment
  my $gene_genomic_seqs_hash = $self->gene_genomic_seqs_hash();
  my $source_genomic_seq_info = $gene_genomic_seqs_hash->{$gene->dbID()};

  unless($source_genomic_seq_info) {
    $self->throw("Could not find the genomic seq info for source gene with dbID: ".$gene->dbID());
  }

  my $source_region_start = ${$source_genomic_seq_info}[0];
  my $source_region_end = ${$source_genomic_seq_info}[1];
  my $source_genome_seq = ${$source_genomic_seq_info}[2];

  say "Unaligned seqs:";
  say $source_genome_seq."\n";
  say $target_region_slice->seq();

  my $target_genome_seq = $target_region_slice->seq();
  if($target_strand != 1) {
    my $revcomp = reverse $target_genome_seq;
    $revcomp =~ tr/ATGCatgc/TACGtacg/;
    $target_genome_seq = $revcomp;
  }

  my ($coverage,$percent_id,$aligned_source_seq,$aligned_target_seq) = align_nucleotide_seqs($source_genome_seq,$target_genome_seq);

  say "Aligned source and target regions with MAFFT: Coverage: ".$coverage.", Percent id: ".$percent_id;
  foreach my $transcript (@$source_transcripts) {
    # If there's a translation we just want the CDS exon for simplicity (though in these cases it's unlikely there would be any UTR)
    my $exons;
    $exons = $transcript->get_all_Exons();

    my $exon = ${$exons}[0];
    my $exon_region_start = $exon->seq_region_start();
    my $exon_region_end = $exon->seq_region_end();
    my $projected_exon = $self->project_feature($transcript,$exon,$source_region_start,$exon_region_start,$exon_region_end,$aligned_source_seq,$aligned_target_seq,$target_region_slice,$target_strand);

    if($projected_exon) {
      my $projected_transcript = Bio::EnsEMBL::Transcript->new(-EXONS => [$projected_exon]);
      $projected_transcript->slice($projected_exon->slice());
      $projected_transcript->analysis($self->analysis());
      $projected_transcript->biotype($transcript->biotype());
      $projected_transcript->stable_id($transcript->dbID());
      $projected_transcript->version($transcript->version());
      $projected_transcript->{'cov'} = $coverage;
      $projected_transcript->{'perc_id'} = $percent_id;
      $projected_transcript->{'source_stable_id'} = $transcript->stable_id();
      $projected_transcript->{'source_biotype_group'} = $transcript->get_Biotype->biotype_group();
      $projected_transcript->{'source_length'} = $transcript->length();

      ###################
      # TEST PRINT
#      say "Projected transcript seq:\n".$projected_transcript->seq->seq();
#      say "Source transcript seq:\n".$transcript->seq->seq();
#      say "Projected exon seq:\n".$projected_exon->seq->seq();
#      say "Source exon seq:\n".$projected_exon->seq->seq();

      my $transcript_seq;
      my $projected_transcript_seq;
      if($transcript->translation()) {
        say "Transcript has a translation, so will apply a CDS across the mapped region";
        $self->set_complete_transcript_cds($projected_transcript);
        $transcript_seq = $transcript->translateable_seq();
        $projected_transcript_seq = $projected_transcript->translateable_seq();
#        say "TEST TRANSLATION: !!!!!!!!";
#        say $transcript->translation->seq();
#        say "\n".$projected_transcript->translation->seq();
      } else {
        $transcript_seq = $transcript->seq->seq();
        $projected_transcript_seq = $projected_transcript->seq->seq();
      }

      my ($coverage,$percent_id,$aligned_source_seq,$aligned_target_seq) = align_nucleotide_seqs($transcript_seq,$projected_transcript_seq);
      say "POST align source seq:\n".$aligned_source_seq;
      say "POST align target seq:\n".$aligned_target_seq;

      my $description_string = "Parent: ".$transcript->stable_id().".".$transcript->version().", Coverage: ".$coverage.", Perc id: ".$percent_id;
      $projected_transcript->description($description_string);
      if($coverage >= $coverage_cutoff and $percent_id >= $perc_id_cutoff) {
        say "Projected transcript passed cut-offs: ".$coverage." cov, ".$percent_id." percent id";
        push(@$good_transcripts,$projected_transcript);
        $good_transcripts_hash->{$transcript->dbID()} = $projected_transcript;
      } else {
        say "Projected transcript failed cut-offs: ".$coverage." cov, ".$percent_id." percent id";
        push(@$best_bad_transcripts,$projected_transcript);
        $best_bad_transcripts_hash->{$transcript->dbID()} = $projected_transcript;
      }
#      say "TEST COVERAGE: !!!";
#      say $aligned_source_seq;
#      say "\n".$aligned_target_seq;

    } # End if($projected_exon)
  } # End foreach my $transcript
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
  say "Aligned seqs:";
  say $aligned_source_seq;
  say $aligned_target_seq;

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
  my $recovered_target_feature_seq = substr($target_seq,$target_seq_start_index,$target_feature_length);

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

#  say "Region slice: ".$region_slice->name();
  my $region_start = $region_slice->seq_region_start();
  my $region_end = $region_slice->seq_region_end();
#  say "Region start: ".$region_start;
#  say "Region end: ".$region_end;

  my $parent_slice = $region_slice->seq_region_Slice();
#  say "Parent slice: ".$parent_slice->name();

  my $exon_start;
  my $exon_end;

  if($target_strand == 1) {
    $exon_start = $region_start + $seq_start_index;
    $exon_end = $region_start + $seq_end_index;
  } else {
    # In this case we need to reverse the coords
    $exon_end = $region_end - $seq_start_index;
    $exon_start = $exon_end - ($seq_end_index - $seq_start_index);
  }


  my $phase = -1;
  my $end_phase = -1;

  my $projected_exon = Bio::EnsEMBL::Exon->new(-start     => $exon_start,
                                               -end       => $exon_end,
                                               -strand    => $target_strand,
                                               -phase     => $phase,
                                               -end_phase => $end_phase,
                                               -analysis  => $self->analysis,
                                               -slice     => $parent_slice);

  if($exon->is_coding($transcript)) {
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
    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts) {
      if ($transcript->translate() and $transcript->translate()->seq()) {
        my $transcript_after_replaced_stops = replace_stops_with_introns($transcript,$max_stops);
        if ($transcript_after_replaced_stops and $transcript_after_replaced_stops->translate()) {
          push(@$minimap_transcripts,$transcript_after_replaced_stops);
        } elsif (!$transcript_after_replaced_stops->translate()) {
          print STDERR "minimap transcript (seq_region_start,seq_region_end,seq_region_strand,seq_region_name) (".$transcript->seq_region_start().",".$transcript->seq_region_end().",".$transcript->seq_region->strand().",".$transcript->seq_region_name().") does not translate after replacing a maximum of $max_stops stops. Discarded.\n";
        } else {
          print STDERR "minimap transcript (seq_region_start,seq_region_end,seq_region_strand,seq_region_name) (".$transcript->seq_region_start().",".$transcript->seq_region_end().",".$transcript->seq_region->strand().",".$transcript->seq_region_name().") has more than the maximum of $max_stops stops. Discarded.\n";
        }
      } else {
        push(@$minimap_transcripts,$transcript);
      }
    }
  }

  return($minimap_transcripts);
}


sub generate_exonerate_transcripts {
  my ($self,$bad_target_transcripts,$source_transcript_id_hash,$target_region_slice,$target_slice_adaptor,$max_intron_size) = @_;

  my $max_stops = 999;
  my $output_transcripts = [];

  foreach my $bad_transcript (@$bad_target_transcripts) {
    my $source_transcript = $source_transcript_id_hash->{$bad_transcript->stable_id()};
    unless($source_transcript) {
      $self->throw("Couldn't find the dbID of the transcript in the source transcript hash when attempting to run exonerate");
    }

    my $exonerate_transcripts = $self->run_exonerate($source_transcript,$target_region_slice,$target_slice_adaptor,$max_intron_size);
    foreach my $transcript (@$exonerate_transcripts) {
      if ($transcript->translate() and $transcript->translate()->seq()) {
        my $transcript_after_replaced_stops = replace_stops_with_introns($transcript,$max_stops);
        if ($transcript_after_replaced_stops and $transcript_after_replaced_stops->translate()) {
          push(@$output_transcripts,$transcript_after_replaced_stops);
        } elsif (!$transcript_after_replaced_stops->translate()) {
          print STDERR "exonerate transcript (seq_region_start,seq_region_end,seq_region_strand,seq_region_name) (".$transcript->seq_region_start().",".$transcript->seq_region_end().",".$transcript->seq_region->strand().",".$transcript->seq_region_name().") does not translate after replacing a maximum of $max_stops stops. Discarded.\n";
        } else {
          print STDERR "exonerate transcript (seq_region_start,seq_region_end,seq_region_strand,seq_region_name) (".$transcript->seq_region_start().",".$transcript->seq_region_end().",".$transcript->seq_region->strand().",".$transcript->seq_region_name().") has more than the maximum of $max_stops stops. Discarded.\n";
        }
      } else {
        push(@$output_transcripts,$transcript);
      }
    }
    #push(@$output_transcripts,@$exonerate_transcripts);
  }

  return($output_transcripts);
}


sub run_exonerate {
  my ($self,$source_transcript,$target_region_slice,$target_slice_adaptor,$max_intron_size) = @_;

  my $output_transcripts = [];
  my $exonerate_length_cutoff = 15000;
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
  if(scalar(@{$runnable->output()})) {
     my $initial_output_transcript = ${$runnable->output()}[0];
     my $output_transcript = $self->update_exonerate_transcript_coords($initial_output_transcript,$target_slice_adaptor);
     if($source_transcript->translation() and $output_transcript->translation()) {
       $self->check_exonerate_translation($source_transcript,$output_transcript);
     }
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

  # Fix the end of the translation if it's set to be beyond the end of the end exon
  # or if it's set to be before the start of the end exon
  my $translation = $output_transcript->translation();
  my $end_exon = $translation->end_Exon(); # 'end_exon_id' exon in 'translation' table
  my $end_exon_length = $end_exon->seq_region_end()-$end_exon->seq_region_start()+1;

  if ($translation->end() > $end_exon_length) {
    print STDERR "Fixing the end of the translation (seq_start,seq_end,start_Exon->seq_region_start,end_Exon->seq_region_end) - (".$translation->start().",".$translation->end().",".$translation->start_Exon()->seq_region_start().",".$translation->end_Exon()->seq_region_end()." from ".$translation->end()." to ".$end_exon_length." because it is beyond the end of the end exon. Setting it to the maximum length of the end exon.\n";
    $translation->end($end_exon_length);
    $output_transcript->translation($translation);
  }

  $translation = $output_transcript->translation();;
  my $end_exon_rank = $end_exon->rank();
  if ($translation->end() <= 0) {
    print STDERR "Fixing the end of the translation (seq_start,seq_end,start_Exon->seq_region_start,end_Exon->seq_region_end) - (".$translation->start().",".$translation->end().",".$translation->start_Exon()->seq_region_start().",".$translation->end_Exon()->seq_region_end()." because it is set to 0. Set it to the end of the previous exon.\n";
    foreach my $exon (@{$output_transcript->get_all_constitutive_Exons()}) {
      if ($exon->rank() == $end_exon_rank-1) {
        $translation->end_Exon($exon);
        $translation->end($exon->length());
	print STDERR "End of the previous exon is ".$exon->length()."\n";
	$output_transcript->translation($translation);
      }
    }
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
    my $transcript_description = "Parent: ".$source_transcript->stable_id().".".$source_transcript->version();
    $transcript->description($transcript_description);
    if($source_transcript->translation()) {
      $source_transcript_seq = $source_transcript->translateable_seq();
      $transcript_seq = $transcript->translateable_seq();

      unless($transcript_seq) {
        compute_translation($transcript);
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

    $transcript_description .= ", Coverage: ".$coverage.", Perc id: ".$percent_id;
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
  my ($self,$good_transcripts_hash,$best_bad_transcripts_hash,$source_transcripts) = @_;

  # This just gets the dbIDs from both hashes and then looks at the dbIDs of the source transcripts to see any ones that are not present
  my $missing_transcripts = [];

  my $db_id_hash = {};
  foreach my $db_id (keys(%{$good_transcripts_hash})) {
    $db_id_hash->{$db_id} = 1;
  }

  foreach my $db_id (keys(%{$best_bad_transcripts_hash})) {
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
  my ($self,$transcripts,$status) = @_;

  # This labels a set of transcripts to help with deciding what to do when examining clusters
  foreach my $transcript (@$transcripts) {
    $transcript->{'status'} = $status;
  }
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
#    my $cov_string = "cov: ".$transcript->{'cov'}." perc_id: ".$transcript->{'perc_id'};
#    $transcript->description($cov_string);
    $gene->add_Transcript($transcript);
  }

  $gene->{'parent_gene_id'} = $parent_gene_id;
  $gene->stable_id($parent_gene_stable_id);
  $gene->version($parent_gene_version);
  $gene->biotype($parent_gene_biotype);
  my $gene_description = "Parent: ".$parent_gene_stable_id.".".$parent_gene_version.", Type: Primary mapping";
  $gene->description($gene_description);

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

1;
