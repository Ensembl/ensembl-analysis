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

use Bio::EnsEMBL::Analysis::Runnable::Minimap2;
use Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(compute_translation);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(set_alignment_supporting_features attach_Slice_to_Transcript attach_Analysis_to_Transcript);
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
  my ($genome_index, $input_file, $paftools_path, $source_adaptor, $target_adaptor, $delete_input_file, $parent_gene_ids, $gene_synteny_hash) = rearrange([qw (GENOME_INDEX INPUT_FILE PAFTOOLS_PATH SOURCE_ADAPTOR TARGET_ADAPTOR DELETE_INPUT_FILE PARENT_GENE_IDS GENE_SYNTENY_HASH)],@args);
  $self->genome_index($genome_index);
  $self->input_file($input_file);
  $self->paftools_path($paftools_path);
  $self->source_adaptor($source_adaptor);
  $self->target_adaptor($target_adaptor);
  $self->delete_input_file($delete_input_file);
  $self->parent_gene_ids($parent_gene_ids);
  $self->gene_synteny_hash($gene_synteny_hash);
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
  foreach my $paf_result (@$paf_results) {
    my @result_cols = split("\t",$paf_result);
    my $gene_id = $result_cols[0];
    if($processed_gene_ids->{$gene_id}) {
      next;
    }
    $high_confidence += $self->process_results(\@result_cols);
    $total_results++;
    $processed_gene_ids->{$gene_id} = 1;
  }
  say "TOTAL RESULTS: ".$total_results;
  say "HIGH CONFIDENCE: ".$high_confidence;
} # End run



sub process_results {
  my ($self, $paf_result) = @_;

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


  my $high_confidence = 0;
  my $source_gene_id = ${$paf_result}[0];
  my $source_genomic_length = ${$paf_result}[1];
  my $source_genomic_start = ${$paf_result}[2];
  my $source_genomic_end = ${$paf_result}[3];
  my $same_strand = ${$paf_result}[4];
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

  my $target_adaptor = $self->target_adaptor();
  my $target_slice_adaptor = $target_adaptor->get_SliceAdaptor();
  my $target_parent_slice = $target_slice_adaptor->fetch_by_region('toplevel',$target_genomic_name);

  unless($target_parent_slice) {
    $self->throw("Could not fetch the slice in the target assembly. Slice name: ".$target_genomic_name);
  }

  my $target_strand = 1;
  my $target_sequence_adaptor = $target_adaptor->get_SequenceAdaptor;
  my $target_genomic_seq = ${ $target_sequence_adaptor->fetch_by_Slice_start_end_strand($target_parent_slice, $target_genomic_start, $target_genomic_end, $target_strand) };
  my $target_region_slice = $target_slice_adaptor->fetch_by_region('toplevel', $target_genomic_name, $target_genomic_start, $target_genomic_end, $target_strand);
  my $target_genomic_fasta = ">".$target_genomic_name."\n".$target_genomic_seq;
  my $target_genome_file = $self->write_input_file([$target_genomic_fasta]);
  my $target_genome_index = $target_genome_file.".mmi";
  my $target_index_command = $self->program()." -d ".$target_genome_index." ".$target_genome_file;
  my $index_result = system($target_index_command);
  if($index_result) {
    $self->throw('The minimap2 index command returned a non-zero exit code. Commandline used:\n'.$target_index_command);
  }

  # Covert the transcripts into fasta records
  my $source_adaptor = $self->source_adaptor();
  my $source_gene_adaptor = $source_adaptor->get_GeneAdaptor();
  my $source_gene = $source_gene_adaptor->fetch_by_dbID($source_gene_id);
  my $source_transcripts = $source_gene->get_all_Transcripts();
  say "Processing a total of ".scalar(@$source_transcripts)." source transcripts";

  my $source_transcript_id_hash = {};
  my $source_transcript_fasta_seqs = [];


  foreach my $source_transcript (@$source_transcripts) {
    my $source_transcript_id = $source_transcript->dbID();
    $source_transcript_id_hash->{$source_transcript_id} = $source_transcript;
    say "STID HASH: ".$source_transcript_id;
    my $source_transcript_sequence = $source_transcript->seq->seq();
    my $fasta_record = ">".$source_transcript->dbID()."\n".$source_transcript_sequence;
    push(@$source_transcript_fasta_seqs,$fasta_record);
  }

  say "Generating initial set of minimap2 mappings";
  my $output_minimap_transcripts = $self->generate_minimap_transcripts($source_transcript_fasta_seqs,$target_genome_index,$target_adaptor,$target_genomic_start);
  my $good_transcripts = [];
  my $bad_minimap_transcripts = $self->check_mapping_quality($output_minimap_transcripts,$source_transcript_id_hash,$good_transcripts);
  say "Number of good transcripts after minimap2: ".scalar(@$good_transcripts);
  say "Number of bad transcripts after minimap2: ".scalar(@$bad_minimap_transcripts);

  say "Running exonerate on bad transcripts";
  my $output_exonerate_transcripts = $self->generate_exonerate_transcripts($bad_minimap_transcripts,$source_transcript_id_hash,$target_region_slice,$target_slice_adaptor);
  my $bad_exonerate_transcripts = $self->check_mapping_quality($output_exonerate_transcripts,$source_transcript_id_hash,$good_transcripts);
  say "Number of good transcripts after exonerate: ".scalar(@$good_transcripts);
  say "Number of bad transcripts after exonerate: ".scalar(@$bad_exonerate_transcripts);

  # Store these in a hash for making selecting the best based on dbID easier
  my $good_transcripts_hash = {};
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

  # Select the best transcript foreach pair across both hashes, and any transcripts that are unique to either hash
  # Then put the corresponding source transcripts into an array to do a global mapping later
  say "Selecting best transcripts out of bad minimap2/exonerate transcripts";
  my $best_bad_transcripts = $self->select_best_transcripts($bad_minimap_transcripts_hash,$bad_exonerate_transcripts_hash);
  my $best_bad_transcripts_hash = {};
  my $bad_source_transcripts = [];
  foreach my $transcript (@$best_bad_transcripts) {
    # If a transcript was bad from minimap and good after exonerate, then we want to skip over it
    if($good_transcripts_hash->{$transcript->stable_id()}) {
      next;
    }

    $best_bad_transcripts_hash->{$transcript->stable_id()} = $transcript;
    push(@$bad_source_transcripts,$source_transcript_id_hash->{$transcript->stable_id()});
  }

  say "Number of bad transcripts after selection: ".scalar(@$bad_source_transcripts);

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
  my $global_mapped_transcripts = $self->generate_minimap_transcripts($transcripts_to_map_fasta_seqs,$target_global_genome_index,$target_adaptor,1);
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

#  my $all_initial_transcripts = [@$good_transcripts,@$good_global_transcripts,@$best_bad_transcripts,@$bad_global_transcripts];
#  my $unique_transcripts_hash = {};
#  foreach my $transcript (@$all_initial_transcripts) {
#    unless($unique_transcripts_hash->{$transcript->stable_id()}) {
#      $unique_transcripts_hash->{$transcript->stable_id()} = $transcript;
#      next;
#    }

#    my $current_selected_transcript = $unique_transcripts_hash->{$transcript->stable_id()};
#    if(($transcript->{'cov'} + $transcript->{'perc_id'}) > ($current_selected_transcript->{'cov'} + $current_selected_transcript->{'perc_id'})) {
#	    $unique_transcripts_hash->{$transcript->stable_id()} = $transcript;
#    }
#  }

#  my $all_transcripts = [];
#  foreach my $transcript_id (keys(%{$unique_transcripts_hash})) {
#    my $transcript = $unique_transcripts_hash->{$transcript_id};
#    push(@$all_transcripts,$transcript);
#  }

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
  my $final_genes = [];
  foreach my $cluster (@$final_clusters) {
    my $gene = $self->create_gene_from_cluster($cluster,$parent_gene_ids,$source_transcript_id_hash);
    say "Created gene: ".$gene->stable_id()." ".$gene->seq_region_name().":".$gene->seq_region_start.":".$gene->seq_region_end.":".$gene->strand();
    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts) {
      say "  Transcript: ".$transcript->stable_id()." ".$transcript->seq_region_start.":".$transcript->seq_region_end.":".$transcript->strand();
#      if($transcript->stable_id() eq 'ENST00000608875') {
#        say $transcript->seq->seq();
#      }
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


sub generate_minimap_transcripts {
  my ($self,$source_transcript_fasta_seqs,$target_genome_index,$target_adaptor,$target_genomic_start) = @_;

  # This will take in the target genomic region (via the index)
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
  );

  $runnable->run();

  my $output_genes = $runnable->output();
  foreach my $gene (@$output_genes) {
    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts) {
      push(@$minimap_transcripts,$transcript);
    }
  }

  return($minimap_transcripts);
}


sub generate_exonerate_transcripts {
  my ($self,$bad_target_transcripts,$source_transcript_id_hash,$target_region_slice,$target_slice_adaptor) = @_;

  my $output_transcripts = [];


 # my $gene = Bio::EnsEMBL::Gene->new(-analysis => $self->analysis);

  my %parameters = ();
  $parameters{-options} = "--model cdna2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE --score 500 ".
                          "--saturatethreshold 100 --dnawordlen 15 --codonwordlen 15 --dnahspthreshold 60 --bestn 1";
  $parameters{-coverage_by_aligned} = 1;

  foreach my $bad_transcript (@$bad_target_transcripts) {
      my $source_transcript = $source_transcript_id_hash->{$bad_transcript->stable_id()};
      unless($source_transcript) {
        $self->throw("Couldn't find the dbID of the transcript in the source transcript hash when attempting to run exonerate");
      }

      say "Running Exonerate on ".$source_transcript->stable_id();
      my $annotation_features;
      if($source_transcript->translation()) {
        $annotation_features = $self->create_annotation_features($source_transcript);
      }

      my $source_transcript_header = $source_transcript->stable_id.'.'.$source_transcript->version;
      my $source_transcript_seq_object = Bio::Seq->new(-display_id => $source_transcript_header, -seq => $source_transcript->seq->seq());

      my $runnable = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript->new(
                       -program  => '/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/opt/exonerate22/bin/exonerate',
                       -analysis => $self->analysis(),
                       -query_type     => "dna",
                       -calculate_coverage_and_pid => 1,
                       -annotation_features => $annotation_features,
                        %parameters,
                     );
      $runnable->target_seqs([$target_region_slice]);
      $runnable->query_seqs([$source_transcript_seq_object]);
      $runnable->run();
      if(scalar(@{$runnable->output()})) {
        my $initial_output_transcript = ${$runnable->output()}[0];
        my $output_transcript = $self->update_exonerate_transcript_coords($initial_output_transcript,$target_slice_adaptor);
        $output_transcript->stable_id($source_transcript->dbID());
#        say "FERGAL TEST NEW SEQ: ".$output_transcript->stable_id;
#        say $output_transcript->seq->seq();

#        say "FERGAL TEST ORIG SEQ: ".$output_transcript->stable_id;
#        say $output_transcript->seq->seq();
#        $output_transcript->start($output_transcript->seq_region_start());
#        $output_transcript->end($output_transcript->seq_region_end());
#        my $exons = $output_transcript->get_all_Exons();
#        foreach my $exon (@$exons) {
#          $exon->start($exon->seq_region_start());
#          $exon->end($exon->seq_region_end());
#        }
        push(@$output_transcripts,$output_transcript);
      }
  }

#  unless(scalar(@$output_transcripts)) {
#    return;
#  }

  #$gene->slice(${$output_transcripts}[0]->slice());
#  say "Adding transcripts to new gene";
#  foreach my $output_transcript (@{$output_transcripts}) {
#    $gene->add_Transcript($output_transcript);
#  }
#  say "Finished adding transcripts";

#  return([$gene]);
  return($output_transcripts);
}


sub update_exonerate_transcript_coords {
  my ($self,$transcript,$target_slice_adaptor) = @_;

  $transcript->flush_supporting_features();
  my $slice_id = $transcript->start_Exon->seqname;
#  primary_assembly:HG02257.pri.mat.f1_v2:JAGYVH010000117.1:461283:475167:1
  $slice_id =~ /[^\:]+\:[^\:]+\:([^\:]+)\:([^\:]+)\:([^\:]+)\:[^\:]+$/;
  my $region_name = $1;
  my $start_offset = $2 - 1;
  my $slice = $target_slice_adaptor->fetch_by_region('toplevel',$region_name);

#  my $slice = $target_slice_adaptor->fetch_by_name($slice_id);
#  my $slice = $initial_output_transcript->start_Exon->slice->seq_region_Slice();
  say "FERGAL SLICE NAME: ".$slice->name();

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
  my $coverage_cutoff = 95;
  my $perc_id_cutoff = 95;
  my $cds_length_diff_cutoff = 0.05;

  foreach my $transcript (@$target_transcripts) {
    my $source_transcript = $source_transcript_id_hash->{$transcript->stable_id()};
    my $source_transcript_seq;
    my $transcript_seq;
    my $cds_length_diff = 0;
    if($source_transcript->translation()) {
      $source_transcript_seq = $source_transcript->translateable_seq();
      $transcript_seq = $transcript->translateable_seq();

      unless($transcript_seq) {
        compute_translation($transcript);
        $transcript_seq = $transcript->translateable_seq();
      }

      unless($transcript_seq) {
        $self->warning("Couldn't find an ORF in transcript mapped from ".$source_transcript->stable_id().". Source transcript has an ORF");
        push(@$bad_transcripts,$transcript);
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

    if($coverage >= $coverage_cutoff and $percent_id >= $perc_id_cutoff and $cds_length_diff <= $cds_length_diff_cutoff) {
      push(@$good_transcripts,$transcript);
    } else {
      say "Transcript ".$source_transcript->stable_id()." failed check: ".$coverage." cov, ".$percent_id." perc_id, ".$cds_length_diff_cutoff." length diff";
      push(@$bad_transcripts,$transcript);
    }
  }

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

    my $cov_string = "cov: ".$transcript->{'cov'}." perc_id: ".$transcript->{'perc_id'};

    $transcript->source($source);
    $transcript->description($cov_string);
    $gene->add_Transcript($transcript);
  }

  $gene->{'parent_gene_id'} = $parent_gene_id;
  $gene->stable_id($parent_gene_stable_id);
  $gene->version($parent_gene_version);
  $gene->biotype($parent_gene_biotype);

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


sub old_code {
  my ($self) = @_;


#  my ($clusters, $unclustered) = cluster_Genes_by_coding_exon_overlap($transcriptomic_genes,$types_hash);
#  my $transcriptomic_clusters = [@$clusters, @$unclustered];



  # At this point there should mostly be a set of good transcripts from minimap2 and some from exonerate
  # Then there are two lists of bad transcripts. Anything in the exonerate list means there was no good transcript from either method,
  # so they need to be assessed in terms of possible methods of recovery
#  my $remaining_bad_transcripts;
#  if(scalar(@$good_transcripts) and scalar(@$bad_exonerate_transcripts)) {
#    say "Found a mix of good and bad transcripts, will attempt to recover bad transcripts in the same region";
#    $remaining_bad_transcripts = $self->access_transcripts_for_recovery($bad_exonerate_transcripts,$source_transcript_id_hash,$good_transcripts);
#  }


#  my $output_genes = [];
#  my $final_gene_hash = {};
#  foreach my $gene (@$output_genes) {
#    my $transcripts = $gene->get_all_Transcripts();
#    foreach my $transcript (@$transcripts) {
#      say "SECOND TRANSCRIPT: ".$transcript->start()." ".$transcript->end()." ".$transcript->strand()." ".$transcript->slice->name();
#      say $transcript->seq->seq();

#      unless($parent_gene_ids->{$transcript->stable_id()}) {
#        $self->throw("The following mapped transcript stable id was not found in the initial list of dbIDs: ".$transcript->stable_id());
#      }

#      my $parent_gene_id = $parent_gene_ids->{$transcript->stable_id()}->{'gene_id'};
#      my $biotype_group = $parent_gene_ids->{$transcript->stable_id()}->{'biotype_group'};
#      my $is_canonical = $parent_gene_ids->{$transcript->stable_id()}->{'is_canonical'};
#      my $source = $parent_gene_ids->{$transcript->stable_id()}->{'source'};

      # Don't see how this would work, but will test
#      if($is_canonical) {
#        $transcript->is_canonical(1);
#      }

#      unless($final_gene_hash->{$parent_gene_id}) {
#        $final_gene_hash->{$parent_gene_id} = [];
#      }

#      my $transcript_id = $transcript->stable_id();
#      my $source_transcript = $source_transcript_id_hash->{$transcript_id};
#      use Data::Dumper;
#      print "SOURCE DUMPER: ".Dumper($source_transcript);

#      if($biotype_group eq 'coding') {
#        compute_translation($transcript);
#        set_alignment_supporting_features($transcript,$source_transcript->translation()->seq(),$transcript->translation()->seq());
#      }

#      my ($coverage,$percent_id) = (0,0);
#      my $aligned_source_seq;
#      my $aligned_target_seq;
#      my $source_transcript_seq = $source_transcript->seq->seq();
#      my $transcript_seq = $transcript->seq->seq();

#      ($coverage,$percent_id,$aligned_source_seq,$aligned_target_seq) = align_nucleotide_seqs($source_transcript_seq,$transcript_seq);
#      $transcript->biotype($source_transcript->biotype());
#      $transcript->created_date($coverage);
#      $transcript->modified_date($percent_id);
#      my $cov_string = "cov: ".$coverage." perc_id: ".$percent_id;
#      if($source_transcript->translation()) {
#        say "Calculating CDS for transcript with stable id: ".$transcript->stable_id();
        # Test
#        compute_translation($transcript);
#        my $source_cds_seq = $source_transcript->translateable_seq();
#        my $computed_target_cds_seq = $transcript->translateable_seq();
#        unless($computed_target_cds_seq) {
#          $self->warning("Existing CDS was not found for transcript ".$transcript->stable_id.", computing translation");
#          compute_translation($transcript);
#          $computed_target_cds_seq = $transcript->translateable_seq();
#        }

#        my ($computed_cds_cov,$computed_cds_perc_id,$computed_aligned_cds_source_seq,$computed_aligned_cds_target_seq) = align_nucleotide_seqs($source_cds_seq,$computed_target_cds_seq);
#        $cov_string .= " cds_cov: ".$computed_cds_cov." cds_perc_id: ".$computed_cds_perc_id;
        # End test

=for comment
        map_cds_location($source_transcript,$transcript,$aligned_source_seq,$aligned_target_seq);
        my $good_cov = 99;
        my $good_ident = 95;
        if($transcript->translation()) {
          my $source_cds_seq = $source_transcript->translateable_seq();
          my $initial_target_cds_seq = $transcript->translateable_seq();
          my $initial_translation = $transcript->translation();
          my ($initial_cds_cov,$initial_cds_perc_id,$initial_aligned_cds_source_seq,$initial_aligned_cds_target_seq) = align_nucleotide_seqs($source_cds_seq,$initial_target_cds_seq);
          if($initial_cds_cov < $good_cov or $initial_cds_perc_id < $good_ident) {
            say "Initial translation failed cut-off thresholds, will try computing translation to compare";
            compute_translation($transcript);
            my $computed_target_cds_seq = $transcript->translateable_seq();
            my ($computed_cds_cov,$computed_cds_perc_id,$computed_aligned_cds_source_seq,$computed_aligned_cds_target_seq) = align_nucleotide_seqs($source_cds_seq,$computed_target_cds_seq);
            # At the moment just look at which combined value is better
            if(($computed_cds_cov + $computed_cds_perc_id) < ($initial_cds_cov + $initial_cds_perc_id)) {
              say "Going with the intial translation over the computed one";
              $transcript->translation($initial_translation);
              $cov_string .= " cds_cov: ".$initial_cds_cov." cds_perc_id: ".$initial_cds_perc_id;
            } else {
              say "Choosing the computed translation over the initial translation";
              $cov_string .= " cds_cov: ".$computed_cds_cov." cds_perc_id: ".$computed_cds_perc_id;
            }
          } else {
            $cov_string .= " cds_cov: ".$initial_cds_cov." cds_perc_id: ".$initial_cds_perc_id;
          }
        } else {
          say "No translation found for transcript: ".$transcript->stable_id();
          say "Will compute translation";
          compute_translation($transcript);
          my $source_cds_seq = $source_transcript->translateable_seq();
          my $computed_target_cds_seq = $transcript->translateable_seq();
          my ($computed_cds_cov,$computed_cds_perc_id,$computed_aligned_cds_source_seq,$computed_aligned_cds_target_seq) = align_nucleotide_seqs($source_cds_seq,$computed_target_cds_seq);
          $cov_string .= " cds_cov: ".$computed_cds_cov." cds_perc_id: ".$computed_cds_perc_id;
        }
=cut

#      } # End if($source_transcript->translation()

#      $transcript->source($source);
#      $transcript->description($cov_string);

#      push(@{$final_gene_hash->{$parent_gene_id}},$transcript);
#    }
#  } # End foreach my $gene

#  my $final_genes = [];
#  foreach my $gene_id (keys(%{$final_gene_hash})) {
#    my $transcripts = $final_gene_hash->{$gene_id};
#    my $gene = Bio::EnsEMBL::Gene->new(-analysis => $self->analysis);

    # Should change this at some point to have a separate has for gene meta data to be tidy
#    my $parent_gene_id = $parent_gene_ids->{${$transcripts}[0]->stable_id()}->{'gene_id'};
#    my $parent_gene_stable_id = $parent_gene_ids->{${$transcripts}[0]->stable_id()}->{'gene_stable_id'};
#    my $parent_gene_version = $parent_gene_ids->{${$transcripts}[0]->stable_id()}->{'gene_version'};
#    my $parent_gene_biotype = $parent_gene_ids->{${$transcripts}[0]->stable_id()}->{'gene_biotype'};

#    foreach my $transcript (@$transcripts) {
#      my $parent_transcript_stable_id = $parent_gene_ids->{$transcript->stable_id()}->{'transcript_stable_id'};
#      my $parent_transcript_version = $parent_gene_ids->{$transcript->stable_id()}->{'transcript_version'};
#      $transcript->stable_id($parent_transcript_stable_id);
#      $transcript->version($parent_transcript_version);
#      $gene->add_Transcript($transcript);
#    }

#    $gene->{'parent_gene_id'} = $parent_gene_id;
#    $gene->stable_id($parent_gene_stable_id);
 #   $gene->version($parent_gene_version);
#    $gene->biotype($parent_gene_biotype);
#    push(@$final_genes, $gene);
#  }

}


1;
