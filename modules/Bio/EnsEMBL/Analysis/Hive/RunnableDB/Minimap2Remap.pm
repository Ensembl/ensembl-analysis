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
use POSIX;
use List::Util qw(min max);

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name align_nucleotide_seqs);
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);
use Bio::EnsEMBL::Analysis::Runnable::Minimap2Remap;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my($self) = @_;

  my $adaptors = Bio::EnsEMBL::DBSQL::DBAdaptor::get_available_adaptors;
  $adaptors->{Sequence} = 'Bio::EnsEMBL::Analysis::Tools::FastaSequenceAdaptor';
  *Bio::EnsEMBL::DBSQL::DBAdaptor::get_available_adaptors = sub {return $adaptors};
  $self->create_analysis;
  my $genome_index = $self->param_required('genome_index');

#  $self->param('region_padding',10000);

  # Define the dna dbs
  my $source_dna_dba = $self->hrdb_get_dba($self->param('source_dna_db'));
  my $target_dna_dba = $self->hrdb_get_dba($self->param('target_dna_db'));
  $self->hrdb_set_con($source_dna_dba,'source_dna_db');
  $self->hrdb_set_con($target_dna_dba,'target_dna_db');
  $source_dna_dba->get_SequenceAdaptor->fasta($self->param_required('source_dna_fasta'));
  my $target_dna_fasta = $self->param_required('genome_index');
  $target_dna_fasta =~ s/(\.fa(sta)?).*$/$1/;
  $target_dna_dba->get_SequenceAdaptor->fasta($target_dna_fasta);

  # Define the source and target gene dbs
  my $source_gene_dba = $self->hrdb_get_dba($self->param('source_gene_db'));
  my $target_gene_dba = $self->hrdb_get_dba($self->param('target_gene_db'));
  $source_gene_dba->dnadb($source_dna_dba);
  $target_gene_dba->dnadb($target_dna_dba);
  $self->hrdb_set_con($source_gene_dba,'source_gene_db');
  $self->hrdb_set_con($target_gene_dba,'target_gene_db');

# [19649] -> this is a gene/batch mapped to the -1 strand in the target, this e.g. has a 0 dist and no issues, so should be interesting
#  my $input_genes = $self->fetch_input_genes_by_id([59266,59269,59270,59293,59271,59274,59295,59276,59278],$source_gene_dba);
#my $input_genes = $self->fetch_input_genes_by_id([59266],$source_gene_dba);

  my $input_id_file = $self->param_required('input_id_file');
  unless(-e $input_id_file) {
    $self->throw("Did not find an input id file containing the list of source genes to be projected/mapped");
  }

  my $initial_gene_ids_hash = {};
  open(IN,$input_id_file);
  while(<IN>) {
    my $line = $_;
    chomp $line;
    if($line =~ /^\d+\tENS/) {
      my @eles = split("\t",$line);
      $initial_gene_ids_hash->{$eles[0]} = $eles[1];
    }
  }
  close IN;

  # TEST
#  my $test_slice = $target_gene_dba->get_SliceAdaptor->fetch_by_region('toplevel','17');
#  my $initial_projected_genes = $target_gene_dba->get_GeneAdaptor->fetch_all_by_Slice($test_slice);
#  my $initial_projected_genes = $target_gene_dba->get_GeneAdaptor->fetch_all();

#  my ($missing_genes,$problematic_transcripts) = $self->analyses_initial_projections($initial_projected_genes,$input_id_file);

  my $input_genes = $self->fetch_source_genes($input_id_file,$source_gene_dba);
#  my $test_slice = $source_gene_dba->get_SliceAdaptor->fetch_by_region('toplevel','7');
#  my $input_genes = $source_gene_dba->get_GeneAdaptor->fetch_all_by_Slice($test_slice);

  my $sequence_adaptor = $source_dna_dba->get_SequenceAdaptor();
  my $target_gene_adaptor = $target_gene_dba->get_GeneAdaptor();


  say "Processing ".scalar(@$input_genes)." genes";


################
# Put in code to sort input genes (might just be pre sorted by api, but just in case)
# Then for each gene record two genes to the left and right. Store these in a gene based hash
# Each key in the hash should be a gene id that points at up to 4 other gene ids that are keys on a hash
# Later on, once the final mapped set has been created and sorted, it will look at the current gene id
# and get the gene ids to the left and right and then compare them to the entries in the hash created now
# If we have two or more matches then the we can be confident in the gene's location
################

  my $sorted_input_genes = [sort { $a->slice->name() cmp $b->slice->name() or
                                   $a->start() <=> $b->start() or
                                   $a->end() <=> $b->end() }  @{$input_genes}];


  my $parent_gene_id_hash = {};

  $self->set_parent_info($sorted_input_genes,$parent_gene_id_hash,$sequence_adaptor);


  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Minimap2Remap->new(
       -analysis          => $self->analysis,
       -program           => $self->param('minimap2_path'),
       -paftools_path     => $self->param('paftools_path'),
       -genome_index      => $genome_index,
       -source_adaptor    => $source_gene_dba,
       -target_adaptor    => $target_gene_dba,
       -parent_genes      => $sorted_input_genes,
       -parent_gene_ids   => $parent_gene_id_hash,
       -no_projection     => $self->param('no_projection'),
  );
  $self->runnable($runnable);
}


sub run {
  my ($self) = @_;
  $self->dbc->disconnect_if_idle() if ($self->param('disconnect_jobs'));
  foreach my $runnable(@{$self->runnable}){
    $runnable->run;
    if ($self->can('filter_results')) {
      $self->output($self->filter_results($runnable->output));
    }
    else {
      $self->output($runnable->output);
    }
  }

  return $self->output;
}


sub write_output {
  my ($self) = @_;
  my $output_dba = $self->hrdb_get_con('target_gene_db');

  my $target_dna_dba = $self->hrdb_get_dba($self->param('target_dna_db'));
  my $target_gene_dba = $self->hrdb_get_dba($self->param('target_gene_db'));
  $target_gene_dba->dnadb($target_dna_dba);

  my $output_gene_adaptor = $target_gene_dba->get_GeneAdaptor;
  my $output_genes = $self->output();
  foreach my $output_gene (@$output_genes) {
#    say "Final gene: ".$output_gene->stable_id()." ".$output_gene->seq_region_start.":".$output_gene->seq_region_end.":".$output_gene->seq_region_strand.":".$output_gene->seq_region_name;
    if($output_gene->{'to_remove'} and !($output_gene->{'to_write'})) {
      $output_gene_adaptor->remove($output_gene);
    } elsif($output_gene->{'to_write'} and !($output_gene->{'to_remove'})) {
      empty_Gene($output_gene);
      $output_gene_adaptor->store($output_gene);
    }
  }

  return 1;
}


sub fetch_source_genes {
  my ($self,$input_id_file,$source_gene_dba) = @_;

  my $source_genes = [];

  open(IN,$input_id_file) or $self->throw("Could not open $input_id_file");
  while(<IN>) {
    my $line = $_;
    my @eles = split("\t",$line);
    my $gene = $source_gene_dba->get_GeneAdaptor->fetch_by_dbID($eles[0]);
    unless($gene) {
      $self->throw("Could not fetch the following gene with dbID ".$eles[0]." from the source gene db");
    }
    push(@$source_genes,$gene);
  }
  close IN or $self->throw("Could not close $input_id_file");

  return($source_genes);
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


sub create_filename {
  my ($self, $stem, $ext, $dir, $no_clean) = @_;
  return create_file_name($stem, $ext, $dir, $no_clean);
}


sub set_parent_info {
  my ($self,$sorted_input_genes,$parent_gene_id_hash,$sequence_adaptor) = @_;

  for(my $i=0; $i<scalar(@$sorted_input_genes); $i++) {
    my $gene = ${$sorted_input_genes}[$i];
    my $gene_id = $gene->dbID();
    my $gene_stable_id = $gene->stable_id();
    my $gene_version = $gene->version();
    my $gene_biotype = $gene->biotype();

    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts) {
      my $transcript_id = $transcript->dbID();
      my $transcript_stable_id = $transcript->stable_id();
      my $transcript_version = $transcript->version();
      my $biotype = $transcript->get_Biotype();
      my $biotype_group = $biotype->biotype_group();
      my $is_canonical = $transcript->is_canonical();
      my $source = $transcript->source();

      $parent_gene_id_hash->{$transcript_id}->{'gene_id'} = $gene_id;
      $parent_gene_id_hash->{$transcript_id}->{'gene_stable_id'} = $gene_stable_id;
      $parent_gene_id_hash->{$transcript_id}->{'gene_version'} = $gene_version;
      $parent_gene_id_hash->{$transcript_id}->{'gene_biotype'} = $gene_biotype;
      $parent_gene_id_hash->{$transcript_id}->{'transcript_stable_id'} = $transcript_stable_id;
      $parent_gene_id_hash->{$transcript_id}->{'transcript_version'} = $transcript_version;
      $parent_gene_id_hash->{$transcript_id}->{'biotype_group'} = $biotype_group;
      $parent_gene_id_hash->{$transcript_id}->{'is_canonical'} = $is_canonical;
      $parent_gene_id_hash->{$transcript_id}->{'source'} = $source;
    }
  }
}

sub batch_input_genes {
  my ($self,$genes,$sequence_adaptor) = @_;

  my $max_batch_span = 100000;
  my $anchor_size = 10000;
  my $anchor_dist = 5000;

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
      $batched_input_genes->{$batch_id}->{'length'} = $current_batch_end - $current_batch_start + 1;
      $batched_input_genes->{$batch_id}->{'midpoint'} = $current_batch_start + ceil($batched_input_genes->{$batch_id}->{'length'}/2);
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
    $batched_input_genes->{$batch_id}->{'length'} = $current_batch_end - $current_batch_start + 1;
    $batched_input_genes->{$batch_id}->{'midpoint'} = $current_batch_start + ceil($batched_input_genes->{$batch_id}->{'length'}/2);
  }

  return($batched_input_genes);
}


sub print_batch_details {
  my ($self,$batched_input_genes) = @_;

  foreach my $id (keys(%$batched_input_genes)) {
    my $batch = $batched_input_genes->{$id};
    say "Batch info: ID: ".$id.", Slice: ".$batch->{'slice_name'}.", Start: ".$batch->{'start'}.", End: ".$batch->{'end'}.", Midpoint: ".$batch->{'midpoint'}.", Length: ".$batch->{'length'};
    my $genes = $batch->{'genes'};
    foreach my $gene (@$genes) {
      say "  Gene: ID: ".$gene->stable_id().", Slice: ".$gene->seq_region_name().", Start: ".$gene->seq_region_start().", End: ".$gene->seq_region_end().", Strand: ".$gene->strand();
    }
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
#    say "ANCHOR LEFT:\n".$anchor_left_seq;

    my $anchor_middle_start = $batch_midpoint - $batch_offset;
    my $anchor_middle_end = $batch_midpoint + $batch_offset;

    if($anchor_middle_start < 1) {
      $anchor_middle_start = 1;
    }

    if($anchor_middle_end > $batch_slice->length()) {
      $anchor_middle_end = $batch_slice->length();
    }

    my $anchor_middle_seq = ${ $sequence_adaptor->fetch_by_Slice_start_end_strand($batch_slice,$anchor_middle_start,$anchor_middle_end,1) };
#    say "ANCHOR MIDDLE:\n".$anchor_middle_seq;

    my $anchor_right_start = $batch_end - $batch_offset;
    my $anchor_right_end = $batch_end + $batch_offset;

    if($anchor_right_start < 1) {
      $anchor_right_start = 1;
    }

    if($anchor_right_end > $batch_slice->length()) {
      $anchor_right_end = $batch_slice->length();
    }

    my $anchor_right_seq = ${ $sequence_adaptor->fetch_by_Slice_start_end_strand($batch_slice,$anchor_right_start,$anchor_right_end,1) };
#    say "ANCHOR RIGHT:\n".$anchor_right_seq;

    $batch->{'anchor_seqs'} = [[$anchor_left_start,$anchor_left_end,$anchor_left_seq],
                               [$anchor_middle_start,$anchor_middle_end,$anchor_middle_seq],
                               [$anchor_right_start,$anchor_right_end,$anchor_right_seq]];
  } # foreach my $id

}


sub map_anchors {
  my ($self,$batched_input_genes,$genome_index) = @_;

  foreach my $id (keys(%$batched_input_genes)) {
    say "Mapping Anchors for batch ID: ".$id;
    my $batch = $batched_input_genes->{$id};
    my $anchors = $batch->{'anchor_seqs'};
    my $batch_file = "batch_".$id.".fa";
    my $batch_output_file = $batch_file.".paf";
    open(OUT,">".$batch_file);
    say OUT ">al\n".${$anchors}[0][2];
    say OUT ">am\n".${$anchors}[1][2];
    say OUT ">ar\n".${$anchors}[2][2];
    close OUT;

    my $minimap2_command = $self->param('minimap2_path')." --cs --secondary=yes -x map-ont -N 10 ".$genome_index." ".$batch_file." > ".$batch_output_file;
#    if($id eq '448') {
#    if($id eq '889') {
    if(-e $batch_output_file) {
#      system($minimap2_command);
      my $target_regions = $self->calculate_target_regions($batch_output_file,$anchors);
    }
  }
}

sub calculate_target_regions {
  my ($self,$batch_output_file,$anchors) = @_;

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

  open(IN,$batch_output_file);
  while(<IN>) {
    my $line = $_;
    my @eles = split("\t",$line);
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
    unless($alignment_identities/$source_hit_length >= 0.99) {
      say "Hit fails the identity cutoff: (".$target_genomic_start."/".$target_genomic_end.")";
      $failed_identity = 1;
#      next;
    }

    unless($source_hit_length >= ($source_length * 0.95)) {
      say "Hit fails the coverage cutoff: (".$target_genomic_start."/".$target_genomic_end."), Hit start: ".$source_hit_start.", Hit end: ".$source_hit_end;
      $failed_coverage = 1;
#      next;
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
  close IN;

#  my $anchor_sets = [];

  my $anchor_sets = $self->get_best_anchor_sets($source_anchor1_midpoint,$source_anchor2_midpoint,$source_anchor3_midpoint,$anchors_left_res,$anchors_middle_res,$anchors_right_res);
  foreach my $anchor_set (@$anchor_sets) {
    my @a1_eles = @{${$anchor_set}[0][1]};
    my @a2_eles = @{${$anchor_set}[1][1]};
    my @a3_eles = @{${$anchor_set}[2][1]};
    my $dist = ${$anchor_set}[3];
    my $issue_count = ${$anchor_set}[4];
    say "Chosen anchor: (".$a1_eles[7].":".$a1_eles[8].") (".$a2_eles[7].":".$a2_eles[8].") (".$a3_eles[7].":".$a3_eles[8]."). Dist: ".$dist.", Issue count: ".$issue_count;
  }

  if(scalar(@$anchor_sets)) {
    return($anchor_sets);
  }

  say "No suitable anchor set found from initial results!";
}


sub get_best_anchor_sets {
  my ($self,$source_anchor1_midpoint,$source_anchor2_midpoint,$source_anchor3_midpoint,$anchors_left_res,$anchors_middle_res,$anchors_right_res) = @_;
  my $best_anchors = [];

  my $source_dist1 = abs($source_anchor2_midpoint - $source_anchor1_midpoint);
  my $source_dist2 = abs($source_anchor3_midpoint - $source_anchor2_midpoint);
  my $source_dist3 = abs($source_anchor3_midpoint - $source_anchor1_midpoint);

#  my $best_anchor_set;
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

#    if($a1_failed_coverage or $a1_failed_identity) {
#      next;
#    }
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
#      if($a2_strand ne $a1_strand or $a1_region ne $a2_region or $a2_failed_coverage or $a2_failed_identity) {
#        next;
#      }

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

#        if($a3_strand ne $a1_strand or $a2_region ne $a3_region or $a3_failed_coverage or $a3_failed_identity) {
#          next;
 #       }
        my $dist2 = abs($a2_mid - $a3_mid);
        my $dist3 = abs($a1_mid - $a3_mid);
        my $final_dist = abs($source_dist1 - $dist1) + abs($source_dist2 - $dist2) + abs($source_dist3 - $dist3);
        my $co_linear = 0;
        if($a1_strand eq '+' and $a2_mid >= $a1_mid and $a3_mid >= $a2_mid) {
          $co_linear = 1;
        } elsif($a1_strand eq '-' and $a2_mid <= $a1_mid and $a3_mid <= $a2_mid) {
          $co_linear = 1;
        } else {
          $issue_tracker++;
        }

 #       unless($co_linear) {
 #         say "Skipping anchor set as they are no co-linear";
 #         next;
 #       }
        say "  Current anchors: (".$a1_eles[7].":".$a1_eles[8].") (".$a2_eles[7].":".$a2_eles[8].") (".$a3_eles[7].":".$a3_eles[8]."). Dist: ".$final_dist.", Co-linear: ".$co_linear.", Issue count: ".$issue_tracker;

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
#        if($best_anchor_set) {
#          my $smallest_current_dist = ${$best_anchor_set}[3];
#          if($final_dist < $smallest_current_dist) {
#            $best_anchor_set = [$a1,$a2,$a3,$final_dist];
#          }
#        } else {
#         say "  Current anchors: (".$a1_eles[7].":".$a1_eles[8].") (".$a2_eles[7].":".$a2_eles[8].") (".$a3_eles[7].":".$a3_eles[8]."). Dist: ".$final_dist;
#         $best_anchor_set = [$a1,$a2,$a3,$final_dist];
 #       }
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
#      push(@$anchor_sets,$anchor_set);
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


1;
