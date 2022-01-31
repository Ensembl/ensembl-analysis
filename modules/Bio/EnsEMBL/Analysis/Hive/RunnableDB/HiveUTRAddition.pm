=head1 LICENSE

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveUTRAddition

=head1 SYNOPSIS


=head1 DESCRIPTION

Fetch a selection of genes which do not have UTR (acceptor) and try to add UTR
using a selection of cDNAs, RNA-seq, IsoSeqs (donor). The donor are ranked based
their correctness, usually cDNAs are considered to be the full sequence. If IsoSeqs
have been 5' and 3' capped they can also fall into the best category. RNA-seq data
the 3' capped data are seen as low quality donors at the moment.

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveUTRAddition;

use strict;
use warnings;
use feature 'say';

use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils qw(make_types_hash_with_genes cluster_Genes);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils qw(
                                                                transfer_supporting_evidence
                                                                validate_Exon_coords
                                                               );
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(
                                                                      calculate_exon_phases
                                                                      is_Transcript_sane
                                                                      all_exons_are_valid
                                                                      intron_lengths_all_less_than_maximum
                                                                      set_start_codon
                                                                      set_stop_codon
                                                                      clone_Transcript
                                                                      has_no_unwanted_evidence
                                                                     );
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(
                                                                       validate_Translation_coords
                                                                       compute_translation
                                                                       contains_internal_stops
                                                                       print_Translation
                                                                       create_Translation
                                                                      );
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene validate_store);
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: Defaults parameters are
               allow_partial_match => 0,
               allowed_input_sets => undef,
               min_size_5prime => 20,
               _max_5prime_utr_size => 200,
               copy_only => 0, # This is for jobs that fail even with lots of mem. Just copy the genes without trying to add UTR
               validate_store => 0, # This will check that the final stored gene is identical to what's in memory
               collapse_redundant_genes => 1,
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    allow_partial_match => 0,
    allowed_input_sets => undef,
    min_size_5prime => 20,
    _max_5prime_utr_size => 200,
    copy_only => 0, # This is for jobs that fail even with lots of mem. Just copy the genes without trying to add UTR
    validate_store => 0, # This will check that the final stored gene is identical to what's in memory
    collapse_redundant_genes => 1,
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Fetch all the genes which can provide UTR (donor) and
              all the genes which do not have UTR or for which the
              UTR could be improved like RNA-seq data sets
 Returntype : None
 Exceptions : Throws if 'dna_db' is not provided
              Throws if 'target_db' is not provided
              Throws if 'source_dbs' is not provided
              Throws if 'iid_type' is not provided

=cut

sub fetch_input {
  my ($self) = @_;

  $self->create_analysis;

  my $target_dba = $self->hrdb_get_dba($self->param_required('target_db'));
  my $dna_dba;
  if($self->param('use_genome_flatfile')) {
    unless($self->param_required('genome_file') && -e $self->param('genome_file')) {
      $self->throw("You selected to use a flatfile to fetch the genome seq, but did not find the flatfile. Path provided:\n".$self->param('genome_file'));
    }
    setup_fasta(
                 -FASTA => $self->param_required('genome_file'),
               );
  } else {
    $dna_dba = $self->hrdb_get_dba($self->param_required('dna_db'));
    $target_dba->dnadb($dna_dba);
  }

  $self->hrdb_set_con($target_dba,'target_db');

  my $input_id_type = $self->param_required('iid_type');

  # Checking that the Hash is set, it will be accessed with $self->biotype_priorities
  if (!exists $self->param_required('utr_biotype_priorities')->{other}) {
    my ($highest_key) = sort { $b <=> $a} values %{$self->param('utr_biotype_priorities')};
    $self->param_required('utr_biotype_priorities')->{other} = $highest_key+1;
  }

  my $input_id = $self->param('iid');

  if($self->param('iid_type') eq 'slice') {
    $self->query($self->fetch_sequence($input_id, $target_dba));
  } else {
    $self->throw("You must specify an input_id type in the config using the 'iid_type' parameter");
  }

  $self->say_with_header("Fetching input genes...");
  my $acceptor_genes = $self->filter_input_genes($self->param_required('acceptor_dbs'), $self->param('allowed_input_sets'));
  $self->complete_early('No genes found in the acceptor databases') unless (@$acceptor_genes);

  my $donor_genes = $self->filter_input_genes($self->param_required('donor_dbs'), $self->param('allowed_input_sets'),1);

  $self->say_with_header("Found ".scalar(@$acceptor_genes)." acceptor genes");
  $self->say_with_header("Found ".scalar(@$donor_genes)." donor genes");

  if($self->param('collapse_redundant_genes')) {
     $self->say_with_header("Removing redundant acceptors (will transfer supporting evidence)");
     $acceptor_genes = $self->remove_redundant_acceptors($acceptor_genes);
     $self->say_with_header("Post filtering the acceptor count is: ".scalar(@$acceptor_genes));

     $self->say_with_header("Removing redundant donors (will consider priority and transfer supporting evidence)");
     $donor_genes = $self->remove_redundant_donors($donor_genes);
     $self->say_with_header("Post filtering the donor count is: ".scalar(@$donor_genes));
  }

  $self->acceptor_genes($acceptor_genes);
  $self->donor_genes($donor_genes);

}


=head2 run

 Arg [1]    : None
 Description: Cluster the acceptor genes with the donor genes. The cluster with
              acceptor and donor genes will be processed
 Returntype : None
 Exceptions : None

=cut

sub run {
  my $self = shift;

  # This is a last resort when jobs start requiring very large mem allocations. Just copy the genes without doing anything in terms of UTR
  # May remove in future if we bring the footprint back down
  if($self->param('copy_only')) {
    $self->warning("Copy only mode selected, no attempts will be made to add UTR. Acceptor genes will be copied unmodified into the output db");
    $self->output($self->acceptor_genes);
    return;
  }

  unless(scalar(@{$self->acceptor_genes})) {
    return;
  }

  my ($type_hash, $genes) = @{make_types_hash_with_genes($self->donor_genes, $self->acceptor_genes, 'donor', 'acceptor')};

  my ($clusters, $unclustered) = cluster_Genes($genes, $type_hash);
  my @genes;
  foreach my $cluster (@$clusters) {
    if (@genes%100 == 0) {
      print scalar(@genes), "\n";
    }
    my $donor_transcripts;
    foreach my $gene (@{$cluster->get_Genes_by_Set('donor')}) {
      push(@$donor_transcripts, @{$gene->get_all_Transcripts});
    }
    foreach my $gene (@{$cluster->get_Genes_by_Set('acceptor')}) {
      my $acceptor_transcripts = $gene->get_all_Transcripts;
      $gene->flush_Transcripts;
      foreach my $acceptor_transcript (@$acceptor_transcripts) {
        if($acceptor_transcript->translate and index($acceptor_transcript->translate->seq, '*') == -1) {
          my $transcript = $self->add_utr($acceptor_transcript, $donor_transcripts);
          if ($transcript) {
            $gene->add_Transcript($transcript);
          }
          else {
            $self->throw("No transcript returned, even if no UTR found the original transcript should always be returned");
          }
        }
        else {
          $self->warning("Transcript with dbID ".$acceptor_transcript->dbID." has a stop codon in the translation. Skipping UTR addition");
        }
      }
      push(@genes, $gene) if (@{$gene->get_all_Transcripts});
    }
  }

  foreach my $single_cluster (@$unclustered) {
    foreach my $single_gene (@{$single_cluster->get_Genes_by_Set('acceptor')}) {
      foreach my $transcript (@{$single_gene->get_all_Transcripts}) {
        calculate_exon_phases($transcript, 0);
      }
      push(@genes, $single_gene);
    }
  }
  $self->output(\@genes);
}


=head2 write_output

 Arg [1]    : None
 Description: Write the genes with or without UTR into the 'target_db'
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my $self = shift;

  my $adaptor = $self->hrdb_get_con('target_db')->get_GeneAdaptor;

  foreach my $gene (@{$self->output}){
    empty_Gene($gene);
    $adaptor->store($gene);

    # We had this because there was an issue with exons not being written. The problem would pop up when genebuilder failed
    # Seems to have been linked to parts of transcripts being off the edge of slices. This has been fixed so this code is
    # disabled by default to save server usage
    if($self->param('validate_store')) {
      my $stored_gene = $adaptor->fetch_by_dbID($gene->dbID);
      my $validate_store = validate_store($gene,$stored_gene);
      unless($validate_store) {
        $self->throw("The stored gene and the gene in memory differed in some key features. Something went wrong on the store");
      }
    }
  }
}


=head2 add_utr

 Arg [1]    : Bio::EnsEMBL::Transcript The transcript to add UTR to
 Arg [2]    : Arrayref of Bio::EnsEMBL::Transcript The transcripts which have UTRs
 Description: Add UTR to Arg[1] if one of the transcript from Arg[2] has the same
              intron structure. The objects from Arg[2] have priorities which means
              that if we add UTR with a transcript from the group 1 we do not try
              to add UTR from group 2.
 Returntype : Bio::EnsEMBL::Transcript
 Exceptions : None

=cut

sub add_utr {
  my ($self, $acceptor_transcript, $donor_transcripts) = @_;

  my $allow_partial_match = $self->param('allow_partial_match');

  my $modified_acceptor_transcript_5prime;
  my $modified_acceptor_transcript_3prime;
  my $transcript;
  my $modified_acceptor_transcript_single_exon;
  if(scalar(@{$acceptor_transcript->get_all_translateable_Exons}) == 1) {
    # I don't like this duplication of code, should find a better way to do it
    $acceptor_transcript->{'5_prime_utr'} = 0;
    $acceptor_transcript->{'3_prime_utr'} = 0;
    foreach my $donor_transcript (@{$donor_transcripts}) {
      my $priority = $self->biotype_priorities($donor_transcript->biotype);
      unless($priority) {
        $self->warning("Transcript biotype was not found in the biotype priorities hash or biotype was set to 0 priority. Skipping.".
                       "Biotype: ".$donor_transcript->biotype);
        next;
      }

      if($acceptor_transcript->{'has_utr'}) {
        # First check if the donor priority is worse (1=best), if it's worse then just skip
        if($priority > $modified_acceptor_transcript_single_exon->{'priority'}) {
          next;
        }
        my $new_transcript_single_exon = $self->add_single_exon_utr($acceptor_transcript,$donor_transcript);
        if($new_transcript_single_exon && ($new_transcript_single_exon->length() > $modified_acceptor_transcript_single_exon->length())) {
          $modified_acceptor_transcript_single_exon = $new_transcript_single_exon;
          $modified_acceptor_transcript_single_exon->{'priority'} = $priority;
          $modified_acceptor_transcript_single_exon->{'donor_5prime_biotype'} = $donor_transcript->biotype;
          $modified_acceptor_transcript_single_exon->{'donor_3prime_biotype'} = $donor_transcript->biotype;
        }
      } else {
        $modified_acceptor_transcript_single_exon = $self->add_single_exon_utr($acceptor_transcript,$donor_transcript);
        if($modified_acceptor_transcript_single_exon) {
          $modified_acceptor_transcript_single_exon->{'priority'} = $priority;
          $modified_acceptor_transcript_single_exon->{'donor_5prime_biotype'} = $donor_transcript->biotype;
          $modified_acceptor_transcript_single_exon->{'donor_3prime_biotype'} = $donor_transcript->biotype;
        }
      }
    }

    if($modified_acceptor_transcript_single_exon) {
      $modified_acceptor_transcript_single_exon->biotype($acceptor_transcript->biotype);
      $self->add_transcript_supporting_features($modified_acceptor_transcript_single_exon,$acceptor_transcript);
      $transcript = $modified_acceptor_transcript_single_exon;
      return($transcript);
    }
  } else {
    $acceptor_transcript->{'5_prime_utr'} = 0;
    $acceptor_transcript->{'3_prime_utr'} = 0;
    my $cds_transcript = $self->create_cds_transcript($acceptor_transcript);
    my $utr_transcript = $self->find_soft_match($cds_transcript,$donor_transcripts);

    if($utr_transcript) {
      unless($acceptor_transcript->translation->seq eq $utr_transcript->translation->seq) {
        $self->throw("The new transcript has a different translation to the original.\nOriginal:\n".
                     $acceptor_transcript->translation->seq."\nNew:\n".$utr_transcript->translation->seq);
      }
      $utr_transcript->biotype($acceptor_transcript->biotype);
      $self->add_transcript_supporting_features($utr_transcript,$acceptor_transcript);
      $self->look_for_both($utr_transcript);
      if($utr_transcript->three_prime_utr && $utr_transcript->{'donor_3prime_biotype'} =~ /^rnaseq/) {
        $utr_transcript = $self->trim_3prime_utr_short_read($utr_transcript);
      }
      if($utr_transcript->five_prime_utr && ($utr_transcript->{'donor_5prime_biotype'} =~ /^rnaseq/ or @{$utr_transcript->get_all_five_prime_UTRs} > 4)) {
        $utr_transcript = $self->trim_5prime_utr_short_read($utr_transcript);
      }
      calculate_exon_phases($utr_transcript, 0);
      return($utr_transcript);
    }
  }

  return($acceptor_transcript);
}


=head2 add_utr_evidence

 Arg [1]    : Arrayref of Bio::EnsEMBL::Exon from the acceptor transcript
 Arg [2]    : Arrayref of Bio::EnsEMBL::Exon from the donor transcript
 Arg [3]    : Bio::EnsEMBL::Transcript donor transcript
 Description: Add the supproting evidences from the donor exons to the
              acceptor exons.
 Returntype : None
 Exceptions : None

=cut

sub add_utr_evidence {
  my ($self, $final_exons, $utr_exons, $utr_transcript) = @_;

  for (my $i = 0; $i < @$final_exons; $i++) {
    my $exon = $final_exons->[$i];
    my @new_sfs;
    for (my $j = $i; $j < @$utr_exons; $j++) {
      my $utr_exon = $utr_exons->[$j];
      next if ($exon->start > $utr_exon->end or $exon->end < $utr_exon->start);
      my @sfs = grep {$_->isa('Bio::EnsEMBL::DnaDnaAlignFeature')} @{$exon->get_all_supporting_features};
      if (scalar(@sfs)) {
        foreach my $sf (@sfs) {
          my $dnaalignfeature = Bio::EnsEMBL::DnaDnaAlignFeature->new(-features => [$sf->ungapped_features], -align_type => 'ensembl');
          push(@new_sfs, $dnaalignfeature);
        }
      }
#      else {
#        my $dnaalignfeature = Bio::EnsEMBL::DnaDnaAlignFeature->new(-cigar_string => $exon->length.'M');
#        $dnaalignfeature->start($exon->start);
#        $dnaalignfeature->end($exon->end);
#        $dnaalignfeature->strand($exon->strand);
#        $dnaalignfeature->analysis($exon->analysis);
#        $dnaalignfeature->slice($exon->slice);
#        $dnaalignfeature->hseqname('rnaseq');
#        $dnaalignfeature->hstrand(1);
#        $dnaalignfeature->hstart($utr_exon->cdna_start($utr_transcript));
#        $dnaalignfeature->hend($utr_exon->cdna_end($utr_transcript));
#        push(@new_sfs, $dnaalignfeature);
#     }
    }
    $exon->add_supporting_features(@new_sfs);
  }
}


=head2 add_five_prime_utr

 Arg [1]    : Bio::EnsEMBL::Transcript acceptor transcript
 Arg [2]    : Bio::EnsEMBL::Transcript donor transcript
 Arg [3]    : Arrayref of Bio::EnsEMBL::Intron from the acceptor transcript
 Arg [4]    : Arrayref of Bio::EnsEMBL::Intron from the donor transcript
 Arg [5]    : String representing the intron structure of the acceptor transcript
 Arg [6]    : String representing the intron structure of the donor transcript
 Description: Try to add UTR to the 5' end of Arg[1] using Arg[2]
 Returntype : Bio::EnsEMBL::Transcript or 0 if UTR was not added
 Exceptions : Throws if it could not get the first non coding exon
              Throws if the translation changed after adding the UTR

=cut

sub add_five_prime_utr {
  my ($self,$transcript_a,$transcript_b,$introns_acceptor,$introns_b,$cds_intron_string_a,$intron_string_b) = @_;

  my $strand = $transcript_a->strand;
  # At this point we have a match, now we need to locate the exon to merge
  my @cds_intron_coords_a = split(":",$cds_intron_string_a);
  my @intron_coords_b = split(":",$intron_string_b);

  my $five_prime_intron_a;
  my $exon_merge_index_b = -1;

  $five_prime_intron_a = $cds_intron_coords_a[0];
  for(my $i=0; $i<scalar(@intron_coords_b); $i++) {
    if($five_prime_intron_a eq $intron_coords_b[$i]) {
      $exon_merge_index_b = $i;
      last;
    }
  }

  # This should not be possible
  if($exon_merge_index_b == -1) {
    $self->throw("The transcript cds structure matched the donor transcript, but something went wrong when trying to match the coords\n".
                 "Intron coords string a:\n".$cds_intron_string_a.
                 "\nIntron coords string b:\n".$intron_string_b."\n");
  }

  my $exons_a = $transcript_a->get_all_Exons();
  my $exons_b = $transcript_b->get_all_Exons();
  foreach my $exon_b (@{$exons_b}) {
    if($exon_b->start == 0) {
      $exon_b->start(1);
    }
  }

  my $merge_exon_candidate_a = ${$exons_a}[0];
  my $merge_exon_candidate_b = ${$exons_b}[$exon_merge_index_b];

  if($strand == 1) {
    if($merge_exon_candidate_b->start >= $merge_exon_candidate_a->start) {
      if ($self->biotype_priorities($transcript_a->biotype) > $self->biotype_priorities($transcript_b->biotype) and
        $transcript_a->coding_region_start > $merge_exon_candidate_b->start and $transcript_a->coding_region_start-$merge_exon_candidate_a->start > $self->min_size_5prime) {
      }
      else {
        return(0);
      }
    } elsif($merge_exon_candidate_b->end != $merge_exon_candidate_a->end) {
      $self->warning("The internal boundry coord of the donor did not match the acceptor, something is wrong.\n".
                     "Donor boundry: ".$merge_exon_candidate_b->end."\nAcceptor boundry: ".$merge_exon_candidate_a->end);
      return(0);
    }
  } else {
    if($merge_exon_candidate_b->end <= $merge_exon_candidate_a->end) {
      if ($self->biotype_priorities($transcript_a->biotype) > $self->biotype_priorities($transcript_b->biotype) and
        $transcript_a->coding_region_end < $merge_exon_candidate_b->end and $merge_exon_candidate_b->end-$transcript_a->coding_region_end > $self->min_size_5prime) {
      }
      else {
        return(0);
      }
    } elsif($merge_exon_candidate_b->start != $merge_exon_candidate_a->start) {
      $self->warning("The internal boundry coord of the donor did not match the acceptor, something is wrong.\n".
                     "Donor boundry: ".$merge_exon_candidate_b->start."\nAcceptor boundry: ".$merge_exon_candidate_a->start);
      return(0);
    }
  }

  my $final_exons = [];
  # Note that for the moment I'm going to convert all exons before the merge candidate to non-coding (this should be true anyway)
  for(my $i=0; $i<$exon_merge_index_b; $i++) {
    my $exon = ${$exons_b}[$i];
    my $out_exon = new Bio::EnsEMBL::Exon(
                                           -START  => $exon->start,
                                           -END       => $exon->end,
                                           -STRAND    => $exon->strand,
                                           -SLICE     => $exon->slice,
                                           -ANALYSIS  => $self->analysis,
                                           -PHASE     => -1,
                                           -END_PHASE => -1);

    push(@{$final_exons},$out_exon);
  }

  my $translation_start = $transcript_a->translation->start;
  my $translation_end = $transcript_a->translation->end;
  # If the donor exon is shorter (or the same length), then just make the acceptor boundry exon the start exon and add
  # all the other ones
  if(($strand == 1 && ($merge_exon_candidate_b->start >= $merge_exon_candidate_a->start)) ||
     ($strand == -1 && ($merge_exon_candidate_b->end <= $merge_exon_candidate_a->end))) {
    my $start_exon = new Bio::EnsEMBL::Exon(
                                           -START     => $merge_exon_candidate_a->start,
                                           -END       => $merge_exon_candidate_a->end,
                                           -STRAND    => $merge_exon_candidate_a->strand,
                                           -SLICE     => $merge_exon_candidate_a->slice,
                                           -ANALYSIS  => $self->analysis,
                                           -PHASE     => $merge_exon_candidate_a->phase,
                                           -END_PHASE => $merge_exon_candidate_a->end_phase);
    if ($transcript_a->{'5_prime_utr'} == 0 or $transcript_a->{'5_prime_utr'} >= $self->biotype_priorities($transcript_b->biotype)) {
      $start_exon->start($merge_exon_candidate_b->start);
      $start_exon->end($merge_exon_candidate_b->end);
      $start_exon->strand($merge_exon_candidate_b->strand);
      $start_exon->slice($merge_exon_candidate_b->slice);
      $start_exon->analysis($self->analysis);
      $start_exon->phase($merge_exon_candidate_b->phase);
      $start_exon->end_phase($merge_exon_candidate_b->end_phase);
    }

    my $supporting_features_a = $merge_exon_candidate_a->get_all_supporting_features();
    $start_exon->add_supporting_features(@{$supporting_features_a});

    push(@{$final_exons},$start_exon);

    # Add the remaining exons
    for(my $i=1; $i<scalar(@{$exons_a}); $i++) {
      my $exon = ${$exons_a}[$i];
      my $out_exon = new Bio::EnsEMBL::Exon(
                                             -START  => $exon->start,
                                             -END       => $exon->end,
                                             -STRAND    => $exon->strand,
                                             -SLICE     => $exon->slice,
                                             -ANALYSIS  => $self->analysis,
                                             -PHASE     => $exon->phase,
                                             -END_PHASE => $exon->end_phase);
      $supporting_features_a = $exon->get_all_supporting_features();
      $out_exon->add_supporting_features(@{$supporting_features_a});
      push(@{$final_exons},$out_exon);
    }

  } else {
    my $merge_exon = new Bio::EnsEMBL::Exon(
                                             -START     => $merge_exon_candidate_b->start,
                                             -END       => $merge_exon_candidate_b->end,
                                             -STRAND    => $merge_exon_candidate_b->strand,
                                             -SLICE     => $merge_exon_candidate_b->slice,
                                             -ANALYSIS  => $self->analysis,
                                             -PHASE     => -1,
                                             -END_PHASE => $merge_exon_candidate_a->end_phase);

    my $supporting_features_a = $merge_exon_candidate_a->get_all_supporting_features();
    $merge_exon->add_supporting_features(@{$supporting_features_a});

    push(@{$final_exons},$merge_exon);

    if ($merge_exon_candidate_a->phase > 0) {
      $translation_start += $merge_exon_candidate_a->phase == 1 ? 2 : 1;
    }

    # Add the remaining exons
    for(my $i=1; $i<scalar(@{$exons_a}); $i++) {
      my $exon = ${$exons_a}[$i];
      my $out_exon = new Bio::EnsEMBL::Exon(
                                             -START  => $exon->start,
                                             -END       => $exon->end,
                                             -STRAND    => $exon->strand,
                                             -SLICE     => $exon->slice,
                                             -ANALYSIS  => $self->analysis,
                                             -PHASE     => $exon->phase,
                                             -END_PHASE => $exon->end_phase);
      push(@{$final_exons},$out_exon);
    }

  }

  my $genomic_start;
  my $genomic_end;
  if($transcript_a->strand == 1) {
    $genomic_start = ($transcript_a->translation->start_Exon->seq_region_start + $translation_start - 1);
    $genomic_end = ($transcript_a->translation->end_Exon->seq_region_start + $translation_end - 1);
  } else {
    $genomic_start = ($transcript_a->translation->start_Exon->seq_region_end - $translation_start + 1);
    $genomic_end = ($transcript_a->translation->end_Exon->seq_region_end - $translation_end + 1);
  }

  my $final_translation = create_Translation($final_exons, $genomic_start, $genomic_end);

  unless ($final_translation) {
    $self->throw("Failed to create a final translation");
  }

  foreach my $seq_edit (@{$transcript_a->translation->get_all_SeqEdits}) {
    $final_translation->add_Attributes($seq_edit->get_Attribute);
  }

  $self->add_utr_evidence($final_exons, $exons_b, $transcript_b);
  my $modified_transcript = Bio::EnsEMBL::Transcript->new(-EXONS => $final_exons);

  # This is a basic sanity check on the UTR itself. First we want to check if the transcript is now abnormally longer (> 100KB longer)
  # If it is then calculate the average 5' intron length of the UTR. If this average is 5' intron length is > 35K then throw it out
  my $utr_intron_count = scalar(@{$final_exons}) - scalar(@{$exons_a});
  my $big_utr_length = 50000;
  my $max_no_intron_extension = 10000;
  my $max_average_5_prime_intron_length = 35000;

  my $modified_transcript_length = $modified_transcript->seq_region_end() - $modified_transcript->seq_region_start();
  my $transcript_a_length = $transcript_a->seq_region_end() - $transcript_a->seq_region_start();
  my $added_length = $modified_transcript_length - $transcript_a_length;  

  if($utr_intron_count == 0 && $added_length > $max_no_intron_extension) {
    return(0);
  }

  $modified_transcript->analysis($transcript_a->analysis);
  $modified_transcript->biotype($transcript_a->biotype);
  $modified_transcript->slice($transcript_a->slice());
  $modified_transcript->translation($final_translation);

  calculate_exon_phases($modified_transcript, 0);


  my $old_translation = $transcript_a->translation->seq;
  my $new_translation = $modified_transcript->translate->seq;

  unless($transcript_a->translation->seq eq $modified_transcript->translate->seq && $modified_transcript->translate->seq eq $modified_transcript->translation->seq) {
    $old_translation =~ s/^X//;
    if ($old_translation eq $new_translation && $modified_transcript->translate->seq eq $new_translation) {
      $transcript_a->translation->start($transcript_a->translation->start+($transcript_a->start_Exon->phase > 0 and $transcript_a->start_Exon->phase == 1 ? 2 : 1));
      $transcript_a->start_Exon->phase(0);
    }
    else {
      $self->throw("There is an issue with the translation after UTR was added. Check above for the sequences, all three should match")
    }
  }

  $transcript_a->{'5_prime_utr'} = $self->biotype_priorities($transcript_b->biotype);
  $modified_transcript->add_supporting_features(grep {$_->isa('Bio::EnsEMBL::DnaDnaAlignFeature')} @{$transcript_b->get_all_supporting_features});

  return($modified_transcript);
}


=head2 add_three_prime_utr

 Arg [1]    : Bio::EnsEMBL::Transcript acceptor transcript
 Arg [2]    : Bio::EnsEMBL::Transcript donor transcript
 Arg [3]    : Arrayref of Bio::EnsEMBL::Intron from the acceptor transcript
 Arg [4]    : Arrayref of Bio::EnsEMBL::Intron from the donor transcript
 Arg [5]    : String representing the intron structure of the acceptor transcript
 Arg [6]    : String representing the intron structure of the donor transcript
 Description: Try to add UTR to the 3' end of Arg[1] using Arg[2]
 Returntype : Bio::EnsEMBL::Transcript or 0 if UTR was not added
 Exceptions : Throws if it could not get the first non coding exon
              Throws if the translation changed after adding the UTR

=cut

sub add_three_prime_utr {
  my ($self,$transcript_a,$transcript_b,$introns_acceptor,$introns_b,$cds_intron_string_a,$intron_string_b) = @_;

  my $strand = $transcript_a->strand;

  # At this point we have a match, now we need to locate the exon to merge
  my @cds_intron_coords_a = split(":",$cds_intron_string_a);
  my @intron_coords_b = split(":",$intron_string_b);

  my $three_prime_intron_a = $cds_intron_coords_a[$#cds_intron_coords_a];
  my $exon_merge_index_b = -1;
  for(my $i=0; $i<scalar(@intron_coords_b); $i++) {
    if($three_prime_intron_a eq $intron_coords_b[$i]) {
        $exon_merge_index_b = $i+1;
        last;
    }
  }

  # This should not be possible
  if($exon_merge_index_b == -1) {
    $self->throw("The transcript cds structure matched the donor transcript, but something went wrong when trying to match the coords");
  }

  my $exons_a = $transcript_a->get_all_Exons();
  my $exons_b = $transcript_b->get_all_Exons();
  my $merge_exon_candidate_a = ${$exons_a}[$#$exons_a];
  my $merge_exon_candidate_b = ${$exons_b}[$exon_merge_index_b];


  if($strand == 1) {
    if($merge_exon_candidate_b->end <= $merge_exon_candidate_a->end) {
      if ($self->biotype_priorities($transcript_a->biotype) > $self->biotype_priorities($transcript_b->biotype) and
        $transcript_a->coding_region_end < $merge_exon_candidate_b->end) {
      }
      else {
        return(0);
      }
    } elsif($merge_exon_candidate_b->start != $merge_exon_candidate_a->start) {
      $self->warning("The internal boundry coord of the donor did not match the acceptor, something is wrong.\n".
                          "Donor boundry: ".$merge_exon_candidate_b->start."\nAcceptor boundry: ".$merge_exon_candidate_a->start);
      return(0);
    }
  } else {
    if($merge_exon_candidate_b->start >= $merge_exon_candidate_a->start) {
      if ($self->biotype_priorities($transcript_a->biotype) > $self->biotype_priorities($transcript_b->biotype) and
        $transcript_a->coding_region_start > $merge_exon_candidate_b->start) {
      }
      else {
        return(0);
      }
    } elsif($merge_exon_candidate_b->end != $merge_exon_candidate_a->end) {
      $self->warning("The internal boundry coord of the donor did not match the acceptor, something is wrong.\n".
                     "Donor boundry: ".$merge_exon_candidate_b->end."\nAcceptor boundry: ".$merge_exon_candidate_a->end);
      return(0);
    }
  }

  my $final_exons = [];

  # First add all the exons from the acceptor
  foreach my $exon (@{$exons_a}) {
    my $out_exon = new Bio::EnsEMBL::Exon(
                                           -START  => $exon->start,
                                           -END       => $exon->end,
                                           -STRAND    => $exon->strand,
                                           -SLICE     => $exon->slice,
                                           -ANALYSIS  => $self->analysis,
                                           -PHASE     => $exon->phase,
                                           -END_PHASE => $exon->end_phase);

    my $supporting_features_a = $exon->get_all_supporting_features();
    $out_exon->add_supporting_features(@{$supporting_features_a});
    push(@{$final_exons},$out_exon);
  }


  # Now look at the merge candidates. There are a few things to check. If the start coord of the donor exon is greater than the start coord of the
  # acceptor then just push the acceptor exon on the final exons array
  if(($strand == 1 && ($merge_exon_candidate_b->end <= $merge_exon_candidate_a->end)) ||
     ($strand == -1 && ($merge_exon_candidate_b->start >= $merge_exon_candidate_a->start))) {
    if ($transcript_a->{'3_prime_utr'} == 0 or $transcript_a->{'3_prime_utr'} >= $self->biotype_priorities($transcript_b->biotype)) {
      my $merge_exon = new Bio::EnsEMBL::Exon(
                                               -START     => $merge_exon_candidate_b->start,
                                               -END       => $merge_exon_candidate_b->end,
                                               -STRAND    => $merge_exon_candidate_b->strand,
                                               -SLICE     => $merge_exon_candidate_b->slice,
                                               -ANALYSIS  => $self->analysis,
                                               -PHASE     => $merge_exon_candidate_a->phase,
                                               -END_PHASE => -1);

      $merge_exon->add_supporting_features(@{$merge_exon_candidate_a->get_all_supporting_features});

      pop(@{$final_exons});
      push(@{$final_exons},$merge_exon);
    }
  }  else {
    my $merge_exon = new Bio::EnsEMBL::Exon(
                                             -START     => $merge_exon_candidate_b->start,
                                             -END       => $merge_exon_candidate_b->end,
                                             -STRAND    => $merge_exon_candidate_b->strand,
                                             -SLICE     => $merge_exon_candidate_b->slice,
                                             -ANALYSIS  => $self->analysis,
                                             -PHASE     => $merge_exon_candidate_a->phase,
                                             -END_PHASE => -1);

    my $supporting_features_a = $merge_exon_candidate_a->get_all_supporting_features();
    $merge_exon->add_supporting_features(@{$supporting_features_a});

    pop(@{$final_exons});
    push(@{$final_exons},$merge_exon);

  }

  for(my $i=$exon_merge_index_b + 1; $i<scalar(@{$exons_b}); $i++) {
    my $exon = ${$exons_b}[$i];
    my $out_exon = new Bio::EnsEMBL::Exon(
                                           -START     => $exon->start,
                                           -END       => $exon->end,
                                           -STRAND    => $exon->strand,
                                           -SLICE     => $exon->slice,
                                           -ANALYSIS  => $self->analysis,
                                           -PHASE     => -1,
                                           -END_PHASE => -1);
    push(@{$final_exons},$out_exon);
  }

  my $genomic_start;
  my $genomic_end;
  if($transcript_a->strand == 1) {
    $genomic_start = ($transcript_a->translation->start_Exon->seq_region_start + $transcript_a->translation->start - 1);
    $genomic_end = ($transcript_a->translation->end_Exon->seq_region_start + $transcript_a->translation->end - 1);
  } else {
    $genomic_start = ($transcript_a->translation->start_Exon->seq_region_end - $transcript_a->translation->start + 1);
    $genomic_end = ($transcript_a->translation->end_Exon->seq_region_end - $transcript_a->translation->end + 1);
  }

  my $final_translation = create_Translation($final_exons, $genomic_start, $genomic_end);
  foreach my $seq_edit (@{$transcript_a->translation->get_all_SeqEdits}) {
    $final_translation->add_Attributes($seq_edit->get_Attribute);
  }
  $self->add_utr_evidence($final_exons, $exons_b, $transcript_b);
  my $modified_transcript = Bio::EnsEMBL::Transcript->new(-EXONS => $final_exons);

  # This is a basic sanity check on the UTR itself. First we want to check if the transcript is now abnormally longer (> 100KB longer)
  # If it is then calculate the average 3' intron length of the UTR. If this average is 3' intron length is > 25K then throw it out
  my $utr_intron_count = scalar(@{$final_exons}) - scalar(@{$exons_a});
  my $big_utr_length = 50000;
  my $max_no_intron_extension = 10000;
  my $max_average_3_prime_intron_length = 35000;

  my $modified_transcript_length = $modified_transcript->seq_region_end() - $modified_transcript->seq_region_start();
  my $transcript_a_length = $transcript_a->seq_region_end() - $transcript_a->seq_region_start();
  my $added_length = $modified_transcript_length - $transcript_a_length;

  if($utr_intron_count == 0 && $added_length > $max_no_intron_extension) {
    return(0);
  }

  $modified_transcript->analysis($transcript_a->analysis);
  $modified_transcript->biotype($transcript_a->biotype);
  $modified_transcript->slice($transcript_a->slice());
  $modified_transcript->translation($final_translation);


  my $modified_translation = $modified_transcript->translation();
  my $old_translation = $transcript_a->translation->seq;
  my $new_translation = $modified_transcript->translate->seq;

  unless($old_translation eq $new_translation && $modified_transcript->translate->seq eq $new_translation) {
    $old_translation =~ s/^X//;
    $self->throw("There is an issue with the translation after UTR was added. Check above for the sequences, all three should match")
      unless ($old_translation eq $new_translation && $modified_transcript->translate->seq eq $new_translation);
  }

  $transcript_a->{'3_prime_utr'} = 1;

  return($modified_transcript);
}




=head2 add_single_exon_utr

 Arg [1]    : Bio::EnsEMBL::Transcript the acceptor transcript
 Arg [2]    : Bio::EnsEMBL::Transcript the donor transcript
 Description: Add UTR to a single exon transcript
 Returntype : Bio::EnsEMBL::Transcript or 0 if could not add UTR
 Exceptions : Throws if the translation changed after adding the UTR

=cut

sub add_single_exon_utr {
  my ($self,$transcript_a,$transcript_b) = @_;

  # The first thing to do is to check if the exon from transcript_a is contained in transcript_b. Contained means
  # that the coordinates could match exactly or reside within the donor exon
  my $exon_a = shift(@{$transcript_a->get_all_CDS});
  my $exons_b = $transcript_b->get_all_Exons();
  return 0 if (scalar(@$exons_b) > 1);

  if(scalar(@{$exons_b}) == 1) {
    my $exon_b = shift(@{$exons_b});
    if($exon_a->start == $exon_b->start && $exon_a->end == $exon_b->end) {
      return 0;
    }
  }

  my $contained = 0;

  my $final_exons = [];
  # First add all the exons from the acceptor
  foreach my $exon_b (@{$exons_b}) {
    # If this is true exon a is contained in exon b and we need to create a merged exon
    if($exon_a->start >= $exon_b->start && $exon_a->end <= $exon_b->end) {
      $contained = 1;

      my $merge_exon = new Bio::EnsEMBL::Exon(
                                             -START     => $exon_b->start,
                                             -END       => $exon_b->end,
                                             -STRAND    => $exon_b->strand,
                                             -SLICE     => $exon_b->slice,
                                             -ANALYSIS  => $self->analysis);

      my $supporting_features_a = $exon_a->get_all_supporting_features();
      $merge_exon->add_supporting_features(@{$supporting_features_a});

      my $start_phase;
      my $end_phase;
      my $translation_shift;
      if($exon_a->strand == 1) {
        $translation_shift = $exon_a->start - $exon_b->start;
        # Work out what the end phase should be
        if($exon_b->end > $exon_a->end) {
          $end_phase = -1;
        } else {
          $end_phase = $exon_a->end_phase();
        }
      } else {
        $translation_shift = $exon_b->end - $exon_a->end;
        if($exon_a->start > $exon_b->start) {
          $end_phase = -1;
        } else {
          $end_phase = $exon_a->end_phase();
        }
      }

      # Set the start phase, if there is a shift we know there is 5' UTR so set to -1
      if($translation_shift) {
        $start_phase = -1;
      } else {
        $start_phase = $exon_a->phase();
      }

      push(@{$final_exons},$merge_exon);
    } else {
      my $out_exon = new Bio::EnsEMBL::Exon(
                                           -START  => $exon_b->start,
                                           -END       => $exon_b->end,
                                           -STRAND    => $exon_b->strand,
                                           -SLICE     => $exon_b->slice,
                                           -ANALYSIS  => $self->analysis,
                                           -PHASE     => -1,
                                           -END_PHASE => -1);

      push(@{$final_exons},$out_exon);
    }
  }

  unless($contained) {
# Single exon acceptor was not contained within a donor exon, no UTR will be added
    return(0);
  }

  my $genomic_start;
  my $genomic_end;
  if($transcript_a->strand == 1) {
    $genomic_start = ($transcript_a->translation->start_Exon->seq_region_start + $transcript_a->translation->start - 1);
    $genomic_end = ($transcript_a->translation->end_Exon->seq_region_start + $transcript_a->translation->end - 1);
  } else {
    $genomic_start = ($transcript_a->translation->start_Exon->seq_region_end - $transcript_a->translation->start + 1);
    $genomic_end = ($transcript_a->translation->end_Exon->seq_region_end - $transcript_a->translation->end + 1);
  }

  my $final_translation = create_Translation($final_exons, $genomic_start, $genomic_end);
  foreach my $seq_edit (@{$transcript_a->translation->get_all_SeqEdits}) {
    $final_translation->add_Attributes($seq_edit->get_Attribute);
  }
  $self->add_utr_evidence($final_exons, $exons_b, $transcript_b);
  my $modified_transcript = Bio::EnsEMBL::Transcript->new(-EXONS => $final_exons);
  $modified_transcript->analysis($transcript_a->analysis);
  $modified_transcript->biotype($transcript_a->biotype);
  $modified_transcript->slice($transcript_a->slice());
  $modified_transcript->translation($final_translation);

  calculate_exon_phases($modified_transcript, 0);

  my $modified_translation = $modified_transcript->translation();

  unless($transcript_a->translation->seq eq $modified_transcript->translate->seq && $modified_transcript->translate->seq eq $modified_transcript->translation->seq) {
    $self->throw("There is an issue with the translation after UTR was added. Check above for the sequences, all three should match");
  }

  $transcript_a->{'has_utr'} = 1;
  $modified_transcript->add_supporting_features(grep {$_->isa('Bio::EnsEMBL::DnaDnaAlignFeature')} @{$transcript_b->get_all_supporting_features});
  return($modified_transcript);
}


=head2 join_transcripts

 Arg [1]    : Bio::EnsEMBL::Transcript 5' modified transcript
 Arg [2]    : Bio::EnsEMBL::Transcript 3' modified transcript
 Description: Create a new transcript from Arg[1] and Arg[2] when
              UTRs could be added on both sides
 Returntype : Bio::EnsEMBL::Transcript
 Exceptions : Throws if the translations for Arg[1] and Arg[2] are different
              Throws if the translation changed after merging both transcripts

=cut

sub join_transcripts {
  my ($self,$transcript_a,$transcript_b) = @_;

  my $joined_exon_set = [];
  unless($transcript_a->translation->seq eq $transcript_b->translation->seq) {
    $self->throw("When attempting to join the modified 5' and 3' transcripts, there was a difference in the translation for each. The translation should be ".
                 "identical at this point between the two.\n5' translation:\n".$transcript_a->translation->seq."\n3' translation:\n".$transcript_b->translation->seq);
  }

# Translations from both partial transcripts are identical. Attemting to merge exon sets
  my $exons_a = $transcript_a->get_all_Exons;
  my $exons_b = $transcript_b->get_all_Exons;

  my $cds_5_prime_index = scalar(@{$exons_a}) - scalar(@{$transcript_a->get_all_CDS});
  my $cds_3_prime_index = scalar(@{$transcript_b->get_all_CDS}) - 1;

  foreach my $exon (@{$exons_a}) {
    my $out_exon = new Bio::EnsEMBL::Exon(
                                           -START  => $exon->start,
                                           -END       => $exon->end,
                                           -STRAND    => $exon->strand,
                                           -SLICE     => $exon->slice,
                                           -ANALYSIS  => $self->analysis,
                                           -PHASE     => $exon->phase,
                                           -END_PHASE => $exon->end_phase);
    $out_exon->add_supporting_features(@{$exon->get_all_supporting_features});
    push(@{$joined_exon_set},$out_exon);
  }

  # Remove end exon as this may be modified in 3' transcript
  pop(@{$joined_exon_set});

  for(my $i=$cds_3_prime_index; $i<scalar(@{$exons_b}); $i++) {
    my $exon = $$exons_b[$i];
    my $out_exon = new Bio::EnsEMBL::Exon(
                                           -START  => $exon->start,
                                           -END       => $exon->end,
                                           -STRAND    => $exon->strand,
                                           -SLICE     => $exon->slice,
                                           -ANALYSIS  => $self->analysis,
                                           -PHASE     => $exon->phase,
                                           -END_PHASE => $exon->end_phase);
    $out_exon->add_supporting_features(@{$exon->get_all_supporting_features});
    push(@{$joined_exon_set},$out_exon);
  }

  my $joined_offset = scalar(@{$joined_exon_set}) - scalar(@{$exons_b});
  $cds_3_prime_index += $joined_offset;

  my %seen;
  my @unique_exons;
  foreach my $exon (@$joined_exon_set) {
    push(@unique_exons, $exon) unless (exists $seen{$exon->seq_region_start.':'.$exon->seq_region_end});
    $seen{$exon->seq_region_start.':'.$exon->seq_region_end} = 1;
  }

  my $genomic_start;
  my $genomic_end;
  if($transcript_a->strand == 1) {
    $genomic_start = ($transcript_a->translation->start_Exon->seq_region_start + $transcript_a->translation->start - 1);
    $genomic_end = ($transcript_a->translation->end_Exon->seq_region_start + $transcript_a->translation->end - 1);
  } else {
    $genomic_start = ($transcript_a->translation->start_Exon->seq_region_end - $transcript_a->translation->start + 1);
    $genomic_end = ($transcript_a->translation->end_Exon->seq_region_end - $transcript_a->translation->end + 1);
  }

  # The translation is the same, but still need to modify the translation so that it has the correct start and end exon
  my $translation = create_Translation(\@unique_exons, $genomic_start, $genomic_end);
  foreach my $seq_edit (@{$transcript_a->translation->get_all_SeqEdits}) {
    $translation->add_Attributes($seq_edit->get_Attribute);
  }

  my $joined_transcript = Bio::EnsEMBL::Transcript->new(-EXONS => \@unique_exons);
  $joined_transcript->analysis($transcript_a->analysis);
  $joined_transcript->biotype($transcript_a->biotype);
  $joined_transcript->slice($transcript_a->slice());
  $joined_transcript->translation($translation);

  calculate_exon_phases($joined_transcript, 0);

  unless($transcript_a->translation->seq eq $joined_transcript->translation->seq) {
    $self->throw("When attempting to join the modified 5' and 3' transcripts, there was a difference in the joined translation. The translation should be identical at ".
                 "this point to the original.\nOriginal translation:\n".$transcript_a->translation->seq."\nJoined translation:\n".$joined_transcript->translation->seq);
  }

  $joined_transcript->add_supporting_features(grep {$_->isa('Bio::EnsEMBL::DnaDnaAlignFeature')} @{$transcript_a->get_all_supporting_features});
  $joined_transcript->add_supporting_features(grep {$_->isa('Bio::EnsEMBL::DnaDnaAlignFeature')} @{$transcript_b->get_all_supporting_features});

  return($joined_transcript);
}


=head2 generate_intron_string

 Arg [1]    : Arayref of Bio::EnsEMBL::Intron
 Arg [2]    : Int start of the acceptor genomic region
 Arg [3]    : Int end of the acceptor genomic region
 Description: Create a string representation of the intron structure.
              XXX..XXX: repeated as many times as there are introns
 Returntype : String
 Exceptions : None

=cut

sub generate_intron_string {
  my ($self,$intron_array, $seq_region_start, $seq_region_end) = @_;

my $intron_string = "";
  my $count = 0;
  foreach my $intron (@{$intron_array}) {
    ++$count if ($intron->seq_region_start < $seq_region_end and $intron->seq_region_end > $seq_region_start);
    my $start = $intron->start();
    my $end = $intron->end();
    $intron_string .= $start."..".$end.":";
  }

  return($intron_string, $count);
}


=head2 add_transcript_supporting_features

 Arg [1]    : Bio::EnsEMBL::Transcript the acceptor transcript
 Arg [2]    : Bio::EnsEMBL::Transcript the donor transcript
 Description: Add transcript supporting evidence to Arg[1] based on Arg[2]
 Returntype : Void
 Exceptions : None

=cut

sub add_transcript_supporting_features {
  my ($self,$transcript_a,$transcript_b) = @_;

  my @supporting_features = ();
  # transcript a is the transcript to add to, transcript b is the one with the sfs
  foreach my $sf (@{$transcript_b->get_all_supporting_features()}) {
    if($sf->overlaps_local($transcript_a)) {
        push(@supporting_features,$sf);
    }
  }

  $transcript_a->add_supporting_features(@supporting_features);

  # This code will make sure that intron supporting features are only added for introns that match
  # This is to accommodate cases with the new soft matching code where some introns are not going to match
  my $introns_a = $transcript_a->get_all_Introns();
  my $introns_a_hash = $self->build_feature_string_hash($introns_a);
  my $intron_support_b = $transcript_b->get_all_IntronSupportingEvidence();
  my $intron_support_b_hash = $self->build_feature_string_hash($intron_support_b);
  foreach my $intron_support (@{$intron_support_b}) {
    my $intron_support_string = $intron_support->seq_region_start."..".$intron_support->seq_region_end;
    if($introns_a_hash->{$intron_support_string}) {
      $transcript_a->add_IntronSupportingEvidence($intron_support);
    }
  }
}


=head2 look_for_both

 Arg [1]    : Bio::EnsEMBL::Transcript
 Description: Check that we have a start codon and a stop codon for our translation
              Try to find them if they are not present and Arg[1] has UTRs
              It updates Arg[1] when needed
 Returntype : Bio::EnsEMBL::Transcript
 Exceptions : Throws if the mapping of the start or stop fails

=cut

sub look_for_both {
  my ($self,$trans) = @_;

  my $time = time;
  my $nupdated_start = 0;
  my $nupdated_end = 0;
  my $metcnt = 1;
  my $maxterdist = 150;

    if ($trans->translation) {
      my $tln = $trans->translation;
      my $coding_start = $trans->cdna_coding_start;
      my $orig_coding_start = $coding_start;
      my $cdna_seq = uc($trans->spliced_seq);
      my @pepgencoords = $trans->pep2genomic(1,1);
      if(scalar(@pepgencoords) > 2) {
        goto TLNEND; # I swore I'd never use this - this code desperately needs a rewrite
      }
      my $pepgenstart = $pepgencoords[0]->start;
      my $pepgenend   = $pepgencoords[$#pepgencoords]->end;

      unless($pepgencoords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
        goto TLNEND; # I swore I'd never use this - this code desperately needs a rewrite
      }
      unless($pepgencoords[$#pepgencoords]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
        goto TLNEND; # I swore I'd never use this - this code desperately needs a rewrite
      }

      my $startseq= substr($cdna_seq,$coding_start-1,3);
      if ($startseq ne "ATG") {
        if ($coding_start > 3) {
            my $had_stop = 0;
            while ($coding_start > 3 && !$had_stop) {
                  my $testseq = substr($cdna_seq,$coding_start-4,3);
                  if ($testseq eq "ATG") {
                          print_Translation($trans) if ($self->debug);
                          my @coords = $trans->cdna2genomic($coding_start-3,$coding_start-1);
                          my $new_start;
                          my $new_end;
                          if(scalar(@coords) > 2) {
                            $self->throw("Shouldn't happen - new coding start maps to >2 locations in genome - I'm out of here\n");
                          } elsif (scalar(@coords) == 2) {
                            if ($trans->strand == 1) {
                                $new_start = $coords[0]->start;
                                  $new_end   = $coords[$#coords]->end;
                              } else {
                                  $new_start = $coords[0]->end;
                                    $new_end   = $coords[$#coords]->start;
                                }
                          } else {
                            $new_start = $coords[0]->start;
                            $new_end   = $coords[0]->end;
                          }

                          unless($coords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
                            # Shouldn't happen - new start maps to gap - I'm out of here
                            next;
                          }
                          unless($coords[$#coords]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
                            # Shouldn't happen - new start (end of) maps to gap - I'm out of here
                            next;
                          }

                          if ($new_end - $new_start == 2) {

                            $nupdated_start++;

                            my $newstartexon;
                            foreach my $exon (@{$trans->get_all_Exons}) {
                              if ($exon->seq_region_end >= $new_start && $exon->seq_region_start <= $new_start) {
                                $newstartexon = $exon;
                                last;
                              }
                            }


                            if ($newstartexon == $tln->start_Exon) {
                              if ($tln->start_Exon->strand == 1) {
                                    $tln->start($new_start - $tln->start_Exon->seq_region_start + 1);
                                  } else {
                                        $tln->start($tln->start_Exon->seq_region_end - $new_end + 1);
                                      }

                                # NAUGHTY, but hey I should have to do this - I've changed the translation after all
                                $trans->{'transcript_mapper'} = undef;
                                $trans->{'coding_region_start'} = undef;
                                $trans->{'coding_region_end'} = undef;
                                $trans->{'cdna_coding_start'} = undef;
                                $trans->{'cdna_coding_end'} = undef;

                            } else {
                                # find exon
                              if (!defined($newstartexon)) {
                                        next;
                                  }
                                # create a copy of if and of current start exon (because of phase change)
                              my $copyexon = new Bio::EnsEMBL::Exon(
                                                                    -start  => $tln->start_Exon->start,
                                                                    -end    => $tln->start_Exon->end,
                                                                    -strand => $trans->strand,
                                                                           );
                              my $copynewstartexon = new Bio::EnsEMBL::Exon(
                                                                            -start  => $newstartexon->start,
                                                                            -end    => $newstartexon->end,
                                                                            -strand => $trans->strand,
                                                                                   );

                                # $copyexon->phase(0);
                                $copyexon->end_phase($tln->start_Exon->end_phase);
                                $copyexon->slice($tln->start_Exon->slice);
                              if ($tln->start_Exon->stable_id) {
                                    $copyexon->stable_id($tln->start_Exon->stable_id . "MET" . $metcnt++);
                                        $copyexon->created($time);
                                        $copyexon->modified($time);
                                        $copyexon->version(1);
                                  }

                                $copynewstartexon->phase($newstartexon->phase);
                                # $copynewstartexon->end_phase(0);
                                $copynewstartexon->slice($newstartexon->slice);
                              if ($newstartexon->stable_id) {
                                    $copynewstartexon->stable_id($newstartexon->stable_id . "MET" . $metcnt++);
                                        $copynewstartexon->created($time);
                                        $copynewstartexon->modified($time);
                                        $copynewstartexon->version(1);
                                  }

                                # TODO evidence

                              if ($copynewstartexon->strand == 1) {
                                    $tln->start($new_start - $copynewstartexon->seq_region_start + 1);
                                  } else {
                                        $tln->start($copynewstartexon->seq_region_end - $new_end + 1);
                                      }

                                # Replace exons in transcript, and fix phases

                                my @newexons;
                                my $inrange = 0;
                              foreach my $exon (@{$trans->get_all_Exons}) {
                                if ($inrange) {
                                        $exon->phase( $newexons[$#newexons]->end_phase );
                                              $exon->end_phase(($exon->length + $exon->phase) % 3);
                                      }
                                if ($exon == $tln->start_Exon) {
                                        $copyexon->phase( $newexons[$#newexons]->end_phase );

                                              push @newexons,$copyexon;
                                              $inrange = 0;
                                      } elsif ($exon == $newstartexon) {
                                              push @newexons,$copynewstartexon;
                                                    $copynewstartexon->end_phase(($exon->length - $tln->start + 1)%3);
                                                    $inrange = 1;
                                            } else {
                                                    push @newexons,$exon;
                                                  }
                              }


                                $trans->flush_Exons;
                              foreach my $exon (@newexons) {
                                    $trans->add_Exon($exon);
                                  }

                                # Reset translation start exon
                              if ($tln->end_Exon == $tln->start_Exon) {
                                    $tln->end_Exon($copyexon);
                                  }
                                $tln->start_Exon($copynewstartexon);

                            }
                          }

                                last;
                        } else {
                          if ($testseq =~ /TAA/ or $testseq =~ /TGA/ or $testseq =~ /TAG/) {
                            $had_stop = 1;
                          } else {
                            $coding_start -= 3;
                          }
                        }
                }
          }
      }

    TLNEND:
      {
        my $coding_end = $trans->cdna_coding_end;
        my $orig_coding_end = $coding_end;


        my $peplen = $trans->translate->length;

        my @pepgencoords = $trans->pep2genomic($peplen,$peplen);

        if(scalar(@pepgencoords) > 2) {
              next;
          }

        my $pepgenstart = $pepgencoords[0]->start;
        my $pepgenend   = $pepgencoords[$#pepgencoords]->end;

        unless($pepgencoords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
          # pep end maps to gap
              next;
          }
        unless($pepgencoords[$#pepgencoords]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
          # pep start (end of) maps to gap\n";
              next;
          }

        my $endseq= substr($cdna_seq,$coding_end-3,3);
        my $cdnalen = length($cdna_seq);

        my $longendseq= substr($cdna_seq,$coding_end-6,12);

        if ($endseq ne "TGA" and $endseq ne "TAA" and $endseq ne "TAG") {
          if (($cdnalen-$coding_end) > 3) {
            while (($cdnalen-$coding_end) > 0 && ($coding_end-$orig_coding_end) <= $maxterdist) {
                    my $testseq = substr($cdna_seq,$coding_end,3);

                    if ($testseq eq "TGA" or $testseq eq "TAA" or $testseq eq "TAG") {

                      my @coords = $trans->cdna2genomic($coding_end+1,$coding_end+3);
                      my $new_start;
                      my $new_end;
                      if(scalar(@coords) > 2) {
                          $self->throw("new end does not map cleanly\n");
                        } elsif (scalar(@coords) == 2) {
                            if ($trans->strand == 1) {
                                  $new_start = $coords[0]->start;
                                      $new_end   = $coords[$#coords]->end;
                                } else {
                                      $new_start = $coords[0]->end;
                                          $new_end   = $coords[$#coords]->start;
                                    }
                          } else {
                              $new_start = $coords[0]->start;
                                $new_end   = $coords[0]->end;
                            }

                      unless($coords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
                         # new start maps to gap
                            next;
                        }
                      unless($coords[$#coords]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
                         # new start (end of) maps to gap\n";
                            next;
                        }

                      if ($new_end - $new_start == 2) {

                          $nupdated_end++;

                            my $newendexon;
                          foreach my $exon (@{$trans->get_all_Exons}) {
                            if ($exon->seq_region_end >= $new_start && $exon->seq_region_start <= $new_start) {
                                    $newendexon = $exon;
                                          last;
                                  }
                          }

                          if ($newendexon == $tln->end_Exon) {
                            if ($tln->end_Exon->strand == 1) {
                                    $tln->end($new_end - $tln->end_Exon->seq_region_start + 1);
                                  } else {
                                          $tln->end($tln->end_Exon->seq_region_end - $new_start + 1);
                                        }

                                # NAUGHTY, but hey I should have to do this - I've changed the translation after all
                                $trans->{'transcript_mapper'} = undef;
                                $trans->{'coding_region_start'} = undef;
                                $trans->{'coding_region_end'} = undef;
                                $trans->{'cdna_coding_start'} = undef;
                                $trans->{'cdna_coding_end'} = undef;

                          } else {
                                # find exon
                            if (!defined($newendexon)) {
                                   # Failed finding new end exon - how can this be?
                                          next;
                                  }
                                # create a copy of if and of current end exon (because of phase change)
                            my $copyexon = new Bio::EnsEMBL::Exon(
                                                                    -start  => $tln->end_Exon->start,
                                                                    -end    => $tln->end_Exon->end,
                                                                    -strand => $trans->strand,
                                                                   );
                            my $copynewendexon = new Bio::EnsEMBL::Exon(
                                                                        -start  => $newendexon->start,
                                                                        -end    => $newendexon->end,
                                                                        -strand => $trans->strand,
                                                                               );

                                $copyexon->phase($tln->end_Exon->phase);
                                $copyexon->end_phase($tln->end_Exon->end_phase);
                                $copyexon->slice($tln->end_Exon->slice);
                            if ($tln->end_Exon->stable_id) {
                                    $copyexon->stable_id($tln->end_Exon->stable_id . "TER" . $metcnt++);
                                          $copyexon->created($time);
                                          $copyexon->modified($time);
                                          $copyexon->version(1);
                                  }

                                $copynewendexon->phase($newendexon->phase);
                                # $copynewendexon->end_phase(0);
                                $copynewendexon->slice($newendexon->slice);
                            if ($newendexon->stable_id) {
                                    $copynewendexon->stable_id($newendexon->stable_id . "TER" . $metcnt++);
                                          $copynewendexon->created($time);
                                          $copynewendexon->modified($time);
                                          $copynewendexon->version(1);
                                  }

                                # TODO evidence

                            if ($copynewendexon->strand == 1) {
                                    $tln->end($new_end - $copynewendexon->seq_region_start + 1);
                                  } else {
                                          $tln->end($copynewendexon->seq_region_end - $new_start + 1 );

                                                my $tercodon = $copynewendexon->seq->subseq($copynewendexon->seq_region_end - $new_start-1, $copynewendexon->seq_region_end - $new_start +1);
                                                #reverse($tercodon);
                                                #$tercodon =~ tr /ACGT/TGCA/;

                                        }

                                # Replace exons in transcript, and fix phases
                                my @newexons;
                                my $inrange = 0;
                            foreach my $exon (@{$trans->get_all_Exons}) {
                              if ($inrange) {
                                $exon->phase( $newexons[$#newexons]->end_phase );
                                $exon->end_phase(($exon->length + $exon->phase) % 3);
                              }
                              if ($exon == $tln->end_Exon) {
                                my $phase = $exon->phase;
                                if ($phase == -1) {
                                    $phase = 0;
                                  }
                                if ($exon == $tln->start_Exon) {
                                    $copyexon->end_phase(($exon->length - $tln->start + 1)%3);
                                  } else {
                                      $copyexon->end_phase(($exon->length + $exon->phase)%3);
                                    }

                                push @newexons,$copyexon;
                                $inrange = 1;
                              } elsif ($exon == $newendexon) {
                                $copynewendexon->phase( $newexons[$#newexons]->end_phase );
                                $copynewendexon->end_phase( -1);

                                push @newexons,$copynewendexon;
                                $inrange = 0;
                              } else {
                                push @newexons,$exon;
                              }
                            }

                                $trans->flush_Exons;
                            foreach my $exon (@newexons) {
                                    $trans->add_Exon($exon);
                                  }

                                # Reset translation start exon
                            if ($tln->end_Exon == $tln->start_Exon) {
                                    $tln->start_Exon($copyexon);
                                  }
                                $tln->end_Exon($copynewendexon);

                          }
                        }
                      last;
                    }
                          $coding_end += 3;
                  }
          }
        }
      }
    }

  return($trans);
}


=head2 filter_input_genes

 Arg [1]    : Arrayref of hashref representing the DB connection details
 Arg [2]    : Hashref, the keys should be logic_names and the values 1 or
                an arrayref of biotypes. If the hashref is empty, fetch all
                the genes
 Arg [3]    : Boolean, update the biotype for any donor transcript based on
              the name of the database, cDNAs and long reads are 'cdna_donor',
              short reads are 'rnaseq_donor'
 Description: Fetch all the genes for each databases in Arg[1] when Arg[2]
              is undef or selection of genes based on Arg[2]
 Returntype : Arrayref of Bio::EnsEMBL::Gene
 Exceptions : Throws if the dbname is not recognised when using Arg[3]

=cut

sub filter_input_genes {
  my ($self, $gene_source_dbs, $allowed_transcript_sets,$standardise_biotypes) = @_;

  my @genes;
  my $slice = $self->query;

  foreach my $db_conn (@$gene_source_dbs) {
    my $db_adaptor = $self->hrdb_get_dba($db_conn);
    my $gene_adaptor = $db_adaptor->get_GeneAdaptor();
    my $dbname = $db_conn->{'-dbname'};

    $self->say_with_header("Fetching from $dbname");
    if($allowed_transcript_sets) {
      foreach my $logic_name (keys %$allowed_transcript_sets) {
        if (ref($allowed_transcript_sets->{$logic_name}) eq 'ARRAY') {
          foreach my $biotype (${$allowed_transcript_sets->{$logic_name}}) {
            push(@genes, @{$gene_adaptor->fetch_all_by_Slice($slice, $logic_name, 1, undef, $biotype)});
          }
        }
        else {
          push(@genes, @{$gene_adaptor->fetch_all_by_Slice($slice, $logic_name, 1)});
        }
      }
    }
    else {
#      my $donor_genes = $gene_adaptor->fetch_all_by_Slice($slice, undef, 1);
      my $donor_genes = $gene_adaptor->fetch_all_by_Slice($slice);
      # This should not have to be done in general, however there is a fix because
      # we have a very large number of assemblies in production and the pipelines
      # will not work without this fix
      if($standardise_biotypes) {
        my $biotype = "";
        if($dbname =~ /\_cdna\_/ || $dbname =~ /\_lrfinal\_/) {
          $biotype = "cdna_donor";
        } elsif($dbname =~ /\_rnalayer\_/ || $dbname =~ /\_rnaseq\_/ || $dbname =~ /\_star\_rs\_layer\_/ || $dbname =~ /\_pcp\_/) {
          $biotype = "rnaseq_donor";
        } else {
          $self->throw("Found an unexpected dbname type for the donor db. Name: ".$dbname);
        }
        foreach my $gene (@$donor_genes) {
          eval {
            $gene->load();
          };
          if ($@) {
            $self->warning("Skipping ".$gene->display_id." from $dbname");
            next;
          };
          $gene->biotype($biotype);
          my $remove_gene = 0;
          my $transcripts = $gene->get_all_Transcripts();
          foreach my $transcript (@$transcripts) {
            $transcript->biotype($biotype);
            eval {
              $transcript->translation;
            };
            if ($@) {
              $remove_gene = 1;
            }
          }
          if ($remove_gene) {
            next;
          }
        }
        push(@genes, @{$donor_genes});
      } else {
        my $direct_output_genes = [];
        foreach my $gene (@$donor_genes) {
          if($gene->biotype eq 'cdna' || $gene->biotype eq 'rnaseq_tissue' || $gene->biotype eq 'rnaseq_merged') {
            push(@$direct_output_genes,$gene);
          } else {
            push(@genes,$gene);
	  }
	}
        if(scalar(@$direct_output_genes)) {
          $self->output($direct_output_genes);
	}
      } # end inner else
    } # end outer else
    $db_adaptor->dbc->disconnect_if_idle();
  }
  return \@genes;
}


=head2 trim_5prime_utr_short_read

 Arg [1]    : Bio::EnsEMBL::Transcript, the transcript with a 5' UTR to trim
 Description: If the 5' UTR of the transcript is longer than the value in '_max_5prime_utr_size'
              we remove all UTR exons until the UTR is equal or smaller than the maximum value.
              If the first UTR exon is more than double the maximum UTR size, it is trimmed down
              to the maximum UTR size. Otherwise the exon is left untouched.
 Returntype : Bio::EnsEMBL::Transcript
 Exceptions : None

=cut

sub trim_5prime_utr_short_read {
  my ($self, $transcript) = @_;

  my $max_utr_size = $self->param('_max_5prime_utr_size');
  my $relative_coding_start = $transcript->cdna_coding_start;
  if ($relative_coding_start > $max_utr_size) {
    $self->say_with_header("Trimming 5' UTR for ".$transcript->display_id.' '.$transcript->biotype.' '.$transcript->source);
    my $exons = $transcript->get_all_Exons;
    $transcript->flush_Exons;
    foreach my $exon (@$exons) {
      if ($relative_coding_start-$exon->length > $max_utr_size) {
        $relative_coding_start -= $exon->length;
        $self->say_with_header("Removing ".$exon->seq_region_start.' '.$exon->seq_region_end);
      }
      else {
        if ($relative_coding_start > $max_utr_size and $exon->length > $max_utr_size*2 ) {
          my $diff = $exon->length+$max_utr_size-$relative_coding_start-1;
          $self->say_with_header("Long exon ".$exon->seq_region_start.' '.$exon->seq_region_end.' '.$exon->length.' '.$diff);
          $relative_coding_start -= $exon->length;
          if ($diff < 0 or $diff > $exon->length) {
            $self->say_with_header('Removing '.$exon->seq_region_start.' '.$exon->seq_region_end.' because it will create a small exon');
            next;
          }
          if ($relative_coding_start < 1) {
            $self->say_with_header("Reset translation ".$exon->seq_region_start.' '.$exon->seq_region_end.' '.$transcript->translation->start);
            my $tln = $transcript->translation->new(
              -start_exon => $transcript->translation->start_Exon,
              -end_exon => $transcript->translation->end_Exon,
              -seq_start => $max_utr_size,
              -seq_end => $transcript->translation->end,
            );
            $transcript->translation($tln);
          }
          if ($exon->strand == 1) {
            $exon->start($exon->end-$diff);
          }
          else {
            $exon->end($exon->start+$diff);
          }
        }
        else {
          $relative_coding_start -= $exon->length;
        }
        $transcript->add_Exon($exon);
      }
    }
    my $ises = $transcript->get_all_IntronSupportingEvidence;
    $transcript->flush_IntronSupportingEvidence;
    foreach my $ise (@$ises) {
      if ($ise->start <= $transcript->end and $ise->end >= $transcript->start) {
        $transcript->add_IntronSupportingEvidence($ise);
      }
    }
  }
  else {
    $self->say_with_header("5' UTR is shorter than $max_utr_size, no need to trim");
  }
  return $transcript;
}


=head2 trim_3prime_utr_short_read

 Arg [1]    : Bio::EnsEMBL::Transcript, the transcript with a 3' UTR to trim
 Description: Looking at the last exon with coding sequence, we try to find the polyA signal
              AATAAA. If it's not present we have a hard limit at 1,000bp.
              When the polyA signla is present we are looking for a cleavage signal around
              15-30 bp downstream. If present we use its position or we are trimming at 15bp.
 Returntype : Bio::EnsEMBL::Transcript
 Exceptions : None

=cut

sub trim_3prime_utr_short_read {
  my ($self,$transcript) = @_;

  $self->say_with_header("Trimming 3' UTR for ".$transcript->display_id.' '.$transcript->biotype.' '.$transcript->source);
  # Only going to use the most common nuclear PAS signal
  my $pas_signal = 'AATAAA';
  my $cleavage_signal = 'CA'; # Not sure about that one
  my $new_end = 1000; # If there is no signal/cleavage, the max 3' end is 1000bp

  my $three_prime_utr = $transcript->three_prime_utr;
  if ($three_prime_utr) {
    my $exons = $transcript->get_all_Exons;
    my $stop_codon_exon = $transcript->translation->end_Exon;
    my $stop_codon_seq = $stop_codon_exon->seq;
    my $translation_end = $transcript->translation->end;
    my $ises = $transcript->get_all_IntronSupportingEvidence;
    $transcript->flush_Exons;
    my $polya_pos = index($stop_codon_seq->seq, $pas_signal);
    if ($polya_pos == -1 or $polya_pos < $translation_end) {
      $self->say_with_header("Did not found $pas_signal in last coding exon");
    }
    else {
      $self->say_with_header("Found $pas_signal at $polya_pos");
      my $cleavage_pos = index(substr($stop_codon_seq->seq, $polya_pos+14, 15), $cleavage_signal);
      if ($cleavage_pos == -1) {
        $self->say_with_header("Did not found $cleavage_signal after $polya_pos");
        $cleavage_pos = 15;
      }
      my $new_end = $polya_pos+$cleavage_pos;
    }
    foreach my $exon (@$exons) {
      if ($exon == $stop_codon_exon) {
        if ($exon->length > $new_end) {
          if ($exon->length-$translation_end > $new_end) {
            if ($exon->strand == 1) {
              $exon->end($exon->start+$translation_end+$new_end);
            }
            else {
              $exon->start($exon->end-$translation_end-$new_end);
            }
          }
        }
        else {
          if ($new_end < $exon->length) {
            if ($exon->strand == 1) {
              $exon->end($exon->start+$new_end);
            }
            else {
              $exon->start($exon->end-$new_end);
            }
          }
        }
        $transcript->add_Exon($exon);
        last;
      }
      $transcript->add_Exon($exon);
    }
    if ($stop_codon_exon != $exons->[-1]) {
      $transcript->flush_IntronSupportingEvidence;
      foreach my $ise (@$ises) {
        if ($ise->start <= $transcript->end and $ise->end >= $transcript->start) {
          $transcript->add_IntronSupportingEvidence($ise);
        }
      }
    }
  }
  else {
    $self->warning("The trim_3prime_utr_short_read was called on a transcript with no 3 prime UTR. Nothing to trim");
  }

  return $transcript;
}

=head2 acceptor_genes

 Arg [1]    : Arrayref of Bio::EnsEMBL::Gene
 Description: Getter/setter for the list of acceptor genes
 Returntype : Arrayref of Bio::EnsEMBL::Gene
 Exceptions : None

=cut

sub acceptor_genes {
  my ($self, $val) = @_;

  if($val) {
    $self->param('_acceptor_genes',$val);
  }

  return $self->param('_acceptor_genes');
}


=head2 donor_genes

 Arg [1]    : Arrayref of Bio::EnsEMBL::Gene
 Description: Getter/setter for the list of acceptor genes
 Returntype : Arrayref of Bio::EnsEMBL::Gene
 Exceptions : None

=cut

sub donor_genes {
  my ($self, $val) = @_;

  if($val) {
    $self->param('_donor_genes',$val);
  }

  return $self->param('_donor_genes');
}


=head2 biotype_priorities

 Arg [1]    : None
 Description: Returns the 'utr_biotype_priorities' hash from the parameters
 Returntype : Hashref
 Exceptions : None

=cut

sub biotype_priorities {
  my ($self,$biotype) = @_;

  # 1 = cdna, 2 = rnaseq, 3 = other
  if (index($biotype, 'cdna') == 0) {
   $biotype = 'cdna';
  } elsif (index($biotype, 'rnaseq') == 0) {
    $biotype = 'rnaseq';
  } else {
    $biotype = 'other';
  }
  return($self->param('utr_biotype_priorities')->{$biotype});
}


=head2 min_size_5prime

 Arg [1]    : None
 Description: Getter for the minimum 5' exon size. It can be set for the
              analysis with 'min_size_5prime'
 Returntype : Int
 Exceptions : None

=cut

sub min_size_5prime {
  my ($self) = @_;

  return $self->param('min_size_5prime');
}



=head2 find_soft_match

 Arg [1]    : A transcript ref for a transcript that you want to add UTR onto
 Arg [2]    : An arrayref of potential UTR donor transcripts
 Description: This method will attempt to add UTR to an acceptor transcript from a list of donor transcripts. It begins by sorting the donor
              transcripts into layers based on their biotype priorities (cdna/IsoSeq being top, with short read data next and other stuff after that)
              For each layer it tries to prioritise donors that match to both the terminal 5' and 3' introns over single end matches. In cases where
              it has to decide between matches, e.g. two donors that match both terminal introns, it considers the following:
              1) Pick one with the most total introns in common with the acceptor transcript
                2) If above is equal, pick the one with the most splicing complexity in the UTR
                  3) If the above are equal pick the one with the longest UTR
                    4) If all the above are tied it just sticks with the first of the pair that it found
              The same logic holds true for 5' and 3' only matches
              In the case where there are single end only matches, it will priortise by the layer first, so if the best match on the 5' end is a cDNA
              this will stay as the best overall match, even if the best 3' end match comes from a lower priority layer such as the short read data.
              This is somewhat important as it means that even if there is a match to both ends in a lower layer, that this will not override a single
              end match in a higher layer

 Returntype : A transcript ref that may or may not have UTR added
 Exceptions : None

=cut

sub find_soft_match {
  my ($self,$acceptor_transcript,$donor_transcripts) = @_;

  # Sort donors by priority first. Putting in to separate layers to begin with
  my $donor_layers = [];
  foreach my $donor_transcript (@$donor_transcripts) {
    my $priority = $self->biotype_priorities($donor_transcript->biotype);
    unless(${$donor_layers}[$priority-1]) {
      ${$donor_layers}[$priority-1] = [];
    }
    push(@{${$donor_layers}[$priority-1]},$donor_transcript);
  }

  # Now donor_layers will have arrayrefs of transcripts that sit at different priorities, with elemenet 0 having
  # the highest priority and the last element having the lowest priority transcripts
  # Loop through donor layers to see if a match can be found on a layer. As soon as a match is found on a layer then stop
  # If multiple matches are found on a layer we must prioritise based on the information content
  # Note the priority is just set to an arbritarily high number, it needs to be something above the levels of priority
  my $best_5_prime_match_transcript;
  my $best_5_prime_match_score = 0;
  my $best_5_prime_match_priority = 777;
  my $best_3_prime_match_transcript;
  my $best_3_prime_match_score = 0;
  my $best_3_prime_match_priority = 777;
  my $best_both_match_transcript;
  my $best_both_match_score = 0;

  # As this is a CDS only transcript, calling get_all_Introns is fine. Calling get_all_CDS_Introns resulted in extreme memory usage as there
  # is some sort of memory leak in the core API
  my $acceptor_introns = $acceptor_transcript->get_all_Introns();

  my $acceptor_introns_count = scalar(@$acceptor_introns);
  unless($acceptor_introns_count) {
    $self->throw('Went into soft intron matching code for single exon acceptor, something went wrong');
  }
  my $acceptor_introns_coord_hash = $self->build_feature_string_hash($acceptor_introns);

  my ($acceptor_cds_intron_string) = $self->generate_intron_string($acceptor_introns,$acceptor_transcript->coding_region_start,$acceptor_transcript->coding_region_end);

  foreach my $layer_transcripts (@$donor_layers) {
    foreach my $donor_transcript (@$layer_transcripts) {
      my $donor_priority = $self->biotype_priorities($donor_transcript->biotype);
      my $donor_introns = $donor_transcript->get_all_Introns();
      my $donor_introns_count = scalar(@$donor_introns);
      my $donor_introns_coord_hash = $self->build_feature_string_hash($donor_introns);
      unless($donor_introns_count) {
        next;
      }

      my ($donor_intron_string) = $self->generate_intron_string($donor_introns,$acceptor_transcript->coding_region_start,$acceptor_transcript->coding_region_end);
      # Want to locate the position of the terminal introns on each side of the acceptor in the donor (if present)
      my $five_prime_match_score = 0;
      my $three_prime_match_score = 0;
      my $both_match_score = 0;

      my $acceptor_terminal_5_prime_intron_start = ${$acceptor_introns}[0]->seq_region_start;
      my $acceptor_terminal_5_prime_intron_end = ${$acceptor_introns}[0]->seq_region_end;
      my $acceptor_terminal_5_prime_string = $acceptor_terminal_5_prime_intron_start."..".$acceptor_terminal_5_prime_intron_end;
      my $acceptor_terminal_3_prime_intron_start = ${$acceptor_introns}[$acceptor_introns_count-1]->seq_region_start;
      my $acceptor_terminal_3_prime_intron_end = ${$acceptor_introns}[$acceptor_introns_count-1]->seq_region_end;
      my $acceptor_terminal_3_prime_string = $acceptor_terminal_3_prime_intron_start."..".$acceptor_terminal_3_prime_intron_end;

      my $total_score = 0;
      my $matched_5_prime = 0;
      my $matched_3_prime = 0;
      foreach my $acceptor_intron_string (keys(%$acceptor_introns_coord_hash)) {
        if($donor_introns_coord_hash->{$acceptor_intron_string}) {
          $total_score++;
          if($acceptor_intron_string eq $acceptor_terminal_5_prime_string) {
            $matched_5_prime = 1;
          }

          if($acceptor_intron_string eq $acceptor_terminal_3_prime_string) {
            $matched_3_prime = 1;
          }
        }
      } # End foreach my $acceptor_intron_string

      my $five_prime_merged_transcript;
      my $three_prime_merged_transcript;
      my $both_ends_merged_transcript;
      if($matched_5_prime) {
        $five_prime_merged_transcript = $self->add_five_prime_utr($acceptor_transcript,$donor_transcript,$acceptor_introns,$donor_introns,$acceptor_cds_intron_string,$donor_intron_string);
      }

      if($matched_3_prime) {
        $three_prime_merged_transcript = $self->add_three_prime_utr($acceptor_transcript,$donor_transcript,$acceptor_introns,$donor_introns,$acceptor_cds_intron_string,$donor_intron_string);
      }

      if($five_prime_merged_transcript && $three_prime_merged_transcript) {
        $both_ends_merged_transcript = $self->join_transcripts($five_prime_merged_transcript,$three_prime_merged_transcript);
      }

      # At this point we know the total score and whether it matched the 5' terminal intron or the 3' terminal intron or both or none
      # If it matched none, then we move on
      # If it matched both then this should get the highest priority, if it is better than the current best
      # If it matched either end then it should be stored as a potential donor, if it is better than the current best end only match
      # First handle the case that it matches both ends, as we want to give preference to scenarios where the donor data comes from a single source
      # Note that this will not override a single end match from a higher priority. So if there is a higher priority existing 5' or 3' match then this bit
      # will not execute. Instead it will go into the else and consider the ends individually. This scenario would crop up if we had a cDNA style support
      # for a 5' UTR but then got a short read match on both ends. The 5' end would keep the cDNA support and the 3' could get the short read support
      if($both_ends_merged_transcript && $donor_priority <= $best_5_prime_match_priority && $donor_priority <= $best_3_prime_match_priority) {
        # If it's valid and the total score is better, then it is the new best
        # Elsif it's valid and has the same score it's better if:
        #   1) It has more combined 5' and 3' UTR features (i.e. higher complexity)
        #   2) It has the same number of features, but the combined length of the 5' and 3' UTR is longer
        # Note that the above rules are just one way of doing this, the choice of complexity over length is debateable
        if($total_score > $best_both_match_score) {
          $best_both_match_score = $total_score;
          $best_both_match_transcript = $both_ends_merged_transcript;
          $best_both_match_transcript->{'donor_5prime_biotype'} = $donor_transcript->biotype;
          $best_both_match_transcript->{'donor_3prime_biotype'} = $donor_transcript->biotype;
        } elsif($total_score == $best_both_match_score) {
          if((scalar(@{$both_ends_merged_transcript->get_all_five_prime_UTRs}) + scalar(@{$both_ends_merged_transcript->get_all_three_prime_UTRs})) >
             (scalar(@{$best_both_match_transcript->get_all_five_prime_UTRs}) + scalar(@{$best_both_match_transcript->get_all_three_prime_UTRs}))) {
            $best_both_match_transcript = $both_ends_merged_transcript;
            $best_both_match_transcript->{'donor_5prime_biotype'} = $donor_transcript->biotype;
            $best_both_match_transcript->{'donor_3prime_biotype'} = $donor_transcript->biotype;
	  } elsif((scalar(@{$both_ends_merged_transcript->get_all_five_prime_UTRs}) + scalar(@{$both_ends_merged_transcript->get_all_three_prime_UTRs})) ==
                  (scalar(@{$best_both_match_transcript->get_all_five_prime_UTRs}) + scalar(@{$best_both_match_transcript->get_all_three_prime_UTRs}))) {
            if((($both_ends_merged_transcript->five_prime_utr->length) + ($both_ends_merged_transcript->three_prime_utr->length)) >
               (($best_both_match_transcript->five_prime_utr->length) + ($best_both_match_transcript->three_prime_utr->length))) {
              $best_both_match_transcript = $both_ends_merged_transcript;
              $best_both_match_transcript->{'donor_5prime_biotype'} = $donor_transcript->biotype;
              $best_both_match_transcript->{'donor_3prime_biotype'} = $donor_transcript->biotype;
	    }
	  }
	} # End elsif($valid && ($total_score == $best_both_match_score))
      } # End if($matched_5_prime && $matched_3_prime)

      # In this case we do not have a match to both ends, we want to identify the best 5' and/or 3' match. It is basically the same as the above
      # but only considering the data from either end when picking the current best. It also only needs to validate the match on one end
      # Makes the same assumption that complexity of splicing is better than total UTR length when choosing UTR
      else {
        if($five_prime_merged_transcript) {
          if($donor_priority > $best_5_prime_match_priority) {
          } elsif($total_score > $best_5_prime_match_score && $donor_priority < $best_5_prime_match_priority ) {
            $best_5_prime_match_score = $total_score;
            $best_5_prime_match_priority = $donor_priority;
            $best_5_prime_match_transcript = $five_prime_merged_transcript;
            $best_5_prime_match_transcript->{'donor_5prime_biotype'} = $donor_transcript->biotype;
          } elsif($total_score == $best_5_prime_match_score) {
            if(scalar(@{$five_prime_merged_transcript->get_all_five_prime_UTRs}) > scalar(@{$best_5_prime_match_transcript->get_all_five_prime_UTRs})) {
              $best_5_prime_match_priority = $donor_priority;
              $best_5_prime_match_transcript = $five_prime_merged_transcript;
              $best_5_prime_match_transcript->{'donor_5prime_biotype'} = $donor_transcript->biotype;
	    } elsif(scalar(@{$five_prime_merged_transcript->get_all_five_prime_UTRs}) == scalar(@{$best_5_prime_match_transcript->get_all_five_prime_UTRs})) {
              if($five_prime_merged_transcript->five_prime_utr->length > $best_5_prime_match_transcript->five_prime_utr->length) {
                $best_5_prime_match_priority = $donor_priority;
                $best_5_prime_match_transcript = $five_prime_merged_transcript;
                $best_5_prime_match_transcript->{'donor_5prime_biotype'} = $donor_transcript->biotype;
              }
            }
          } # End elsif($total_score == $best_5_prime_match_score)
        }

        if($three_prime_merged_transcript) {
          if($donor_priority > $best_3_prime_match_priority) {
	  } elsif($total_score > $best_3_prime_match_score) {
            $best_3_prime_match_score = $total_score;
            $best_3_prime_match_priority = $donor_priority;
            $best_3_prime_match_transcript = $three_prime_merged_transcript;
            $best_3_prime_match_transcript->{'donor_3prime_biotype'} = $donor_transcript->biotype;
          } elsif($total_score == $best_3_prime_match_score) {
            if(scalar(@{$three_prime_merged_transcript->get_all_three_prime_UTRs}) > scalar(@{$best_3_prime_match_transcript->get_all_three_prime_UTRs})) {
              $best_3_prime_match_priority = $donor_priority;
              $best_3_prime_match_transcript = $three_prime_merged_transcript;
              $best_3_prime_match_transcript->{'donor_3prime_biotype'} = $donor_transcript->biotype;
	    } elsif(scalar(@{$three_prime_merged_transcript->get_all_three_prime_UTRs}) == scalar(@{$best_3_prime_match_transcript->get_all_three_prime_UTRs})) {
              if($three_prime_merged_transcript->three_prime_utr->length > $best_3_prime_match_transcript->three_prime_utr->length) {
                $best_3_prime_match_priority = $donor_priority;
                $best_3_prime_match_transcript = $three_prime_merged_transcript;
                $best_3_prime_match_transcript->{'donor_3prime_biotype'} = $donor_transcript->biotype;
              }
            }
          } # End elsif($total_score == $best_5_prime_match_score)
        } # End if($three_prime_merged_transcript)
      } # End else
    } # End foreach my $donor_transcript

    # At this point we look to see if we had either one best match, two single end matches, one single end match or no match
    # Depending on what we had we either do 3 merges (from two separate single end, then a final joint merge), one merge (from one best match or
    # one single end) or we don't merge anything cos we had no match
    my $final_transcript;
    if($best_both_match_transcript) {
      $final_transcript = $best_both_match_transcript;
      return($final_transcript);
    } elsif($best_5_prime_match_transcript && $best_3_prime_match_transcript) {
 #     $final_transcript = $best_5_prime_match_transcript;
      $final_transcript = $self->join_transcripts($best_5_prime_match_transcript,$best_3_prime_match_transcript);
      $final_transcript->{'donor_5prime_biotype'} = $best_5_prime_match_transcript->{'donor_5prime_biotype'};
      $final_transcript->{'donor_3prime_biotype'} = $best_3_prime_match_transcript->{'donor_3prime_biotype'};
      return($final_transcript);
    }
  } # End foreach my $layer_transcripts

  # If we get to this point without anything having been returned it means we have either 5' or 3' or no UTR
  if($best_5_prime_match_transcript) {
    return($best_5_prime_match_transcript);
  } elsif($best_3_prime_match_transcript) {
    return($best_3_prime_match_transcript);
  } else {
    return(0);
  }
}


sub build_feature_string_hash {
  my ($self,$features) = @_;

  my $feature_string_hash = {};
  foreach my $feature (@$features) {
    my $start = $feature->seq_region_start;
    my $end = $feature->seq_region_end;
    $feature_string_hash->{$start."..".$end} = 1;
  }

  return($feature_string_hash);
}


sub print_gene_details {
  my ($self,$gene) = @_;

  say "Gene:\n".$gene->seq_region_start."..".$gene->seq_region_end.":".$gene->strand;
  my $transcripts = $gene->get_all_Transcripts;
  foreach my $transcript (@$transcripts) {
    say "  Transcript:\n  ".$transcript->seq_region_start."..".$transcript->seq_region_end.":".$transcript->strand;
    my $exons = $transcript->get_all_Exons();
    my $exon_string = "";
    foreach my $exon (@$exons) {
      $exon_string .= '('.$exon->seq_region_start."..".$exon->seq_region_end.")..";
    }
    chop($exon_string);
    chop($exon_string);
    say "    Exons:\n    ".$exon_string;
  }
}

sub print_transcript_details {
  my ($self,$transcript) = @_;

  say "Transcript:\n".$transcript->seq_region_start."..".$transcript->seq_region_end.":".$transcript->strand;
  my $exons = $transcript->get_all_Exons();
  my $exon_string = "";
  foreach my $exon (@$exons) {
    $exon_string .= '('.$exon->seq_region_start."..".$exon->seq_region_end.")..";
  }
  chop($exon_string);
  chop($exon_string);
  say "    Exons:\n    ".$exon_string;
}


sub create_cds_transcript {
  my ($self,$transcript) = @_;

  my $cds_transcript = Bio::EnsEMBL::Transcript->new();
  my $cds_exons = $transcript->get_all_translateable_Exons;
  foreach my $cds_exon (@$cds_exons) {
    $cds_transcript->add_Exon($cds_exon);
  }

  $cds_transcript->biotype($transcript->biotype);
  $cds_transcript->analysis($transcript->analysis);
  $cds_transcript->slice($transcript->slice);
  my $start_exon = ${$cds_exons}[0];
  my $end_exon = ${$cds_exons}[scalar(@$cds_exons)-1];
  my $translation = Bio::EnsEMBL::Translation->new();
  $translation->start_Exon($start_exon);
  $translation->start(1);
  $translation->end_Exon($end_exon);
  $translation->end($end_exon->length());
  foreach my $seq_edit (@{$transcript->translation->get_all_SeqEdits}) {
    $translation->add_Attributes($seq_edit->get_Attribute);
  }
  $cds_transcript->translation($translation);
  calculate_exon_phases($cds_transcript, 0);

  $cds_transcript->{'5_prime_utr'} = 0;
  $cds_transcript->{'3_prime_utr'} = 0;

  return($cds_transcript);
}


sub remove_redundant_acceptors {
  my ($self,$genes) = @_;

  my $output_genes = [];
  my $transcript_strings = {};
  foreach my $gene (@$genes) {
    my $kept_transcript_count = 0;
    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts) {
      my $transcript_string = $self->transcript_string($transcript);
      if($transcript_strings->{$transcript_string}) {
#        $self->transfer_supporting_features($transcript,$transcript_strings->{$transcript_string})
      } else {
        $transcript_strings->{$transcript_string} = $transcript;
        $kept_transcript_count++;
      }
    } # end foreach my $transcript
    if($kept_transcript_count != 0) {
      push(@$output_genes,$gene)
    }
  } # end foreach my $gene

  return($output_genes);
}


sub transfer_supporting_features {
  my ($self,$transcript1,$transcript2) = @_;

  # This will transfer supporting evidence from transcript1 to transcript2
  my $exons1 = $transcript1->get_all_Exons();
  my $exons2 = $transcript2->get_all_Exons();

  unless(scalar(@$exons1) == scalar(@$exons2)) {
    $self->throw("Exon count differs between transcripts. Something is wrong as they should be identical");
  }

  for(my $i=0; $i<scalar(@$exons1); $i++) {
    my $exon1 = $$exons1[$i];
    my $exon2 = $$exons2[$i];
    unless($exon1->start == $exon2->start && $exon1->end == $exon2->end && $exon1->strand == $exon2->strand) {
      $self->throw("The exon coordinates do not match, something has gone wrong");
    }
    transfer_supporting_evidence($exon1,$exon2);
  }

  my $intron_supporting_features1 = $transcript1->get_all_IntronSupportingEvidence();
  foreach my $intron_support (@{$intron_supporting_features1}) {
    $transcript2->add_IntronSupportingEvidence($intron_support);
  }
}


sub transcript_string {
  my ($self,$transcript,$ignore_cds) = @_;

  my $transcript_string = $transcript->start.":".$transcript->end.":".$transcript->strand.":";
  unless($ignore_cds) {
    $transcript_string .= $transcript->cdna_coding_end.":".$transcript->cdna_coding_end.":";
  }
  my $exons = $transcript->get_all_Exons();
  my $exon_string = "";
  foreach my $exon (@$exons) {
    $exon_string .= ":".$exon->seq_region_start.":".$exon->seq_region_end.":";
  }
  $transcript_string .= $exon_string;
  return($transcript_string);
}


sub remove_redundant_donors {
  my ($self,$genes) = @_;

  # Note the below code is not ideal as it will not remove redundant genes if they happen to appear later in the list than a
  # identical gene with a worse priority. Most implementations I can think of to solve this are less than ideal so I am leaving
  # for the moment
  my $output_genes = [];
  my $transcript_strings = {};
  foreach my $gene (@$genes) {
    my $kept_transcript_count = 0;
    my $transcripts = $gene->get_all_Transcripts();
    foreach my $transcript (@$transcripts) {
      my $transcript_string = $self->transcript_string($transcript,1);
      if($transcript_strings->{$transcript_string}) {
        my $biotype1 = $transcript_strings->{$transcript_string}->biotype;
        my $biotype2 = $transcript->biotype;
        my $priority1 = $self->biotype_priorities($biotype1);
        my $priority2 = $self->biotype_priorities($biotype2);
        if($priority1 <= $priority2) {
#          $self->transfer_supporting_features($transcript,$transcript_strings->{$transcript_string});
	} else {
          # In this case we have a better priority for the current transcript, so we should transfer supporting evidence from
          # the one in the hash to it, and then make it the one in the hash
 #         $self->transfer_supporting_features($transcript_strings->{$transcript_string},$transcript);
          $transcript_strings->{$transcript_string} = $transcript;
          $kept_transcript_count++;
        }
      } else {
        $transcript_strings->{$transcript_string} = $transcript;
        $kept_transcript_count++;
      }
    } # end foreach my $transcript
    if($kept_transcript_count != 0) {
      push(@$output_genes,$gene)
    }
  } # end foreach my $gene
  return($output_genes);
}

1;
