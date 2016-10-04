#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveUTRAddition;

use strict;
use warnings;
use feature 'say';
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene attach_Slice_to_Gene);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(calculate_exon_phases);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');
use Data::Dumper;


sub fetch_input {
  my $self = shift;
  my $test_case = 0;

  my $analysis = Bio::EnsEMBL::Analysis->new(
                                              -logic_name => $self->param('logic_name'),
                                              -module => $self->param('module'),
                                            );
  $self->analysis($analysis);

  my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
  $self->hrdb_set_con($dna_dba,'dna_db');


  my $utr_output_dba = $self->hrdb_get_dba($self->param('utr_output_db'));
  $utr_output_dba->dnadb($dna_dba);
  $self->hrdb_set_con($utr_output_dba,'utr_output_db');

  my $input_id = $self->param('iid');
  my $input_id_type = $self->param('iid_type');

  my $utr_biotype_priorities = $self->param('utr_biotype_priorities');
  $self->biotype_priorities($utr_biotype_priorities);

  if($test_case) {
    $self->donor_transcripts($self->donor_test_cases());
    $self->acceptor_transcripts($self->acceptor_test_cases());
  } else {
    # If the input is an array of gene ids then loop through them and fetch them from the input db
    if($input_id_type eq 'slice') {
       $self->load_cluster($input_id);
    } else {
      $self->throw("You have selected an input id type using iid_type that is not supported by the module. Offending type: ".$input_id_type);
    }
  }

  return 1;
}

sub run {
  my $self = shift;

  my $donor_transcripts = $self->donor_transcripts();
  my $acceptor_transcripts = $self->acceptor_transcripts();

  say "Attempting to add utr to non utr transcripts";

  my $utr_transcripts = $self->add_utr($donor_transcripts,$acceptor_transcripts);

  $self->output_transcripts($utr_transcripts);

  return 1;
}

sub write_output {
  my $self = shift;

  my $adaptor = $self->hrdb_get_con('utr_output_db')->get_GeneAdaptor;
  my $slice = $self->cluster_slice();

#  my @output = @{$self->output_genes};
  my @output = @{$self->output_transcripts};
  say "Writing genes to output db";

  foreach my $transcript (@output){
    say "Attempting to store gene on: ".$slice->name();
    my $output_gene = Bio::EnsEMBL::Gene->new();
    $output_gene->analysis($transcript->analysis);
    $output_gene->biotype($transcript->biotype);
    $output_gene->add_Transcript($transcript);
    attach_Slice_to_Gene($output_gene,$slice);
    empty_Gene($output_gene);
    $adaptor->store($output_gene);
  }
  say "...finished writing genes to output db";

  return 1;
}


sub load_cluster {
  my ($self,$cluster) = @_;

  my $gene = new Bio::EnsEMBL::Gene();
  my $dna_dba = $self->hrdb_get_con('dna_db');
  my $input_dbs = $self->param('input_gene_dbs');

  my $donor_transcripts = [];
  my $acceptor_transcripts = [];

  # Loop through the adaptors listed in the cluster and retrieve the transcripts from the corresponding db using the adaptor
  # Transcripts are pushed onto either the donor or acceptor array. The key thing to note is that the acceptor array is expecting
  # an adaptor tagged as 'no_utr_db'. Every transcript from this db is an acceptor. Everything else is considered a donor
  my $internal_transcript_id = 0;
  foreach my $adaptor_name (keys(%{$cluster})) {
    unless($input_dbs->{$adaptor_name}) {
      $self->throw("You are using a cluster type input id, but one of the the adaptor names in the input id did not match a corresponding db hash\n".
                   "All adaptor names must have a correspondingly named db hash passed in. Offending adaptor name:\n".$adaptor_name);
    }
    my $dba = $self->hrdb_get_dba($input_dbs->{$adaptor_name});
    unless($dba) {
      $self->throw("You are using a cluster type input id, but one of the the adaptor names in the input id did not match a corresponding db hash".
                   "All adaptor names must have a correspondingly named db hash passed in. Offending adaptor name:\n".$adaptor_name);
    }
    $dba->dnadb($dna_dba);
    $self->hrdb_set_con($dba,$adaptor_name);
    my $transcript_adaptor = $dba->get_TranscriptAdaptor();
    my $transcript_id_array = $cluster->{$adaptor_name};
    foreach my $transcript_id (@{$transcript_id_array}) {
      my $transcript = $transcript_adaptor->fetch_by_dbID($transcript_id);

      # This is not being used at the moment (it's is used in similar code to be able to track transcripts using a unique id, since dbID may not be unique)
      $transcript->{'internal_transcript_id'} = $internal_transcript_id;
      $internal_transcript_id++;

      if($adaptor_name eq 'no_utr_db') {
        push(@{$acceptor_transcripts},$transcript);
      } else {
        push(@{$donor_transcripts},$transcript);
      }
    }

  }


  my $donor_count = scalar(@{$donor_transcripts});
  my $acceptor_count = scalar(@{$acceptor_transcripts});

  say "utr transcript count: ".$donor_count;
  $self->donor_transcripts($donor_transcripts);

  say "no utr transcript count: ".$acceptor_count;
  $self->acceptor_transcripts($acceptor_transcripts);

  # If we have no acceptor transcripts in a region then just mark the job as finished
  # Don't do the same for no donor transcripts as we want all acceptors to be output
  unless($acceptor_count) {
    $self->input_job->autoflow(0);
    $self->complete_early('No acceptor transcripts to process');
  }

   # This bit is slightly awkward, basically we want a slice, all the transcripts should have the same slice
   # so take if from the first acceptor (since at this point we know there is at least one acceptor) and then
   # get the toplevel equivalent. This doesn't really need to be done, but seems more correct
   # The important thing is that a slice is set at this point and can be retrieved later by cluster_slice
   my $slice_name = ${$acceptor_transcripts}[0]->slice()->seq_region_name;
   my $slice_adaptor = $dna_dba->get_SliceAdaptor();
   my $slice = $slice_adaptor->fetch_by_region('toplevel',$slice_name);
   say "Setting slice to: ".$slice->name();
   $self->cluster_slice($slice);
}


sub add_utr {
  my ($self,$donor_transcripts,$acceptor_transcripts) = @_;

  my $final_transcripts = [];
  foreach my $acceptor_transcript (@{$acceptor_transcripts}) {
    if($acceptor_transcript->biotype() =~ /^rnaseq_/ || scalar(@{$donor_transcripts}) == 0) {
      push(@{$final_transcripts},$acceptor_transcript);
      next;
    }

    my $cds_introns_a = $acceptor_transcript->get_all_CDS_Introns();
    my $cds_intron_string_a = $self->generate_intron_string($cds_introns_a);

    unless($cds_intron_string_a) {
      say "Acceptor has no CDS intorns and is therefore single exon. Not adding UTR";
      push(@{$final_transcripts},$acceptor_transcript);
      next;
    }

    $acceptor_transcript->{'5_prime_utr'} = 0;
    $acceptor_transcript->{'3_prime_utr'} = 0;
    say "Checking transcript ".$acceptor_transcript->dbID()." for potential UTR transcript match:";
    my $modified_acceptor_transcript_5prime;
    my $modified_acceptor_transcript_3prime;

    foreach my $donor_transcript (@{$donor_transcripts}) {

     my $priority = $self->biotype_priorities($donor_transcript->biotype);
     unless($priority) {
       $self->warning("Transcript biotype was not found in the biotype priorities hash or biotype was set to 0 priority. Skipping.".
                      "Biotype: ".$donor_transcript->biotype);
       next;
     }

     ########################
     # Add in some code for checking if the donor transcript has a CDS or not
     # If it does the behaviour should be changed from get_all_Introns to get_all_CDS_Introns
     ########################


     my $introns_b = $donor_transcript->get_all_Introns();
     say "\nCDS intron coords (A):";
     foreach my $intron (@{$cds_introns_a}) {
       print "(".$intron->start."..".$intron->end.")";
     }

     say "\nIntron coords (B):";
     foreach my $intron (@{$introns_b}) {
       print "(".$intron->start."..".$intron->end.")";
     }

     print "\n";

     if(scalar(@{$cds_introns_a}) > scalar(@{$introns_b})) {
       say "Acceptor has more introns than donor, so will not add UTR";
       next;
     }

     my $intron_string_b = $self->generate_intron_string($introns_b);

     # Unless we have a match of the cds intron coords of the target to the introns coords of the donor, return 0
     unless($intron_string_b =~ $cds_intron_string_a) {
       say "\n-----------------------------------------------------------------------";
       say "Acceptor CDS introns coords do not match a set in the donor transcript:";
       say $cds_intron_string_a." (acceptor intron coords)";
       say $intron_string_b." (donor intron coords)";
       say "-----------------------------------------------------------------------";
       next;
     }


     say "\n-----------------------------------------------------------------------------------------------------";
     say "Acceptor CDS introns coords match a set in the donor transcript, attempting to add UTR!!!!!!!!!!!";
     say "-----------------------------------------------------------------------------------------------------";

     if($acceptor_transcript->{'5_prime_utr'}) {
       say "5' UTR has been attached to the acceptor transcript already";
       # First check if the donor priority is worse (1=best), if it's worse then just skip
       if($priority > $modified_acceptor_transcript_5prime->{'priority'}) {
         say "No adding UTR as there is already 5' donor UTR from a biotype with a better priority";
         next;
       }
       my $new_transcript_5prime = $self->add_five_prime_utr($acceptor_transcript,$donor_transcript,$cds_introns_a,$introns_b,$cds_intron_string_a,$intron_string_b);
       if($new_transcript_5prime && ($new_transcript_5prime->length() > $modified_acceptor_transcript_5prime->length())) {
         say "A longer or higher UTR donor has been found, selecting as current 5' UTR";
         $modified_acceptor_transcript_5prime = $new_transcript_5prime;
         $modified_acceptor_transcript_5prime->{'priority'} = $priority;
       }
     } else {
       $modified_acceptor_transcript_5prime = $self->add_five_prime_utr($acceptor_transcript,$donor_transcript,$cds_introns_a,$introns_b,$cds_intron_string_a,$intron_string_b);
       if($modified_acceptor_transcript_5prime) {
         $modified_acceptor_transcript_5prime->{'priority'} = $priority;
       }
     }

     if($acceptor_transcript->{'3_prime_utr'}) {
       say "3' UTR has been attached to the acceptor transcript already";
       # First check if the donor priority is worse (1=best), if it's worse then just skip
       if($priority > $modified_acceptor_transcript_3prime->{'priority'}) {
         say "No adding UTR as there is already 3' donor UTR from a biotype with a better priority";
         next;
       }
       my $new_transcript_3prime = $self->add_three_prime_utr($acceptor_transcript,$donor_transcript,$cds_introns_a,$introns_b,$cds_intron_string_a,$intron_string_b);
       if($new_transcript_3prime && ($new_transcript_3prime->length > $modified_acceptor_transcript_3prime->length())) {
         say "A longer UTR donor has been found, selecting as current 3' UTR";
         $modified_acceptor_transcript_3prime = $new_transcript_3prime;
         $modified_acceptor_transcript_3prime->{'priority'} = $priority;
      }
     } else {
       $modified_acceptor_transcript_3prime = $self->add_three_prime_utr($acceptor_transcript,$donor_transcript,$cds_introns_a,$introns_b,$cds_intron_string_a,$intron_string_b);
       if($modified_acceptor_transcript_3prime) {
         $modified_acceptor_transcript_3prime->{'priority'} = $priority;
       }
     }

    } # End foreach my $donor_transcript

    # At this point we either have the final UTR on both ends, the final on one end only or no UTR. The tricky situation is when UTR has been added
    # to both ends. In this case we have to merge both transcripts into a single transcript
    if($modified_acceptor_transcript_5prime && $modified_acceptor_transcript_3prime) {
      say "Added both 5' and 3' UTR. Creating final joined transcript: ";
      my $joined_transcript = $self->join_transcripts($modified_acceptor_transcript_5prime,$modified_acceptor_transcript_3prime);
      $joined_transcript->biotype($acceptor_transcript->biotype);
      $self->add_transcript_supporting_features($joined_transcript,$acceptor_transcript);
      push(@{$final_transcripts},$joined_transcript);
    } elsif($modified_acceptor_transcript_5prime) {
      say "Added 5' UTR only";
      $modified_acceptor_transcript_5prime->biotype($acceptor_transcript->biotype);
      $self->add_transcript_supporting_features($modified_acceptor_transcript_5prime,$acceptor_transcript);
      push(@{$final_transcripts},$modified_acceptor_transcript_5prime);
    } elsif($modified_acceptor_transcript_3prime) {
      say "Added 3' UTR only";
      $modified_acceptor_transcript_3prime->biotype($acceptor_transcript->biotype);
      $self->add_transcript_supporting_features($modified_acceptor_transcript_3prime,$acceptor_transcript);
      push(@{$final_transcripts},$modified_acceptor_transcript_3prime);
    } else {
      say "No UTR added to transcript";
      push(@{$final_transcripts},$acceptor_transcript);
    }
  }

  return($final_transcripts);
}


sub add_five_prime_utr {
  my ($self,$transcript_a,$transcript_b,$cds_introns_a,$introns_b,$cds_intron_string_a,$intron_string_b) = @_;

  my $strand = $transcript_a->strand;
#  my $modified_transcript;

  say "\nAttempting to add 5' UTR";

  say "CDS INTRON STRING: ".$cds_intron_string_a;

  # At this point we have a match, now we need to locate the exon to merge
  my @cds_intron_coords_a = split(":",$cds_intron_string_a);
  my @intron_coords_b = split(":",$intron_string_b);

  my $five_prime_intron_a;
  my $exon_merge_index_b = -1;

  $five_prime_intron_a = $cds_intron_coords_a[0];
  for(my $i=0; $i<scalar(@intron_coords_b); $i++) {
    if($five_prime_intron_a eq $intron_coords_b[$i]) {
      say "Index of exon to the 5' side of terminal 5' intron found in donor at exon index: ".$i;
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
  my $merge_exon_candidate_a = ${$exons_a}[0];
  my $merge_exon_candidate_b = ${$exons_b}[$exon_merge_index_b];


  say "Merge candidate exon acceptor: ".$merge_exon_candidate_a->start."..".$merge_exon_candidate_a->end;
  say "Merge candidate exon donor: ".$merge_exon_candidate_b->start."..".$merge_exon_candidate_b->end;

  if($strand == 1) {
    if($merge_exon_candidate_b->start >= $merge_exon_candidate_a->start) {
      say "Merge candidate exon from donor is first exon and has a start that is >= acceptor first exon start, therefore not adding UTR";
      return(0);
    } elsif($merge_exon_candidate_b->end != $merge_exon_candidate_a->end) {
      $self->warning("The internal boundry coord of the donor did not match the acceptor, something is wrong.\n".
                     "Donor boundry: ".$merge_exon_candidate_b->end."\nAcceptor boundry: ".$merge_exon_candidate_a->end);
      return(0);
    }
  } else {
    if($merge_exon_candidate_b->end <= $merge_exon_candidate_a->end) {
      say "Merge candidate exon from donor is first exon and has a start that is <= acceptor first exon end (- strand) , therefore not adding UTR";
      return(0);
    } elsif($merge_exon_candidate_b->start != $merge_exon_candidate_a->start) {
      $self->warning("The internal boundry coord of the donor did not match the acceptor, something is wrong.\n".
                     "Donor boundry: ".$merge_exon_candidate_b->start."\nAcceptor boundry: ".$merge_exon_candidate_a->start);
      return(0);
    }
  }

  my $final_exons = [];
  my $final_translation;
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

  # If the donor exon is shorter (or the same length), then just make the acceptor boundry exon the start exon and add
  # all the other ones
  if(($strand == 1 && ($merge_exon_candidate_b->start >= $merge_exon_candidate_a->start)) ||
     ($strand == -1 && ($merge_exon_candidate_b->end <= $merge_exon_candidate_a->end))) {
    say "Donor boundry exon is shorter than candidate, therefore no merge of boundry exon data will occur";

    my $start_exon = new Bio::EnsEMBL::Exon(
                                           -START     => $merge_exon_candidate_a->start,
                                           -END       => $merge_exon_candidate_a->end,
                                           -STRAND    => $merge_exon_candidate_a->strand,
                                           -SLICE     => $merge_exon_candidate_a->slice,
                                           -ANALYSIS  => $self->analysis,
                                           -PHASE     => $merge_exon_candidate_a->phase,
                                           -END_PHASE => $merge_exon_candidate_a->end_phase);
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
      $supporting_features_a = $merge_exon_candidate_a->get_all_supporting_features();
      $out_exon->add_supporting_features(@{$supporting_features_a});
      push(@{$final_exons},$out_exon);
    }

    my $translation = Bio::EnsEMBL::Translation->new();
    $translation->start_Exon($start_exon);
    $translation->start($transcript_a->translation()->start());
    $translation->end_Exon($$final_exons[$#$final_exons]);
    $translation->end($transcript_a->translation()->end());
    $final_translation = $translation;
  } else {
    say "Donor boundry exon is longer than acceptor, therefore merge of exon data will occur";
    my $merge_exon = new Bio::EnsEMBL::Exon(
                                             -START     => $merge_exon_candidate_b->start,
                                             -END       => $merge_exon_candidate_b->end,
                                             -STRAND    => $merge_exon_candidate_b->strand,
                                             -SLICE     => $merge_exon_candidate_b->slice,
                                             -ANALYSIS  => $self->analysis,
                                             -PHASE     => -1,
                                             -END_PHASE => $merge_exon_candidate_a->phase);

    my $supporting_features_a = $merge_exon_candidate_a->get_all_supporting_features();
    $merge_exon->add_supporting_features(@{$supporting_features_a});

    push(@{$final_exons},$merge_exon);


    my $new_translation_start;
    if($strand == 1) {
      $new_translation_start = $transcript_a->translation->start + ($merge_exon_candidate_a->start - $merge_exon_candidate_b->start);
    } else {
      $new_translation_start = $transcript_a->translation->start + ($merge_exon_candidate_b->end - $merge_exon_candidate_a->end);
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

 #   # Add all the supporting features from the donor transcript
 #   for(my $i=0; $i<scalar(@{$exons_a}); $i++) {
 #     my $exon_b = ${$exons_b}[$i];
 #     my $supporting_features_b = $merge_exon_candidate_b->get_all_supporting_features();
 #     $$final_exons[$i]->add_supporting_features(@{$supporting_features_b});
 #   }

    my $translation = Bio::EnsEMBL::Translation->new();
    $translation->start_Exon($merge_exon);
    $translation->start($new_translation_start);
    $translation->end_Exon($$final_exons[$#$final_exons]);
    $translation->end($transcript_a->translation()->end());
    $final_translation = $translation;

  }

  say "Old translation start: ".$transcript_a->translation->start;
  say "New translation start: ".$final_translation->start;

  say "\nOriginal exon coords (A):";
  foreach my $exon (@{$exons_a}) {
    print "(".$exon->start."..".$exon->end.")";
  }

  say "\nOriginal exon coords (B):";
  foreach my $exon (@{$exons_b}) {
    print "(".$exon->start."..".$exon->end.")";
  }

  say "\nModified exon coords:";
  foreach my $exon (@{$final_exons}) {
    print "(".$exon->start."..".$exon->end.")";
  }

  my $modified_transcript = Bio::EnsEMBL::Transcript->new(-EXONS => $final_exons);

  # This is a basic sanity check on the UTR itself. First we want to check if the transcript is now abnormally longer (> 100KB longer)
  # If it is then calculate the average 5' intron length of the UTR. If this average is 5' intron length is > 35K then throw it out
  my $utr_intron_count = scalar(@{$final_exons}) - scalar(@{$exons_a});
  my $big_utr_length = 50000;
  my $max_no_intron_extension = 5000;
  my $max_average_5_prime_intron_length = 35000;
  
  my $modified_transcript_length = $modified_transcript->seq_region_end() - $modified_transcript->seq_region_start();
  my $transcript_a_length = $transcript_a->seq_region_end() - $transcript_a->seq_region_start();
  my $added_length = $modified_transcript_length - $transcript_a_length;  

  if($utr_intron_count == 0 && $added_length > $max_no_intron_extension) {
    say "\nNot adding UTR as no introns were present but transcript length was extended past max allowed value.";
    say "Max allowed genomic extension (for UTR with no introns): ".$max_no_intron_extension;
    say "Observed extension length: ".$added_length;
    return(0);
  } elsif ($added_length > $big_utr_length) {
      return 0; 
  }
    
#  } elsif($utr_intron_count >= 0 && $added_length >= $big_utr_length) {
#    say "Modified transcript is significantly longer due to 5' UTR addition (length added: ".$added_length."), calculating average 5' intron length.";
#    my $total_length = 0;
#    for(my $i=0; $i<$utr_intron_count; $i++) {
#      my $exon_1 = ${$exons_b}[$i];
#      my $exon_2 = ${$exons_b}[$i+1];
#      if($strand == 1) {
#        $total_length += $exon_2->start - $exon_1->end + 1;
#      } else {
#        $total_length += $exon_1->start - $exon_2->end + 1;
#      }
#    }
#    my $average_5_prime_intron_length = $total_length / $utr_intron_count;
#    if($average_5_prime_intron_length > $max_average_5_prime_intron_length) {
#      say "\nAverage intron length exceeds the max allowed value for 5' UTR, not adding UTR";
#      say "Allowed max average 5' UTR intron size: ".$max_average_5_prime_intron_length;
#      say "Observed average 5' UTR intron size: ".$average_5_prime_intron_length;
#      return 0;
#    }
# }

  $modified_transcript->analysis($transcript_a->analysis);
  $modified_transcript->biotype($transcript_a->biotype);
  $modified_transcript->slice($transcript_a->slice());
  $modified_transcript->translation($final_translation);


  say "\n";
  say "Acceptor original sequence:\n".$transcript_a->seq->seq;
  say "Acceptor current sequence:\n".$modified_transcript->seq->seq;
  say "Acceptor original translation:\n".$transcript_a->translation->seq;
  say "Acceptor current translation (from translateable seq):\n".$modified_transcript->translate->seq;
  say "Acceptor current translation (from translation object string):\n".$modified_transcript->translation->seq;

  unless($transcript_a->translation->seq eq $modified_transcript->translate->seq && $modified_transcript->translate->seq eq $modified_transcript->translation->seq) {
    $self->throw("There is an issue with the translation after UTR was added. Check above for the sequences, all three should match");
  }

  $transcript_a->{'5_prime_utr'} = 1;

  return($modified_transcript);
}


sub add_three_prime_utr {
  my ($self,$transcript_a,$transcript_b,$cds_introns_a,$introns_b,$cds_intron_string_a,$intron_string_b) = @_;

  my $strand = $transcript_a->strand;
#  my $modified_transcript;

  say "\nAttempting to add 3' UTR";

  say "CDS INTRON STRING: ".$cds_intron_string_a;

  # At this point we have a match, now we need to locate the exon to merge
  my @cds_intron_coords_a = split(":",$cds_intron_string_a);
  my @intron_coords_b = split(":",$intron_string_b);

  my $three_prime_intron_a = $cds_intron_coords_a[$#cds_intron_coords_a];
  my $exon_merge_index_b = -1;
  for(my $i=0; $i<scalar(@intron_coords_b); $i++) {
    if($three_prime_intron_a eq $intron_coords_b[$i]) {
      say "Index of exon to the 3' side of terminal 3' intron found in donor at exon index: ".$i;
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


  say "Merge candidate exon acceptor: ".$merge_exon_candidate_a->start."..".$merge_exon_candidate_a->end;
  say "Merge candidate exon donor: ".$merge_exon_candidate_b->start."..".$merge_exon_candidate_b->end;

  if($strand == 1) {
    if($merge_exon_candidate_b->end <= $merge_exon_candidate_a->end) {
      say "Merge candidate exon from donor is last exon and has an end that is <= acceptor last exon end, therefore not adding UTR";
      return(0);
    } elsif($merge_exon_candidate_b->start != $merge_exon_candidate_a->start) {
      $self->warning("The internal boundry coord of the donor did not match the acceptor, something is wrong.\n".
                          "Donor boundry: ".$merge_exon_candidate_b->start."\nAcceptor boundry: ".$merge_exon_candidate_a->start);
      return(0);
    }
  } else {
    if($merge_exon_candidate_b->start >= $merge_exon_candidate_a->start) {
      say "Merge candidate exon from donor is last exon and has a start that is >= acceptor first exon end (- strand) , therefore not addign UTR";
      return(0);
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

    my $supporting_features_a = $merge_exon_candidate_a->get_all_supporting_features();
    $out_exon->add_supporting_features(@{$supporting_features_a});
    push(@{$final_exons},$out_exon);
  }

  my $final_translation;

  # Now look at the merge candidates. There are a few things to check. If the start coord of the donor exon is greater than the start coord of the
  # acceptor then just push the acceptor exon on the final exons array
  if(($strand == 1 && ($merge_exon_candidate_b->end <= $merge_exon_candidate_a->end)) ||
     ($strand == -1 && ($merge_exon_candidate_b->start >= $merge_exon_candidate_a->start))) {
    say "Donor boundry exon is shorter than candidate, therefore no merge of boundry exon data will occur";
    my $translation = Bio::EnsEMBL::Translation->new();
    $translation->start_Exon($$final_exons[0]);
    $translation->start($transcript_a->translation()->start());
    $translation->end_Exon($$final_exons[$#$final_exons]);
    $translation->end($transcript_a->translation()->end());
    $final_translation = $translation;
  }  else {
    say "Donor boundry exon is longer than acceptor, therefore merge of boundry exon data will occur";
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

    my $translation = Bio::EnsEMBL::Translation->new();
    $translation->start_Exon($$final_exons[0]);
    $translation->start($transcript_a->translation()->start());
    $translation->end_Exon($merge_exon);
    $translation->end($transcript_a->translation()->end());
    $final_translation = $translation;

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

  say "Old translation start: ".$transcript_a->translation->start;
  say "New translation start: ".$final_translation->start;


  say "\nOriginal exon coords:";
  foreach my $exon (@{$exons_a}) {
    print "(".$exon->start."..".$exon->end.")";
  }

  say "\nModified exon coords:";
  foreach my $exon (@{$final_exons}) {
    print "(".$exon->start."..".$exon->end.")";
  }

  # Add all the supporting features from the donor transcript
#  for(my $i=0; $i<scalar(@{$exons_a}); $i++) {
#    my $exon_b = ${$exons_b}[$i];
#    my $supporting_features_b = $merge_exon_candidate_b->get_all_supporting_features();
#    $$final_exons[$i]->add_supporting_features(@{$supporting_features_b});
#  }

  my $modified_transcript = Bio::EnsEMBL::Transcript->new(-EXONS => $final_exons);

  # This is a basic sanity check on the UTR itself. First we want to check if the transcript is now abnormally longer (> 100KB longer)
  # If it is then calculate the average 3' intron length of the UTR. If this average is 3' intron length is > 25K then throw it out
  my $utr_intron_count = scalar(@{$final_exons}) - scalar(@{$exons_a});
  my $big_utr_length = 50000;
  my $max_no_intron_extension = 5000;
  my $max_average_3_prime_intron_length = 35000;

  my $modified_transcript_length = $modified_transcript->seq_region_end() - $modified_transcript->seq_region_start();
  my $transcript_a_length = $transcript_a->seq_region_end() - $transcript_a->seq_region_start();
  my $added_length = $modified_transcript_length - $transcript_a_length;

  if($utr_intron_count == 0 && $added_length > $max_no_intron_extension) {
    say "\nNot adding UTR as no introns were present but transcript length was extended past max allowed value.";
    say "Max allowed genomic extension (for UTR with no introns): ".$max_no_intron_extension;
    say "Observed extension length: ".$added_length;
    return(0);
    
  } elsif ($added_length > $big_utr_length) {
      return 0; 
  }
    
#  } elsif($utr_intron_count >= 0 && $added_length >= $big_utr_length) {
#    say "Modified transcript is significantly longer due to 3' UTR addition (length added: ".$added_length."), calculating average 3' intron length.";
#    my $total_length = 0;
#    for(my $i=$exon_merge_index_b; $i<$utr_intron_count; $i++) {
#      my $exon_1 = ${$exons_b}[$i];
#      my $exon_2 = ${$exons_b}[$i+1];
#      if($strand == 1) {
#        $total_length += $exon_2->start - $exon_1->end + 1;
#      } else {
#        $total_length += $exon_1->start - $exon_2->end + 1;
#      }
#    }
#    my $average_3_prime_intron_length = $total_length / $utr_intron_count;
#    if($average_3_prime_intron_length > $max_average_3_prime_intron_length) {
#      say "\nAverage intron length exceeds the max allowed value for 3' UTR, not adding UTR";
#      say "Allowed max average 3' UTR intron size: ".$max_average_3_prime_intron_length;
#      say "Observed average 3' UTR intron size: ".$average_3_prime_intron_length;
#      return 0;
#    }
#  }

  $modified_transcript->analysis($transcript_a->analysis);
  $modified_transcript->biotype($transcript_a->biotype);
  $modified_transcript->slice($transcript_a->slice());
  $modified_transcript->translation($final_translation);


  my $modified_translation = $modified_transcript->translation();
  say "\n";
  say "Acceptor original sequence:\n".$transcript_a->seq->seq;
  say "Acceptor original translateable seq:\n".$transcript_a->translateable_seq();
  say "Acceptor current sequence:\n".$modified_transcript->seq->seq;
  say "Acceptor current translateable seq:\n".$modified_transcript->translateable_seq();
  say "Acceptor original translation:\n".$transcript_a->translation->seq;
  say "Acceptor current translation (from translateable seq):\n".$modified_transcript->translate->seq;
  say "Acceptor current translation (from translation object string):\n".$modified_transcript->translation->seq;

  unless($transcript_a->translation->seq eq $modified_transcript->translate->seq && $modified_transcript->translate->seq eq $modified_transcript->translation->seq) {
    $self->throw("There is an issue with the translation after UTR was added. Check above for the sequences, all three should match");
  }

  $transcript_a->{'3_prime_utr'} = 1;

  return($modified_transcript);
}


sub join_transcripts {
  my ($self,$transcript_a,$transcript_b) = @_;

#  my $joined_transcript;
  my $joined_exon_set = [];
  unless($transcript_a->translation->seq eq $transcript_b->translation->seq) {
    $self->throw("When attempting to join the modified 5' and 3' transcripts, there was a difference in the translation for each. The translation should be ".
                 "identical at this point between the two.\n5' translation:\n".$transcript_a->translation->seq."\n3' translation:\n".$transcript_b->translation->seq);
  }

  say "Translations from both partial transcripts are identical. Attemting to merge exon sets";
  my $exons_a = $transcript_a->get_all_Exons;
  my $exons_b = $transcript_b->get_all_Exons;

  my $cds_5_prime_index = scalar(@{$exons_a}) - scalar(@{$transcript_a->get_all_CDS});
  my $cds_3_prime_index = scalar(@{$transcript_b->get_all_CDS}) - 1;

  say "CDS 5' exon index: ".$cds_5_prime_index;
  say "CDS 3' exon index: ".$cds_3_prime_index;

  foreach my $exon (@{$exons_a}) {
    my $out_exon = new Bio::EnsEMBL::Exon(
                                           -START  => $exon->start,
                                           -END       => $exon->end,
                                           -STRAND    => $exon->strand,
                                           -SLICE     => $exon->slice,
                                           -ANALYSIS  => $self->analysis,
                                           -PHASE     => $exon->phase,
                                           -END_PHASE => $exon->end_phase);
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
    push(@{$joined_exon_set},$out_exon);
  }

  my $joined_offset = scalar(@{$joined_exon_set}) - scalar(@{$exons_b});
  $cds_3_prime_index += $joined_offset;

  # The translation is the same, but still need to modify the translation so that it has the correct start and end exon
  my $translation = Bio::EnsEMBL::Translation->new();
  $translation->start_Exon($$joined_exon_set[$cds_5_prime_index]);
  $translation->start($transcript_a->translation()->start());
  $translation->end_Exon($$joined_exon_set[$cds_3_prime_index]);
  $translation->end($transcript_a->translation()->end());

  my $joined_transcript = Bio::EnsEMBL::Transcript->new(-EXONS => $joined_exon_set);
  $joined_transcript->analysis($transcript_a->analysis);
  $joined_transcript->biotype($transcript_a->biotype);
  $joined_transcript->slice($transcript_a->slice());
  $joined_transcript->translation($translation);

  unless($transcript_a->translation->seq eq $joined_transcript->translation->seq) {
    $self->throw("When attempting to join the modified 5' and 3' transcripts, there was a difference in the joined translation. The translation should be identical at ".
                 "this point to the original.\nOriginal translation:\n".$transcript_a->translation->seq."\nJoined translation:\n".$joined_transcript->translation->seq);
  }

  say "Joined translation: ".$joined_transcript->translation->seq();

  say "Joined exon coords:";
  foreach my $exon (@{$joined_exon_set}) {
    print "(".$exon->start."..".$exon->end.")";
  }
  print "\n";

  return($joined_transcript);
}


sub generate_intron_string {
  my ($self,$intron_array) = @_;

  say "FM2 GEN INT: ";
  my $intron_string = "";
  foreach my $intron (@{$intron_array}) {
    my $start = $intron->start();
    my $end = $intron->end();
    $intron_string .= $start."..".$end.":";
    print "(".$start."..".$end.")";
  }

  print "\n";

  return($intron_string);
}


sub check_cds_introns {
  my ($self,$cds_a,$cds_b) = @_;

  my $exon_count_a = scalar(@{$cds_a});
  my $exon_count_b = scalar(@{$cds_b});

  unless($exon_count_a && $exon_count_b) {
    return 0;
  }

  unless($exon_count_a == $exon_count_b) {
    return 0;
  }

  if($exon_count_a == 1) {
    my $exon_a = shift(@{$cds_a});
    my $exon_b = shift(@{$cds_b});
    if($self->features_overlap($exon_a,$exon_b)) {
      return 1;
    } else {
      return 0;
    }
  } else {
    my $intron_coords_a = $self->intron_coords($cds_a);
    my $intron_coords_b = $self->intron_coords($cds_b);

    if($intron_coords_a eq $intron_coords_b) {
      return 1;
    } else {
      return 0;
    }
  }

}


sub add_transcript_supporting_features {
  my ($self,$transcript_a,$transcript_b) = @_;

  my @supporting_features = ();
  # transcript a is the transcript to add to, transcript b is the one with the sfs
  foreach my $sf (@{$transcript_b->get_all_supporting_features()}) {
    if(features_overlap($sf,$transcript_a)) {
        push(@supporting_features,$sf);
    }
  }

  $transcript_a->add_supporting_features(@supporting_features);

  my @intron_support = @{$transcript_b->get_all_IntronSupportingEvidence()};

  foreach my $intron_support (@intron_support) {
    $transcript_a->add_IntronSupportingEvidence($intron_support);
  }

}

sub features_overlap {
  my ( $featureA, $featureB ) = @_;

  if (
     ( $featureA->seq_region_start() <= $featureB->seq_region_end() ) &&
     ( $featureA->seq_region_end() >= $featureB->seq_region_start() ) )
    {
    return 1;
  }

  return 0;
}

sub features_overlap_new {
  my ($self,$feature_a,$feature_b) = @_;

  if(($feature_a->seq_region_start() <= $feature_b->seq_region_end() ) && ($feature_a->seq_region_end() >= $feature_b->seq_region_start())) {
    return 1;
  }

  return 0;
}

sub intron_coords {
  my ($self,$exons) = @_;

  my $intron_coords = "";
  for(my $i=0; $i<(scalar@{$exons})-1; $i++) {
    my $exon_a = $$exons[$i];
    my $exon_b = $$exons[$i+1];
    $intron_coords .= $exon_a->end()."..".$exon_b->start().":";
  }

  unless($intron_coords) {
    $self->throw("Issue with finding the intron coordinates!");
  }

  return($intron_coords);
}



sub recluster_genes {
  my ($self,$processed_genes) = @_;

  my $single_transcript_genes = [];
  my $output_genes = [];
  foreach my $processed_gene (@{$processed_genes}) {
    my $transcripts = $processed_gene->get_all_Transcripts();

    # If the gene has a single transcript then there is no need to recluster
    if(scalar(@{$transcripts}) == 1) {
      push(@{$single_transcript_genes},$processed_gene);
      next;
    }

    # If there is more than one transcript, then make a new gene for each
    elsif(scalar(@{$transcripts}) > 1) {
      foreach my $transcript (@{$transcripts}) {
        my $single_transcript_gene = Bio::EnsEMBL::Gene->new();
        $single_transcript_gene->add_Transcript($transcript);
        push(@{$single_transcript_genes},$single_transcript_gene);
      }
    }

    else {
      $self->throw("Got a gene without any transcripts");
    }
  }


  my $biotypes_hash = $self->get_all_biotypes($single_transcript_genes);
  my $biotypes_array = [keys(%$biotypes_hash)];
  my $types_hash;
  $types_hash->{genes} = $biotypes_array;

  say "Reclustering gene models based on final transcript set";
  my ($clusters, $unclustered) = cluster_Genes($single_transcript_genes,$types_hash);
  say "...finished reclustering genes models";
  say "Found clustered sets: ".scalar(@{$clusters});
  say "Found unclustered sets: ".scalar(@{$unclustered});

  say "Processing gene clusters into single gene models...";
  $output_genes = $self->process_clusters($clusters,$unclustered);
  say "...finished processing gene clusters into single gene models";
  say "Made output genes: ".scalar(@{$output_genes});

  $self->output_genes($output_genes);

}


sub donor_transcripts {
  my ($self,$val) = @_;

  if($val) {
    $self->param('_donor_transcripts',$val);
  }

  return($self->param('_donor_transcripts'));
}

sub acceptor_transcripts {
  my ($self,$val) = @_;

  if($val) {
    $self->param('_acceptor_transcripts',$val);
  }

  return($self->param('_acceptor_transcripts'));
}


sub output_genes {
  my ($self,$val) = @_;

  if($val) {
    $self->param('_output_genes',$val);
  }

  return($self->param('_output_genes'));
}

sub output_transcripts {
  my ($self,$val) = @_;

  if($val) {
    $self->param('_output_transcripts',$val);
  }

  return($self->param('_output_transcripts'));
}


sub biotype_priorities {
  my ($self,$biotype_priorities_hash) = @_;
  if($biotype_priorities_hash) {
    $self->param('_biotype_priorities',$biotype_priorities_hash);
  }
  return($self->param('_biotype_priorities'));
}

sub get_all_biotypes {
  my ($self,$master_genes_array) = @_;

  my $master_biotypes_hash = {};

  foreach my $gene (@{$master_genes_array}) {
    unless($master_biotypes_hash->{$gene->biotype}) {
      $master_biotypes_hash->{$gene->biotype} = 1;
    }
  }
  return($master_biotypes_hash);
}


sub process_clusters {

  my ($self,$clustered,$unclustered) = @_;

  my $output_genes = [];
  foreach my $single_cluster (@{$unclustered}) {
    my $cluster_genes = $single_cluster->get_Genes();
    foreach my $single_gene (@{$cluster_genes}) {
      my $output_gene = Bio::EnsEMBL::Gene->new();
      $output_gene->slice($self->query());
      $output_gene->analysis($self->analysis);
      $output_gene->biotype($self->analysis->logic_name);
      my $transcripts = $single_gene->get_all_Transcripts();
      my $single_transcript = shift(@{$transcripts});
      $output_gene->add_Transcript($single_transcript);
      push(@{$output_genes},$output_gene);
    }
  }

  foreach my $single_cluster (@{$clustered}) {
    my $combined_gene = Bio::EnsEMBL::Gene->new();
    $combined_gene->slice($self->query());
    $combined_gene->analysis($self->analysis);
    $combined_gene->biotype($self->analysis->logic_name);

    my $cluster_genes = $single_cluster->get_Genes();

    foreach my $single_gene (@{$cluster_genes}) {
      my $transcripts = $single_gene->get_all_Transcripts();
      my $single_transcript = shift(@{$transcripts});
      $combined_gene->add_Transcript($single_transcript);
    }
    push(@{$output_genes},$combined_gene);
  }

  return($output_genes);

}

sub cluster_slice {
  my ($self,$slice) = @_;
  if($slice) {
    $self->param('_cluster_slice',$slice);
  }

  return($self->param('_cluster_slice'));
}

sub donor_test_cases {
  my ($self) = @_;

  my $dba = $self->hrdb_get_con('utr_output_db');
  my $slice_adaptor = $dba->get_SliceAdaptor();
  my $slice = $slice_adaptor->fetch_by_name('chromosome:Mmul_8.0.1:9:1:129882849:1');

  my $exon_1 = new Bio::EnsEMBL::Exon(
                                       -START     => 100,
                                       -END       => 199,
                                       -STRAND    => 1,
                                       -SLICE     => $slice,
                                       -ANALYSIS  => $self->analysis,
                                       -PHASE     => -1);

  my $exon_2 = new Bio::EnsEMBL::Exon(
                                       -START     => 300,
                                       -END       => 399,
                                       -STRAND    => 1,
                                       -SLICE     => $slice,
                                       -ANALYSIS  => $self->analysis,
                                       -PHASE     => -1);

  my $exon_3 = new Bio::EnsEMBL::Exon(
                                       -START     => 500,
                                       -END       => 599,
                                       -STRAND    => 1,
                                       -SLICE     => $slice,
                                       -ANALYSIS  => $self->analysis,
                                       -PHASE     => -1);

  my $exon_4 = new Bio::EnsEMBL::Exon(
                                       -START     => 700,
                                       -END       => 799,
                                       -STRAND    => 1,
                                       -SLICE     => $slice,
                                       -ANALYSIS  => $self->analysis,
                                       -PHASE     => -1);

  my $exon_5 = new Bio::EnsEMBL::Exon(
                                       -START     => 900,
                                       -END       => 999,
                                       -STRAND    => 1,
                                       -SLICE     => $slice,
                                       -ANALYSIS  => $self->analysis,
                                       -PHASE     => -1);

  my  @exons_set_1 = ($exon_1,$exon_2,$exon_3,$exon_4,$exon_5);

  my $transcript_1 = new Bio::EnsEMBL::Transcript(
                                                   -EXONS => \@exons_set_1,
                                                   -STRAND    => 1,
                                                   -SLICE     => $slice,
                                                   -ANALYSIS  => $self->analysis);


  my $exon_6 = new Bio::EnsEMBL::Exon(
                                       -START     => 900,
                                       -END       => 999,
                                       -STRAND    => -1,
                                       -SLICE     => $slice,
                                       -ANALYSIS  => $self->analysis,
                                       -PHASE     => -1);

  my $exon_7 = new Bio::EnsEMBL::Exon(
                                       -START     => 1100,
                                       -END       => 1199,
                                       -STRAND    => -1,
                                       -SLICE     => $slice,
                                       -ANALYSIS  => $self->analysis,
                                       -PHASE     => -1);

  my $exon_8 = new Bio::EnsEMBL::Exon(
                                       -START     => 1300,
                                       -END       => 1399,
                                       -STRAND    => -1,
                                       -SLICE     => $slice,
                                       -ANALYSIS  => $self->analysis,
                                       -PHASE     => -1);

  my $exon_9 = new Bio::EnsEMBL::Exon(
                                       -START     => 1500,
                                       -END       => 1599,
                                       -STRAND    => -1,
                                       -SLICE     => $slice,
                                       -ANALYSIS  => $self->analysis,
                                       -PHASE     => -1);

  my $exon_10 = new Bio::EnsEMBL::Exon(
                                       -START     => 1700,
                                       -END       => 1799,
                                       -STRAND    => -1,
                                       -SLICE     => $slice,
                                       -ANALYSIS  => $self->analysis,
                                       -PHASE     => -1);

  my $exon_11 = new Bio::EnsEMBL::Exon(
                                       -START     => 1900,
                                       -END       => 1999,
                                       -STRAND    => -1,
                                       -SLICE     => $slice,
                                       -ANALYSIS  => $self->analysis,
                                       -PHASE     => -1);

  my $exon_12 = new Bio::EnsEMBL::Exon(
                                       -START     => 2100,
                                       -END       => 2199,
                                       -STRAND    => -1,
                                       -SLICE     => $slice,
                                       -ANALYSIS  => $self->analysis,
                                       -PHASE     => -1);

  my $exon_13 = new Bio::EnsEMBL::Exon(
                                       -START     => 2300,
                                       -END       => 2399,
                                       -STRAND    => -1,
                                       -SLICE     => $slice,
                                       -ANALYSIS  => $self->analysis,
                                       -PHASE     => -1);

  my  @exons_set_2 = ($exon_12,$exon_11,$exon_10,$exon_9,$exon_8,$exon_7,$exon_6);

  my $transcript_2 = new Bio::EnsEMBL::Transcript(
                                                   -EXONS => \@exons_set_2,
                                                   -STRAND    => -1,
                                                   -SLICE     => $slice,
                                                   -ANALYSIS  => $self->analysis);

  my  @exons_set_3 = ($exon_13,$exon_12,$exon_11,$exon_10,$exon_9,$exon_8,$exon_7,$exon_6);

  my $transcript_3 = new Bio::EnsEMBL::Transcript(
                                                   -EXONS => \@exons_set_3,
                                                   -STRAND    => -1,
                                                   -SLICE     => $slice,
                                                   -ANALYSIS  => $self->analysis);

  say "Created the following test donor transcripts: ";
  say "DONOR T1: (".$transcript_1->start.":".$transcript_1->end.":".$transcript_1->strand.")";
  my $exons = $transcript_1->get_all_Exons();
  foreach my $exon (@{$exons}) {
    print "(".$exon->start."..".$exon->end.")";
  }
  print "\n";

  say "DONOR T2: (".$transcript_2->start.":".$transcript_2->end.":".$transcript_2->strand.")";
  $exons = $transcript_2->get_all_Exons();
  foreach my $exon (@{$exons}) {
    print "(".$exon->start."..".$exon->end.")";
  }
  print "\n";

  say "DONOR T3: (".$transcript_3->start.":".$transcript_3->end.":".$transcript_3->strand.")";
  $exons = $transcript_3->get_all_Exons();
  foreach my $exon (@{$exons}) {
    print "(".$exon->start."..".$exon->end.")";
  }
  print "\n";

  return([$transcript_1,$transcript_2,$transcript_3]);
}

sub acceptor_test_cases {
  my ($self) = @_;

  my $dba = $self->hrdb_get_con('utr_output_db');
  my $slice_adaptor = $dba->get_SliceAdaptor();
  my $slice = $slice_adaptor->fetch_by_name('chromosome:Mmul_8.0.1:9:1:129882849:1');

  my $exon_1 = new Bio::EnsEMBL::Exon(
                                       -START     => 550,
                                       -END       => 599,
                                       -STRAND    => 1,
                                       -SLICE     => $slice,
                                       -ANALYSIS  => $self->analysis,
                                       -PHASE     => 0,
                                       -END_PHASE => 0);

  my $exon_2 = new Bio::EnsEMBL::Exon(
                                      -START     => 700,
                                      -END       => 759,
                                      -STRAND    => 1,
                                      -SLICE     => $slice,
                                      -ANALYSIS  => $self->analysis,
                                      -PHASE     => 0,
                                      -END_PHASE => 0);

  my $translation_1 = new Bio::EnsEMBL::Translation(
                                                     -START_EXON => $exon_1,
                                                     -END_EXON => $exon_2,
                                                     -SEQ_START => 1,
                                                     -SEQ_END => 49,
                                                   );
  my  @exons_set_1 = ($exon_1,$exon_2);

  my $transcript_1 = new Bio::EnsEMBL::Transcript( -DBID  => 1,
                                                   -EXONS => \@exons_set_1,
                                                   -STRAND    => 1,
                                                   -SLICE     => $slice,
                                                   -ANALYSIS  => $self->analysis);

  $transcript_1->translation($translation_1);

  my $exon_3 = new Bio::EnsEMBL::Exon(
                                       -START     => 1350,
                                       -END       => 1399,
                                       -STRAND    => -1,
                                       -SLICE     => $slice,
                                       -ANALYSIS  => $self->analysis,
                                       -PHASE     => 0,
                                       -END_PHASE => 0);

  my $exon_4 = new Bio::EnsEMBL::Exon(
                                      -START     => 1500,
                                      -END       => 1599,
                                      -STRAND    => -1,
                                      -SLICE     => $slice,
                                      -ANALYSIS  => $self->analysis,
                                      -PHASE     => 0,
                                      -END_PHASE => 0);

  my $exon_5 = new Bio::EnsEMBL::Exon(
                                      -START     => 1700,
                                      -END       => 1759,
                                      -STRAND    => -1,
                                      -SLICE     => $slice,
                                      -ANALYSIS  => $self->analysis,
                                      -PHASE     => 0,
                                      -END_PHASE => 0);

  my $translation_2 = new Bio::EnsEMBL::Translation(
                                                     -START_EXON => $exon_5,
                                                     -END_EXON => $exon_3,
                                                     -SEQ_START => 1,
                                                     -SEQ_END => 49,
                                                   );
  my  @exons_set_2 = ($exon_5,$exon_4,$exon_3);

  my $transcript_2 = new Bio::EnsEMBL::Transcript( -DBID  => 2,
                                                   -EXONS => \@exons_set_2,
                                                   -STRAND    => -1,
                                                   -SLICE     => $slice,
                                                   -ANALYSIS  => $self->analysis);

  $transcript_2->translation($translation_2);

  say "ACCEPTOR T1: (".$transcript_1->start.":".$transcript_1->end.":".$transcript_1->strand.")";
  my $exons = $transcript_1->get_all_Exons();
  foreach my $exon (@{$exons}) {
    print "(".$exon->start."..".$exon->end.")";
  }
  print "\n";
  say "ACCEPTOR TN1: ".$transcript_1->translation()->seq();


  say "ACCEPTOR T2: (".$transcript_2->start.":".$transcript_2->end.":".$transcript_2->strand.")";
  $exons = $transcript_2->get_all_Exons();
  foreach my $exon (@{$exons}) {
    print "(".$exon->start."..".$exon->end.")";
  }
  print "\n";
  say "ACCEPTOR TN2: ".$transcript_2->translation()->seq();

  return([$transcript_1,$transcript_2]);
}

1;
