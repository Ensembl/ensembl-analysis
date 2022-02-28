=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCleanGeneset

=head1 SYNOPSIS


=head1 DESCRIPTION

Clean te geneset to remove really bad transcripts which have small introns
or single exon models in an intron

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCleanGeneset;

use strict;
use warnings;
use feature 'say';

use File::Spec::Functions qw(catfile);
use File::Path qw(make_path);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub fetch_input {
  my $self = shift;

  if($self->param('skip_analysis')) {
    $self->complete_early('Skip analysis flag is enabled, so no cleaning will occur');
  }

  $self->create_analysis;

  my $dna_dba = $self->get_database_by_name('dna_db');

  my $input_dba = $self->get_database_by_name('input_db', $dna_dba);
  $self->hrdb_set_con($input_dba,'input_db');

  $self->load_genes();

  my $output_path = $self->param_required('output_path');
  unless(-e $output_path) {
    make_path($output_path);
  }

  my $blessed_biotypes = $self->param('blessed_biotypes');
  $self->blessed_biotypes($blessed_biotypes);

  return 1;
}

sub run {
  my ($self) = @_;


  my $genes = $self->genes_to_qc();
  my $transcript_ids_to_remove = $self->clean_genes($genes);
  $self->transcript_ids_to_remove($transcript_ids_to_remove);
  return 1;
}

sub write_output {
  my $self = shift;

  my $transcript_ids_to_remove = $self->transcript_ids_to_remove();
  my $output_file = catfile($self->param('output_path'), 'transcript_ids_to_remove.txt');
  open(OUT, ">$output_file") || $self->throw("Could not open $output_file for writing");
  foreach my $id (@{$transcript_ids_to_remove}) {
    say OUT $id;
  }
  close(OUT) || $self->throw("Could not open $output_file for writing");

  return 1;
}

sub load_genes {
  my ($self) = @_;

  my $dba = $self->hrdb_get_con('input_db');
  my $gene_adaptor = $dba->get_GeneAdaptor();
  my $genes = $gene_adaptor->fetch_all();#_by_Slice($slice);
  $self->genes_to_qc($genes);
}


sub genes_to_qc {
  my ($self,$genes) = @_;
  if($genes) {
    $self->param('_genes_to_qc',$genes);
  }

  return($self->param('_genes_to_qc'));
}


sub blessed_biotypes {
  my ($self,$blessed_biotypes) = @_;
  if($blessed_biotypes) {
    $self->param('_blessed_biotypes',$blessed_biotypes);
  }

  return($self->param('_blessed_biotypes'));
}


sub clean_genes {
  my ($self,$genes) = @_;
  my $transcript_ids_to_remove = [];
  my $tiny_gene_size = 500; # Note that this should be made a param. Usually there are very few such genes in a geneset
  say "Have ".scalar(@{$genes})." genes to process";

  my $gene_strings;
  foreach my $gene (@{$genes}) {
    my $gene_string = $gene->start.":".$gene->end.":".$gene->seq_region_name;
    my $transcripts = $gene->get_all_Transcripts();
    my ($normal_transcripts,$flagged_transcripts) = @{$self->flag_transcripts($transcripts,$gene_string)};
    say "Checking gene with dbID: ".$gene->dbID;
    say "Total transcript count: ".scalar(@{$transcripts});
    say "Normal transcript count: ".scalar(@{$normal_transcripts});
    say "Flagged transcript count: ".scalar(@{$flagged_transcripts});

    # If there are some normal transcripts then assess the flagged ones for removal
    if(scalar(@{$normal_transcripts})) {
      my $db_ids_for_removal = $self->assess_flagged_transcripts($gene,$normal_transcripts,$flagged_transcripts);
      foreach my $db_id (@{$db_ids_for_removal}) {
        push(@{$transcript_ids_to_remove},$db_id);
      }
    }

# NOTE: I'm disabling this code because between genebuilder and the above sub, it shouldn't pick up anything, however
# it could be quite useful in other scenarios, so we can revive it if needed
    # Check the transcripts for redundancy in terms of the terminal exons. This is checked for all transcripts,
    # not just the flagged ones
#    if(scalar(@{$transcripts}) > 1) {
#      my $db_ids_for_removal = $self->assess_transcripts_for_redundancy($transcripts);
#      foreach my $db_id (@{$db_ids_for_removal}) {
#        push(@{$transcript_ids_to_remove},$db_id);
#      }
#    }
#  }

    # This section is for transcripts that are single exon or only frameshift introns
    if(scalar(@{$transcripts}) == 1) {
      my $transcript = shift @{$transcripts};
      if($self->blessed_biotypes()->{$transcript->biotype}) {
        next;
      }

#      if(scalar(@{$transcript->get_all_Exons}) == 1 || $transcript->{'_contains_only_frameshift_introns'}) {
#        # There are a few things we can consider to either keep or delete single exon genes. They are:
#        # 1) biotype, at this point we have the biotypes classified based on hit coverage and percent id, so only take the best ones
#        # 2) utr, if there is utr present then this gives some extra confidence in an model and could override the biotype decision
#        # 3) location, if the model is in the intron of a gene on the same strand then we will almost always want to delete it. The
#        #    only case we wouldn't would be if the intron of the other gene looked dodgy (but this can be tricky to assess)
#
 #       if($self->single_exon_within_intron($transcript)) {
 #         say "Found single exon or frameshift intron only gene within another gene, will remove: ".$transcript->dbID.", ".$transcript->biotype;
 #         push(@{$transcript_ids_to_remove},$transcript->dbID);
 #       }
 #     }
 #   } elsif(($gene->end - $gene->start + 1) <= $tiny_gene_size) {
 #     foreach my $transcript (@{$transcripts}) {
 #       unless(scalar(@{$transcript->get_all_five_prime_UTRs}) || scalar(@{$transcript->get_all_three_prime_UTRs})) {
 #         say "Found tiny gene that does not have a 95/95 alignment: ".$transcript->dbID.", ".$transcript->biotype;
 #         push(@{$transcript_ids_to_remove},$transcript->dbID);
 #       }
 #     }

      my $logic_name = $transcript->analysis->logic_name();
      if($logic_name =~ /^genblast/ && scalar(@{$transcript->get_all_Exons()}) <= 3) {
        if($self->check_protein_models($transcript)) {
          say "Found protein model with weak supporting evidence: ".$transcript->dbID.", ".$transcript->biotype;
          push(@{$transcript_ids_to_remove},$transcript->dbID);
        }
      }
    }

    if($gene_strings->{$gene_string}) {
      say "Found a duplicate gene string:";
      say "Duplicate start: ".$gene->start;
      say "Duplicate end: ".$gene->end;
      say "Duplicate strand: ".$gene->strand;
      say "Duplicate name: ".$gene->seq_region_name;
    } else {
      $gene_strings->{$gene_string} = 1;
    }

  } # End foreach my $gene (@{$genes})
  return($transcript_ids_to_remove);
}

sub flag_transcripts {
  my ($self,$transcripts,$gene_string) = @_;
  my $flagged_transcripts = [];
  my $normal_transcripts = [];
  foreach my $transcript (@{$transcripts}) {
    my $exons = $transcript->get_all_Exons();
    $gene_string .= ":".$self->generate_exon_string($exons);
    if($transcript->translation) {
      $gene_string .= $transcript->translation->seq;
    }

    if($self->blessed_biotypes()->{$transcript->biotype}) {
      next;
    }


    my $introns = $transcript->get_all_Introns();
    my $translateable_seq = $transcript->translateable_seq();
    my $start_codon  = uc( substr( $translateable_seq, 0, 3 ) );
    my $end_codon  = uc( substr( $translateable_seq, -3 ) );
    my $found_frameshift_intron = 0;
    my $found_non_canonical_splice = 0;
    my $frame_shift_intron_count = 0;
    foreach my $intron (@{$introns}) {
      if($intron->length < 10) {
        $found_frameshift_intron = 1;
        $frame_shift_intron_count++;
      } elsif($intron->is_splice_canonical() == 0) {
        $found_non_canonical_splice = 1;
      }
    }

    if($found_frameshift_intron) {
      $transcript->{'_is_frameshift'} = 1;
    }

    if($frame_shift_intron_count && ($frame_shift_intron_count == scalar(@{$introns}))) {
      $transcript->{'_contains_only_frameshift_introns'} = 1;
    }

    if($found_non_canonical_splice) {
      $transcript->{'_is_non_canonical'} = 1;
    }

    unless($start_codon eq "ATG" && ($end_codon eq 'TAG') || ($end_codon eq 'TAA') || ($end_codon eq 'TGA')) {
      $transcript->{'_non_start_stop_complete'} = 1;
    }


    if($found_frameshift_intron || $found_non_canonical_splice) {
      push(@{$flagged_transcripts},$transcript);
    } else {
      push(@{$normal_transcripts},$transcript);
    }
  }

  return([$normal_transcripts,$flagged_transcripts]);

}

sub assess_flagged_transcripts {
  my ($self,$gene,$normal_transcripts,$flagged_transcripts) = @_;

  my $redundancy_cutoff = $self->param('flagged_redundancy_coverage_threshold');
  unless(defined($redundancy_cutoff)) {
    $redundancy_cutoff = 95;
  }

  my $db_ids_to_remove = [];
  foreach my $flagged_transcript (@{$flagged_transcripts}) {
    # If there's a frameshift issue first check if any normal transcripts have the same cds intron
    # set if we remove the frameshift introns, if this is the case just remove it
    if($flagged_transcript->{'_is_frameshift'}) {
      my $redundant_ignoring_frameshifts = $self->check_frameshift_redundancy($flagged_transcript,$normal_transcripts);
      if($redundant_ignoring_frameshifts) {
        say "Found a flagged transcript whose CDS is redundant if you ignore the frameshift introns, dbID: ".$flagged_transcript->dbID;
        $flagged_transcript->{'_marked_for_removal'} = 1;
        push(@{$db_ids_to_remove},$flagged_transcript->dbID);
        next;
      }
    }

    # At this point we want to check how much info would be lost by removing the flagged transcript
    # There are a couple of ways to look at it, the first is the contribution of the transcript to
    # the pool of exons in the normal transcripts. The second way is to look to see if it has an
    # effect on the structure. On one hand we don't want things that are effectively joining two
    # separate genes. But on the other hand we don't want to remove things that repreaent the dominant
    # transcript. The easiest call at the moment might be to look to see how many exons overlap with
    # an ordered set of the ordered transcripts, as soon as you find one with enough of an overlap
    # mark the transcript for deletion
    my @sorted_normal_transcripts = sort { $b->length <=> $a->length } @{$normal_transcripts};
    foreach my $normal_transcript (@sorted_normal_transcripts) {
      my $coverage =  $self->calculate_coverage($normal_transcript,$flagged_transcript);
      if($coverage >= $redundancy_cutoff) {
        say "Found a flagged transcript (".$flagged_transcript->dbID.") that is covered by a normal transcript (".$normal_transcript->dbID.")";
        say "Coverage: ".$coverage." (coverage threshold: ".$redundancy_cutoff.")";
        $flagged_transcript->{'_marked_for_removal'} = 1;
        push(@{$db_ids_to_remove},$flagged_transcript->dbID);
        last;
      }
    }
  }
  return($db_ids_to_remove);
}


sub assess_transcripts_for_redundancy {
  my ($self,$transcripts) = @_;

  my $db_ids_to_remove = [];
  for(my $i=0; $i<scalar(@{$transcripts}) - 1; $i++) {
    my $transcript_a = ${$transcripts}[$i];
    # Skip if this is already marked for removal
    if($transcript_a->{'_marked_for_removal'}) {
      next;
    }

    my $introns_a = $transcript_a->get_all_Introns();
    unless($introns_a) {
      next;
    }

    # Check if the transcript has a cleaned intron string in for transcripts with frameshift introns
    my $intron_string_a;
    if($transcript_a->{'_cleaned_intron_string'}) {
      $intron_string_a = $transcript_a->{'_cleaned_intron_string'};
    } else {
      $intron_string_a = $self->generate_intron_string($introns_a);
    }

    for(my $j=$i+1; $j<scalar(@{$transcripts}); $j++) {
      my $transcript_b = ${$transcripts}[$j];
      # Skip if this is already marked for removal
      if($transcript_b->{'_marked_for_removal'}) {
        next;
      }

      my $introns_b = $transcript_b->get_all_Introns();
      unless($introns_b) {
        next;
      }

      # Check if the transcript has a cleaned intron string in for transcripts with frameshift introns
      my $intron_string_b;
      if($transcript_b->{'_cleaned_intron_string'}) {
        $intron_string_b = $transcript_b->{'_cleaned_intron_string'};
      } else {
        $intron_string_b = $self->generate_intron_string($introns_b);
      }
      if($intron_string_a eq $intron_string_b) {
        my $db_id_to_remove = $self->assess_transcript_with_identical_introns($transcript_a,$transcript_b);
        if($db_id_to_remove) {
          push(@{$db_ids_to_remove},$db_id_to_remove);
          if($db_id_to_remove == $transcript_a->dbID) {
            last;
          }
        }
      }
    }
  }

  return($db_ids_to_remove);
}


sub assess_transcript_with_identical_introns {
  my ($self,$transcript_a,$transcript_b) = @_;
  my $selected_id_for_removal = 0;
  my $redundancy_cutoff = $self->param('general_redundancy_coverage_threshold');
  my $translation_a_length = $transcript_a->translation->length;
  my $translation_b_length = $transcript_b->translation->length;

  my $translation_coverage = 0;
  if($translation_a_length > $translation_b_length) {
    $translation_coverage = ($translation_b_length / $translation_a_length) * 100;
  } else {
    $translation_coverage = ($translation_a_length / $translation_b_length) * 100;
  }

  if($translation_coverage >= $redundancy_cutoff) {
    my $transcript_a_length = $transcript_a->length;
    my $transcript_b_length = $transcript_a->length;
    my $transcript_coverage = 0;
    my $initial_db_id;
    if($transcript_a_length > $transcript_b_length) {
      $transcript_coverage = ($transcript_b_length / $transcript_a_length) * 100;
      $initial_db_id = $transcript_a->dbID;
    } else {
      $transcript_coverage = ($transcript_a_length / $transcript_b_length) * 100;
      $initial_db_id = $transcript_b->dbID;
    }

    if($transcript_coverage >= $redundancy_cutoff) {
      my $transcript_a_incomplete = $transcript_a->{'_non_start_stop_complete'};
      my $transcript_b_incomplete = $transcript_b->{'_non_start_stop_complete'};

      if($initial_db_id == $transcript_a->dbID) {
        unless($transcript_a_incomplete && !$transcript_b_incomplete) {
          $selected_id_for_removal = $initial_db_id;
        } else {
          $selected_id_for_removal = $transcript_b->dbID;
        }
      } else {
        unless($transcript_b_incomplete && !$transcript_a_incomplete) {
          $selected_id_for_removal = $initial_db_id;
       } else {
          $selected_id_for_removal = $transcript_a->dbID;
        }
      }
    }
  }

  return($selected_id_for_removal);
}


sub single_exon_within_intron {
  my ($self,$transcript) = @_;

  my $slice = $transcript->feature_Slice();
  my $genes = $slice->get_all_Genes();

  my $overlap_on_same_strand_count = 0;
  foreach my $gene (@{$genes}) {
    # Note that as this is a feature slice, any genes with strand 1 are on the same strand as the transcript
    if($gene->strand == 1) {
      $overlap_on_same_strand_count++;
    }
  }

  say "FM2 overlap_on_same_strand_count: ".$overlap_on_same_strand_count;

  # If there is more than just the parent gene then this single exon model must lie in an intron of another gene
  # Otherwise it would have been part of the overlapping gene
  if($overlap_on_same_strand_count > 1) {
    return(1);
  } else {
    return(0);
  }

}
sub transcript_ids_to_remove {
  my ($self,$transcript_ids_to_remove) = @_;

  if($transcript_ids_to_remove) {
    $self->param('_transcript_ids_to_remove',$transcript_ids_to_remove);
  }

  return($self->param('_transcript_ids_to_remove'));
}


sub check_frameshift_redundancy {
  my ($self,$frameshift_transcript,$normal_transcripts) = @_;
  my @frameshift_transcript_cds_introns = @{$frameshift_transcript->get_all_CDS_Introns()};
  my @cleaned_cds_introns = ();
  my $clean_count = 0;
  foreach my $intron (@frameshift_transcript_cds_introns) {
    if($intron->length >= 10) {
      $cleaned_cds_introns[$clean_count] = $intron;
      $clean_count++;
    }
  }

  my $cleaned_intron_string = $self->generate_intron_string(\@cleaned_cds_introns);
  # Store the string for the transcript redundancy comparisons later
  $frameshift_transcript->{'_cleaned_intron_string'} = $cleaned_intron_string;
  foreach my $normal_transcript (@{$normal_transcripts}) {
    my $normal_cds_introns = $normal_transcript->get_all_CDS_Introns();
    my $normal_intron_string = $self->generate_intron_string($normal_cds_introns);
    if($normal_intron_string eq $cleaned_intron_string) {
      # The CDS structure is redundant except the frameshift introns. Drop this transcript
      return 1;
    }
  }

  return 0;
}


sub calculate_coverage {
  my ($self, $transcript_a, $transcript_b) = @_;

  my $coverage = 0;
  my $overlap = 0;
  foreach my $exon_a (@{$transcript_a->get_all_Exons()}) {
    foreach my $exon_b (@{$transcript_b->get_all_Exons()}) {
      $overlap += overlap_length($exon_a, $exon_b);
    }
  }

  $coverage = ($overlap / $transcript_b->length) * 100;

  return $coverage;
}


sub overlap_length {
  my ($feature_a, $feature_b) = @_;

  if (!features_overlap($feature_a, $feature_b)) {
    return 0;
  }

  my $min_end = $feature_a->seq_region_end();
  if ($feature_b->seq_region_end() < $min_end) {
    $min_end = $feature_b->seq_region_end();
  }

  my $max_start = $feature_a->seq_region_start();
  if ($feature_b->seq_region_start() > $max_start) {
    $max_start = $feature_b->seq_region_start();
  }

  return $min_end - $max_start + 1;
}


sub features_overlap {
  my ($feature_a, $feature_b) = @_;

  if(($feature_a->seq_region_start() <= $feature_b->seq_region_end() ) && ($feature_a->seq_region_end() >= $feature_b->seq_region_start())) {
    return 1;
  }

  return 0;
}


sub generate_intron_string {
  my ($self,$intron_array) = @_;

  my $intron_string = "";
  foreach my $intron (@{$intron_array}) {
    my $start = $intron->start();
    my $end = $intron->end();
    $intron_string .= $start."..".$end.":";
  }

  return($intron_string);
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


sub generate_exon_string {
  my ($self,$exon_array) = @_;

  my $exon_string = "";
  foreach my $exon (@{$exon_array}) {
    my $start = $exon->start();
    my $end = $exon->end();
    $exon_string .= $start."..".$end.":";
    print "(".$start."..".$end.")";
  }

  print "\n";

  return($exon_string);
}


sub check_protein_models {
  my ($self,$transcript) = @_;

  if(scalar(@{$transcript->get_all_five_prime_UTRs}) || scalar(@{$transcript->get_all_three_prime_UTRs})) {
    return(0);
  }

  my $supporting_features = $transcript->get_all_supporting_features();
  unless(scalar(@$supporting_features)) {
    say "Transcript has no supporting features to assess, looking at exons";
    my $translation = $transcript->translation->seq;
    my $exons = $transcript->get_all_Exons;
    my $all_exon_supporting_features = [];

    my $coverage;
    my $hit_name;
    # Loop till an exon supporting feature is found
    foreach my $exon (@$exons) {
      my $exon_supporting_features = $exon->get_all_supporting_features;
      if(scalar(@$exon_supporting_features)) {
        $coverage = ${$exon_supporting_features}[0]->hcoverage;
        $hit_name = ${$exon_supporting_features}[0]->hseqname;
        last;
      }
    }

    # If an exon supporting feature is found, base decision on that
    if($coverage) {
      if($coverage < 75 && $translation !~ /^M/) {
       say "Removing transcript ".$transcript->seq_region_name.":".$transcript->seq_region_start.":".$transcript->seq_region_end.":".$hit_name.":".$coverage;
       say "Translation:\n".$translation;
       return(1);
      } else {
        return(0);
      }
    }

    my $biotype = $transcript->biotype;
    $biotype =~ /\_(\d+)$/;
    my $classification = $1;
    if($classification <= 4) {
      return(0);
    } else {
      say "Removing based on classification: ".$biotype;
      return(1);
    }
  }

  my $supporting_feature = shift(@$supporting_features);
  my $coverage = $supporting_feature->hcoverage;
  my $translation = $transcript->translation->seq;
  my $hit_name = $supporting_feature->hseqname;
  if($coverage < 75 && $translation !~ /^M/) {
    say "Removing transcript ".$transcript->seq_region_name.":".$transcript->seq_region_start.":".$transcript->seq_region_end.":".$hit_name.":".$coverage;
    say "Translation:\n".$translation;
    return(1);
  }

  return(0);

}

1;
