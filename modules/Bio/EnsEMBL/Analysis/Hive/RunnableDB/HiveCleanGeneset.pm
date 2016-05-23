#!/usr/bin/env perl

# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCleanGeneset;

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

  my $input_dba = $self->hrdb_get_dba($self->param('input_db'));
  $self->hrdb_set_con($input_dba,'input_db');
  $input_dba->dnadb($dna_dba);


  my $input_id = $self->param('iid');
  my $input_id_type = $self->param('iid_type');

  # If the input is an array of gene ids then loop through them and fetch them from the input db
  if($input_id_type eq 'slice') {
    $self->load_genes($input_id);
  } else {
    $self->throw("You have selected an input id type using iid_type that is not supported by the module. Offending type: ".$input_id_type);
  }

  my $output_path = $self->param('output_path');
  unless(-e $output_path) {
    system("mkdir -p ".$output_path);
  }

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
  my $output_file = $self->param('output_path')."/ids_to_remove.".$$.".txt";
  open(OUT,">".$output_file);
  foreach my $id (@{$transcript_ids_to_remove}) {
    say OUT $id;
  }
  close OUT;

  return 1;
}

sub load_genes {
  my ($self,$input_id) = @_;

  my $dba = $self->hrdb_get_con('input_db');
  my $slice_adaptor = $dba->get_SliceAdaptor();
  my $slice = $slice_adaptor->fetch_by_name($input_id);
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

sub clean_genes {
  my ($self,$genes) = @_;
  my $transcript_ids_to_remove = [];

  say "Have ".scalar(@{$genes})." genes to process";
  foreach my $gene (@{$genes}) {
    my $transcripts = $gene->get_all_Transcripts();
    my ($normal_transcripts,$flagged_transcripts) = @{$self->flag_transcripts($transcripts)};
    say "Checking gene with dbID: ".$gene->dbID;
    say "Total transcript count: ".scalar(@{$transcripts});
    say "Normal transcript count: ".scalar(@{$normal_transcripts});
    say "Flagged transcript count: ".scalar(@{$flagged_transcripts});

    if(scalar(@{$transcripts}) == scalar(@{$flagged_transcripts})) {
      say "All transcripts are flagged so no cleaning will take place";
      next;
    } else {
      my $db_ids_for_removal = $self->access_flagged_transcripts($gene,$normal_transcripts,$flagged_transcripts);
      foreach my $db_id (@{$db_ids_for_removal}) {
        push(@{$transcript_ids_to_remove},$db_id);
      }
    }
  }

  return($transcript_ids_to_remove);
}

sub flag_transcripts {
  my ($self,$transcripts) = @_;
  my $flagged_transcripts = [];
  my $normal_transcripts = [];
  foreach my $transcript (@{$transcripts}) {
    my $introns = $transcript->get_all_Introns();
    my $found_frameshift_intron = 0;
    my $found_non_canonical_splice = 0;
    foreach my $intron (@{$introns}) {
      if($intron->length < 10) {
        $found_frameshift_intron = 1;
      } elsif($intron->is_splice_canonical() == 0) {
        $found_non_canonical_splice = 1;
      }
    }

    if($found_frameshift_intron) {
      $transcript->{'_is_frameshift'} = 1;
    }

    if($found_non_canonical_splice) {
      $transcript->{'_is_non_canonical'} = 1;
    }

    if($found_frameshift_intron || $found_non_canonical_splice) {
      push(@{$flagged_transcripts},$transcript);
    } else {
      push(@{$normal_transcripts},$transcript);
    }
  }

  return([$normal_transcripts,$flagged_transcripts]);

}

sub access_flagged_transcripts {
  my ($self,$gene,$normal_transcripts,$flagged_transcripts) = @_;

  my $redundancy_cutoff = $self->param('redundancy_coverage_threshold');
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
        push(@{$db_ids_to_remove},$flagged_transcript->dbID);
        last;
      }
    }
  }
  return($db_ids_to_remove);
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


1;
