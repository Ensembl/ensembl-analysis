#!/usr/bin/env perl

# Copyright [2017-2022] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck;

use strict;
use warnings;
use feature 'say';

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    min_size_empty_slice => 1000000,
    min_allowed_feature_counts => undef,
    skip_rnaseq => 0,
    skip_projection => 0,
  }
}

=head2 fetch_input

 Arg [1]    : None
 Description: Use skip_analysis to avoid running the sanity checks
 Returntype : None
 Exceptions : None

=cut

sub fetch_input {
  my $self = shift;

  if ($self->param('skip_analysis')) {
    $self->complete_early('Skip check flag is enabled, so no check will be carried out');
  }
  else {
    $self->param_required('sanity_check_type');
    my $target_db = $self->get_database_by_name('target_db');
    $self->hrdb_set_con($target_db, 'target_db');
  }

  return 1;
}


sub run {
  my $self = shift;

  $self->run_sanity_checks($self->param('sanity_check_type'));

  return 1;
}


sub write_output {
  my ($self) = @_;

  return 1;
}


sub run_sanity_checks {
  my ($self,$type) = @_;

  if($type eq 'genome_preparation_checks') {
    $self->post_genome_preparation();
  } elsif($type eq 'gene_db_checks') {
    $self->gene_db_checks();
  } elsif($type eq 'final_core_checks') {
    $self->final_core_checks();
  } else {
    $self->throw('The sanity check type you specified is not recognised. Type: '.$type);
  }
}


sub post_genome_preparation {
  my ($self) = @_;

  my $logic_names = $self->param('min_allowed_feature_counts');
  if ($logic_names) {
    my $genome_prep_dba = $self->hrdb_get_con('target_db');

    my $feature_adaptors = {
                             'repeat' => $genome_prep_dba->get_RepeatFeatureAdaptor(),
                             'simple' => $genome_prep_dba->get_SimpleFeatureAdaptor(),
                             'prediction transcript' => $genome_prep_dba->get_PredictionTranscriptAdaptor(),
                             'dna align' => $genome_prep_dba->get_DnaAlignFeatureAdaptor(),
                             'protein align' => $genome_prep_dba->get_ProteinAlignFeatureAdaptor(),
                           };

    my $previous_feature_type = "";
    my $error_string = '';
    foreach my $logic_name (sort {$logic_names->{$a}->[1] cmp $logic_names->{$b}->[1]} keys %{$logic_names}) {
      my ($min_count,$feature_type) = @{$logic_names->{$logic_name}};
      if($feature_type ne $previous_feature_type) {
        say "\nChecking ".$feature_type." features";
        $previous_feature_type = $feature_type;
      }
      my $observed_count = $self->count_features($logic_name,$feature_adaptors->{$feature_type});
      say "Count for ".$logic_name.": ".$observed_count;
      unless($observed_count >= $min_count) {
        $error_string .= "Observed value was too low for $logic_name, min allowed: $min_count, observed: $observed_count\n";
      }
    }
    if ($error_string) {
      $self->throw($error_string);
    }
  }
}

sub gene_db_checks {
  my ($self) = @_;

  my $genome_prep_dba = $self->hrdb_get_con('target_db');
  my $min_allowed_feature_counts = $self->param('min_allowed_feature_counts');
  my $slices = $genome_prep_dba->get_SliceAdaptor->fetch_all('toplevel');
  if ($min_allowed_feature_counts) {

    my $logic_names = $min_allowed_feature_counts->{'logic_names'};
    my $biotypes = $min_allowed_feature_counts->{'biotypes'};
    unless($biotypes) {
      $self->throw("No entries for biotypes where found in the min_allowed_feature_counts hash");
    }
    unless($logic_names) {
        $self->throw("No entries for logic_names where found in the min_allowed_feature_counts hash");
    }

    my $observed_biotype_counts = {};
    my $observed_logic_name_counts = {};
    foreach my $slice (@$slices) {
      foreach my $gene (@{$slice->get_all_Genes}) {
        foreach my $transcript (@{$gene->get_all_Transcripts}) {
          my $biotype = $transcript->biotype;
          if (exists $observed_biotype_counts->{$biotype}) {
            $observed_biotype_counts->{$biotype}++;
          }
          else {
            $observed_biotype_counts->{$biotype} = 1;
          }
        }
        my $logic_name = $gene->analysis->logic_name;
        if (exists $observed_logic_name_counts->{$logic_name}) {
          $observed_logic_name_counts->{$logic_name}++;
        }
        else {
          $observed_logic_name_counts->{$logic_name} = 1;
        }
      }
    }

    $self->say_with_header("Checking logic names:");
    foreach my $logic_name (keys(%{$logic_names})) {

      if ($self->param('skip_rnaseq') && $logic_name =~ /_rnaseq/) {
        $self->warning("Skipping count of $logic_name as rnaseq db flagged as absent");
      }
      elsif ($self->param('skip_projection') && $logic_name =~ /^project_/) {
        $self->warning("Skipping count of $logic_name as projection db flagged as absent");
      }
      else {
        my $min_count = $logic_names->{$logic_name};
        my $observed_count = $observed_logic_name_counts->{$logic_name};
        $self->say_with_header("Observed gene/transcript count for $logic_name : $observed_count");
        if ($observed_count < $min_count) {
          $self->throw("Gene/transcript count is below the min value of $min_count for $logic_name, observed count: $observed_count");
        }
      }

    }

    # Count the transcripts per biotype
    # Not particularly happy with this, I think it is good to have a way to specify generic parts of a biotype, but it comes at the
    # expense of having to pattern match, which is dangerous
    $self->say_with_header("Checking biotypes:");
    my $observed_biotype_keys = qr/${\ join('|', map quotemeta, keys %$observed_biotype_counts) }/;
    foreach my $biotype (keys(%{$biotypes})) {
      if ($self->param('skip_rnaseq') && $biotype =~ /^rnaseq\_/) {
        $self->warning("Skipping count of $biotype as rnaseq db flagged as absent");
      }
      elsif ($self->param('skip_projection') && $biotype =~ /^realign\_/) {
        $self->warning("Skipping count of $biotype as projection db flagged as absent");
      }
      else {
        my $min_count = $biotypes->{$biotype};
        my $observed_count = 0;
        while($observed_biotype_keys =~ s/$biotype[^|\)]*//) {
          my $matching_biotype = $&;
          my $matching_biotype_count = $observed_biotype_counts->{$matching_biotype};
          $observed_count += $matching_biotype_count;
        }

        $self->say_with_header("Observed gene/transcript count for $biotype : $observed_count");
        unless(defined $observed_count && ($observed_count >= $min_count)) {
          $self->throw("Gene/transcript count is below the min value of ".$min_count." for ".$biotype.", observed count: ".$observed_count);
        }
      }
    }
  }
  else {
    my $gene_adaptor = $genome_prep_dba->get_GeneAdaptor;
    my $min_size_empty_slice = $self->param('min_size_empty_slice');
    $self->warning('No count provided, will fail if no gene is found or if no gene is found in region bigger than '.$min_size_empty_slice);
    foreach my $slice (@$slices) {
      my $count = $gene_adaptor->count_all_by_Slice($slice);
      if ($count == 0 and $slice->length > $min_size_empty_slice) {
        $self->throw($slice->name.' has no genes');
      }
    }
  }
}

sub final_core_checks {
  my ($self) = @_;
  # To be implemented
  # Should cover most of the common healthchecks. Should check that all the logic names and external db ids make sense across the different
  # table types
  # After this point we should be certain that it will pass the set of critical healtchecks that production will require for handover
}

sub count_features {
  my ($self,$logic_name,$adaptor) = @_;

  my $feature_set = $adaptor->fetch_all_by_logic_name($logic_name);
  my $feature_count = scalar(@{$feature_set});
  $feature_set = [];
  return($feature_count);
}

1;

