#!/usr/bin/env perl

# Copyright [2017] EMBL-European Bioinformatics Institute
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

sub fetch_input {
  my $self = shift;

  $self->param_required('sanity_check_type');

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

  if($type eq 'post_genome_preparation') {
    $self->post_genome_preparation();
  } else {
    $self->throw('The sanity check type you specified is not recognised. Type: '.$type);
  }
}


sub post_genome_preparation {
  my ($self) = @_;

  # ToDo, add in a hash containing min counts based on average counts for other species in the clade
  my $logic_names = $self->param_required('feature_logic_names');
#  my $simple_logic_names = $self->param_required('simple_feature_logic_names');
#  my $daf_logic_names = $self->param_required('dna_align_feature_logic_names');
#  my $paf_logic_names = $self->param_required('protein_align_feature_logic_names');
#  my $prediction_logic_names = $self->param_required('prediction_transcript_logic_names');

  my $target_db = $self->param_required('target_db');
  my $genome_prep_dba = $self->hrdb_get_dba($target_db);

  my $feature_adaptors = {
                           'repeat' => $genome_prep_dba->get_RepeatFeatureAdaptor(),
                           'simple' => $genome_prep_dba->get_SimpleFeatureAdaptor(),
                           'prediction' => $genome_prep_dba->get_PredictionTranscriptAdaptor(),
                           'dna align' => $genome_prep_dba->get_DnaAlignFeatureAdaptor(),
                           'protein align' => $genome_prep_dba->get_ProteinAlignFeatureAdaptor(),
                         };
#  my $repeat_feature_adaptor = $genome_prep_dba->get_RepeatFeatureAdaptor();
#  my $simple_feature_adaptor = $genome_prep_dba->get_SimpleFeatureAdaptor();
#  my $daf_adaptor = $genome_prep_dba->get_DnaAlignFeatureAdaptor();
#  my $paf_adaptor = $genome_prep_dba->get_ProteinAlignFeatureAdaptor();
#  my $prediction_transcript_adaptor = $genome_prep_dba->get_PredictionTranscriptAdaptor();

  my $feature_counts = {};
  my $previous_feature_type = "";
  foreach my $logic_name (sort {$logic_names->{$a} cmp $logic_names->{$b}} keys %{$logic_names}) {
    my $feature_type = $logic_names->{$logic_name};
    my $feature_adaptor = $feature_adaptors->{$feature_type};
    if($feature_type ne $previous_feature_type) {
      say "\nChecking ".$feature_type." features";
      $previous_feature_type = $feature_type;
    }
    $self->count_features($logic_name,$feature_adaptor,$feature_counts);
  }

#  say "Checking repeat features";
#  my $repeat_counts = $self->count_features($repeat_logic_names,$repeat_feature_adaptor);
#  say "\nChecking simple features";
#  my $simple_feature_counts = $self->count_features($simple_logic_names,$simple_feature_adaptor);
#  say "\nChecking dna align features";
#  my $daf_counts = $self->count_features($daf_logic_names,$daf_adaptor);
#  say "\nChecking protein align features";
#  my $paf_counts = $self->count_features($paf_logic_names,$paf_adaptor);
#  say "\nChecking prediction transcript features";
#  my $prediction_transcript_counts = $self->count_features($prediction_logic_names,$prediction_transcript_adaptor);

  if($self->param('min_allowed_feature_counts')) {
#    my $all_counts = {%$repeat_counts,%$simple_feature_counts,%$daf_counts,%$paf_counts,%$prediction_transcript_counts};
#    $self->check_min_counts($all_counts);
    $self->check_min_counts($feature_counts);
  }
}

sub count_features {
  my ($self,$logic_name,$adaptor,$feature_counts) = @_;

  my $feature_set = $adaptor->fetch_all_by_logic_name($logic_name);
  my $feature_count = scalar(@{$feature_set});
  unless($feature_count) {
    $self->throw("Feature count is 0 for ".$logic_name);
  } else {
    say "Feature count for ".$logic_name.": ".$feature_count;
  }

  $feature_counts->{$logic_name} = $feature_count;
}

#sub count_features {
#  my ($self,$logic_names,$adaptor) = @_;
#  my $feature_counts = {};
#  foreach my $logic_name (keys(%{$logic_names})) {
#    my $feature_set = $adaptor->fetch_all_by_logic_name($logic_name);
#    my $feature_count = scalar(@{$feature_set});
#    unless($feature_count) {
#      $self->throw("Feature count is 0 for ".$logic_name);
#    } else {
#      say "Feature count for ".$logic_name.": ".$feature_count;
#    }
#    $feature_counts->{$logic_name} = $feature_count;
#  }
#  return($feature_counts);
#}


sub check_min_counts {
  my ($self,$observed_counts) = @_;
  my $min_allowed_counts = $self->param('min_allowed_feature_counts');
  foreach my $logic_name (keys(%$min_allowed_counts)) {
    my $expected_value = $min_allowed_counts->{$logic_name};
    unless(exists $observed_counts->{$logic_name}) {
      $self->throw("A logic name from the min counts hash was not found in the observed hash. Logic name: ".$logic_name);
    }
    my $observed_value = $observed_counts->{$logic_name};
    unless($observed_value >= $expected_value) {
      $self->throw("Observed value was too low for ".$logic_name.", min allowed: ".$expected_value.", observed: ".$observed_value);
    }
  }
}


1;
