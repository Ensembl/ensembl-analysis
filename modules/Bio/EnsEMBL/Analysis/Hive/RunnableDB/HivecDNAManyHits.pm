#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivecDNAManyHits;

use strict;
use warnings;
use feature 'say';


use Bio::EnsEMBL::Analysis::Tools::Utilities qw(hrdb_get_dba);
use parent ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');


sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    threshold => 20,
    column_names => ['iid'],
    many_hits_process_threshod => .90,
  }
}

sub fetch_input {
  my $self = shift;

  my $db = hrdb_get_dba($self->param_required('target_db'));
  my $slice_adaptor = $db->get_SliceAdaptor;
  my %hit_names;
  foreach my $slice (@{$slice_adaptor->fetch_all('toplevel', undef, 1)}) {
    foreach my $transcript (@{$slice->get_all_Transcripts}) {
      ++$hit_names{$transcript->get_all_supporting_features->[0]->hseqname};
    }
  }
  my @many_hits;
  my $threshold = $self->param('threshold');
  foreach my $key (keys %hit_names) {
    push(@many_hits, $key) if ($hit_names{$key} > $threshold);
  }
  if (@many_hits) {
    if ($self->param_is_defined('old_db')) {
      my $old_db = hrdb_get_dba($self->param_required('old_db'));
      my $transcript_adaptor = $old_db->get_TranscriptAdaptor;
      my @to_process;
      $threshold *= $self->param('many_hits_process_threshod');
      foreach my $hitname (@many_hits) {
        my $transcripts = $transcript_adaptor->fetch_all_by_transcript_supporting_evidence($hitname, 'dna_align_feature');
        push(@to_process, $hitname) unless (scalar(@$transcripts) > $threshold);
      }
      if (@to_process) {
        $self->param('inputlist', \@to_process);
      }
      else {
        $self->complete_early(scalar(@many_hits).' cDNAs had more than '.$self->param('threshold').' hits but were already in the previous database');
      }
    }
    else {
      $self->param('inputlist', \@many_hits);
    }
  }
  else {
    $self->complete_early("No cDNAs had more than $threshold hits");
  }
}


1;
