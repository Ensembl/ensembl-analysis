# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Tools::IgTranscriptFilter;

use strict;
use warnings;

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );



sub new{
  my ($class, @args) = @_;
  my $self = bless {},$class;

  my ($percent_range,
      $min_percent,
      $min_score) = 
          rearrange([
                     'PERCENTRANGE',
                     'MINPERCENT',
                     'MINSCORE'
                     ], @args);

  $self->percent_range($percent_range) if defined $percent_range;
  $self->min_percent($min_percent) if defined $min_percent;
  $self->min_score($min_score) if defined $min_score;

  return $self;
}

sub filter_results{
  my ($self, $transcripts) = @_;

  #############
  # assumption: transcripts are sorted by score, best to worst.
  # This is how ExonerateTranscript returns them.
  #############

  my (@keep, %by_hit);

  foreach my $t (@$transcripts) {
    my ($sf) = @{$t->get_all_supporting_features};

    next if $sf->percent_id < $self->min_percent;
    next if $sf->score < $self->min_score;

    push @{$by_hit{$sf->hseqname}}, [$t, $sf];
  }

  foreach my $hid (keys %by_hit) {
    my @entries = @{$by_hit{$hid}};

    # only consider hits that have coverage at least as good as first,
    # and with percent_id within percent_range of the first

    my $best_cov = $entries[0]->[1]->score;
    my $best_percent  = $entries[0]->[1]->percent_id;

    foreach my $hit (@entries) {
      if ($hit->[1]->score >= $best_cov and
          ($hit->[1]->percent_id >= $best_percent or
           $best_percent - $hit->[-1]->percent_id <= $self->percent_range)) {
        push @keep, $hit->[0];
      }
    }

  }

  return \@keep;
}




sub percent_range {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_percent_range} = $val;
  }

  if (exists $self->{_percent_range}) {
    return $self->{_percent_range};
  } else {
    return 0;
  }
}


sub min_percent {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_min_percent} = $val;
  }

  if (exists $self->{_min_percent}) {
    return $self->{_min_percent};
  } else {
    return 0;
  }
}


sub min_score {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_min_score} = $val;
  }

  if (exists $self->{_min_score}) {
    return $self->{_min_score};
  } else {
    return 0;
  }
}



1;
