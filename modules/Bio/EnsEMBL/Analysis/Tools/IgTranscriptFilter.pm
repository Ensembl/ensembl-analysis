# Ensembl module for Bio::EnsEMBL::Analysis::Tools::IgTranscriptFilter
#
# Copyright (c) 2004 Ensembl
#


# $Source: /tmp/ENSCOPY-ENSEMBL-ANALYSIS/modules/Bio/EnsEMBL/Analysis/Tools/IgTranscriptFilter.pm,v $
# $Revision: 1.3 $
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
