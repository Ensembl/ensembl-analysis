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

Bio::EnsEMBL::Analysis::Tools::FeatureFilterOnGenome

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Tools::FeatureFilterOnGenome;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils qw(cluster_AlignFeatures);

use parent 'Bio::EnsEMBL::Analysis::Tools::FeatureFilter';


=head2 filter_results

 Arg [1]    : Arrayref of Bio::EnsEMBL::BaseAlignFeature, Bio::EnsEMBL::DnaPepAlignFeature or Bio::EnsEMBL::DnaDnaAlignFeature
 Description: Cluster the hits, filter the unclustered hits based on the parameters given to the filter, no pruning necessary.
              Then filter the clustered hits and prune if needed
 Returntype : Arrayref Bio::EnsEMBL::BaseAlignFeature
 Exceptions : Throws if Arg[1] is not an arrayref
              Throws if Arg[1] is emtpy

=cut

sub filter_results {
  my ($self, $features) = @_;
  if (!$features || ref($features) ne 'ARRAY') {
    throw("Must pass filter_results an arrayref not ".$features.
          " FeatureFilter::filter_results");
  }
  my %types = map {$_->analysis->logic_name => 1} @$features;
  my ($clusters, $unclustered) = cluster_AlignFeatures($features, {hit => [keys %types]});
  my @features;
  foreach my $cluster (@$unclustered) {
    my ($feature) = @{$cluster->get_AlignFeatures};
    if ($feature->score > $self->min_score) {
      if (!$feature->can('p_value') || !defined $feature->p_value || $feature->p_value < $self->max_pvalue) {
        push(@features, @{$cluster->get_AlignFeatures});
      }
    }
  }
  foreach my $cluster (@$clusters) {
    my %validhit;
    my %hitarray;
    my %totalscore;
    my $minstart;
    my $maxend = 0;

    #filtering by score
    #sorting by score so we use the highest scoring feature first
    #The score filter basically takes all features belonging to
    #one hit id provided that at least one of its features has a score
    #greater than the min score and a pvalue less than the max pvalue
    foreach my $f (sort { $b->score <=> $a->score } @{$cluster->get_AlignFeatures}) {
      if ($f->score > $self->min_score) {
        if (!exists $validhit{$f->hseqname}) {
          $validhit{$f->hseqname} = 0;
          $totalscore{$f->hseqname} = 0;
        }
        $totalscore{$f->hseqname} += $f->score;
        if (!$f->can('p_value') || !defined $f->p_value || $f->p_value < $self->max_pvalue) {
          if ($validhit{$f->hseqname} < $f->score) {
            $validhit{$f->hseqname} = $f->score;
          }
          if ($f->end > $maxend) {
            $maxend = $f->end;
          }
          if (!$minstart or $minstart > $f->start) {
            $minstart = $f->start;
          }
        }
      }
      #if a hit doesn't pass the min score or max pvalue threshold
      #but the hit id has features which do this ensures all features
      #for that hit id are kept regardless of score and pvalue
      if ($validhit{$f->hseqname}) {
        if (!exists($hitarray{$f->hseqname})) {
          $hitarray{$f->hseqname} = [];
        }
        push(@{$hitarray{$f->hseqname}}, $f);
      }
    }

    my @hit_ids = sort {$validhit{$b} <=> $validhit{$a}
                           or $totalscore{$b} <=> $totalscore{$a}
                             or $a cmp $b
                           } keys %validhit;

    #this sorts the hseqnames of the valid hits first by the max score
    #then by the highest total score
    #then on alphabetical order before feeding the names to the
    #coverage based filter
    my %accepted_hit_ids;
    #This coverage filter works on principle if one hit
    #for a particular id is lower than the coverage it takes all the hits
    #It only considers hits who have scores and p values with in the limits
    #First it generates an array which is the length of the max end of all
    #the features is. Then it considers each strand separately.
    #Then for each hit name it goes though each features if the feature is
    #on the appropriate strand it looks at the array index for each base pair
    #in that feature and checks if the coverage is too high. Provided one
    #basepair has less coverage than required the feature is kept
    #otherwise that hit id is marked to be thrown away.

    if ($self->filter_on_coverage) {
      foreach my $strand(1, -1) {
        my @list;
        $list[$maxend-$minstart] = 0; #perl will automatically extend this array
        foreach my $name(@hit_ids) {
          my $hole = 0;
          foreach my $f (@{$hitarray{$name}}) {
            next if ($f->strand != $strand);
            if ($f->score > $self->min_score and $f->can('p_value') and defined $f->p_value and $f->p_value < $self->max_pvalue) {
              foreach my $i ( ($f->start-$minstart) .. ($f->end-$minstart) ) {
                unless( $list[$i] ) {
                  $list[$i] = 0;
                }
                if ( $list[$i] < $self->coverage ) {
                  # accept!
                  $hole = 1;
                  last;
                }
              }
            }
          }
          if ($hole == 0) {
            next ;
          }
          $accepted_hit_ids{$name} = 1;
          foreach my $f ( @{$hitarray{$name}} ) {
            if ($f->strand == $strand) {
              for my $i ( ($f->start-$minstart) .. ($f->end-$minstart) ) {
                $list[$i]++;
              }
            }
          }
        }
      }
    }
    else{
      foreach my $name(keys(%hitarray)) {
        $accepted_hit_ids{$name} = 1;
      }
    }

    my @pruned_features;
    foreach my $name(keys(%accepted_hit_ids)) {
      my @tmp = @{$hitarray{$name}};

      my $remaining;
      if ($self->prune) {
        $remaining = $self->prune_features(\@tmp);
      }else{
        $remaining = \@tmp;
      }

      push(@pruned_features, @$remaining);
    }

    if ($self->hard_prune) { #this pruning is done across all features
      @pruned_features = @{$self->prune_features(\@features)};
    }
    push(@features, @pruned_features);

  }

  return \@features;
}


=head2 prune_features_by_strand

 Arg [1]    : Bio::EnsEMBL::Analysis::Tools::FeatureFilter
 Arg [2]    : int, strand (must be 1 or -1)
 Arg [3]    : arrayref of features
 Description: Works out coverage of each base pair and throws out features of a base pair is covered by to many hits.
              The lowest scoring features are thrown out first
              1. first create an array of all features of the appropriate strand

              2. Then create a separate array whose last index is the last base pair
              covered by any features

              3. go through that array and increment each index for each feature which
              covers that base pair

              4.Then it cycles through the bases covered array collecting a list of
              base pairs whose coverage is above the limit

              5. the over covered list is sorted so the bases which have the greatest
              excess coverage are considered first

              6. for each over covered base pair all all of the features are looked at
              each time a feature covers this base pair it is thrown away and the
              coverage of each base the rejected feature covers is decremented.

              7. Each over covered base is considered untill the coverage is decreased
              to acceptable limits

              8. the remaining feature set is returned
 Returntype: arrayref of features
 Exceptions: none

=cut


sub prune_features_by_strand {
  my ($self, $strand, $in) = @_;

  my @input_for_strand;
  foreach my $f (@$in) {
    push @input_for_strand, $f if $f->strand eq $strand;
  }

  return if (!@input_for_strand);

  my ($first_object) = sort{ $a->start <=> $b->start } @input_for_strand;
  my $first_base = $first_object->start;
  my ($last_object) = sort{ $b->end <=> $a->end } @input_for_strand;
  my $length = $last_object->end-$first_base;

  # fs_per_base: set element i to the number of features covering base i
  my @fs_per_base;
  foreach  my $base (0..$length) {
    $fs_per_base[$base] = 0;	# initialise
  }
  foreach my $f (@input_for_strand) {
    foreach my $covered_base (($f->start-$first_base)..($f->end-$first_base)) {
      $fs_per_base[$covered_base]++;
    }
  }

  # put the worst features first, so they get removed with priority
  my @sorted_fs = sort {$a->score <=> $b->score || $a->start <=> $b->start || $b->end <=> $a->end} @input_for_strand;

  # over_covered_bases: list of base numbers where coverage must be
  # reduced, listed worst-case-first
  my $max_coverage = $self->coverage;
  my @over_covered_bases;
  foreach  my $base (0..$length) {
    if (($fs_per_base[$base] - $max_coverage) > 0) {
      push @over_covered_bases, $base;
    }
  }

  foreach my $base (sort { $fs_per_base[$b] <=> $fs_per_base[$a] } @over_covered_bases) {
    my $f_no = 0;
    while ($fs_per_base[$base] > $max_coverage) {
      my $start = $sorted_fs[$f_no]->start-$first_base;
      my $end = $sorted_fs[$f_no]->end-$first_base;
      if ($start <= $base and $end >= $base) {	# cut this feature
        splice @sorted_fs, $f_no, 1;	# same index will give next feature
          foreach my $was_covered ($start..$end) {
            $fs_per_base[$was_covered]--;
          }
      } else {	# didn't overlap this base, move on to next feature
        $f_no++;
      }
    }
  }
  return \@sorted_fs;
}

1;
