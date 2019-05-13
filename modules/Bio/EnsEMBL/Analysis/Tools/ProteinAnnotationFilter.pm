# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

=head1 NAME

  Bio::EnsEMBL::Analysis::Tools::ProteinAnnotationFilter

=head1 SYNOPSIS

  my $filter = new Bio::EnsEMBL::Analysis::Tools::ProteinAnnotationFilter
  new->()
  my @filtered_results = @{$filter->filter_results(\@results)};

=head1 DESCRIPTION

This can be used as a default filter for ProteinAnnotation. It
removes overlaps in a "sensible" way. 

=cut


package Bio::EnsEMBL::Analysis::Tools::ProteinAnnotationFilter;

use warnings ;
use strict;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);


sub new {
  my ($caller, @args) = @_;
  my $class = ref($caller) || $caller;

  my $self = bless({}, $class);
  
  my ($sort_func,
      $cut) = rearrange(['sort',
                         'cut',
                               ], @args);

  if (defined $sort_func) {
    $self->sort_function($sort_func);
  }
  if (defined $cut) {
    $self->cut($cut);
  }

  return $self;
}


=head2 filter_results

 Title   : filter_results
 Description: Filters the given features set by basic
  overlap. The scores are used to determine priority 
  (i.e. high-scoring features over-rule low-scoring ones)
  but this can be over-ridden by supplying a sort function
  reference
                                                      
=cut

sub filter_results {
  my ($self, $res) = @_;

  my $sort_routine = $self->sort_function;

  my (%res_by_name, @all_fps);

  foreach my $f (@$res) {
    push @{$res_by_name{$f->seqname}->{$f->hseqname}}, $f;
  }

  foreach my $seqname (keys %res_by_name) {
    # stage one: remove overlapping hits to the same family
    my @kept_hits;
    foreach my $fam (keys %{$res_by_name{$seqname}}) {
      my @hits_so_far;
      my @hits = sort $sort_routine @{$res_by_name{$seqname}->{$fam}};

      foreach my $h (@hits) {
        my $overlap = 0;
        foreach my $kh (@hits_so_far) {
          if ($h->start <= $kh->end and
              $h->end   >= $kh->start) {
            $overlap = 1;
            last;
          }
        }
        if (not $overlap) {
          push @hits_so_far, $h;
        }      
      }

      push @kept_hits, @hits_so_far;
    }
    
    # stage 2: truncate or remove remaining overlaps between families

    @kept_hits = sort $sort_routine @kept_hits;

    if (not $self->cut) {
      my @hits_so_far;

      foreach my $h (@kept_hits) {
        my $overlap = 0;
        foreach my $kh (@hits_so_far) {
          if ($h->start <= $kh->end and
              $h->end   >= $kh->start) {
            $overlap = 1;
            last;
          }
        }
        if (not $overlap) {
          push @hits_so_far, $h;
        }
      }

      push @all_fps, @hits_so_far;
    } else {
      my @hit_string;
      
      foreach my $h (reverse @kept_hits) {
        foreach my $idx ($h->start .. $h->end) {
          $hit_string[$idx] = $h;
        }
      }
      
      my $max = scalar(@hit_string) - 1;
      my @domains;
      
      foreach my $idx(1 .. $max) {
        if (defined($hit_string[$idx])) {
          if (not @domains or
              $domains[-1]->{end} < $idx - 1 or
              $hit_string[$idx] != $domains[-1]->{domain}) {
            push @domains, {
              start => $idx,
              end   => $idx,
              domain => $hit_string[$idx],
            };
          } elsif ($hit_string[$idx] == $domains[-1]->{domain}) {
            $domains[-1]->{end} = $idx;
          }
        } 
      }
      
      my @fps;
      foreach my $dom (@domains) {
        my $f = $dom->{domain};
        my $st = $dom->{start};
        my $en = $dom->{end};
        
        if ($st == $f->start and
            $en   == $f->end) {
          # the feature did not have to be truncated
          push @fps, $f;
        } else {
          # truncation; create a new version of the feature
          my $newfp = Bio::EnsEMBL::ProteinFeature->new(-seqname => $f->seqname,
                                                        -start   => $st,
                                                        -end     => $en,
                                                        -hseqname => $f->hseqname,
                                                        -hstart  => $f->hstart,
                                                        -hend    => $f->hend,
                                                        -score   => $f->score,
                                                        -percent_id => $f->percent_id,
                                                        -analysis => $f->analysis);
          push @fps, $newfp;
          
        }
      }
      
      push @all_fps, @fps;
    }
  }

  return \@all_fps;
}


sub cut {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_cut} = $val;
  }

  return $self->{_cut};
}


sub sort_function {
  my ($self, $func) = @_;

  if (defined $func) {
    $self->{_sort_func} = $func;
  }
  if (defined $self->{_sort_func}) {
    return $self->{_sort_func};
  } else {
    my $ref = sub { $b->score <=> $a->score };
    return $ref;
  }
}

1;
