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

=head1 NAME

  Bio::EnsEMBL::Analysis::Tools::ClusterFilter

=head1 SYNOPSIS

  my $filter = new Bio::EnsEMBL::Analysis::Tools::ClusterFilter
  new->(
        -coverage => 80,
        -percent_id => 90,
        -best_in_genome => 1,
       );

  my @filtered_results = @{$filter->filter_results(\@results)};

=head1 DESCRIPTION

This is the standard module used for filtering WGA2GenesDirect transcripts. It takes all the transcript alignments riginated from one gene
and cluster them together by location, it then selects the best cluster based in coverage and percent id of the transcripts within each cluster.
It was originally used to avoid transcript from one gene to be allocated in different chromosomes in the projection process (which would have no sense).

=cut


package Bio::EnsEMBL::Analysis::Tools::ClusterFilter;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );



=head2 new

  Returntype: Bio::EnsEMBL::Analysis::Tools::ClusterFilter
  Exceptions: none
  Example   : 

=cut



sub new{
  my ($class,@args) = @_;
  my $self = bless {}, $class;
  &verbose('WARNING');
  my ($min_coverage,
      $min_percent,
      $max_stops,
      $best_in_genome,
      ) = 
        rearrange(['coverage',
                   'percent_id',
                   'max_stops',
                   'best_in_genome'], @args); 

  ######################
  #SETTING THE DEFAULTS#
  ######################

  $self->min_coverage($min_coverage) if defined $min_coverage;
  $self->min_percent($min_percent) if defined $min_percent;
  $self->max_stops($max_stops) if defined $max_stops;
  $self->best_in_genome($best_in_genome) if $best_in_genome;

  return $self;
}


#filter methods



=head2 filter_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::ClusterFilter
  Arg [2]   : arrayref of Trancripts
  Function  : filter the given Transcruipts in the tried and trusted manner
  Returntype: arrayref
  Exceptions: throws if passed nothing or not an arrayref
  Example   : 

=cut


sub filter_results{
  my ($self, $transcripts) = @_;
  
  # results are Bio::EnsEMBL::Transcripts with exons and supp_features
  
  my @good_matches;


  TRANSCRIPT:
  foreach my $transcript (@$transcripts ){

    my $coverage  = $self->_get_transcript_coverage($transcript);
    my $percent_id  = $self->_get_transcript_percent_id($transcript);
    my $pep = $transcript->translate->seq;
    my $num_stops = $pep =~ tr/\*/\*/;

    if (defined $self->min_percent and $percent_id < $self->min_percent) {
      print "Transcript REJECTED with: perc_ident: $percent_id\n";
      next TRANSCRIPT;
    }
    if (defined $self->min_coverage and $coverage < $self->min_coverage){
      print "Transcript REJECTED with: coverage: $coverage\n";
      next TRANSCRIPT;
    }
    if (defined $self->max_stops and $num_stops > $self->max_stop){
      print "Rejecting proj due to too many STOPS ",$num_stops,"\n";
      next TRANSCRIPT;
    }
    
    print "Transcripts passed filter (cov=$coverage, per=$percent_id, stops=$num_stops\n";
    push @good_matches, $transcript;
  }

  my @best_matches;

  if ($self->best_in_genome == 1){

    @good_matches = sort { $a->slice->seq_region_name cmp $b->slice->seq_region_name or
                             $a->start <=> $b->start } @good_matches;

    my @c;

    foreach my $t (@good_matches) {

      if (not @c or
        $c[-1]->{name} ne $t->slice->seq_region_name or
        $c[-1]->{end} < $t->start) {

        push @c, { 
          name  => $t->slice->seq_region_name,
          start => $t->start,
          end   => $t->end,
          trans => [$t],
        };

      } else {

        push @{$c[-1]->{trans}}, $t;

        if ($t->{end} > $c[-1]->{end}) {

          $c[-1]->{end} = $t->end;

        }
      }
    }

    #####################

    my ($best_c, $best_cov);
    foreach my $c (@c) {

      my @trans = @{$c->{trans}};
      my ($total_cov, $total_pid);

      foreach my $t (@trans) {

        my ($sf) = @{$t->get_all_supporting_features};
        my $cov = $sf->score;
        my $pid = $sf->percent_id;
      
        $total_cov += $cov;
        $total_pid += $pid;
      }

      $total_cov /= scalar(@trans);
      $total_pid /= scalar(@trans);
      $c->{average_coverage} = $total_cov;
      $c->{average_percent_id} = $total_pid;
    }

    @c = sort { $b->{average_coverage} <=> $a->{average_coverage} or
                 $b->{average_percent_id} <=> $a->{average_percent_id} } @c;

    if (@c && @{$c[0]->{trans}}){
      my $best_c = shift @c;
      push @best_matches, @{$best_c->{trans}};
      foreach my $oc (@c) {
        if ($oc->{average_coverage} >= $best_c->{average_coverage} and
            $oc->{average_percent_id} >= $best_c->{average_percent_id}) {
          push @best_matches, @{$oc->{trans}};
        }
      }
    }
  }

  return \@best_matches;

}

############################################################

sub _get_transcript_coverage{
  my ($self,$tran) = @_;

  if (@{$tran->get_all_supporting_features} and
      defined $tran->get_all_supporting_features->[0]->hcoverage) {
    my ($evi) = @{$tran->get_all_supporting_features};
    return $evi->hcoverage;
  } else {
    my @exons = @{$tran->get_all_Exons};
    my ($evi) = @{$exons[0]->get_all_supporting_features};
    return $evi->score;
  }
}

############################################################

sub _get_transcript_percent_id{
  my ($self,$tran) = @_;

  my ($sf);

  if (@{$tran->get_all_supporting_features}) {
    ($sf) = @{$tran->get_all_supporting_features};
  } else {
    my @exons = @{$tran->get_all_Exons};
    ($sf) = @{$exons[0]->get_all_supporting_features};    
  }

  return $sf->percent_id;
}

############################################################

sub _get_transcript_evidence_id{
  my ($self,$tran) = @_;

  my ($sf);

  if (@{$tran->get_all_supporting_features}) {
    ($sf) = @{$tran->get_all_supporting_features};
  } else {
    my @exons = @{$tran->get_all_Exons};
    ($sf) = @{$exons[0]->get_all_supporting_features};    
  }
  
  return $sf->hseqname;
}

############################################################

# containers

sub min_coverage{
  my $self = shift;
  $self->{'_min_coverage'} = shift if(@_);

  return exists($self->{'_min_coverage'}) ? $self->{'_min_coverage'} : undef;
}

sub min_percent{
  my $self = shift;
  $self->{'_min_percent'} = shift if(@_);

  return exists($self->{'_min_percent'}) ? $self->{'_min_percent'} : undef;
}

sub max_stops{
  my $self = shift;
  $self->{'_max_stops'} = shift if(@_);

  return exists($self->{'_max_stops'}) ? $self->{'_max_stops'} : undef;
}


sub best_in_genome{
  my $self = shift;
  $self->{'_best_in_genome'} = shift if(@_);

  return exists($self->{'_best_in_genome'}) ? $self->{'_best_in_genome'} : undef;
}


1;
