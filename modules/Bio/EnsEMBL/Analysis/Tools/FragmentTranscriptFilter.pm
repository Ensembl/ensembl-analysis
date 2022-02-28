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

  Bio::EnsEMBL::Analysis::Tools::FragmentTranscriptFilter

=head1 SYNOPSIS

  my $filter = new Bio::EnsEMBL::Analysis::Tools::FragmentTranscriptFilter->
  new->(
        -coverage => 80,
        -percent_id => 90,
       );

  my @filtered_results = @{$filter->filter_results(\@results)};

=head1 DESCRIPTION

This is a best-in-genome filter is designed for mapping proteins/cDNAs/ESTs to
a low-coverage, fragmented genome, where different parts a single transcribed sequence 
may validly map to different top-level sequences in the target.

=cut


package Bio::EnsEMBL::Analysis::Tools::FragmentTranscriptFilter;

use strict;
use warnings;

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );


use vars qw (@ISA);

@ISA = qw(Bio::EnsEMBL::Root);



=head2 new

  Returntype: Bio::EnsEMBL::Analysis::Tools::FragmentTranscriptFilter
  Exceptions: none
  Example   : 

=cut



sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  &verbose('WARNING');
  my (
      $min_percent,
      $min_score,
      ) = 
        rearrange([
                   'PERCENT_ID',
                   'SCORE',
                   ], @args); 

  ######################
  #SETTING THE DEFAULTS#
  ######################

  $self->min_percent($min_percent) if defined $min_percent;
  $self->min_score($min_score) if defined $min_score;

  return $self;
}


#filter methods



=head2 filter_results

  Arg [1]   : self
  Arg [2]   : arrayref of Trancripts
  Function  : filter the given Transcruipts in the tried and trusted manner
  Returntype: arrayref
  Exceptions: throws if passed nothing or not an arrayref
  Example   : 

=cut



sub filter_results{
  my ($self, $transcripts) = @_;

  my %trans_by_hid;

  foreach my $tran (@$transcripts) {
    # transcript will only have one supporting feature for use cases of this filter
    my ($sf) = @{$tran->get_all_supporting_features};

    next if defined($self->min_score) and $sf->score < $self->min_score;
    next if defined($self->min_percent) and $sf->percent_id < $self->min_percent;

    push @{$trans_by_hid{$sf->hseqname}}, { 
      tran => $tran,
      score => $sf->score,
    };
  }

  my @good_transcripts;

  foreach my $hid (keys %trans_by_hid) {
    # sort transcripts by score
    my @trans = map { $_->{tran} } sort { $b->{score} <=> $a->{score} } @{$trans_by_hid{$hid}};
    
    my (@all_t_sfs);

    TRANSCRIPT:
    foreach my $tran (@trans) {
      my @t_sfs;

      foreach my $exon (@{$tran->get_all_Exons}) {
        my ($sf) = @{$exon->get_all_supporting_features};

        # check that this does not overlap with any of the previous overlap features
        foreach my $f (@all_t_sfs) {
          if ($sf->hstart <= $f->hend and $sf->hend >= $f->hstart) {            
            next TRANSCRIPT;
          }
        }
       
        push @t_sfs, $sf;
      }

      # if we get here, the transcript has passed;
      push @all_t_sfs, @t_sfs;
      push @good_transcripts, $tran;
    }
  }

  return \@good_transcripts;
}

# containers


sub min_percent{
  my $self = shift;
  $self->{'_min_percent'} = shift if(@_);

  return exists($self->{'_min_percent'}) ? $self->{'_min_percent'} : undef;
}

sub min_score{
  my $self = shift;
  $self->{'_min_score'} = shift if(@_);

  return exists($self->{'_min_score'}) ? $self->{'_min_score'} : undef;
}



1;
