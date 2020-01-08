# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

  Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter

=head1 SYNOPSIS

  my $filter = new Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter
  new->(
        -best_in_genome => 1,
        -reject_processed_pseudos => 1,
        -coverage => 80,
        -percent_id => 90,
       );

  my @filtered_results = @{$filter->filter_results(\@results)};

=head1 DESCRIPTION

This is the standard module used for filtering Exonerate transcripts

=cut


package Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use Data::Dumper;

use vars qw (@ISA);



=head2 new

  Returntype: Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter
  Exceptions: none
  Example   : 

=cut



sub new{
  my ($caller,@args) = @_;

  my $class = ref($caller) || $caller;
  my $self = bless({},$class);

  &verbose('WARNING');
  my ($min_coverage,
      $min_percent,
      $best_in_genome,
      $rpp) = 
        rearrange([
                   'COVERAGE',
                   'PERCENT_ID', 
                   'BEST_IN_GENOME',
                   'REJECT_PROCESSED_PSEUDOS',], @args); 

  ######################
  #SETTING THE DEFAULTS#
  ######################

  if (defined ($min_coverage)) {
    $self->min_coverage($min_coverage);
  } else {
    warn("\n\tmin_coverage not set, setting it to zero (0)!\n\n$!");
    $self->min_coverage(0);
  }
  if (defined ($min_percent)) {
    $self->min_percent($min_percent);
  } else {
    warn("\n\tmin_percent not set, setting it to zero (0)!\n\n$!");
    $self->min_percent(0);
  }
  if (defined ($best_in_genome)) {
    $self->best_in_genome($best_in_genome);
  } else {
    warn("\n\tbest_in_genome not set, setting it to one (1)!\n\n$!");
    $self->best_in_genome(1);
  }
  if (defined ($rpp)) {
    $self->reject_processed_pseudos($rpp);
  } else {
    warn("\n\treject_processed_pseudos, setting it to one (1)!\n\n$!");
    $self->reject_processed_pseudos(1);
  }

  return $self;
}


#filter methods



=head2 filter_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::DefaultExonerateFilter
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

  my %matches;

TRAN:
  foreach my $transcript (@$transcripts ){
    my $coverage  = $self->_get_transcript_coverage($transcript);
    my $percent_id  = $self->_get_transcript_percent_id($transcript);


    my $id = $self->_get_transcript_evidence_id($transcript);
    push @{$matches{$id}}, {
      transcript  => $transcript,
      coverage    => $coverage,
      percent_id  => $percent_id,
      num_exons   => scalar(@{$transcript->get_all_Exons}),
      is_spliced  => $self->_transcript_is_spliced($transcript),
    };
  }
  
  my %matches_sorted_by_coverage;
  my %selected_matches;
 
 QUERY:
  foreach my $query_id ( keys( %matches ) ){
    @{$matches_sorted_by_coverage{$query_id}} = 
        sort { $b->{coverage}   <=> $a->{coverage} or
               $b->{num_exons}  <=> $a->{num_exons} or
               $b->{percent_id} <=> $a->{percent_id} } @{$matches{$query_id}};

    my $max_coverage;
    my $perc_id_of_best;
    my $count = 0;
    my $splices_elsewhere = 0;
    my $best_has_been_seen = 0;
    
    #print STDERR "####################\n";
    #print STDERR "Matches for $query_id:\n";

  #get the slice_name, start, and end of the best hit and store them to be used
  #in checking of any other good hit may be overlaping the best prediction
  # This is to avoid a match that spans over a region where two or more proteins of
  # a same family are close together and they get merged by a wrong "good quallity" alignment
    my $best_transcript = ${$matches_sorted_by_coverage{$query_id}}[0]->{transcript};
    my $best_start = $best_transcript->start;
    my $best_end = $best_transcript->end;
    my $best_slice;
    if (!defined $best_transcript->slice) {
      $best_slice = $best_transcript->start_Exon->seqname;
    } else {
      $best_slice = $best_transcript->slice;
    }
    
  TRANSCRIPT:
    foreach my $hit ( @{$matches_sorted_by_coverage{$query_id}} ){
      $count++;

      my ($accept, $label);      
      my $transcript = $hit->{transcript};
      my $strand = $transcript->strand;
      my $coverage = $hit->{coverage};
      my $percent_id = $hit->{percent_id};
      my $is_spliced = $hit->{is_spliced};

      my $transcript_slice;
      if (!defined $transcript->slice) {
        $transcript_slice = $transcript->start_Exon->seqname;
      } else {
        $transcript_slice = $transcript->slice; 
      }

      unless ($max_coverage){
        $max_coverage = $coverage;
      }
      unless ( $perc_id_of_best ){
	$perc_id_of_best = $percent_id;
      }

      if ( $count == 1 ){
	$label = 'best_match';
      } elsif ( $count > 1 && 
                $splices_elsewhere && 
                ! $is_spliced) {
	$label = 'potential_processed_pseudogene';
      } else{
	$label = $count;
      }

      if ( $count == 1 && $is_spliced ){
	$splices_elsewhere = 1;
      }

      if ( $self->best_in_genome ){
        # we keep the hit with the best coverage...
	if ($coverage == $max_coverage &&
            # as long as it has coverage/percent_id above limits or...
            (($coverage >= $self->min_coverage && 
              $percent_id >= $self->min_percent)
             ||
             # ...if coverage is significanly greater than the
             # specified minimum, then we are willing to accept
             # hits that have a percent_id just below the specified
             # minimum
             ($coverage   >= (1 + 5/100) * $self->min_coverage &&
              $percent_id >= (1 - 3/100) * $self->min_percent))) { 
	  if ( $self->reject_processed_pseudos
	       && $count > 1 
	       && $splices_elsewhere 
	       && ! $is_spliced) {
	    $accept = 'NO';
	  }
          # ... if one transcript with lower quality completely overlaps
          # the best one don't accept the lower quality one.
          elsif ($best_slice eq $transcript_slice &&
                 $best_start > $transcript->start &&
                 $best_end   < $transcript->end){
            $accept = 'NO';
          }
	  else {
	    $accept = 'YES';
	    push( @good_matches, $transcript);
	  }
	}
	else{
	  $accept = 'NO';
	}
	
      }
      else{
        # we keep anything which is within the 2% of the best score...
	if ($coverage >= (0.98 * $max_coverage) && 
            # as long as it has coverage/percent_id above limits or...
            (($coverage >= $self->min_coverage && 
              $percent_id >= $self->min_percent)
             ||              
             # ...if coverage is significanly greater than the
             # specified minimum, then we are willing to accept
             # hits that have a percent_id just below the specified
             # minimum
             ($coverage   >= (1 + 5/100) * $self->min_coverage &&
              $percent_id >= (1 - 3/100) * $self->min_percent))) {
          
          
	  ############################################################
	  # non-best matches are kept only if they are not unspliced with the
	  # best match being spliced - otherwise they could be processed pseudogenes
	  if ( $self->reject_processed_pseudos &&
	       $count > 1 &&
	       $splices_elsewhere &&
	       ! $is_spliced) {
	    $accept = 'NO';
	  }
	  else{
	    $accept = 'YES';
	    push( @good_matches, $transcript);
	  }
	}
	else{
	  $accept = 'NO';
	}
      }
    }
  }
  
  return \@good_matches;

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

sub _transcript_is_spliced {
  my ($self, $tran) = @_;

  my @exons = sort { $a->start <=> $b->start } @{$tran->get_all_Exons};

  if ( scalar (@exons) > 1 ){    
    # check that there are non "frameshift" introns
    for(my $i=0; $i < @exons - 1; $i++){
      my $intron_len = $exons[$i+1]->start - $exons[$i]->end - 1;
      if ( $intron_len > 9 ){
        return 1;
      }
    }
  }

  return 0;
}


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


sub best_in_genome{
  my $self = shift;
  $self->{'_best_in_genome'} = shift if(@_);

  return exists($self->{'_best_in_genome'}) ? $self->{'_best_in_genome'} : 0;
}

sub reject_processed_pseudos {
  my $self = shift;
  $self->{'_reject_processed_pseudos'} = shift if(@_);

  return exists($self->{'_reject_processed_pseudos'}) ? $self->{'_reject_processed_pseudos'} : 0;
}



1;
