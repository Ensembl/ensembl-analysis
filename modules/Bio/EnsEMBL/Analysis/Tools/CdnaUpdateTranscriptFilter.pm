# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

  Bio::EnsEMBL::Analysis::Tools::CdnaUpdateTranscriptFilter

=head1 SYNOPSIS

  my $filter = new Bio::EnsEMBL::Analysis::Tools::CdnaUpdateTranscriptFilter
  new->(
        -best_in_genome => 1,
        -reject_processed_pseudos => 1,
        -coverage => 80,
        -percent_id => 90,
		-verbosity => 0,
       );

  my @filtered_results = @{$filter->filter_results(\@results)};

=head1 DESCRIPTION

This is a module used for filtering Exonerate transcripts

=cut


package Bio::EnsEMBL::Analysis::Tools::CdnaUpdateTranscriptFilter;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::UnmappedObject;




=head2 new

  Returntype: Bio::EnsEMBL::Analysis::Tools::CdnaUpdateTranscriptFilter
  Exceptions: none
  Example   : 

=cut



sub new {
  my ($class,@args) = @_;

  my $self = bless {},$class;
  &verbose('WARNING');

  my ($min_coverage,
      $min_percent,
      $best_in_genome,
      $rpp,
	  $verbosity) = 
        rearrange([
                   'COVERAGE',
                   'PERCENT_ID', 
                   'BEST_IN_GENOME',
                   'REJECT_PROCESSED_PSEUDOS',
				   'VERBOSITY',], @args); 

  ######################
  #SETTING THE DEFAULTS#
  ######################

  $self->min_coverage($min_coverage) if defined $min_coverage;
  $self->min_percent($min_percent) if defined $min_percent;
  $self->best_in_genome($best_in_genome) if defined $best_in_genome;
  $self->reject_processed_pseudos($rpp) if defined $rpp;
  $self->verbosity($verbosity) if defined $verbosity;
  return $self;
}


#filter methods



=head2 filter_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::DefaultExonerateFilter
  Arg [2]   : arrayref of Trancripts
  Function  : filter the given Transcripts in the tried and trusted manner
  Returntype: arrayref
  Exceptions: throws if passed nothing or not an arrayref
  Example   : 

=cut



sub filter_results {
  my ($self, $transcripts) = @_;
  # results are Bio::EnsEMBL::Transcripts with exons and supp_features
  my @good_matches;  
  my %matches;
  my $printing = $self->verbosity;
  my @rejected;
  my $analysis = $transcripts->[0]->analysis;

TRAN:
  
  foreach my $transcript (@$transcripts ){
    my $coverage  = $self->_get_transcript_coverage($transcript);
    my $percent_id  = $self->_get_transcript_percent_id($transcript);

	#identify the longest intron
	my @sorted_exons = sort {$a->start <=> $b->start} @{$transcript->get_all_Exons()};
	my $num_exons = scalar @sorted_exons;
	my $max_intron = 0;
	my $short_intron = 0; #flag to check that not all introns are long
	
	for my $a (0 .. $#sorted_exons - 1){
		my $intron_len = ($sorted_exons[$a+1]->start -  $sorted_exons[$a]->end) - 1;
		if ($intron_len > $max_intron){
			$max_intron = $intron_len;
		}
		if ($intron_len < 250000){
			$short_intron = 1;
		}
	}

    my $id = $self->_get_transcript_evidence_id($transcript);
    push @{$matches{$id}}, {
      transcript  => $transcript,
      coverage    => $coverage,
      percent_id  => $percent_id,
      num_exons   => $num_exons,
      is_spliced  => $self->_transcript_is_spliced($transcript),
	  max_intron  => $max_intron,
	  short_intron => $short_intron,
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
    
  TRANSCRIPT:
    
    foreach my $hit ( @{$matches_sorted_by_coverage{$query_id}} ){
	  my $coverage = $hit->{coverage};
      my $percent_id = $hit->{percent_id};
      my $is_spliced = $hit->{is_spliced};
	
      unless ($max_coverage){
        $max_coverage = $coverage;
      }
      unless ( $perc_id_of_best ){
		$perc_id_of_best = $percent_id;
      }
	
      #sd3
      #single exon genes (ie mouse olfactory receptors) being thrown out in favour of 
      #low percentage id multi-exon genes of equivalent coverage
      #this loop checks for high quality matches in multi-exon hits before the splices_elsewhere flag is set
      
      #so if have good coverage & percentage id:
      if ( (($coverage >= $self->min_coverage && $percent_id >= $self->min_percent) 
      
        #or have high coverage & a slightly lower percentage id
        || ($coverage   >= (1 + 5/100) * $self->min_coverage &&
            $percent_id >= (1 - 3/100) * $self->min_percent)) 
	    && $is_spliced){
        
	$splices_elsewhere = 1;
	
	last;
      }
      
    }	
    
    foreach my $hit ( @{$matches_sorted_by_coverage{$query_id}} ){
      $count++;

      my ($accept, $label);      
      my $transcript = $hit->{transcript};
      my $strand = $transcript->strand;
      my $coverage = $hit->{coverage};
      my $percent_id = $hit->{percent_id};
      my $is_spliced = $hit->{is_spliced};
	  my $max_intron = $hit->{max_intron};
	  my $short_intron = $hit->{short_intron};
	  my $num_exons = $hit->{num_exons};
     
      if ( $count == 1 ){
		$label = 'best_match';
      } elsif ( $count > 1 && 
        	$splices_elsewhere && 
            ! $is_spliced) {
		$label = 'potential_processed_pseudogene';
      } else{
		$label = $count;
      }
      
#      if ( $count == 1 && $is_spliced ){ #old way of doing it
#	     $splices_elsewhere = 1;
#      }
	
     
  if ( $self->best_in_genome ){
    # we keep the hit with the best coverage...
	  if ($coverage == $max_coverage &&
              # as long as it has coverage/percent_id above limits or...
              (($coverage >= $self->min_coverage && 
            	$percent_id >= $self->min_percent)
               ||
               # ...if coverage is significantly greater than the
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
	      if ($printing){
		      print  "rpp $query_id\n";
          if ($printing == 2) {
            push(@rejected, Bio::EnsEMBL::UnmappedObject->new(
             -type => 'cDNA',
             -identifier => $query_id,
             -summary => 'Processed pseudogene',
             -full_desc => 'Rejected as a processed pseudogene because there are multiple-exon hits with the same coverage which have been rejected for other reasons',
             -target_score => $count,
             -analysis => $analysis,
            ));
          }
		  }
		}
		elsif (($short_intron == 0) && ($num_exons > 1)){
			#all long introns
			$accept = 'NO';
			if ($printing){
			  print "only long introns $query_id\n";
        if ($printing == 2) {
          push(@rejected, Bio::EnsEMBL::UnmappedObject->new(
           -type => 'cDNA',
           -identifier => $query_id,
           -summary => 'All long introns',
           -full_desc => 'Every intron in these hits is of length 250000-400000bp, we require at least one intron to be shorter than 250000bp',
           -target_score => $count,
           -analysis => $analysis,
          )) if ($count == 1);
        }
			}
		}
		#Only accept long intron hits with very high coverage and percent_id
		elsif($max_intron > 250000 ){ 
			if (($coverage >= 98) && ($percent_id >= 98)){
				$accept = 'YES';
				push( @good_matches, $transcript);
				
				#print "accept: intron $max_intron coverage $coverage \%id $percent_id $query_id\n"; 
				
				#find out which introns are long - ie are they the first and last introns?
#				my @exons = @{$transcript->get_all_Exons()};
#				my @sorted_exons = sort {$a->start <=> $b->start} @exons;
#				my $num_introns = scalar @exons - 1;
#
#				for my $a (0 .. $#sorted_exons - 1){
#					my $intron_len = ($sorted_exons[$a+1]->start -  $sorted_exons[$a]->end) - 1;
#					if ($intron_len > 250000){
#						print "intron position ".($a + 1)." out of $num_introns\n"; 		
#					}
#				}
				
			}else{
				$accept = 'NO';
				if ($printing){
				  print "reject: intron $max_intron coverage $coverage \%id $percent_id $query_id\n";
          if ($printing == 2) {
            my $summary = 'Low coverage with long intron';
            my $full_desc = 'Hits containing introns longer than 250000bp are rejected if coverage is less than 98% - see query_score for coverage';
            my $score = $coverage;
            if ($percent_id < 98) {
              $summary = 'Low percent_id with long intron';
              $full_desc = 'Hits containing introns longer than 250000bp are rejected if percentage identity is less than 98% - see query_score for percent_id';
              $score = $percent_id;
            }
            push(@rejected, Bio::EnsEMBL::UnmappedObject->new(
             -type => 'cDNA',
             -identifier => $query_id,
             -summary => $summary,
             -full_desc => $full_desc,
             -target_score => $count,
             -analysis => $analysis,
             -query_score => $score,
            )) if ($count == 1);
          }
				}
			}	
		}else{
	      $accept = 'YES';
	      push( @good_matches, $transcript);

		}
	  }
	  else{
		$accept = 'NO';
 		if ($printing){
		  print "max_coverage $max_coverage coverage $coverage \%id $percent_id $query_id\n"; 
      if ($printing == 2) {
        my $summary = 'Low coverage';
        my $full_desc = 'Coverage of the best alignment is less than 90% - see query_score for coverage';
        my $score = $coverage;
        if ($percent_id < $self->min_percent) {
          $summary = 'Low percent_id';
          $full_desc = 'Percentage identity of the best alignment is less than 97% - see query_score for percent_id';
          $score = $percent_id;
        }
        push(@rejected, Bio::EnsEMBL::UnmappedObject->new(
         -type => 'cDNA',
         -identifier => $query_id,
         -summary => $summary,
         -full_desc => $full_desc,
         -target_score => $count,
         -query_score => $score,
         -analysis => $analysis,
        )) if ($count == 1);
      }
		}

	  }

  }
  else{
        # we keep anything which is within the 2% of the best score...
	if ($coverage >= (0.98 * $max_coverage) && 
            # as long as it has coverage/percent_id above limits or...
            (($coverage >= $self->min_coverage && 
              $percent_id >= $self->min_percent)
             ||              
             # ...if coverage is significantly greater than the
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
  return(\@good_matches, \@rejected);

}

############################################################

sub _get_transcript_coverage {
  my ($self,$tran) = @_;
  if (@{$tran->get_all_supporting_features} and
      defined $tran->get_all_supporting_features->[0]->hcoverage) {
    return $tran->get_all_supporting_features->[0]->hcoverage;
  } else {
    return $tran->get_all_Exons->[0]->get_all_supporting_features->score;
  }
}

############################################################

sub _get_transcript_percent_id {
  my ($self,$tran) = @_;

  my $sf;

  if (@{$tran->get_all_supporting_features}) {
    $sf = $tran->get_all_supporting_features->[0];
  } else {
    $sf = $tran->get_all_Exons->[0]->get_all_supporting_features->[0];
  }

  return $sf->percent_id;
}

############################################################

sub _get_transcript_evidence_id {
  my ($self,$tran) = @_;

  my $sf;

  if (@{$tran->get_all_supporting_features}) {
    $sf = $tran->get_all_supporting_features->[0];
  } else {
    $sf = $tran->get_all_Exons->[0]->get_all_supporting_features->[0];
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

sub min_coverage {
  my $self = shift;
  $self->{'_min_coverage'} = shift if(@_);

  return exists($self->{'_min_coverage'}) ? $self->{'_min_coverage'} : undef;
}

sub min_percent {
  my $self = shift;
  $self->{'_min_percent'} = shift if(@_);

  return exists($self->{'_min_percent'}) ? $self->{'_min_percent'} : undef;
}


sub best_in_genome {
  my $self = shift;
  $self->{'_best_in_genome'} = shift if(@_);

  return exists($self->{'_best_in_genome'}) ? $self->{'_best_in_genome'} : 0;
}

sub reject_processed_pseudos {
  my $self = shift;
  $self->{'_reject_processed_pseudos'} = shift if(@_);

  return exists($self->{'_reject_processed_pseudos'}) ? $self->{'_reject_processed_pseudos'} : 0;
}

sub verbosity {
  my $self = shift;
  $self->{'_verbosity'} = shift if(@_);

  return exists($self->{'_verbosity'}) ? $self->{'_verbosity'} : 0;
}


1;
