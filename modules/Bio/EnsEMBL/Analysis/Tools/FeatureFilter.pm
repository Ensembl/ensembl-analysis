# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

  Bio::EnsEMBL::Analysis::Tools::FeatureFilter

=head1 SYNOPSIS

  my $featurefilter = new Bio::EnsEMBL::Analysis::Tools::FeatureFilter
  new->(
        -min_score => 200,
        -max_pvalue => 0.01,
        -coverage => 5,
        -prune => 1,
       );

  my @filtered_results = @{$featurefilter->filter_results(\@results)};

=head1 DESCRIPTION

This is the standard module used for filtering blast results
It filters on 4 separate things, It will throw out features whose score
is below a certain value or whose pvalue is above a certain value. Coverage
keeps all hit ids where at least of its features base pairs is covered by 
less than its value in features (10 as standard). Prune and hard prune 
function in the same manner but prune only considers sets of features with
the same hit id but hard prune considers all of them. The feature set
which is passed to prune_features is first split on the basis of strand
then each set of features on one strand are considered and the number of
features which cover each base pair in the range of basepairs covered is
counted. Then starting with the lowest scoring features and the highest
covered base pairs first a single feature is thrown out if at least one
of its base pairs covered the query sequence too much (this level is
also set by coverage)

here is an asci drawing which should hopefully explain. In this system
there are no min scores and coverage is set to 5 but nothing is thrown
away by the initial coverage filter as every feature has every feature has
at least on base pair whose coverage is below 5. The features with bits 
scores of 50 and 100 though would be thrown away by prune as base 21 is
covered by 7 features and these are the lowest scoring of those 7 features
if normal prune was used these features would all need to share the same 
hit id but if hard prune was used that would not be true

                                 bit score
                ---------------- 50
 ---------------------           100
               ---------         250
           ---------------       500
    -------------------          750
                    -----------  1000
 --------------------------      1250
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
1   5    1    1    2    2    3
         0    5    0    5    0

=head1 CONTACT

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

=cut


package Bio::EnsEMBL::Analysis::Tools::FeatureFilter;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use vars qw (@ISA);

@ISA = qw();



=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::FeatureFilter
  Arg [2]   : int, minimum score
  Arg [3]   : int, maximum p value
  Arg [4]   : int, maximun coverage
  Arg [5]   : int, toggle whether to prune features
  Arg [6]   : int, toggle whether to use hard prune
  Arg [7]   : int, toggle whether to filter on coverage or not
  Function  : create a new FeatureFilter object
  Returntype: Bio::EnsEMBL::Analysis::Tools::FeatureFilter
  Exceptions: none
  Example   : 

=cut



sub new{
  my ($class,@args) = @_;
  my $self = bless {},$class;
  &verbose('WARNING');
  my ($min_score, $max_pvalue, $coverage,
      $prune, $hard_prune, $filter_on_coverage) = 
        rearrange(['MIN_SCORE', 'MAX_PVALUE', 'COVERAGE',
                   'PRUNE', 'HARD_PRUNE', 'FILTER_ON_COVERAGE'], @args);
  ######################
  #SETTING THE DEFAULTS#
  ######################
  $self->min_score(-100000);
  $self->max_pvalue(0.1);
  $self->coverage(10);
  $self->prune(0);
  $self->hard_prune(0);
  $self->filter_on_coverage(1);
  ######################
  $self->min_score($min_score) if(defined $min_score);
  $self->max_pvalue($max_pvalue) if(defined $max_pvalue);
  $self->coverage($coverage) if(defined $coverage);
  $self->prune($prune) if(defined $prune);
  $self->hard_prune($hard_prune) if(defined $hard_prune);
  $self->filter_on_coverage($filter_on_coverage) 
    if(defined($filter_on_coverage));
  return $self;
}


#containers

=head2 min_score 

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::FeatureFilter
  Arg [2]   : variable, generally int or string
  Function  : This describes the 6 container methods below
  min_score, max_pvalue, coverage, prune, hard_prune and 
  filter on coverage. The all take, store and return their give
  variable
  Returntype: int/string 
  Exceptions: none
  Example   : none

=cut


sub min_score{
  my $self = shift;
  $self->{'min_score'} = shift if(@_);
  return $self->{'min_score'};
}

sub max_pvalue{
  my $self = shift;
  $self->{'max_pvalue'} = shift if(@_);
  return $self->{'max_pvalue'};
}

sub coverage{
  my $self = shift;
  $self->{'coverage'} = shift if(@_);
  return $self->{'coverage'};
}

sub prune{
  my $self = shift;
  $self->{'prune'} = shift if(@_);
  return $self->{'prune'};
}

sub hard_prune{
  my $self = shift;
  $self->{'hard_prune'} = shift if(@_);
  return $self->{'hard_prune'};
}

sub filter_on_coverage{
  my $self = shift;
  $self->{'filter_on_coverage'} = shift if(@_);
  return $self->{'filter_on_coverage'};
}



#filter methods



=head2 filter_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::FeatureFilter
  Arg [2]   : arrayref of features, there is not type checking but
  it is expected that the features have start, end, score and hseqname
  methods (Bio::EnsEMBL::FeaturePairs is what it was written against)
  Function  : filter the given features by score, pvalue coverage
  Returntype: arrayref
  Exceptions: throws if passed nothing or not an arrayref
  Example   : 

=cut



sub filter_results{
  my ($self, $features) = @_;
  if(!$features || ref($features) ne 'ARRAY'){
    throw("Must pass filter_results an arrayref not ".$features.
          " FeatureFilter::filter_results");
  }
  my %validhit;
  my %hitarray;
  my %totalscore;
  my $maxend = 0;
 
  #filtering by score
  #sorting by score so we use the highest scoring feature first
  #The score filter basically takes all features belonging to
  #one hit id provided that at least one of its features has a score
  #greater than the min score and a pvalue less than the max pvalue
  #print "Passed ".@$features." results\n";
  @$features = sort { $b->score <=> $a->score } @$features;
  
  foreach my $f(@$features){
    if($f->score > $self->min_score){ 
      if(!exists $validhit{$f->hseqname}){
        $validhit{$f->hseqname} = 0;
        $totalscore{$f->hseqname} = 0;
      }
      $totalscore{$f->hseqname} += $f->score;
      if(!$f->can('p_value') || !defined $f->p_value 
         || $f->p_value < $self->max_pvalue){
        if($validhit{$f->hseqname} < $f->score){
          $validhit{$f->hseqname} = $f->score;
        }
        if($f->end > $maxend){
          $maxend = $f->end;
        }
      }
    }
    #if a hit doesn't pass the min score or max pvalue threshold
    #but the hit id has features which do this ensures all features
    #for that hit id are kept regardless of score and pvalue
    if($validhit{$f->hseqname}){
      if(!exists($hitarray{$f->hseqname})){
        $hitarray{$f->hseqname} = [];
      }
      push(@{$hitarray{$f->hseqname}}, $f);
    }
  }
  
  my $total_count = 0;
  foreach my $id(keys(%hitarray)){
    $total_count += @{$hitarray{$id}};
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
  
  if($self->filter_on_coverage){
    my @strands = (1, -1);
    foreach my $strand(@strands){
      my @list;
      $list[$maxend] = 0; #perl will automatically extend this array
    NAME:foreach my $name(@hit_ids){

        my $hole = 0;
      FEATURE:foreach my $f(@{$hitarray{$name}}){
          next FEATURE if $f->strand != $strand;
          if( (($f->score > $self->min_score) && ($f->can('p_value'))
             && defined $f->p_value && $f->p_value < $self->max_pvalue) ){
          INDEX:foreach my $i ( $f->start .. $f->end ) {
              unless( $list[$i] ){
                $list[$i] = 0;
              }
              if( $list[$i] < $self->coverage ) {
                # accept!
                $hole = 1;
                last;
              }
            }
          }
        }
        if($hole == 0){
          next ;
        }
        $accepted_hit_ids{$name} = 1;
        foreach my $f ( @{$hitarray{$name}} ) {
          if ($f->strand == $strand) {
            for my $i ( $f->start .. $f->end ) {
              $list[$i]++; 
            }
          }
        }
      }
    }
  }else{
    foreach my $name(keys(%hitarray)){
      $accepted_hit_ids{$name} = 1;
    }
  }
  $total_count=0;
  foreach my $id(keys(%accepted_hit_ids)){
    $total_count += @{$hitarray{$id}};
  }
 
  my @features;
 
  foreach my $name(keys(%accepted_hit_ids)){
    my @tmp = @{$hitarray{$name}};
   
    my $remaining;
    if($self->prune){
      $remaining = $self->prune_features(\@tmp);
    }else{
      $remaining = \@tmp;
    }
    
    push(@features, @$remaining);
  }
  
  if($self->hard_prune){ #this pruning is done across all features
    @features = @{$self->prune_features(\@features)};
  }
  
  return \@features;
}



=head2 prune_features

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::FeatureFilter
  Arg [2]   : arrayref of features
  Function  : calls to prune_features_by_strand with each strand
  Returntype: arrayref of features
  Exceptions: none
  Example   : 

=cut



sub prune_features{
  my ($self, $features) = @_;
  my @output;
  my $forward = $self->prune_features_by_strand(1, $features);
  push(@output, @$forward) if($forward);
  my $reverse = $self->prune_features_by_strand(-1, $features);
  push(@output, @$reverse) if($reverse);
  return \@output;
}




=head2 prune_features_by_strand

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::FeatureFilter
  Arg [2]   : int, strand (must be 1 or -1)
  Arg [3]   : arrayref of features
  Function  : works out coverage of each base pair and throws out
  features of a base pair is covered by to many hits. The lowest scoring
  features are thrown out first
  Returntype: arrayref of features
  Exceptions: none
  Example   : 

=cut


#How prune_features_by_strand  works

#1. first create an array of all features of the appropriate strand

#2. Then create a separate array whose last index is the last base pair
#covered by any features

#3. go through that array and increment each index for each feature which
#covers that base pair

#4.Then it cycles through the bases covered array collecting a list of
#base pairs whose coverage is above the limit

#5. the over covered list is sorted so the bases which have the greatest
#excess coverage are considered first

#6. for each over covered base pair all all of the features are looked at
#each time a feature covers this base pair it is thrown away and the
#coverage of each base the rejected feature covers is decremented.

#7. Each over covered base is considered untill the coverage is decreased 
#to acceptable limits 

#8. the remaining feature set is returned





sub prune_features_by_strand {
   my ($self, $strand, $in) = @_;
   
   my @input_for_strand;
   foreach my $f (@$in) {
     push @input_for_strand, $f if $f->strand eq $strand;
   }
   
   return () if !@input_for_strand;
  
   my @sorted_fs = sort{ $a->start <=> $b->start } @input_for_strand;
   my $first_base = $sorted_fs[0]->start;
   @sorted_fs = sort{ $a->end <=> $b->end } @input_for_strand;
   my $last_base = $sorted_fs[$#sorted_fs]->end;

   # fs_per_base: set element i to the number of features covering base i
   my @fs_per_base;
   foreach  my $base ($first_base..$last_base) {
     $fs_per_base[$base] = 0;	# initialise
   }
   foreach my $f (@input_for_strand) {
     foreach my $covered_base ($f->start..$f->end) {
       $fs_per_base[$covered_base]++;
     }
   }

   # put the worst features first, so they get removed with priority
   @sorted_fs = sort { $a->score <=> $b->score } @input_for_strand;
   #note if you have two features with the same score you may not always
   #throw away the same one the code below eliminates this problem but it
   #isn't used as standard
   #@sorted_fs = sort { $a->score <=> $b->score || 
   #                      $a->start <=> $b->start } @input_for_strand;
 
   
   @input_for_strand = ();	# free some memory?

   # over_covered_bases: list of base numbers where coverage must be
   # reduced, listed worst-case-first
   my $max_coverage = $self->coverage;
   my @over_covered_bases;
   #print "Order of features in prune\n";
   foreach my $base ($first_base..$last_base) {
     my $excess_fs = $fs_per_base[$base] - $max_coverage;
     if ($excess_fs > 0) {
       push @over_covered_bases, $base;
     }
   }
   @over_covered_bases = sort { $fs_per_base[$b] <=> $fs_per_base[$a] }
     @over_covered_bases;

   foreach my $base (@over_covered_bases) {
     my $f_no = 0;
     while ($fs_per_base[$base] > $max_coverage) {
       my $start = $sorted_fs[$f_no]->start;
       my $end = $sorted_fs[$f_no]->end;
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
