package Bio::EnsEMBL::Analysis::Tools::FeatureFilter;

use strict;
use warnings;

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use vars qw (@ISA);

@ISA = qw(Bio::EnsEMBL::Root);


sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  &verbose('WARNING');
  my ($features, $min_score, $max_pvalue, $coverage,
      $prune, $hard_prune, $filter_on_coverage) = 
        rearrange(['FEATURES', 'MIN_SCORE', 
                   'MAX_PVALUE', 'COVERAGE',
                   'PRUNE', 'HARD_PRUNE', 
                   'FILTER_ON_COVERAGE'], @args);
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
  $self->features($features);
  $self->min_score($min_score) if(defined $min_score);
  $self->max_pvalue($max_pvalue) if(defined $max_pvalue);
  $self->coverage($coverage) if(defined $coverage);
  $self->prune($prune) if(defined $prune);
  $self->hard_prune($hard_prune) if(defined $hard_prune);
  $self->filter_on_coverage($filter_on_coverage) 
    if(defined($filter_on_coverage));
  #print "Have min score of ".$self->min_score."\n";
  return $self;
}


#containers

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
sub features{
  my ($self, $features) = @_;
  if(!$self->{'features'}){
    $self->{'features'} = [];
  }
  if($features){
    throw("Must pass features an arrayref not ".$features.
          "FeatureFilter:features") unless(ref($features) eq 'ARRAY');
    if(@{$self->{'features'}} == 0){
      $self->{'features'} = $features;
    }else{
      push(@{$self->{'features'}}, $features);
    }
  }
  return $self->{'features'};
}

sub clear_features{
  my ($self) = @_;
  $self->{'features'} = [];
}

#filter methods


sub filter_results{
  my ($self, $features) = @_;
  if(!$features){
    $features = $self->features;
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
  @$features = sort { $b->score <=> $a->score } @$features;
  #print "Have ".@$features. " features to filter\n";
  foreach my $f(@$features){
    #print "Score ".$f->score." compared to min ".$self->min_score."\n";
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
  #print STDERR "After filtering on score have ".$total_count." features ".
  #  " on ".keys(%hitarray)." hit ids\n";
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
              #print "bases covered ".$list[$i]." max coverage ".
              #  $self->coverage."\n" if($verbose);
              if( $list[$i] < $self->coverage ) {
                # accept!
                #print "Accepting ".$name."\n" if($verbose);
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
    $total_count = 0;
    foreach my $id(keys(%accepted_hit_ids)){
      $total_count += @{$hitarray{$id}};
    }
  }else{
    foreach my $name(keys(%hitarray)){
      $accepted_hit_ids{$name} = 1;
    }
  }
  my @features;
  if($self->prune){
    foreach my $name(keys(%accepted_hit_ids)){
      my @tmp = @{$self->prune_features(\@{$hitarray{$name}})};
      push(@features, @tmp);
    }
  }else{
    foreach my $name(keys(%hitarray)){
     push(@features, @{$hitarray{$name}});
    }
  }
  if($self->hard_prune){
    @features = @{$self->prune_features(\@features)};
  }
  return \@features;
}


sub prune_features{
  my ($self, $features) = @_;
  my @output;
  my $forward = $self->prune_features_by_strand(1, $features);
  push(@output, @$forward) if($forward);
  my $reverse = $self->prune_features_by_strand(-1, $features);
  push(@output, @$reverse) if($reverse);
  return \@output;
}


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

   @input_for_strand = ();	# free some memory?

   # over_covered_bases: list of base numbers where coverage must be
   # reduced, listed worst-case-first
   my $max_coverage = $self->coverage;
   my @over_covered_bases;
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
