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

Bio::EnsEMBL::Analysis::Tools::Algorithms::AlignFeatureCluster

=head1 SYNOPSIS


=head1 DESCRIPTION

This object holds one or more features which has been clustered according to 
comparison criteria external to this class (for instance, in the 
methods compare and _compare_AlignFeatures methods of the class AlignFeatureComparison).
Each AlignFeatureCluster object holds the IDs of the features clustered and the beginning and end coordinates
of each one (taken from the start and end coordinates of the first and last exon in the correspondig
get_all_Exons array)

=head1 CONTACT

eae@sanger.ac.uk

=cut

# Let the code begin ...

package Bio::EnsEMBL::Analysis::Tools::Algorithms::AlignFeatureCluster;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::BaseAlignFeature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning );

=head1 METHODS

=cut

#########################################################################


=head2 new

new() initializes the attributes:

$self->{'_benchmark_features'}
$self->{'_prediction_features'}

=cut

sub new {
  my ($class, $ignore_strand, $whatever) =@_ ;

  if (ref($class)){
    $class = ref($class);
  }
  my $self = {};
  bless($self,$class);

  $self->{_ignore_strand} = $ignore_strand;

  if ($whatever){
    throw( "Can't pass an object to new() method. Use put_AlignFeatures() to include Bio::EnsEMBL::AlignFeature in cluster");
  }

  $self->{_cached_start}  = undef;
  $self->{_cached_end}    = undef;
  $self->{_cached_strand} = undef;
  
  $self->{v} = 0 ; # verbosity 
  return $self;
}


=head1 Range-like methods

Methods start and end are typical for a range. We also implement the boolean
and geometrical methods for a range.

=head2 start

  Title   : start
  Usage   : $start = $feature_cluster->end();
  Function: get/set the start of the range covered by the cluster. This is re-calculated and set everytime
            a new feature is added to the cluster
  Returns : a number
  Args    : optionaly allows the start to be set

=cut

# method to get the start of the cluster, which we take to be the left_most exon_coordinate
# i.e. the start coordinate of the first exon ordered as { $a->start <=> $b->start }, regardless of the strand

sub start {
  my ($self, $start) = @_ ;

  if ($start) {
    throw( "$start is not an integer") unless $start =~/^[-+]?\d+$/;
    $self->{_cached_start} = $start;
  }

  if (!defined($self->{_cached_start})) {
    my $start;

    foreach my $feature (@{$self->get_AlignFeatures}) {
      my $this_start = $feature->start;
      unless ( $start ){
        $start = $this_start;
      }
      if ( $this_start < $start ){
        $start = $this_start;
      }
    }
    $self->{_cached_start} = $start;
  }
  return $self->{_cached_start};
}


############################################################

=head2 end

  Title   : end
  Usage   : $end = $feature_cluster->end();
  Function: get/set the end of the range covered by the cluster. This is re-calculated and set everytime
            a new feature is added to the cluster
  Returns : the end of this range
  Args    : optionaly allows the end to be set
          : using $range->end($end
=cut


# method to get the end of the cluster, which we take to be the right_most exon_coordinate
# this being the end coordinate of the first exon ordered as { $b->end <=> $a->end }, regardless of the strand

sub end {
  my ($self, $end) = @_ ;

  if ($end) {
    throw( "$end is not an integer") unless $end =~/^[-+]?\d+$/;
    $self->{_cached_end} = $end;
  }

  if (!defined($self->{_cached_end})) {
    my $end;

    foreach my $feature (@{$self->get_AlignFeatures}) {
      my $this_end = $feature->end;
      unless ( $end ){
        $end = $this_end;
      }
      if ( $this_end > $end ){
        $end = $this_end;
      }
    }
    $self->{_cached_end} = $end;
  }
  return $self->{_cached_end};
}


############################################################

=head2 length

  Title   : length
  Usage   : $length = $range->length();
  Function: get/set the length of this range
  Returns : the length of this range
  Args    : optionaly allows the length to be set
          : using $range->length($length)

=cut

sub length {
  my $self = shift @_ ;
  if (@_) {
    $self->confess( ref($self)."->length() is read-only") ;
  }
  return ( $self->{_cached_end} - $self->{_cached_start} + 1 ) ;
}


#########################################################################

=head2 strand

  Title   : strand
  Usage   : $strand = $feature->strand();
  Function: get/set the strand of the features in the cluster.
            The strand is set in put_AlignFeatures when the first feature is added to the cluster
            in this method there is also a check for strand consistency everytime a new feature is added
            Returns 0/undef when strand not set
  Returns : the strandidness (-1, 0, +1)
  Args    : optionaly allows the strand to be set

=cut


sub strand {
  my ($self, $strand) = shift;

  if ($strand) {
    if ( $self->{_cached_strand} ) {
      print "Strand called, but ignoring\n";
      return 0;
    }

    $self->{_cached_strand} = $strand ;
  }

  if (!defined($self->{_cached_strand})) {
    my @features = @{ $self->get_AlignFeatures } ;
    unless (@features) {
      $self->warning("cannot retrieve the strand in a cluster with no features");
    }
    my $strand;
    foreach my $feature (@features) {
      if (ref($feature) =~ m/AlignFeature/) {
        $feature->strand;
      }
    }
    if (!defined $strand) {
      throw("Strand not defined");
    }
    $self->{_cached_strand} = $strand ;
  }
  return $self->{_cached_strand};
}




#########################################################################

=head2 put_AlignFeatures

  function to include one or more features in the cluster.
  Useful when creating a cluster. It takes as argument an array of features, it returns nothing.

=cut

sub put_AlignFeatures {
  my ($self, $features, $ignore_strand)= @_ ;

  if ( !defined( $self->{'_types_sets'} ) ){
    throw( "Cluster lacks references to feature-types, unable to put the feature");
  } 

  unless ( ref($features) =~ m/ARRAY/ ) {   
    throw("Only take array ref !\n") ; 
  }

# Adjust cluster boundaries with new added features

  foreach my $feature (@$features) {
    if ( !defined ($self->{_cached_start}) || $feature->start < $self->start ) {
      $self->start ( $feature->start ) ;
    }
    if ( !defined ($self->{_cached_end}) || $feature->end > $self->end ) {
      $self->end ( $feature->end );
    }
  }

# Check strand consistency

  foreach my $feature (@$features) {
    if (!$ignore_strand) {
      if ( defined ($self->{_cached_strand} ) ) {
        if ( $self->strand != $feature->strand ) {
          warning( "You're trying to put $feature in a cluster of opposite strand");
        }
      }
    } else {
  # we can ignore the strand, and do nothing
    }
  }
 

 FEATURE:
  foreach my $feature (@$features) {
    throw("undef for feature. Cannot put_AlignFeatures") if (!$feature) ;
    my $feature_logicname = $feature->analysis->logic_name ;
    foreach my $set_name ( keys %{$self->{'_types_sets'}}) { 
      my $set = $self->{'_types_sets'}{$set_name} ; 
      foreach my $type ( @{$set} ) {
        if ($feature_logicname eq $type) {
          push ( @{ $self->{'_feature_sets'}{$set_name} }, $feature ) ;
          next FEATURE; 
        }
      }
    }
    throw("Failed putting feature of type " . $feature->analysis->logic_name . "\n");
  }
}



sub get_sets_included {
  my $self = shift;
  my @included_sets;

  foreach my $set_name ( keys %{$self->{'_types_sets'}}) {
    if (defined( $self->{'_feature_sets'}{$set_name})) {
      push @included_sets,$set_name;
    }
  }
  return \@included_sets;
}


#########################################################################

=head2 get_AlignFeatures

  it returns the array of features in the AlignFeatureCluster object

=cut

sub get_AlignFeatures {
  my $self = shift @_ ;

  my @features ;
  if (!defined( $self->{'_feature_sets'} ) ) {
    $self->warning("The feature array you try to retrieve is empty") ;
    @features = () ;
  }

  foreach my $set_name (keys %{$self->{'_feature_sets'}}) {
    push( @features, @{ $self->{'_feature_sets'}{$set_name} } ) ;
  }

  return \@features;
}


#########################################################################


=head2 get_AlignFeature_Count

  it returns the number of genes in the AlignFeatureCluster object

=cut

sub get_AlignFeature_Count {
  my $self = shift @_;

  my @dafs = @{$self->get_AlignFeatures} ;
  return scalar(@dafs);
}




sub feature_Types {
  my ($self, $set_name, $types) = @_;
  $self->{'_types_sets'}{$set_name} = $types;
  return $self->{'_types_sets'}{$set_name};
}

#########################################################################

sub get_AlignFeatures_by_Set() {
  my ($self,$set) = @_;

  unless ($set){
    throw( "must provide a set");
  }

  my @selected_features;
  if ($self->{v}){ 
    for (keys %{ $self->{_feature_sets} } ) { 
      print " i know the following sets : $_\n" ; 
    }
  } 
  if (!defined($self->{'_feature_sets'}{$set})) {
    # throw("No features of set name $set");
    #warning("No features of set name $set in cluster");
  }else{
    push @selected_features, @{$self->{'_feature_sets'}{$set}};
  }
  return \@selected_features;
}
#########################################################################

1;
