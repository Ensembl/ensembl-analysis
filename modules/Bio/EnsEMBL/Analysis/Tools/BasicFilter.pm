# Ensembl module for Bio::EnsEMBL::Analysis::Tools::FeatureFilter
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

  Bio::EnsEMBL::Analysis::Tools::BasicFilter

=head1 SYNOPSIS

  my $filter = Bio::EnsEMBL::Analysis::Tools::BasicFilter
  ->new(
        -methods => {
                     score => 'greaterthan 500',
                     percent_id => 90,
                     }
        );

  my $filtered_results = $filter->filter_results($features);

=head1 DESCRIPTION

This module will take a hash keyed on method name with the 
value either being a straight cut off which it is then assumed
all features should have a value greater than or if the cut off
number is prefaced with either lessthan or greaterthan you can
vary what comparison is made. It is important to note that 
features passed in must beable to call the method specified in 
as the key and all cut offs must be numeric

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::Tools::BasicFilter;

use strict;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );





sub new {
  my ($class, @args) = @_;
  my $self = bless {},$class;

  my ($methods) = rearrange(['METHODS'], @args);
  $self->methods($methods);
  return $self;
}



sub methods{
  my ($self, $methods) = @_;
  if($methods){
    throw("The value passed into BasicFilter::methods ".
          "must be a hash ref not a ".$methods) 
      if(ref($methods) ne 'HASH');
    $self->{'methods'} = $methods;
  }
  return $self->{'methods'};
}



sub filter_results{
  my ($self, $features) = @_;
  if(!$features || ref($features) ne 'ARRAY'){
    throw("Must pass filter_results an arrayref not ".
          $features." BasicFilter::filter_results");
  }
  my $methods = $self->methods;
  my @filtered_features;
  FEATURE:foreach my $feature(@$features){
    foreach my $method(keys(%$methods)){
      throw("The features passed in ".$feature->[0].
            " must have the method specified ".$method)
        unless($feature->can($method));
      my $value = $methods->{$method};
      my @values = split /\s+/, $value;
      my ($checkmethod, $cutoff);
      if(@values == 1){
        $cutoff = $values[0];
        $checkmethod = 'greaterthan';
      }elsif(@values == 2){
        $cutoff = $values[1];
        $checkmethod = $values[0];
      }else{
        throw("String ".$value." from ".$method." key isn't ".
              "in the expected format\n");
      }
      next FEATURE unless($self->$checkmethod($feature, 
                                              $method, 
                                              $cutoff));
    }
    push(@filtered_features, $feature);
  }
  return \@filtered_features;
}



sub lessthan{
  my ($self, $feature, $method, $cutoff) = @_;
  return $feature if($feature->$method < $cutoff);
}

sub greaterthan{
  my ($self, $feature, $method, $cutoff) = @_;
  return $feature if($feature->$method > $cutoff);
}

1;
