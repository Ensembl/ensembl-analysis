# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

=cut

package Bio::EnsEMBL::Analysis::Tools::BasicFilter;

use warnings ;
use strict;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );





=head2 new

  Arg [1]   : hashref, containing list of method names and
  cut off values
  Function  : create a Bio::EnsEMBL::Analysis::Tools::
  BasicFilter object
 Returntype: Bio::EnsEMBL::Analysis::Tools::BasicFilter
  Exceptions: none
  Example   : see docs above

=cut


sub new {
  my ($class, @args) = @_;
  my $self = bless {},$class;

  my ($methods) = rearrange(['METHODS'], @args);
  $self->methods($methods);
  return $self;
}



=head2 methods

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::BasicFilter
  Arg [2]   : hashref containing a list of method names and cut
  off values
  Function  : container method
  Returntype: hashref
  Exceptions: throw if not passed a hashref
  Example   : 

=cut


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



=head2 filter_results

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::BasicFilter
  Arg [2]   : arryref of objects
  Function  : to filter objects on specific criteria
  Returntype: arrayref of objects
  Exceptions: throws if not passed an arrayref or if objects
  cant do the specified methods
  Example   : 

=cut


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
      throw("Cut off ".$cutoff." must be a numeric value from".
            " ".$method." ")
        unless($cutoff =~ /\d+/);
      next FEATURE unless($self->$checkmethod($feature, 
                                              $method, 
                                              $cutoff));
    }
    push(@filtered_features, $feature);
  }
  return \@filtered_features;
}




=head2 lessthan

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::BasicFilter
  Arg [2]   : object
  Arg [3]   : string which is method to call on object
  Arg [4]   : cutoff object needs to have value greater or
  lesser than
  Function  : returns object if method returns value less than specific cutoff
  Returntype: object
  Exceptions: none
  Example   : 

=cut



sub lessthan{
  my ($self, $feature, $method, $cutoff) = @_;
  return $feature if($feature->$method < $cutoff);
}

=head2 greaterthan

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::BasicFilter
  Arg [2]   : object
  Arg [3]   : string which is method to call on object
  Arg [4]   : cutoff object needs to have value greater or
  lesser than
  Function  : returns object if method returns value greater
  than specific cutoff
  Returntype: object
  Exceptions: none
  Example   : 

=cut
sub greaterthan{
  my ($self, $feature, $method, $cutoff) = @_;
  return $feature if($feature->$method > $cutoff);
}

1;
