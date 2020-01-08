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

Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils - base class for
genebuild utility methods

=head1 SYNOPSIS

  use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(coord_string id 
                                                       empty_Object);

  or 

  use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils 

  to get all methods

=head1 DESCRIPTION

This is a base class for Utility modules for genebuilding and
other gene manupulation code. This module probably wont be used
direct outside of these classes but it can been if needed

It provides some simple functionality that all of the utility
modules need like id, to provide a sensible id string or
coord string to provide a basic, start, end, strand, 
seq_region_name string for printing or removing databases connections
from objects

These modules are heavily based on the Utils modules which
can be found in EnsEMBL::Analysis::Pipeline::Tools

=head1 CONTACT

please send any questions to http://lists.ensembl.org/mailman/listinfo/dev

=head1 METHODS

the rest of the documention details the exported static
class methods

=cut


package Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils;

use strict;
use warnings;

use Exporter qw(import);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

our @EXPORT_OK = qw(
                    coord_string
                    seq_region_coord_string
                    id
                    empty_Object
                    hashkey_Object
                    lies_inside_of_slice
                  );


=head2 coord_string

  Arg [1]   : Bio::EnsEMBL::Feature
  Function  : Returns a string with the start, end, strand
              and slice name that the feature is on, delimited by spaces.
              Coord positions printed will be relative to the start of
              the slice
  Returntype: string
  Exceptions: throws if no feature is passed in
  Example   : my $coord_string = coord_string($transcript);

=cut



sub coord_string{
  my $feature = shift;
  my ($p, $f, $l) = caller;
  throw("Must be passed a feature") if(!$feature);
  my $string = $feature->start."\t".$feature->end."\t".$feature->strand."\t".$feature->slice->seq_region_name;
  return $string;
}


=head2 seq_region_coord_string

  Arg [1]   : Bio::EnsEMBL::Feature
  Function  : Returns a string with the seq_region_start, seq_region_end, strand
              and slice name that the feature is on, delimited by spaces
  Returntype: string
  Exceptions: throws if no feature is passed in
  Example   : my $coord_string = seq_region_coord_string($transcript);

=cut



sub seq_region_coord_string{
  my $feature = shift;
  my ($p, $f, $l) = caller;
  throw("Must be passed a feature") if(!$feature);
  my $string = $feature->seq_region_start."\t".$feature->seq_region_end."\t".$feature->strand."\t".$feature->slice->seq_region_name;
  return $string;
}


=head2 id

  Arg [1]   : Bio::EnsEMBL::Feature
  Function  : Returns a string containing an appropriate label
              for the feature
  Returntype: string
  Exceptions: none
  Example   : 

=cut



sub id {
  my $feature = shift;
  my $id;

  if($feature->can('stable_id') && $feature->stable_id){
    $id = $feature->stable_id;
  }elsif($feature->can('dbID') && $feature->dbID) {
    $id = $feature->dbID;
  }else{ 
    $id = 'no-id';
  }
  if($feature->can('biotype') && $feature->biotype){
    $id .= "_".$feature->biotype;
  }
  return $id;
}



=head2 empty_Object

  Arg [1]   : Bio::EnsEMBL::Storeable or an object which inherits from it
  Arg [2]   : Boolean, whether to remove the stable id from the given object
  Function  : remove the dbID, adaptor and if appropriate the stable id
  Returntype: Bio::EnsEMBL::Storeable
  Exceptions: n/a
  Example   : empty_Object($object);

=cut



sub empty_Object{
  my ($object, $include_stable_id, $new_analysis) = @_;
  $object->adaptor(undef);
  $object->dbID(undef);
  $object->stable_id(undef) if($include_stable_id &&
                               $object->can("stable_id"));
  if ($new_analysis and $object->can('analysis')) {
      $object->analysis($new_analysis);
  }
  return $object;
}



=head2 lies_inside_of_slice

  Arg [1]   : Bio::EnsEMBL::Feature
  Arg [2]   : Bio::EnsEMBL::Slice
  Function  : Ensures the transcript within the slice, completely on
              the lower end, it can overhang the upper end.
  Returntype: Boolean, 1 for pass, 0 for fail, i.e. lies outside of slice
              or across lower boundary
  Exceptions: none
  Example   :

=cut


sub lies_inside_of_slice{
  my ($feature, $slice) = @_;
  if($feature->start > $slice->length || 
     $feature->end < 1){
    warning(id($feature)." lies off edge if slice ".
            $slice->name);
    return 0;
  }
  if($feature->start < 1 && $feature->end > 1){
    warning(id($feature)." lies over lower boundary".
            " of slice ".$slice->name);
    return 0;
  }
  return 1;
}


=head2 hashkey_Object

 Arg [1]    : ArrayRef of features
 Description: Creates a string which would be unique to a certain point
              by concatenating start, end, strand of the objects with ':'
              as a separator. This is mostly used for exons and introns
 Returntype : String
 Exceptions : None

=cut

sub hashkey_Object {
  my ($objects) = @_;

  my $hashkey = '';
  foreach my $elm (@$objects) {
    $hashkey .= $elm->start.':'.$elm->end.':'.$elm->strand.':';
  }
  return substr($hashkey, 0, -1);
}

1;
