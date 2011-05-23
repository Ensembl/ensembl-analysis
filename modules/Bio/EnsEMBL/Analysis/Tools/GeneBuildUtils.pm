=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

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

=head1 METHODS

the rest of the documention details the exported static
class methods

=cut


package Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils;

use strict;
use warnings;
use Exporter;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning
                                      stack_trace_dump);
use vars qw (@ISA  @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(coord_string seq_region_coord_string id empty_Object lies_inside_of_slice);


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
  my ($object, $include_stable_id) = @_;
  $object->adaptor(undef);
  $object->dbID(undef);
  $object->stable_id(undef) if($object->can("stable_id") && 
                               $include_stable_id);
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


1;
