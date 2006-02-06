
=head1 NAME

Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils - base class for
genebuild utility methods

=head1 SYNOPSIS

  use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(coord_string id);

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
seq_region_name string for printing

These modules are heavily based on the Utils modules which
can be found in EnsEMBL::Analysis::Pipeline::Tools

=head1 CONTACT

please send any questions to ensembl-dev@ebi.ac.uk

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
@EXPORT = qw(coord_string id create_filename);


=head2 coord_string

  Arg [1]   : Bio::EnsEMBL::Feature
  Function  : returns a string with the start, end, strand
  and slice name that the feature is on delimited by spaces
  Returntype: string
  Exceptions: throws if no feature is passed in
  Example   : 

=cut



sub coord_string{
  my $feature = shift;
  throw("Must be passed a feature") if(!$feature);
  my $string = "start ".$feature->start." end ".$feature->end." strand ".$feature->strand." seq region name ".$feature->slice->seq_region_name;
  return $string;
}




=head2 id

  Arg [1]   : Bio::EnsEMBL::Feature
  Function  : returns a string containing an appropriate label
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
  }elsif($feature->can('dbID')) {
    $id = $feature->dbID;
  }else{ 
    $id = 'no-id';
  }
  if($feature->can('biotype') && $feature->biotype){
    $id .= "_".$feature->biotype
  }
  return $id;
}



=head2 create_filename

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable
  Arg [2]   : string, stem of filename
  Arg [3]   : string, extension of filename
  Arg [4]   : directory file should live in
  Function  : create a filename containing the PID and a random number
  with the specified directory, stem and extension
  Returntype: string, filename
  Exceptions: throw if directory specifed doesnt exist
  Example   : my $queryfile = $self->create_filename('seq', 'fa');

=cut



sub create_filename{
  my ($stem, $ext, $dir) = @_;
  if(!$dir){
    $dir = '/tmp/';
  }
  $stem = '' if(!$stem);
  $ext = '' if(!$ext);
  throw($dir." doesn't exist Runnable:create_filename") unless(-d $dir);
  my $num = int(rand(100000));
  my $file = $dir."/".$stem.".".$$.".".$num.".".$ext;
  while(-e $file){
    $num = int(rand(100000));
    $file = $dir."/".$stem.".".$$.".".$num.".".$ext;
  }
  return $file;
}

1;
