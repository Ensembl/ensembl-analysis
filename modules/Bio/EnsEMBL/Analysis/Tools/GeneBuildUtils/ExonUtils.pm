=head1 NAME

Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils - utilities for transcript objects

=head1 SYNOPSIS

  use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils qw(Exon_info);

  or 

  use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils 

  to get all methods

=head1 DESCRIPTION

All methods in this class should take a Bio::EnsEMBL::Exon
object as their first argument.

The methods provided should carry out some standard 
functionality for said objects such as printing info or
transfering evidence

=head1 CONTACT

please send any questions to ensembl-dev@ebi.ac.uk

=head1 METHODS

the rest of the documention details the exported static
class methods

=cut


package Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils;

use strict;
use warnings;
use Exporter;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning
                                      stack_trace_dump);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(coord_string id);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::EvidenceUtils qw(print_Evidence clone_Evidence);
use Bio::EnsEMBL::Exon;

use vars qw (@ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
             print_Exon
             clone_Exon
             Exon_info
             exon_length_less_than_maximum
            );




=head2 print_Exon

  Arg [1]   : Bio::EnsEMBL::Exon
  Arg [2]   : string, indent (\t) to tag infront of printed
  string
  Function  : print information about given exon and
  associated evidence
  Returntype: none
  Exceptions: throws if not passed a first argument
  Example   : 

=cut



sub print_Exon{
  my ($exon, $indent) = @_;
  throw("Must be passed an exon") if(!$exon);
  print Exon_info($exon, $indent)."\n";
  foreach my $evidence(@{$exon->get_all_supporting_features}){
    my $evidence_indent = $indent."\t";
    print_Evidence($evidence, $evidence_indent);
  }
}



=head2 Exon_info

  Arg [1]   : Bio::EnsEMBL::Exon
  Arg [2]   : indent, string
  Function  : returns a string with info about exon in it
  Returntype: string
  Exceptions: throws if not passed a first argument
  Example   : 

=cut



sub Exon_info{
  my ($exon, $indent) = @_;
  throw("Must be passed an exon") if(!$exon);
  $indent = '' if(!$indent);
  my $coord_string = coord_string($exon);
  my $id = id($exon);
  return $indent."EXON: ".$id." ".$coord_string." ".
    $exon->phase." ".$exon->end_phase;
}


=head2 clone_Exon

  Arg [1]   : Bio::EnsEMBL::Exon
  Function  : produces a whole new exon object identical
  to the object passed
  Returntype: Bio::EnsEMBL::Exon
  Exceptions: none
  Example   : 

=cut



sub clone_Exon{
  my ($exon) = @_;
  my $newexon = Bio::EnsEMBL::Exon->new();
  $newexon->start      ($exon->start);
  $newexon->end        ($exon->end);
  $newexon->phase      ($exon->phase);
  $newexon->end_phase  ($exon->end_phase);
  $newexon->strand     ($exon->strand);
  $newexon->dbID       ($exon->dbID);
  $newexon->slice      ($exon->slice);
  $newexon->analysis   ($exon->analysis);
  foreach my $sf(@{$exon->get_all_supporting_features}){
    my $newsf = clone_Evidence($sf);
    $newexon->add_supporting_features($sf);
  }
  return $newexon;
}



=head2 exon_length_less_than_maximum

  Arg [1]   : Bio::EnsEMBL::Exon
  Arg [2]   : int, max length
  Function  : checks if exon is longer than specified length
  Returntype: boolean, 0 for longer, 1 for shorter
  Exceptions: none
  Example   : 

=cut



sub exon_length_less_than_maximum{
  my ($exon, $max_length) = @_;
  if($exon->length >= $max_length){
    logger_warning(id($exon)." is longer than max length ".
                   $max_length);
    return 0;
  }
  return 1;
}


#METHODS

#transfer supporting evidence, a method to transfer the 
#supporting evidence between 2 exons

1;
