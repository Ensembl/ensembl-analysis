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
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_verbosity
                                             logger_info
                                             logger_warning);
use vars qw (@ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
             print_Exon
             clone_Exon
             Exon_info
             exon_length_less_than_maximum
             transfer_supporting_evidence
            );




=head2 print_Exon

  Arg [1]   : Bio::EnsEMBL::Exon
  Arg [2]   : string, indent (\t) to tag infront of printed
  Arg [3]   : boolean (1 if evidence attached to exon should be printend) 
  Function  : print information about given exon and
  associated evidence
  Returntype: none
  Exceptions: throws if not passed a first argument
  Example   : 

=cut



sub print_Exon{
  my ($exon, $indent, $print_ev) = @_;
  $indent = "" if(!$indent);
  throw("Must be passed an exon") if(!$exon);
  $print_ev = 1 unless defined $print_ev ;  
  print Exon_info($exon, $indent)."\n";
  if ($print_ev) { 
    foreach my $evidence(@{$exon->get_all_supporting_features}){
      my $evidence_indent = $indent."\t";
      print_Evidence($evidence, $evidence_indent);
    }
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
  return $indent."EXON: ".$id." ".$coord_string." phase ".
    $exon->phase." end_phase ".$exon->end_phase." length ".
      $exon->length;
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
  $newexon->stable_id  ($exon->stable_id);
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



=head2 transfer_supporting_evidence

  Arg [1]   : Bio::EnsEMBL::Exon
  Arg [2]   : Bio::EnsEMBL::Exon
  Function  : transfer evidence from the first exon to the 
  second exon
  Returntype: Bio::EnsEMBL::Exon
  Exceptions: 
  Example   : 

=cut

#Note this method modified the target exon which is 
#passed in

sub transfer_supporting_evidence{
  my ($source_exon, $target_exon) = @_;
  my %target_evidence;
  my %source_evidence;
  foreach my $sf
    (@{$target_exon->get_all_supporting_features}){
      my $unique_id = $sf->start."-".$sf->end."-".
        $sf->strand."-".$sf->hseqname."-".$sf->cigar_string;
      $target_evidence{$unique_id} = $sf;
    }
  foreach my $sf
    (@{$source_exon->get_all_supporting_features}){
      my $unique_id = $sf->start."-".$sf->end."-".
        $sf->strand."-".$sf->hseqname."-".$sf->cigar_string;
      if(!$target_evidence{$unique_id}){
        $source_evidence{$unique_id} = $sf;
      }
    }

  foreach my $sf(values(%source_evidence)){
    logger_info("Adding ".$sf->hseqname." to the target exon");
    $target_exon->add_supporting_feature($sf);
  }
  return $target_exon;
}


1;
