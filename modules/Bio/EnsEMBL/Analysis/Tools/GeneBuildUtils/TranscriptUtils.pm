=head1 NAME

Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils - utilities for gene objects

=head1 SYNOPSIS

  use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(clone_Gene);

  or 

  use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils 

  to get all methods

=head1 DESCRIPTION

All methods in this class should take a Bio::EnsEMBL::Transcript
object as their first argument.

The methods provided should carry out some standard 
functionality for said objects such as printing info, and 
cloning and checking phase consistency or splice sites etc

=head1 CONTACT

please send any questions to ensembl-dev@ebi.ac.uk

=head1 METHODS

the rest of the documention details the exported static
class methods

=cut

package Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils;

use strict;
use warnings;
use Exporter;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning
                                      stack_trace_dump);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(print_Translation clone_Translation);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils qw(print_Exon clone_Exon);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(coord_string id);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::EvidenceUtils qw (print_Evidence clone_Evidence);
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_verbosity
                                             logger_info
                                             logger_warning);
use vars qw (@ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(print_Transcript clone_Transcript 
             print_Transcript_evidence print_just_Transcript
             is_strand_consistent lies_inside_of_slice
             exons_shorter_than);




=head2 print_Transcript

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : string, this should be a string or spaces or tabs
  to indent the printed string
  Function  : print information about the transcript and its
  children objects, using indent to make the format readable
  Returntype: n/a
  Exceptions: none
  Example   : print_Transcript($transcript);

=cut


sub print_Transcript{
  my ($transcript, $indent) = @_;
  $indent = '' if(!$indent);
  print_just_Transcript($transcript, $indent);
  my $translation_indent = $indent."\t";
  print_Translation($transcript, $translation_indent);
  foreach my $exon(@{$transcript->get_all_Exons}){
    my $exon_indent = $translation_indent."\t";
    print_Exon($exon, $exon_indent);
  }
}


=head2 print_just_Transcript

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : string, indent
  Function  : print info about just the transcript
  Returntype: n/a
  Exceptions: none
  Example   : print_just_Transcript($transcript);

=cut



sub print_just_Transcript{
  my ($transcript, $indent) = @_;
  $indent = '' if(!$indent);
  my $coord_string = coord_string($transcript);
  my $id = id($transcript);
  print $indent."TRANSCRIPT: ".$id." ".$coord_string."\n";
}

=head2 print_Transcript_evidence

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : string, an indent
  Function  : print the transcripts supporting evidence
  Returntype: n/a
  Exceptions: none
  Example   : print_Transcript_evidence($transcript);

=cut



sub print_Transcript_evidence{
  my ($transcript, $indent) = @_;
  print $indent."TRANSCRIPT EVIDENCE:\n";
  foreach my $evidence(@{$transcript->get_all_supportingfeatures}){
    my $evidence_indent = $indent."\t";
    print_Evidence($evidence, $evidence_indent);
  }
}



=head2 clone_Transcript

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : produce a new copy of the transcript object passed
  in so it can be altered without impact on the original objects
  the only bit it doesnt keep is the adaptor so the cloned 
  object can be stored
  Returntype: Bio::EnsEMBL::Transcript
  Exceptions: none 
  Example   : my $newtranscript = clone_Transcript($transcript);

=cut



sub clone_Transcript{
  my ($transcript) = @_;
  my $newtranscript = new Bio::EnsEMBL::Transcript();
  foreach my $exon(@{$transcript->get_all_Exons}){
    my $newexon = clone_Exon($exon);
    $newtranscript->add_Exon($newexon);
  }
  foreach my $sf(@{$transcript->get_all_supporting_features}){
    my $newsf = clone_Evidence($sf);
    $newtranscript->add_supporting_features($newsf);
  }
  my $newtranslation = clone_Translation($transcript, 
                                         $newtranscript);
  $newtranscript->translation($newtranslation);
  my $attribs = $transcript->get_all_Attributes();
  $newtranscript->add_Attributes(@$attribs);
  $newtranscript->slice($transcript->slice);
  $newtranscript->dbID($transcript->dbID);
  return $newtranscript;
}



=head2 is_strand_consistent

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : checks if strand is consistent between 
  transcript and first exon and in multiexon genes between 
  all exons
  Returntype: boolean, 1 if true undef if not
  Exceptions: none
  Example   : throw("Strands not consistent") 
  if(!is_strand_consistent($transcript));

=cut



sub is_strand_consistent{
  my ($transcript) = @_;
  my $exons = $transcript->get_all_Exons;
  if($exons->[0]->strand != $transcript->strand){
    logger_warning("Strands are inconsistent between the ".
                   "first exon and the transcript for ".
                   id($transcript));
    return undef;
  }
  if(@$exons >= 2){
    for(my $i = 1;$i < @$exons;$i++){
      if($exons->[$i]->strand != $exons->[$i-1]->strand){
        logger_warning("Strands are inconsistent between ".
                       "exon $i exon and exon ".($i-1)." for ".
                       id($transcript));
        return undef;
      }
    }
  }
  return 1;
}



=head2 lies_inside_of_slice

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : Bio::EnsEMBL::Slice
  Function  : ensures the transcript within the slice, 
  completely on the lower end, it can overhang the upper end
  Returntype: boolean, 1 for true undef for false
  Exceptions: none
  Example   : 

=cut


sub lies_inside_of_slice{
  my ($transcript, $slice) = @_;
  if($transcript->start > $slice->length || 
     $transcript->end < 1){
    logger_warning(id($transcript)." lies off edge if slice ".
                   $slice->name);
    return undef;
  }
  if($transcript->start < 1 && $transcript->end > 1){
    logger_warning(id($transcript)." lies over lower boundary".
                   " of slice ".$slice->name);
  }
  return 1;
}



=head2 exons_shorter_than

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : int, max length
  Function  : checks if any of the exons of given transcript
  are longer than specified length
  Returntype: boolean, 1 for true, undef for false
  Exceptions: 
  Example   : 

=cut



sub exons_shorter_than{
  my ($transcript, $max_length) = @_;
  foreach my $exon(@{$transcript->get_all_Exons}){
    if($exon->length > $max_length){
      logger_warning(id($exon)." from ".id($transcript).
                     " is longer than max length ".
                     $max_length);
      return undef;
    }
  }
  return 1;
}


##METHODS NEEDED

#checks

#phase consistency (account of -1/[012] exception)
#folded transcript, exon starting before previous ends!
#suspect evidence (is this still needed to avoid NGs)
#low complexity
#is_spliced, does contain at least one exon longer than 9bps and how many?
#canonical splice sites, does it have all canoical splice sites?
#perhaps also what percentage is useful?
#orf coverage, how does the spliced length compare to the genomic extent?

#utilitys
#split trancripts, split transcripts on long introns
#replace stops with introns, frameshift around stop codons!
#list_evidence, a list of ids that support the gene

1;
