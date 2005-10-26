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

use vars qw (@ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(print_Transcript clone_Transcript 
             print_Transcript_evidence);




sub print_Transcript{
  my ($transcript, $indent) = @_;
  my $coord_string = coord_string($transcript);
  my $id = id($transcript);
  $indent = '' if(!$indent);
  print $indent."TRANSCRIPT: ".$id." ".$coord_string."\n";
  my $translation_indent = $indent."\t";
  print_Translation($transcript, $translation_indent);
  foreach my $exon(@{$transcript->get_all_Exons}){
    my $exon_indent = $translation_indent."\t";
    print_Exon($exon, $exon_indent);
  }
}


sub print_Transcript_evidence{
  my ($transcript, $indent) = @_;
  print $indent."TRANSCRIPT EVIDENCE:\n";
  foreach my $evidence(@{$transcript->get_all_supportingfeatures}){
    my $evidence_indent = $indent."\t";
    print_Evidence($evidence, $evidence_indent);
  }
}

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





1;
