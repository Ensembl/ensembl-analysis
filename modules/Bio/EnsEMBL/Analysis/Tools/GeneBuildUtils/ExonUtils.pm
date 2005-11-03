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
@EXPORT = qw(print_Exon clone_Exon print_just_Exon);



sub print_Exon{
  my ($exon, $indent) = @_;
  throw("Must be passed an exon") if(!$exon);
  print_just_Exon($exon, $indent);
  foreach my $evidence(@{$exon->get_all_supporting_features}){
    my $evidence_indent = $indent."\t";
    print_Evidence($evidence, $evidence_indent);
  }
}


sub print_just_Exon{
  my ($exon, $indent) = @_;
  throw("Must be passed an exon") if(!$exon);
  $indent = '' if(!$indent);
  my $coord_string = coord_string($exon);
  my $id = id($exon);
  print $indent."EXON: ".$id." ".$exon->start." ".$exon->end." ".$exon->strand." ".$exon->phase."\n";
}

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


#METHODS

#transfer supporting evidence, a method to transfer the 
#supporting evidence between 2 exons

1;
