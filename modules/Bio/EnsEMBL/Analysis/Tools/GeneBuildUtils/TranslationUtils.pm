package Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils;

use strict;
use warnings;
use Exporter;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning
                                      stack_trace_dump);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(id);
use Bio::EnsEMBL::Translation;
use vars qw (@ISA  @EXPORT);


@ISA = qw(Exporter);
@EXPORT = qw(print_Translation clone_Translation 
             print_peptide print_just_Translation 
             starts_with_met ends_with_stop 
             contains_internal_stops);



sub print_Translation{
  my ($transcript,$indent) = @_;
  $indent = '' if(!$indent);
  print_just_Translation($transcript, $indent);
  print_Translation_genomic_coords($transcript, $indent);
  my $met = starts_with_met($transcript);
  print $indent."peptide starts with a methionine\n" if($met);
  my $stop = ends_with_stop($transcript);
  print $indent."peptide ends with a stop codon\n" if($stop);
  my $num_stops = contains_internal_stops($transcript);
  print $indent."peptide contains ".$num_stops.
    " internal stop codons\n" if($num_stops);
}



sub print_peptide{
  my ($transcript) = @_;
  print ">";
  print_just_Translation($transcript, '');
  print $transcript->translate->seq."\n";
}




sub print_just_Translation{
  my ($transcript, $indent) = @_;
  my $translation = $transcript->translation;
  my $id = id($translation);
  my $start_id = id($translation->start_Exon);
  my $end_id = id($translation->end_Exon);
  my $cdna_length = length($transcript->translateable_seq);
  print $indent."TRANSLATION: ".$id." ".$start_id." ".
    $translation->start." ".$end_id." ".$translation->end." ".
      $cdna_length."\n";
}

sub print_Translation_genomic_coords{
  my ($transcript, $indent) = @_;
  print $indent."genomic coords ".
    $transcript->coding_region_start." ".
      $transcript->coding_region_end."\n";
}

sub starts_with_met{
  my ($transcript) = @_;
  my $pep = $transcript->translate->seq;
  return 1 if($pep =~ /^M/);
  return 0;
}

sub ends_with_stop{
  my ($transcript) = @_;
  my $cdna_seq = uc($transcript->translateable_seq);
  return 1 if($cdna_seq =~ /(TAA|TGA|TAG)$/);
  return 0;
}

sub contains_internal_stops{
  my ($transcript) = @_;
  my $pep = $transcript->translate->seq;
  my $num = $pep =~ /\*/;
  return $num;
}

sub clone_Translation{
  my ($transcript, $newtranscript) = @_;
  my $newtranslation = Bio::EnsEMBL::Translation->new();
  my $translation = $transcript->translation;
  my %new_exons;
  foreach my $exon(@{$newtranscript->get_all_Exons}){
    my $id_string = $exon->start."-".$exon->end."-".
      $exon->strand;
    $new_exons{$id_string} = $exon;
  }
  my $start_exon = $translation->start_Exon;
  my $old_start_id_string = $start_exon->start."-".
    $start_exon->end."-".$start_exon->strand;
  my $new_start_Exon = $new_exons{$old_start_id_string};
  $newtranslation->start_Exon;
  $newtranslation->start($translation->start);
  my $end_exon = $translation->end_Exon;
  my $old_end_id_string = $end_exon->end."-".
    $end_exon->end."-".$end_exon->strand;
  my $new_end_Exon = $new_exons{$old_end_id_string};
  $newtranslation->end_Exon;
  $newtranslation->end($translation->end);
  my $attribs = $translation->get_all_Attributes();
  $newtranslation->add_Attributes(@$attribs);
  return $newtranslation;
}


1;
