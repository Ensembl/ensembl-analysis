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


=head2 print_Translation

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : string, indent
  Function  : prints info about translation of the given 
  transcript including its coords, if starts with a met etc
  Returntype: n/a
  Exceptions: none
  Example   : print_Translation($transcript);

=cut



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




=head2 print_peptide

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : prints the peptide sequence of the given 
  transcript
  Returntype: n/a
  Exceptions: none
  Example   : 

=cut



sub print_peptide{
  my ($transcript) = @_;
  print ">";
  print_just_Translation($transcript, '');
  print $transcript->translate->seq."\n";
}



=head2 print_just_Translation

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : string, indent
  Function  : prints just the exon and translation start
  and end info
  Returntype: n/a
  Exceptions: none
  Example   : 

=cut



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


=head2 print_Translation_genomic_coords

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : string, indent
  Function  : print the genomic coordinates of the translation
  Returntype: n/a
  Exceptions: none
  Example   : print_Translation_genomic_coords($transcript);

=cut



sub print_Translation_genomic_coords{
  my ($transcript, $indent) = @_;
  print $indent."genomic coords ".
    $transcript->coding_region_start." ".
      $transcript->coding_region_end."\n";
}



=head2 starts_with_met

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : checks if peptide starts with a methoinine
  Returntype: boolean
  Exceptions: none
  Example   : 

=cut



sub starts_with_met{
  my ($transcript) = @_;
  my $pep = $transcript->translate->seq;
  return 1 if($pep =~ /^M/);
  return 0;
}


=head2 ends_with_stop

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : checks if transcript ends with a stop codon
  Returntype: boolean
  Exceptions: none
  Example   : 

=cut



sub ends_with_stop{
  my ($transcript) = @_;
  my $cdna_seq = uc($transcript->translateable_seq);
  return 1 if($cdna_seq =~ /(TAA|TGA|TAG)$/);
  return 0;
}


=head2 contains_internal_stops

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : counts how many internal stops there are
  Returntype: int
  Exceptions: none
  Example   : 

=cut



sub contains_internal_stops{
  my ($transcript) = @_;
  my $pep = $transcript->translate->seq;
  my $num = $pep =~ /\*/;
  return $num;
}



=head2 clone_Translation

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : Bio::EnsEMBL::Transcript
  Function  : copies the translation of transcript 1 to 
  transcript2
  Returntype: 
  Exceptions: 
  Example   : 

=cut



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
  $newtranslation->start_Exon($new_start_Exon);
  $newtranslation->start($translation->start);
  my $end_exon = $translation->end_Exon;
  my $old_end_id_string = $end_exon->start."-".
    $end_exon->end."-".$end_exon->strand;
  my $new_end_Exon = $new_exons{$old_end_id_string};
  throw($old_end_id_string." failed to get exon") 
   if(!$new_end_Exon);
  $newtranslation->end_Exon($new_end_Exon);
  $newtranslation->end($translation->end);
  my $attribs = $translation->get_all_Attributes();
  $newtranslation->add_Attributes(@$attribs);
  $newtranslation->dbID($translation->dbID);
  return $newtranslation;
}



#METHODS NEEDED

#run_translate, a method to run the orf finder program 

#compute translation, a method which finds a new translation for
#a transcript using run_translate

#return translation, a method if you just want the translation
#object returned and not a whole new transcript

1;
