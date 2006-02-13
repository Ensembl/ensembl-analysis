package Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils;

use strict;
use warnings;
use Exporter;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning
                                      stack_trace_dump);
use Bio::EnsEMBL::Analysis::Tools::Logger;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(id);
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::SequenceUtils qw(write_seqfile);

use vars qw (@ISA  @EXPORT);


@ISA = qw(Exporter);
@EXPORT = qw(print_Translation
             clone_Translation
             print_peptide
             Translation_info
             starts_with_met
             ends_with_stop
             contains_internal_stops
             print_Translation_genomic_coords
             run_translate
             compute_translation
             return_translation
            );


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
  print Translation_info($transcript, $indent)."\n";
  print_Translation_genomic_coords($transcript, $indent);
  logger_warning("Transcript is less than 3bp can't print ".
                 "any more info") if($transcript->length < 3);
  return if($transcript->length < 3);
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
  logger_warning("Transcript is less than 3bp can't print ".
                 "Its peptide") if($transcript->length < 3);
  return if($transcript->length < 3);
  print ">";
  print Translation_info($transcript, '')."\n";
  print $transcript->translate->seq."\n";
}



=head2 Translation_info

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : string, indent
  Function  : 
  Returntype: returns a string of translation start exon,
  start, end exon, end, cdna length
  Exceptions: none
  Example   : 

=cut



sub Translation_info{
  my ($transcript, $indent) = @_;
  my $translation = $transcript->translation;
  my $id = id($translation);
  my $start_id = id($translation->start_Exon);
  my $end_id = id($translation->end_Exon);
  my $cdna_length = length($transcript->translateable_seq);
  return $indent."TRANSLATION: ".$id." ".$start_id." ".
    $translation->start." ".$end_id." ".$translation->end." ".
      $cdna_length." ";
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
  throw("Failed to find exon ".$old_start_id_string." for ".
        id($transcript)) if(!$new_start_Exon);
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


=head2 return_translation

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : calculates translation for transcript
  Returntype: Bio::EnsEMBL::Translation
  Exceptions: 
  Example   : 

=cut



sub return_translation{
  my ($trans) = @_;
  $trans = compute_translation($trans);
  return $trans->translation;
}


=head2 compute_translation

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : run run_translate and give transcript new
  translation based on those results
  Returntype: Bio::EnsEMBL::Transcript
  Exceptions: Warns in unable to create translation
  Example   : 

=cut



sub compute_translation{
  my ($transcript) = @_;

  my @met_predictions = @{run_translate 
                            ($transcript, 1)};
  my @nomet_predictions = @{run_translate 
                              ($transcript)};

  my $orf;
  if(@met_predictions && @nomet_predictions){
    my $met_best = $met_predictions[0];
    my $no_met_best = $nomet_predictions[0];
    if($no_met_best->[0] > (2*$met_best->[0])){
      $orf = $no_met_best;
    }else{
      $orf = $met_best;
    }
  }elsif(@met_predictions){
    $orf = $met_predictions[0];
  }elsif(@nomet_predictions){
    $orf = $nomet_predictions[0];
  }else{
    logger_warning(id($transcript)." has no translations ");
    return $transcript;
  }
  #Here we take the best prediction with a methionine unless 
  #there aren't any of the best prediction without a 
  #methoinine is more than twice the length

  my $orf_start = $orf->[1];
  my $orf_end = $orf->[2];
  my $translation = Bio::EnsEMBL::Translation->new();
  logger_info("Best orf for ".id($transcript)." ".$orf->[0].
              " long ".$orf_start." to ".$orf_end);
  my ($translation_start, $translation_end, 
      $translation_start_Exon, $translation_end_Exon);
  my $exon_count = 0;
  my $pos = 1;
  foreach my $exon(@{$transcript->get_all_Exons}){
    $exon_count++;
    logger_info("exon:$exon_count exon_length:".$exon->length." pos:$pos orf_start:$orf_start orf_end:$orf_end pos+:".($pos + $exon->length - 1));
    if ( $orf_start >= $pos && $orf_start <= $pos 
         + $exon->length - 1 ){
      $translation_start = $orf_start - $pos+1;
      $translation_start_Exon = $exon;
    }
    if($orf_end >= $pos && $orf_end <= $pos 
       + $exon->length - 1){
      $translation_end = $orf_end - $pos + 1;
      $translation_end_Exon = $exon;
    }
    $pos += $exon->length;
  }
  if(!$translation_start || !$translation_end 
     || !$translation_start_Exon
     || !$translation_end_Exon){
    logger_warning("problems making the translation ".
                   "for ".id($transcript));
    return $transcript;
  }else{
    $translation->start($translation_start);
    $translation->end($translation_end);
    $translation->start_Exon($translation_start_Exon);
    $translation->end_Exon($translation_end_Exon);
    $transcript->translation($translation);
  }

  my $found_start = 0;
  my $found_end = 0;
  my $last_end_phase;
  foreach my $exon(@{$transcript->get_all_Exons}){
    $exon->phase(-1);
    $exon->end_phase(-1);
    if($translation->start_Exon == $exon){
      if($translation->start == 1){
        $exon->phase(0);
      }
      $found_start = 1;
    }elsif($found_start and not $found_end){
      $exon->phase($last_end_phase);
    }
    my $end_phase;
    if($exon == $translation->start_Exon){
      $end_phase = ($exon->end - ($exon->start + 
                                  $translation->start 
                                  - 1) +1 ) %3;
    }else{
      $end_phase = (($exon->length + $exon->phase) %3)
    }
    $exon->end_phase($end_phase) 
      if(($exon == $translation->end_Exon
         && $exon->length == $translation->end) ||
         ($found_start and not $found_end));

    $found_end = 1 if($exon == $translation->end_Exon);
    $last_end_phase = $exon->end_phase;
  }
  return $transcript;
}

=head2 run_translate

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : boolean, do want predictions with met or 
  not
  Function  : run the program translate and find coordinates
  of translations in transcripts cdna
  Returntype: arrayref of an array or arrays, each element 
  containing an array of length, start and end
  Exceptions: 
  Example   : 

=cut



sub run_translate{
  my ($trans,$met) = @_;

  my $trans_id = id($trans);
  my $seq = $trans->seq;
  $seq->display_id($trans_id);

  my $file = write_seqfile($seq);
  my $command = "/usr/local/ensembl/bin/translate";
  $command .= " -m " if($met);
  $command .= " ".$file." | ";
  logger_info($command);
  open ( ORF, $command ) || throw( "Error running translate" );
  
  my @orf_predictions;
 ORF:
  while ( <ORF> ){
    chomp;
    next ORF unless /\>/;
    my @entries = split;
    next ORF unless ( $entries[3] && $entries[5] );
    my $id = $entries[1];
    my $orf_length = $entries[3];
    $orf_length =~s/\,//;
    $entries[5] =~/(\d+)\.\.(\d+)/;
    my $orf_start = $1;
    my $orf_end   = $2;
    next ORF if $orf_start>=$orf_end;
    my @prediction = ($orf_length,$orf_start,$orf_end);
    push( @orf_predictions, \@prediction );
  }
  my @sorted_predictions = 
    map { $_->[1] } sort { $b->[0] <=> $a->[0] } map { [$_->[0], $_] } @orf_predictions;
  unlink $file;
  return \@sorted_predictions;
}




1;
