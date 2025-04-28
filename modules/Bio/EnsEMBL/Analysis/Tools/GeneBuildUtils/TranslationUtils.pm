# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

package Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils;

use strict;
use warnings;
use feature 'say';
use Exporter;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning
                                      stack_trace_dump);
use Bio::EnsEMBL::Analysis::Tools::Logger;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(id);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::EvidenceUtils qw(clone_Evidence);
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(write_seqfile);
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);

use vars qw (@ISA  @EXPORT);


@ISA = qw(Exporter);
@EXPORT = qw(
             print_Translation
             clone_Translation
             print_peptide
             Translation_info
             starts_with_met
             ends_with_stop
             contains_internal_stops
             print_Translation_genomic_coords
             run_translate
             orfs_prediction
             compute_best_translation
             compute_translation
             return_translation
             validate_Translation_coords
             add_ORF_to_transcript
             compute_6frame_translations_for_transcript 
             create_Translation
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
  if ($transcript->translation) {
    print Translation_info($transcript, $indent)."\n";
    print_Translation_genomic_coords($transcript, $indent);
    warning("Transcript is less than 3bp can't print ".
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
  else {
    print $indent, "NO translation\n";
  }
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
  warning("Transcript is less than 3bp can't print ".
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
  $indent = '' if(!$indent);
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
  my $num = ($pep =~ s/\*/\*/g);   
  $num = 0 unless $num > 0;
  logger_info(Translation_info($transcript)." contains internal stops ".$pep)
    if($num);
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
  my ($transcript, $newtranscript, $clone_xrefs) = @_;
  $clone_xrefs = 1 if(!defined($clone_xrefs));
  #print "\n**CLONING TRANSLATION**\n\n";
  my $newtranslation = Bio::EnsEMBL::Translation->new();
  my $translation = $transcript->translation;
  
  my %new_exons;
  foreach my $exon(@{$newtranscript->get_all_Exons}){
    my $id_string = $exon->start."-".$exon->end."-".
      $exon->strand;
    #print "Adding ".$id_string." exon to hash\n";
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
  $newtranslation->stable_id($translation->stable_id); 

  $newtranslation->version($translation->version);
  $newtranslation->created_date($translation->created_date); 
  $newtranslation->modified_date($translation->modified_date); 
  if ($clone_xrefs){
    foreach my $DBEntry (@{$translation->get_all_DBEntries}){
      $newtranslation->add_DBEntry($DBEntry);
    }
  }
  return $newtranslation;
  #print "\n\n**\n\n";
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


=head2 orfs_prediction

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : boolean, do want predictions with met or not
  Function  : Find coordinates of translations in transcripts cdna
  Returntype: arrayref of an array or arrays, each element
  containing an array of length, start and end
  Exceptions: None

=cut

sub orfs_prediction {
  my ($transcript, $met_only) = @_;

  my $seq = $transcript->seq;

  my @orf_predictions;
  foreach my $frame (0, 1, 2) {
    my $translation = $seq->translate(-frame => $frame);
    my $start = $frame+1;
    foreach my $subseq (split('\*', $translation->seq)) {
      my $length_subseq = length($subseq);
      if ($length_subseq) {
        my $first_met = index($subseq, 'M');
        if ($met_only and $first_met != -1) {
          my $length = $length_subseq-$first_met;
          $first_met *= 3;
          push(@orf_predictions, [$length, $start+$first_met, $start+$first_met+(($length+1)*3)-1]);
        }
        elsif (!$met_only and $first_met != 0) {
          push(@orf_predictions, [$length_subseq, $start, $start+(($length_subseq+1)*3)-1]);
        }
      }
      $start += ($length_subseq*3)+3;
    }
    # If the end is longer than the sequence, we probably don't have a stop codon
    if (@orf_predictions and $orf_predictions[-1]->[2] > $seq->length) {
      $orf_predictions[-1]->[2] -= 3;
    }
  }
  my @sorted_predictions = sort { $b->[0] <=> $a->[0] } @orf_predictions;
  return \@sorted_predictions;
}


=head2 compute_best_translation

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : run run_translate and give transcript new
  translation based on those results
  Returntype: Bio::EnsEMBL::Transcript
  Exceptions: Warns in unable to create translation

=cut



sub compute_best_translation{
  my ($transcript) = @_;

  my $met_predictions = orfs_prediction($transcript, 1);
  my $nomet_predictions = orfs_prediction($transcript);

  # choosing the best ORF
  my $orfs;
  if (@$met_predictions && @$nomet_predictions) {
    if ($nomet_predictions->[0]->[0] > (2*$met_predictions->[0]->[0])) {
      $orfs = $nomet_predictions;
    }
    else {
      $orfs = $met_predictions;
    }
  }
  elsif (@$met_predictions) {
    $orfs = $met_predictions;
  }
  elsif (@$nomet_predictions) {
    $orfs = $nomet_predictions;
  }
  else {
    warning(id($transcript)." has no translations ");
    return $transcript;
  }

  # add ORF to transcript
  #Here we take the best prediction with a methionine unless
  #there aren't any of the best prediction without a
  #methoinine is more than twice the length
  my $min_length = $orfs->[0]->[0]*.9;
  my $orf_start = $orfs->[0]->[1];
  my $orf_end = $orfs->[0]->[2];
  my ($translation_start, $translation_end,
      $translation_start_Exon, $translation_end_Exon);
  my $exon_count = 0;
  my $pos = 1;
  my $exons = $transcript->get_all_Exons;
  foreach my $exon (@$exons) {
    $exon_count++;
    if ( $orf_start >= $pos && $orf_start <= $pos + $exon->length - 1 ) {
      $translation_start = $orf_start - $pos+1;
      $translation_start_Exon = $exon;
    }
    if ($orf_end >= $pos && $orf_end <= $pos + $exon->length - 1) {
      $translation_end = $orf_end - $pos + 1;
      $translation_end_Exon = $exon;
    }
    $pos += $exon->length;
  }
  if(!$translation_start || !$translation_end || !$translation_start_Exon || !$translation_end_Exon) {
    throw('problems making the translation for '.id($transcript));
  }
  if ($translation_start_Exon == $translation_end_Exon and @$exons > 1 and @$orfs > 1) {
    my ($better_translation_start, $better_translation_end, $better_translation_start_Exon, $better_translation_end_Exon);
    for (my $index = 1; $index < @$orfs; $index++) {
      last if ($min_length > $orfs->[$index]->[0]);
      my $orf_start = $orfs->[$index]->[1];
      my $orf_end = $orfs->[$index]->[2];
      my $exon_count = 0;
      my $pos = 1;
      foreach my $exon (@$exons) {
        $exon_count++;
        if ( $orf_start >= $pos && $orf_start <= $pos + $exon->length - 1 ) {
          $better_translation_start = $orf_start - $pos+1;
          $better_translation_start_Exon = $exon;
        }
        if ($orf_end >= $pos && $orf_end <= $pos + $exon->length - 1) {
          $better_translation_end = $orf_end - $pos + 1;
          $better_translation_end_Exon = $exon;
        }
        $pos += $exon->length;
      }
      if(!$better_translation_start || !$better_translation_end || !$better_translation_start_Exon || !$better_translation_end_Exon) {
        throw('problems making the translation for '.id($transcript));
      }
      if ($better_translation_start_Exon != $better_translation_end_Exon) {
        $translation_start = $better_translation_start;
        $translation_end = $better_translation_end;
        $translation_start_Exon = $better_translation_start_Exon;
        $translation_end_Exon = $better_translation_end_Exon;
        last;
      }
    }
  }
  my $translation = Bio::EnsEMBL::Translation->new();
  $translation->start($translation_start);
  $translation->end($translation_end);
  $translation->start_Exon($translation_start_Exon);
  $translation->end_Exon($translation_end_Exon);
  $transcript->translation($translation);
  my $found_start = 0;
  my $found_end = 0;
  my $last_end_phase;
  my $first_exon = 1;
  foreach my $exon(@{$transcript->get_all_Exons}){
    $exon->phase(-1);
    $exon->end_phase(-1);

    if($translation->start_Exon == $exon){
      if($translation->start == 1 && $first_exon){
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
      $end_phase = (($exon->length + $exon->phase) %3);
    }
    if(($exon == $translation->end_Exon && $exon->length == $translation->end)){
      $exon->end_phase($end_phase);
    }

    $found_end = 1 if($exon == $translation->end_Exon);

    if (($found_start and !$found_end)) {
      $exon->end_phase($end_phase);
    }

    $last_end_phase = $exon->end_phase;
    $first_exon = 0;
  }
  return $transcript;
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

  # choosing the best ORF 
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
    warning(id($transcript)." has no translations ");
    return $transcript;
  } 

  # add ORF to transcript  
  #Here we take the best prediction with a methionine unless 
  #there aren't any of the best prediction without a 
  #methoinine is more than twice the length 
  
  $transcript = add_ORF_to_transcript($orf,$transcript) ;  

  return $transcript;
}



sub add_ORF_to_transcript{
  my  ($orf,$transcript) = @_; 

  my $orf_start = $orf->[1];
  my $orf_end = $orf->[2];
  my $translation = Bio::EnsEMBL::Translation->new();
  logger_info ("Best orf for ".id($transcript)." ".$orf->[0] ." long ".$orf_start." to ".$orf_end);
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
    warning("problems making the translation ".
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
  my $first_exon = 1;
  # print "Setting phases on transcript after adding translation\n";
  foreach my $exon(@{$transcript->get_all_Exons}){
    $exon->phase(-1);
    $exon->end_phase(-1);

    # print "  Have exon " . $exon->start . " " . $exon->end . "\n";
    if($translation->start_Exon == $exon){
      if($translation->start == 1 && $first_exon){
        $exon->phase(0);
        # print "   setting start phase on it to 0 (tstart = 1 and is start_Exon)\n";
      }
      $found_start = 1;
    }elsif($found_start and not $found_end){
      $exon->phase($last_end_phase);
      # print "   setting start phase on it to last_end_phase ($last_end_phase)\n";
    }
    my $end_phase;
    if($exon == $translation->start_Exon){
      $end_phase = ($exon->end - ($exon->start + 
                                  $translation->start 
                                  - 1) +1 ) %3;
      # print "   start_Exon end phase calculated ($end_phase)\n";
    }else{
      $end_phase = (($exon->length + $exon->phase) %3);
      # print "   end phase calculated ($end_phase)\n";
    }
    if(($exon == $translation->end_Exon && $exon->length == $translation->end)){
      # print "   setting end phase to $end_phase (end exon condition)\n";
      $exon->end_phase($end_phase);
    }

    $found_end = 1 if($exon == $translation->end_Exon);

    if (($found_start and !$found_end)) {
      # print "   setting end phase to $end_phase (found_start and not found_end condition)\n";
      $exon->end_phase($end_phase);
    }

    $last_end_phase = $exon->end_phase;
    $first_exon = 0;
  }
  return $transcript;
} 



=head2 compute_6frame_translations_for_transcript

  Arg [1]   : Bio::EnsEMBL::Transcript 
  Function  : computes all possible 6-frame-translations for a transcript 
              and returns the reference to an array of new Bio::EnsEMBL::Transcript objects. 
              The array contains one "new" transcript created for each unique translation (peptide) sequence
  Returntype: An array of Bio::EnsEMBL::Transcript objects
  Example   : @new_transcripts  = @{ compute_6frame_translations_for_transcript($transcript) } ; 

=cut



sub compute_6frame_translations_for_transcript{
  my ($transcript) = @_;
  my @met_predictions = @{run_translate ($transcript, 1)};
  my @nomet_predictions = @{run_translate ($transcript)}; 
  my %translations; 
  foreach my $orf ( @met_predictions, @nomet_predictions ) {  
    my $nt = new Bio::EnsEMBL::Transcript( -EXONS => $transcript->get_all_Exons  ) ; 
    $nt->biotype($transcript->biotype);  
    $nt= add_ORF_to_transcript($orf,$nt);
    $translations{$nt->translate->seq} = $nt;  
    foreach my $tsf(@{$transcript->get_all_supporting_features}) {
      my $cloned_tsf = clone_Evidence($tsf);
      $nt->add_supporting_features($cloned_tsf);
    }
  }   
  return [values %translations] ;
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

  my $seq = $trans->seq; 
  my $trans_id = id($trans);
  $seq->display_id($trans_id);

  my $file = write_seqfile($seq);
  my $command = "/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/bin/translate";
  $command .= " -m " if($met);
  $command .= " ".$file." | ";
  logger_info($command);
  open ( ORF, $command ) || throw( "Could not run command '$command': $!" );
  
  my @orf_predictions;
 ORF:
  while ( <ORF> ){
    chomp;
    next ORF unless /\>\s*(\S+).*length\s+(\d+),.*\s+(\d+)\.\.(\d+)/;
    my $id = $1;
    my $orf_length = $2;
    my $orf_start = $3;
    my $orf_end   = $4;
    if ($orf_start>=$orf_end ) {  
      # print "can't compute translation for this transcript as translation start >= translation end : $orf_start >= $orf_end \n " ;
      next ORF  ; 
    } 

    # Check if there's a stop codon and add it
    if ($orf_end + 3 <= $seq->length) {
      my $codon = uc($seq->subseq($orf_end+1,$orf_end+3));
      if ($codon eq 'TAA' || $codon eq 'TAG' || $codon eq 'TGA' ) {
        $orf_end+=3;
      }
    }
    my @prediction = ($orf_length,$orf_start,$orf_end);
    push( @orf_predictions, \@prediction );
  }
  close(ORF) || throw("Error running '$command': $!");
  my @sorted_predictions = 
    map { $_->[1] } sort { $b->[0] <=> $a->[0] } map { [$_->[0], $_] } @orf_predictions;
  unlink $file;
  return \@sorted_predictions;
}


=head2 validate_Translation_coords

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : (optional) boolean, allow negative start coordinates
  Function  : check coordinates of translation of transcript
  Returntype: boolean, 1 if translation is ok

=cut

sub validate_Translation_coords{
  my ($transcript, $allow_negative_start) = @_;
  my $translation = $transcript->translation;
  if(!$allow_negative_start && ($translation->start < 1)){
    warning(Translation_info($transcript)." has a start  ".$translation->start.
            " less than 1");
    return 0;
  }
  if($translation->end < 0 || $translation->end > $translation->end_Exon->length){
    warning(Translation_info($transcript)." has an end ".$translation->end.
            " less than 1 or greater than ".$translation->end_Exon->length);
    return 0;
  }
  return 1;
}


=head2 create_Translation

 Arg [1]    : Arrayref of Bio::EnsEMBL::Exons
 Arg [2]    : Int start, genomic start of the cds
 Arg [3]    : Int end, genomic end of the cds
 Description: Create a Bio::EnsEMBL::Translation object based on an array of Bio::EnsEMBL::Exons
 Returntype : Bio::EnsEMBL::Translation
 Exceptions : None

=cut

sub create_Translation {
  my ($exons, $genomic_start, $genomic_end) = @_;

  unless(scalar(@{$exons})) {
    throw("No exons passed in");
  }

  unless($genomic_start && $genomic_start >= 1) {
    throw("Genomic start needs to be a value >= 1");
  }

  unless($genomic_end && $genomic_end >= 1) {
    throw("Genomic end needs to be a value >= 1");
  }

  my $translation = Bio::EnsEMBL::Translation->new();
  foreach my $exon (@$exons) {
    if ($genomic_start >= $exon->seq_region_start and $genomic_start <= $exon->seq_region_end) {
      if ($exon->strand == 1) {
        $translation->start_Exon($exon);
        $translation->start($genomic_start-$exon->seq_region_start+1);
      }
      else {
        $translation->start_Exon($exon);
        $translation->start($exon->seq_region_end-$genomic_start+1);
      }
    }
    if ($genomic_end >= $exon->seq_region_start and $genomic_end <= $exon->seq_region_end) {
      if ($exon->strand == 1) {
        $translation->end_Exon($exon);
        $translation->end($genomic_end-$exon->seq_region_start+1);
      }
      else {
        $translation->end_Exon($exon);
        $translation->end($exon->seq_region_end-$genomic_end+1);
      }
    }
  }

  unless($translation->start) {
    throw("Failed to set a translation start");
  }

  unless($translation->end) {
    throw("Failed to set a translation end");
  }

  unless($translation->start_Exon) {
    throw("Failed to set a start exon");
  }

  unless($translation->end_Exon) {
    throw("Failed to set a end exon");
  }

  return $translation;
}

1;
