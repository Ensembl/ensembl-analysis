# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

=head1 NAME

Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils - utilities for transcript objects

=head1 SYNOPSIS

  use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(clone_Transcript);

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

please send any questions to http://lists.ensembl.org/mailman/listinfo/dev

=head1 METHODS

the rest of the documention details the exported static
class methods

=cut

package Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils;

use strict;
use warnings;
use Exporter qw(import);
use feature 'say';

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning stack_trace_dump);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(print_Translation clone_Translation print_peptide);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils qw(print_Exon clone_Exon Exon_info exon_length_less_than_maximum Exon_info merge_exons get_upstream_Intron get_downstream_Intron get_upstream_splice_sites get_downstream_splice_sites validate_Exon_coords);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::IntronUtils qw(intron_length_less_than_maximum get_splice_sites);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(seq_region_coord_string id empty_Object);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::EvidenceUtils qw (print_Evidence clone_Evidence);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(write_seqfile align_proteins_with_alignment);
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);

use Bio::EnsEMBL::Analysis;
use Bio::SeqIO;

use POSIX qw/ceil/;

our @EXPORT_OK = qw(
             all_exons_are_valid
             are_phases_consistent
             are_splice_sites_canonical
             are_strands_consistent
             attach_Slice_to_Transcript
             attach_Analysis_to_Transcript
             attach_Analysis_to_Transcript_no_support
             attach_Analysis_to_Transcript_no_overwrite
             calculate_exon_phases
             clone_Transcript
             coding_coverage
             convert_to_genes
             convert_translateable_exons_to_exon_extended_objects
             count_non_canonical_splice_sites
             count_real_introns
             dump_cDNA_file
             empty_Transcript
             evidence_coverage
             evidence_coverage_greater_than_minimum
             exon_lengths_all_less_than_maximum
             exon_overlap
             coding_exon_overlap
             exonic_proportion
             features_overlap
             fully_load_Transcript
             get_downstream_Intron_from_Exon
             get_downstream_splice_sites
             get_evidence_ids
             get_upstream_Intron_from_Exon
             get_upstream_splice_sites
             has_no_unwanted_evidence
             identical_Transcripts
             intron_lengths_all_less_than_maximum
             is_not_folded
             is_Transcript_sane
             is_spliced
             list_evidence
             overlap_length
             print_Transcript
             print_Transcript_and_Exons
             print_Transcript_evidence
             remove_initial_or_terminal_short_exons
             remove_short_frameshift_introns
             replace_stops_with_introns
             set_start_codon
             set_stop_codon
             split_Transcript
             tidy_split_transcripts
             trim_cds_to_whole_codons
             Transcript_info
             has_polyA_signal
             set_alignment_supporting_features
            );



=head2 convert_translateable_exons_to_exon_extended_objects 

  Arg [0]   : Bio::EnsEMBL::Transcript 
  Function  : converts all translateable Exons of a Bio::EnsEMBL::Transcript object 
              into  ExonExtended objects to access addtional attributes, like previouus exon etc 
  Returntype: Arrayref of Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended objects 
  Example   : 

=cut


sub convert_translateable_exons_to_exon_extended_objects {
  my ( $transcript ) = @_ ;

  my $exons = $transcript->get_all_translateable_Exons;
  my @conv_exons ;

   for (my $i=0 ; $i < scalar ( @$exons) ; $i++) {
      my $pte = ${$exons}[$i] ;
      bless $pte,"Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended" ;

      if ( ($i+1) <scalar(@$exons)) {
        $pte->next_exon(${$exons}[$i+1]) ;
      }
      if ( ($i-1)>=0) {
        $pte->prev_exon(${$exons}[$i-1]) ;
      }
      $pte->transcript($transcript) ;
      push @conv_exons , $pte ;
   }
   return \@conv_exons ;
}






=head2 print_Transcript

  Arg [1]   : Bio::EnsEMBL::Transcript or Aref of Bio::EnsEMBL::Transcript-objects
  Arg [2]   : string, this should be a string or spaces or tabs to indent the 
              printed string

  Function  : print information about the transcript and its
  children objects, using indent to make the format readable
  Returntype: none
  Exceptions: none
  Example   : print_Transcript($transcript);

=cut


sub print_Transcript{
  my ($tref, $indent) = @_;

  $indent = '' if(!$indent);
  my @transcripts ; 
  if (ref($tref)=~m/ARRAY/){
    @transcripts = @$tref; 
  }else{
    push @transcripts, $tref; 
  }
  for my $transcript ( @transcripts ) { 
    print Transcript_info($transcript, $indent)."\n";
    print_Transcript_evidence($transcript);
    my $translation_indent = $indent."\t";
    print_Translation($transcript, $translation_indent) ;
    foreach my $exon(@{$transcript->get_all_Exons}){
      my $exon_indent = $translation_indent."\t";
      print_Exon($exon, $exon_indent);
    }
  }
}


=head2 print_Transcript_and_Exons

  Arg [1]   : Bio::EnsEMBL::Transcript or Aref of Bio::EnsEMBL::Transcript-objects
  Arg [2]   : string, this should be a string or spaces or tabs to indent the 
              printed string

  Function  : print information about the transcript and its Exons
              using indent to make the format readable
  Returntype: none
  Exceptions: none
  Examples  : print_Transcript_and_Exons(\@transcript) 
              print_Transcript_and_Exons($transcript) 

=cut


sub print_Transcript_and_Exons{
  my ($tref, $indent) = @_;
  
  $indent = '' if(!$indent);

  my @tr = (ref($tref)=~m/ARRAY/) ? @$tref : $tref ; 

  for my $transcript ( @tr )  { 
    print "\n" .  Transcript_info($transcript, $indent)."\n"; 
    #map { print Exon_info($_, $indent."\t")."\n"  } @{$transcript->get_all_Exons} ; 
    my $count = 1;
    foreach my $exon(@{$transcript->get_all_Exons}){
      print $indent."\t".$count." ".Exon_info($exon)."\n";
      $count++;
    }
  }
}



=head2 Transcript_info

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : string, indent
  Function  : return string of info about the transcript
  Returntype: String
  Exceptions: none
  Example   : print_just_Transcript($transcript);

=cut



sub Transcript_info{
  my ($transcript, $indent) = @_;
  $indent = '' if(!$indent);
  my $coord_string = seq_region_coord_string($transcript); 
  my $id = $transcript->display_id;
  my $logic_name = $transcript->analysis ? $transcript->analysis->logic_name : 'NO_ANALYSIS';
  return $indent."TRANSCRIPT: ".$id." ".$coord_string.' '.$logic_name;
}

=head2 print_Transcript_evidence

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : string, an indent
  Function  : print the transcripts supporting evidence
  Returntype: none
  Exceptions: none
  Example   : print_Transcript_evidence($transcript);

=cut



sub print_Transcript_evidence{
  my ($transcript, $indent) = @_;
  $indent = '' if(!$indent);
  print $indent."TRANSCRIPT EVIDENCE:\n";
  foreach my $evidence(@{$transcript->get_all_supporting_features}){
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
  my ($transcript, $clone_xrefs) = @_; 

  $clone_xrefs = 1 if(!defined($clone_xrefs));
  my $newtranscript = new Bio::EnsEMBL::Transcript();
  foreach my $exon(@{$transcript->get_all_Exons}){
    
    my $newexon = clone_Exon($exon);
    $newtranscript->add_Exon($newexon);
  }
  foreach my $sf(@{$transcript->get_all_supporting_features}){
    my $newsf = clone_Evidence($sf);
    $newtranscript->add_supporting_features($newsf);
  }
  if ( $transcript->get_all_IntronSupportingEvidence) {
    foreach my $ise(@{$transcript->get_all_IntronSupportingEvidence}){
      $newtranscript->add_IntronSupportingEvidence($ise);
    }
  }
  my $newtranslation;
  if($transcript->translation){
    $newtranslation = clone_Translation($transcript, 
                                        $newtranscript, $clone_xrefs);
  }
  $newtranscript->translation($newtranslation);
  my $attribs = $transcript->get_all_Attributes();
  $newtranscript->add_Attributes(@$attribs);
  $newtranscript->slice($transcript->slice);
  $newtranscript->biotype($transcript->biotype);
  $newtranscript->dbID($transcript->dbID);
  $newtranscript->stable_id($transcript->stable_id);
  $newtranscript->version($transcript->version);
  $newtranscript->source($transcript->source);
  $newtranscript->analysis($transcript->analysis);
  if ($clone_xrefs){
    foreach my $DBEntry (@{$transcript->get_all_DBEntries}){
      $newtranscript->add_DBEntry($DBEntry);
    } 
   $newtranscript->display_xref($transcript->display_xref); 
  } 
  
  return $newtranscript;
}



=head2 are_strands_consistent

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : checks if strand is consistent between 
  transcript and first exon and in multiexon genes between 
  all exons
  Returntype: boolean, 1 if pass 0 if fail (ie strands 
                                                inconsistent)
  Exceptions: none
  Example   : throw("Strands not consistent") 
  if(!are_strands_consistent($transcript));

=cut



sub are_strands_consistent{
  my ($transcript) = @_;

  my $exons = $transcript->get_all_Exons;
  if($exons->[0]->strand != $transcript->strand){
    warn("Strands are inconsistent between the ".
            "first exon and the transcript for ".
            $transcript->display_id);
    return 0;
  }
  if(@$exons >= 2){
    for(my $i = 1;$i < @$exons;$i++){
      if($exons->[$i]->strand != $exons->[$i-1]->strand){
        warn("Strands are inconsistent between ".
                "exon $i exon and exon ".($i-1)." for ".
                $transcript->display_id);
        return 0;
      }
    }
  }
  return 1;
}





=head2 exon_lengths_all_less_than_maximum

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : int, max length
  Function  : checks if any of the exons of given transcript
  are longer than specified max length
  Returntype: boolean, 1 for pass, 0 for fail (ie exon beyond
                                                   max length)
  Exceptions: none
  Example   : 

=cut



sub exon_lengths_all_less_than_maximum{
  my ($transcript, $max_length) = @_;

  foreach my $exon(@{$transcript->get_all_Exons}){
    if(!exon_length_less_than_maximum($exon, $max_length)){
      warn("Transcript ".$transcript->display_id." has ".
                  "exon longer than ".$max_length);
      return 0;
    }
  }
  return 1;
}

=head2 intron_lengths_all_less_than_maximum

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : int, max length
  Function  : checks if any of the introns of given transcript
  are longer than specified max length
  Returntype: boolean, 1 for pass, 0 for fail 
  (ie intron beyond max length)
  Exceptions: none
  Example   : 

=cut



sub intron_lengths_all_less_than_maximum{
  my ($transcript, $max_length) = @_;

  foreach my $intron(@{$transcript->get_all_Introns}){
    if(!intron_length_less_than_maximum($intron, 
                                        $max_length)){
      #warning("Transcript ".id($transcript)." has ".
      #        "intron ".$intron->length." longer than ".$max_length);
      return 0;
    }
  }
  if(@{$transcript->get_all_Introns} == 0){
    my $warn = "intron_lengths_all_less_than_maximum is an ".
      "inappropriate test for a single exon gene";
    logger_info($warn);
  }
  return 1;
}


=head2 are_phases_consistent

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : to check that end phase of one exon is always
  the same as start phase of the next. The only acceptable
  exception to this is an end phase of -1 and a start phase
  of 0/1/2 or vice versa as UTR exons always have end phase
  of -1 but the whole of the next exon may be coding and
  to give it a start phase of -1 is considered misleading
  Returntype: boolean, 1 for pass 0 for fail (ie phase 
                                                  inconsistency)
  Exceptions: none
  Example   : 

=cut



sub are_phases_consistent{
  my ($transcript) = @_;

  my $exons = $transcript->get_all_Exons;

  if (not $transcript->translation) {
    # all phases should be -1
    foreach my $e (@$exons) {
      if ($e->phase != -1 or $e->end_phase != -1) {
        warn("Non-coding transcript does not have -1 phases");
#        return 0;
      }
    }
  } else {
    my $tr = $transcript->translation;

    for(my $i=1;$i < @$exons;$i++){
      my $prev_end_phase = $exons->[$i-1]->end_phase;
      my $this_phase = $exons->[$i]->phase;

      if($prev_end_phase != $this_phase) {
        # exception: allow a -1<->0 transition if exon is
        # translation start exon
        if ($prev_end_phase == -1 and 
            $this_phase == 0 and
            $tr->start_Exon == $exons->[$i]) {
          # okay
          next;
        } elsif ($prev_end_phase == 0 and
                 $this_phase == -1 and
                 $tr->end_Exon == $exons->[$i-1]) {
          # okay
          next;
        } else {
          warn("Coding transcript has inconsistent phases");
          return 0;
        }
      }
    }
  }

  return 1;
}



=head2 calculate_exon_phases

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : Given a transcript, calculates and sets
    exon phases according to translation and given
    start phase

=cut

sub calculate_exon_phases {
  my ($transcript, $start_phase) = @_;

  foreach my $e (@{$transcript->get_all_Exons}) {
    $e->phase(-1);
    $e->end_phase(-1);
  }

  if ($transcript->translation) {
    my $tr = $transcript->translation;

    my @exons = @{$transcript->get_all_Exons};

    while($exons[0] != $tr->start_Exon) {
      shift @exons;
    }
    while($exons[-1] != $tr->end_Exon) {
      pop @exons;
    }

    # set phase of for first coding exon
    my $cds_len = $exons[0]->length;
    if ($tr->start == 1) {
      $exons[0]->phase($start_phase);
      if ($start_phase > 0) {
        $cds_len += $start_phase;
      }
    } else {
      $cds_len -= ($tr->start - 1);
    }
    $exons[0]->end_phase($cds_len % 3);
    # set phase for internal coding exons      
    for(my $i=1; $i < @exons; $i++) {
      $exons[$i]->phase($exons[$i-1]->end_phase);
      $exons[$i]->end_phase(($exons[$i]->length + $exons[$i]->phase) % 3);
    }
        
    # set phase for last coding exon
    if ($exons[-1]->length > $tr->end) {
      $exons[-1]->end_phase(-1);
    }
  }
}

=head2 is_not_folded

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : check if any exons start before the previous exon
  finished
  Returntype: Boolean, 1 for pass, 0 for fail 
  Exceptions: 
  Example   : 

=cut



sub is_not_folded{
  my ($transcript) = @_;

  my $exons = $transcript->get_all_Exons;
  if(@$exons == 1){
    my $warn = "is_not_folded ".
      "is an inappropriate test for a single ".
        "exon gene";
    logger_info($warn);
    return 1;
  }
  for(my $i = 1;$i < @$exons;$i++){
    $exons->[$i]->stable_id('');
    if($exons->[$i]->strand == 1){
      if($exons->[$i]->start < $exons->[$i-1]->end){
        warning($transcript->display_id." is folded");
        warn($i." ".id($exons->[$i])." has a start which ".
                    "is less than ".($i-1)." ".id($exons->[$i-1]).
                    " end");
        return 0;
      }
    }else{
      if($exons->[$i]->end > $exons->[$i-1]->start){
        warning($transcript->display_id." is folded");
        warn($i." ".id($exons->[$i])." has a end which ".
                    "is greater than ".($i-1)." ".id($exons->[$i-1]).
                    " start");
        return 0;
      }
    }
  }
  return 1;
}

=head2 has_no_unwanted_evidence

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : hashref, a hash of ids which are unwanted
  Function  : To ensure the given transcript is not supported
  by any of the unwanted evidence
  Returntype: boolean, 1 for no unwanted evidence, 0 for
  unwanted evidence
  Exceptions: 
  Example   : 

=cut



sub has_no_unwanted_evidence{
  my ($transcript, $ids) = @_;
  $ids = {} if(!$ids);
  $ids->{'NG_'} = 1;

  my $evidence = get_evidence_ids($transcript);
  foreach my $evi(keys(%$evidence)){
    foreach my $unwanted (keys(%$ids)) {
      if($evi =~ /$unwanted/){
        warn($transcript->display_id." has ".$evi.
                " unwanted evidence");
        return 0;
      }
    }
  }
  return 1;
}



=head2 get_evidence_ids

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : gets a hashref of all the evidence supporting
  the given transcript
  Returntype: hashref
  Exceptions: 
  Example   : 

=cut



sub get_evidence_ids{
  my ($transcript) = @_;
  my %hash;
  foreach my $sf(@{$transcript->get_all_supporting_features}){
    $hash{$sf->hseqname} = 1;
  }
  foreach my $exon(@{$transcript->get_all_Exons}){
    foreach my $sf(@{$exon->get_all_supporting_features}){
      $hash{$sf->hseqname} = 1;
    }
  }
  return \%hash;
}



=head2 is_spliced

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : int, intron size
  Function  : calculates is the transcript contains at least 
  one "real" intron
  Returntype: boolean 1 if it does, 0 if not
  Exceptions: 
  Example   : 

=cut



sub is_spliced{
  my ($transcript, $intron_size) = @_;
  my $count = count_real_introns($transcript, $intron_size);
  warning($transcript->display_id." has no introns ".
          "longer than $intron_size bps") if(!$count);
  return 0 if(!$count);
  return 1;
}



=head2 count_real_introns

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : int, intron size
  Function  : counts the number of introns longer than
  the specified size, by default this is 9bp
  Returntype: int, the number of real introns
  Exceptions: 
  Example   : 

=cut


sub count_real_introns{
  my ($transcript, $intron_size) = @_;
  $intron_size = 9 if(!$intron_size);
  my $real_count = 0 ;
  my @introns = @{$transcript->get_all_Introns};
  my $warn = $transcript->display_id." has no introns. count_real".
    "_introns makes no sense in those terms";
  logger_info($warn) if(@introns == 0);
  foreach my $intron(@introns){
    $real_count++ if($intron->length > $intron_size);
  }
  logger_info($transcript->display_id." has ".$real_count." introns ".
              "longer than ".$intron_size." out of ".
              @introns." introns");
  return $real_count;
}



=head2 are_splice_sites_canonical

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : check if all splice sites are canonical GT/AG,
  GC/AG or AT/AC pairings
  Returntype: boolean, 1 for yes 0 for no
  Exceptions: 
  Example   : 

=cut



sub are_splice_sites_canonical{
  my ($transcript) = @_;
  my $introns = $transcript->get_all_Introns;
  my $non_canonical_count = 
    count_non_canonical_splice_sites($transcript);
  if($non_canonical_count){
    warn(Transcript_info($transcript)." contains ".
                $non_canonical_count." non canonical ".
                "splice sites out of ".@$introns.
                " introns");
    return 0;
  }
  if(@$introns == 0){
    my $warn ="are_splice_sites_canonical is an ".
      "inappropriate test for a single exon gene"; 
    logger_info($warn);
  }
  return 1;
}



=head2 count_non_canonical_splice_site

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : count the number of non canonical splice sites
  Returntype: int, the number of non canonical splice sites
  Exceptions: 
  Example   : 

=cut



sub count_non_canonical_splice_sites{
  my ($transcript) = @_;
  my $slice = $transcript->slice;
  my $none = 0;
  foreach my $intron(@{$transcript->get_all_Introns}){
    my ($upstream_site, $downstream_site) =
      get_splice_sites($intron);
    $none++ unless(($upstream_site eq 'GT' && 
                    $downstream_site eq 'AG') || 
                   ($upstream_site eq 'AT' && 
                    $downstream_site eq 'AC') ||
                   ($upstream_site eq 'GC' && 
                    $downstream_site eq 'AG') )
  }
  return $none;
}



=head2 exonic_proportion

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : calculates what proportion of the transcript
  extent is made up of exonic sequence
  Returntype: int
  Exceptions: 
  Example   : 

=cut



sub exonic_proportion{
  my ($transcript) = @_;
  my $genomic_extent = ($transcript->end -
                        $transcript->start) + 1;
  my $cdna_length = $transcript->length;
  my $value = ($cdna_length/$genomic_extent) * 100;
  return $value;
}



=head2 coding_coverage

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : calculate the ratio of the translateable seq
  of the transcript to the full length sequence
  Returntype: int
  Exceptions: 
  Example   : 

=cut



sub coding_coverage{
  my ($transcript) = @_;
  my $cdna_length = $transcript->length;
  my $coding_seq_length = 
    length($transcript->translateable_seq);
  my $value = ($coding_seq_length/$cdna_length) * 100;
  return $value;
}



=head2 list_evidence

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : produce a unique list of evidence to support
  a gene
  Returntype: listref 
  Exceptions: 
  Example   : 

=cut


sub list_evidence{
  my ($transcript) = @_;
  my $hash = get_evidence_ids($transcript);
  my @ids =  keys(%$hash);
  return \@ids;
}



=head2 split_Transcript

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : int, max intron length
  Arg [3]   : boolean, 1 to carry out checks, 0 to not be default
  is 1
  Function  : to split transcripts on introns which are
  too long, discard any single exon transcripts left behind 
  and return the remaining
  Returntype: arrayref of Bio::EnsEMBL::Transcript objects
  Exceptions: throws if not passed a Bio::EnsEMBL::Transcript
  object

=cut

#This method is designed to split transcripts on
#long introns and hopefully maintain clean and sensible
#transcripts

sub split_Transcript{
  my ($transcript, $max_intron_length, $intron_numbers) = @_;

  throw("TranscriptUtils:split_Transcript will not work on single exon transcripts")
    if(@{$transcript->get_all_Exons} == 0);
  #cloning the transcript to ensure if all else fails
  #the original is fine
  my %intron_numbers;
  if (defined $intron_numbers) {
    map { $intron_numbers{$_} = 1 } @$intron_numbers;
  }
  my ($cds_start, $cds_end);
  if ($transcript->translation) {
    my ($c) = $transcript->cdna2genomic($transcript->cdna_coding_start, 
                                        $transcript->cdna_coding_start);
    $cds_start = $c->start;
    
    ($c) = $transcript->cdna2genomic($transcript->cdna_coding_end,
                                     $transcript->cdna_coding_end);
    $cds_end   = $c->start;

    if ($cds_start > $cds_end) {
      ($cds_start, $cds_end) = ($cds_end, $cds_start);
    }
  }
  throw("cds_start ".$cds_start." or cds_end.".$cds_end." is not defined") if(!$cds_start || !$cds_end);

  my $cloned_transcript = clone_Transcript($transcript);

  #creating first new transcript
  my $curr_transcript = new Bio::EnsEMBL::Transcript;
  $curr_transcript->biotype($transcript->biotype);
  $curr_transcript->analysis($transcript->analysis);
  my @split_transcripts;
  
  my $first_exon = 1;
  my $intron_count = 0;
  
  #adding first new transcript to list
  push(@split_transcripts, $curr_transcript);
  my $last_exon;
  #checking each intron
 INTRON:
  foreach my $intron(@{$cloned_transcript->get_all_Introns}){
    $intron_count++;

    my $prev_exon = $intron->prev_Exon;
    my $next_exon = $intron->next_Exon;
    
    #If considering first intron and prev exon
    #is the first exon then need to add it to the 
    #transcript and check if it is the start of translation
    if($first_exon){
      $curr_transcript->add_Exon($prev_exon);
      $first_exon = 0;
    }
    
    #If the intron length is less than the specified max
    #then you can add it to the current transcript and
    #check if it is the start of translation or end of
    #translation
    if((not defined $max_intron_length or 
       intron_length_less_than_maximum($intron, $max_intron_length)) and
       not exists $intron_numbers{$intron_count}) {
      $curr_transcript->add_Exon($next_exon);
      next INTRON;
    }else{
      #If the intron is longer than maximum length a new
      #transcript must be started. This means the last exon 
      #needs to be set as the end of the previous translation
      #and a new transcript must be created

      my $t = Bio::EnsEMBL::Transcript->new;
      $t->add_Exon($next_exon);
      $curr_transcript = $t;
      $curr_transcript->biotype($transcript->biotype);
      $curr_transcript->analysis($transcript->analysis);
      $last_exon = $next_exon;
      push(@split_transcripts, $curr_transcript);
    }
  }
 
  # now add the translations, and trim back transcripts that
  # do not end in whole codons
  my @processed;
  
  foreach my $stran (@split_transcripts) {
    if ($transcript->translation) {
      if ($stran->seq_region_end >= $cds_start and
          $stran->seq_region_start <= $cds_end) {
        # at least part of this transcript is coding
        my $tr = Bio::EnsEMBL::Translation->new;
        $stran->translation($tr);

        my @exons = @{$stran->get_all_Exons};
        foreach my $e (@exons) {
          if ($cds_start >= $e->seq_region_start and $cds_start < $e->seq_region_end) {
            # start of translation is in this exon
            if ($stran->strand > 0) {
              $tr->start_Exon($e);
              $tr->start( $cds_start - $e->seq_region_start + 1);
            } else {
              $tr->end_Exon($e);
              $tr->end( $e->seq_region_end - $cds_start + 1);
            }
          }
          if ($cds_end >= $e->seq_region_start and $cds_end <= $e->seq_region_end) {
            if ($stran->strand > 0) {
              $tr->end_Exon($e);
              $tr->end( $cds_end - $e->seq_region_start + 1);
            } else {
              $tr->start_Exon($e);
              $tr->start( $e->seq_region_end - $cds_end + 1);
            }
          }
        }

        if (not $tr->start_Exon) {
          $tr->start_Exon($exons[0]);
          $tr->start(1);
        }
        if (not $tr->end_Exon) {
          $tr->end_Exon($exons[-1]);
          $tr->end($exons[-1]->length);
        }
      }
    }
    if ($stran->translation) {
      #$stran = trim_cds_to_whole_codons($stran);
    }
    foreach my $sf (@{$transcript->get_all_supporting_features}) {
      my @ugs;
      foreach my $ug ($sf->ungapped_features) {
        if ($ug->start >= $stran->start and 
            $ug->end <= $stran->end) {
          push @ugs, $ug;
        }
      }
      if (@ugs) {
        my $newf = $sf->new(-features => \@ugs) if(@ugs);
        $stran->add_supporting_features($newf);
      }
    }
    push @processed, $stran;
  }

  my $info = "split_Transcript Returning " . scalar(@processed) . " transcripts";
  logger_info($info);
  return \@processed;
}


=head2 tidy_split_transcripts

  Arg [1]   : arrayref Bio::EnsEMBL::Transcript
  Arg [2]   : Bio::EnsEMBL::Transcript
  Function  : Intended to be called after split_Transcripts
    Performs a series of standard gene-buildy checks on the
    given transcripts,
  Args      : The original transcripts, and the results of split_Transcript
  Returns   : arrayref of transcripts that pass checks

=cut

sub tidy_split_transcripts{
  my ($orig_tran, $split_trans) = @_;

  my @keep;

  foreach my $stran (@{$split_trans}) {
    if(@{$stran->get_all_Exons} == 1){
      my ($exon) = @{$stran->get_all_Exons};
      my $warn = $stran->display_id." only has one exon\n".
          Exon_info($exon);
      warning($warn);
      next;
    }

    if($stran->translate->seq =~ /\*/){
      my $warn = Transcript_info($stran).
          " does not translate\n";
      warning($warn);
      next;
    }

    my $initial_peptide = $orig_tran->translate->seq;    

    if(!$initial_peptide =~ /$stran->translate->seq/){
      my $warn = Transcript_info($stran)." translation ".
          "does not appear to be a subset of the original ".
          "translation";
      warning($warn);
      next;
    }

    #and if the strands or phases are not consistent
    if(!are_strands_consistent($stran)){
      my $warn = Transcript_info($stran)." has ".
          "inconsistent strands";
      warning($warn);
      next;
    }

    if(!are_phases_consistent($stran)){
      my $warn = Transcript_info($stran)." has ".
          "inconsistent phases";
      warning($warn);
      next;
    }

    push @keep, $stran;
  }

  return(\@keep);
}


=head2 trim_cds_to_whole_codons

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : when start/end of translation is flush to
    start/end of transcript, and phase is non-zero, trim
    back to whole codon. If we have UTR adjacent to a non-zero
    phase codon, then something is wrong, so best to leave it
    alone    
  Returntype: Bio::EnsEMBL::Transcript 
  Exceptions: 
  Example   : 

=cut

sub trim_cds_to_whole_codons {
  my ($transcript) = @_;
  my $remove5 = 0;
  my $remove3 = 0;  

  if ($transcript->translation) {
    my @exons = @{$transcript->get_all_Exons};
    my $tr = $transcript->translation;

    if ($tr->start_Exon == $exons[0] and 
        $tr->start_Exon->phase > 0 and
        $tr->start == 1) {
      $remove5 = (3 - $tr->start_Exon->phase) % 3;
    }
    if ($tr->end_Exon == $exons[-1] and 
        $tr->end_Exon->end_phase > 0 and
        $tr->end == $tr->end_Exon->length) {
      $remove3 = $tr->end_Exon->end_phase;
    }

    if ($remove5 or $remove3) {
      my $cloned_transcript = clone_Transcript($transcript);
      my @exons = @{$cloned_transcript->get_all_Exons};
      my $cloned_tr = $cloned_transcript->translation;
      if ($remove5) {
        while(@exons and $exons[0]->length <= $remove5) {
          $remove5 -= $exons[0]->length;
          shift @exons;
        }
        $cloned_transcript->flush_Exons;
        map { $cloned_transcript->add_Exon($_) } @exons;

        if ($cloned_transcript->strand > 0) {
          $exons[0]->start($exons[0]->start + $remove5);
        } else {
          $exons[0]->end($exons[0]->end - $remove5);
        }
        $exons[0]->phase(0);
        $cloned_tr->start_Exon($exons[0]);
        $cloned_tr->start(1);
      }
      if ($remove3) {
        while(@exons and $exons[-1]->length <= $remove3) {
          $remove3 -= $exons[-1]->length;
          pop @exons;
        }
        if ($cloned_transcript->strand > 0) {
          $exons[-1]->end($exons[-1]->end - $remove3);
        } else {
          $exons[-1]->start($exons[-1]->start + $remove3);
        }
        $exons[-1]->end_phase(0);
        $cloned_tr->end_Exon($exons[-1]);
        $cloned_tr->end($exons[-1]->length);
      }
      return $cloned_transcript;
    }
  } 
  return $transcript;
}


=head2 replace_stops_with_introns

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : int, max stop count
  Function  : replace any inframe stops with
  introns
  Returntype: Bio::EnsEMBL::Transcript 
  Exceptions: method returns zero if translation contains stop codon adjacent to gap. 
  Example   :  
      my $new_tr = replace_stops_with_introns($tr);
      if ($new_tr == 0) {
        print STDERR "Skipping transcript with internal stop codon next to gapped sequence: ".$tr->dbID."\n";
        next;
      }
  Notes     : Needs extension to deal with transcript_supporting_feats

=cut

sub replace_stops_with_introns{
  my ($transcript,$max_stops) = @_;

  my $translation_start_shift = 0; # in number of bases
  my $translation_end_shift = 0;   # in number of bases
  my $end_exon_shift = 0;          # in number of exons

  #foreach my $exon (@{$transcript->get_all_Exons}) {
  #  print "DEBUG: Exon ".$exon->start."-".$exon->end.":".$exon->strand.' ^ '.$exon->rank($transcript).' @ '.$exon->length.' sp '.$exon->phase. ' ep '.$exon->end_phase."\n";
  #}

  my $newtranscript = clone_Transcript($transcript);
  my @exons = @{$newtranscript->get_all_Exons};
  my $pep = $newtranscript->translate->seq;
  #print 'DEBUG: peptide: ', $pep, "\n";
  my $removed_exon_count = 0;
  # gaps adjacent to internal stop codons - skip
# I'm not sure this happens often but I had a peptide equal to '*'
  return 0 if ($pep eq '*' or $pep =~ /X\*/ || $pep =~ /\*X/);

  my $num_stops = $pep =~ s/\*/\*/g;

  # The next few bits of code are to do some checks in terms of the allowed number of stops and to throw
  # if there is an issue. There is at least one redundant check later in the code, however it is faster
  # to have these here

  # This is the default behaviour, throw if > 1 stop
  if ($num_stops > 1 && !$max_stops) {
    throw("Transcript does not have exactly one stop codon; it has $num_stops stops. Multiple stops replacement is ".
          "experimental and requires a value to be passed to max_stops");
  }

  # If max_stops is not a sensible value then throw
  if(defined($max_stops) && $max_stops <= 0) {
    throw("You have passed a value to max_stops but this value is <= 0. The value passed to max_stops should be ".
          ">= 1");
  }

  # If max_stops is defined and there are more stops than the value then throw
  if(defined($max_stops) && ($num_stops > $max_stops)) {
    warning("You have set max_stops to ".$max_stops." however the number of stops in the translation is ".$num_stops);
    return(0);
  }

  # Warn that there are internal stops
  if ($num_stops > 0) {
    warning("Transcript has ".$num_stops." internal stops\n");
  }

  while($pep =~ /\*/g) {
    # find the position of the stop codon within the peptide
    my $position = pos($pep);
    print "Replacing stop at pos: ".$position."\n";

    # and find out the genomic start position of this stop codon
    my @coords = $newtranscript->pep2genomic($position, $position);

    foreach my $stop (@coords) {
      print "Found stop at position start ".$stop->start." end ".$stop->end." on strand ".$transcript->strand."\n";
      # locate the exon that this stop lies in
      my @new_exons;
      foreach my $exon (@exons) {
        #print 'DEBUG: ', $exon->rank($newtranscript), ' $$ ', $exon->seq->seq, "\n";
        # NOTE that at this point the stop will always lie on a translateable exon
        if ($stop->start > $exon->start and $stop->end < $exon->end) {
          # This stop lies _completely_ within an exon and not on its
          # boundary. We therefore can split the exon into two UNLESS
          # (the stop starts at the start of the translation OR
          #  the stop ends at the end of the translation)

          if ( ($transcript->translation->start_Exon->start == $exon->start and
                $transcript->translation->genomic_start() == $stop->start and
                $transcript->strand == 1)
               or
               ($transcript->translation->start_Exon->start == $exon->start and
                $transcript->translation->genomic_end() == $stop->end and
                $transcript->strand == -1) ) {
            # if the stop starts at the start of the translation
            # the translation start is shifted and the exon is not divided

            $translation_start_shift += 3; # translation start and end are "stranded", no need to look at strand
            print("The stop starts at the start of the translation within an exon, not boundary.\n");

            push @new_exons,$exon; # exon not changed but translation start will be shifted
            next;

          } elsif ( ($transcript->translation->end_Exon->start == $exon->start and
                     $transcript->translation->genomic_end() == $stop->end and
                     $transcript->strand == 1)
                    or
                    ($transcript->translation->end_Exon->start == $exon->start and
                     $transcript->translation->genomic_start() == $stop->start and
                     $transcript->strand == -1) )
          {

            # if the stop ends at the end of the translation
            # the translation end is shifted and the exon is not divided

            $translation_end_shift -= 3; # translation start and end are "stranded", no need to look      at strand
            print("The stop ends at the end of the translation within an exon, not boundary.\n");

            push @new_exons,$exon; # exon not changed but translation start will be shifted
            next;

          } else {
            # the stop end DOES NOT match the translation end
            $end_exon_shift += 1;
            if ($transcript->translation->end_Exon->start == $exon->start) {
              # and this is the last translateable exon
              if ($transcript->strand == 1) {
                $translation_end_shift -= ($stop->end-$exon->start+1); # translation start and end are "stranded"
              } else {
                $translation_end_shift -= ($exon->end-$stop->start+1);
              }
            }
          }

          print("---I am NOT a boundary stop\n");

          my $exon_left = Bio::EnsEMBL::Exon->
              new(-slice     => $exon->slice,
                  -start     => $exon->start,
                  -end       => $stop->start - 1,
                  -strand    => $exon->strand,
                  -phase     => $exon->strand < 0 ? 0 : $exon->phase,
                  -end_phase => $exon->strand < 0 ? $exon->end_phase  :0);
          my $exon_right = Bio::EnsEMBL::Exon->
              new(-slice     => $exon->slice,
                  -start     => $stop->end + 1,
                  -end       => $exon->end,
                  -strand    => $exon->strand,
                  -phase     => $exon->strand < 0 ? $exon->phase : 0,
                  -end_phase => $exon->strand < 0 ? 0 : $exon->end_phase);

          my @sfs = @{$exon->get_all_supporting_features}; 
          my (@ug_left, @ug_right);

          my $feature_isa = undef;
          my $feature_unit_length = undef;

          foreach my $f (@sfs) {
          	if ($f->isa("Bio::EnsEMBL::DnaDnaAlignFeature")) {

          	  if ($feature_isa and ($feature_isa ne "Bio::EnsEMBL::DnaDnaAlignFeature")) {
          	  	throw("All the supporting features must be of the same type.");
          	  }

          	  $feature_unit_length = 1; # ie feature is a cDNA alignment
          	  $feature_isa = "Bio::EnsEMBL::DnaDnaAlignFeature";

                } elsif ($f->isa("Bio::EnsEMBL::DnaPepAlignFeature")) {

         	  if ($feature_isa and ($feature_isa ne "Bio::EnsEMBL::DnaPepAlignFeature")) {
                throw("All the supporting features must be of the same type.");
              }

         	  $feature_unit_length = 3; # ie feature is a protein alignment
          	  $feature_isa = "Bio::EnsEMBL::DnaPepAlignFeature";

         	} else {
          	  throw("Feature ".$f->dbID()." is not Bio::EnsEMBL::DnaDnaAlignFeature nor Bio::EnsEMBL::DnaPepAlignFeature.");
          	}

            foreach my $ug ($f->ungapped_features) {
              $ug->analysis($newtranscript->analysis);

              my $orignial_analysis = $ug->analysis;
              if ($ug->start >= $exon_left->start &&
                  $ug->end <= $exon_left->end+3) {
                # completely within the left-side of the split
                # (including the stop length +3)
                push @ug_left, $ug;

              } elsif($ug->end >= $exon_left->start &&
                      $ug->end <= $exon_left->end+3) {
                # This case might crop up if the feature went over the edge of the end of the
                # left exon. Possibly because of a previous stop removal. I'm keeping this as
                # a separate case to draw attention to the possibility
                warning("Feature only partially overlaps with left exon. Will add anyway.");
                push @ug_left, $ug;

              } elsif ($ug->start >= $exon_right->start-3 && 
                       $ug->end <= $exon_right->end) {
                # completely within the right-side of the split
                # (including the stop length -3)
                push @ug_right, $ug;

              } elsif($ug->start >= $exon_right->start-3 &&
                      $ug->start <= $exon_right->end) {
                # This case might crop up if the feature went over the edge of the end of the
                # right exon. Possibly because of a previous stop removal. I'm keeping this as
                # a separate case to draw attention to the possibility
                warning("Feature only partially overlaps with right exon. Will add anyway.");
                push @ug_right, $ug;

              } elsif ($ug->start >= ($exon_left->start-3) && $ug->end <= ($exon_right->end+3)) {

                # this ug must span the split
                my $fp_left = Bio::EnsEMBL::FeaturePair->new();
                if ($ug->slice) {
                  $fp_left->slice($ug->slice);
              }
              $fp_left->seqname   ($ug->seqname);
              $fp_left->strand    ($ug->strand);
              $fp_left->hstrand    ($ug->hstrand);
              $fp_left->hseqname  ($ug->hseqname);
              $fp_left->score     ($ug->score);
              $fp_left->percent_id($ug->percent_id);
              $fp_left->start     ($ug->start);
              $fp_left->end       ($stop->start - 1);
              $fp_left->external_db_id($ug->external_db_id);
              $fp_left->hcoverage($ug->hcoverage);
              $fp_left->analysis($orignial_analysis) ;           

              my $fp_right = Bio::EnsEMBL::FeaturePair->new();
              if ($ug->slice) {
                $fp_right->slice($ug->slice);
              }
              $fp_right->seqname   ($ug->seqname);
              $fp_right->strand    ($ug->strand);
              $fp_right->hstrand    ($ug->hstrand);
              $fp_right->hseqname  ($ug->hseqname);
              $fp_right->score     ($ug->score);
              $fp_right->percent_id($ug->percent_id);
              $fp_right->start     ($stop->end + 1);
              $fp_right->end       ($ug->end);
              $fp_right->external_db_id($ug->external_db_id);
              $fp_right->hcoverage($ug->hcoverage);
              $fp_right->analysis($orignial_analysis) ;           

                  # here's the state of play:
                  #
                  #                        fp_left          fp_right
                  #                         s    e        s          e
                  #                         ======        ============
                  #          1 >---------------------------------------------> strand
                  # 
                  #                         s                        e
                  #             un gapped   ==========================
                  #
                  #         -1 <---------------------------------------------< strand
                  #                         ======        ============
                  #                         s    e        s          e
                  #                        fp_right         fp_left
                  #

              if ($exon->strand > 0) {
                $fp_left->hstart($ug->hstart);
                $fp_left->hend($ug->hstart +
                               ceil($fp_left->length / $feature_unit_length) -
                               1);
                $fp_right->hend ($ug->hend);
                $fp_right->hstart($ug->hend -
                                  ceil($fp_right->length / $feature_unit_length) +
                               1);
              } else {
                $fp_right->hstart($ug->hstart);
                $fp_right->hend($ug->hstart +
                               ceil($fp_right->length / $feature_unit_length) -
                               1);
                $fp_left->hend ($ug->hend);
                $fp_left->hstart($ug->hend -
                                  ceil($fp_left->length / $feature_unit_length) +
                               1);
              }
               
#              if ($exon->strand < 0) {
#                # if we are on the reverse strand
#                # we swap the right and the left
#                my $tmp_fp = $fp_left;
#                $fp_left = $fp_right;
#                $fp_right = $tmp_fp;
#              }
               
              if ($fp_left->end >= $fp_left->start) {
                push @ug_left, $fp_left;
              }
              if ($fp_right->end >= $fp_right->start) {
                push @ug_right, $fp_right;
              }
            }

            elsif($ug->start < $exon_left->start && $ug->end < $exon_left->start ||
                  $ug->start > $exon_right->end && $ug->end > $exon_right->end
                 ) {
              warning("Feature is present but lies fully outside the left and right exons, not adding");
            }

            else {
              throw("Something about this feature has not been covered in the conditionals, edit code");
            }

          } # foreach my $ug ($f->ungapped_features) {
        } # foreach my $f (@sfs) {

          $exon_left = add_dna_align_features_by_hitname_and_analysis(\@ug_left,$exon_left,$feature_isa) ;
          $exon_right =add_dna_align_features_by_hitname_and_analysis(\@ug_right,$exon_right,$feature_isa) ;

          if ($exon->strand < 0) {
            if ($exon_right->end >= $exon_right->start) {
              push @new_exons, $exon_right;
            }
            if ($exon_left->end >= $exon_left->start) {
              push @new_exons, $exon_left;
            }
          } else {
            if ($exon_left->end >= $exon_left->start) {
              push @new_exons, $exon_left;
            }
            if ($exon_right->end >= $exon_right->start) {
              push @new_exons, $exon_right;
            }
          }
        } elsif($stop->start == $exon->start && $stop->end == $exon->end) {

          warning("Exon is a stop codon, removing the exon");

          # This will later be added to the end_exon_index to account for the removed
          # exon or exons. This works on the test case and seems sensible, but is
          # difficult to thoroughly test.
          $end_exon_shift -= 1;
        } elsif ($stop->start == $exon->start) {
          # stop lies at the start of the exon
          print("---stop lies at the start of the exon\n");
          # note that +3 has been replaced with $stop->end-$stop->start+1 to
          # fix the rare cases where stops lie on two consecutive exons
          my $orig_exon_start = $exon->start;
          $exon->start($exon->start + ($stop->end-$stop->start+1));
          # Because the stop length may not now be 3 bases long now we need to fix the phase
          if ( $transcript->strand == -1 ) {
            $exon->end_phase(0);
          } else {
            $exon->phase(0);
          }
          
          #print "ttees: " .$transcript->translation->end_Exon->start . "\n";
          #print "exonstart: " . $exon->start . "\n";
          if ($transcript->translation->end_Exon->start == $orig_exon_start and
                 $transcript->strand == 1) {
            # this is the last translateable exon on the forward strand
            $translation_end_shift -= $stop->end-$stop->start+1;
          }

          push @new_exons, $exon;
          # I'm removing the call to the truncate sub because it still needs work, sometimes it truncates things
          # it doesn't need to. As we are short on time the simplist thing is to just use the old way, which is
          # the above line of code, where no modification occurs to the supporting features. Might leave slight
          # overhang, but this is no big deal
#          my $new_exon = truncate_exon_features($exon,$newtranscript->analysis,0);
#          push @new_exons, $new_exon;

        } elsif ($stop->end == $exon->end ) {
          # stop lies at the end of the exon
          print("---stop lies at the end of the exon\n");
          # note that +3 has been replaced with $stop->end-$stop->start+1 to
          # fix the rare cases where stops lie on two consecutive exons
          #print "DB8 e end: ". $exon->end. " s sta: ". $stop->start. " s end: " .$stop->end. " s len: " .$stop->length."\n"; 
          $exon->end($exon->end - ($stop->end-$stop->start+1));  
          
          # Because the stop length may not now be 3 bases long now we need to fix the phase
          if ( $transcript->strand == -1 ) {
            $exon->phase(0);
          } else {
            $exon->end_phase(0);
          }

          if ($transcript->translation->end_Exon->start == $exon->start and
              $transcript->strand == -1) {
            # this is the last translateable exon on the reverse strand
            $translation_end_shift -= $stop->end-$stop->start+1;
          }

          push @new_exons, $exon;
          # I'm removing the call to the truncate sub because it still needs work, sometimes it truncates things
          # it doesn't need to. As we are short on time the simplist thing is to just use the old way, which is
          # the above line of code, where no modification occurs to the supporting features. Might leave slight
          # overhang, but this is no big deal
#          my $new_exon = truncate_exon_features($exon,$newtranscript->analysis,1);
#          push @new_exons, $new_exon;

        } else {
          # this exon is unaffected by this stop
          print("---this exon is unaffected by this stop\n");
          push @new_exons, $exon;
        }

      }

      @exons = @new_exons;
    } # foreach

  } #end of while loop;  by this time, we hope there are not stop codons in the peptide

  # Useful for debugging: keep.
  #foreach my $exon (@exons) {
  #  print "DEBUG: NEW Exon ".$exon->start."-".$exon->end.":".$exon->strand.' ^ '.$exon->rank($transcript).' @ '.$exon->length.' sp '.$exon->phase. ' ep '.$exon->end_phase."\n";
  #}


  # this removes the old exons and replaces with new exon
  # by first cloning the old transcript and then replacing the exon
  # we should be keeping the info attached to the transcript eg. xrefs, sequence edits (atttribs)
  $newtranscript->flush_Exons;
  foreach my $exon (@exons) {
    $newtranscript->add_Exon($exon);
  }

  # assign translation values
  my $translation = Bio::EnsEMBL::Translation->new();

  # get new start exon
  my $start_exon_index = 0; # I'll use 1-based exon count.
  foreach my $exon (@{$transcript->get_all_Exons}) {
    $start_exon_index++;
    if ($transcript->translation->start_Exon->start == $exon->start) {
      last;
    }
  }
  my $new_start_exon;
  foreach my $exon (@{$newtranscript->get_all_Exons}) {
    $start_exon_index--;
    if ($start_exon_index <= 0) {
      $new_start_exon = $exon;
      last;
    }
  }
  $translation->start_Exon($new_start_exon);

  # get new end exon
  my $end_exon_index = 0; # I'll use 1-based exon count.
  foreach my $exon (@{$transcript->get_all_Exons}) {
    $end_exon_index++;
    if ($transcript->translation->end_Exon->start == $exon->start) {
      last;
    }
  }

  #print("DEBUG old transcript number of exons:".@{$transcript->get_all_Exons}."\n");
  #print("DEBUG new transcript number of exons:".@{$newtranscript->get_all_Exons}."\n");

  #print("DEBUG old end_exon_index is $end_exon_index\n");
  my $new_end_exon_index = $end_exon_index+$end_exon_shift; # I'll use 1-based exon count.
 #print("DEBUG new end_exon_index is $new_end_exon_index\n");

  my $new_end_exon;
  foreach my $exon (@{$newtranscript->get_all_Exons}) {
    #print("DEBUG new end exon ".$exon->start." ".$exon->end."\n");
    $new_end_exon_index--;
    if ($new_end_exon_index <= 0) {
      $new_end_exon = $exon;
      last;
    }
  }

  $translation->end_Exon($new_end_exon);
  $translation->start($transcript->translation->start + $translation_start_shift);
  $translation->end($transcript->translation->end + $translation_end_shift);
  $newtranscript->translation($translation);

  my $old_transl_len = $transcript->translate->length();
  my $new_transl_len = $newtranscript->translate->length();

  # The first case to test for is if the edited translation is longer than the original, this shouldn't happen
  if($new_transl_len > $old_transl_len) {
    throw("The edited transcript has a longer translation than the original, something has gone wrong.".
          " Original translation has length ".$old_transl_len.", edited translation has length ".$new_transl_len.
          "\n>original\n".$transcript->translate->seq."\n>edited\n".$newtranscript->translate->seq);
  }

  # As long as the new translation is less than or equal to the original length do a few more checks.
  # If the max_stops parameter is active, make sure the edited translation has had a length of at
  # least the original length minus the value of max_stops
  elsif($max_stops && (($old_transl_len-$max_stops) > $new_transl_len)) {
    throw("The edited transcript is shorter than allowed by the max_stops parameter. Currently max_stops is ".
          "set as: ".$max_stops."\nThe original translation has length ".$old_transl_len.", edited translation has length ".
          $new_transl_len."\n>original\n".$transcript->translate->seq."\n>edited\n".$newtranscript->translate->seq);
  }

  # The next issue is if max_stops is not present and the edited translation has had more than one stop removed.
  # By default this should expect to only remove one.
  elsif(!$max_stops && (($old_transl_len-1) > $new_transl_len)) {
    throw("The edited transcript had more than one internal stop removed. The default is to allow one removal. If you ".
          "want to remove more than one interal stop you can pass in the max_stops parameter. "."\nThe original translation has length ".
          $old_transl_len.", edited translation has length ".$new_transl_len."\n>original\n".$transcript->translate->seq.
          "\n>edited\n".$newtranscript->translate->seq);
  }

  # Hopefully at this point the removal of the stops is okay.
  else {
    print "Original translation has length ".$old_transl_len." but edited translation has length ".$new_transl_len.
          "\n>original\n".$transcript->translate->seq."\n>edited\n".$newtranscript->translate->seq."\n";
  }

  # add xrefs from old translation to new translation
  my $old_translation = $transcript->translation  ;  
  foreach my $DBEntry (@{$old_translation->get_all_DBEntries}){
     $translation->add_DBEntry($DBEntry);
  }

  return $newtranscript;
}


=head2 truncate_exon_features

  Arg [1]  : reference to an array of ungapped features
  Arg [2]  : 0 (truncate from start) or 1 (truncate from end)
  Function : this is called when the start or end of the exon is a stop codon
             In this case we want to recalculate the supporting features

=cut

sub truncate_exon_features {

  my ($exon,$analysis,$start_or_end) = @_;

  my $new_exon   = Bio::EnsEMBL::Exon->
  new(-slice     => $exon->slice,
    -start     => $exon->start,
    -end       => $exon->end,
    -strand    => $exon->strand,
    -phase     => $exon->phase,
    -end_phase => $exon->end_phase );

  # This is currently making the assumption that the feature spans over
  # the length of the exon. If for some reason it didn't then perhaps a
  # check could go here to make sure the processing is actually needed

  my @sfs = @{$exon->get_all_supporting_features};
  my (@ug);

  my $feature_isa = undef;
  my $feature_unit_length = undef;

  foreach my $f (@sfs) {
    if ($f->isa("Bio::EnsEMBL::DnaDnaAlignFeature")) {

      if ($feature_isa and ($feature_isa ne "Bio::EnsEMBL::DnaDnaAlignFeature")) {
        throw("All the supporting features must be of the same type.");
      }

      $feature_unit_length = 1; # ie feature is a cDNA alignment
      $feature_isa = "Bio::EnsEMBL::DnaDnaAlignFeature";

    } elsif ($f->isa("Bio::EnsEMBL::DnaPepAlignFeature")) {

      if ($feature_isa and ($feature_isa ne "Bio::EnsEMBL::DnaPepAlignFeature")) {
        throw("All the supporting features must be of the same type.");
      }

      $feature_unit_length = 3; # ie feature is a protein alignment
      $feature_isa = "Bio::EnsEMBL::DnaPepAlignFeature";

    } else {
      throw("Feature ".$f->dbID()." is not Bio::EnsEMBL::DnaDnaAlignFeature nor Bio::EnsEMBL::DnaPepAlignFeature.");
    }

    foreach my $ug ($f->ungapped_features) {
      $ug->analysis($analysis);
      my $orignial_analysis = $ug->analysis;
      my $trunc_feature = Bio::EnsEMBL::FeaturePair->new();
      if ($ug->slice) {
        $trunc_feature->slice($ug->slice);
      }
      $trunc_feature->seqname   ($ug->seqname);
      $trunc_feature->strand    ($ug->strand);
      $trunc_feature->hseqname  ($ug->hseqname);
      $trunc_feature->score     ($ug->score);
      $trunc_feature->percent_id($ug->percent_id);
      if($start_or_end == 0) {
        $trunc_feature->start     ($ug->start+3);
        $trunc_feature->end       ($ug->end);
      }
      else {
        $trunc_feature->start     ($ug->start);
        $trunc_feature->end       ($ug->end-3);
      }
      $trunc_feature->external_db_id($ug->external_db_id);
      $trunc_feature->hcoverage($ug->hcoverage);
      $trunc_feature->analysis($orignial_analysis) ;

      $trunc_feature->hend ($ug->hend);
      $trunc_feature->hstart($ug->hend -
        ceil($trunc_feature->length / $feature_unit_length) +
        1);
      if ($trunc_feature->end >= $trunc_feature->start) {
        push @ug, $trunc_feature;
      }
    }}

  $new_exon = add_dna_align_features_by_hitname_and_analysis(\@ug,$new_exon,$feature_isa) ;

return $new_exon;
}

=head2 add_dna_align_features_by_hitname_and_analysis

  Arg [1]  : reference to an array of ungapped features
  Arg [1]  : Bio::EnsEMBL::Exon
  Function : to add supporting features to an exon 

=cut

sub add_dna_align_features_by_hitname_and_analysis {
  my ($ug_ref,$exon,$feature_isa) = @_ ;
  my %group_features_by_hitname_and_analysis ;
  for my $ug ( @$ug_ref ) {
    push @{$group_features_by_hitname_and_analysis{$ug->analysis->logic_name}{$ug->hseqname }} , $ug ;
  }
  for my $logic_name ( keys %group_features_by_hitname_and_analysis  ) {
    for my $hseqname  ( keys  %{$group_features_by_hitname_and_analysis{$logic_name}}) {
      my @features = @{$group_features_by_hitname_and_analysis{$logic_name}{$hseqname}};
      my $f;
      if ($feature_isa eq "Bio::EnsEMBL::DnaDnaAlignFeature") {
        $f = Bio::EnsEMBL::DnaDnaAlignFeature->new(-features => \@features, -align_type => 'ensembl');
      } elsif ($feature_isa eq "Bio::EnsEMBL::DnaPepAlignFeature") {
        $f = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@features, -align_type => 'ensembl');
      } else {
      	throw("Cannot create feature. Unknown feature isa.");
      }
      $exon->add_supporting_features($f);
    }
  }
  return $exon ;
}


=head2 remove_initial_or_terminal_short_exons

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : int, minimum length
  Function  : to trim transcripts of initial or terminal exons
  which are consider too short, by default less than 3bp
  Returntype: Bio::EnsEMBL::Transcript
  Exceptions: 
  Example   : 

=cut


sub remove_initial_or_terminal_short_exons{
  my ($transcript, $min_length) = @_;

  $min_length = 3 unless($min_length);
  throw("TranscriptUtils::remove_initial_or_terminal_short_exons will not work ".
        " if ".Transcript_info($transcript)." has no translation") if(!$transcript->translation);
  my %translateable = %{_identify_translateable_exons
                          ($transcript)};

  my $cloned_transcript = clone_Transcript($transcript);
  my $start_exon = $cloned_transcript->translation->
    start_Exon;
  my $start_exon_id = $start_exon->start.":".
    $start_exon->end.":".$start_exon->strand.":".
      $start_exon->phase;
  my $end_exon = $cloned_transcript->translation->end_Exon;
  #If the first exon is shorter than min length
  if($cloned_transcript->start_Exon->length < $min_length){
    my @exons = @{$cloned_transcript->get_all_Exons};
    FIVE_PRIME_END: while(scalar(@exons)){
      my $exon = shift(@exons);
      #throw away each exon which is shorter than the min lenght
      #and recreate the transcript from the point where the
      #exons become longer again
      #This follows pretty much the same princples are translation
      #creation in the spliting of transcripts
      if($exon->length > $min_length){
        my $unique_string = $exon->start.":".$exon->end.":".
          $exon->strand.":".$exon->phase;
        if($translateable{$unique_string} && $exon ne 
           $start_exon){
          $cloned_transcript->translation->start_Exon
            ($exon);
          $cloned_transcript->translation->start(1);
          $exon = _trim_translation_start($exon);
          $cloned_transcript->flush_Exons;
          $cloned_transcript->add_Exon($exon);
          foreach my $e(@exons){
            $cloned_transcript->add_Exon($e);
          }
          last FIVE_PRIME_END;
        }
      }else {
        print id($exon)." is being trimmed  - it is too short\n"; 
        my $warn = id($exon)." is being trimmed ".
          "it is too short\n".Exon_info($exon);
        logger_info($warn);
      }
    }
  }
  #This is as above but it starts from the last exon and works
  #backwards
  if($cloned_transcript->end_Exon->length < $min_length){
    my @exons = @{$cloned_transcript->get_all_Exons};
  THREE_PRIME_END: while(@exons){
      my $exon = pop @exons;
      if($exon->length > $min_length){
        my $unique_string = $exon->start.":".$exon->end.":".
          $exon->strand.":".$exon->phase;
        if($translateable{$unique_string} && $exon !=
           $cloned_transcript->translation->end_Exon){
          $cloned_transcript->translation->end_Exon($exon);
          $exon = _trim_translation_end($exon);
        }
        $cloned_transcript->flush_Exons;
        $cloned_transcript->add_Exon($exon);
        my $translation_end = $exon->end - $exon->start + 1;
        $cloned_transcript->translation->end($translation_end);
        foreach my $e(@exons){
          $cloned_transcript->add_Exon($e);
        }
        last THREE_PRIME_END;
      }else{
        print id($exon)." is being trimmed  - it is too short\n"; 
        my $warn = id($exon)." is being trimmed ".
          "it is too short\n".Exon_info($exon);
        logger_info($warn);
      }
    }
  }
  return $cloned_transcript;
}



=head2 _identify_translateable_exons

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : get a list of exons in the translation
  Returntype: hashref
  Exceptions: 
  Example   : 

=cut


#This method assumes the exons are in the correct 5'-3' order if
#this ever changes this will fail to work without adding a sort

sub _identify_translateable_exons{
  my ($transcript) = @_;

  my $start_exon = 
    $transcript->translation->start_Exon;
  my $end_exon = 
    $transcript->translation->end_Exon;
  my %translateable;
 EXON:foreach my $ex(@{$transcript->get_all_Exons}){
    my $unique_string = $ex->start.":".$ex->end.":".
      $ex->strand.":".$ex->phase;
    if($ex ne $start_exon && !%translateable){
      next EXON;
    }
    if($ex eq $end_exon){
      $translateable{$unique_string} = $ex;
      last EXON;
    }
    $translateable{$unique_string} = $ex;
  }
  return \%translateable;
}





=head2 _trim_translation_start/end

  Arg [1]   : Bio::EnsEMBL::Exon
  Function  : trim untranslated bases from the start or end of
  an exon
  Returntype: Bio::EnsEMBL::Exon 
  Exceptions: 
  Example   : 

=cut


#Should these be here, they are not exported but they take exons not transcripts
#so don't really follow convention, if they are ever made public they should

#This adjusts the appropriate end of an exon to ensure the
#translation has no untranslated basepairs preceeding it


sub _trim_translation_start{
  my ($exon) = @_;
  my $addition = 3 - $exon->phase;
  $addition = 0 if ($exon->phase == 0);
  my $tmp_start;
  if($exon->strand == 1){
    $tmp_start = $exon->start + $addition;
    $exon->start($tmp_start);
  }else{
    $tmp_start = $exon->end - $addition;
    $exon->end($tmp_start);
  }
  $exon->phase(0);
  return $exon;
}

#as above but for translation ends

sub _trim_translation_end{
  my ($exon) = @_;
  my $tmp_end;
  if($exon->strand == 1){
    $tmp_end = $exon->end - $exon->end_phase;
    $exon->end($tmp_end);
  }else{
    $tmp_end = $exon->start + $exon->end_phase;
    $exon->start($tmp_end);
  }
  return $exon;
}



=head2 dump_cDNA_file

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : string, filename
  Arg [3]   : string, format (optional)
  Function  : dump file using Bio::SeqIO
  Returntype: filename
  Exceptions: throws if fails to write or close file
  Example   : 

=cut


sub dump_cDNA_file{
  my ($transcript, $filename, $format) = @_;
  $format = 'fasta' if(!$format);
  logger_info("You are going to dump the cdna of ".
              $transcript->display_id." into ".$filename." format ".
              $format);
  my $seq = $transcript->seq;
  $seq->display_id($transcript->display_id);
  $filename = write_seqfile($seq, $filename, $format);
  return $filename;
}

=head2 convert_to_genes( $tref, $analysis )

Arg[0]     : Arrayref of Bio::EnsEMBL::Transcript objects
Arg[1]     : Bio::EnsEMBL::Analysis object (opt)
Name       : convert_to_genes
Function   : converts all transcripts to an arrayref of single transcript genes (biotype of transcript becomes biotype of gene)
Returntype : Arrayref of Bio::EnsEMBL::Gene objects

=cut

sub convert_to_genes {
  my ($tref, $analysis, $biotype ) = @_;
  my @genes ;

  my @tr = (ref($tref)=~m/ARRAY/) ? @$tref : ($tref) ; 

  for my $t (@tr) {
    throw("Problems with ".Transcript_info($t)." undef coords") if(!$t->start || !$t->end);
    fully_load_Transcript($t);
    $analysis = $t->analysis unless $analysis ;
    if($biotype){ $t->biotype($biotype); }
    my $g = Bio::EnsEMBL::Gene->new() ;
    $g->biotype( $t->biotype ) ;
    $g->add_Transcript($t);
    $g->analysis($analysis) ;
    $g->dbID($t->dbID);
    push @genes, $g ;
  }
  for my $g(@genes) {
    throw("there are no transcripts for this gene\n") if scalar(@{$g->get_all_Transcripts}) == 0 ;
    for my $tr ( @{$g->get_all_Transcripts} ) {
      throw("there are no exons  for this transcript \n") if scalar(@{$tr->get_all_Exons}) == 0 ;
    }
  }
  foreach my $gene(@genes){
    throw("Problems with ".Gene_info($gene)." undef coords") if(!$gene->start || !$gene->end);
  }
  return \@genes ;
}


=head2 evidence_coverage

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : Bio::Seq 
  Function  : 
  Returntype: 
  Exceptions: 
  Example   : 

=cut

#note this method is designed to work with protein evidence, it should work with
#dna evidence but I am not sure

sub evidence_coverage{
  my ($transcript, $evidence) = @_;

  throw("Can't work out evidence coverage for ".$transcript->display_id.
        " without any evidence") if(!$evidence);
  my $matches = 0;
  my $pstart  = 0;
  my $pend    = 0;
  my $evidence_name = $evidence->id;
  my $plength;

  foreach my $exon(@{$transcript->get_all_Exons}) {
    $pstart = 0;
    $pend   = 0;
    
    foreach my $f(@{$exon->get_all_supporting_features}){
      
      if($evidence_name ne $f->hseqname){
        warning("$evidence_name ne " . $f->hseqname . "\n");
      }
      
      if((!$pstart) || $pstart > $f->hstart){
        $pstart = $f->hstart;
      }
      
      if((!$pend) || $pend < $f->hend){
        $pend= $f->hend;
      }
    }
    $matches += ($pend - $pstart + 1);
  }

  $plength = $evidence->length;
  if(!defined($plength) || $plength == 0){
    warning("TranscriptUtils: no sensible length for ".$evidence_name." assuming 0% ".
            "coverage");
    return 0;
  }
  my $coverage = int(100 * $matches/$plength);
  foreach my $sf(@{$transcript->get_all_supporting_features}){
    $sf->hcoverage($coverage) if($sf->hseqname eq $evidence_name);
  }
  return $coverage;
}


sub get_upstream_Intron_from_Exon{
  my ($transcript, $exon) = @_;
  return get_upstream_Intron($exon, $transcript);
}

sub get_upstream_splice_sites_from_Exon{
  my ($transcript, $exon) = @_;
  return get_upstream_splice_sites($exon, $transcript);
}

sub get_downstream_Intron_from_Exon{
  my ($transcript, $exon) = @_;
  return get_downstream_Intron($exon, $transcript);
}

sub get_downstream_splice_sites_from_Exon{
  my ($transcript, $exon) = @_;
  return get_downstream_splice_sites($exon, $transcript);
}


sub attach_Slice_to_Transcript{
  my ($transcript, $slice) = @_;
  throw('You need a Bio::EnsEMBL::Slice object not a "'.ref($slice).'"')
    unless ($slice and ref($slice) eq 'Bio::EnsEMBL::Slice');
  $transcript->slice($slice);
  foreach my $sf(@{$transcript->get_all_supporting_features}){
    $sf->slice($slice);
  }

  foreach my $ise (@{$transcript->get_all_IntronSupportingEvidence}){
    $ise->slice($slice);
  }

  foreach my $exon(@{$transcript->get_all_Exons}){
    $exon->slice($slice);
    foreach my $sf(@{$exon->get_all_supporting_features}){
      $sf->slice($slice);
    }
  }
}

sub attach_Analysis_to_Transcript{
  my ($transcript, $analysis) = @_;
  throw('You need a Bio::EnsEMBL::Analysis object not a "'.ref($analysis).'"')
    unless ($analysis and ref($analysis) eq 'Bio::EnsEMBL::Analysis');
  $transcript->analysis($analysis);
  foreach my $sf(@{$transcript->get_all_supporting_features}){
    $sf->analysis($analysis);
  }
  foreach my $sf (@{$transcript->get_all_IntronSupportingEvidence}){
    $sf->analysis($analysis);
  }
  foreach my $exon(@{$transcript->get_all_Exons}){
    $exon->analysis($analysis);
    foreach my $sf(@{$exon->get_all_supporting_features}){
      $sf->analysis($analysis);
    }
  }
}
sub attach_Analysis_to_Transcript_no_support{
  my ($transcript, $analysis) = @_;
  $transcript->analysis($analysis);
}

=head2 attach_Analysis_to_Transcript_no_overwrite

 Arg [1]    : Bio::EnsEMBL::Transcript
 Arg [2]    : Bio::EnsEMBL::Analysis
 Description: Attach Arg[2] to Arg[1] and all supporting evidence and exons
              unless the object already has an analysis
              This is usefull if your databases analysis table are not synchronised
 Returntype : None
 Exceptions : Throws if Arg[2] is not a Bio::EnsEMBL::Analysis

=cut

sub attach_Analysis_to_Transcript_no_overwrite {
  my ($transcript, $analysis) = @_;

  throw('You need a Bio::EnsEMBL::Analysis object not a "'.ref($analysis).'"')
    unless ($analysis and ref($analysis) eq 'Bio::EnsEMBL::Analysis');
  $transcript->analysis($analysis) unless ($transcript->analysis);
  foreach my $sf (@{$transcript->get_all_supporting_features}) {
    $sf->analysis($analysis) unless ($sf->analysis);
  }
  foreach my $sf (@{$transcript->get_all_IntronSupportingEvidence}){
    $sf->analysis($analysis) unless ($sf->analysis);
  }
  foreach my $exon (@{$transcript->get_all_Exons}) {
    $exon->analysis($analysis) unless ($exon->analysis);
    foreach my $sf (@{$exon->get_all_supporting_features}) {
      $sf->analysis($analysis) unless ($sf->analysis);
    }
  }
}

sub fully_load_Transcript{
  my ($transcript, $keep_xrefs) = @_;

  $keep_xrefs = 1 if(!defined($keep_xrefs));

  if ($transcript->translation) {
    $transcript->get_all_alternative_translations();
    $transcript->translate;
    $transcript->translation->get_all_Attributes;
    $transcript->translation->get_all_DBEntries if($keep_xrefs);
    $transcript->translation->get_all_ProteinFeatures;
    $transcript->translation->get_all_SeqEdits;
  }

 EXON:foreach my $e(@{$transcript->get_all_Exons}){
    $e->analysis;
    $e->stable_id;
    $e->get_all_supporting_features;
  }

  $transcript->stable_id;
  $transcript->analysis;
  $transcript->get_all_Attributes;
  $transcript->get_all_DBEntries if($keep_xrefs);
  $transcript->get_all_supporting_features;
  $transcript->get_all_alternative_translations();
  $transcript->get_all_SeqEdits();
  $transcript->get_all_IntronSupportingEvidence();

  return $transcript;
}

sub empty_Transcript{
  my ($transcript, $remove_stable_id, $remove_xrefs) = @_;
  fully_load_Transcript($transcript);
  foreach my $sf(@{$transcript->get_all_supporting_features}){
    empty_Object($sf);
  }
  foreach my $ise(@{$transcript->get_all_IntronSupportingEvidence}){
    empty_Object($ise);
  }
  empty_Object($transcript->translation, $remove_stable_id) 
    if($transcript->translation);;
 EXON:foreach my $e(@{$transcript->get_all_Exons}){
  SF:foreach my $sf(@{$e->get_all_supporting_features}){
      empty_Object($sf, $remove_stable_id);
    }
    empty_Object($e, $remove_stable_id);
  }
  if($remove_xrefs) {
    $transcript->display_xref(undef);
    # It is naughty to go into the transcript & translation object
    # but we have no API method to do this:
    $transcript->{'dbentries'} = [];
    if ($transcript->translation) {
      $transcript->translation->{'dbentries'} = [];
    }
  }
  empty_Object($transcript, $remove_stable_id);
  return $transcript;
}

sub all_exons_are_valid{
  my ($transcript, $max_length, $allow_negative_start) = @_;

  foreach my $exon(@{$transcript->get_all_Exons}){
    throw(Transcript_info($transcript)." seems to contain an undefined exon") 
      if(!$exon); 
    if(!exon_length_less_than_maximum($exon, $max_length)){
      warn("Transcript ".Transcript_info($transcript)." has ".
                  "exon longer than ".$max_length);
      return 0;
    }
    if(!validate_Exon_coords($exon, $allow_negative_start)){
      warn("Transcript ".Transcript_info($transcript)." has ".
                  "invalid exon coords ".Exon_info($exon));
      return 0;
    }
  }
  return 1;
}


sub evidence_coverage_greater_than_minimum{
  my ($transcript, $evidence, $min_coverage) = @_;
  if ($evidence) {
    my $coverage = evidence_coverage($transcript, $evidence);
    if (!($coverage > $min_coverage)) {
      warning(Transcript_info($transcript)." based on ".$evidence->id." has too ".
            "low coverage ".$coverage." of the evidence");
    return 0;
    }
  }
  else {
    warning ("There has been no evidence passed in for ".Transcript_info($transcript).
             "Assuming coverage is fine");
  }
  return 1;
}


=head2 identical_Transcripts

 Title   : identical_Transcripts
 Usage   : $identical = identical_Transcripts($transcript1, $transcript2);
 Function: compares 2 Transcripts. DOES NOT CHECK TO SEE WHETHER PHASES ARE IDENTICAL
           OR WHETHER EVIDENCE IS IDENTICAL.
           Transcripts are compared on the basis of (i) number of Exons,
           (ii) start and end coordinates for each Transcript, (iii) strand
 Example :
 Returns : 1 if identical, 0 if not indentical
 Args    : Transcript, Transcript

=cut

sub identical_Transcripts {
  my ($transcript1, $transcript2) = @_;

  my @exons1 = @{$transcript1->get_all_Exons()};
  my @exons2 = @{$transcript2->get_all_Exons()};

  # compare no. of Exons
  if (scalar(@exons1) != scalar(@exons2)) {
    return 0;
  }
  # compare Exon coordinates and strand
  foreach ( my $num = 0; $num < scalar(@exons1); $num++ ) {

    if ($exons1[$num]->strand != $exons2[$num]->strand) {
      return 0;
    }
    if ($exons1[$num]->start != $exons2[$num]->start) {
      return 0;
    }
    if ($exons1[$num]->end != $exons2[$num]->end) {
      return 0;
    }
  }

  # you may want to check the evidence and phase at some stage by adding a wrapper method
  return 1;
}

sub set_start_codon{
  my ($transcript) = @_;

  # check transcript has a translation
  if(!$transcript->translation || !$transcript->translation->start_Exon){
    warning("Transcript has no translation, or no start exon - maybe a pseudogene?");
    return $transcript;
  }
  my $cloned_transcript = clone_Transcript($transcript);
  # useful info in genomic coordinates
  my $strand = @{$cloned_transcript->get_all_Exons}[0]->strand;
  my $translation       = $cloned_transcript->translation;
  my $start_exon        = $translation->start_Exon;
  my $cdna_coding_start = $cloned_transcript->cdna_coding_start;
  my $cdna_seq          = uc($cloned_transcript->spliced_seq);
  my @pepgencoords      = $cloned_transcript->pep2genomic(1,1);

  if(scalar(@pepgencoords) > 2) {
    logger_info("peptide start does not map cleanly - not modifying transcript");
    return $cloned_transcript;
  }

  my $pepgenstart = $pepgencoords[0]->start;
  my $pepgenend   = $pepgencoords[$#pepgencoords]->end;

  unless($pepgencoords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate') && 
	 $pepgencoords[$#pepgencoords]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
    logger_info("peptide coordinate(s)  maps to gap - not modifying transcript");
    return $cloned_transcript;
  }

  ############################################################
  # first see whether the transcript already begins with ATG
  my $first_codon = substr($cdna_seq, $cdna_coding_start-1, 3);

  if ( uc($first_codon) eq 'ATG' ){
    logger_info("transcript already starts with ATG - no need to modify");
    return $cloned_transcript;
  }
  
  ############################################################
  # now look at the previous codon
  ############################################################
  # first the simplest cases
  if($cdna_coding_start>3){
    # the previous codon is in the cdna
    $first_codon = substr($cdna_seq, $cdna_coding_start-4, 3);
    if ($first_codon ne 'ATG'){
      logger_info("Upstream codon is not an ATG - not modifying transcript");
      return $cloned_transcript;
    }else{
      # save current coords, just in case we need to revert
      my $current_translation_start = $cloned_transcript->translation->start;
      my $current_start_exon        = $cloned_transcript->translation->start_Exon;
      my $current_start_exon_start  = $current_start_exon->start;
      my $current_start_exon_end  = $current_start_exon->end;
      my $current_start_exon_phase = $current_start_exon->phase;
      my $newstartexon;
      my $current_newstartexon_endphase;

      my @coords = $cloned_transcript->cdna2genomic($cdna_coding_start-3,$cdna_coding_start-1,$strand);
      my $new_start;
      my $new_end;

      # check not mapping to gaps
      unless($coords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate') &&
	     $coords[$#coords]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
	logger_info("new coordinate(s) maps to gap - not modifying transcript");
	return $cloned_transcript;
      }

      if (scalar(@coords) > 2){
	print STDERR "coordinate mapping not done cleanly - not modifying transcript!\n";
	return $cloned_transcript;
      }elsif(scalar(@coords) == 2){
	logger_info("new start codon split across intron");
	logger_info("coord[0] = " . $coords[0]->start . " " . $coords[0]->end);
	logger_info("coord[1] = " . $coords[1]->start . " " . $coords[1]->end);
	
	if($strand == 1){
	  $new_start = $coords[0]->start;
	  $new_end   = $coords[$#coords]->end;}
	else{
	  $new_start = $coords[0]->end;
	  $new_end   = $coords[$#coords]->start; 
	}
	
	# find exon
	my $newstartexon = get_previous_Exon($cloned_transcript, $start_exon);

	if (!$newstartexon) {
	  logger_info("Failed finding new start exon - not modifying transcript");
	  return $cloned_transcript;
	}

	# save in case we need to revert
	my $current_newstartexon_endphase = $newstartexon->end_phase;
	
	my $newphase;
	if ($strand == 1) {
	  $newphase = $newstartexon->end - $new_start + 1;
	} else {
	  $newphase = $new_start - $newstartexon->start + 1;
	}
	
	$start_exon->phase($newphase);
	$newstartexon->end_phase($newphase);
	$translation->start_Exon($newstartexon);
	$translation->start($newstartexon->length-$newphase+1);

	# make sure it still translates, and revert if necessary
	eval{
	  $cloned_transcript->translate;
	};
	if($@){
	  logger_info("problem with modified transcript - reverting coordinates");
	  $cloned_transcript->start_Exon($current_start_exon);
	  $cloned_transcript->start_Exon->start($current_start_exon_start);
	  $cloned_transcript->start_Exon->end($current_start_exon_end);
	  $translation->start($current_translation_start);
	  $cloned_transcript->start_Exon->phase($current_start_exon_phase);
	  if ($newstartexon){
	    $newstartexon->end_phase($current_newstartexon_endphase);
	  }
	}	
	
	$cloned_transcript->recalculate_coordinates;
	return $cloned_transcript;
      }else{
	logger_info("New start codon doesn't split across introns - but which exon is it in");
	$new_start = $coords[0]->start;
	$new_end   = $coords[0]->end;
	if (($strand == 1  && $new_end == $pepgenstart-1) ||
	    ($strand == -1 && $new_start == $pepgenend+1)) { 
	  $translation->start($translation->start-3);
	} else{
          # find exon
	  my $newstartexon = get_previous_Exon($cloned_transcript, $start_exon);
	  if (!$newstartexon) {
	    logger_info("Failed finding new start exon - how can this be?");
	    return $cloned_transcript;
	  }

	  $current_newstartexon_endphase = $newstartexon->end_phase;
           
	  # make the boundary phases 0 - the ATG is the last codon of $newstartexon 
	  # as we know it doesn't cross the intron
	  $start_exon->phase(0);
	  $newstartexon->end_phase(0);

	  # Reset translation start exon
	  $translation->start_Exon($newstartexon);
	  $translation->start($newstartexon->length-2);
	}
	
	# make sure it still translates, and revert if necessary
	eval{
	  $cloned_transcript->translate;
	};

	if($@){
          logger_info("problem with modified transcript - reverting coordinates");
	  $cloned_transcript->start_Exon($current_start_exon);
	  $cloned_transcript->start_Exon->start($current_start_exon_start);
	  $cloned_transcript->start_Exon->end($current_start_exon_end);
	  $translation->start($current_translation_start);
	  $cloned_transcript->start_Exon->phase($current_start_exon_phase);
	  $newstartexon->end_phase($current_newstartexon_endphase);
	}
	$cloned_transcript->recalculate_coordinates;
	return $cloned_transcript;
      } 
    }
  }

  ############################################################
  # more complex cases: the previous codon falls off the cdna
  else{
    my $codon_start;
    my $codon_end; 
    
    if ($strand == 1) {
      $codon_start = $pepgenstart - 3;
      $codon_end   = $pepgenstart - 1;
    } else {
      $codon_start = $pepgenend + 1;
      $codon_end   = $pepgenend + 3;
    }
    
    my $seq_adaptor = $start_exon->slice->adaptor->db->get_SequenceAdaptor;
    my $codonseq      = uc(${$seq_adaptor->fetch_by_Slice_start_end_strand
                           ($start_exon->slice, $codon_start,$codon_end, 
                            $strand)});
    
    if ($codonseq ne "ATG") {
      logger_info("upstream codon (faling off the slice) is not ATG - not modifying transcript");
      return $cloned_transcript;
    }
    else{
      # fun fun fun

      # save current coordinates in case we need to revert
      my $current_start_exon         = $start_exon;
      my $current_start_exon_start    = $start_exon->start;
      my $current_start_exon_end      = $start_exon->end;
      my $current_start_exon_phase    = $start_exon->phase;
      my $current_start_exon_endphase = $start_exon->end_phase;
      my $current_translation_start  = $translation->start;
      my $current_translation_end    = $translation->end;

      if($strand == 1){
	$start_exon->start($codon_start)
      }
      else{
	$start_exon->end($codon_end)
      }
      $start_exon->phase(0);
      if ($translation->end_Exon == $start_exon){
	$translation->end($translation->end + (4-$translation->start));
      }
      $translation->start(1);

      # make sure it still translates, and revert if necessary
      eval{
	$cloned_transcript->translate;
      };
      if($@){
	logger_info("problem with modified transcript - reverting coordinates");
	$cloned_transcript->start_Exon($current_start_exon);
	$cloned_transcript->start_Exon->start($current_start_exon_start);
	$cloned_transcript->start_Exon->end($current_start_exon_end);
	$translation->start($current_translation_start);
	$translation->end($current_translation_end);
	$cloned_transcript->start_Exon->phase($current_start_exon_phase);
	$cloned_transcript->start_Exon->end_phase($current_start_exon_endphase);
      }

      $cloned_transcript->recalculate_coordinates;
	return $cloned_transcript;
    }
  }
}

sub get_previous_Exon{
  my ($transcript, $exon ) = @_;
    
  # this order the exons 5' to 3'
  
  my @exons = @{$transcript->get_all_Exons};
  
  for (my $i=0; $i<=$#exons; $i++ ){
    if ( $exons[$i]->start == $exon->start 
	 && 
	 $exons[$i]->end   == $exon->end
	 &&
	 $exons[$i]->strand == $exon->strand 
	 &&
	 $i > 0 
       ){
      return $exons[$i-1];
    }
  }
  return undef;
}

sub get_next_Exon{
  my ($transcript, $exon ) = @_;
    
  # this order the exons 5' to 3'

  my @exons = @{$transcript->get_all_Exons};
  for (my $i=0; $i<=$#exons; $i++ ){
    if ( $exons[$i]->start == $exon->start 
	 && 
	 $exons[$i]->end   == $exon->end
	 &&
	 $exons[$i]->strand == $exon->strand 
	 &&
	 ($i+1) <= $#exons
       ){
      return $exons[$i+1];
    }
  }
  return undef;
}

sub set_stop_codon{
  my ($transcript) = @_;
  my  $verbose = 0;
  unless ( $transcript->translation ){
    logger_info("transcript has no translation - cannot put the stops");
    return $transcript;
  }
  my $cloned_transcript = clone_Transcript($transcript);
  
  my $end      = $cloned_transcript->translation->end;
  my $end_exon = $cloned_transcript->translation->end_Exon;
  
  ############################################################
  # first see whether the transcript already include the stop:  taa/tag/tga
  # this gives you the sequence 5' to 3'

  my $bioseq = $end_exon->seq; 
  my $last_codon;
  if ( $end > 2 ){
    $last_codon = $bioseq->subseq( $end - 2, $end );
  }else{
    my $donor    = 3 - $end;
    my $acceptor = $end;

    my $previous_exon = get_previous_Exon( $cloned_transcript, $end_exon );
    if ($previous_exon ){
      my $subseq_start = 
        $previous_exon->end - $previous_exon->start + 1 - $donor + 1;
      my $subseq_end =  $previous_exon->end - $previous_exon->start + 1;
      my $donor_seq = $previous_exon->seq->subseq
        ($subseq_start, $subseq_end);
      my $acceptor_seq = $end_exon->seq->subseq( 1, $end );
      $last_codon = $donor_seq.$acceptor_seq;
    }
  }
  if ( uc($last_codon) eq 'TAA' || uc($last_codon) eq 'TAG' 
       || uc($last_codon) eq 'TGA' ){ 
    logger_info("transcript already has a stop at the end - no need ".
                "to modify");
    return $cloned_transcript;
  }
  ############################################################
  # now look at the next codon
  ############################################################
  # first the simplest case
  if ( $end + 3 <= ($end_exon->end - $end_exon->start + 1) ){
    my $next_codon = $bioseq->subseq( $end+1, $end+3 );      
    if ( uc($next_codon) eq 'TAA' || uc($next_codon) eq 'TAG' || 
         uc($next_codon) eq 'TGA'){ 
      logger_info("simple-case: next codon is a stop - ".
                  "extending translation\n");
      $cloned_transcript->translation->end( $end + 3 );
      $cloned_transcript->recalculate_coordinates;
      return $cloned_transcript;
    }else{
      logger_info("next codon is not a stop - not modifying translation");
      return $cloned_transcript;
    }
  }
  
  ############################################################
  # more complex cases we need to know if there is a next exon:
  my $next_exon = get_next_Exon( $cloned_transcript, $end_exon );
  
  if ( $next_exon ){
    ############################################################
    # how many bases of the next codon sit in $end_exon?
    my $donor_bases_count = ( $end_exon->end - $end_exon->start + 1 ) 
      - $end;
    my $acceptor_bases_count = 3 - $donor_bases_count;
    
    ############################################################
    # get the next codon
    my $next_bioseq = $next_exon->seq;
    my $donor;
    if ( $donor_bases_count == 0 ){
      $donor = '';
    }
    else{
      $donor = $bioseq->subseq( $end+1, ($end_exon->end - 
                                         $end_exon->start + 1 ));
    }
    my $acceptor = $next_bioseq->subseq( 1, $acceptor_bases_count );

    my $next_codon = $donor.$acceptor;
    if ( uc($next_codon) eq 'TAA' || uc($next_codon) eq 'TAG' 
         || uc($next_codon) eq 'TGA'){ 
      logger_info("shared-codon: next codon is a stop - ".
                  "extending translation");

      $cloned_transcript->translation->end_Exon( $next_exon );
      $cloned_transcript->translation->end( $acceptor_bases_count );
      ############################################################
      # re-set the phases:
      $end_exon->end_phase($donor_bases_count%3);
      $next_exon->phase( $donor_bases_count%3 );
      $cloned_transcript->recalculate_coordinates;
      return $cloned_transcript;
    } else{
      logger_info("next codon is not a stop - not modifying translation");
      return $cloned_transcript;
    }
  }elsif( $end + 3 > ($end_exon->end - $end_exon->start + 1) ){
    # there is no next exon and the next codon would fall off the end 
    # of the exon 

    # need to get the slice sequence
    my $adaptor =  $end_exon->slice->adaptor;
    if ( $adaptor ){
      my $donor_bases_count = ( $end_exon->end - $end_exon->start + 1 ) 
        - $end;
      my $acceptor_bases_count = 3 - $donor_bases_count;

      # the sequence from the current end exon is:
      my $donor;
      if ( $donor_bases_count == 0 ){
        $donor = '';
      } else {
        $donor = $bioseq->subseq( $end+1, ( $end_exon->end - 
                                            $end_exon->start + 1 ));
      }

      ############################################################
      # here we distinguish the strands
      if ( $end_exon->strand == 1 ){
        my $slice_start = $end_exon->slice->start;

        ############################################################
        # calculate the next codon start/end in chr coordinates 

        my $codon_start = $slice_start + ( $end_exon->start + $end - 1 );
        my $codon_end   = $codon_start + 2;

        my $codon_slice = $adaptor->fetch_by_region
          ($end_exon->slice->coord_system->name, 
           $end_exon->slice->seq_region_name, $codon_start, $codon_end );
        my $codon = $codon_slice->seq;

        ############################################################
        if ( uc($codon) eq 'TAA' || uc($codon) eq 'TAG' || uc($codon) eq 'TGA'){ 
          logger_info("forward-strand:next codon (falling off the exon) ".
                      "is a stop - extending translation");
          $end_exon->end( $end_exon->end + $acceptor_bases_count );
          $cloned_transcript->translation->end( $end + 3 );
          
          ############################################################
          # update the exon sequence:	    	    
          my $seq_string = $end_exon->slice->subseq( $end_exon->start, $end_exon->end, $end_exon->strand );
      
          $cloned_transcript->translation->end_Exon($end_exon);
      
          $cloned_transcript->recalculate_coordinates;
          return $cloned_transcript;
        }
        else{
          logger_info("next codon (falling off the exon) is not a stop - not modifying");
          return $cloned_transcript;
        }
      } else {
        my $slice_start = $end_exon->slice->start;
        
        ############################################################
        # calculate the next codon start/end in chr coordinates 

        my $codon_end   = $slice_start + $end_exon->end - $end - 1;
        my $codon_start = $codon_end - 2;
       
        if($codon_start <= 0){
          logger_info("Can't extend the transcript off the end of a ".
                  $end_exon->slice->coord_system->name);
          return $cloned_transcript;
        }
        my $codon_slice = $adaptor->fetch_by_region
          ($end_exon->slice->coord_system->name, 
           $end_exon->slice->seq_region_name, $codon_start, $codon_end );
        my $pre_codon = $codon_slice->seq;
        
        # need to reverse and complement:
        my $codon;
        ( $codon = reverse $pre_codon ) =~tr/gatcGATC/ctagCTAG/; 
       ;
        if ( uc($codon) eq 'TAA' || uc($codon) eq 'TAG' || uc($codon) eq 'TGA'){ 
          logger_info("reverse-strand: next codon (falling off the exon) is a stop - extending translation\n");

          $end_exon->start( $end_exon->start - $acceptor_bases_count);
          $cloned_transcript->translation->end( $end + 3 );

          ############################################################
          # update the exon sequence:	    	    
          my $seq_string = $end_exon->slice->subseq( $end_exon->start, $end_exon->end, $end_exon->strand );

          $cloned_transcript->translation->end_Exon( $end_exon );
          $cloned_transcript->recalculate_coordinates;
          return $cloned_transcript;
        } else {
          logger_info("next codon (falling off the exon) is not a stop - not modifying");
          return $cloned_transcript;
        }
      }
    } else {
      logger_info("cannot get an adaptor to get the sequence - not modifying the translation");
      return $cloned_transcript;
    }
  } else {
    logger_info("There is no downstream exon - and no stop codon beyond the last exon - not modifying");
    return $cloned_transcript;
  }
}




=head2 is_Transcript_sane

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : checks some simple facts about a transcript structure to ensure
  its sane
  Returntype: boolean 
  Exceptions: none
  Example   : 

=cut


sub is_Transcript_sane {
  my ($transcript) = @_;

  # Exon coord sanity.
  my $is_sane = 1;

  foreach my $exon ( @{ $transcript->get_all_Exons() } ) {
    if ( $exon->start() > $exon->end() ) {
      $is_sane = 0;
      last;
    }
  }

  if ( $is_sane && !are_strands_consistent($transcript) ) {
    $is_sane = 0;
  }
  if ( $is_sane && !are_phases_consistent($transcript) ) {
    $is_sane = 0;
  }
  if ( $is_sane && !is_not_folded($transcript) ) {
    $is_sane = 0;
  }

  return $is_sane;
} ## end sub is_Transcript_sane


sub features_overlap {
# return 1 if there featureA overlaps feature B
# return 0 otherwise
  my ($featureA,$featureB) = @_;

  if (($featureA->seq_region_start() <= $featureB->seq_region_end()) and
      ($featureA->seq_region_end() >= $featureB->seq_region_start())) {
    return 1;
  }
  return 0;
}

sub overlap_length {
# return the length of the overlap between featureA and featureB 
  my ($featureA,$featureB) = @_;

  if (!features_overlap($featureA,$featureB)) {
    return 0;
  }

  my $min_end = $featureA->seq_region_end();
  if ($featureB->seq_region_end() < $min_end) {
    $min_end = $featureB->seq_region_end();
  }

  my $max_start = $featureA->seq_region_start();
  if ($featureB->seq_region_start() > $max_start) {
    $max_start = $featureB->seq_region_start();
  }
  return $min_end-$max_start+1;
}

sub exon_overlap {
# return the length of the exonic overlap between
# transcriptA and transcriptB
  my ($trancriptA,$trancriptB) = @_;

  my $overlap = 0;

  foreach my $exonA (@{$trancriptA->get_all_Exons()}) {
    foreach my $exonB (@{$trancriptB->get_all_Exons()}) {
      $overlap += overlap_length($exonA,$exonB);
    }
  }
  return $overlap;
}

sub coding_exon_overlap {
# return the length of the coding exonic overlap between
# transcriptA and transcriptB
  my ($trancriptA,$trancriptB) = @_;

  my $overlap = 0;

  foreach my $exonA (@{$trancriptA->get_all_translateable_Exons()}) {
    foreach my $exonB (@{$trancriptB->get_all_translateable_Exons()}) {
      $overlap += overlap_length($exonA,$exonB);
    }
  }
  return $overlap;
}

sub remove_short_frameshift_introns {
# return the transcript after removing the short
# frameshift introns found ie merging the
# exons separated by the short frameshift introns
  my $transcript = shift;
  my $transcript_no_frameshift = clone_Transcript($transcript);
  
  $transcript_no_frameshift->flush_Exons();
    
  my @frameshift_attributes = @{$transcript->get_all_Attributes('Frameshift')};
  
  # get all the short frameshift intron numbers
  my @intron_numbers = ();
  foreach my $frameshift_attribute (@frameshift_attributes) {
  	push(@intron_numbers,$frameshift_attribute->value());
  }
  	
  # make sure they are sorted
  my @sorted_intron_numbers = sort {$a <=> $b} @intron_numbers;
  	
  my $removed_introns = 0;
  my $current_exon;
  my $next_exon;
  my $current_exon_index = 1;
  my $merged_exon;
  my @exons = @{$transcript->get_all_Exons()};
  
  SHORT_INTRON: foreach my $intron_number (@sorted_intron_numbers) {

    while ($current_exon = shift(@exons)) {
    
      if ($current_exon_index < $intron_number-$removed_introns) {
      	
      	$transcript_no_frameshift->add_Exon(clone_Exon($current_exon));
      	$current_exon_index++;

      } else {

        $next_exon = shift(@exons);
        $merged_exon = merge_exons($current_exon,$next_exon);	

      	$removed_introns++;
      	      	
      	# the new merged exon has to be processed in case there are
      	# consecutive short frameshift introns so the current_exon_index is not
      	# changed and the new merged exon is added to the exons array 
      	unshift(@exons,$merged_exon);
      	
      	next SHORT_INTRON;
      	
      } # else
    } # while current exon 
  } # foreach short intron
  
  # add the remaining exons
  while ($current_exon = pop(@exons)) {
  	$transcript_no_frameshift->add_Exon(clone_Exon($current_exon));
  }
  return $transcript_no_frameshift;
}


=head2 has_polyA_signal

 Arg [1]    : Bio::EnsEMBL::Transcript object
 Arg [2]    : Boolean $lenient (optional), DNA sequence
 Description: Checks for the presence of AATAAA between 10 and 30nt from the end
              When Arg[2] is set to true, it checks for A[AG]TAAA between 0 and 30nt
              from the end
 Returntype : Boolean true when found, false when not found
 Exceptions : Throws if Arg[1] is not a Bio::EnsEMBL::Transcript object

=cut

sub has_polyA_signal {
  my ($transcript, $lenient) = @_;

  throw('You need to pass a Bio::EnsEMBL::Transcript object not a "'.ref($transcript).'"')
    unless ($transcript->isa('Bio::EnsEMBL::Transcript'));

  my $regex = 'A{2}TA{3}.{10,30}$';
  $regex =~ 'A[AG]TA{3}.{0,30}$' if ($lenient);

  return $transcript->end_Exon->seq->seq =~ /$regex/;
}

=head2 set_alignment_supporting_features

 Arg [1]    : Bio::EnsEMBL::Transcript object
 Arg [2]    : string. Alignment query protein sequence.
 Arg [3]    : string. Alignment target protein sequence.
 Description: 
  This sets supporting features for all exons. It does this by comparing all non-split codons in the exons to their corresponding positions in the alignment.
  For each non-split codon, it makes a feature pair where the start and end are the codon start and end and the
  hit start and end is the position of the corresponding query sequence amino acid in the UNGAPPED version of the
  query sequence. So it 1) Finds the codon in the alignment 2) Finds the index of the corresponding query amino acid in the unaligned query.
  Once codons that have support have had the corresponding hit start and end worked out, it then joins then together
  It does this in a simple way, as the proto feature pairs are already in codon by codon order it just looks at consecutive ones and checks that they are from adjacent codons. If they are it then checks if the hit end of the first is directly beside the hit start of the second. If they are then the left proto-feature pair is deleted and the right one is extended to encompass the left one. This repeats for as long as possible to build maximum lenght feature pairs
  In the cases where there is a gap in the query, there will be no proto SF for the corresponding codon.
  In the case where there is a gap in the target (genome), the hit start and hit ends will not be side by side and thus not joined together.
  By following this an exon might have several feature pairs, some might be side by side and unjoined, some might be on either side of codons with no feature pairs. Often there is just one long feature pair covering all non-split codons in the exon.
  All of the feature pairs on an exon are then passing into a single object (DnaPepAlignFeature). By doing this the API will automatically work out the cigar strings for each exon. By looking at gaps between the feature pairs and at feature pairs that are adjacent but the hit end/start are not adjacent it will model the indels across the exon.
  For the transcript supporting feature you need to add every individual feature pair into an array ref and add to a DnaPepAlignFeature.
 Returntype : N/A
 Exceptions : Throws if Arg[1] is not a Bio::EnsEMBL::Transcript object

=cut

sub set_alignment_supporting_features {
  my ($transcript,$query_seq,$target_seq) = @_;

  if (!($transcript->isa('Bio::EnsEMBL::Transcript'))) {
    throw('You need to pass a Bio::EnsEMBL::Transcript object not a "'.ref($transcript).'"')
  }

  # This will store all feature pairs in the end so that they can be added as a transcript supporting feature
  my $all_exon_supporting_features = [];
  
  # Keep track of the previous feature pair so features not in order can be discarded
  my $prev_feature_pair;

  my $coverage;
  my $percent_id;
  ($query_seq,$target_seq,$coverage,$percent_id) = align_proteins_with_alignment($query_seq,$target_seq);

  # Add a stop to the alignment seqs. Basically this will allow me to ignore the final codon (which is a stop)
  $query_seq .= '*';
  $target_seq .= '*';

  say "";
  say "------------------------------------------------";
  say "NEW TRANSCRIPT";
  say "------------------------------------------------";

  my $codon_index= 0;
  my $exons = $transcript->get_all_Exons();
  my $i=0;
  for($i=0; $i<scalar(@{$exons}); $i++) {

    my $proto_supporting_features = [];
    my $exon = $$exons[$i];
    my $exon_seq = $exon->seq->seq();
    my @nucleotide_array = split('',$exon_seq);    

    my $phase = $exon->phase();
    my $end_phase = $exon->end_phase();
    my $start_index = 0;

    say "";
    say "------------------------------------------------";
    say "NEW EXON";
    say "------------------------------------------------";
    say "FM2 ESTART: ".$exon->start;
    say "FM2 EEND: ".$exon->end;
    say "FM2 ESTRAND: ".$exon->strand;
    say "FM2 EPHASE: ".$exon->phase;
    say "FM2 EENDPHASE: ".$exon->end_phase;

    # If the phase is not 0 then the first codon is a split one. The phase is then number of bases missing from
    # the codon. so if you take the phase from 3 you get the number of bases in the split codon
    if($phase) {
      $start_index += (3 - $phase);
    }

    for(my $k=$start_index; $k<scalar(@nucleotide_array); $k+=3) {

     # Ending on a split codon, so increase the index and finish the loop
     if($k+2 >= scalar(@nucleotide_array)) {
        $codon_index++;        
        say "Skipping split codon";

        last;
      }

      my $target_start;
      my $target_end;
      if($exon->strand == 1) {
         $target_start = $exon->start + $k;
         $target_end = $target_start + 2;
      } else {
        $target_start = $exon->end - $k - 2;
        $target_end = $target_start + 2;
      }

      # Now need to look at the alignment index for this codon. If there is a gap in the query sequence at that index
      # then the codon should be skipped. If it isn't a gap at that codon position then
      my $codon_alignment_index = find_codon_alignment_index($codon_index,$target_seq);

      my $query_char = substr($query_seq,$codon_alignment_index,1);
      my $target_char = substr($target_seq,$codon_alignment_index,1);

      # check whether the selenocysteine exon attribute was added by the HiveCesar2 module
      if (exists($exon->{'selenocysteine'})) {

        if (substr($query_seq,$codon_alignment_index,3) eq 'NNN' and
            substr($target_seq,$codon_alignment_index,3) eq 'TGA') {

          my $seq_edit = Bio::EnsEMBL::SeqEdit->new(
                                       -CODE    => '_selenocysteine',
                                       -NAME    => 'Selenocysteine',
                                       -DESC    => 'Selenocysteine',
                                       -START   => $k/3,
                                       -END     => $k/3,
                                       -ALT_SEQ => 'U' # 1-letter code for selenocysteine amino acid Sec
                                       );
          my $translation = $transcript->translation();
          $translation->add_Attributes($seq_edit->get_Attribute());
        }

      }

      if($query_char eq '-') {
        say "FM2 TSTART: ".$target_start;
        say "FM2 TEND: ".$target_end;
        say "FM2 CODON INDEX: ".$codon_index;
        say "FM2 CODON ALIGNMENT INDEX: ".$codon_alignment_index;
        say "FM2 CODON CHARS: '".$nucleotide_array[$k].$nucleotide_array[$k+1].$nucleotide_array[$k+2]."'";
        say "FM2 ALIGNMENT CHARS: '".$query_char."'='".$target_char."'";
        say "FM2 SKIPPING CODON BECAUSE OF ALIGMENT GAP";
        $codon_index++;
        next;
      } elsif(($codon_alignment_index >= length($query_seq)-1) && $query_char eq '*' && $target_char eq '*') {
        say "FM2 LAST CODON IS STOP SO SKIPPING";
        last;
      }

      my $hit_start = find_hit_start($codon_alignment_index,$query_seq);
      my $hit_end = $hit_start;
      say "FM2 TSTART: ".$target_start;
      say "FM2 TEND: ".$target_end;
      say "FM2 HSTART: ".$hit_start;
      say "FM2 HEND: ".$hit_end;
      say "FM2 CODON INDEX: ".$codon_index;
      say "FM2 CODON ALIGNMENT INDEX: ".$codon_alignment_index;
      say "FM2 CODON CHARS: '".$nucleotide_array[$k].$nucleotide_array[$k+1].$nucleotide_array[$k+2]."'";
      say "FM2 ALIGNMENT CHARS: '".$query_char."'---'".$target_char."'";
      $codon_index++;

      push(@{$proto_supporting_features},{'tstart' => $target_start,
                                          'tend'   => $target_end,
                                          'hstart' => $hit_start,
                                          'hend'   =>$hit_end});
    }

    my $joined_supporting_features = join_supporting_features($proto_supporting_features,$exon->strand);
    my $exon_feature_pairs = [];
    foreach my $joined_feature (@{$joined_supporting_features}) {
      my $feature_pair = Bio::EnsEMBL::FeaturePair->new(
                                                        -start      => $joined_feature->{'tstart'},
                                                        -end        => $joined_feature->{'tend'},
                                                        -strand     => $exon->strand,
                                                        -hseqname   => $transcript->stable_id(),#$transcript->{'accession'},
                                                        -hstart     => $joined_feature->{'hstart'},
                                                        -hend       => $joined_feature->{'hend'},
                                                        -hcoverage  => $coverage,
                                                        -percent_id => $percent_id,
                                                        -slice      => $exon->slice,
                                                        -analysis   => $transcript->analysis);
      say "FM2 ADD SUPPORTING EVIDENCE START: ".$feature_pair->start;
      say "FM2 ADD SUPPORTING EVIDENCE END: ".$feature_pair->end;
      say "ADD SUPPORTING EVIDENCE hSTART: ".$feature_pair->hstart;
      say "ADD SUPPORTING EVIDENCE hEND: ".$feature_pair->hend;

      push(@{$exon_feature_pairs},$feature_pair);
      if (scalar(@{$exon_feature_pairs})) {
        my $final_exon_supporting_features = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => $exon_feature_pairs);
        $exon->add_supporting_features($final_exon_supporting_features);
      } else {
        warning("No supporting features added for exon.\nExon start: ".$exon->start."\nExon end: ".$exon->end);
      }
      
      # exon feature pairs which are not in order compared to the previous feature in the previous exon
      # will not be added to all exon supporting features to make the transcript supporting evidence
      # because the Core API assumes that both all seq region coordinates and hit sequence coordinates will
      # be in order but this would not be the case here
      if ($prev_feature_pair) {
        if (($feature_pair->start() < $prev_feature_pair->end() and $feature_pair->strand() == 1) or
            ($feature_pair->end() > $prev_feature_pair->start() and $feature_pair->strand() == -1)) {
          # if features are not sorted, do not add them
          warning("Feature pair not in order. Added to exon supporting features but not added to transcript supporting features.");
          next;
        } else {
          push(@{$all_exon_supporting_features},$feature_pair);
        }
      } else {
        push(@{$all_exon_supporting_features},$feature_pair);
      }
      $prev_feature_pair = $feature_pair;
    }
  }

  if (scalar(@$all_exon_supporting_features) > 0) {
    my $transcript_supporting_features = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => $all_exon_supporting_features);

    if ($transcript_supporting_features) {
      $transcript->add_supporting_features($transcript_supporting_features);
    } else {
      warning("There are some all_exon_supporting_features but no $transcript_supporting_features for transcript ".$transcript->dbID()." ".$transcript->stable_id());
    }
  } else {
    warning("There are no all_exon_supporting_features and, therefore, no transcript_supporting_features for transcript ".$transcript->dbID()." ".$transcript->stable_id());
  }
}

sub find_codon_alignment_index {
  my ($codon_index,$align_seq) = @_;

  my $align_index = -1;
  my $char_count = 0;
  for(my $j=0; $j<length($align_seq); $j++) {

    my $char = substr($align_seq,$j,1);
    if($char eq '-') {
      next;
    }

    if($char_count == $codon_index) {
      $align_index = $j;
      last;
    }

    $char_count++;

  }

  # For the last codon
  if($align_index == -1 && $char_count == $codon_index) {
    $align_index = length($align_seq) - 1;
  }
  unless($align_index >= 0) {
    throw("Did not find the alignment index for the codon");
  }

  return($align_index);
}

sub find_hit_start {
  my ($alignment_index,$align_seq) = @_;

  my $sub_seq = substr($align_seq,0,$alignment_index + 1);
  $sub_seq =~ s/\-//g;

  my $hit_start = length($sub_seq);

  return($hit_start);
}

sub join_supporting_features {
  my ($proto_supporting_features,$strand) = @_;

  # At this point all supporting features for the exon have been built on a per codon basis. Before making
  # a feature pair we want to combine adjacent supporting features. If both the genomic end of the left
  # feature matches the genomic start - 1 of the next then they can be joined if the hit end of the first
  # is hit start - 1 of the next

  my $joined_supporting_features = [];
  if($strand == 1) {
    for(my $i=0; $i<scalar(@{$proto_supporting_features})-1; $i++) {
      my $left_proto = $$proto_supporting_features[$i];
      my $right_proto = $$proto_supporting_features[$i+1];

      if($left_proto->{'tend'} == ($right_proto->{'tstart'} - 1)) {
        if($left_proto->{'hend'} == ($right_proto->{'hstart'} - 1)) {
          # If this is the case then the codons and hits are contiguous and so they can be joined
          $right_proto->{'tstart'} = $left_proto->{'tstart'};
          $right_proto->{'hstart'} = $left_proto->{'hstart'};
          $$proto_supporting_features[$i] = 0;
          $$proto_supporting_features[$i+1] = $right_proto;
        }
      }
    }
  } else {
    for(my $i=scalar(@{$proto_supporting_features})-1; $i>0; $i--) {
      my $left_proto = $$proto_supporting_features[$i];
      my $right_proto = $$proto_supporting_features[$i-1];

      say "FM2 PROTO LTS: ".$left_proto->{'tstart'};
      say "FM2 PROTO LTE: ".$left_proto->{'tend'}; 
      say "FM2 PROTO LHS: ".$left_proto->{'hstart'};
      say "FM2 PROTO LHE: ".$left_proto->{'hend'};
      say "FM2 PROTO RTS: ".$right_proto->{'tstart'};
      say "FM2 PROTO RTE: ".$right_proto->{'tend'};
      say "FM2 PROTO RHS: ".$right_proto->{'hstart'};
      say "FM2 PROTO RHE: ".$right_proto->{'hend'};

      if($left_proto->{'tend'} == ($right_proto->{'tstart'} - 1)) {
        if($left_proto->{'hend'} == ($right_proto->{'hstart'} + 1)) {
          # If this is the case then the codons and hits are contiguous and so they can be joined
          $right_proto->{'tstart'} = $left_proto->{'tstart'};
          $right_proto->{'hstart'} = $left_proto->{'hstart'};
          $$proto_supporting_features[$i] = 0;
          $$proto_supporting_features[$i-1] = $right_proto;
        }
      }
    }
  }
  foreach my $proto_sf (@{$proto_supporting_features}) {
    if($proto_sf) {
      say "FM2 JOINED TSTART: ".$proto_sf->{'tstart'};
      say "FM2 JOINED TEND: ".$proto_sf->{'tend'};
      say "FM2 JOINED HSTART: ".$proto_sf->{'hstart'};
      say "FM2 JOINED HEND: ".$proto_sf->{'hend'};

      # If it's the negative strand then swap the start and end of the hit
      if($strand == -1) {
        my $temp = $proto_sf->{'hstart'};
        $proto_sf->{'hstart'} = $proto_sf->{'hend'};
        $proto_sf->{'hend'} = $temp;
      }
      push(@{$joined_supporting_features},$proto_sf);
    }
  }

  return($joined_supporting_features);

}

1;
