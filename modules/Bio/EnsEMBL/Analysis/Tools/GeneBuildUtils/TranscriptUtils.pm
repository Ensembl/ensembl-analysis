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

please send any questions to dev@ensembl.org

=head1 METHODS

the rest of the documention details the exported static
class methods

=cut

# $Source: /tmp/ENSCOPY-ENSEMBL-ANALYSIS/modules/Bio/EnsEMBL/Analysis/Tools/GeneBuildUtils/TranscriptUtils.pm,v $
# $Revision: 1.85 $
package Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils;

use strict;
use warnings;
use Exporter;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning stack_trace_dump);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(print_Translation clone_Translation print_peptide);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils qw(print_Exon clone_Exon Exon_info exon_length_less_than_maximum Exon_info get_upstream_Intron get_downstream_Intron get_upstream_splice_sites get_downstream_splice_sites validate_Exon_coords);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::IntronUtils qw(intron_length_less_than_maximum get_splice_sites);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(seq_region_coord_string id empty_Object);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::EvidenceUtils qw (print_Evidence clone_Evidence);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(write_seqfile);
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);

use Bio::EnsEMBL::Analysis;
use Bio::SeqIO;

use vars qw (@ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
             all_exons_are_valid
             are_phases_consistent
             are_splice_sites_canonical
             are_strands_consistent
             attach_Slice_to_Transcript
             attach_Analysis_to_Transcript
             attach_Analysis_to_Transcript_no_support
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
             exonic_proportion
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
             print_Transcript
             print_Transcript_and_Exons
             print_Transcript_evidence
             remove_initial_or_terminal_short_exons
             replace_stops_with_introns
             set_start_codon
             set_stop_codon
             split_Transcript
             tidy_split_transcripts
             trim_cds_to_whole_codons
             Transcript_info
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
  return $indent."TRANSCRIPT: ".$id." ".$coord_string;
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

  $transcript->sort;
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

  #$transcript->sort;
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

  #$transcript->sort; Commented out as it has been deprecated. Transcript are now sorted by default.
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
      $cds_len += $start_phase;
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

  #$transcript->sort;
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
      if ($stran->end >= $cds_start and
          $stran->start <= $cds_end) {
        # at least part of this transcript is coding
        my $tr = Bio::EnsEMBL::Translation->new;
        $stran->translation($tr);

        my @exons = @{$stran->get_all_Exons};
        foreach my $e (@exons) {
          if ($cds_start >= $e->start and $cds_start < $e->end) {
            # start of translation is in this exon
            if ($stran->strand > 0) {
              $tr->start_Exon($e);
              $tr->start( $cds_start - $e->start + 1);
            } else {
              $tr->end_Exon($e);
              $tr->end( $e->end - $cds_start + 1);
            }
          }
          if ($cds_end >= $e->start and $cds_end <= $e->end) {
            if ($stran->strand > 0) {
              $tr->end_Exon($e);
              $tr->end( $cds_end - $e->start + 1);
            } else {
              $tr->start_Exon($e);
              $tr->start( $e->end - $cds_end + 1);
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
  my ($transcript) = @_;

  $transcript->sort;
  my $newtranscript = clone_Transcript($transcript);

  my @exons = @{$newtranscript->get_all_Exons};
  my $pep = $newtranscript->translate->seq;

  # gaps adjacent to internal stop codons - skip
  return 0 if ($pep =~ /X\*/ || $pep =~ /\*X/);

  while($pep =~ /\*/g) {
    my $position = pos($pep);

    my @coords = $newtranscript->pep2genomic($position, $position);

    foreach my $stop (@coords) {
      # locate the exon that this stop lies in
      my @new_exons;
      foreach my $exon (@exons) {
        
        if (($stop->start > ($exon->start + 2)) && ($stop->end < $exon->end - 2)) {
          # this stop lies _completely_ within an exon and not on it's
          # boundary. We therefore can split the exon into two 
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

          # if either of the exons are exactly 3bp long, do not
          # span the gap with an intron, instead prune the tiny exon and intron

          # check we don't have a triplet
          if (($exon_left->length >= 3) && ($exon_right->length >= 3)) {
            # both exons are >3bp in length
            # need to split the supporting features between the two
            my @sfs = @{$exon->get_all_supporting_features}; 
            my (@ug_left, @ug_right);

            foreach my $f (@sfs) { 
              foreach my $ug ($f->ungapped_features) {  
                $ug->analysis($newtranscript->analysis);
                my $orignial_analysis = $ug->analysis;
                if (($ug->start == $exon_left->start && $ug->end == $exon_left->end) ||
                  ($ug->start == $exon_right->start && $ug->end == $exon_right->end)) {
                  print STDERR "There's one base in it - cannot split due to out of phase error\n";
                  last;
                } elsif (($ug->start + 2 > $exon_left->start) || ($ug->end < $exon_left->end -2)) {
                 print STDERR "There's two bases in it - cannot split due to out of phase error\n";
                 last;
                } elsif ($ug->start >= $exon_left->start && 
                    $ug->end <= $exon_left->end) {
                  # completely within the left-side of the split
                  push @ug_left, $ug;

                } elsif ($ug->start >= $exon_right->start && 
                         $ug->end <= $exon_right->end) {
                 
                  # completely within the right-side of the split
                  push @ug_right, $ug;
                } else {
                  # this ug must span the split
                  my $fp_left = Bio::EnsEMBL::FeaturePair->new();
                  if ($ug->slice) {
                    $fp_left->slice($ug->slice);
                  }
                 $fp_left->seqname   ($ug->seqname);
                 $fp_left->strand    ($ug->strand);
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
                    $fp_left->hend($fp_left->hstart +
                                   ($fp_left->length / 3) - 
                                   1);
                  
                    $fp_right->hend ($ug->hend);
                    $fp_right->hstart($ug->hend - 
                                      ($fp_right->length / 3) + 
                                      1);
                  } else {
                    # we are on the other strand
                    $fp_left->hend ($ug->hend);
                    $fp_left->hstart($ug->hend - 
                                     ($fp_left->length / 3) + 
                                     1);
                 
                    $fp_right->hstart($ug->hstart);
                    $fp_right->hend($fp_right->hstart +
                                   ($fp_right->length / 3) - 
                                   1);
                 }
                
                 if ($fp_left->end >= $fp_left->start) { 
                    push @ug_left, $fp_left;
                 }
                 if ($fp_right->end >= $fp_right->start) {
                    push @ug_right, $fp_right;
                 }
               }
             }
           } 

           sub add_dna_align_features_by_hitname_and_analysis {   
              my ( $ug_ref, $exon ) = @_ ;  
              my %group_features_by_hitname_and_analysis ; 
              for my $ug ( @$ug_ref ) {  
                 push @{$group_features_by_hitname_and_analysis{$ug->analysis->logic_name}{$ug->hseqname }} , $ug ; 
              }  
              for my $logic_name ( keys %group_features_by_hitname_and_analysis  ) {  
                 for my $hseqname  ( keys  %{$group_features_by_hitname_and_analysis{$logic_name}}) {  
                    my @features = @{$group_features_by_hitname_and_analysis{$logic_name}{$hseqname}}; 
                    my $f = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@features); 
                    $exon->add_supporting_features($f);  
                 } 
              }  
              return $exon ; 
            }   

           $exon_left = add_dna_align_features_by_hitname_and_analysis(\@ug_left,$exon_left) ; 
           $exon_right =add_dna_align_features_by_hitname_and_analysis(\@ug_right,$exon_right) ;  

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
          } elsif ($exon_left->length == 3) {
            # Need to prune a tiny left exon
            push @new_exons, $exon_right;

          } elsif ($exon_right->length == 3) {
            # Need to prune a tiny right exon
            push @new_exons, $exon_left;
          } else {
            verbose("new Exon must have length that is =3bp. There's nothing we can do to fix this stop codon.");
            push @new_exons, $exon;
          }
        } elsif ($stop->start == $exon->start -1) {
          # stop lies at the start of the exon
          $exon->start($exon->start + 3); 
          push @new_exons, $exon;
          die;
        } elsif ($stop->end == $exon->end ) {
          # stop lies at the end of the exon
          $exon->end($exon->end - 3); 
          push @new_exons, $exon;
        } else {
          # this exon is unaffected by this stop
          push @new_exons, $exon;
        }
      }
      
      @exons = @new_exons;
    }
  }
  
  $newtranscript->flush_Exons; 
  foreach my $exon (@exons) {
    $newtranscript->add_Exon($exon);
  } 

  my $translation = Bio::EnsEMBL::Translation->new();
  $translation->start_Exon($exons[0]);
  $translation->end_Exon($exons[-1]);
  $translation->start(1);
  $translation->end($exons[-1]->end - $exons[-1]->start + 1);
  $newtranscript->translation($translation); 

  my $old_translation = $transcript->translation  ;  

  foreach my $DBEntry (@{$old_translation->get_all_DBEntries}){
     $translation->add_DBEntry($DBEntry);
  }
  return $newtranscript;
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

  $transcript->sort;
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
  $transcript->slice($slice);
  foreach my $sf(@{$transcript->get_all_supporting_features}){
    $sf->slice($slice);
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
  $transcript->analysis($analysis);
  foreach my $sf(@{$transcript->get_all_supporting_features}){
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

  $transcript->sort;
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
  warning ("There has been no evidence passed in for ".Transcript_info($transcript).
           "Assumung coverage is fine") if(!$evidence);
  return 1 if(!$evidence);
  my $coverage = evidence_coverage($transcript, $evidence);
  warning(Transcript_info($transcript)." based on ".$evidence->id." has too ".
          "low coverage ".$coverage." of the evidence") 
    unless($coverage > $min_coverage);
  return 0 unless($coverage > $min_coverage);
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

  $transcript->sort();

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

1;
