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
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(print_Translation clone_Translation print_peptide);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils qw(print_Exon clone_Exon Exon_info exon_length_less_than_maximum Exon_info get_upstream_Intron get_downstream_Intron get_upstream_splice_sites get_downstream_splice_sites validate_Exon_coords);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::IntronUtils qw(intron_length_less_than_maximum get_splice_sites);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(coord_string id empty_Object);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::EvidenceUtils qw (print_Evidence clone_Evidence);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(write_seqfile);
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);

use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Seg;
use Bio::SeqIO;

use vars qw (@ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
             are_phases_consistent
             are_splice_sites_canonical
             are_strands_consistent
             attach_Slice_to_Transcript
             calculate_exon_phases
             clone_Transcript
             coding_coverage
             convert_to_genes
             count_non_canonical_splice_sites
             dump_cDNA_file
             empty_Transcript
             evidence_coverage
             exon_lengths_all_less_than_maximum
             exonic_proportion
             get_downstream_Intron_from_Exon
             get_downstream_splice_sites
             get_evidence_ids
             get_upstream_splice_sites
             get_upstream_Intron_from_Exon
             has_no_unwanted_evidence
             intron_lengths_all_less_than_maximum
             is_not_folded
             is_spliced 
             list_evidence
             low_complexity_less_than_maximum
             replace_stops_with_introns
             remove_initial_or_terminal_short_exons
             get_evidence_ids
             dump_cDNA_file
             convert_to_genes
             evidence_coverage
             get_upstream_Intron_from_Exon
             get_downstream_Intron_from_Exon
             get_upstream_splice_sites
             get_downstream_splice_sites
             attach_Slice_to_Transcript
             empty_Transcript
             fully_load_Transcript
             all_exons_are_valid
             print_Transcript
             print_Transcript_and_Exons
             print_Transcript_evidence
             split_Transcript
             tidy_split_transcripts
             trim_cds_to_whole_codons
             Transcript_info
             evidence_coverage_greater_than_minimum
             count_real_introns
             identical_Transcripts
            );


=head2 print_Transcript

  Arg [1]   : Bio::EnsEMBL::Transcript or Aref of Bio::EnsEMBL::Transcript-objects
  Arg [2]   : string, this should be a string or spaces or tabs to indent the 
              printed string

  Function  : print information about the transcript and its
  children objects, using indent to make the format readable
  Returntype: n/a
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
  Returntype: n/a
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
    map { print Exon_info($_, $indent."\t")."\n"  } @{$transcript->get_all_Exons} ; 
  }
}


=head2 Transcript_info

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : string, indent
  Function  : return string of info about the transcript
  Returntype: 
  Exceptions: none
  Example   : print_just_Transcript($transcript);

=cut



sub Transcript_info{
  my ($transcript, $indent) = @_;
  $indent = '' if(!$indent);
  my $coord_string = coord_string($transcript); 
  my $id = id($transcript);
  return $indent."TRANSCRIPT: ".$id." ".$coord_string;
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
  my $newtranslation;
  if($transcript->translation){
    $newtranslation = clone_Translation($transcript, 
                                        $newtranscript);
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
    warning("Strands are inconsistent between the ".
            "first exon and the transcript for ".
            id($transcript));
    return 0;
  }
  if(@$exons >= 2){
    for(my $i = 1;$i < @$exons;$i++){
      if($exons->[$i]->strand != $exons->[$i-1]->strand){
        warning("Strands are inconsistent between ".
                "exon $i exon and exon ".($i-1)." for ".
                id($transcript));
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
      logger_info("Transcript ".id($transcript)." has ".
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
      warning("Transcript ".id($transcript)." has ".
              "intron longer than ".$max_length);
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
        logger_info("Non-coding transcript does not have -1 phases");
        return 0;
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
          logger_info("Coding transcript has inconsistent phases");
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
  Returntype:boolean, 1 for pass, 0 for fail 
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
        warning(id($transcript)." is folded");
        logger_info($i." ".id($exons->[$i])." has a start which ".
                    "is less than ".($i-1)." ".id($exons->[$i-1]).
                    " end");
        return 0;
      }
    }else{
      if($exons->[$i]->end > $exons->[$i-1]->start){
        warning(id($transcript)." is folded");
        logger_info($i." ".id($exons->[$i])." has a end which ".
                    "is greater than ".($i-1)." ".id($exons->[$i-1]).
                    " start");
        return 0;
      }
    }
  }
  return 1;
}



=head2 low_complexity_less_than_maximum

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : int, maximum low complexity
  Function  : calculate how much low complexity a
  transcripts translation has and check if it is above
  the specificed threshold
  Returntype: boolean, 1 for pass 0 for fail
  Exceptions: none
  Example   : 

=cut



sub low_complexity_less_than_maximum{
  my ($transcript, $complexity_threshold) = @_;
  my $peptide = $transcript->translate;
  my $seg = Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Seg->new
    (
     -query => $peptide,
     -analysis => Bio::EnsEMBL::Analysis->new
     (
      -logic_name => 'seg',
      -program_file => 'seg',
     )
    );
  $seg->run;
  my $low_complexity = $seg->get_low_complexity_length;
  logger_info(id($transcript)." has ".$low_complexity.
              " low complexity sequence");
  if($low_complexity > $complexity_threshold){
    warning(id($transcript)."'s low ".
            "complexity is above ".
            "the threshold ".$complexity_threshold.
            "\n");
    return 0;
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
        warning(id($transcript)." has ".$evi.
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
  warning(id($transcript)." has no introns ".
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
  my $warn = id($transcript)." has no introns. count_real".
    "_introns makes no sense in those terms";
  logger_info($warn) if(@introns == 0);
  foreach my $intron(@introns){
    $real_count++ if($intron->length > $intron_size);
  }
  logger_info(id($transcript)." has ".$real_count." introns ".
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
    logger_info(id($transcript)." contains ".
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
  Example   : 

=cut

#This method is designed to split transcripts on
#long introns and hopefully maintain clean and sensible
#transcripts

sub split_Transcript{
  my ($transcript, $max_intron_length, $intron_numbers) = @_;
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
      $stran = trim_cds_to_whole_codons($stran);
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
      my $warn = id($stran)." only has one exon\n".
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
  Exceptions: 
  Example   : 
  Notes     : Needs extension to deal with transcript_supporting_feats
=cut

sub replace_stops_with_introns{
  my ($transcript) = @_;

  my $newtranscript = clone_Transcript($transcript);
  my @exons = @{$newtranscript->get_all_Exons};
  my $pep = $newtranscript->translate->seq;

  while($pep =~ /\*/g) {
    my $position = pos($pep);

    my @coords = $newtranscript->pep2genomic($position, $position);

    foreach my $stop (@coords) {
      # locate the exon that this stop lies in
      my @new_exons;
      foreach my $exon (@exons) {
        if ($stop->start >= $exon->start and $stop->end <= $exon->end) {
          # this stop lies completely within an exon. We split the exon
          # into two
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
          # need to split the supporting features between the two
          my @sfs = @{$exon->get_all_supporting_features};
          my (@ug_left, @ug_right);
          foreach my $f (@sfs) {
            foreach my $ug ($f->ungapped_features) {
              if ($ug->start >= $exon_left->start and 
                  $ug->end <= $exon_left->end) {
                #completely within the left-side of the split
                push @ug_left, $ug;
              } elsif ($ug->start >= $exon_right->start and 
                       $ug->end <= $exon_right->end) {
                #completely within the right-side of the split
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
          
          if (@ug_left) {
            my $f = Bio::EnsEMBL::DnaPepAlignFeature->
                new(-features => \@ug_left);
            $exon_left->add_supporting_features($f);
          }
          if (@ug_right) {
            my $f = Bio::EnsEMBL::DnaPepAlignFeature->
                new(-features => \@ug_right);
            $exon_right->add_supporting_features($f);
          }
          
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
  $min_length = 3 unless($min_length);
  throw("TranscriptUtils::remove_initial_or_terminal_short_exons will not work ".
        " if ".id($transcript)." has no translation") if(!$transcript->translation);
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
      }else{
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
              id($transcript)." into ".$filename." format ".
              $format);
  my $seq = $transcript->seq;
  $seq->display_id(id($transcript));
  $filename = write_seqfile($seq, $filename, $format);
  return $filename;
}

=head2 convert_to_genes( $tref, $analysis )

Arg[0]     : Arrayref of Bio::EnsEMBL::Transcript objects
Arg[1]     : Bio::EnsEMBL::Analysis object (opt)
Name       : convert_to_genes
Function   : converts all transcripts to Genes (biotype of transcript becomes biotype of gene)
Returntype : Arrayref of Bio::EnsEMBL::Gene objects

=cut

sub convert_to_genes {
  my ($tref, $analysis, $biotype ) = @_;
  my @genes ;

  my @tr = (ref($tref)=~m/ARRAY/) ? @$tref : ($tref) ; 

  for my $t (@tr) {
    $analysis = $t->analysis unless $analysis ;
    $t->biotype($biotype);
    my $g = Bio::EnsEMBL::Gene->new() ;
    $g->biotype( $t->biotype ) ;
    $g->add_Transcript($t);
    $g->analysis($analysis) ;
    push @genes, $g ;
  }
  for my $g(@genes) {
    throw("there are no transcripts for this gene\n") if scalar(@{$g->get_all_Transcripts}) == 0 ;
    for my $tr ( @{$g->get_all_Transcripts} ) {
      throw("there are no exons  for this transcript \n") if scalar(@{$tr->get_all_Exons}) == 0 ;
    }
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

  throw("Can't work out evidence coverage for ".id($transcript).
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


sub fully_load_Transcript{
  my ($transcript, $keep_xrefs) = @_;

  $keep_xrefs = 1 if(!defined($keep_xrefs));

  if ($transcript->translation) {
    $transcript->translate;
    $transcript->translation->get_all_Attributes;
    $transcript->translation->get_all_DBEntries if($keep_xrefs);
    $transcript->translation->get_all_ProteinFeatures;
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

  return $transcript;
}

sub empty_Transcript{
  my ($transcript, $remove_stable_id, $remove_xrefs) = @_;
  fully_load_Transcript($transcript);
  foreach my $sf(@{$transcript->get_all_supporting_features}){
    empty_Object($sf);
  }
  empty_Object($transcript->translation, $remove_stable_id) 
    if($transcript->translation);;
 EXON:foreach my $e(@{$transcript->get_all_Exons}){
  SF:foreach my $sf(@{$e->get_all_supporting_features}){
      empty_Object($sf, $remove_stable_id);
    }
    empty_Object($e, $remove_stable_id);
  }
  $transcript->display_xref(undef) if($remove_xrefs);
  empty_Object($transcript, $remove_stable_id);
  return $transcript;
}

sub all_exons_are_valid{
  my ($transcript, $max_length) = @_;

  foreach my $exon(@{$transcript->get_all_Exons}){
    throw(Transcript_info($transcript)." seems to contain an undefined exon") 
      if(!$exon); 
    if(!exon_length_less_than_maximum($exon, $max_length)){
      logger_info("Transcript ".id($transcript)." has ".
                  "exon longer than ".$max_length);
      return 0;
    }
    if(!validate_Exon_coords($exon)){
      logger_info("Transcript ".id($transcript)." has ".
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
           Transcript, Transcript
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
1;
