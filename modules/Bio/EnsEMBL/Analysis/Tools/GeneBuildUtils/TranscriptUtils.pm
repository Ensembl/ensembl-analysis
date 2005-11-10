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
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(print_Translation clone_Translation);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils qw(print_Exon clone_Exon Exon_info 
                                                                exon_length_less_than_maximum);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::IntronUtils qw(intron_length_less_than_maximum get_splice_sites);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(coord_string id);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::EvidenceUtils qw (print_Evidence clone_Evidence);
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_verbosity
                                             logger_info
                                             logger_warning);
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Seg;
use vars qw (@ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
             print_Transcript
             clone_Transcript
             print_Transcript_evidence
             Transcript_info
             are_strands_consistent
             lies_inside_of_slice
             exon_lengths_all_less_than_maximum
             intron_lengths_all_less_than_maximum
             are_phases_consistent
             is_not_folded
             low_complexity_less_than_maximum
             has_no_unwanted_evidence
             is_spliced 
             are_splice_sites_canonical
             count_non_canonical_splice_sites
            );




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
  print Transcript_info($transcript, $indent)."\n";
  my $translation_indent = $indent."\t";
  print_Translation($transcript, $translation_indent);
  foreach my $exon(@{$transcript->get_all_Exons}){
    my $exon_indent = $translation_indent."\t";
    print_Exon($exon, $exon_indent);
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
  return $indent."TRANSCRIPT: ".$id." ".$coord_string."\n";
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
  my $newtranslation = clone_Translation($transcript, 
                                         $newtranscript);
  $newtranscript->translation($newtranslation);
  my $attribs = $transcript->get_all_Attributes();
  $newtranscript->add_Attributes(@$attribs);
  $newtranscript->slice($transcript->slice);
  $newtranscript->dbID($transcript->dbID);
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
    logger_warning("Strands are inconsistent between the ".
                   "first exon and the transcript for ".
                   id($transcript));
    return 0;
  }
  if(@$exons >= 2){
    for(my $i = 1;$i < @$exons;$i++){
      if($exons->[$i]->strand != $exons->[$i-1]->strand){
        logger_warning("Strands are inconsistent between ".
                       "exon $i exon and exon ".($i-1)." for ".
                       id($transcript));
        return 0;
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
  Returntype: boolean, 1 for pass 0 for fail ie(lies outside
                                                    of slice or
                                                    across lower 
                                                    boundary)
  Exceptions: none
  Example   : 

=cut


sub lies_inside_of_slice{
  my ($transcript, $slice) = @_;
  if($transcript->start > $slice->length || 
     $transcript->end < 1){
    logger_warning(id($transcript)." lies off edge if slice ".
                   $slice->name);
    return 0;
  }
  if($transcript->start < 1 && $transcript->end > 1){
    logger_warning(id($transcript)." lies over lower boundary".
                   " of slice ".$slice->name);
    return 0;
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
      logger_info("Transcript ".id($transcript)." has ".
                  "intron longer than ".$max_length);
      return 0;
    }
  }
  if(@{$transcript->get_all_Introns} == 0){
    my $warn = "intron_lengths_less_than_maximum is an ".
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
  my $end_phase = $exons->[0]->end_phase;
  if(@$exons == 1){
    my $warn = "are_phases_consistent ".
      "is an inappropriate test for a single ".
        "exon gene";
    logger_info($warn);
    return 1;
  }
  for(my $i=1;$i < @$exons;$i++){
    if($exons->[$i]->phase != $end_phase){
      my $warn = (id($transcript)." has inconsistent ".
        "phases between exon ".id($exons->[$i])."-".
          $i." and ".id($exons->[$i-1])."-".($i-1));
      logger_warning($warn)
        unless($end_phase == -1 &&
               $exons->[-1]->phase != -1 || 
               $end_phase != -1 &&
               $exons->[-1]->phase == -1);;
      return 0 unless($end_phase == -1 &&
                          $exons->[-1]->phase != -1 || 
                          $end_phase != -1 &&
                          $exons->[-1]->phase == -1);
    }
    $end_phase = $exons->[$i]->end_phase;
  }
  return 1;
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
        logger_warning(id($transcript)." is folded");
        logger_info($i." ".id($exons->[$i])." has a start which ".
                    "is less than ".($i-1)." ".id($exons->[$i-1]).
                    " end");
        return 0;
      }
    }else{
      if($exons->[$i]->end > $exons->[$i-1]->start){
        logger_warning(id($transcript)." is folded");
        logger_info($i." ".id($exons->[$i])." has a end which ".
                    "is greater than ".($i-1)." ".id($exons->[$i-1]).
                    " start");
        return 0;
      }
    }
  }
  return 1;
}



=head2 low_complexity_below_maximum

  Arg [1]   : Bio::EnsEMBL::Transcript
  Arg [2]   : int, maximum low complexity
  Function  : calculate how much low complexity a
  transcripts translation has and check if it is above
  the specificed threshold
  Returntype: boolean, 1 for pass 0 for fail
  Exceptions: none
  Example   : 

=cut



sub low_complexity_below_maximum{
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
    logger_warning(id($transcript)."'s low ".
                   "complexity is above ".
                   "the threshold ".$complexity_threshold.
                   "\n");
    return 0;
  }
  return 1;
}


sub has_no_unwanted_evidence{
  my ($transcript, $ids) = @_;
  $ids = {} if(!$ids);
  $ids->{'NG_'} = 1;
  my $return = 1;
  my $evidence = _get_evidence_ids($transcript);
  foreach my $evidence(keys(%$evidence)){
    foreach my $unwanted(keys(%$ids)){
      if($evidence =~ /$unwanted/){
        logger_warning(id($transcript)." has ".$evidence.
                       " unwanted evidence");
        return 0;
      }
    }
  }
  return 1;
}



=head2 _get_evidence_ids

  Arg [1]   : Bio::EnsEMBL::Transcript
  Function  : get a unique hash of evidence ids from the given
  transcript. Note this is curretnly an internal method and is
  not exported from the module
  Returntype: hashref of unique ids 
  Exceptions: none
  Example   : 

=cut


sub _get_evidence_ids{
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



sub is_spliced{
  my ($transcript, $intron_size) = @_;
  my $count = count_real_introns($transcript, $intron_size);
  logger_warning(id($transcript)." has no introns ".
                 "longer than $intron_size bps") if(!$count);
  return 0 if(!$count);
  return 1;
}

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


sub are_splice_sites_canonical{
  my ($transcript) = @_;
  my $introns = $transcript->get_all_Introns;
  my $non_canonical_count = 
    count_non_canonical_splice_sites($transcript);
  if($non_canonical_count){
    logger_warning(id($transcript)." contains ".
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



##METHODS NEEDED

#checks



#perhaps also what percentage is useful?
#orf coverage, how does the spliced length compare to the genomic extent?

#utilitys
#split trancripts, split transcripts on long introns
#replace stops with introns, frameshift around stop codons!
#list_evidence, a list of ids that support the gene

1;
