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

Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils - utilities for transcript objects

=head1 SYNOPSIS

  use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils qw(Exon_info);

  or 

  use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils 

  to get all methods

=head1 DESCRIPTION

All methods in this class should take a Bio::EnsEMBL::Exon
object as their first argument.

The methods provided should carry out some standard 
functionality for said objects such as printing info or
transfering evidence

=head1 CONTACT

please send any questions to http://lists.ensembl.org/mailman/listinfo/dev

=head1 METHODS

the rest of the documention details the exported static
class methods

=cut


package Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils;

use strict;
use warnings;
use Exporter;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning
                                      stack_trace_dump);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(seq_region_coord_string id);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::EvidenceUtils qw(print_Evidence clone_Evidence);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::IntronUtils qw(get_splice_sites);
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);
use vars qw (@ISA @EXPORT);

@ISA = qw(Exporter);
@EXPORT = qw(
             print_Exon
             clone_Exon
             Exon_info
             exon_length_less_than_maximum
             transfer_supporting_evidence
             merge_exons
             get_upstream_Intron
             get_downstream_Intron
             get_upstream_splice_sites
             get_downstream_splice_sites
             create_Exon
             validate_Exon_coords
            );




=head2 print_Exon

  Arg [1]   : Bio::EnsEMBL::Exon
  Arg [2]   : string, indent (\t) to tag infront of printed

  Function  : print information about given exon and
  associated evidence
  Returntype: none
  Exceptions: throws if not passed a first argument
  Example   : 

=cut



sub print_Exon{
  my ($exon, $indent ) = @_;
  $indent = "" if(!$indent);
  throw("Must be passed an exon") if(!$exon);
  print Exon_info($exon, $indent)."\n";
  foreach my $evidence(@{$exon->get_all_supporting_features}){
    my $evidence_indent = $indent."\t";
    print_Evidence($evidence, $evidence_indent);
  }
}



=head2 Exon_info

  Arg [1]   : Bio::EnsEMBL::Exon
  Arg [2]   : indent, string
  Function  : returns a string with info about exon in it
  Returntype: string
  Exceptions: throws if not passed a first argument
  Example   : 

=cut



sub Exon_info{
  my ($exon, $indent) = @_;
  throw("Must be passed an exon") if(!$exon);
  $indent = '' if(!$indent);
  my $coord_string = seq_region_coord_string($exon);
  my $id = id($exon);
  return $indent."EXON: ".$id." ".$coord_string."  ".
    $exon->phase." ".$exon->end_phase." ".
      $exon->length;
}


=head2 clone_Exon

  Arg [1]   : Bio::EnsEMBL::Exon
  Function  : produces a whole new exon object identical
  to the object passed
  Returntype: Bio::EnsEMBL::Exon
  Exceptions: none
  Example   : 

=cut



sub clone_Exon{
  my ($exon) = @_;

  my @sfs;

  foreach my $sf(@{$exon->get_all_supporting_features}){
    my $newsf = clone_Evidence($sf);
    push(@sfs, $newsf);
  }
  
    
  my $newexon = create_Exon($exon->start, $exon->end, $exon->phase, $exon->end_phase,
                            $exon->strand, $exon->analysis, \@sfs, $exon->dbID,
                            $exon->slice, $exon->stable_id, $exon->version);
  return $newexon;
}


sub create_Exon{
  my ($start, $end, $phase, $end_phase, $strand, $analysis, $sfs, $dbID, $slice, $stable_id, $version) = @_;
  
  my $newexon = Bio::EnsEMBL::Exon->new();
  $newexon->start      ($start);
  $newexon->end        ($end);
  $newexon->phase      ($phase);
  $newexon->end_phase  ($end_phase);
  $newexon->strand     ($strand);
  $newexon->analysis   ($analysis);
  $newexon->dbID       ($dbID);
  $newexon->slice      ($slice);
  $newexon->stable_id  ($stable_id);
  $newexon->version    ($version);
  foreach my $sf(@$sfs){
    $newexon->add_supporting_features($sf);
  }
  return $newexon;
}


=head2 exon_length_less_than_maximum

  Arg [1]   : Bio::EnsEMBL::Exon
  Arg [2]   : int, max length
  Function  : checks if exon is longer than specified length
  Returntype: boolean, 0 for longer, 1 for shorter
  Exceptions: none
  Example   : 

=cut



sub exon_length_less_than_maximum{
  my ($exon, $max_length) = @_;
  if($exon->length >= $max_length){
    warning(id($exon)." is longer than max length ".
            $max_length);
    return 0;
  }
  return 1;
}



=head2 transfer_supporting_evidence

  Arg [1]   : Bio::EnsEMBL::Exon
  Arg [2]   : Bio::EnsEMBL::Exon
  Function  : transfer evidence from the first exon to the 
  second exon
  Returntype: Bio::EnsEMBL::Exon
  Exceptions: 
  Example   : 

=cut

#Note this method modifies the target exon which is_
#passed in

sub transfer_supporting_evidence{
  my ($source_exon, $target_exon) = @_;
  my %target_evidence;
  my %source_evidence;
  foreach my $sf
    (@{$target_exon->get_all_supporting_features}){
      my $unique_id = $sf->start."-".$sf->end."-".
        $sf->strand."-".$sf->hseqname."-".$sf->cigar_string;
      $target_evidence{$unique_id} = $sf;
    }
  foreach my $sf
    (@{$source_exon->get_all_supporting_features}){
      my $unique_id = $sf->start."-".$sf->end."-".
        $sf->strand."-".$sf->hseqname."-".$sf->cigar_string;
      if(!$target_evidence{$unique_id}){
        $source_evidence{$unique_id} = $sf;
      }
    }

  foreach my $sf(values(%source_evidence)){
  	if ($sf->overlaps($target_exon)) {
      logger_info("Adding ".$sf->hseqname." to the target exon");
      $target_exon->add_supporting_features($sf);
  	} else {
      logger_info("Not adding ".$sf->hseqname." to the target exon because there is no overlap");
  	}
  }
  return $target_exon;
}


sub validate_Exon_coords{
  my ($exon, $allow_negative_start) = @_;
  throw("Must pass ExonUtils::validate_Exon_coords an exon ") if(!$exon);
  if(!$allow_negative_start && ($exon->start < 0)){
    warning(id($exon)." ".$exon->start." is less than zero");
    return 0;
  }
  elsif($exon->start > $exon->end){
    warning(id($exon)." ".$exon->start." is greater than ".$exon->end);
    return 0;
  }
  return 1;
}

=head2 merge_exons

  Arg [1]   : Bio::EnsEMBL::Exon
  Arg [2]   : Bio::EnsEMBL::Exon
  Function  : merges the two exons by creating a new one
              setting the start to the first exon start and
              the end to the second exon end, and copying all
              supporting evidence to the new merged exon
  Returntype: Bio::EnsEMBL::Exon
  Exceptions: none
  Example   : 

=cut

sub merge_exons{
  my ($exon1,$exon2) = @_;

  my @sfs;

  foreach my $sf (@{$exon1->get_all_supporting_features()},@{$exon2->get_all_supporting_features()}){
    my $newsf = clone_Evidence($sf);
    push(@sfs,$newsf);
  }

  my $merged_start;
  my $merged_end;
  my $merged_phase;
  my $merged_end_phase;
  
  if ($exon1->strand() == 1) {
  	$merged_start = $exon1->start();
  	$merged_end = $exon2->end();
  	$merged_phase = $exon1->phase();
  	$merged_end_phase = $exon2->end_phase();
  } else {
  	$merged_start = $exon2->start();
    $merged_end = $exon1->end();
    $merged_phase = $exon1->phase();
    $merged_end_phase = $exon2->end_phase();
  }

  return create_Exon($merged_start,$merged_end,$merged_phase,$merged_end_phase,
                     $exon1->strand(),$exon1->analysis(),\@sfs,$exon1->dbID(),
                     $exon1->slice(),$exon1->stable_id(),$exon1->version());
}

=head2 Intron methods

  Arg [1]   : Bio::EnsEMBL::Exon
  Arg [2]   : Bio::EnsEMBL::Transcript
  Function  : these methods get either intron objects or splice site pairs
  on the basis of a give exon and transcript
  Returntype: either Bio::EnsEMBL::Intron or 2 strings. Note the 5" most exon will give
  no upstream intron/splice sites and the 3" most exon will give no downstream
  Exceptions: 
  Example   : 

=cut



sub get_upstream_Intron{
  my ($exon, $transcript) = @_;
  foreach my $intron(@{$transcript->get_Introns}){
    if($exon eq $intron->next_Exon){
      return $intron;
    }
  }
  return undef;
}

sub get_downstream_Intron{
  my ($exon, $transcript) = @_;
  foreach my $intron(@{$transcript->get_Introns}){
    if($exon eq $intron->prev_Exon){
      return $intron;
    }
  }
  return undef;
}

sub get_upstream_splice_sites{
  my ($exon, $transcript) = @_;
  my $upstream_intron = get_upstream_Intron($exon, $transcript);
  if($upstream_intron){
    return get_splice_sites($upstream_intron);
  }else{
    return undef;
  }
}


sub get_downstream_splice_sites{
  my ($exon, $transcript) = @_;
  my $downstream_intron = get_downstream_Intron($exon, $transcript);
  if($downstream_intron){
    return get_splice_sites($downstream_intron);
  }else{
    return undef;
  }
}

1;
