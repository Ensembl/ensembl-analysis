# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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
#
# EnsEMBL module for TranscriptChecker
#
# Cared for by Steve Searle <searle@sanger.ac.uk>
#

# POD documentation - main docs before the code

=pod

=head1 NAME

Bio::EnsEMBL::Test::TranscriptChecker


=head1 SYNOPSIS
Module to check the validity of a transcript

=head1 DESCRIPTION
Performs various checks on a Bio::EnsEMBL::Transcript object. These include:
  1. Transcripts with no exons
  2. Transcripts with more than a specified maximum number of exons
  3. Exon rank ordering is correct
  4. Exons don't overlap or contain other exons
  5. Transcript translates
  6. Translation is longer than a specified minimum length
  7. Short or long introns
  8. Introns with non GT .. AG splice sites (warnings)
  9. Short or long exons
 10. ATG immediately after 5' UTR (warning)
 11. There are supporting features for every exon (warning)
 12. The translation is not all X
 13. Transcript with no transcript_supporting_features

This class does not use stable_ids but instead dbIDs because stable_ids are
not set until after the gene build.

=head1 CONTACT

  Steve Searle <searle@sanger.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package TranscriptChecker;
use vars qw(@ISA);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use strict;

use Checker;
use Bio::EnsEMBL::Intron;
use Bio::Seq;


@ISA = qw( Checker );


sub new {
  my ($class, @args) = @_;

  my $self = bless {},$class;

  my ( $transcript, $minshortintronlen, $maxshortintronlen, $minlongintronlen,
       $minshortexonlen, $maxshortexonlen, $minlongexonlen, $maxexonstranscript,
       $mintranslationlen, $ignorewarnings, $slice, $adaptor ) = rearrange
	 ( [ qw ( TRANSCRIPT
		  MINSHORTINTRONLEN
		  MAXSHORTINTRONLEN
		  MINLONGINTRONLEN
		  MINSHORTEXONLEN
		  MAXSHORTEXONLEN
		  MINLONGEXONLEN
                  MAXEXONTRANSCRIPT
                  MINTRANSLATIONLEN
                  IGNOREWARNINGS
                  SLICE
                  ADAPTOR
	      )], @args );


  if( !defined $transcript ) {
    $self->throw("Transcript must be set in new for TranscriptChecker")
  }
  $self->transcript($transcript);

  if( defined $minshortintronlen ) {
    $self->minshortintronlen( $minshortintronlen );
  } else {
    $self->minshortintronlen(3);
  }
  if( defined $maxshortintronlen ) {
    $self->maxshortintronlen( $maxshortintronlen );
  } else {
    $self->maxshortintronlen(50);
  }
  if( defined $minlongintronlen ) {
    $self->minlongintronlen( $minlongintronlen );
  } else {
    $self->minlongintronlen(100000);
  }

  if( defined $minshortexonlen ) {
    $self->minshortexonlen( $minshortexonlen );
  } else {
    $self->minshortexonlen(3);
  }
  if( defined $maxshortexonlen ) {
    $self->maxshortexonlen( $maxshortexonlen );
  } else {
    $self->maxshortexonlen(10);
  }
  if( defined $minlongexonlen ) {
    $self->minlongexonlen( $minlongexonlen );
  } else {
    $self->minlongexonlen(50000);
  }

  if( defined $mintranslationlen ) {
    $self->mintranslationlen( $mintranslationlen );
  } else {
    $self->mintranslationlen(10);
  }

  if( defined $maxexonstranscript) {
    $self->maxexonstranscript( $maxexonstranscript );
  } else {
    $self->maxexonstranscript(150);
  }
  if( defined $ignorewarnings) { $self->ignorewarnings($ignorewarnings); }
  if( defined $slice) { $self->slice($slice); }
  if( defined $adaptor ) { $self->adaptor( $adaptor )}

  $self->{_errors} = [];
  $self->{_warnings} = [];

  return $self;
}


sub mintranslationlen {
  my ( $self, $arg ) = @_;
  if (defined $arg) {
    $self->{_mintranslationlen} = $arg;
  }
  return $self->{_mintranslationlen};
}

sub maxshortintronlen {
  my ( $self, $arg ) = @_;
  if (defined $arg) {
    $self->{_maxshortintronlen} = $arg;
  }
  return $self->{_maxshortintronlen};
}

sub minshortintronlen {
  my ( $self, $arg ) = @_;
  if (defined $arg) {
    $self->{_minshortintronlen} = $arg;
  }
  return $self->{_minshortintronlen};
}

sub minlongintronlen {
  my ( $self, $arg ) = @_;
  if (defined $arg) {
    $self->{_minlongintronlen} = $arg;
  }
  return $self->{_minlongintronlen};
}

sub maxexonstranscript {
  my ( $self, $arg ) = @_;
  if (defined $arg) {
    $self->{_maxexonstranscript} = $arg;
  }
  return $self->{_maxexonstranscript};
}

sub maxshortexonlen {
  my ( $self, $arg ) = @_;
  if (defined $arg) {
    $self->{_maxshortexonlen} = $arg;
  }
  return $self->{_maxshortexonlen};
}

sub minshortexonlen {
  my ( $self, $arg ) = @_;
  if (defined $arg) {
    $self->{_minshortexonlen} = $arg;
  }
  return $self->{_minshortexonlen};
}

sub minlongexonlen {
  my ( $self, $arg ) = @_;
  if (defined $arg) {
    $self->{_minlongexonlen} = $arg;
  }
  return $self->{_minlongexonlen};
}



=head2 transcript

 Title   : transcript
 Usage   : $obj->transcript($newval)
 Function:
 Returns : value of transcript
 Args    : newvalue (optional)

=cut

sub transcript {
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      if (!($value->isa('Bio::EnsEMBL::Transcript'))) {
         $self->throw("transcript passed a non Bio::EnsEMBL::Transcript object\n");
      }
      $self->{_transcript} = $value;
    }
    return $self->{_transcript};

}

sub slice {
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      if (!($value->isa('Bio::EnsEMBL::Slice'))) {
         $self->throw("vc passed a non Bio::EnsEMBL::Slice object\n");
      }
      $self->{_slice} = $value;
    }
    return $self->{_slice};

}

sub adaptor {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{_adaptor} = $arg;
  }
  return $self->{_adaptor};
}


sub output {
  my $self = shift;

  my $transcript = $self->transcript;

  print "\n===\n";
  print "Transcript " . $transcript->dbID . " " .
        (defined($transcript->stable_id) ? $transcript->stable_id : "") .  "\n";

  my $vcoffset = 0;
  if ($self->slice) {
    $vcoffset = $self->slice->start - 1;
  }

  my $translation = $transcript->translation;

  printf "       %9s %9s %9s %9s %9s %9s %9s\n", "dbID","Start","End","Strand","Length","Phase","EndPhase";
  foreach my $exon (@{$transcript->get_all_Exons()}) {
    printf "  Exon %9d %9d %9d %9d %9d %9d %9d", $exon->dbID, $exon->start + $vcoffset, $exon->end + $vcoffset,
           $exon->strand, ($exon->end-$exon->start+1), $exon->phase, $exon->end_phase;
    if (defined($translation)) {
      if ($exon == $translation->start_Exon) {
        print " S (" . $translation->start . " (" .
                 (($exon->strand == 1 ? $exon->start+$translation->start-1 : $exon->end-$translation->start+1) + $vcoffset) ."))";
      }
      if ($exon == $translation->end_Exon) {
        print " E (" . $translation->end . " (" .
                 (($exon->strand == 1 ? $exon->start+$translation->end-1 : $exon->end-$translation->end+1) + $vcoffset) ."))";
      }
    }
    print "\n";
  }

  if ($translation) {
    my $pepstr = $transcript->translate->seq;
    $pepstr =~ s/(.{0,80})/  $1\n/g;
    print "\nTranslation (" . $translation->dbID . "):\n";
    print $pepstr;
  }

  $self->SUPER::output;

  if ($self->adaptor) {
    $self->print_raw_data;
  }
}

sub print_raw_data {
#  my $self = shift;
#
#  my $dbId = $self->transcript->dbID;
#  my $sgp_type = $self->adaptor->static_golden_path_type;
#
#  my $q = qq(
#    SELECT e.exon_id,
#           if(sgp.raw_ori=1,(e.seq_start-sgp.raw_start+sgp.chr_start),
#              (sgp.chr_start+sgp.raw_end-e.seq_end)) as start,
#           if(sgp.raw_ori=1,(e.seq_end-sgp.raw_start+sgp.chr_start),
#              (sgp.chr_start+sgp.raw_end-e.seq_start)) as end,
#           if (sgp.raw_ori=1,e.strand,(-e.strand)) as strand,
#           sgp.chr_name,
#           abs(e.seq_end-e.seq_start)+1 as length,
#           et.rank,
#           if(e.exon_id=tl.start_exon_id, (concat(tl.seq_start," (start)",
#              if(e.exon_id=tl.end_exon_id,(concat(" ",tl.seq_end," (end)")),
#              ("")))),if (e.exon_id=tl.end_exon_id,
#              (concat(tl.seq_end," (end)")),(""))) as transcoord,
#           if(e.sticky_rank>1,(concat("sticky (rank = ",e.sticky_rank,")")),
#              ("")) as sticky
#     FROM  translation tl, exon e, transcript tr, exon_transcript et,
#           static_golden_path sgp
#     WHERE e.exon_id=et.exon_id AND
#           et.transcript_id=tr.transcript_id AND
#           sgp.raw_id=e.contig_id AND sgp.type = '$sgp_type' AND
#           tr.transcript_id = $dbId AND
#           tr.translation_id=tl.translation_id
#     ORDER BY et.rank
#  );
#  my $sth = $self->adaptor->prepare($q) || $self->throw("can't prepare: $q");
#  my $res = $sth->execute || $self->throw("can't execute: $q");
#
#  print "\nTranscript data from SQL query:\n";
#  printf "%-9s: %-9s %-9s %-15s %-6s %-6s %-4s %-16s %-11s\n",
#         "Exon ID", "Start", "End", "Chromosome", "Strand", "Length",
#         "Rank", "Translation_Info", "Sticky_Info";
#  while( my ($id, $start, $end, $strand, $chrname, $length,
#             $rank, $transstr, $stickystr)= $sth->fetchrow_array) {
#    printf "%9d: %9d %9d %-15s %6d %6d %4d %16s %11s\n",$id, $start, $end,
#           $chrname, $strand, $length, $rank, $transstr, $stickystr;
#  }
}


sub check {
  my $self = shift;

  my $transcript = $self->transcript;
  #print "Transcript " . $transcript->dbID . " ";
  my @exons = @{$transcript->get_all_Exons()};

  my $numexon = scalar(@exons);

  if ($numexon == 0) {
    $self->add_Error("No exons\n", 'noexons');
    return;

  } elsif ($numexon > $self->maxexonstranscript) {
    $self->add_Error("Unusually large number of exons (" .
                     $numexon . ")\n", 'manyexons');
  }


# Sorting is needed because rank is flacky
# This allows other tests to be performed even if ranks are wrong
  my @sortedexons;

  my $transstart;
  my $transend;
  my $vcoffset = 0;
  if ($self->slice) {
    $vcoffset = $self->slice->start - 1;
  }

  if ($exons[0]->strand == 1) {
    @sortedexons = sort {$a->start <=> $b->start} @exons ;
    $transstart = ($sortedexons[0]->start+$vcoffset);
    $transend   = ($sortedexons[$#sortedexons]->end+$vcoffset);
  } else {
    @sortedexons = sort {$b->start <=> $a->start} @exons;
    $transstart = ($sortedexons[$#sortedexons]->start+$vcoffset);
    $transend   = ($sortedexons[0]->end+$vcoffset);
  }

# Check for rank errors
  my $exnum = 0;
  EXON: foreach my $exon (@sortedexons) {
    if ($exon != $exons[$exnum++]) {
      $self->add_Error("Incorrect exon ranks (first at exon " .
                       $exons[$exnum++]->dbID . ")\n".
                       "NOTE: Further error checking (except translations) " .
                       "done with resorted exons.\n" ,'exonrank');
      last EXON;
    }
  }


  $self->check_Structure(\@sortedexons, $transstart, $transend);

  $self->check_Phases;

  if (defined($self->transcript->translation)) {
    $self->check_Translation;
    $self->check_UTRs(\@sortedexons);
  }

  $self->check_Supporting_Evidence(\@sortedexons);

  # check transcript supporting evidence
  if (!scalar(@{$transcript->get_all_supporting_features})) {
    $self->add_Error("No supporting evidence for transcript ".$transcript->dbID ."\n",'nosupport');
  }
}

sub check_Phases {
  my $self = shift;

  if (defined($self->transcript->translation)) {
    my @exons = @{$self->transcript->get_all_translateable_Exons};
    my $prev_exon = $exons[0];
    my $start_phase = $prev_exon->phase;

    if ($start_phase == -1) {
      $start_phase = 0;
    }

    my $zero_phase_count = 0;
    if ($start_phase == 0) {
      $zero_phase_count++;
    }

    my @calc_phases;
    my $phase = $start_phase;
    for (my $i=0; $i<scalar(@exons);$i++) {
      push @calc_phases, $phase;
      $phase = (($exons[$i]->length + $phase) % 3);
    }

    for (my $i=1; $i<scalar(@exons);$i++) {
      my $exon=$exons[$i];
      if ($exon->phase == 0) {
        $zero_phase_count++;
      }
      if ($exon->phase != $prev_exon->end_phase) {
        $self->add_Error("EndPhase/Phase mismatch (" . $prev_exon->end_phase . " and " . $exon->phase .
                         ") between exons " .  $prev_exon->dbID . " and " . $exon->dbID .
                         "\n", 'endphase_phase');
      }

      if ($exon->phase != $calc_phases[$i]) {
        $self->add_Error("Phase wrong for " . $exon->dbID .  " is " . $exon->phase . " should be " .
                         $calc_phases[$i] . "\n", 'startphase');
      }

      my $calc_endphase = (($exon->length + $calc_phases[$i]) % 3);
      #print $self->transcript->stable_id . " " . $exon->end_phase . " " . $calc_endphase . "\n";
      if ($exon->end_phase != $calc_endphase && $i != $#exons) {
        $self->add_Error("EndPhase wrong for " . $exon->dbID .  " is " . $exon->end_phase . " should be " .
                         $calc_endphase . "\n", 'endphase');
      }
      $prev_exon = $exon;
    }
    if ($zero_phase_count/scalar(@exons) < 0.20 && scalar(@exons) > 8) {
      $self->add_Warning("Suspiciously high number of non phase zero exon (num phase zero = $zero_phase_count nexon = " . scalar(@exons).")\n");
    }
  }
}

sub check_Translation {
  my $self = shift;
  my $pepseq = undef;

  return if (!defined($self->transcript->translation));

  eval {
    $pepseq = $self->transcript->translate;
  };
  if (defined($pepseq)) {
    my $pepseqstr = $pepseq->seq;
    my $peplen = length($pepseqstr);
    # print "Pep seq = $pepseqstr\n";
    if ($pepseqstr =~ /\*/) {
      $self->add_Error("Translation failed - Translation contains stop codons\n",'transstop');
      $self->get_initial_frame;
      return 1;
    } elsif ($peplen == 0) {
      $self->add_Error("Translation failed - Translation has zero length\n",'transzerolen');
      return 1;
    } elsif ($pepseqstr =~ /^X*$/) {
      $self->add_Error("Translation failed - Translation is all X (probably on sequence gap)\n",'transallx');
      return 1;
    } elsif ($peplen < $self->mintranslationlen) {
      $self->add_Error("Short (" . $peplen . " residue) translation\n",'transminlen');
    }
    return 0;
  } else {
    $self->add_Error("Translation failed.\n",'transfail');
    return 1;
  }
}

sub check_Structure {
  my $self = shift;
  my $sortedexons = shift;
  my $transstart = shift;
  my $transend = shift;

  my $prev_exon = undef;
  my $totalexonlen  = 0;
  foreach my $exon (@$sortedexons) {
    my $exlen = $exon->length;
    $totalexonlen+=$exlen;

    if ($exlen >= $self->minshortexonlen &&
        $exlen <= $self->maxshortexonlen) {
      $self->add_Error("Short exon (" . $exlen .
                       " bases) for exon " . $exon->dbID . "\n", 'shortexon');
    } elsif ($exon->length >= $self->minlongexonlen) {
      $self->add_Error("Long exon (" . $exlen .
                       " bases) for exon " . $exon->dbID . "\n", 'longexon');
    }
    if (defined $prev_exon) {
      if ($exon->strand != $prev_exon->strand) {
        $self->add_Error("Exons on different strands\n",'mixedstrands');

      } elsif (($exon->strand == 1 && $exon->start < $prev_exon->end) ||
               ($exon->strand == -1 && $exon->end > $prev_exon->start)) {
        $self->add_Error("Incorrect exon ordering (maybe duplicate exons or embedded exons)\n", 'exonorder');

      } else {
        $self->check_Intron($prev_exon, $exon);
      }
    }
    $prev_exon = $exon;
  }
  my $exondensity = $totalexonlen/($transend-$transstart+1);
  if (scalar(@$sortedexons) > 1) {
    if ($exondensity > 0.8) {
      $self->add_Warning("Unusually high exon density (exon len/genomic len) (value $exondensity)\n");
    } elsif ($exondensity < 0.005) {
      $self->add_Warning("Unusually low exon density (exon len/genomic len) (value $exondensity)\n");
    }
  }
}

#NOTE: This method doesn't use the UTR methods in Transcript because:
#  1. They rely on stable_ids
#  2. They will fetch unnecessary (for this purpose) sequence
#
sub check_UTRs {
  my $self = shift;
  my $exons = shift;

  my $translation = $self->transcript->translation;

  return if (!defined($translation));
  my $rank = 0;
  my $trans_start_exon = $translation->start_Exon;
  my $trans_end_exon = $translation->end_Exon;
  my $foundstart = 0;
  my $foundend = 0;

  my $start_is_atg = 0;
  my $stop_is_term = 0;

  EXON: foreach my $exon (@$exons) {
    if ($exon == $trans_start_exon) {
      $foundstart = 1;
      if ($translation->start > 3 || $rank > 0) {
        my $startcodon = substr($exon->seq->seq,$translation->start-1,3);
        if ($startcodon ne "ATG") {
          $self->add_Warning("No ATG at five prime of transcript CDS with UTR (has $startcodon)\n");
        } else {
          $start_is_atg = 1;
        }

        my $len = $exon->length;
        my $end_phase = $exon->end_phase;
        if ($exon == $trans_end_exon) { $len = $translation->end; }
        if ($end_phase == -1) { $end_phase = 0; }

        if (($len - $translation->start + 1 + ((3-$end_phase)%3)) % 3) {
          $self->add_Warning("Translation start not on codon boundary (or end_phase error on first coding exon)\n");
          if ($start_is_atg) {
            $self->add_Error("ATG as first three bases in case where first codon supposedly incomplete\n");
          }
        }
      }
    }
    if ($exon == $trans_end_exon) {
      $foundend = 1;
      if ($translation->end < ($exon->length-2) || $rank < scalar(@$exons) - 1) {
        my $stopcodon = substr($exon->seq->seq,$translation->end-3,3);
        if ($stopcodon !~ /TAA|TAG|TGA/) {
          $self->add_Warning("No TAA, TAG or TGA at three prime end of transcript CDS with UTR (has $stopcodon)\n");
        } else {
          $stop_is_term = 1;
        }
      }
    }
    $rank++;
  }
  if (!$foundstart) {
    $self->add_Error("Didn't find translation->start_exon (" .
                     $trans_start_exon->dbID .  ")\n");
  }
  if (!$foundend) {
    $self->add_Error("Didn't find translation->end_exon (" .
                     $trans_end_exon->dbID .  ")\n");
  }

  my $cdslen = 0;
  foreach my $exon (@{$self->transcript->get_all_translateable_Exons}) {
    $cdslen += $exon->length;
  }
  if ($cdslen%3) {
    if ($start_is_atg && $stop_is_term) {
      $self->add_Error("CDS length not multiple of 3 in transcript with ATG->Stop\n");
    } else {
      $self->add_Warning("CDS length not multiple of 3\n");
    }
  }
}

sub check_Supporting_Evidence {
  my $self = shift;
  my $exons = shift;

  EXON: foreach my $exon (@$exons) {
    if (!scalar(@{$exon->get_all_supporting_features})) {
      # $self->add_Error("No supporting evidence for exon ".$exon->dbID ."\n",'nosupport');
    }
  }
}

sub check_Intron {
  my $self = shift;
  my $prev_exon = shift || $self->throw("Prev_exon must be passed");
  my $exon = shift || $self->throw("Exon must be passed");

  if ($prev_exon->overlaps($exon)) {
    $self->add_Error("Overlapping exons in transcript for exons " . $prev_exon->dbID .
                     " and " . $exon->dbID . "\n", 'overlapexon');
    return;
  }

  my $intron;
  if (Bio::EnsEMBL::Intron->can('upstream_Exon')) {
    $intron = new Bio::EnsEMBL::Intron;
    $intron->upstream_Exon($prev_exon);
    $intron->downstream_Exon($exon);
  } else {
    $intron = new Bio::EnsEMBL::Intron($prev_exon, $exon);
  }

  my $intlen = $intron->length;

  if ($intlen >= $self->minshortintronlen &&
      $intlen <= $self->maxshortintronlen) {
    $self->add_Error("Short intron (" . $intlen .
                     " bases) between exons " . $prev_exon->dbID .
                     " and " . $exon->dbID . "\n", 'shortintron');
  } elsif ($intlen >= $self->minlongintronlen) {
    $self->add_Error("Long intron (" . $intlen . ") between exons " .
                     $prev_exon->dbID . " and " . $exon->dbID .
                     "\n", 'longintron');
  }

  # if ($intlen >= $self->minshortintronlen) {
  if ($intlen >= 2) {
    my $intseq = $intron->seq;
    if (ref($intseq)) {
      $intseq = $intseq->seq;
    }
    #print "Checking intron\n";
    if (substr($intseq,0,2) ne "GT") {
      $self->add_Warning("Non consensus 5' intron splice site sequence (".
                         substr($intseq,0,2) . ") after exon " .
                         $prev_exon->dbID . " (intron length = " . $intlen .
                         ")\n", 'noncons5');
    }
    if (substr($intseq,-2,2) ne "AG") {
      $self->add_Warning("Non consensus 3' intron splice site sequence (".
                         substr($intseq,-2,2) .
                         ") before exon " .  $exon->dbID . " (intron length = ".
                         $intlen . ")\n", 'noncons3');
    }
  }
}

sub get_initial_frame {
  my ($self) = @_;

  my $trans = $self->transcript;

  $trans->sort;

  # print " Trying to find translating frame\n";
  my $cdna = new Bio::Seq(-seq => $trans->translateable_seq);

  foreach my $frame (0 .. 2) {
    my $pepseq = $cdna->translate(undef, undef, $frame)->seq;
    chop $pepseq;
    if ($pepseq !~ /\*/) {
      $self->add_Error("Frame $frame does translate\n",'wrongframe');
      return $frame;
    }
  }
#  print "Failed to translate transcript in any phase\n";
#  foreach my $frame (0 .. 2) {
#    my $pepseq = $cdna->translate(undef, undef, $frame)->seq;
#    chop $pepseq;
#    print "Pepseq frame $frame = " . $pepseq . "\n";
#  }
  return 0;
}
