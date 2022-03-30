# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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
=pod

=head1 NAME

polA-clipping.pl

=head1 DESCRIPTION

script to parse a fasta file and identify sequences with polyA/T tails/heads
these are then clipped and stored in a file 

Copied directly from sd3's /ensembl-pipeline/scripts/EST/new_polyA_clipping.pl
and made into a rough module 

Can also return which end of a sequence was clipped (head or tail) and how many bases were removed

clipping:
  the non-A/T sequences at the ends must be <=10bp (set by $buffer)  
  the polyA/T string must be >4bp to be removed
  it only clips polyA tails or polyT heads using a sliding window of 3 bp
  the clipping is recursive but only clips one end of a sequence
  the head/tail is only clipped if the polyA/T string is longer than the non-polyA/T string at the end of the sequence

perl new_polyA_clipping.pl sequences.fasta polyat_clipped.out

=cut

package Bio::EnsEMBL::Analysis::Tools::PolyAClipping;

use strict;
use warnings;
use Bio::Seq;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Exporter qw(import);

our @EXPORT_OK = qw(clip_if_necessary prepare_seq);


=head2 clip_if_necessary 

  Arg [1]    : Bio::Seq  
  Arg [2]    : Int - buffer size (optional)
  Arg [3]    : Int - window_size  (optional)
  Example    : my ($clipped, $clip_end, $num_bases_removed) = clip_if_necessary($cdna,$buffer,$window) 
  Description: Decides whether a cDNA needs to be clipped
  Returntype : String - the clipped sequence - or undef
  Exceptions : none

=cut
sub clip_if_necessary {
  my ($unprepared, $buffer, $window_size) = @_;
  my $clipped = undef;
  my $clipped_seq = "";
  my $clip_end;

  if (!defined $buffer) {
    #length of region allowed to mismatch at ends of sequence
    $buffer = 10;
  }
  if (!defined $window_size) {
    $window_size = 3;
  }
  my $end_region = $buffer + $window_size;

  my $unclipped = prepare_seq($unprepared);
  #print ">".$unclipped->display_id."\n".$unclipped->seq."\n";

  # # #
  #  If the sequence is long enough, clip it
  # # #
  my $seq = $unclipped->seq;
  if (length $seq >= $end_region) {
    #check for poly-A/T tails/heads - which should have mostly been removed...
    #decide whether looking for polyT head or polyA tail:
    my $head = substr($seq, 0, $end_region);
    my $tail = substr($seq, -$end_region);

    my $t_count = 0;
    my $a_count = 0;

    while ($head=~/([Tt]+)/g) { #will match greedily
      if (length $1 > $t_count) {
        $t_count = length $1;
      }
    }
    while ($tail=~/([Aa]+)/g) { #will match greedily
      if (length $1 > $a_count) {
        $a_count = length $1;
      }
    }

    #decide whether to clip heads, tails or unsure... set appropriate flag
    if (($a_count > $t_count) && ($a_count > 4)) { #call the subroutine to trim the polyA tail:
      $clipped_seq = clip_polya($seq, "tail", $buffer, $window_size, $end_region);
      $clip_end = "tail";

    } elsif (($a_count < $t_count) && ($t_count > 4)) { #call the subroutine to trim the polyT head:
      $clipped_seq = clip_polya($seq, "head", $buffer, $window_size, $end_region);
      $clip_end = "head";

    } else {
      if ($a_count > 4) { #only do for ones which appear to have a head/tail:
        #tied - not sure which to do -try both and choose afterwards
        my $clipped_head = clip_polya($seq, "head", $buffer, $window_size, $end_region);

        my $clipped_tail = clip_polya($seq, "tail", $buffer, $window_size, $end_region);

        #choose one which clipped the most:
        if (length $clipped_head < length $clipped_tail) {
          $clipped_seq = $clipped_head;
          $clip_end = "head";
        } elsif (length $clipped_tail < length $clipped_head) {
          $clipped_seq = $clipped_tail;
          $clip_end = "tail";
        } else { #still can't tell, leave as original seq
          $clipped_seq = $seq;
        }

      } else { # $a_count <= 4
        #not going to be clipped
        $clipped_seq = $seq;
      }
    }
  } else { # length $seq < $end_region so no need to clip
    $clipped_seq = $seq;
  }

  # # #
  # Check that we have a clipped sequence
  # # #
  my $num_bases_removed;
  if (defined $clipped_seq) {
    $clipped = $unclipped;
    eval {
      $clipped->seq($clipped_seq);
      $clipped->desc("");
    };
    #print "$unclipped->display_id: $seq\n$cdna->display_id: $clipped_seq\n\n";
    $num_bases_removed = (length $seq) - (length $clipped_seq);
  } else {
    #the entire sequence seems to be polyA/T tail/head
    warning("Sequence ".$unclipped->display_id." has been removed:\n$seq\n");
  }
  # We are now finished reading in all the sequences
  return ($clipped, $clip_end, $num_bases_removed);
}

=head2 prepare_seq

  Arg [1]    : Bio::Seq
  Example    : my $prepared = prepare_seq($unprepared)
  Description: Checks the display_id is OK
               Makes all sequence uppercase
  Returntype : Bio::Seq
  Exceptions : Can't read id

=cut
sub prepare_seq {
  my ($unclipped) = @_;
  my $prepared = new Bio::Seq;
    # # #
    # just want the id on the first line:
    # # #
    my $tmp = $unclipped->display_id;
    $tmp =~ s/\>//;
    my $id;
    if ($tmp =~ m/^gi\|\d+\|ref\|(N[MR]_.+)\|/) {   #RefSeq partially curated entries
      $id = $1;
    } elsif ($tmp =~ m/^gi\|\d+\|gb\|(\w+\.\d)\|/) {  #Genbank entries
      $id = $1;
    } elsif ($tmp =~ m/^gi\|\d+\|dbj\|(\w+\.\d)\|/) {  #DDBJ entries
      $id = $1;
    } elsif ($tmp =~ m/^gi\|\d+\|emb\|(\w+\.\d)\|/) {  #EMBL entries
      $id = $1;
    } elsif ($tmp =~ m/^[\w\d]+\s([\w\.\d]+)\s.+\n{1}/) {
      $id = $1;
    } elsif ($tmp =~m/^[\w\.\d]+\s.+\n{1}/) {
        #already in correct format - do nothing
      $id = $tmp;
    } elsif ($tmp =~ m/[\w\.\d]+/) {
      $id = $tmp;
    } else {
      throw("Cannot extract the input id from: \n".$unclipped->id."\nTry making your file so that the first line of each entry ".
            "only contains the id, eg: \n>BC000830.1\n");
    }

    # make a new clipped cdna
    $prepared->display_id($id);
    $prepared->accession_number($unclipped->accession_number);
    $prepared->version($unclipped->version);

    # change seq to uppercase (as pattern matching later on is case-sensitive)
    my $seq = $unclipped->seq;
    $seq =~ tr/a-z/A-Z/;
    $prepared->seq($seq);

    return $prepared;
}


=head2 clip_polya 

  Arg [1]    : String - the cDNA sequence 
  Arg [2]    : String - "head" or "tail" - end of sequence to be clipped
  Arg [3]    : Int - buffer
  Arg [4]    : Int - window
  Arg [5]    : Int - end region
  Example    : $clipped_seq = clip_polya($seq, "head"); 
  Description: Clips polyA or polyT sequences (case-sensitive) 
  Returntype : String - the clipped sequence - or undef
  Exceptions : none

=cut

sub clip_polya {
  my ($seq, $end, $buffer, $window_size, $end_region) = @_;  

  my @seq = split//, $seq;
  my $length = length $seq;
  my $a_count = 0;
  my $t_count = 0;

  my $clipped_seq = $seq;

  # # #
  # Count the number of consecutive A's or T's
  # # #
  if ($end eq "tail") {
    my $tail = substr($seq, -$end_region); 
    while ($tail =~ m/([Aa]+)/g){ #will match greedily
      if (length $1 > $a_count){
        $a_count = length $1;
      }
    }
  } elsif ($end eq "head") {
    my $head = substr($seq, 0, $end_region); 
    while ($head =~ m/([Tt]+)/g){ #will match greedily
      if (length $1 > $t_count){
        $t_count = length $1;
      }
    }
  }

  # # #
  # Treat the tail (high a_count) and head (high t_count) differently 
  # first look at tails
  # looking only for polyA tail - use moving window looking for strings of 2/3 A's
  # # #
  if ($a_count > 4) {
    #moving through seq starting from end - allow for buffer region:
    for (my $i = ($length - 1); $i > ($length - $buffer); $i--) {
      my $match = 0;
      for (my $j = $i; $j > ($i - $window_size); $j--) { #check a window
        if ($seq[$j] eq 'A') { 
          $match++;
        }
      }

      if ($match > 1) { #if (2+)/3 = A - looks like a polyA tail
        #in a polyA region - want to see how far this extends:
        my $pos = $i;
        while ($pos > ($window_size - 1)) {
          #move the window one position:
          if (ord($seq[$pos]) == 65) { #65 = decimal for 'A'
            $match--;
          }
          if (ord($seq[$pos - $window_size]) == 65) {
            $match++;
          }

          if ($match < 2) { 
            #at end of the polyA region:
            #find length of polyA region: polya_len = (seq length - non-tail length - post-tail buffer length)
            my $polya_len = $length - ($pos - ($window_size - 1)) - ($length  - ($i + 1)); 
            #test to see if polyA string > post polyA region:
            if ($polya_len > ($length - $i)) { 
              #we now want to clip end of sequence
              #identify last non-A in this window:
              my $len;
              for (my $j = ($pos - $window_size); $j <= $pos; $j++) {
                if (ord($seq[$j]) != 65) {
                  $len = ($j + 1);
                } else {
                  last;
                }
              }
              $clipped_seq = substr($seq, 0, $len);

              #now, it might be that the sequence look something like ....AAAAACGAAAAA
              #in which case the newly clipped ....AAAAACG can be reexamined
              $clipped_seq = clip_polya($clipped_seq, $end, $buffer, $window_size, $end_region);

            } # end of  if ($polya_len > ($length - $i)) {
            $pos = 0; #break out of while loop
          } # end of if ($match < 2) {
          $pos--;
        } # end of while ($pos > ($window_size - 1)) {

        if ($match > 1) {
          #then the while loop was finished whilst still on a polyA-tail
          $clipped_seq = undef;
        }

        last; #move onto a new sequence
      } # end of if ($match > 1) {
    } # close the forlopp [i]

  # # #
  # now look at the heads
  # looking only for polyT head - use moving window looking for strings of 2/3 T's
  # # #
  } elsif ($t_count > 4) {
    #moving through seq from front:
    for (my $i = 0; $i <= $buffer; $i++) {
      my $match = 0;
      for (my $j = $i; $j < ($i + $window_size); $j++) { #check a window 
        if ($seq[$j] eq 'T') { 
          $match++;
        }
      }

      if ($match > 1) { #if (2+)/3 = T - looks like a polyT head 
        my $pos = $i;
        #in a polyT region - want to see how far this extends:
        while ($pos < ($length - $window_size)) {
          #move the window one position
          if (ord($seq[$pos]) == 84) { #eq 'T'
            $match--;
          }
          if (ord($seq[$pos + $window_size]) == 84) {
            $match++;
          }

          if ($match < 2) { 
            #at end of polyT region:
            #find length of polyT region: polyt_len = (head length - pre-head buffer length)
            my $polyt_len = ($pos + ($window_size - 1)) - ($i - 1);  
            #test to see if polyT string > pre polyT region:
            if ($polyt_len > $i) { 
              #we now want to clip front of sequence:
              #identify first non-T in this window:
              my $len;
              for (my $j = ($pos + ($window_size)); $j >= $pos; $j--) {
                if (ord($seq[$j]) != 84) {
                  $len = $j;
                } else {
                  last;
                }
              }
              $clipped_seq = substr($seq, $len);

              #now, it might be that the sequence look something like TTTTTCGTTTTT...
              #in which case the newly clipped CGTTTTT... can be reexamined
              $clipped_seq = clip_polya($clipped_seq, $end, $buffer, $window_size, $end_region);
            } # end of if ($match < 2) {             
            $pos = $length; #break out of while loop
          } # end of while ($pos < ($length - $window_size)) {
          $pos++;
        } # end of while ($pos < ($length - $window_size)) { 

        if ($match > 1) {
          #then the while loop was finished whilst still on a polyA-tail
          $clipped_seq = undef;
        }
        last; #move onto a new sequence
      } # end of if ($match > 1) { 
    } # end of for (my $i = 0; $i <= $buffer; $i++) {
  } else {
    # no polyA or polyT
    # do nothing
  } # end of } elsif ($t_count > 4) {
  return $clipped_seq;
}
1;
