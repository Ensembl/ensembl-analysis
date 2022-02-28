=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Analysis::Tools::Filter::cDNA2GenomeTranscriptFilter

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Tools::Filter::cDNA2GenomeTranscriptFilter;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw);
use Data::Dumper;
use parent ('Bio::EnsEMBL::Analysis::Tools::Filter');


=head2 filter_results

 Arg [1]    : Arrayref of Bio::EnsEMBL::Transcript
 Description: Check if the alignment has a protein. We only want the models
              with protein and we replace C in the cigar string to M
 Returntype : Arrayref of Bio::EnsEMBL::Transcript
 Exceptions : None

=cut

sub filter_results{
  my ($self, $transcripts) = @_;

  throw('You should give an arrayref of Bio::EnsEMBL::Transcript') unless (ref($transcripts) eq 'ARRAY');
  my @good_matches;
  my %matches;

  foreach my $transcript (@$transcripts ) {
    next unless ($transcript->translation and $transcript->translation->length > 3);
    my $coverage  = $self->_get_transcript_coverage($transcript);
    my $percent_id  = $self->_get_transcript_percent_id($transcript);
    my $id = $self->_get_transcript_evidence_id($transcript);

    push @{$matches{$id}}, {
      transcript  => $transcript,
      coverage    => $coverage,
      percent_id  => $percent_id,
      num_exons   => scalar(@{$transcript->get_all_Exons}),
      is_spliced  => $self->_transcript_is_spliced($transcript),
    };
  }

  my %matches_sorted_by_coverage;

  foreach my $query_id ( keys %matches ) {
    @{$matches_sorted_by_coverage{$query_id}} =
        sort { $b->{coverage}   <=> $a->{coverage} or
               $b->{num_exons}  <=> $a->{num_exons} or
               $b->{percent_id} <=> $a->{percent_id} } @{$matches{$query_id}};

    my $max_coverage;
    my $perc_id_of_best;
    my $count = 1;
    my $splices_elsewhere = 0;
    my $best_has_been_seen = 0;

  #get the slice_name, start, and end of the best hit and store them to be used
  #in checking of any other good hit may be overlaping the best prediction
  # This is to avoid a match that spans over a region where two or more proteins of
  # a same family are close together and they get merged by a wrong "good quallity" alignment
    my $best_transcript = ${$matches_sorted_by_coverage{$query_id}}[0]->{transcript};
    my $best_start = $best_transcript->start;
    my $best_end = $best_transcript->end;
    my $best_slice;
    if (defined $best_transcript->slice) {
      $best_slice = $best_transcript->slice;
    }
    else {
      $best_slice = $best_transcript->start_Exon->seqname;
    }

    foreach my $hit (@{$matches_sorted_by_coverage{$query_id}}) {
      my $accept;
      my $label = $count++;
      my $transcript = $hit->{transcript};
      my $strand = $transcript->strand;
      my $coverage = $hit->{coverage};
      my $percent_id = $hit->{percent_id};
      my $is_spliced = $hit->{is_spliced};

      my $transcript_slice;
      if (defined $transcript->slice) {
        $transcript_slice = $transcript->slice;
      }
      else {
        $transcript_slice = $transcript->start_Exon->seqname;
      }

      $max_coverage = $coverage unless ($max_coverage);
      $perc_id_of_best = $percent_id unless ($perc_id_of_best);

      if ($count == 1) {
        $label = 'best_match';
        $splices_elsewhere = 1 if ($is_spliced);
      }
      elsif ($splices_elsewhere && ! $is_spliced) {
        $label = 'potential_processed_pseudogene';
      }

      if ($self->best_in_genome) {
        # we keep the hit with the best coverage...
        if ($coverage == $max_coverage &&
            # as long as it has coverage/percent_id above limits or...
            (($coverage >= $self->min_coverage &&
              $percent_id >= $self->min_percent)
             ||
             # ...if coverage is significanly greater than the
             # specified minimum, then we are willing to accept
             # hits that have a percent_id just below the specified
             # minimum
             ($coverage   >= (1 + 5/100) * $self->min_coverage &&
              $percent_id >= (1 - 3/100) * $self->min_percent))) {
          if ( $self->reject_processed_pseudos
               && $count > 1
               && $splices_elsewhere
               && ! $is_spliced) {
            $accept = 'NO';
          }
          # ... if one transcript with lower quality completely overlaps
          # the best one don't accept the lower quality one.
          elsif ($best_slice eq $transcript_slice &&
                 $best_start > $transcript->start &&
                 $best_end   < $transcript->end) {
            $accept = 'NO';
          }
          else {
            $accept = 'YES';
            push( @good_matches, $transcript);
          }
        }
        else{
          $accept = 'NO';
        }
      }
      else{
        # we keep anything which is within the 2% of the best score...
        if ($coverage >= (0.98 * $max_coverage) &&
            # as long as it has coverage/percent_id above limits or...
            (($coverage >= $self->min_coverage &&
              $percent_id >= $self->min_percent)
             ||
             # ...if coverage is significanly greater than the
             # specified minimum, then we are willing to accept
             # hits that have a percent_id just below the specified
             # minimum
             ($coverage   >= (1 + 5/100) * $self->min_coverage &&
              $percent_id >= (1 - 3/100) * $self->min_percent))) {


	        ############################################################
	        # non-best matches are kept only if they are not unspliced with the
	        # best match being spliced - otherwise they could be processed pseudogenes
          if ( $self->reject_processed_pseudos &&
               $count > 1 &&
               $splices_elsewhere &&
               ! $is_spliced) {
            $accept = 'NO';
          }
          else{
            $accept = 'YES';
            push(@good_matches, $transcript);
          }
        }
        else{
          $accept = 'NO';
        }
      }
    }
  }

  return \@good_matches;
}


=head2 _get_transcript_coverage

 Arg [1]    : Bio::EnsEMBL::Transcript
 Description: Retrieve the coverage of the alignment
 Returntype : Int
 Exceptions : None

=cut

sub _get_transcript_coverage{
  my ($self, $tran) = @_;

  if (@{$tran->get_all_supporting_features} and
      defined $tran->get_all_supporting_features->[0]->hcoverage) {
    return $tran->get_all_supporting_features->[0]->hcoverage;
  }
  else {
    return $tran->get_all_Exons->[0]->get_all_supporting_features->[0]->score;
  }
}


=head2 _get_transcript_percent_id

 Arg [1]    : Bio::EnsEMBL::Transcript
 Description: Retrieve the percentage of identity of the alignment
 Returntype : Int
 Exceptions : None

=cut

sub _get_transcript_percent_id{
  my ($self, $tran) = @_;

  my $sf;

  if (@{$tran->get_all_supporting_features}) {
    $sf = $tran->get_all_supporting_features->[0];
  } else {
    $sf = $tran->get_all_Exons->[0]->get_all_supporting_features->[0];
  }

  return $sf->percent_id;
}


=head2 _get_transcript_evidence_id

 Arg [1]    : Bio::EnsEMBL::Transcript
 Description: Retrieve the id of the supporting evidence
 Returntype : String
 Exceptions : None

=cut

sub _get_transcript_evidence_id{
  my ($self, $tran) = @_;

  my $sf;

  if (@{$tran->get_all_supporting_features}) {
    $sf = $tran->get_all_supporting_features->[0];
  }
  else {
    $sf = $tran->get_all_Exons->[0]->get_all_supporting_features->[0];
  }

  return $sf->hseqname;
}


=head2 _transcript_is_spliced

 Arg [1]    : Bio::EnsEMBL::Transcript
 Description: Check if the transcript is spliced and not just a framshift
              It returns 1 if at least one of the intron is longer than
              9 bases.
 Returntype : Boolean
 Exceptions : None

=cut

sub _transcript_is_spliced {
  my ($self, $tran) = @_;

  my @exons = sort { $a->start <=> $b->start } @{$tran->get_all_Exons};

  if (scalar(@exons) > 1) {
    # check that there are non "frameshift" introns
    for(my $i = 0; $i < @exons-1; $i++) {
      my $intron_len = $exons[$i+1]->start-$exons[$i]->end-1;
      return 1 if ($intron_len > 9);
    }
  }

  return 0;
}


sub _update_transcript_cigars {
  my ($self, $transcript) = @_;
  my $cigar_line = $transcript->get_all_supporting_features->[0]->cigar_string;
  foreach my $exon (@{$transcript->get_all_Exons}) {
    foreach my $sf (@{$exon->get_all_supporting_features}) {
      $sf->cigar_line($self->_update_cigar_string($sf->cigar_line));
    }
  }
  $transcript->get_all_supporting_features->[0]->cigar_string($self->_update_cigar_string($cigar_line));
}


=head2 _update_cigar_string

 Arg [1]    : String $cigar
 Description: Change Arg[1] from a cdna2genome cigar string which contains C where
              it match a codon. We only want M when it matches.
 Returntype : String
 Exceptions : Throw if Arg[1] is not set

=cut

sub _update_cigar_string {
  my ($self, $cigar_line) = @_;

  throw('You should give a cigar string to modify') unless ($cigar_line);
  while ($cigar_line =~ /((\d+)[MC](\d+)[MC])/gc) {
    my $regex = $1;
    my $new_cigar = $2+$3;
    $new_cigar .= 'M';
    $cigar_line =~ s/$regex/$new_cigar/;
  }
  $cigar_line =~ tr/C/M/;
  return $cigar_line;
}

1;
