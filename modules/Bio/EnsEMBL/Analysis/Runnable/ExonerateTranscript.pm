=head1 LICENSE

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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript - 

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript->new(
								 -query_seqs     => \@q_seqs,
                                                             [or -query_file     => $q_file]   
								 -query_type     => 'dna',
								 -target_seqs    => \@t_seqs,
                                                             [or -target_file    => $t_file]   
                                                                 -program      => $exonerate,
								 -options        => $options,
								);

 $runnable->run; #create and fill Bio::Seq object
 my @transcripts = @{$runnable->output};
 
=head1 DESCRIPTION

This module handles a specific use of the Exonerate (G. Slater) program, namely 
the prediction of the transcript structure in a piece of genomic DNA by the alignment 
of a 'transcribed' sequence (EST, cDNA or protein). The results is a set of 
Bio::EnsEMBL::Transcript objects


=head1 METHODS

=cut

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::Runnable::BaseExonerate;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(calculate_exon_phases);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::BaseExonerate);


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($coverage_aligned) =
      rearrange([qw(
                    COVERAGE_BY_ALIGNED
                    )
                 ], @args);


  if (defined($coverage_aligned)) {
    $self->coverage_as_proportion_of_aligned_residues($coverage_aligned);
  } else {
    $self->coverage_as_proportion_of_aligned_residues(1);
  }


  return $self;
}



############################################################
#
# Analysis methods
#
############################################################


sub parse_results {
  my ($self, $fh) = @_;

  my %strand_lookup = ('+' => 1, '-' => -1, '.' => 1);

  # Each alignment will be stored as a transcript with 
  # exons and supporting features.  Initialise our
  # transcript.

  my @transcripts;

  # Parse output looking for lines beginning with 'RESULT:'.
  # Each line represents a distinct match to one sequence
  # containing multiple 'exons'.

 TRANSCRIPT:
  while (<$fh>){
    print STDERR $_ if $self->_verbose;

    next unless /^RESULT:/;

    chomp;

    my ($tag, $q_id, $q_start, $q_end, $q_strand, $t_id, $t_start, $t_end,
	$t_strand, $score, $perc_id, $q_length, $t_length, $gene_orientation,
	@align_components) = split;

    $t_strand = $strand_lookup{$t_strand};
    $q_strand = $strand_lookup{$q_strand};
    $gene_orientation = $strand_lookup{$gene_orientation};

    # Read vulgar information and extract exon regions.
    my ($exons, %flag) = $self->_parse_vulgar_block($t_start,
                                           $t_end,
                                           $t_strand,
                                           $t_length,
                                           $q_start, 
                                           $q_end,
                                           $q_strand,
                                           $q_length,
                                           \@align_components);

    # now we have extracted the exons and the coordinates are with 
    # reference to the forward strand of the query and target, we can 
    # use the gene_orienation to flip the strands if necessary
    if ($gene_orientation == -1 and $t_strand == 1) {
      $t_strand *= -1;
      $q_strand *= -1;
    }

    my $covered_count = 0;
    if ($self->coverage_as_proportion_of_aligned_residues) {
      foreach my $exon (@$exons) {
        foreach my $sf (@{$exon->{sf}}) {
          $covered_count += $sf->{query_end} - $sf->{query_start} + 1;
        }
      }
    } else {
      $covered_count = abs($q_end - $q_start);
    }

    my $coverage = sprintf("%.2f", 100 * $covered_count / $q_length);

    # Build FeaturePairs for each region of query aligned to a single
    # Exon.  Create a DnaDnaAlignFeature from these FeaturePairs and then
    # attach this to our Exon.
    my $transcript = Bio::EnsEMBL::Transcript->new();
    $transcript->{'accession'} = $q_id;
    $transcript->{'query_length'} = $q_length;
    $transcript->{'score'} = $score;

    my (@tran_feature_pairs,
        $cds_start_exon, $cds_start,
        $cds_end_exon, $cds_end);

    foreach my $proto_exon (@$exons){
      if ($flag{$proto_exon}) {
        my $link_attrib = Bio::EnsEMBL::Attribute->new
            ( -CODE => 'lost_frameshift',
              -NAME => 'lost_frameshift',
              -DESCRIPTION => 'Frameshift on the query sequence is lost in the target sequence',
              -VALUE => "internal frameshift exon $flag{$proto_exon}",
            ) ;
        $transcript->add_Attributes($link_attrib) ;
      }
 
      # Build our exon and set its key values.
      my $exon = Bio::EnsEMBL::Exon->new();
      
      $exon->seqname($t_id);
      $exon->start($proto_exon->{exon_start});
      $exon->end($proto_exon->{exon_end});
      $exon->phase($proto_exon->{phase});
      $exon->end_phase($proto_exon->{end_phase});
      $exon->strand($t_strand);
            
      my @exon_feature_pairs;
      foreach my $sf (@{$proto_exon->{sf}}){
        my $feature_pair = Bio::EnsEMBL::FeaturePair->new(-seqname    => $t_id,
                                                          -start      => $sf->{target_start},
                                                          -end        => $sf->{target_end},
                                                          -strand     => $t_strand,
                                                          -hseqname   => $q_id,
                                                          -hstart     => $sf->{query_start},
                                                          -hend       => $sf->{query_end},
                                                          -hstrand    => $q_strand,
                                                          -hcoverage  => $coverage,
                                                          -score      => $coverage,
                                                          -percent_id => $perc_id);


	push @exon_feature_pairs, $feature_pair;
        push @tran_feature_pairs, $feature_pair;
      }

      if (@exon_feature_pairs) {
        # Use our feature pairs for this exon to create a single 
        # supporting feature (with cigar line).
        my $supp_feature;
        
        eval{
          if ($self->query_type eq 'protein') {
            $supp_feature =
                Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@exon_feature_pairs, -align_type => 'ensembl');
          } else {
            $supp_feature = 
                Bio::EnsEMBL::DnaDnaAlignFeature->new(-features => \@exon_feature_pairs, -align_type => 'ensembl');
          }
        };
        if ($@){
          warning($@);
          next TRANSCRIPT;
        }      
        $exon->add_supporting_features($supp_feature);
      }      

      if (exists $proto_exon->{coding_start}) {
        if (not defined $cds_start_exon) {
          $cds_start_exon = $exon;
          $cds_start = $proto_exon->{coding_start};
        }
        $cds_end_exon = $exon;
        $cds_end = $proto_exon->{coding_end};
      }


      $transcript->add_Exon($exon);
    }

    # Create a single supporting feature for the whole transcript
    my $t_supp_feat;
    eval{
      if ($self->query_type eq 'protein') {
        $t_supp_feat =
            Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@tran_feature_pairs, -align_type => 'ensembl');
      } else {
        $t_supp_feat = 
            Bio::EnsEMBL::DnaDnaAlignFeature->new(-features => \@tran_feature_pairs, -align_type => 'ensembl');
      }
    };
    if ($@) {
      warning("Could not create Transcript supporting feature");
    } else {
      $transcript->add_supporting_features($t_supp_feat);
    }

    my @exons = @{$transcript->get_all_Exons};

    if (scalar(@exons)) {
      if (defined $cds_start_exon) {
        my $translation = Bio::EnsEMBL::Translation->new();
        $translation->start_Exon($cds_start_exon);
        $translation->start($cds_start);      
        $translation->end_Exon($cds_end_exon);
        $translation->end($cds_end);

        $transcript->translation($translation);
      }
      
      calculate_exon_phases($transcript, 0);


      $transcript->analysis($self->analysis);
      push @transcripts, [$score, $transcript];
    }
  }

  @transcripts = sort { $b->[0] <=> $a->[0] } @transcripts;
  @transcripts = map { $_->[1] } @transcripts;

  return \@transcripts;
}

sub _parse_vulgar_block {
  my ($self, 
      $target_start, $target_end, $target_strand, $target_length,
      $query_start, $query_end,  $query_strand, $query_length,
      $vulgar_components) = @_;

  # This method works along the length of a vulgar line 
  # exon-by-exon.  Matches that comprise an exon are 
  # grouped and an array of 'proto-exons' is returned.
  # Coordinates from the vulgar line are extrapolated 
  # to actual genomic/query coordinates.

  my @exons;
  my $exon_number = 0;
  my %flag ;

  # We sometimes need to increment all our start coordinates. Exonerate 
  # has a coordinate scheme that counts _between_ nucleotides at the start.
  # However, for reverse strand matches 
  
  my ($query_in_forward_coords, $target_in_forward_coords);
  my ($cumulative_query_coord, $cumulative_target_coord);

  if ($target_start > $target_end) {
    warn("For target, start and end are in thew wrong order for a reverse strand match")
        if $target_strand != -1;
    $cumulative_target_coord = $target_start;
    $target_in_forward_coords = 1;
  } else {
    $cumulative_target_coord = $target_start + 1;
    $target_in_forward_coords = 0;
  }
  if ($query_start > $query_end) {
    warn("For query, start and end are in thew wrong order for a reverse strand match")
        if $query_strand != -1;
    $cumulative_query_coord = $query_start;
    $query_in_forward_coords = 1;
  } else {
    $cumulative_query_coord = $query_start + 1;
    $query_in_forward_coords = 0;
  }


  while (@$vulgar_components){
    throw("Something funny has happened to the input vulgar string." .
		 "  Expecting components in multiples of three, but only have [" .
		 scalar @$vulgar_components . "] items left to process.")
      unless scalar @$vulgar_components >= 3;

    my $type                = shift @$vulgar_components;
    my $query_match_length  = shift @$vulgar_components;
    my $target_match_length = shift @$vulgar_components;

    throw("Vulgar string does not start with a match.  Was not " . 
		 "expecting this. (Have type $type)")
      if (scalar @exons == 0) && $type ne 'M' && $type ne 'S' && $type ne 'C';

    if ($type eq 'M' or $type eq 'S' or $type eq 'C'){
      my %hash = (type => $type);

      if ($target_strand == -1) {
        if ($target_in_forward_coords) {
          $hash{target_start} = $cumulative_target_coord - ($target_match_length - 1);
          $hash{target_end}   = $cumulative_target_coord;
        } else {
          $hash{target_end}   = $target_length - ($cumulative_target_coord - 1);
          $hash{target_start} = $hash{target_end} - ($target_match_length - 1);
        }
      } else {
        $hash{target_start} = $cumulative_target_coord;
        $hash{target_end}   = $cumulative_target_coord + ($target_match_length - 1);
      }

      if ($query_strand == -1) {
        if ($query_in_forward_coords) {
          $hash{query_start} = $cumulative_query_coord - ($query_match_length - 1);
          $hash{query_end}   = $cumulative_query_coord;
        } else {
          $hash{query_end}   = $query_length - ($cumulative_query_coord - 1);
          $hash{query_start} = $hash{query_end} - ($query_match_length - 1);
        }
      } else {
        $hash{query_start} = $cumulative_query_coord;
        $hash{query_end}   = $cumulative_query_coord + ($query_match_length - 1);
      }

      # there is nothing to add if this is the last state of the exon
      $exons[$exon_number]->{gap_end}   = 0;
      push @{$exons[$exon_number]->{sf}}, \%hash;
    }
    elsif ($type eq "G") {
      if (exists($exons[$exon_number]->{sf})) {
        # this is the gap in the middle of an exon, or at the end. Assume it is 
        # at the end, and then reset if we see another match state in this exon
        $exons[$exon_number]->{gap_end}   = $target_match_length;
      } else {
        # this is a gap at the start of an exon; 
        $exons[$exon_number]->{gap_start} = $target_match_length;
      }
    }
    elsif ($type eq "I" or
           $type eq "F") {


      # if the intron or frameshift is only on the query sequence, 
      # the target sequence will have only one exon for the region
      # instead of two on the query sequence
      if ($target_match_length > 0) {

      # in protein mode, any insertion on the genomic side should be treated as 
      # an intron to ensure that the result translates. However, we allow for
      # codon insertions in the genomic sequence with respect to the protein. 
      # This introduces the possibility of in-frame stops, but I don't
      # think "introning over" these insertions is appropriate here. 

      # if we see a gap/intron immediately after an intron, the current exon is "empty"
        if ($exons[$exon_number]) {
          $exon_number++;
        }
      } elsif ($type eq "I") {
          throw("Really odd case here, with lost intron !!\n") ;
      } else {
          $flag{$exons[$exon_number]} = $exon_number ;
      }
    }

    if ($target_in_forward_coords and $target_strand == -1) {
      $cumulative_target_coord -= $target_match_length;
    } else {
      $cumulative_target_coord += $target_match_length;
    }
    if ($query_in_forward_coords and $query_strand == -1) {
      $cumulative_query_coord  -= $query_match_length;
    }
    else {
      $cumulative_query_coord  += $query_match_length;
    }

  }

  for(my $i = 0; $i < @exons; $i++) {
    my $ex = $exons[$i];
    my @ex_sf = @{$ex->{sf}};

    if ($target_strand == -1) {
      $ex->{exon_start} = $ex_sf[-1]->{target_start};
      $ex->{exon_end}   = $ex_sf[0]->{target_end};

      if (exists $ex->{gap_start}) {
        $ex->{exon_end} += $ex->{gap_start};
      }
      if (exists $ex->{gap_end}) {
        $ex->{exon_start} -= $ex->{gap_end};
      }

    } else {
      $ex->{exon_start} = $ex_sf[0]->{target_start};
      $ex->{exon_end}   = $ex_sf[-1]->{target_end};

      if (exists $ex->{gap_start}) {
        $ex->{exon_start} -= $ex->{gap_start};
      }
      if (exists $ex->{gap_end}) {
        $ex->{exon_end} += $ex->{gap_end};
      }
    }

    # split codons are a pain. If the query is dna, we must be in the 
    # cdna2genome model so they must be part of the supporting feature. 
    # If query is protein, they need to be removed from the supporting feature

    if ($self->query_type eq 'dna') {
      map { $_->{type} = 'C' if $_->{type} eq 'S' } @ex_sf;
      if (my @cod = grep { $_->{type} eq 'C' } @ex_sf) {
        # at least part of this exon is coding

        @cod = sort { $a->{target_start} <=> $b->{target_start} } @cod;
        my $cod_start = $cod[0]->{target_start};
        my $cod_end   = $cod[-1]->{target_end};
        my $cod_len   = $cod_end - $cod_start + 1;

        if ($target_strand == -1) {
          $ex->{coding_start} = $ex->{exon_end} - $cod_end + 1;
        } else {
          $ex->{coding_start} = $cod_start - $ex->{exon_start} + 1;
        }

        $ex->{coding_end}   = $ex->{coding_start} + $cod_len - 1;
      }
      # merge together abutting ungapped features
      my @nr_sf;
      foreach my $sf (@ex_sf) {
        my $merged = 0;
        if (@nr_sf and
            $sf->{type} eq 'C' and
            $nr_sf[-1]->{type} eq 'C' and
            $sf->{query_start} == $nr_sf[-1]->{query_end} + 1) {
          if ($target_strand == -1) {
            if ($sf->{target_end} == $nr_sf[-1]->{target_start} - 1) {
              $nr_sf[-1]->{target_start} = $sf->{target_start};
              $nr_sf[-1]->{query_end}  = $sf->{query_end};
              $merged = 1;
            }
          } else {
            if ($sf->{target_start} == $nr_sf[-1]->{target_end} + 1) {
              $nr_sf[-1]->{target_end} = $sf->{target_end};
              $nr_sf[-1]->{query_end}  = $sf->{query_end};
              $merged = 1;
            }
          }         
        }
        if (not $merged) {
          push @nr_sf, $sf;
        }
      }
      @ex_sf = @nr_sf;
    } else {
      # query type is protein
      @ex_sf = grep { $_->{type} ne 'S' } @ex_sf;
      $ex->{coding_start} = 1;
      $ex->{coding_end} = $ex->{exon_end} - $ex->{exon_start} + 1;
    }

    $ex->{sf} = \@ex_sf;
  }

  return \@exons, %flag ;
}




############################################################
#
# get/set methods
#
############################################################

sub coverage_as_proportion_of_aligned_residues {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_coverage_aligned} = $val;
  }
  return $self->{_coverage_aligned};
}


1;

