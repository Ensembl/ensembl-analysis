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

Bio::EnsEMBL::Analysis::Runnable::ExonerateAlignFeature - 

=head1 SYNOPSIS

  my $runnable = 
    Bio::EnsEMBL::Analysis::Runnable::ExonerateAlignFeature->new(
     -query_seqs     => \@q_seqs,
     -query_type     => 'dna',
     -target_seqs    => \@t_seqs,
     -options        => $options,
    );

 $runnable->run; #create and fill Bio::Seq object
 my @results = $runnable->output;
 
=head1 DESCRIPTION

This module handles a specific use of the Exonerate (G. Slater) program, 
to align clone sequences with genomic sequences.

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Runnable::ExonerateAlignFeature;

use warnings ;
use strict;
use feature 'say';

use parent ('Bio::EnsEMBL::Analysis::Runnable::BaseExonerate');
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );


=head2 new

 Arg [REDUNDANT_IDS]: Hashref of arrays, the keys for the hashes are accessions
                       and the arrays contains accessions of sequences which are
                       the same as the sequence of the key
 Description        : Constructor
 Returntype         : Bio::EnsEMBL::Analysis::Runnable::ExonerateAlignFeature
 Exceptions         : None

=cut

sub new {
  my ( $class, @args ) = @_;

  my $self = $class->SUPER::new(@args);
  my ($redundant_ids) = rearrange( [qw(REDUNDANT_IDS)],@args );
  $self->redundant_ids($redundant_ids);
  return $self;
}

#
# Implementation of method in abstract superclass
#

=head2 parse_results

 Arg [1]    : FileHandle
 Description: Parse the result from the program executed
 Returntype : None
 Exceptions : Throws if failed to process the vulgar block

=cut

sub parse_results {
  my ( $self, $input, $write_results ) = @_;

  my @features;
  my $redundant_ids = $self->redundant_ids;

  my $fh;
  if($write_results) {
    unless(open($fh,$input)) {
      throw("Could not open file: ".$input);
     };
  } else {
    $fh = $input;
  }

  while (<$fh>){

    next unless /^RESULT:/;

    chomp;

    my (
      $tag, $q_id, $q_start, $q_end, $q_strand, 
      $t_id, $t_start, $t_end, $t_strand, $score, 
      $perc_id, $q_length, $t_length, $gene_orientation,
      @vulgar_blocks, $total_match_length
    ) = split;


    my $cigar_string='';  
    while (@vulgar_blocks){
      throw("Something funny has happened to the input vulgar string." .
		 "  Expecting components in multiples of three, but only have [" .
		 scalar @vulgar_blocks . "] items left to process.")
      unless scalar @vulgar_blocks >= 3;

      # check for introns if we are using exonerate2genes model
      my $str = join(',',@vulgar_blocks);
#      print "$str\n";
      if ( $str =~ /^3,(\d+),(\d+),I,(\d+),(\d+),5,(\d+),(\d+)/ or
	   $str =~ /^5,(\d+),(\d+),I,(\d+),(\d+),3,(\d+),(\d+)/ ) {
	
	splice(@vulgar_blocks,0,9);
	my $query_match_length = $2+$4+$6 - ( $1+$3+$5 );
	my $match_type = 'I';
	$cigar_string .= $query_match_length.$match_type;
#	print "INTRON $1 $2 $3 $4 $5 $6\n";
      }

      $str = join(',',@vulgar_blocks);
#      print "$str\n";
      my $match_type          = shift @vulgar_blocks;
      my $query_match_length  = shift @vulgar_blocks;
      my $target_match_length = shift @vulgar_blocks;
      $total_match_length +=  $query_match_length;

      if ($match_type eq "G"){
	if ($query_match_length == 0){
          $match_type="I";
          $query_match_length = $target_match_length;
        }elsif ($target_match_length == 0){
          $match_type="D";
	}
      }    
      $cigar_string .= $query_match_length.$match_type;     
    }
    my $hcoverage = $total_match_length / $q_length * 100;

    my $feature = 
      $self->make_feature(
        $q_id, $q_length, $q_start, $q_end, $q_strand, 
        $t_id, $t_length, $t_start, $t_end, $t_strand, 
        $score, $perc_id, $cigar_string, $hcoverage , $gene_orientation
      );

    if($feature){
      push @features, $feature;
    }else{
      warn "Clone end feature from probe :$q_id doesnt match well enough\n";
    }

    if ($redundant_ids) {
      my $redundant_ids = $redundant_ids->{$q_id};
      if(scalar(@$redundant_ids)) {
        foreach my $id (@$redundant_ids) {
          my $feature = $self->make_feature(
            $id, $q_length, $q_start, $q_end, $q_strand,
            $t_id, $t_length, $t_start, $t_end, $t_strand,
            $score, $perc_id, $cigar_string, $hcoverage , $gene_orientation);

          if($feature){
            push @features, $feature;
          } else{
            warn "Clone end feature from probe :$q_id doesnt match well enough\n";
          }
        } # end foreach my $id
      } # end if(scalar(@$redundant_ids
    }
  }

  if($write_results) {
    close $fh;
  }
  return \@features;
}

#
# Create dna align feature objects: 
#

=head2 make_feature

 Arg [1]    : String
 Arg [2]    : Integer
 Arg [3]    : Integer
 Arg [4]    : Integer
 Arg [5]    : Integer
 Arg [6]    : String
 Arg [7]    : Integer
 Arg [8]    : Integer
 Arg [9]    : Integer
 Arg [10]   : Integer
 Arg [11]   : Integer
 Arg [12]   : Integer
 Arg [13]   : String
 Arg [14]   : Integer
 Arg [15]   : Integer
 Description: Create a (Dna|Pep)AlignFeature based on the result from exonerate
 Returntype : Bio::EnsEMBL::DnaDnaAlignFeature or Bio::EnsEMBL::DnaPepAlignFeature
 Exceptions : Throws if the strand is wrong

=cut

sub make_feature{
  my ($self,
      $q_id, $q_len, $q_start, $q_end, $q_strand,
      $t_id, $t_len, $t_start, $t_end, $t_strand,
      $score, $perc_id, $cigar_string, $hcoverage, $gene_orientation) = @_;
 
  if($q_strand eq '+'){
    $q_strand = 1;
  } elsif ($q_strand eq '-') {
    $q_strand = -1;
  } else {
    $q_strand = $q_start <= $q_end ? 1 : -1;
  }
 
  if($t_strand eq '+'){
    $t_strand = 1;
  } elsif ($t_strand eq '-') {
    $t_strand = -1;
  } else {
    $t_strand = $t_start <= $t_end ? 1 : -1;
  }

  if ($t_start > $t_end) {
    ($t_start, $t_end) = ($t_end, $t_start);
  } elsif ($t_strand < 0) {
    ($t_start, $t_end) = ($t_len - $t_end, $t_len - $t_start);
  }

  # for reverse strand matches, Exonerate reports end => start 
  if ($q_start > $q_end) {
    ($q_start, $q_end) = ($q_end, $q_start);
  } elsif ($q_strand < 0) {
    # coordinates are in strand of reverse complemented sequence. Need to flip
    ($q_start, $q_end) = ($q_len - $q_end, $q_len - $q_start);
  }

  if ($gene_orientation eq '-') {
    $t_strand *= -1;
    $q_strand *= -1;
    # need to reverse the cigar string
    my @numbers = split(/\D+/,$cigar_string);
    my @letters = split(/\d+/,$cigar_string);
    $cigar_string = '';
    for (my $i = $#numbers ; $i >= 0 ; $i-- ) {
      $cigar_string .= $numbers[$i] . $letters[$i+1];
    }
  }

  # correct for exonerate half-open coords
  $q_start+=1;
  $t_start+=1;

  my $obj_name = $self->query_type =~ /dna/i 
      ? "Bio::EnsEMBL::DnaDnaAlignFeature"
      : "Bio::EnsEMBL::DnaPepAlignFeature";
  
  my $feature =
    $obj_name->new(
      -seqname      => $t_id,
      -start        => $t_start,
      -end          => $t_end,
      -strand       => $t_strand,
      -hseqname     => $q_id,
      -hstart       => $q_start,
      -hend         => $q_end,
      -hstrand      => $q_strand,
      -score        => $score,
      -percent_id   => $perc_id,
      -cigar_string => $cigar_string,
      -hcoverage    => $hcoverage,
    );

  $feature->{"_intron"} = 1 unless $gene_orientation eq '.';
  return $feature;
}


=head2 redundant_ids

 Arg [1]    : Hashref of arrays, the keys for the hashes are accessions
               and the arrays contains accessions of sequences which are
               the same as the sequence of the key
 Description: Getter/Setter
 Returntype : Hashref
 Exceptions : None

=cut

sub redundant_ids {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_redundant_ids'} = $value;
  }

  return($self->{'_redundant_ids'});
}

1;

