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

Bio::EnsEMBL::Analysis::Runnable::ExonerateCloneEnds - 

=head1 SYNOPSIS

  my $runnable = 
    Bio::EnsEMBL::Analysis::Runnable::ExonerateCloneEnds->new(
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

package Bio::EnsEMBL::Analysis::Runnable::ExonerateCloneEnds;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Analysis::Runnable::BaseExonerate;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::BaseExonerate);

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);

  return $self;
}

#
# Implementation of method in abstract superclass
#
sub parse_results {
  my ( $self, $fh ) = @_;
  
  my @features;
  
  while (<$fh>){

    next unless /^RESULT:/;

    chomp;
    
    my (
      $tag, $q_id, $q_start, $q_end, $q_strand, 
      $t_id, $t_start, $t_end, $t_strand, $score, 
      $perc_id, $q_length, $t_length, $gene_orientation,
      @vulgar_blocks
    ) = split;

    my $cigar_string='';  
    while (@vulgar_blocks){
      throw("Something funny has happened to the input vulgar string." .
		 "  Expecting components in multiples of three, but only have [" .
		 scalar @vulgar_blocks . "] items left to process.")
      unless scalar @vulgar_blocks >= 3;

      my $match_type          = shift @vulgar_blocks;
      my $query_match_length  = shift @vulgar_blocks;
      my $target_match_length = shift @vulgar_blocks;


      if ($match_type eq "G"){
	if ($query_match_length == 0){
	    $match_type="D";
            $query_match_length = $target_match_length;
        }elsif ($target_match_length == 0){
            $match_type="I";
	}

      }
    
      $cigar_string .= $query_match_length.$match_type;
     
    }

    my $feature = 
      $self->make_feature(
        $q_id, $q_start, $q_end, $q_strand, 
        $t_id, $t_start, $t_end, $t_strand, $score, 
        $perc_id, $q_length, $cigar_string
      );

    if($feature){
      push @features, $feature;
    }else{
      warn "Clone end feature from probe :$q_id doesnt match well enough\n";
    }
  }

  return \@features;
}

#
# Create dna align feature objects: 
#
sub make_feature{
  my ($self, @args) = @_;
  
  my (
    $tag, $q_id, $q_start, $q_end, $q_strand, 
    $t_id, $t_start, $t_end, $t_strand, $score, 
    $perc_id, $q_length, $cigar_string 
  ) = @_;
 
  if($q_strand eq '+'){
    $q_strand = 1;
    if($t_strand eq '+'){
      $t_strand = 1;
    }elsif($t_strand eq '-'){
      $t_strand = -1;
    }else{
      throw "unrecognised target strand symbol: $t_strand\n";
    }
  }elsif($q_strand eq '-'){
    $q_strand = -1;
    if($t_strand eq '-'){
      $t_strand = -1;
    }elsif($t_strand eq '+'){
      $t_strand = 1;
    }else{
      throw "unrecognised target strand symbol: $t_strand\n";
    }
  }else{    
      throw "unrecognised query strand symbol: $q_strand\n";
  }
  
  # Exonerate reports query start -1 so in some cases we get alignments with hit start 0
  # we add 1 to avoid this situation.
  $q_start+=1;
  $t_start+=1;

  # for reverse strand matches, Exonerate reports end => start 
  if ($q_start > $q_end) {
    ($q_start, $q_end) = ($q_end, $q_start);
  }
  if ($t_start > $t_end) {
    ($t_start, $t_end) = ($t_end, $t_start);
  }

  my $feature =
    new Bio::EnsEMBL::DnaDnaAlignFeature(
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
    );
  
  return $feature;
}

1;

