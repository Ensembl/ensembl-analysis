
=pod

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::ExonerateCloneEnds

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

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::Runnable::ExonerateCloneEnds;

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
#  print $cigar_string,"\n";

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

  #because of the 'in-between' coordinates.
  $t_start += 1;
  

  # Everything is 'flipped' into the forward strand of the probe -
  # so a hit on the reverse strand of the probe (q_strand = '-1')
  # is altered:  q_strand = '+1' and t_strand => -1 x t_strand. 
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
      $t_strand = 1;
    }elsif($t_strand eq '+'){
      $t_strand = -1;
    }else{
      throw "unrecognised target strand symbol: $t_strand\n";
    }
  }else{    
      throw "unrecognised query strand symbol: $q_strand\n";
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
#  print $t_id," Percent_id: ",$perc_id," Score: ", $score,"\n";
  
  return $feature;
}

1;

