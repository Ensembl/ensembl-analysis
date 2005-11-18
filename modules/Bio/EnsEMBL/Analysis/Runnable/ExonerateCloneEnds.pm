
=pod

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript

=head1 SYNOPSIS

  my $runnable = 
    Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript->new(
     -query_seqs     => \@q_seqs,
     -query_type     => 'dna',
     -target_seqs    => \@t_seqs,
     -options        => $options,
    );

 $runnable->run; #create and fill Bio::Seq object
 my @results = $runnable->output;
 
=head1 DESCRIPTION

This module handles a specific use of the Exonerate (G. Slater) program, to

=head1 CONTACT

ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::Runnable::ExonerateCloneEnds;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Analysis::Runnable::BaseExonerate;
use Bio::EnsEMBL::Feature;
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
    #print STDERR $_ if $self->_verbose;

    next unless /^RESULT:/;

    chomp;
    
    my (
      $tag, $q_id, $q_start, $q_end, $q_strand, 
      $t_id, $t_start, $t_end, $t_strand, $score, 
      $perc_id, $q_length, $t_length, $gene_orientation,
 #     $cigar_string, $analysis,
      @vulgar_blocks
    ) = split;
    
    my($match_type, $query_match_length, $target_match_length,@rest_of_vulgar) = @vulgar_blocks;
    
    if(@rest_of_vulgar){
      throw (
        "There is more than a simple match ('M') vulgar output: @vulgar_blocks for tag $q_id mapped to region $t_id \n"
      );
    }
    
    if(!($match_type eq 'M')){
      throw "I have received a starting Vulgar symbol $match_type which is not a match!"; 
    }
    
    my $feature = 
      $self->make_feature(
        $q_id, $q_start, $q_end, $q_strand, 
        $t_id, $t_start, $t_end, $t_strand, $score, 
        $perc_id, $q_length, $query_match_length,
        #, $cigar_string, $percent_id, $analysis
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
    $perc_id, $q_length, $query_match_length
    #, $cigar_string, $percent_id, $analysis
  ) = @_;

  #because of the 'in-between' coordinates.
  $t_start += 1;
  
  my $mismatch_count;
  
  # If we miss by more than one, don't create the feature.
  if($query_match_length == $q_length){
    if($score == 125){
      $mismatch_count = 0;
    }else{
      $mismatch_count = 1;
    }
  }elsif($query_match_length == $q_length -1){
    $mismatch_count = 1;
  }else{
    return undef;
  }

  # Everything is 'flipped' into the forward strand of the probe -
  # so a hit on the reverse strand of the probe (q_strand = '-1')
  # is altered:  q_strand = '+1' and t_strand => -1 x t_strand. 
  if($q_strand eq '+'){
    if($t_strand eq '+'){
      $t_strand = 1;
    }elsif($t_strand eq '-'){
      $t_strand = -1;
    }else{
      throw "unrecognised target strand symbol: $t_strand\n";
    }
  }elsif($q_strand eq '-'){
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
    new Bio::EnsEMBL::DnaAlignFeature(
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
  #    -analysis     => $analysis
    );

  # attach the slice name onto the feature: let the runnabledb
  # sort out whether it's valid.
  $feature->seqname($t_id);
  
  return $feature;
}

1;

