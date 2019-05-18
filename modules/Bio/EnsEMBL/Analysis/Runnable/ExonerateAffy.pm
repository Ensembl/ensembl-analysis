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

Bio::EnsEMBL::Analysis::Runnable::ExonerateAffy - 

=head1 SYNOPSIS

  my $runnable = 
    Bio::EnsEMBL::Analysis::Runnable::ExonerateAffy->new(
     -query_seqs     => \@q_seqs,
     -query_type     => 'dna',
     -target_seqs    => \@t_seqs,
     -options        => $options,
    );

 $runnable->run; #create and fill Bio::Seq object
 my @results = $runnable->output;
 
=head1 DESCRIPTION

This module handles a specific use of the Exonerate (G. Slater) program, to
align Affymetrix probes to a target genome. (The resulting alignments will
probably be stored in an ensembl db as Bio::EnsEMBL::AffyFeature objects.)

NOTE: the AffyFeature objects refer to AffyProbe id's, and they in turn 
refer to AffyArray id's. Both the AffyProbes and the AffyArray's should
be pre-loaded into the ensembl db: there are separate RunnableDB
/RunnableDB's to do this from the Affymetrix data sets. This runnable
just creates fake affy probes in order to create reasonable-looking
affy features.

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Runnable::ExonerateAffy;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Analysis::Runnable::BaseExonerate;
use Bio::EnsEMBL::AffyArray;
use Bio::EnsEMBL::AffyProbe;
use Bio::EnsEMBL::AffyFeature;
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
  
  my @affy_features;
  
  while (<$fh>){
    #print STDERR $_ if $self->_verbose;

    next unless /^RESULT:/;

    chomp;
    
    my (
      $tag, $q_id, $q_start, $q_end, $q_strand, 
      $t_id, $t_start, $t_end, $t_strand, $score, 
      $perc_id, $q_length, $t_length, $gene_orientation,
      @vulgar_blocks
    ) = split;
    
    my($match_type, $query_match_length, $target_match_length,@rest_of_vulgar) = @vulgar_blocks;
    
    if(@rest_of_vulgar){
      throw (
        "There is more than a simple match ('M') vulgar output: @vulgar_blocks for 
        affy tag $q_id mapped to region $t_id \n"
      );
    }
    
    if(!($match_type eq 'M')){
      throw "I have received a starting Vulgar symbol $match_type which is not a match!"; 
    }
    
    my $affy_feature = 
      $self->make_affy_feature(
        $q_id, $q_start, $q_end, $q_strand, 
        $t_id, $t_start, $t_end, $t_strand, $score, 
        $q_length, $query_match_length
      );

    if($affy_feature){
      push @affy_features, $affy_feature;
    }else{
      warn "Affy feature from probe :$q_id doesnt match well enough\n";
    }
  }

  return \@affy_features;
}

#
# Create affy features attached to 'fake' affy array's and affy probe
# objects: these are identified by name only, and are then persisted by
# runnabledb: the runnabledb can sort out whether the attached probe/array
# is valid (exists in the db).
#
sub make_affy_feature{
  my ($self, @args) = @_;
  
  my (
    $tag, $q_id, $q_start, $q_end, $q_strand, 
    $t_id, $t_start, $t_end, $t_strand, $score, 
    $q_length, $query_match_length
  ) = @_;

  #because of the 'in-between' coordinates.
  $t_start += 1;
  
  my $probe_internal_id = $q_id;
  my $mismatch_count;
  
  if(!($probe_internal_id =~ /\d+/)){
    throw "Affy probe headers MUST be the internal db ids of the affy probes for this parser to work!\n";
  }
  
  # If we miss by more than one, don't create the feature.
  if($query_match_length == $q_length){
    if($score == 125){
      $mismatch_count = 0;
    }elsif ($score == 116) {
      $mismatch_count = 1;
    } else {
      return undef;
    }
  }elsif($query_match_length == $q_length -1){
    if ($score == 120) {
      $mismatch_count = 1;
    } else {
      return undef;
    }
  }else{
    return undef;
  }

  #Create a dummy probe object to pass out. It must be checked and
  #replaced by the parent runnabledb before being stored.
  my $affy_probe = 
    new Bio::EnsEMBL::AffyProbe(
      -dbID =>$q_id,
      -name => "replace this probe name",
      -arrayname => "replace this array name",
    );
   
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
  
  my $affy_feature =
    new Bio::EnsEMBL::AffyFeature(
      -probe => $affy_probe,
      -start => $t_start,
      -end => $t_end,
      -strand => $t_strand,
      -mismatchcount => $mismatch_count
    );

  # attach the slice name onto the feature: let the runnabledb
  # sort out whether it's valid.
  $affy_feature->seqname($t_id);
  
  return $affy_feature;
}

1;

