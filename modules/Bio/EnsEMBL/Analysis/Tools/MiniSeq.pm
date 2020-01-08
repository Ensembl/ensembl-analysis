# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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
# Spangle version of the ensembl Transcript object

# POD documentation - main docs before the code

=head1 NAME

MiniSeq - an artificially constructed cDNA on a genomic sequence

=head1 SYNOPSIS

This module is used when we only want to run an analysis over
a part or multiple parts of a sequence. 

=head1 DESCRIPTION

Contains details of coordinates of all exons that make
up a gene transcript.

Creation:
   
     my $mini = new Bio::EnsEMBL::Analysis::Tools::MiniSeq(-id      => $id,
						    -pairaln => $pairaln);



   $pairaln is a Bio::EnsEMBL::Analysis::PairAlign object containing
   one or more feature pairs that represent the mapping of
   the genomic coords to the cDNA coords

Manipulation:

    my $align   = $mini->get_PairAlign;
    my $cdnaseq = $mini->get_cDNA_sequence;

# We now do some analysis on the cdnaseq that puts sequence
# features on it.  We don't want the coordsin the cDNA frame so
# we now pass them back to the miniseq to convert them into 
# genomic coordinates.

    my  @newfeatures = $mini->convert_FeaturePairs(@featurepairs);


=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::Analysis::Tools::MiniSeq;
use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis::Tools::PairAlign;
use Bio::PrimarySeq;
use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump verbose);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

@ISA = qw();

sub new {
  my($class,@args) = @_;
  my $self = {};
  bless $self, $class;

  my ($id,$pairalign) = rearrange([qw(ID PAIRALIGN)],@args);
  verbose('warning');
  #No defaults

  $self->id($id);
  $self->pairAlign($pairalign);

  throw("No input id for MiniSeq")        unless($self->id);
  throw("No input pairalign for MiniSeq") unless($self->pairAlign);

  
  return $self;
}

=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: 
 Returns : value of id
 Args    : newvalue (optional)


=cut

sub id {
  my ($self,$arg) = @_;

  if(defined($arg)) {
    $self->{'_id'} = $arg;
  }

  return $self->{'_id'};
}

=head2 pairAlign

 Title   : pairAlign
 Usage   : $self->pairAlign($pair)
 Function: Get/set method for the pairalign object
           that stores the cDNA-genomic exon mapping
 Returns : Bio::EnsEMBL::Analysis::PairAlign
 Args    : 


=cut

sub pairAlign {
  my ($self,$pair) = @_;

  if ($pair) {
    if( ! $pair->isa("Bio::EnsEMBL::Analysis::Tools::PairAlign") ) {
      throw("$pair is not a Bio::EnsEMBL::Analysis::Tools::PairAlign!");
    }
    foreach my $p (@{$pair->eachFeaturePair}) {
      if ($p->strand != 1) {
        throw("Can't have a PairAlign object where the strand of the first ".
              "object is reversed");
      }
    }
    $self->{'_pair'} = $pair;
  }
  return $self->{'_pair'};
}


=head2 get_cDNA_sequence

 Title   : get_cDNA_sequence
 Usage   : my $seq = $self->get_cDNA_sequence
 Function: Returns the cdna sequence corresponding
           to the cDNA in the pairAlign object
 Example : 
 Returns : Bio::PrimarySeq
 Args    : none


=cut

sub get_cDNA_sequence {
  my ($self) = @_;

  my $seqstr = "";

  my @exons = @{$self->pairAlign->eachFeaturePair};
  return unless (scalar @exons > 0);

  foreach my $exon (@exons) {
    $seqstr .= $exon->seq;
  }
  return new Bio::PrimarySeq('-id' => "genomic" ,
                             -seq => $seqstr);

}

=head2 convert_FeaturePair

 Title   : convert_FeaturePair
 Usage   : my @newfeatures = $self->convert_FeaturePairs($feature)
 Function: Converts feature pair coordinates on the cDNA sequence
           into an array of feature pairs on the genomic sequence
 Example : 
 Returns : Bio::EnsEMBL::FeaturePair
 Args    : Array of Bio::EnsEMBL::FeaturePair


=cut

sub convert_FeaturePair {
    my ($self,$feature) = @_;

    my @newfeatures;

    my @tmp = @{$self->pairAlign->convert_FeaturePair($feature)};
    push(@newfeatures,@tmp);

    return \@newfeatures;
}

=head2 convert_SeqFeature

 Title   : convert_FeaturePair
 Usage   : my @newfeatures = $self->convert_FeaturePairs($feature)
 Function: Converts feature coordinates on the cDNA sequence
           into an array of features on the genomic sequence
 Example : 
 Returns : Bio::EnsEMBL::FeaturePair
 Args    : Array of Bio::EnsEMBL::FeaturePair


=cut

sub convert_SeqFeature {
    my ($self,$feature) = @_;

    my @newfeatures;

    my @tmp = @{$self->pairAlign->convert_cDNA_feature($feature)};
    push(@newfeatures,@tmp);

    return \@newfeatures;
}

=head2 convert_PepFeaturePair

 Title   : convert_PepFeaturePair
 Usage   : my @newfeatures = $self->convert_PepFeaturePair($feature)
 Function: Converts feature pair coordinates on the cDNA sequence
           into an array of feature pairs on the genomic sequence
           Peptide coordinates are maintained.
 Example : 
 Returns : Array of Bio::EnsEMBL::FeaturePair
 Args    : Bio::EnsEMBL::FeaturePair


=cut

sub convert_PepFeaturePair {
    my ($self,$feature) = @_;

    my @newfeatures;

    my @tmp = @{$self->pairAlign->convert_FeaturePair($feature)};

    # replace protein coordinates
    $tmp[0]->hstart($feature->hstart);
    $tmp[0]->hend($feature->hend);
# SMJS strand of peptide should be positive?
    $tmp[0]->hstrand(1);

    push(@newfeatures,@tmp);

    return \@newfeatures;
}

1;
