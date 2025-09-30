# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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
# BioPerl module for PairAlign object
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

PairAlign - Dna pairwise alignment module

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Contains list of sub alignments making up a dna-dna alignment

Creation:
   
    my $pair = new Bio::EnsEMBL::FeaturePair(-start  => $qstart,
                                             -end    => $qend,
                                             -strand => $qstrand, 
                                             -hstart => $hstart,
                                             -hend   => $hend, 
                                              -hend   => $hstrand,        
    );

    my $pairaln   = new Bio::EnsEMBL::Analysis::Tools::PairAlign;
    $pairaln->addFeaturePair($pair);

Any number of pair alignments can be added to the PairAlign object


Manipulation:

To convert between coordinates : 

    my $cdna_coord = $pair->genomic2cDNA($gen_coord);
    my $gen_coord  = $pair->cDNA2genomic($cdna_coord);

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::Analysis::Tools::PairAlign;

use warnings ;
use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
use vars qw(@ISA);
use strict;


@ISA = qw();

sub new {
  my($class,@args) = @_;
  my $self = {};
  bless $self, $class;

  return $self; # success - we hope!
}

sub addFeaturePair {
  my ($self,$pair) = @_;
  push(@{$self->{'_pairs'}},$pair);
}


=head2 eachFeaturePair

 Title   : eachFeaturePait
 Usage   : my @pairs = $pair->eachFeaturePair
 Function: 
 Example : 
 Returns : Array of Bio::SeqFeature::FeaturePair
 Args    : none


=cut

sub eachFeaturePair {
  my ($self) = @_;

  if (defined($self->{'_pairs'})) {
    return $self->{'_pairs'};
  }
}


=head2 get_hstrand

  Arg [1]   : Bio::EnsEMBL::Analysis::Tools::PairAlign
  Function  : gets the hit strand of the first feature in the array
  Returntype: 1/-1
  Exceptions: 
  Example   : 

=cut



sub get_hstrand {
  my ($self) = @_;

  my $features = $self->eachFeaturePair;

  return $features->[0]->hstrand;
}

=head2 genomic2cDNA

 Title   : genomic2cDNA
 Usage   : my $cdna_coord = $pair->genomic2cDNA($gen_coord)
 Function: Converts a genomic coordinate to a cdna coordinate
 Example : 
 Returns : int
 Args    : int


=cut

sub genomic2cDNA {
  my ($self,$coord) = @_;
  my @pairs = @{$self->eachFeaturePair};

  @pairs = sort {$a->start <=> $b->start} @pairs;

  my $newcoord;

 HOMOL: while (my $sf1 = shift(@pairs)) {
    next HOMOL unless ($coord >= $sf1->start && $coord <= $sf1->end);

    if ($sf1->strand == 1 && $sf1->hstrand == 1) {
      $newcoord = $sf1->hstart + ($coord - $sf1->start);
      last HOMOL;
    } elsif ($sf1->strand == 1 && $sf1->hstrand == -1) {
      $newcoord = $sf1->hend   - ($coord - $sf1->start);
      last HOMOL;
    } elsif ($sf1->strand == -1 && $sf1->hstrand == 1) {
      $newcoord = $sf1->hstart + ($sf1->end - $coord);
      last HOMOL;
    } elsif ($sf1->strand == -1 && $sf1->hstrand == -1) {
      $newcoord = $sf1->hend   - ($sf1->end - $coord);
      last HOMOL;
    } else {
      throw("ERROR: Wrong strand value in FeaturePair (" . $sf1->strand . 
            "/" . $sf1->hstrand . "\n");
    }
  }

  if (defined($newcoord)) {
    return $newcoord;
  } else {
    throw("Couldn't convert $coord");
  }
}

=head2 cDNA2genomic

 Title   : cDNA2genomic
 Usage   : my $gen_coord = $pair->genomic2cDNA($cdna_coord)
 Function: Converts a cdna coordinate to a genomic coordinate
 Example : 
 Returns : int
 Args    : int


=cut

sub cDNA2genomic {
  my ($self,$coord) = @_;

  my @pairs = @{$self->eachFeaturePair};

  my $newcoord;

 HOMOL: while (my $sf1 = shift(@pairs)) {
    next HOMOL unless ($coord >= $sf1->hstart && $coord <= $sf1->hend);

    if ($sf1->strand == 1 && $sf1->hstrand == 1) {
      $newcoord = $sf1->start + ($coord - $sf1->hstart);
      last HOMOL;
    } elsif ($sf1->strand == 1 && $sf1->hstrand == -1) {
      $newcoord = $sf1->start  +($sf1->hend - $coord);
      last HOMOL;
    } elsif ($sf1->strand == -1 && $sf1->hstrand == 1) {
      $newcoord = $sf1->end   - ($coord - $sf1->hstart);
      last HOMOL;
    } elsif ($sf1->strand == -1 && $sf1->hstrand == -1) {
      $newcoord = $sf1->end   - ($sf1->hend - $coord);
      last HOMOL; 
    } else {
      throw("ERROR: Wrong strand value in homol (" . $sf1->strand . "/" . 
            $sf1->hstrand . "\n");
    }
  }

  if (defined ($newcoord)) {
    return $newcoord;
  } else {
    throw("Couldn't convert $coord\n");
  }
}

sub find_Pair {
  my ($self,$coord) = @_;

  foreach my $p (@{$self->eachFeaturePair}) {
    if ($coord >= $p->hstart && $coord <= $p->hend) {
      return $p;
    }
  }
}

=head2 convert_cDNA_feature

 Title   : convert_cDNA_feature
 Usage   : my @newfeatures = $self->convert_cDNA_feature($f);
 Function: Converts a feature on the cDNA into an array of 
           features on the genomic (for features that span across introns);
 Example :
 Returns : @Bio::EnsEMBL::FeaturePair
 Args    : Bio::EnsEMBL::FeaturePair

=cut

sub convert_cDNA_feature {
  my ($self,$feature) = @_;

  my $foundstart = 0;
  my $foundend   = 0;

  my @pairs = @{$self->eachFeaturePair};
  my @newfeatures;

 HOMOL:while (my $sf1 = shift(@pairs)) {
    my $skip = 0;

    $skip = 1 unless ($feature->start >= $sf1->hstart 
                      && $feature->start <= $sf1->hend);

    if($skip){
      next HOMOL;
    }
    if ($feature->end >= $sf1->hstart && $feature->end <= $sf1->hend) {
      $foundend = 1;
    }

    my $startcoord = $self->cDNA2genomic($feature->start);
    my $endcoord;

    if ($sf1->hstrand == 1) {
      $endcoord   = $sf1->end;
    } else {
      $endcoord   = $sf1->start;
    }

    if ($foundend) {
      $endcoord = $self->cDNA2genomic($feature->end);
    }

    my $tmpf = new Bio::EnsEMBL::Feature(-seqname => $feature->seqname,
                                         -start   => $startcoord,
                                         -end     => $endcoord,
                                         -strand  => $feature->strand);
    push(@newfeatures,$tmpf);
    last;
  }

  # Now the rest of the pairs until we find the endcoord

  while ((my $sf1 = shift(@pairs)) && ($foundend == 0)) {

    if ($feature->end >= $sf1->hstart && $feature->end <= $sf1->hend) {
      $foundend = 1;
    }

    my $startcoord;
    my $endcoord;

    if ($sf1->hstrand == 1) {
      $startcoord = $sf1->start;
      $endcoord   = $sf1->end;
    } else {
      $startcoord = $sf1->end;
      $endcoord   = $sf1->start;
    }

    if ($foundend) {
      $endcoord = $self->cDNA2genomic($feature->end);
    }

    my $tmpf = new Bio::EnsEMBL::Feature(-seqname => $feature->seqname,
                                         -start   => $startcoord,
                                         -end     => $endcoord,
                                         -strand  => $feature->strand);
    push(@newfeatures,$tmpf);
  }

  return \@newfeatures;
}


sub convert_FeaturePair {
  my ($self,$pair) = @_;

  my $hstrand = $self->get_hstrand;

  my $feat = $self->create_Feature($pair->start, $pair->end, $pair->strand,
                                   $pair->slice);
  my @newfeatures = @{$self->convert_cDNA_feature($feat)};
  my @newpairs;

  my $hitpairaln  = new Bio::EnsEMBL::Analysis::Tools::PairAlign;
  $hitpairaln->addFeaturePair($pair);

  foreach my $new (@newfeatures) {

    # Now we want to convert these cDNA coords into hit coords

    my $hstart1 = $self->genomic2cDNA($new->start);
    my $hend1   = $self->genomic2cDNA($new->end);

    my $hstart2 = $hitpairaln->genomic2cDNA($hstart1);
    my $hend2   = $hitpairaln->genomic2cDNA($hend1);

    # We can now put the final feature together

    my $finalstrand = $hstrand * $pair->strand * $pair->hstrand;

    if ($hstart2 > $hend2) {
      my $tmp = $hstart2;
      $hstart2 = $hend2;
      $hend2   = $tmp;
    }

    my $finalpair = $self->create_FeaturePair($new->start, $new->end, 
                                              $new->strand,
                                              $hstart2, $hend2, 
                                              $finalstrand, $pair->score);

    push(@newpairs,$finalpair);

  }

  return \@newpairs;
}

sub create_FeaturePair {
    my ($self, $start, $end, $strand, $hstart, $hend, 
        $hstrand, $score) = @_;
   
    my $fp = Bio::EnsEMBL::FeaturePair->new(
                                            -start    => $start,
                                            -end      => $end,
                                            -strand   => $strand,
                                            -hstart   => $hstart,
                                            -hend     => $hend,
                                            -hstrand  => $hstrand,
                                            -score    => $score,
                                           );

   
    return $fp;
}

sub create_Feature{
  my ($self, $start, $end, $strand,  $slice) = @_;

  my $feat = new Bio::EnsEMBL::Feature(-start   => $start,
                                       -end     => $end,
                                       -strand  => $strand,
                                       -slice   => $slice,
                                      );
  return $feat;
}

1;







