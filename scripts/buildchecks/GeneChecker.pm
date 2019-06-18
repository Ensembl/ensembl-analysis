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
# limitations under the License.#
# EnsEMBL module for GeneChecker
#
# Cared for by Steve Searle <searle@sanger.ac.uk>
#

# POD documentation - main docs before the code

=pod

=head1 NAME

Bio::EnsEMBL::Test::GeneChecker


=head1 SYNOPSIS
Module to check the validity of a gene

=head1 DESCRIPTION
Performs various checks on a Bio::EnsEMBL::Gene object. These include:
  1. Genes with long extents
  2. Genes with no transcripts
  2. Genes with more than a certain number of transccripts
  3. Genes with transcripts on both strands

This class does not use stable_ids but instead dbIDs because stable_ids are
not set until after the gene build.

=head1 CONTACT

  Steve Searle <searle@sanger.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package GeneChecker;
use vars qw(@ISA);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use strict;

use Checker;
use Bio::EnsEMBL::Intron;
use Bio::Seq;


@ISA = qw( Checker );


sub new {
  my ($class, @args) = @_;

  my $self = bless {},$class;

  my ( $gene, $maxgenelen, $maxtransgene, $ignorewarnings, $slice, $adaptor ) = rearrange
	 ( [ qw ( GENE
		  MAXGENELEN
		  MAXTRANSGENE
                  IGNOREWARNINGS
                  SLICE
                  ADAPTOR
	      )], @args );


  if( !defined $gene ) {
    $self->throw("Gene must be set in new for GeneChecker")
  }
  $self->gene($gene);

  if( defined $maxgenelen ) {
    $self->maxgenelen( $maxgenelen );
  } else {
    $self->maxgenelen(2_000_000);
  }
  if( defined $maxtransgene ) {
    $self->maxtransgene( $maxtransgene );
  } else {
    $self->maxtransgene(10);
  }

  if( defined $ignorewarnings) { $self->ignorewarnings($ignorewarnings); }
  if( defined $slice) { $self->slice($slice); }
  if( defined $adaptor ) { $self->adaptor( $adaptor )}

  $self->{_errors} = [];
  $self->{_warnings} = [];

  return $self;
}


sub maxtransgene {
  my ( $self, $arg ) = @_;
  if (defined $arg) {
    $self->{_maxtransgene} = $arg;
  }
  return $self->{_maxtransgene};
}

sub maxgenelen {
  my ( $self, $arg ) = @_;
  if (defined $arg) {
    $self->{_maxgenelen} = $arg;
  }
  return $self->{_maxgenelen};
}


=head2 gene

 Title   : gene
 Usage   : $obj->gene($newval)
 Function:
 Returns : value of gene
 Args    : newvalue (optional)

=cut

sub gene {
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      if (!($value->isa('Bio::EnsEMBL::Gene'))) {
         $self->throw("gene passed a non Bio::EnsEMBL::Gene object\n");
      }
      $self->{_gene} = $value;
    }
    return $self->{_gene};

}

sub slice {
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      if (!($value->isa('Bio::EnsEMBL::Slice'))) {
         $self->throw("vc passed a non Bio::EnsEMBL::Slice object\n");
      }
      $self->{_slice} = $value;
    }
    return $self->{_slice};

}

sub adaptor {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{_adaptor} = $arg;
  }
  return $self->{_adaptor};
}


sub output {
  my $self = shift;

  my $gene = $self->gene;

  print "\n++++++++++++++++++\n";

  print "Gene " . $self->gene->dbID . " " .
        (defined($self->gene->stable_id) ? $self->gene->stable_id : "") . " type " . $self->gene->biotype . "\n";


  my $vcoffset = 0;
  if ($self->slice) {
    $vcoffset = $self->slice->start - 1;
  }

  print " Extents:         " . ($self->gene->start+$vcoffset) . " " . ($self->gene->end+$vcoffset) . " " . $self->gene->strand . "\n";
  print " Length:          " . $self->gene->length . "\n";
  print " Num Transcripts: " . scalar(@{$self->gene->get_all_Transcripts}) . "\n";
  foreach my $trans (@{$self->gene->get_all_Transcripts}) {
    print "   Transcript " . $trans->dbID . " " . ($trans->start+$vcoffset) . " " . ($trans->end+$vcoffset) . " " . $trans->length . "\n";
  }

  $self->SUPER::output;
}


sub check {
  my $self = shift;

  $self->check_Structure;

}

sub check_Structure {
  my $self = shift;

  if ($self->gene->length > $self->maxgenelen) {
    $self->add_Error("Gene too long\n", 'genelen');
  }

  if (scalar(@{$self->gene->get_all_Transcripts}) == 0) {
    $self->add_Error("No transcripts in gene\n", 'notransgene');
  }

  if (scalar(@{$self->gene->get_all_Transcripts}) > $self->maxtransgene) {
    $self->add_Error("Too many transcripts in gene\n", 'numtransgene');
  }

  my $firststrand = undef;
  foreach my $trans (@{$self->gene->get_all_Transcripts}) {
    if (!defined($firststrand)) {
      $firststrand = $trans->strand;
    }
    if ($trans->strand != $firststrand) {
      $self->add_Error("Transcripts on both strands in gene\n", 'bothstrandgene');
    }
  }
}

1;
