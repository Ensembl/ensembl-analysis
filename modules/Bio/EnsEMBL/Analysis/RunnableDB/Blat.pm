# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::Blat
#
# Copyright (c) 2009 WormBase
#

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Blat

=head1 SYNOPSIS

  my $blat = Bio::EnsEMBL::Analysis::RunnableDB::Blat->
  new(
      -input_id => 'file_name',
      -db => $db,
      -analysis => $analysis,
     );
  $blat->fetch_input;
  $blat->run;
  $blat->write_output;

=head1 DESCRIPTION

This module provides an interface between the ensembl database and
the Runnable Blat which wraps the program Blat

This module can fetch appropriate input from the database
pass it to the runnable then write the results back to the database
in the dna_align_feature  tables

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Blat;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Runnable::Blat;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB);

=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Blat
  Function  : fetch data out of database and create runnable
  Returntype: 1
  Exceptions: none
  Example   :

=cut

sub fetch_input {
  my ($self) = @_;
  my %parameters;
  if ( $self->parameters_hash ) {
    %parameters = %{ $self->parameters_hash };
  }
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Blat->new( -query     => $self->input_id,
                                                              -program   => $self->analysis->program_file,
                                                              -analysis  => $self->analysis,
                                                              -database  => $self->analysis->db_file,
                                                              -gapped    => 1,
                                                              -queryhseq => 1,
                                                              %parameters, );
  $self->runnable($runnable);
  return 1;
}

sub run {
  my ($self) = @_;

  my ( @raw_features, @filtered_features, %feats_by_hseq );

  foreach my $runnable ( @{ $self->runnable } ) {
    $runnable->run;
    push @raw_features, @{ $runnable->output };
  }

  foreach my $feat (@raw_features) {
    # remove feature where coverage is less than 25%
    next if $feat->score < 25;
    # remove features with an implausibly long span
    next if $feat->end - $feat->start + 1 > 100000;

    push @{ $feats_by_hseq{ $feat->hseqname } }, $feat;
  }

  foreach my $hseq ( keys %feats_by_hseq ) {
    my @hits = sort { $b->score <=> $a->score } @{ $feats_by_hseq{$hseq} };

    my $cutoff = $hits[0]->score*0.75;
    foreach my $hit (@hits) {
      if ( $hit->score >= $cutoff ) {
        push @filtered_features, $hit;
      }
    }
  }

  $self->output( \@filtered_features );
} ## end sub run

sub write_output {
  my ($self) = @_;

  foreach my $feat ( @{ $self->output } ) {
    my $slice = $self->db->get_SliceAdaptor->fetch_by_region( 'toplevel', $feat->seqname );
    if ( not defined $slice ) {
      throw( "Could not fetch slice from the db: " . $feat->hseqname . "\n" );
    }
    $feat->slice($slice);
  }

  $self->SUPER::write_output();
}

=head2 get_adaptor

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Blat
  Function  : get dna_align_feature adaptor
  Returntype: Bio::EnsEMBL::DBSQL::DnaAlignFeatureAdaptor
  Exceptions: none
  Example   :

=cut

sub get_adaptor {
  my ($self) = @_;
  return $self->db->get_DnaAlignFeatureAdaptor;
}

1;
