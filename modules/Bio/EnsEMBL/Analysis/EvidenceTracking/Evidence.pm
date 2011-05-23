=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Storable;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
#use Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::EvidenceAdaptor;

@ISA = qw(Bio::EnsEMBL::Storable);


sub new {
  my($class,@args) = @_;

  my $self = bless {},$class;

  my ($id, $hit_name, $is_aligned, $seq_region_id, $seq_region_start,$seq_region_end, $seq_region_strand, $adaptor) =
          rearrange([qw(DBID
                        HIT_NAME       
                        IS_ALIGNED
                        SEQ_REGION_NAME
                        SEQ_REGION_START
                        SEQ_REGION_END
                        SEQ_REGION_STRAND
                        ADAPTOR
                        )],@args);

  $self->dbID( $id );
  $self->is_aligned( $is_aligned || 0 );
#  $self->input_seq_id( $inputseq_id );
  $self->hit_name( $hit_name );
  $self->seq_region_name( $seq_region_id );
  $self->seq_region_start( $seq_region_start );
  $self->seq_region_end( $seq_region_end );
  $self->seq_region_strand( $seq_region_strand );
  $self->adaptor( $adaptor );
  return $self; # success - we hope!
}


sub add_track {
  my $self = shift;
  my $track = shift;

  if ($track and !$track->isa('Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack')) {
      throw('add_track is waiting for a Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack object!');
  }
  if ($track->is_last) {
      unshift(@{$self->{'_tracks'}}, $track);
  }
  else {
      push(@{$self->{'_tracks'}}, $track);
  }
}

#sub update_current_track {
#    my $self = shift;
#    my ($key, $value) = @_;
#
#    $self->get_track->$key($value);
#}

sub get_tracks {
    my $self = shift;
    return $self->{'_tracks'};
}

sub has_tracks {
    my $self = shift;
    return exists $self->{'_tracks'};
}


sub fetch_track_by_logic_name {
    my $self = shift;
    my $logic_name = shift;

    my @a_tracks;
    foreach my $track (@{$self->{'_tracks'}}) {
        push(@a_tracks, $track) if ($track->logic_name eq $logic_name);
    }
    return \@a_tracks;
}


sub input_seq_id {
  my $self = shift;
  $self->{'inputseq_id'} = shift if ( @_ );
  return $self->{'inputseq_id'};
}

sub hit_name {
  my $self = shift;
  $self->{'hit_name'} = shift if ( @_ );
  return $self->{'hit_name'};
}

sub is_aligned {
  my $self = shift;
  $self->{'is_aligned'} = shift if ( @_ );
  return $self->{'is_aligned'};
}

sub seq_region_name {
  my $self = shift;
  $self->{'seq_region_id'} = shift if ( @_ );
  return $self->{'seq_region_id'};
}

sub seq_region_start {
  my $self = shift;
  $self->{'seq_region_start'} = shift if ( @_ );
  return $self->{'seq_region_start'};
}

sub seq_region_end {
  my $self = shift;
  $self->{'seq_region_end'} = shift if ( @_ );
  return $self->{'seq_region_end'};
}

sub seq_region_strand {
  my $self = shift;
  $self->{'seq_region_strand'} = shift if ( @_ );
  return $self->{'seq_region_strand'};
}

sub is_stored {
  my $self = shift;
  my $db = shift;

  # uniquely defined by the evidence_id
  # and the location on the genome
  my $stored_objs = $db->get_EvidenceAdaptor->fetch_all_by_hit_name($self->hit_name);

  foreach my $stored (@$stored_objs) {
#    print STDERR "seq_region_id:".$self->seq_region_id.", seq_region_start:".$self->seq_region_start.
#                 " seq_region_end:".$self->seq_region_end." seq_region_strand:".$self->seq_region_end."\n".
#                 "seq_region_id:".$stored->seq_region_id.", seq_region_start:".$stored->seq_region_start.
#                 " seq_region_end:".$stored->seq_region_end." seq_region_strand:".$stored->seq_region_end."\n";
    if ($self->is_aligned eq $stored->is_aligned) {
        if ($self->is_aligned eq 'y') {
            next unless ($self->seq_region_name eq $stored->seq_region_name &&
                        $self->seq_region_start == $stored->seq_region_start &&
                        $self->seq_region_end == $stored->seq_region_end &&
                        $self->seq_region_strand == $stored->seq_region_strand );
        }
#      warning("Evidence is_stored: Looking for exact match coordinates. Would you rather look for overlaps?");
        $self->dbID($stored->dbID);
      return $stored->dbID;
    }
  }

  return 0;
}

1;


