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

Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence - Store Evidence informations

=head1 SYNOPSIS

use Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence;

my $evidence = Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence->new(
    -hit_name => $name,
    -is_aligned => 'n'
    );

my $evidence2 = Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence->new(
    -hit_name => $other_name,
    -is_aligned => 'y',
    -seq_region_id => $seq_region_id,
    -seq_region_start => $seq_region_start,
    -seq_region_end => $seq_region_end,
    -seq_region_strand => $seq_region_strand
    );

=head1 DESCRIPTION

  Object to store the information about the evidence used. At the moment an evidence
  is unique by its name and position. If the same sequence has 2 alignements, one
  overlaping the other, there will be 2 evidences. They need to have a perfect match
  to be considered the same.
  If a sequence was used in an analysis but did not align it will have the flag is_aligned
  set to 'n'.
  If it was aligned, is_aligned is set to 'y', and the position is provided.
  At the beginning of the analysis all is_aligned flags are set to 'u'.

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Storable;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

@ISA = qw(Bio::EnsEMBL::Storable);

=head2 new

 Arg [1]    : $dbid, int
 Arg [2]    : $hit_name, string, name of the evidence
 Arg [3]    : $is_aligned, enum, 'y', 'n' or 'u'
 Arg [4]    : $seq_region_id, int
 Arg [5]    : $seq_region_start, int
 Arg [6]    : $seq_region_end, int
 Arg [7]    : $seq_region_strand, int
 Arg [8]    : $adaptor, Bio::EnsEMBL::Analysis::DBSQL::EvidenceAdaptor object
 Example    : $evidence = Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence->new(
    -hit_name => $other_name,
    -is_aligned => 'y',
    -seq_region_id => $seq_region_id,
    -seq_region_start => $seq_region_start,
    -seq_region_end => $seq_region_end,
    -seq_region_strand => $seq_region_strand
    );
 Description: Constructor
 Returntype : Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence
 Exceptions : 


=cut

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

=head2 add_track

 Arg [1]    : $evidencetrack, Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack object
 Example    : $evidence->add_track($evidencetrack);
 Description: Add the track evidence object to a list. If the flag is_last is set, so if this
              is the evidence track for the current analysis it is put to the index 0.
 Returntype : 
 Exceptions : if it's not a EvidenceTrack object


=cut

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

=head2 get_tracks

 Example    : $ra_tracks = $evidence->get_tracks;
 Description: Getter for the list of evidence tracks attached to the evidence.
              There should be one track for each analysis ran.
 Returntype : listref of Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack
 Exceptions : 


=cut

sub get_tracks {
    my $self = shift;
    return $self->{'_tracks'};
}

=head2 has_tracks

 Example    : if ($evidence->has_tracks) {...};
 Description: Return 1 if this evidence has evidence tracks attached
 Returntype : boolean
 Exceptions : 


=cut

sub has_tracks {
    my $self = shift;
    return exists $self->{'_tracks'};
}

=head2 fetch_track_by_logic_name

 Arg [1]    : $logic_name, string
 Example    : $ra_evidencetracks = $evidence->fetch_track_by_logic_name($logic_name);
 Description: Get the evidence tracks attached for the analysis given
 Returntype : listref of Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack
 Exceptions : 


=cut

sub fetch_track_by_logic_name {
    my $self = shift;
    my $logic_name = shift;

    my @a_tracks;
    foreach my $track (@{$self->{'_tracks'}}) {
        push(@a_tracks, $track) if ($track->logic_name eq $logic_name);
    }
    return \@a_tracks;
}


=head2 input_seq_id

 Arg [1]    : $input_seq_id, int [optional]
 Example    : $evidence->input_seq_id($input_seq_id);
 Description: Getter/Setter for the input sequence id
 Returntype : integer, the input sequence id
 Exceptions : 


=cut

sub input_seq_id {
  my $self = shift;
  $self->{'inputseq_id'} = shift if ( @_ );
  return $self->{'inputseq_id'};
}

=head2 hit_name

 Arg [1]    : $hit_name, int [optional]
 Example    : $evidence->hit_name($hit_name);
 Description: Getter/Setter for the hit name
 Returntype : integer, the hit name
 Exceptions : 


=cut

sub hit_name {
  my $self = shift;
  $self->{'hit_name'} = shift if ( @_ );
  return $self->{'hit_name'};
}

=head2 is_aligned

 Arg [1]    : $is_aligned, int [optional]
 Example    : $evidence->is_aligned($is_aligned);
 Description: Getter/Setter if the sequence is aligned
 Returntype : boolean
 Exceptions : 


=cut

sub is_aligned {
  my $self = shift;
  $self->{'is_aligned'} = shift if ( @_ );
  return $self->{'is_aligned'};
}

=head2 seq_region_name

 Arg [1]    : $seq_region_name, string [optional]
 Example    : $evidence->seq_region_name($seq_region_name);
 Description: Getter/Setter for the seq_region name
 Returntype : integer, the seq_region name
 Exceptions : 


=cut

sub seq_region_name {
  my $self = shift;
  $self->{'seq_region_id'} = shift if ( @_ );
  return $self->{'seq_region_id'};
}

=head2 seq_region_start

 Arg [1]    : $seq_region_start, int [optional]
 Example    : $evidence->seq_region_start($seq_region_start);
 Description: Getter/Setter for the seq_region start
 Returntype : integer, the seq_region start
 Exceptions : 


=cut

sub seq_region_start {
  my $self = shift;
  $self->{'seq_region_start'} = shift if ( @_ );
  return $self->{'seq_region_start'};
}

=head2 seq_region_end

 Arg [1]    : $seq_region_end, int [optional]
 Example    : $evidence->seq_region_end($seq_region_end);
 Description: Getter/Setter for the seq_region end
 Returntype : integer, the seq_region end
 Exceptions : 


=cut

sub seq_region_end {
  my $self = shift;
  $self->{'seq_region_end'} = shift if ( @_ );
  return $self->{'seq_region_end'};
}

=head2 seq_region_strand

 Arg [1]    : $seq_region_strand, int [optional]
 Example    : $evidence->seq_region_strand($seq_region_strand);
 Description: Getter/Setter for the seq_region strand
 Returntype : integer, the seq_region strand
 Exceptions : 


=cut

sub seq_region_strand {
  my $self = shift;
  $self->{'seq_region_strand'} = shift if ( @_ );
  return $self->{'seq_region_strand'};
}

=head2 is_stored

 Arg [1]    : $evidence_adaptor, a Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::EvidenceAdaptor object
 Example    : $evidence->is_stored($evidence_adaptor);
 Description: Test if the object is alreday stored
 Returntype : boolean
 Exceptions : 


=cut

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


