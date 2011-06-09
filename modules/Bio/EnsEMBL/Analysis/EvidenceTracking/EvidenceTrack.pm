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

Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack - Keep track information for one evidence in one run

=head1 SYNOPSIS

use Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack;

$evidencetrack = Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack->new(
    -evidence_id     => $evidence_id,
    -analysis_run_id => $analysis_run_id,
    -analysis_id     => $self->anlysis->dbID,
    -reason_id       => $reason_id,
    -input_id        => $self->input_id,
    -is_last         => 1
    );

$evidence->add($evidencetrack);

=head1 DESCRIPTION

  Object to store the information about the run for each evidence.
  It is the central object, it connects everything together.

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

@ISA = qw(Bio::EnsEMBL::Storable);


=head2 new

 Arg [1]    : $evidence_id, int
 Arg [2]    : $analysis_run_id, int
 Arg [3]    : $analysis_id, int
 Arg [4]    : $reason_id, int
 Arg [5]    : $input_id, int
 Arg [6]    : $last, int
 Arg [7]    : $adaptor, Bio::EnsEMBL::Analysis::DBSQL::EvidenceTrackAdaptor object
 Example    : $evidencetrack = Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack->new(
    -evidence_id     => $evidence_id,
    -analysis_run_id => $analysis_run_id,
    -analysis_id     => $self->anlysis->dbID,
    -reason_id       => $reason_id,
    -input_id        => $self->input_id,
    -is_last         => 1
    );
 Description: Constructor
 Returntype : Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack
 Exceptions : 


=cut

sub new {
  my($class,@args) = @_;

  my $self = bless {},$class;

  my ($evidence_id, $analysis_run_id, $analysis_id, $reason_id, $input_id, $is_last, $adaptor) =
          rearrange([qw(EVIDENCE_ID
                        ANALYSIS_RUN_ID       
                        ANALYSIS_ID       
                        REASON_ID
                        INPUT_ID
                        IS_LAST
                        ADAPTOR
                        )],@args);

  $self->evidence_id( $evidence_id );
  $self->analysis_run_id ( $analysis_run_id );
  $self->analysis_id ( $analysis_id );
  $self->reason_id ( $reason_id );
  $self->input_id ( $input_id );
  $self->is_last ( $is_last );
  $self->adaptor ( $adaptor );
  return $self; # success - we hope!
}

=head2 clone

 Example    : $evidencetrack_clone = $evidencetrack->clone;
 Description: Clone the object
 Returntype : Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack
 Exceptions : 


=cut

sub clone {
    my $self = shift;

    my $clone = Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack->new(
        -evidence_id => $self->evidence_id,
        -analysis_run_id => $self->analysis_run_id,
        -analysis_id => $self->analysis_id,
        -reason_id => $self->reason_id,
        -input_id => $self->input_id,
        -is_last => $self->is_last,
        -adaptor => $self->adaptor
        );
    return $clone;
}


=head2 evidence_id

 Arg [1]    : $evidence_id, int [optional]
 Example    : $evidence->evidence_id($evidence_id);
 Description: Getter/Setter for the evidence id
 Returntype : integer, the evidence id
 Exceptions : 


=cut

sub evidence_id {
  my $self = shift;
  $self->{'evidence_id'} = shift if ( @_ );
  return $self->{'evidence_id'};
}

=head2 analysis_run_id

 Arg [1]    : $analysis_run_id, int [optional]
 Example    : $evidence->analysis_run_id($analysis_run_id);
 Description: Getter/Setter for the analysis run id
 Returntype : integer, the analysis run id
 Exceptions : 


=cut

sub analysis_run_id {
  my $self = shift;
  $self->{'analysis_run_id'} = shift if ( @_ );
  return $self->{'analysis_run_id'};
}

=head2 analysis_id

 Arg [1]    : $analysis_id, int [optional]
 Example    : $evidence->analysis_id($analysis_id);
 Description: Getter/Setter for the analysis id
 Returntype : integer, the analysis id
 Exceptions : 


=cut

sub analysis_id {
  my $self = shift;
  $self->{'analysis_id'} = shift if ( @_ );
  return $self->{'analysis_id'};
}

=head2 reason_id

 Arg [1]    : $reason_id, int [optional]
 Example    : $evidence->reason_id($reason_id);
 Description: Getter/Setter for the reason id
 Returntype : integer, the reason id
 Exceptions : 


=cut

sub reason_id {
  my $self = shift;
  $self->{'reason_id'} = shift if ( @_ );
  return $self->{'reason_id'};
}

=head2 input_id

 Arg [1]    : $input_id, int [optional]
 Example    : $evidence->input_id($input_id);
 Description: Getter/Setter for the input id
 Returntype : integer, the input id
 Exceptions : 


=cut

sub input_id {
  my $self = shift;
  $self->{'input_id'} = shift if ( @_ );
  return $self->{'input_id'};
}

=head2 is_last

 Arg [1]    : $is_last, int [optional]
 Example    : $evidence->is_last($is_last);
 Description: Getter/Setter if the evidence track is for the current run
 Returntype : boolean
 Exceptions : 


=cut

sub is_last {
  my $self = shift;
  $self->{'is_last'} = shift if ( @_ );
  return $self->{'is_last'};
}


1;
