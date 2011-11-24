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

 Arg [1]    : $evidence, Bio::EnsEMBL::Analysis::DBSQL::Evidence object
 Arg [2]    : $analysis_run, Bio::EnsEMBL::Analysis::DBSQL::AnalysisRun object
 Arg [3]    : $reason, Bio::EnsEMBL::Analysis::DBSQL::Reason object
 Arg [4]    : $input_id, int
 Arg [5]    : $adaptor, Bio::EnsEMBL::Analysis::DBSQL::EvidenceTrackAdaptor object
 Example    : $evidencetrack = Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack->new(
    -evidence     => $evidence,
    -analysis_run => $analysis_run,
    -reason       => $reason,
    -input_id     => $input_id,
    );
 Description: Constructor
 Returntype : Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack
 Exceptions : 


=cut

sub new {
  my($class,@args) = @_;

  my $self = bless {},$class;

  my ($evidence, $analysis_run, $reason, $input_id, $adaptor) =
          rearrange([qw(EVIDENCE
                        ANALYSIS_RUN       
                        REASON
                        INPUT_ID
                        ADAPTOR
                        )],@args);

  $self->evidence( $evidence);
  $self->analysis_run( $analysis_run);
  $self->reason( $reason);
  $self->input_id ( $input_id );
  $self->adaptor ( $adaptor );
  return $self; # success - we hope!
}


=head2 evidence_id

 Arg [1]    : $evidence_id, int [optional]
 Example    : $evidence->evidence_id($evidence_id);
 Description: Getter/Setter for the evidence id
 Returntype : integer, the evidence id
 Exceptions : 


=cut

sub evidence {
  my $self = shift;
  my $evidence = shift;
  if ($evidence) {
      throw('Need to pass an Evidence object') unless $evidence->isa('Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence');
      $self->{'evidence'} = $evidence;
  }
  return $self->{'evidence'};
}

=head2 analysis_run_id

 Arg [1]    : $analysis_run_id, int [optional]
 Example    : $evidence->analysis_run_id($analysis_run_id);
 Description: Getter/Setter for the analysis run id
 Returntype : integer, the analysis run id
 Exceptions : 


=cut

sub analysis_run {
  my $self = shift;
  my $analysis_run = shift;
  if ($analysis_run) {
      throw('Need to pass an AnalysisRun object') unless $analysis_run->isa('Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun');
      $self->{'analysis_run'} = $analysis_run;
  }
  return $self->{'analysis_run'};
}

=head2 reason_id

 Arg [1]    : $reason_id, int [optional]
 Example    : $evidence->reason_id($reason_id);
 Description: Getter/Setter for the reason id
 Returntype : integer, the reason id
 Exceptions : 


=cut

sub reason {
  my $self = shift;
  $self->{'reason'} = shift if ( @_ );
  return $self->{'reason'};
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

=head2 name

 Example    : $name = $evidencetrack->name;
 Description: Get the name of the sequence used as supporting evidence
 Returntype : String
 Exceptions : 


=cut

sub name {
    my $self = shift;

    return $self->evidence->name;
}

=head2 code

 Example    : $code = $evidencetrack->code;
 Description: Get the the reason code for the track
 Returntype : String
 Exceptions : 


=cut

sub code {
    my $self = shift;

    return $self->reason->code;
}


1;
