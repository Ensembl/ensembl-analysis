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

Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Storable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

@ISA = qw(Bio::EnsEMBL::Storable);


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


sub evidence_id {
  my $self = shift;
  $self->{'evidence_id'} = shift if ( @_ );
  return $self->{'evidence_id'};
}

sub analysis_run_id {
  my $self = shift;
  $self->{'analysis_run_id'} = shift if ( @_ );
  return $self->{'analysis_run_id'};
}

sub analysis_id {
  my $self = shift;
  $self->{'analysis_id'} = shift if ( @_ );
  return $self->{'analysis_id'};
}

sub reason_id {
  my $self = shift;
  $self->{'reason_id'} = shift if ( @_ );
  return $self->{'reason_id'};
}

sub input_id {
  my $self = shift;
  $self->{'input_id'} = shift if ( @_ );
  return $self->{'input_id'};
}

sub is_last {
  my $self = shift;
  $self->{'is_last'} = shift if ( @_ );
  return $self->{'is_last'};
}

#sub is_stored {
#  my $self = shift;
#  my $db = shift;
#
#  if($db and $db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
#    if ($db->get_reason_by_evidence_and_run_id($self->evidence_id, $self->analysis_run_id)) {
#        return 1;
#    }
#    return 0;
#  }
#  else {
#    throw('db argument must be a Bio::EnsEMBL::DBSQL::DBAdaptor');
#  }
#}

1;
