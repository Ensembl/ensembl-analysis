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

Bio::EnsEMBL::Analysis::EvidenceTracking::InputSeq - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::EvidenceTracking::InputSeq;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Storable;
use Data::Dumper;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
#use Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::InputSeqAdaptor;

@ISA = qw(Bio::EnsEMBL::Storable);


sub new {
  my($class,@args) = @_;

  my $self = bless {},$class;

  my ($id, $hit_name, $external_db_id, 
      $molecule_type, $submission_date, $adaptor) =
          rearrange([qw(DBID
                        HIT_NAME       
                        EXTERNAL_DB_ID
                        MOLECULE_TYPE
                        SUBMISSION_DATE
                        ADAPTOR
                        )],@args);

  $self->dbID      ( $id ) if (defined $id);
  $self->hit_name ( $hit_name );
  $self->external_db_id ( $external_db_id );
  $self->molecule_type ( $molecule_type );
  $self->submission_date ( $submission_date );
  $self->adaptor   ( $adaptor );
  return $self; # success - we hope!
}


sub hit_name {
  my $self = shift;
  $self->{'hit_name'} = shift if ( @_ );
  return $self->{'hit_name'};
}

sub external_db_id {
  my $self = shift;
  $self->{'external_db_id'} = shift if ( @_ );
  return $self->{'external_db_id'};
}

sub molecule_type {
  my $self = shift;
  $self->{'molecule_type'} = shift if ( @_ );
  return $self->{'molecule_type'};
}

sub submission_date {
  my $self = shift;
  $self->{'submission_date'} = shift if ( @_ );
  return $self->{'submission_date'};
}

sub is_stored {
  my $self = shift;
  my $db = shift;

  if($db and $db->isa('Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::InputSeqAdaptor')) {
    if ($db->fetch_by_hit_name($self->hit_name)) {
        print STDERR $self->hit_name, ' is stored!',"\n";
        return 1;
    }
    return 0;
  }
  else {
    throw('db argument must be a Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::InputSeqAdaptor');
  }
}

1;
