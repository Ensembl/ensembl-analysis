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

Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::ReasonAdaptor - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::ReasonAdaptor; 

use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Analysis::EvidenceTracking::Reason;
use Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( throw );

use vars '@ISA';
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 store

 Arg [1]    : $reason, Bio::EnsEMBL::EvidenceTracking::Reason 
 Example    : $reason_adaptor->store($reason);
 Description: Store the new reason
 Returntype : an integer, the dbID
 Exceptions : if not successful


=cut

sub store {
  my ($self, $reason) = @_;

  if (!ref $reason || !$reason->isa('Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence') ) {
    throw("Must store a Evidence object, not a $reason");
  }

  my $db = $self->db();

  my $sth = $self->prepare('INSERT into reasons ( reason_id, code, reason ) VALUES ( ?, ?, ? )');

  $sth->bind_param( 1, $reason->dbID, SQL_INTEGER );
  $sth->bind_param( 2, $reason->dbID, SQL_VARCHAR );
  $sth->bind_param( 3, $reason->info, SQL_VARCHAR );

  $sth->execute();
  $sth->finish();
  my $reason_dbID = $sth->{'mysql_insertid'};

  # set the adaptor and dbID on the original passed in reason not the
  # transfered copy
  $reason->adaptor($self);
  print STDERR "Stored reason object ".$reason->dbID."\n";
  return $reason_dbID;
}

=head2 fetch_by_dbID

 Arg [1]    : $id, integer
 Example    : $reason_adaptor->fetch_by_dbID($id);
 Description: Fetch the reason by dbID, unique element
 Returntype : Bio::EnsEMBL::EvidenceTracking::Reason
 Exceptions : 


=cut

sub fetch_by_dbID {
  my $self = shift;
  my $reason_dbID = shift;
  my $constraint = 'r.reason_id = '.$reason_dbID;
  my ($reason) = @{ $self->generic_fetch($constraint) };
  return $reason;
}

=head2 fetch_by_code

 Arg [1]    : $code, string
 Example    : $reason_adaptor->fetch_by_code($code);
 Description: Fetch the reason by code, unique element
 Returntype : Bio::EnsEMBL::EvidenceTracking::Reason
 Exceptions : 


=cut

sub fetch_by_code {
  my $self = shift;
  my $code = shift;
  my $constraint = 'r.code = '.$code;
  my ($reason) = @{ $self->generic_fetch($constraint) };
  return $reason;
}


=head2 fetch_all

 Example    : $reason_adaptor->fetch_all();
 Description: Fetch all reasons
 Returntype : listref of Bio::EnsEMBL::EvidenceTracking::Reason object
 Exceptions : 


=cut

sub fetch_all {
  my $self = shift;
  my $reason_dbID = shift;
  return $self->generic_fetch();
}



###################
# Private methods #
###################

=head2 _tables

 Example    : $self->_tables;
 Description: Return the table and its abbreviation
 Returntype : a listref of string
 Exceptions : 


=cut

sub _tables {
  return ( ['reasons' , 'r']);
}

=head2 _columns

 Example    : $self->_columns;
 Description: Return a list of columns
 Returntype : a listref of string
 Exceptions : 


=cut

sub _columns {
  return ( 'r.reason_id', 'r.code', 'r.reason');
}

=head2 _objs_from_sth

 Arg [1]    : $sth
 Example    : $self->_objs_from_sth($sth);
 Description: Put the result of the query in Bio::EnsEMBL::EvidenceTracking::Reason objects
 Returntype : listref of Bio::EnsEMBL::EvidenceTracking::Reason
 Exceptions : 


=cut


sub _objs_from_sth {
  my ($self, $sth) = @_;

  my @out;
  my ( $reason_id, $code, $reason);
  $sth->bind_columns( \$reason_id, \$code, \$reason);

  while($sth->fetch()) {
    push @out, Bio::EnsEMBL::Analysis::EvidenceTracking::Reason->new(
              -dbID    => $reason_id,
              -adaptor => $self,
              -code    => $code,
              -info    => $reason
              );
  }
  return \@out;
}


1;
