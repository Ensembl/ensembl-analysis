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

Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DatabaseAdaptor - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DatabaseAdaptor; 

use strict;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning stack_trace_dump );
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Analysis::EvidenceTracking::Database;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 store

 Arg [1]    : $analysis_run, Bio::EnsEMBL::Pipeline::Database
 Example    : $analysisrun_adaptor->store($analysis_run);
 Description: Store the new run of the analysis
 Returntype : an integer, the dbID
 Exceptions : if not successful


=cut

sub store {
  my $self = shift;
  my $database = shift;

  if (!ref $database || !$database->isa('Bio::EnsEMBL::Analysis::EvidenceTracking::Database') ) {
    throw('Must store a Bio::EnsEMBL::Analysis::EvidenceTracking::Database object, not a '.ref($database));
  }

  my $db = $self->db();

  my $sth = $self->prepare("INSERT into dbs ( db_name, instance ) 
                            VALUES ( ?,? )");

  $sth->bind_param( 1, $database->db_name, SQL_VARCHAR );
  $sth->bind_param( 2, $database->instance, SQL_VARCHAR );

  $sth->execute();
  $sth->finish();
  my $database_dbID = $sth->{'mysql_insertid'};

  $database->adaptor($self);
  $database->dbID($database_dbID);
  print STDERR "Stored inputseq object ".$database->dbID."\n";
  return $database_dbID;
}

=head2 fetch_by_dbID

 Arg [1]    : $id, integer
 Example    : $analysisrun_adaptor->fetch_by_dbID($id);
 Description: Fetch the evidence by dbID, unique element
 Returntype : Bio::EnsEMBL::Pipeline::Database
 Exceptions : 


=cut

sub fetch_by_dbID {
  my $self = shift;
  my $database_id = shift;
  my $constraint = "d.db_id = '$database_id'";
  my ($database) = @{ $self->generic_fetch($constraint) };
  return $database;
}

=head2 fetch_DB

 Arg [1]    : $id, integer
 Example    : $analysisrun_adaptor->fetch_by_dbID($id);
 Description: Fetch the evidence by dbID, unique element
 Returntype : Bio::EnsEMBL::Pipeline::Database
 Exceptions : 


=cut

sub fetch_DB {
  my $self = shift;
  my ($db_name, $instance) = @_;

  my $constraint = "d.db_name = '$db_name' AND d.instance = '$instance'";
  my ($database) = @{ $self->generic_fetch($constraint) };
  return $database;
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
  my $self = shift;
  return (['dbs' , 'd']);
}

=head2 _columns

 Example    : $self->_columns;
 Description: Return a list of columns
 Returntype : a listref of string
 Exceptions : 


=cut

sub _columns {
  my $self = shift;
  return ( 'd.db_id', 'd.db_name', 'd.instance' );
}

=head2 _objs_from_sth

 Arg [1]    : $sth
 Example    : $self->_objs_from_sth($sth);
 Description: Put the result of the query in Bio::EnsEMBL::Pipeline::InputSeq objects
 Returntype : listref of Bio::EnsEMBL::Pipeline::InputSeq
 Exceptions : 


=cut

sub _objs_from_sth {
  my ($self, $sth) = @_;

  my @out;
  my ( $db_id, $db_name, $instance );
  $sth->bind_columns( \$db_id, \$db_name, \$instance);

  while($sth->fetch()) {
    push @out, Bio::EnsEMBL::Analysis::EvidenceTracking::Database->new(
              -dbID     => $db_id,
              -adaptor  => $self,
              -db_name  => $db_name,
              -instance => $instance,
              );
  }
  return \@out;
}


1;
