package Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DatabaseAdaptor; 

use strict;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning stack_trace_dump );
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Analysis::EvidenceTracking::Database;
use Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 new

  Args       : Bio::EnsEMBL::DBSQL::DBAdaptor
  Example    : my $aa = new Bio::EnsEMBL::Pipeline::DBSQL::PipelineAdaptor();
  Description: Creates a new Bio::EnsEMBL::Pipeline::DBSQL::PipelineAdaptor object and
               internally loads and caches all the Pipeline objects from the 
               database.
  Returntype : Bio::EnsEMBL::Pipeline::DBSQL::PipelineAdaptor
  Exceptions : none
  Caller     : Bio::EnsEMBL::DBSQL::DBAdaptor

=cut

sub new {
  my ($class, $db) = @_;
 
  my $self = $class->SUPER::new($db);
 
  #load and cache all of the Pipeline objects
#  $self->fetch_all;

  return $self;
}

=head2 store

 Arg [1]    : $analysis_run, Bio::EnsEMBL::Pipeline::Database
 Example    : $analysisrun_adaptor->store($analysis_run);
 Description: Store the new run of the analysis
 Returntype : an integer, the dbID
 Exceptions : if not successful


=cut

sub store {
  my ($self, $database) = @_;

  if (!ref $database || !$database->isa('Bio::EnsEMBL::Analysis::EvidenceTracking::Database') ) {
    throw("Must store a Database object, not a $database");
  }

  my $db = $self->db();

  # make an sql statement
  my $original = $database;

  my $sth = $self->prepare("INSERT into dbs ( db_name, instance ) 
                            VALUES ( ?,? )");

  $sth->bind_param( 1, $database->db_name, SQL_VARCHAR );
  $sth->bind_param( 2, $database->instance, SQL_VARCHAR );

  $sth->execute();
  $sth->finish();
  my $database_dbID = $sth->{'mysql_insertid'};

  # set the adaptor and dbID on the original passed in analysis_run not the
  # transfered copy
#  $original->adaptor($self);
#  $original->dbID($database_dbID);
  $database->adaptor($self);
  $database->dbID($database_dbID);
  print STDERR "Stored inputseq object ".$original->dbID."\n";
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



#@@@@@@@
# Done @
#@@@@@@@


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
