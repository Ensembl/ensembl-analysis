package Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::InputSeqAdaptor; 

use strict;
use Bio::EnsEMBL::Analysis::EvidenceTracking::InputSeq;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning stack_trace_dump );
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Data::Dumper;

use vars '@ISA';
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 fetch_by_dbID

 Arg [1]    : $dbID
 Example    : $inputseq_adaptor->fetch_by_dbID($dbID);
 Description: Return a sequence
 Returntype : Bio::EnsEMBL::EvidenceTracking::InputSeq
 Exceptions : 


=cut

sub fetch_by_dbID {
  my $self = shift;
  my $inputseq_id = shift;
  my $constraint = "i.input_seq_id = '$inputseq_id'";
  my ($inputseq_obj) = @{ $self->generic_fetch($constraint) };
  return $inputseq_obj;
}

=head2 fetch_by_hit_name

 Arg [1]    : $hit_name, a string
 Example    : $inputseq_adaptor->fetch_by_hit_name($hit_name);
 Description: Fetch a sequence by its name
 Returntype : Bio::EnsEMBL::EvidenceTracking::InputSeq
 Exceptions : 


=cut

sub fetch_by_hit_name {
  my $self = shift;
  my $hit_name = shift;
  my $constraint = "i.hit_name = '$hit_name'";
  my ($inputseq_obj) = @{ $self->generic_fetch($constraint) };
  return $inputseq_obj;
}

=head2 fetch_all_by_molecule_type

 Arg [1]    : $type, a string between 'protein', 'cdna', 'est'
 Example    : $inputseq_adaptor->fetch_all_by_molecule_type($type);
 Description: Fetch all the molecules from a type
 Returntype : listref of Bio::EnsEMBL::EvidenceTracking::InputSeq
 Exceptions : 


=cut

sub fetch_all_by_molecule_type {
  my $self = shift;
  my $molecule_type = shift;
  throw("molecule_type is required") unless ($molecule_type);

  # fetch list of inputseq_ids
  my $sth = $self->prepare("SELECT i.input_seq_id 
                            FROM input_seq i
                            WHERE molecule_type = ?");
  $sth->bind_param(1, $molecule_type, SQL_VARCHAR);
  $sth->execute();
  my @array = @{$sth->fetchall_arrayref()};
  $sth->finish();

  # fetch all inputseq objects
  my @inputseq_ids = map {$_->[0]} @array;
  my $inputseq_objs = $self->fetch_all_by_dbID_list(\@inputseq_ids);

  return $inputseq_objs;
}

=head2 fetch_all_by_external_db_id

 Arg [1]    : $external_db_id, an int
 Example    : $inputseq_adaptor->fetch_all_by_external_db_id($external_db_id);
 Description: Fetch all molecule from an external database (Uniprot,...)
 Returntype : listref of Bio::EnsEMBL::EvidenceTracking::InputSeq
 Exceptions : 


=cut

sub fetch_all_by_external_db_id {
  my $self = shift;
  my $external_db_id = shift;
  throw("external_db_id is required") unless ($external_db_id);

  # fetch list of inputseq_ids
  my $sth = $self->prepare("SELECT i.input_seq_id 
                            FROM input_seq i
                            WHERE external_db_id = ?");
  $sth->bind_param(1, $external_db_id, SQL_INTEGER);
  $sth->execute();
  my @array = @{$sth->fetchall_arrayref()};
  $sth->finish();

  # fetch all inputseq objects
  my @inputseq_ids = map {$_->[0]} @array;
  my $inputseq_objs = $self->fetch_all_by_dbID_list(\@inputseq_ids);

  return $inputseq_objs;
}

=head2 fetch_all_by_dbID_list

 Arg [1]    : $ra_dbIDs, a listref of dbID
 Example    : $inputseq_adaptor->fetch_all_by_dbID_list($ra_dbIDs);
 Description: Fetch a list of molecule via their dbID
 Returntype : listref of Bio::EnsEMBL::EvidenceTracking::InputSeq
 Exceptions : 


=cut


sub fetch_all_by_dbID_list {
  my ($self,$id_list_ref) = @_;

  if(!defined($id_list_ref) || ref($id_list_ref) ne 'ARRAY') {
    throw("inputseq id list reference argument is required");
  }

  return [] if(!@$id_list_ref);

  my @out;
  #construct a constraint like 't1.table1_id = 123'
  my @tabs = $self->_tables;
  my ($name, $syn) = @{$tabs[0]};

  # mysql is faster and we ensure that we do not exceed the max query size by
  # splitting large queries into smaller queries of 200 ids
  my $max_size = 200;
  my @id_list = @$id_list_ref;

  while(@id_list) {
    my @ids;
    if(@id_list > $max_size) {
      @ids = splice(@id_list, 0, $max_size);
    } else {
      @ids = splice(@id_list, 0);
    }

    my $id_str;
    if(@ids > 1)  {
      $id_str = " IN (" . join(',', @ids). ")";
    } else {
      $id_str = " = " . $ids[0];
    }

    my $constraint = "${syn}.${name}_id $id_str";
    push @out, @{$self->generic_fetch($constraint)};
  }
  return \@out;
}

=head2 store

 Arg [1]    : $input_seq, a Bio::EnsEMBL::EvidenceTracking::InputSeq
 Example    : $inputseq_adaptor->store($input_seq);
 Description: Store a molecule
 Returntype : int, the dbID of the molecule
 Exceptions : 


=cut


sub store {
  my ($self, $inputseq_obj) = @_;

  if (!ref $inputseq_obj || !$inputseq_obj->isa('Bio::EnsEMBL::Analysis::EvidenceTracking::InputSeq') ) {
    throw("Must store a InputSeq object, not a $inputseq_obj");
  }

  if ($inputseq_obj->is_stored($self)) {
    print STDERR "Already stored this InputSeq ".$inputseq_obj->hit_name."\n"; 
    return $self->fetch_by_hit_name($inputseq_obj->hit_name)->dbID;
  }

  # make an sql statement
  my $original = $inputseq_obj;

  my $sth = $self->prepare("INSERT into input_seq ( hit_name,
                            external_db_id,molecule_type) 
                            VALUES ( ?,?,?)");

  $sth->bind_param( 1, $inputseq_obj->hit_name, SQL_VARCHAR );
  $sth->bind_param( 2, $inputseq_obj->external_db_id, SQL_INTEGER );
  $sth->bind_param( 3, $inputseq_obj->molecule_type, SQL_VARCHAR );

  $sth->execute();
  $sth->finish();
  my $inputseq_obj_dbID = $sth->{'mysql_insertid'};

  # We may not have a submission date
  my $date_sth = $self->prepare("UPDATE input_seq set submission_date = ? 
                                WHERE input_seq_id = ?");

  $date_sth->bind_param( 1, $inputseq_obj->submission_date, SQL_DATETIME);
  $date_sth->bind_param( 2, $inputseq_obj_dbID, SQL_INTEGER );
  $date_sth->execute();
  $date_sth->finish();

  # set the adaptor and dbID on the original passed in inputseq_obj not the
  # transfered copy
  $original->adaptor($self);
  $original->dbID($inputseq_obj_dbID);
  print STDERR "Stored inputseq object ".$original->dbID."\n";
  return $inputseq_obj_dbID;
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
  return (['input_seq' , 'i']);
}

=head2 _columns

 Example    : $self->_columns;
 Description: Return a list of columns
 Returntype : a listref of string
 Exceptions : 


=cut

sub _columns {
  my $self = shift;
  return ( 'i.input_seq_id', 'i.hit_name', 'i.external_db_id',
           'i.molecule_type', 'i.submission_date');
}

=head2 _objs_from_sth

 Arg [1]    : $sth
 Example    : $self->_objs_from_sth($sth);
 Description: Put the result of the query in Bio::EnsEMBL::EvidenceTracking::InputSeq objects
 Returntype : listref of Bio::EnsEMBL::EvidenceTracking::InputSeq
 Exceptions : 


=cut


sub _objs_from_sth {
  my ($self, $sth) = @_;

  my @out;
  my ( $inputseq_id, $hit_name, $external_db_id, $molecule_type, $submission_date);
  $sth->bind_columns( \$inputseq_id, \$hit_name, \$external_db_id, \$molecule_type, \$submission_date);

  while($sth->fetch()) {
    push @out, Bio::EnsEMBL::Analysis::EvidenceTracking::InputSeq->new(
              -dbID             => $inputseq_id,
              -adaptor          => $self,
              -hit_name         => $hit_name,
              -external_db_id   => $external_db_id,
              -molecule_type    => $molecule_type, 
              -submission_date  => $submission_date
              );
  }
  return \@out;
}

1;
