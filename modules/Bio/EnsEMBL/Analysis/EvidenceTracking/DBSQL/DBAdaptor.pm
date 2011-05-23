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

Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

#
# Object for storing the connection to the analysis database
#
# Written by Simon Potter <scp@sanger.ac.uk>
# Based on Michele Clamp's Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::Obj
#
# Copyright GRL/EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code


=pod

=head1 NAME

Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor -
adapter class for EnsEMBL Analysis DB

=head1 SYNOPSIS

    my $dbobj = new Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor;
    $dbobj->do_funky_db_stuff;

=head1 DESCRIPTION

Interface for the connection to the analysis database

=head1 CONTACT

Post general queries to B<ensembl-dev@ebi.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor;


use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Root;

# Inherits from the base bioperl object

@ISA = qw(Bio::EnsEMBL::DBSQL::DBAdaptor);


# new() inherited from Bio::EnsEMBL::DBSQL::BaseAdaptor

sub get_available_adaptors {
  my ($self) = @_;

  my $pairs = $self->SUPER::get_available_adaptors();

  $pairs->{'InputSeq'}       = 'Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::InputSeqAdaptor';
  $pairs->{'AnalysisRun'}       = 'Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::AnalysisRunAdaptor';
  $pairs->{'EvidenceTrack'}  = 'Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::EvidenceTrackAdaptor';
  $pairs->{'Evidence'}       = 'Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::EvidenceAdaptor';
  $pairs->{'Database'}       = 'Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DatabaseAdaptor';

  return $pairs; 
}

=head2 get_InputSeqAdaptor

 Title   : get_InputSeqAdaptor
 Usage   : $db->get_InputSeqAdaptor
 Function: The Adaptor for InputSeq objects in this db
 Example :
 Returns : Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::InputSeqAdaptor
 Args    : nothing

=cut

=head2 get_EvidenceTrackAdaptor

 Title   : get_EvidenceTrackAdaptor
 Usage   : $db->get_EvidenceTrackAdaptor
 Function: The Adaptor for EvidenceTrack objects in this db
 Example :
 Returns : Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::EvidenceTrackAdaptor
 Args    : nothing


=cut

=head2 get_EvidenceAdaptor

 Title   : get_EvidenceAdaptor
 Usage   : $db->get_EvidenceAdaptor
 Function: The Adaptor for Evidence objects in this db
 Example :
 Returns : Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::EvidenceAdaptor
 Args    : nothing


=cut


=head2 _db_handle

 Title   : _db_handle
 Usage   : $sth = $dbobj->_db_handle($dbh);
 Function: Get/set method for the database handle
 Example :
 Returns : A database handle object
 Args    : A database handle object

=cut

sub _db_handle {
    my ($self,$arg) = @_;

    if (defined($arg)) {
	$self->{_db_handle} = $arg;
    }
    return $self->{_db_handle};
}


sub prepare {
  my ($self, @args) = @_;

  $self->dbc->prepare(@args);
}


=head2 DESTROY

 Title   : DESTROY
 Usage   :
 Function:
 Example :
 Returns :
 Args    :

=cut

sub DESTROY {
   my ($obj) = @_;

   if( $obj->{'_db_handle'} ) {
       $obj->{'_db_handle'}->disconnect;
       $obj->{'_db_handle'} = undef;
   }
}

=head2 Some Utility Stuff

 Access to the meta table of the schema including
 
 get_meta_value_by_key     - retrieve value by key
 store_meta_key_value      - write key, value pair
 remove_meta_key           - delete by key name
 make_meta_value_from_hash - flatten a hash to a string (uses dbi->quote to escape)
 make_hash_from_meta_value - returns a hash from a previously flattened string
 
 
    $dbA->store_meta_key_value('my_key', 'the value');

    my $value = $dbA->get_meta_value_by_key('my_key');
    
    $dbA->remove_meta_key('my_key');

    my %hash = ('-host' => 'pfam', '-port' => '3306'...);
    my $flat = $dbA->make_meta_value_from_hash(\%hash);
    my %retrieved = %{$dbA->make_hash_from_meta_value($flat)};

=cut

sub get_meta_value_by_key{
    my ($self, $key, $value) = @_;
    $self->throw("No key to get supplied") unless $key;
    my $sth = $self->prepare(qq{
        SELECT meta_value
        FROM   meta
        WHERE  meta_key = ? LIMIT 1
    }); # ONLY RETRIEVES FIRST ENTRY,
    # SHOULD THERE BE A UNIQUE KEY ON meta_key??
    $sth->execute($key);
    my $row = $sth->fetchrow_arrayref();
    $value = $row->[0] if $row;
    $sth->finish();
    return $value;
}
sub store_meta_key_value{
    my ($self, $key, $value) = @_;
    $self->throw("No key|value to insert supplied") unless $key && $value;
    my $sth = $self->prepare(qq{
        INSERT INTO meta (meta_key, meta_value) VALUES (?, ?)
    });
    $sth->execute($key, $value);
    $sth->finish();
    return undef;
}
sub remove_meta_key{
    my ($self, $key) = @_;
    $self->throw("No key to remove supplied") unless $key;
    my $sth = $self->prepare(qq{
	DELETE
	FROM   meta
	WHERE  meta_key = ?
    });
    $sth->execute($key);
    $sth->finish;
    return undef;
}
sub make_meta_value_from_hash{
    my ($self, $hash) = @_;
    my $dbh = $self->_db_handle();
    return join(",\n", map{ $dbh->quote($_)." => ".$dbh->quote($hash->{$_}) } keys(%$hash));
}
sub make_hash_from_meta_value{
    my ($self,$string) = @_;
    if($string){
        my $hash = { eval $string };
        $@ ? die "error evaluating $string" : return $hash || {};
    }
    return {};
}
1;
