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

    my $dbobj = new Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor;
    $dbobj->do_funky_db_stuff;

=head1 DESCRIPTION


=head1 METHODS

=cut

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor;


use vars qw(@ISA);
use strict;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Root;

# Inherits from the base bioperl object

@ISA = qw(Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor);


# new() inherited from Bio::EnsEMBL::DBSQL::BaseAdaptor

sub get_available_adaptors {
  my ($self) = @_;

  my $pairs = $self->SUPER::get_available_adaptors();

  $pairs->{'InputSeq'}       = 'Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::InputSeqAdaptor';
  $pairs->{'AnalysisRun'}    = 'Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::AnalysisRunAdaptor';
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

=head2 get_AnalysisRunAdaptor

 Title   : get_AnalysisRunAdaptor
 Usage   : $db->get_AnalysisRunAdaptor
 Function: The Adaptor for AnalysisRun objects in this db
 Example :
 Returns : Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::AnalysisRunAdaptor
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

=head2 get_DatabaseAdaptor

 Title   : get_DatabaseAdaptor
 Usage   : $db->get_DatabaseAdaptor
 Function: The Adaptor for Database objects in this db
 Example :
 Returns : Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DatabaseAdaptor
 Args    : nothing

=cut



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

sub get_meta_values_by_key{
    my ($self, $key, $value) = @_;
    $self->throw("No key to get supplied") unless $key;
    my $sth = $self->prepare(qq{
        SELECT meta_value
        FROM   meta
        WHERE  meta_key = ? 
    }); # ONLY RETRIEVES FIRST ENTRY,
    # SHOULD THERE BE A UNIQUE KEY ON meta_key??
    $sth->execute($key);
    my $row = $sth->fetchrow_arrayref();
    $sth->finish();
    return @{$row};
}
sub has_meta_value_by_key {
    my $self = shift;
    my ($key, $value) = @_;

    $self->throw("No key to get supplied") unless $key;
    my $sth = $self->prepare(qq{
        SELECT meta_value
        FROM   meta
        WHERE  meta_key = ? 
        AND    meta_value = ?
    }); 
    $sth->execute($key, $value);
    my $row = $sth->fetchrow_arrayref();
    $sth->finish();
    return 1 if $row;
    return 0;
}

sub update_meta_key_by_value {
    my $self = shift;
    my ($old_key, $new_key, $value) = @_;

    my $sth = $self->prepare(qq{
        UPDATE meta SET meta_key = ?
        WHERE meta_key = ?
        AND meta_value = ?
    });
    $sth->execute($new_key, $old_key, $value);
    $sth->finish();
}
1;
