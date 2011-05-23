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

Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::EvidenceAdaptor - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::EvidenceAdaptor; 

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning stack_trace_dump );
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

use vars '@ISA';
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 store

 Arg [1]    : $evidence, Bio::EnsEMBL::EvidenceTracking::Evidence 
 Example    : $evidence_adaptor->store($evidence);
 Description: Store the new evidence
 Returntype : an integer, the dbID
 Exceptions : if not successful


=cut

sub store {
  my ($self, $evidence) = @_;

  if (!ref $evidence || !$evidence->isa('Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence') ) {
    throw("Must store a Evidence object, not a $evidence");
  }

  my $db = $self->db();

  # make an sql statement
  my $original = $evidence;

  my $sth = $self->prepare('INSERT into evidence ( input_seq_id, is_aligned ) VALUES ( (SELECT input_seq_id FROM input_seq where hit_name = ?), ? )');

  $sth->bind_param( 1, $evidence->hit_name, SQL_VARCHAR );
  $sth->bind_param( 2, $evidence->is_aligned, SQL_VARCHAR );

  $sth->execute();
  $sth->finish();
  my $evidence_dbID = $sth->{'mysql_insertid'};

  # set the adaptor and dbID on the original passed in evidence not the
  # transfered copy
  $original->adaptor($self);
  $original->dbID($evidence_dbID);
  if ($evidence->is_aligned eq 'y') {
      $sth = $self->prepare('INSERT INTO evidence_coord (evidence_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand) VALUES ( ?, (SELECT seq_region_id FROM seq_region where name = ?), ?, ?, ? )');
      $sth->bind_param( 1, $evidence_dbID, SQL_INTEGER );
      $sth->bind_param( 2, $evidence->seq_region_name, SQL_VARCHAR );
      $sth->bind_param( 3, $evidence->seq_region_start, SQL_INTEGER );
      $sth->bind_param( 4, $evidence->seq_region_end, SQL_INTEGER );
      $sth->bind_param( 5, $evidence->seq_region_strand, SQL_INTEGER );
      $sth->execute();
      $sth->finish();
  }
  print STDERR "Stored evidence object ".$original->dbID."\n";
  return $evidence_dbID;
}

=head2 fetch_by_dbID

 Arg [1]    : $id, integer
 Example    : $evidence_adaptor->fetch_by_dbID($id);
 Description: Fetch the evidence by dbID, unique element
 Returntype : Bio::EnsEMBL::EvidenceTracking::Evidence
 Exceptions : 


=cut

sub fetch_by_dbID {
  my $self = shift;
  my $evidence_dbID = shift;
  my $constraint = 'e.evidence_id = "'.$evidence_dbID.'"';
  my ($evidence) = @{ $self->generic_fetch($constraint) };
  return $evidence;
}

=head2 fetch_all_by_slice

 Arg [1]    : $slice, a slice
 Example    : $evidence_adaptor->fetch_all_by_slice($slice);
 Description: Return a list of evidence contained in the slice
 Returntype : listref of Bio::EnsEMBL::EvidenceTracking::Evidence
 Exceptions : 


=cut

sub fetch_all_by_slice {
  my $self = shift;
  my $slice = shift;

  my $constraint = 'e.seq_region_id = '.$slice->get_seq_region_id.' e.seq_region_start > '.($slice->start-1).' e.seq_region_end < '.($slice->end+1);
  return $self->generic_fetch($constraint);
}

=head2 fetch_all_by_inputseq_id

 Arg [1]    : $input_id
 Example    : $evidence_adaptor->fetch_all_by_inputseq_id($input_id);
 Description: Return a list of evidence get from the input id $input_id,
              if there was multiple run for the analysis, it returns the
              evidences for every run
 Returntype : listref of Bio::EnsEMBL::EvidenceTracking::Evidence
 Exceptions : 


=cut

sub fetch_all_by_hit_name {
  my $self = shift;
  my $hit_name = shift;

#  sub _tables { return (['evidence',' e'], ['input_seq', 'i']); }
#  sub _columns { return ( 'i.hit_name'); }
#  sub _left_join { return 'e.input_seq_id = i.input_seq_id'; }
#  my $constraint = 'i.hit_name = '.$hit_name;
  my $constraint = 'i.hit_name = "'.$hit_name.'"';
  return $self->generic_fetch($constraint);
}


#@@@@@@@
# Done @
#@@@@@@@

=head2 fetch_all_by_analysis

 Arg [1]    : $analysis_id
 Example    : $evidence_adaptor->fetch_all_by_analysis($analysis_id);
 Description: Return a list of evidence get with the analysis $analysis_id,
              if there was multiple run for the analysis, it returns the
              evidences for every run
 Returntype : listref of Bio::EnsEMBL::EvidenceTracking::Evidence
 Exceptions : 


=cut

sub fetch_all_by_analysis {
  my $self = shift;
  my $analysis_id = shift;
#  my $constraint = 'e.seq_region_id = '.$slice->get_seq_region_id.' e.seq_region_start > '.($slice->start-1).' e.seq_region_end < '.($slice->end+1);
#  return $self->generic_fetch($constraint);
}

=head2 fetch_by_name

 Arg [1]    : $hit_name, sequence name as know in the database
 Example    : $evidence_adaptor->fetch_by_name($hit_name);
 Description: Fetch evidences by hit name, can return multiple evidence
 Returntype : listref of Bio::EnsEMBL::EvidenceTracking::Evidence
 Exceptions : 


=cut

sub fetch_by_name {
}

=head2 fetch_all_by_analysis_run

 Arg [1]    : $analysis_run_id
 Example    : $evidence_adaptor->fetch_all_by_analysis_run($analysis_run_id);
 Description: Return a list of evidence get for the run $analysis_run_id of
              the analysis
 Returntype : listref of Bio::EnsEMBL::EvidenceTracking::Evidence
 Exceptions : 


=cut

sub fetch_all_by_analysis_run {
}

=head2 fetch_all_after_analysis_run

 Arg [1]    : $analysis_run_id
 Example    : $evidence_adaptor->fetch_all_after_analysis_run($analysis_run_id);
 Description: Return a list of evidence which are kept for the next analsyis
 Returntype : listref of Bio::EnsEMBL::EvidenceTracking::Evidence
 Exceptions : 


=cut

sub fetch_all_after_analysis_run {
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
  return ( ['evidence' , 'e'], ['input_seq', 'i'] );
}

=head2 _columns

 Example    : $self->_columns;
 Description: Return a list of columns
 Returntype : a listref of string
 Exceptions : 


=cut

sub _columns {
  return ( 'e.evidence_id', 'i.hit_name', 'e.is_aligned' );
}

sub _left_join {
    return (['evidence', 'e.input_seq_id = i.input_seq_id']);
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
  my ( $evidence_id, $hit_name, $is_aligned);
  $sth->bind_columns( \$evidence_id, \$hit_name, \$is_aligned);

  while($sth->fetch()) {
    push @out, Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence->new(
              -dbID              => $evidence_id,
              -adaptor           => $self,
              -hit_name          => $hit_name,
              -is_aligned        => $is_aligned
              );
    if ($out[-1]->is_aligned eq 'y') {
        my $query = 'SELECT sr.name, ec.seq_region_start, ec.seq_region_end, ec.seq_region_strand FROM evidence_coord ec LEFT JOIN seq_region sr ON sr.seq_region_id = ec.seq_region_id WHERE ec.evidence_id = '.$evidence_id;
        my $q2_sth = $self->dbc->prepare($query);
        $q2_sth->execute();
        my ($seq_region_name, $seq_region_start, $seq_region_end, $seq_region_strand) = @{$q2_sth->fetch()};
        $out[-1]->seq_region_name($seq_region_name);
        $out[-1]->seq_region_start($seq_region_start);
        $out[-1]->seq_region_end($seq_region_end);
        $out[-1]->seq_region_strand($seq_region_strand);
    }
  }
  return \@out;
}
1;
