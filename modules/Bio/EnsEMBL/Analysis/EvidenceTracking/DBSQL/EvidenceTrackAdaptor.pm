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

Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::EvidenceTrackAdaptor - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::EvidenceTrackAdaptor; 

use strict;
use Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning stack_trace_dump );
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

use vars '@ISA';
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 store

 Arg [1]    : $analysis_run, Bio::EnsEMBL::EvidenceTracking::EvidenceTrack
 Example    : $analysisrun_adaptor->store($analysis_run);
 Description: Store the new run of the analysis
 Returntype : an integer, the dbID
 Exceptions : if not successful


=cut

sub store {
  my ($self, $evidence_track) = @_;

  if (!ref $evidence_track || !$evidence_track->isa('Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack') ) {
    throw('Must store a EvidenceTrack object, not a '.$evidence_track);
  }

  my $db = $self->db();

  # make an sql statement
  my $original = $evidence_track;

  my $sth = $self->prepare('INSERT into track_evidence ( evidence_id, analysis_run_id,
                            reason_id,input_id) 
                            VALUES ( ?,?,?,?)');

  $sth->bind_param( 1, $evidence_track->evidence_id, SQL_INTEGER );
  $sth->bind_param( 2, $evidence_track->analysis_run_id, SQL_INTEGER );
  $sth->bind_param( 3, $evidence_track->reason_id, SQL_INTEGER );
  $sth->bind_param( 4, $evidence_track->input_id, SQL_VARCHAR );

  $sth->execute();
  $sth->finish();
  my $evidence_track_dbID = $sth->{'mysql_insertid'};

  # set the adaptor and dbID on the original passed in evidence_track not the
  # transfered copy
  $original->adaptor($self);
  $original->dbID($evidence_track_dbID);
  print STDERR 'Stored evidencetrack object ',$evidence_track->evidence_id, "\n";
  return $evidence_track_dbID;
}

=head2 fetch_all_by_evidence

 Arg [1]    : $evidence
 Example    : $analysisrun_adaptor->fetch_all_by_evidence($evidence);
 Description: Return a list of track for this evidence
 Returntype : listref of Bio::EnsEMBL::EvidenceTracking::EvidenceTrack
 Exceptions : 


=cut

sub fetch_all_by_evidence {
  my $self = shift;
  my $evidence = shift;
  my $constraint = 'te.evidence_id = "'.$evidence->dbID.'"';
  return $self->generic_fetch($constraint);
}

=head2 fetch_all_by_reason

 Arg [1]    : $reason
 Example    : $analysisrun_adaptor->fetch_all_by_reason($reason);
 Description: Return a list of track for this reason
 Returntype : listref of Bio::EnsEMBL::EvidenceTracking::EvidenceTrack
 Exceptions : 


=cut

sub fetch_all_by_reason {
  my $self = shift;
  my $reason = shift;
  my $constraint = 'te.reason_id = "'.$reason->dbID.'"';
  return $self->generic_fetch($constraint);
}

=head2 fetch_all_by_input_id

 Arg [1]    : $input_id
 Example    : $analysisrun_adaptor->fetch_all_by_input_id($input_id);
 Description: Return a list of track for this input_id
 Returntype : listref of Bio::EnsEMBL::EvidenceTracking::EvidenceTrack
 Exceptions : 


=cut

sub fetch_all_by_input_id {
  my $self = shift;
  my $input_id = shift;
  my $constraint = 'te.input_id = "'.$input_id.'"';
  return $self->generic_fetch($constraint);
}

=head2 fetch_by_analysis

 Arg [1]    : $analysis_id
 Example    : $analysisrun_adaptor->fetch_by_analysis($analysis_id);
 Description: Return all run for the analysis $analysis_id
 Returntype : listref of Bio::EnsEMBL::EvidenceTracking::EvidenceTrack
 Exceptions : 


=cut

sub fetch_all_by_analysis {
  my $self = shift;
  my $analysis_id = shift;
      
  my $constraint = 'ar.analysis_id = "'.$analysis_id.'"';
  return $self->generic_fetch($constraint);
}

sub fetch_all_by_analysis {
  my $self = shift;
  my ($analysis_id, $reason_id) = @_;
      
  my $constraint = 'ar.analysis_id = "'.$analysis_id.'" AND te.reason_id = '.$reason_id;
  return $self->generic_fetch($constraint);
}

=head2 fetch_by_analysis_run_id

 Arg [1]    : $analysis_run_id
 Example    : $analysisrun_adaptor->fetch_by_analysis_run_id($analysis_run_id);
 Description: Return all run for the analysis $analysis_run_id
 Returntype : listref of Bio::EnsEMBL::EvidenceTracking::EvidenceTrack
 Exceptions : 


=cut

sub fetch_by_analysis_run_id {
  my $self = shift;
  my $analysis_run_id = shift;
  my $constraint = 'te.analysis_run_id = '.$analysis_run_id;
  return $self->generic_fetch($constraint);
}

#@@@@@@@
# Done @
#@@@@@@@

=head2 fetch_all_unmapped_evidence

 Example    : $evidence_adaptor->fetch_all_unmapped_evidence;
 Description: Return a list of sequence that have not been mapped to the genome
 Returntype : listref of Bio::EnsEMBL::EvidenceTracking::Evidence
 Exceptions : 


=cut

sub fetch_all_unmapped_evidence {
}

=head2 fetch_all_unmapped_evidence_by_analysis_run

 Arg [1]    : $analysis_run_id, int
 Example    : $evidence_adaptor->fetch_all_unmapped_evidence_by_analysis_run($analysis_run_id);
 Description: Return a list of sequence that have not been mapped to the genome
              for a specific run
 Returntype : listref of Bio::EnsEMBL::EvidenceTracking::Evidence
 Exceptions : 


=cut

sub fetch_all_unmapped_evidence_by_analysis_run {
}

=head2 fetch_all_unmapped_evidence_by_reason

 Arg [1]    : $reason_id
 Example    : $evidence_adaptor->fetch_all_unmapped_evidence_by_reason($reason_id);
 Description: Return a list of evidence that have not been mapped for a reason
 Returntype : listref of Bio::EnsEMBL::EvidenceTracking::Evidence
 Exceptions : 


=cut

sub fetch_all_unmapped_evidence_by_reason {
}



###################
# Private methods #
###################

=head2 _tables

 Example    : $self->_tables;
 Description: Return the table and its abbreviation
 Returntype : a listref of string


=cut

sub _tables {
  my $self = shift;
  return (['track_evidence' , 'te'], ['analysis_run', 'ar']);
}

=head2 _columns

 Example    : $self->_columns;
 Description: Return a list of columns
 Returntype : a listref of string


=cut

sub _columns {
  my $self = shift;
  return ( 'te.evidence_id', 'te.analysis_run_id', 'te.reason_id',
           'te.input_id', 'ar.analysis_id');
}

=head2 _left_join

   Example    : $self->_left_join;
   Description: Return the left join
   Returntype : a string


=cut

sub _left_join {
    return ['analysis_run', 'te.analysis_run_id = ar.analysis_run_id'];
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

  my @tmp;
  my ( $evidence_id, $analysis_run_id, $reason_id, $input_id, $analysis_id);
  $sth->bind_columns( \$evidence_id, \$analysis_run_id, \$reason_id, \$input_id, \$analysis_id);

  while($sth->fetch()) {
    push @tmp, Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack->new(
              -adaptor         => $self,
              -evidence_id     => $evidence_id,
              -analysis_run_id => $analysis_run_id,
              -reason_id       => $reason_id, 
              -input_id        => $input_id,
              -analysis_id     => $analysis_id,
              -is_last         => 0,
              );
  }
  my @out = sort {$a->analysis_run_id <=> $b->analysis_run_id} @tmp;
  $out[-1]->is_last(1) unless (scalar(@out) == 0);
  return \@out;
}
1;
