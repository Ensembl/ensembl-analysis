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
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning stack_trace_dump );
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack;

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

  my $evidence = $evidence_track->evidence;
  if ($evidence->dbID eq '') {
      my $evidence_adaptor = $db->get_EvidenceAdaptor;
      $evidence_adaptor->store($evidence) unless ($evidence->is_stored($db));
  }


  my $sth = $self->prepare('INSERT into track_evidence ( evidence_id, analysis_run_id,
                            reason_id,input_id) 
                            VALUES ( ?,?,?,?)');

  $sth->bind_param( 1, $evidence->dbID, SQL_INTEGER );
  $sth->bind_param( 2, $evidence_track->analysis_run->dbID, SQL_INTEGER );
  $sth->bind_param( 3, $evidence_track->reason->dbID, SQL_INTEGER );
  $sth->bind_param( 4, $evidence_track->input_id, SQL_VARCHAR );

  $sth->execute();
  $sth->finish();
  my $evidence_track_dbID = $sth->{'mysql_insertid'};

  # set the adaptor and dbID on the original passed in evidence_track not the
  # transfered copy
#  $original->adaptor($self);
#  $original->dbID($evidence_track_dbID);
  print STDERR 'Stored evidencetrack object ',$evidence_track->evidence->dbID, "\n";
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
      
  my $constraint = 'ar.analysis_id = '.$analysis_id;
  return $self->generic_fetch($constraint);
}

=head2 fetch_by_logic_name

 Arg [1]    : $logic_name
 Example    : $analysisrun_adaptor->fetch_by_logic_name($logic_name);
 Description: Return all run for the analysis $logic_name
 Returntype : listref of Bio::EnsEMBL::EvidenceTracking::EvidenceTrack
 Exceptions : 


=cut

sub fetch_by_logic_name {
  my $self = shift;
  my $logic_name = shift;
  sub _tables {
    my $self = shift;
    return (['track_evidence' , 'te'], ['analysis_run', 'ar'], ['analysis', 'a']);
  }

  sub _left_join {
      return ['analysis_run', 'te.analysis_run_id = ar.analysis_run_id', 'analysis', 'ar.analysis_id = a.analysis_id'];
  }
      
  my $constraint = 'a.logic_name = "'.$logic_name.'"';
  return $self->generic_fetch($constraint);
}

=head2 fetch_all

 Example    : $evidencetrack_adaptor->fetch_all();
 Description: Fetch all the evidence tracks.
             WARNING!! It can take a looooong time...
 Returntype : listref of Bio::EnsEMBL::EvidenceTracking::EvidenceTrack
 Exceptions : 


=cut

sub fetch_all {
  my $self = shift;
      
  return $self->generic_fetch();
}

=head2 fetch_all_by_analysis_and_reason

 Arg [1]    : $analysis_id, int
 Arg [2]    : $reason_id, int
 Example    : @trackevidence = @{$evidence_adaptor->fetch_all_by_analysis_and_reason($analysis_id, $reason_id)};
 Description: Fetch all track evidence for an analysis with some kind of reason
 Returntype : listref of Bio::EnsEMBL::EvidenceTracking::EvidenceTrack objects
 Exceptions : 


=cut

sub fetch_all_by_analysis_and_reason {
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
           'te.input_id');
}

=head2 _left_join

   Example    : $self->_left_join;
   Description: Return a list of left joins
   Returntype : a listref of string


=cut

sub _left_join {
    return ['analysis_run', 'te.analysis_run_id = ar.analysis_run_id'];
}

=head2 _objs_from_sth

 Arg [1]    : $sth
 Example    : $self->_objs_from_sth($sth);
 Description: Put the result of the query in Bio::EnsEMBL::EvidenceTracking::EvidenceTrack objects
 Returntype : listref of Bio::EnsEMBL::EvidenceTracking::EvidenceTrack
 Exceptions : 


=cut


sub _objs_from_sth {
  my ($self, $sth) = @_;

  my @out;
  my ( $evidence_id, $analysis_run_id, $reason_id, $input_id);
  $sth->bind_columns( \$evidence_id, \$analysis_run_id, \$reason_id, \$input_id);

  my $reason_adaptor = $self->db->get_ReasonAdaptor;
  my %h_reasons;
  foreach my $reason (@{$reason_adaptor->fetch_all}) {
    $h_reasons{$reason->dbID} = $reason;
  }
  while($sth->fetch()) {
    my $evidence = $self->_cache_by_evidence_id($evidence_id);
    my $analysis_run = $self->_cache_by_analysis_run_id($analysis_run_id);
    
    push @out, Bio::EnsEMBL::Analysis::EvidenceTracking::EvidenceTrack->new(
              -adaptor      => $self,
              -evidence     => $evidence,
              -analysis_run => $analysis_run,
              -reason       => $h_reasons{$reason_id},
              -input_id     => $input_id,
              );
  }
#  my @out = sort {$a->analysis_run->dbID <=> $b->analysis_run->dbID} @tmp;
  return \@out;
}

sub _cache_by_evidence_id {
    my $self = shift;
    my $evidence_id = shift;

    $self->{'_evidence_cache'}->{$evidence_id} = $self->db->get_EvidenceAdaptor->fetch_by_dbID($evidence_id)
        unless (exists $self->{'_evidence_cache'}->{$evidence_id});
    return $self->{'_evidence_cache'}->{$evidence_id};
}

sub _cache_by_analysis_run_id {
    my $self = shift;
    my $analysis_run_id = shift;

    $self->{'_analysis_run_cache'}->{$analysis_run_id} = $self->db->get_AnalysisRunAdaptor->fetch_by_dbID($analysis_run_id)
        unless (exists $self->{'_analysis_run_cache'}->{$analysis_run_id});
    return $self->{'_analysis_run_cache'}->{$analysis_run_id};
}


1;
