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

  # I can also lock the table but I'm afraid it will just make the analysis longer for not a big deal... So I'm limiting
  # People will have to remove the duplicates
  my $sth = $self->prepare('INSERT into evidence ( input_seq_id, is_aligned ) VALUES ( (SELECT input_seq_id FROM input_seq where hit_name = ? LIMIT 1), ? )');

  $sth->bind_param( 1, $evidence->input_seq->hit_name, SQL_VARCHAR );
  $sth->bind_param( 2, $evidence->is_aligned, SQL_VARCHAR );

  $sth->execute();
  $sth->finish();
  my $evidence_dbID = $sth->{'mysql_insertid'};

  # set the adaptor and dbID on the original passed in evidence not the
  # transfered copy
  $original->adaptor($self);
  $original->dbID($evidence_dbID);
  if ($evidence->is_aligned eq 'y') {
      print STDERR 'INSERT: ', $evidence_dbID, ' ', $evidence->seq_region_name, ' ',$evidence->seq_region_start, ' ',$evidence->seq_region_end, ' ',$evidence->seq_region_strand, "\n";
      $sth = $self->prepare('INSERT INTO evidence_coord (evidence_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand) VALUES ( ?, (SELECT sr.coord_system_id FROM seq_region sr, seq_region_attrib sra LEFT JOIN attrib_type at ON at.attrib_type_id = sra.attrib_type_id WHERE sr.seq_region_id = sra.seq_region_id AND at.code = "toplevel" AND sr.name = ?), ?, ?, ? )');
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

=head2 fetch_exact_evidence

 Arg [1]    : $id, integer
 Example    : $evidence_adaptor->fetch_by_dbID($id);
 Description: Fetch the evidence by dbID, unique element
 Returntype : Bio::EnsEMBL::EvidenceTracking::Evidence
 Exceptions : 


=cut

sub is_evidence_exists {
  my $self = shift;
  my $evidence = shift;
  my $query;
  if ($evidence->is_aligned eq 'y') {
    $query = 'SELECT e.evidence_id FROM evidence e, input_seq i, evidence_coord ec LEFT JOIN seq_region sr ON sr.seq_region_id = ec.seq_region_id'.
        ' WHERE ec.evidence_id = e.evidence_id AND ec.seq_region_start = '.$evidence->seq_region_start.
        ' AND ec.seq_region_end = '.$evidence->seq_region_end.
        ' AND ec.seq_region_strand = '.$evidence->seq_region_strand.
        ' AND sr.name = "'.$evidence->seq_region_name.
        '" AND e.input_seq_id = i.input_seq_id AND i.hit_name = "'.$evidence->input_seq->hit_name.'"';
  }
  else {
    $query = 'SELECT e.evidence_id FROM evidence e LEFT JOIN input_seq i ON i.input_seq_id = e.input_seq_id'.
        ' WHERE e.is_aligned = "'.$evidence->is_aligned.'"'.
        ' AND i.hit_name = "'.$evidence->input_seq->hit_name.'"';
  }
  my $sth = $self->dbc->prepare($query);
  $sth->execute();
  my $dbID = $sth->fetch();
  return @{$dbID}[0] if (defined $dbID);
  return 0;
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

  my $constraint = 'i.hit_name = "'.$hit_name.'"';
  return $self->generic_fetch($constraint);
}

=head2 fetch_all_unmapped_evidence

 Example    : $evidence_adaptor->fetch_all_unmapped_evidence;
 Description: Fetch all the evidence that did not aligned on the genome,
              regardless of the reason or the analysis
 Returntype : listref of Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence object
 Exceptions : 


=cut

sub fetch_all_unmapped_evidence {
    my $self = shift;

    my $constraint = 'e.is_aligned = "n"';
    return $self->generic_fetch($constraint);
}

=head2 fetch_all_unmapped_evidence_by_analysis

 Arg [1]    : $analysis, an Bio::EnsEMBL::Analysis object
 Example    : $evidence_adaptor->fetch_all_unmapped_evidence_by_analysis($analysis);
 Description: Fetch all evidences that did not align on the genome for an analysis
 Returntype : listref of Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence object
 Exceptions : throw if $analysis is not an Bio::EnsEMBL::Analysis object


=cut

sub fetch_all_unmapped_evidence_by_analysis {
    my $self = shift;
    my $analysis = shift;

    throw('Need a Bio::EnsEMBL::Analysis object for method fetch_all_unmapped_evidence_by_analysis, not a '.ref($analysis))
        unless $analysis->isa('Bio::EnsEMBL::Analysis');
    sub _tables {
      return ( ['evidence' , 'e'], ['input_seq', 'i'], ['evidence_coord', 'ec'], ['seq_region', 'sr'], ['analysis_run', 'ar'], ['track_evidence', 'te'] );
    }
    sub _left_join {
        return (['evidence_coord', 'e.evidence_id = ec.evidence_id'],
                ['input_seq', 'e.input_seq_id = i.input_seq_id'],
                ['seq_region', 'sr.seq_region_id = ec.seq_region_id'],
                ['track_evidence', 'e.evidence_id = te.evidence_id'],
                ['analysis_run', 'te.analysis_run_id = ar.analysis_run_id']);
    }
    my $constraint = 'e.is_aligned = "n" AND ar.analysis_id = '.$analysis->dbID;
    return $self->generic_fetch($constraint);
}

=head2 fetch_all_unmapped_evidence_by_reason

 Arg [1]    : $reason, an Bio::EnsEMBL::Analysis::EvidenceTracking::Reason object
 Example    : $evidence_adaptor->fetch_all_unmapped_evidence_by_reason($reason);
 Description: Fetch all evidences that did not align on the genome because of a reason
 Returntype : listref of Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence object
 Exceptions : throw if $reason is not an Bio::EnsEMBL::Analysis::EvidenceTracking::Reason object


=cut

sub fetch_all_unmapped_evidence_by_reason {
    my $self = shift;
    my $reason = shift;

    throw('Need a Bio::EnsEMBL::Analysis::EvidenceTracking::Reason object for method fetch_all_unmapped_evidence_by_reason, not a '.ref($reason))
        unless $reason->isa('Bio::EnsEMBL::Analysis::EvidenceTracking::Reason');
    sub _tables {
      return ( ['evidence' , 'e'], ['input_seq', 'i'], ['evidence_coord', 'ec'], ['seq_region', 'sr'], ['track_evidence', 'te'] );
    }
    sub _left_join {
        return (['evidence_coord', 'e.evidence_id = ec.evidence_id'],
                ['input_seq', 'e.input_seq_id = i.input_seq_id'],
                ['seq_region', 'sr.seq_region_id = ec.seq_region_id'],
                ['track_evidence', 'e.evidence_id = te.evidence_id']);
    }
    my $constraint = 'e.is_aligned = "n" AND te.reason_id = '.$reason->dbID;
    return $self->generic_fetch($constraint);
}

=head2 fetch_all_unmapped_evidence_by_analysis_and_reason

 Arg [1]    : $analysis, an Bio::EnsEMBL::Analysis object
 Arg [2]    : $reason, an Bio::EnsEMBL::Analysis::EvidenceTracking::Reason object
 Example    : $evidence_adaptor->fetch_all_unmapped_evidence_by_analysis_and_reason($analysis, $reason);
 Description: Fetch all evidence that did not align on the genome for an analysis and because of
              a reason
 Returntype : listref of Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence object
 Exceptions : throw if $reason is not an Bio::EnsEMBL::Analysis::EvidenceTracking::Reason object
              or if $analysis is not an Bio::EnsEMBL::Analysis object


=cut

sub fetch_all_unmapped_evidence_by_analysis_and_reason {
    my $self = shift;
    my ($analysis, $reason) = shift;

    throw('Need a Bio::EnsEMBL::Analysis::EvidenceTracking::Reason object for method fetch_all_unmapped_evidence_by_reason, not a '.ref($reason))
        unless $reason->isa('Bio::EnsEMBL::Analysis::EvidenceTracking::Reason');
    throw('Need a Bio::EnsEMBL::Analysis object for method fetch_all_unmapped_evidence_by_analysis, not a '.ref($analysis))
        unless $analysis->isa('Bio::EnsEMBL::Analysis');
    sub _tables {
      return ( ['evidence' , 'e'], ['input_seq', 'i'], ['evidence_coord', 'ec'], ['seq_region', 'sr'], ['track_evidence', 'te'] );
    }
    sub _left_join {
        return (['evidence_coord', 'e.evidence_id = ec.evidence_id'],
                ['input_seq', 'e.input_seq_id = i.input_seq_id'],
                ['seq_region', 'sr.seq_region_id = ec.seq_region_id'],
                ['track_evidence', 'e.evidence_id = te.evidence_id']);
    }
    my $constraint = 'e.is_aligned = "n" AND e.evidence_id = te.evidence_id AND te.reason_id = '.$reason->dbID
                    .' AND te.analysis_run_id = ar.analysis_run_id AND ar.analysis_id = '.$analysis->dbID;
    return $self->generic_fetch($constraint);
}

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
    my $analysis = shift;

    throw('Need a Bio::EnsEMBL::Analysis object for method fetch_all_unmapped_evidence_by_analysis, not a '.ref($analysis))
        unless $analysis->isa('Bio::EnsEMBL::Analysis');
    sub _tables {
      return ( ['evidence' , 'e'], ['input_seq', 'i'], ['evidence_coord', 'ec'], ['seq_region', 'sr'], ['analysis_run', 'ar'], ['track_evidence', 'te'] );
    }
    sub _left_join {
        return (['evidence_coord', 'e.evidence_id = ec.evidence_id'],
                ['input_seq', 'e.input_seq_id = i.input_seq_id'],
                ['seq_region', 'sr.seq_region_id = ec.seq_region_id'],
                ['track_evidence', 'e.evidence_id = te.evidence_id'],
                ['analysis_run', 'te.analysis_run_id = ar.analysis_run_id']);
    }
    my $constraint = 'e.evidence_id = te.evidence_id AND te.analysis_run_id = ar.analysis_run_id AND ar.analysis_id = '.$analysis->dbID;
    return $self->generic_fetch($constraint);
}

=head2 fetch_by_name

 Arg [1]    : $hit_name, sequence name as know in the database
 Example    : $evidence_adaptor->fetch_by_name($hit_name);
 Description: Fetch evidences by hit name, can return multiple evidence
 Returntype : listref of Bio::EnsEMBL::EvidenceTracking::Evidence
 Exceptions : 


=cut

sub fetch_by_name {
    my $self = shift;
    my $name = shift;

    my $constraint = 'i.hit_name = '.$name;
    return $self->generic_fetch($constraint);
}

=head2 fetch_all_after_analysis_run

 Arg [1]    : $analysis_run_id
 Example    : $evidence_adaptor->fetch_all_after_analysis_run($analysis_run_id);
 Description: Return a list of evidence which are kept for the next analsyis
 Returntype : listref of Bio::EnsEMBL::EvidenceTracking::Evidence
 Exceptions : 


=cut

sub fetch_all_after_analysis_run {
    my $self = shift;
    my $analysis = shift;

    throw('Need a Bio::EnsEMBL::Analysis object for method fetch_all_unmapped_evidence_by_analysis, not a '.ref($analysis))
        unless $analysis->isa('Bio::EnsEMBL::Analysis');
    sub _tables {
      return ( ['evidence' , 'e'],
               ['input_seq', 'i'],
               ['evidence_coord', 'ec'],
               ['seq_region', 'sr'],
               ['analysis_run', 'ar'],
               ['track_evidence', 'te'] );
    }
    sub _left_join {
        return (['evidence_coord', 'e.evidence_id = ec.evidence_id'],
                ['input_seq', 'e.input_seq_id = i.input_seq_id'],
                ['seq_region', 'sr.seq_region_id = ec.seq_region_id'],
                ['track_evidence', 'e.evidence_id = te.evidence_id'],
                ['analysis_run', 'te.analysis_run_id = ar.analysis_run_id']);
    }
    my $constraint = 'ar.analysis_id = '.$analysis->dbID.' AND te.reason_id < 100';
    return $self->generic_fetch($constraint);
}


#@@@@@@@
# Done @
#@@@@@@@

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
  return ( ['evidence' , 'e'], ['input_seq', 'i'], ['evidence_coord', 'ec'], ['seq_region', 'sr'] );
}

=head2 _columns

 Example    : $self->_columns;
 Description: Return a list of columns
 Returntype : a listref of string
 Exceptions : 


=cut

sub _columns {
  return ( 'e.evidence_id', 'i.input_seq_id', 'e.is_aligned', 'sr.name', 'ec.seq_region_start', 'ec.seq_region_end', 'ec.seq_region_strand' );
}

=head2 

 Example    : $self->_left_join;
 Description: Return a list of left joins
 Returntype : a listref of string
 Exceptions : 


=cut

sub _left_join {
    return (['evidence_coord', 'e.evidence_id = ec.evidence_id'], ['input_seq', 'e.input_seq_id = i.input_seq_id'], ['seq_region', 'sr.seq_region_id = ec.seq_region_id']);
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
  my ( $evidence_id, $input_seq_id, $is_aligned, $seq_region_name, $seq_region_start, $seq_region_end, $seq_region_strand);
  $sth->bind_columns( \$evidence_id, \$input_seq_id, \$is_aligned, \$seq_region_name, \$seq_region_start, \$seq_region_end, \$seq_region_strand);

  while($sth->fetch()) {
    my $input_seq = $self->db->get_InputSeqAdaptor->fetch_by_dbID($input_seq_id);
    push @out, Bio::EnsEMBL::Analysis::EvidenceTracking::Evidence->new(
              -dbID       => $evidence_id,
              -adaptor    => $self,
              -input_seq  => $input_seq,
              -is_aligned => $is_aligned
              );
    if ($out[-1]->is_aligned eq 'y') {
        $out[-1]->seq_region_name($seq_region_name);
        $out[-1]->seq_region_start($seq_region_start);
        $out[-1]->seq_region_end($seq_region_end);
        $out[-1]->seq_region_strand($seq_region_strand);
    }
  }
  return \@out;
}
1;
