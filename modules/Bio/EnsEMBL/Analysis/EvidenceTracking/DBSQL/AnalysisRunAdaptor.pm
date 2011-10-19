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

Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::AnalysisRunAdaptor - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::AnalysisRunAdaptor; 

use strict;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( deprecate throw warning stack_trace_dump );
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun;
use Data::Dumper;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);


=head2 store

 Arg [1]    : $analysis_run, Bio::EnsEMBL::Pipeline::AnalysisRun
 Example    : $analysisrun_adaptor->store($analysis_run);
 Description: Store the new run of the analysis
 Returntype : an integer, the dbID
 Exceptions : if not successful


=cut

sub store {
  my ($self, $analysis_run) = @_;

  if (!ref $analysis_run || !$analysis_run->isa('Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun') ) {
    throw("Must store a AnalysisRun object, not a $analysis_run");
  }

  my $db = $self->db();
  my @a_input_dbs = @{$analysis_run->input_dbs};
  my @a_output_dbs = @{$analysis_run->output_dbs};
  my $input_db_id = shift @a_input_dbs;
  my $output_db_id = shift @a_output_dbs;
#  print STDERR Dumper($analysis_run);
  while (my $db = shift @a_input_dbs) {
      $input_db_id .= ':'.$db->dbID;
  }
  while (my $db = shift @a_output_dbs) {
      $output_db_id .= ':'.$db->dbID;
  }

  # make an sql statement
  my $sth = $self->prepare("INSERT INTO analysis_run ( analysis_id,
                            run_date,input_db_id, output_db_id) 
                            VALUES ( ?,NOW(),?,?)");

  $sth->bind_param( 1, $analysis_run->analysis_id, SQL_INTEGER );
  $sth->bind_param( 2, $input_db_id, SQL_VARCHAR );
  $sth->bind_param( 3, $output_db_id, SQL_VARCHAR );

  $sth->execute();
  $sth->finish();
  my $analysis_run_dbID = $sth->{'mysql_insertid'};

  # set the adaptor and dbID on the original passed in analysis_run not the
  # transfered copy
  $analysis_run->adaptor($self);
  $analysis_run->dbID($analysis_run_dbID);
  print STDERR "Stored AnanlysisRun object ".$analysis_run->dbID."\n";
  return $analysis_run_dbID;
}

=head2 update

 Arg [1]    : a Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun object
 Example    : $analysisrun_adaptor->update($analysis_run);
 Description: Update the analysis_run table
 Returntype : void
 Exceptions : throw if not of the type Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun


=cut

sub update {
  my ($self, $analysis_run) = @_;

  if (!ref $analysis_run || !$analysis_run->isa('Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun') ) {
    throw("Must store a AnalysisRun object, not a $analysis_run");
  }

  my $db = $self->db();
  my $input_db_id;
  my $output_db_id;
  for my $tuple ([$analysis_run->input_dbs, $input_db_id], [$analysis_run->output_dbs, $output_db_id]) {
      my @tmp = @{$tuple->[0]};
      my $tuple->[1] = shift(@tmp);
      while (my $db = shift(@tmp)) {
          $tuple->[1] .= ':'.$db->dbID;
      }
  }

  # make an sql statement
  my $sth = $self->prepare("UPDATE analysis_run SET analysis_id = ?,
                            input_db_id = ?, output_db_id = ?) 
                            WHERE analysis_run_id = ?");

  $sth->bind_param( 1, $analysis_run->analysis_id, SQL_INTEGER );
  $sth->bind_param( 2, $analysis_run->input_db_id, SQL_VARCHAR );
  $sth->bind_param( 3, $analysis_run->output_db_id, SQL_VARCHAR );
  $sth->bind_param( 4, $analysis_run->dbID, SQL_INTEGER );

  $sth->execute();
  $sth->finish();
}

=head2 fetch_by_dbID

 Arg [1]    : $id, integer
 Example    : $analysisrun_adaptor->fetch_by_dbID($id);
 Description: Fetch the evidence by dbID, unique element
 Returntype : Bio::EnsEMBL::Pipeline::AnalysisRun
 Exceptions : 


=cut

sub fetch_by_dbID {
  my $self = shift;
  my $analysis_run_id = shift;
  my $constraint = "ar.analysis_run_id = '$analysis_run_id'";
  my ($analysis_run) = @{ $self->generic_fetch($constraint) };
  return $analysis_run;
}

=head2 fetch_by_analysis

 Arg [1]    : $analysis_id
 Example    : $analysisrun_adaptor->fetch_by_analysis($analysis_id);
 Description: Return all run for the analysis $analysis_id
 Returntype : listref of Bio::EnsEMBL::Pipeline::AnalysisRun
 Exceptions : 


=cut

sub fetch_by_analysis {
  my $self = shift;
  my $analysis_id = shift;
  my $constraint = 'ar.analysis_id = '.$analysis_id;
  return $self->generic_fetch($constraint);
}

=head2 get_current_analysis_run

 Example    : $analysisrun_adaptor->get_current_analysis_run;
 Description: Return the current analysis that has been run, the last one
 Returntype : Bio::EnsEMBL::Pipeline::AnalysisRun object
 Exceptions : 


=cut

sub get_current_analysis_run {
    my $self = shift;

    my $constraint = 'ar.analysis_run_id IS NOT null ORDER BY analysis_run_id DESC LIMIT 1';
    my ($current_analysis) = $self->generic_fetch($constraint);
    return $current_analysis->[0];
}

=head2 get_current_analysis_run_by_analysis_id

 Arg [1]    : $analysis_id, int
 Example    : $analysisrun_adaptor->get_current_analysis_run_by_analysis_id($analysis_id);
 Description: Return the current analysis that has been run for a specific analysis $analysis_id
 Returntype : Bio::EnsEMBL::Pipeline::AnalysisRun object
 Exceptions : 


=cut

sub get_current_analysis_run_by_analysis_id {
    my $self = shift;
    my $analysis_id = shift;

    my $constraint = 'ar.analysis_id = '.$analysis_id.' ORDER BY analysis_run_id DESC LIMIT 1';
    my ($current_analysis) = $self->generic_fetch($constraint);
    return $current_analysis->[0];
}

=head2 fetch_all_by_evidence

 Arg [1]    : $evidence
 Example    : $analysisrun_adaptor->fetch_all_by_evidence($evidence);
 Description: Return a list of analysis run on this evidence
 Returntype : listref of Bio::EnsEMBL::Pipeline::AnalysisRun
 Exceptions : 


=cut

sub fetch_all_by_evidence {
  my $self = shift;
  my $evidence = shift;
  sub _tables {
    return (['analysis_run' , 'ar'], ['track_evidence', 'te']);
  }
  sub _left_join {
      return ( ['track_evidence', 'te.analysis_run_id = ar.analysis_run_id'] );
  }
  my $constraint = 'te.evidence_id = '.$evidence->dbID;
  return $self->generic_fetch($constraint);
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
  return (['analysis_run' , 'ar']);
}

=head2 _columns

 Example    : $self->_columns;
 Description: Return a list of columns
 Returntype : a listref of string
 Exceptions : 


=cut

sub _columns {
  return ( 'ar.analysis_run_id', 'ar.analysis_id', 'ar.run_date',
           'ar.input_db_id', 'ar.output_db_id');
}

=head2 _objs_from_sth

 Arg [1]    : $sth
 Example    : $self->_objs_from_sth($sth);
 Description: Put the result of the query in Bio::EnsEMBL::Pipeline::InputSeq objects
 Returntype : listref of Bio::EnsEMBL::Pipeline::AnalysisRun object
 Exceptions : 


=cut


sub _objs_from_sth {
  my ($self, $sth) = @_;

  my @out;
  my ( $analysis_run_id, $analysis_id, $run_date, $input_db_id, $output_db_id);
  $sth->bind_columns( \$analysis_run_id, \$analysis_id, \$run_date, \$input_db_id, \$output_db_id);

  my @a_db_in;
  my @a_db_out;
  while($sth->fetch()) {
    foreach my $db (split(':', $input_db_id)) {
        push(@a_db_in, $self->db->get_DatabaseAdaptor->fetch_by_dbID($db));
    }
    foreach my $db (split(':', $output_db_id)) {
        push(@a_db_out, $self->db->get_DatabaseAdaptor->fetch_by_dbID($db));
    }
    push @out, Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun->new(
              -dbID         => $analysis_run_id,
              -adaptor      => $self,
              -analysis_id  => $analysis_id,
              -run_date     => $run_date,
              -input_db_id  => \@a_db_in,
              -output_db_id => \@a_db_out
              );
  }
  return \@out;
}


1;
