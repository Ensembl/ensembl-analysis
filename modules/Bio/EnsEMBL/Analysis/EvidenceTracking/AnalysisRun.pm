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

Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun - Object with all the information about the run

=head1 SYNOPSIS

  use Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun;

  my $analysis_run = Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun->new(
     -analysis_id  => $analysis->dbID,
     -input_db_id  => $in_databases,
     -output_db_id => $out_databases
    );

=head1 DESCRIPTION

  This module stores the current of the analysis, its goal is to provide a way
  to differenciate between several runs of one analysis.

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Storable;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

@ISA = qw(Bio::EnsEMBL::Storable);


=head2 new

 Arg [1]    : $dbid, int
 Arg [2]    : $analysis_id, int
 Arg [3]    : $input_db_id, string. An integer or a string of number separated with :, 2:4
 Arg [4]    : $output_db_id, string. An integer or a string of number separated with :, 5
 Arg [5]    : $run_date, when the object is stored it's set with the MySQL now()
 Arg [6]    : $adaptor, Bio::EnsEMBL::Analysis::DBSQL::AnalysisRunAdaptor object
 Example    : $analysis_run = Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun->new(
              -analysis_id  => $analysis->dbID,
              -input_db_id  => $in_databases,
              -output_db_id => $out_databases
            );
 Description: Constructor
 Returntype : Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun
 Exceptions : 


=cut

sub new {
  my($class,@args) = @_;

  my $self = bless {},$class;

  my ($id, $analysis_id, $input_db_id, $output_db_id, $run_date, $adaptor) =
          rearrange([qw(DBID
                        ANALYSIS_ID
                        INPUT_DB_ID
                        OUTPUT_DB_ID
                        RUN_DATE
                        ADAPTOR
                        )],@args);

  $self->dbID      ( $id ) if (defined $id);
  $self->analysis_id ( $analysis_id );
  $self->input_db_id ( $input_db_id );
  $self->output_db_id ( $output_db_id );
  $self->run_date ( $run_date );
  $self->adaptor   ( $adaptor );
  return $self; # success - we hope!
}

=head2 analysis_id

 Arg [1]    : $analysis_id, int [optional]
 Example    : $analysis_run->analysis_id($analysis_id);
 Description: Getter/Setter for the analysis id
 Returntype : integer, the analysis id
 Exceptions : 


=cut

sub analysis_id {
  my $self = shift;
  $self->{'analysis_id'} = shift if ( @_ );
  return $self->{'analysis_id'};
}

=head2 input_db_id

 Arg [1]    : $input_db_id, string [optional]
 Example    : $analysis_run->input_db_id($input_db_id);
 Description: Getter/Setter for the databases used as input
              It's either a number or a strin of number like 1:2
 Returntype : string, the input db ids
 Exceptions : 


=cut

sub input_db_id {
  my $self = shift;
  $self->{'input_db_id'} = shift if ( @_ );
  return $self->{'input_db_id'};
}

=head2 output_db_id

 Arg [1]    : $output_db_id, string [optional]
 Example    : $analysis_run->output_db_id($output_db_id);
 Description: Getter/Setter for the databases used as output
              It's either a number or a strin of number like 1:2
 Returntype : string, the output db ids
 Exceptions : 


=cut

sub output_db_id {
  my $self = shift;
  $self->{'output_db_id'} = shift if ( @_ );
  return $self->{'output_db_id'};
}

=head2 run_date

 Arg [1]    : $run_date, int [optional]
 Example    : $analysis_run->run_date($run_date);
 Description: Getter/Setter for the analysis id
 Returntype : integer, the analysis id
 Exceptions : 


=cut

sub run_date {
  my $self = shift;
  $self->{'run_date'} = shift if ( @_ );
  return $self->{'run_date'};
}

=head2 is_stored

 Arg [1]    : $analysisrun_adaptor, a Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::AnalysisRunAdaptor object
 Example    : $analysis_run->is_stored($analysisrun_adaptor);
 Description: Test if the object is alreday stored
 Returntype : boolean
 Exceptions : 


=cut

sub is_stored {
  my $self = shift;
  my $db = shift;

  if($db and $db->isa('Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::AnalysisRunAdaptor')) {
    my $current_analysisrun = $db->get_current_analysis_run_by_analysis_id($self->analysis_id);
    return 0 unless (defined $current_analysisrun);
    if ($current_analysisrun->run_date eq $self->run_date) {
        return 1;
    }
    return 0;
  }
  else {
    throw('db argument must be a Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::AnalysisRunAdaptor not '.ref($db));
  }
}

1;
