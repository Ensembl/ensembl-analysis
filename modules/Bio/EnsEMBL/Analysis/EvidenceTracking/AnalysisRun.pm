package Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Storable;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

@ISA = qw(Bio::EnsEMBL::Storable);


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


sub analysis_id {
  my $self = shift;
  $self->{'analysis_id'} = shift if ( @_ );
  return $self->{'analysis_id'};
}

sub input_db_id {
  my $self = shift;
  $self->{'input_db_id'} = shift if ( @_ );
  return $self->{'input_db_id'};
}

sub output_db_id {
  my $self = shift;
  $self->{'output_db_id'} = shift if ( @_ );
  return $self->{'output_db_id'};
}

sub run_date {
  my $self = shift;
  $self->{'run_date'} = shift if ( @_ );
  return $self->{'run_date'};
}

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
