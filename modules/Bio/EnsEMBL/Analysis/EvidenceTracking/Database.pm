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

Bio::EnsEMBL::Analysis::EvidenceTracking::Database - Keep information about used databases

=head1 SYNOPSIS

  use Bio::EnsEMBL::Analysis::EvidenceTracking::Database;

  my $output_db = Bio::EnsEMBL::Analysis::EvidenceTracking::Database->new(
    -dbname => $db_name,
    -instance => $instance
    );

=head1 DESCRIPTION

  Get all the information about the databases used to store the analysis.
  It allows us to know which one where use for each run

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::EvidenceTracking::Database;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Storable;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

@ISA = qw(Bio::EnsEMBL::Storable);

=head2 new

 Arg [1]    : $dbid, int
 Arg [2]    : $db_name, string, name of the database
 Arg [3]    : $instance, string, name of the instance
 Arg [4]    : $adaptor, Bio::EnsEMBL::Analysis::DBSQL::AnalysisRunAdaptor object
 Example    : $analysis_run = Bio::EnsEMBL::Analysis::EvidenceTracking::AnalysisRun->new(
              -db_name  => $db_name,
              -instance => $instance
            );
 Description: Constructor
 Returntype : Bio::EnsEMBL::Analysis::EvidenceTracking::Database
 Exceptions : 


=cut

sub new {
  my($class,@args) = @_;

  my $self = bless {},$class;

  my ($id, $db_name, $instance, $adaptor) =
          rearrange([qw(DBID
                        DB_NAME
                        INSTANCE
                        ADAPTOR
                        )],@args);

  $self->dbID      ( $id ) if (defined $id);
  $self->db_name ( $db_name );
  $self->instance ( $instance );
  $self->adaptor   ( $adaptor );
  return $self; # success - we hope!
}


=head2 db_name

 Arg [1]    : $db_name, int [optional]
 Example    : $database->db_name($db_name);
 Description: Getter/Setter for the database name
 Returntype : integer, the database name
 Exceptions : 


=cut

sub db_name {
  my $self = shift;
  $self->{'db_name'} = shift if ( @_ );
  return $self->{'db_name'};
}

=head2 instance

 Arg [1]    : $instance, int [optional]
 Example    : $database->instance($instance);
 Description: Getter/Setter for the Mysql instance name
 Returntype : integer, the instance name
 Exceptions : 


=cut

sub instance {
  my $self = shift;
  $self->{'instance'} = shift if ( @_ );
  return $self->{'instance'};
}

=head2 is_stored

 Arg [1]    : $database_adaptor, a Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::Database object
 Example    : $database->is_stored($database_adaptor);
 Description: Test if the object is alreday stored
 Returntype : boolean
 Exceptions : 


=cut

sub is_stored {
  my $self = shift;
  my $db = shift;

  if($db and $db->isa('Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor')) {
    my $database = $db->get_DatabaseAdaptor->fetch_DB($self->db_name, $self->instance);
    if ($database) {
        $self->dbID($database->dbID);
        return $self->dbID;
    }
    return 0;
  }
  else {
    throw('db argument must be a Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::DBAdaptor');
  }
}

1;
