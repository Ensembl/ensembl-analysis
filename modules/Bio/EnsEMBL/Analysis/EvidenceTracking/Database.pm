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

Bio::EnsEMBL::Analysis::EvidenceTracking::Database - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::EvidenceTracking::Database;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Storable;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

@ISA = qw(Bio::EnsEMBL::Storable);


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


sub db_name {
  my $self = shift;
  $self->{'db_name'} = shift if ( @_ );
  return $self->{'db_name'};
}

sub instance {
  my $self = shift;
  $self->{'instance'} = shift if ( @_ );
  return $self->{'instance'};
}

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
