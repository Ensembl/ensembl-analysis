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
