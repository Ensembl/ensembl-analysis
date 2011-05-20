package Bio::EnsEMBL::Analysis::EvidenceTracking::Reason;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Storable;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Analysis::EvidenceTracking::DBSQL::ReasonAdaptor;

@ISA = qw(Bio::EnsEMBL::Storable);


sub new {
  my($class,@args) = @_;

  my $self = bless {},$class;

  my ($id, $info, $adaptor) = rearrange([qw(DBID INFO ADAPTOR)],@args);

  $self->dbID    ( $id );
  $self->info  ( $info );
  return $self; # success - we hope!
}


sub info {
  my $self = shift;
  $self->{'info'} = shift if ( @_ );
  return $self->{'info'};
}

1;
