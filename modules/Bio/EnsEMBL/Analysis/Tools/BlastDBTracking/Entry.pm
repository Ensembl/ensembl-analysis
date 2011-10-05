package Bio::EnsEMBL::Analysis::Tools::BlastDBTracking::Entry;

use namespace::autoclean;
use Moose;

has filename       => ( is => 'ro', isa => 'Str', required => 1 );
has version        => ( is => 'ro', isa => 'Str', required => 1 );
has sanger_version => ( is => 'ro', isa => 'Int',               );
has installation   => ( is => 'ro', isa => 'Int',               );
has count          => ( is => 'ro', isa => 'Int',               );
has checksum       => ( is => 'ro', isa => 'Str',               );
has from_file      => ( is => 'ro', isa => 'Bool',              );

1;
