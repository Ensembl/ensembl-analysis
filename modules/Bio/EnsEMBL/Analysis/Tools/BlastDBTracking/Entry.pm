# $Source: /tmp/ENSCOPY-ENSEMBL-ANALYSIS/modules/Bio/EnsEMBL/Analysis/Tools/BlastDBTracking/Entry.pm,v $
# $Revision: 1.4 $
package Bio::EnsEMBL::Analysis::Tools::BlastDBTracking::Entry;

use warnings ;
use strict ;
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
