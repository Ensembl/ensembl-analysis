#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;
use Bio::EnsEMBL::Utils::ConfigRegistry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor;

my $curr_release = $ENV{ENSEMBL_RELEASE};
#my $curr_release = 97;
my $prev_release = $curr_release - 1;

# ---------------------- CURRENT COMPARA DATABASE ---------------------------------

my $compara_dbs = {
#    'compara_curr'   => [ 'mysql-ens-compara-prod-1', 'ensembl_compara_'.$curr_release ],
    'compara_curr'   => [ 'mysql-ens-mirror-1', 'ensembl_compara_'.$curr_release ],
};

foreach my $alias_name ( keys %$compara_dbs ) {
  my ( $host, $db_name ) = @{ $compara_dbs->{$alias_name} };
  my ( $user, $pass ) = ( 'ensadmin', $ENV{'ENSADMIN_PSW'} );
  ( $user, $pass ) = ( 'ensro', '' ) if ( $alias_name  =~ /_prev/ );
  Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
            -host => $host,
            -user => $user,
            -pass => $pass,
            -port => get_port($host),
            -species => $alias_name,
            -dbname  => $db_name,
        );
}

# ---------------------- CURRENT CORE DATABASES ---------------------------------

# The majority of core databases live on staging servers:
Bio::EnsEMBL::Registry->load_registry_from_url(
   "mysql://ensro\@mysql-ens-sta-1.ebi.ac.uk:4519/$curr_release");

sub get_port {
    my $host = shift;
    my $port = `$host port`;
    chomp $port;
    return $port;
  }

1;
