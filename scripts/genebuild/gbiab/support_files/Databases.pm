package Reg;
use warnings;
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Production::DBSQL::DBAdaptor;
use Bio::EnsEMBL::MetaData::DBSQL::MetaDataDBAdaptor;
use Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor;
{
Bio::EnsEMBL::Production::DBSQL::DBAdaptor->new(
-host => 'mysql-ens-meta-prod-1',
-port => '4483',
-dbname => 'ensembl_production',
-user => 'ensro',
-species => 'multi',
-group => 'production',
);
Bio::EnsEMBL::MetaData::DBSQL::MetaDataDBAdaptor->new(
-host => 'mysql-ens-meta-prod-1',
-port => '4483',
-dbname => 'ensembl_metadata',
-user => 'ensro',
-species => 'multi',
-group => 'metadata',
);
Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor->new(
-host => 'mysql-ens-meta-prod-1',
-port => '4483',
-dbname => 'ncbi_taxonomy',
-user => 'ensro',
-species => 'multi',
-group => 'taxonomy',
);
}
