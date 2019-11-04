use strict;
use warnings;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::OntologyDBAdaptor;
use Bio::EnsEMBL::Production::DBSQL::DBAdaptor;

{

  Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host    => '',
    -port    => ,
    -user    => 'ensro',
    -dbname  => '',
    -species => '',
    -group   => 'core',
  );

  my $db_suffix = $ENV{'ENSEMBL_RELEASE'} - 1;

  Bio::EnsEMBL::DBSQL::OntologyDBAdaptor->new(
      -HOST   => 'mysql-ens-sta-1',
      -PORT   => 4519,
      -USER   => 'ensro',
      -DBNAME => 'ensembl_ontology_'.$db_suffix,
      -SPECIES => 'multi',
      -GROUP   => 'ontology'
  );

  Bio::EnsEMBL::Production::DBSQL::DBAdaptor->new(
      '-species' => 'multi',
      '-group'   => 'production',
      '-host'    => 'mysql-eg-pan-prod.ebi.ac.uk',
      '-port'    => 4276,
      '-user'    => 'ensro',
      '-dbname'  => 'ensembl_production'
  );

}
1;
