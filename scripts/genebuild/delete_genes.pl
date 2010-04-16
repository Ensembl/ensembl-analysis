#!/usr/local/ensembl/bin/perl -w

=head1 NAME

  delete_genes.pl

=head1 DESCRIPTION

  Given a list of gene_ids, deletes the genes from the specified
  database. If config_dbname is provided reads the database details
  from the Bio::EnsEMBL::Analysis::Config::GeneBuild::Databases
  configuration.

=head1 OPTIONS

  -dbhost   host name for database (gets put as host= in locator)
  -dbport   For RDBs, what port to connect to (port= in locator)
  -dbname   For RDBs, what name to connect to (dbname= in locator)
  -dbuser   For RDBs, what username to connect as (dbuser= in locator)
  -dbpass   For RDBs, what password to use (dbpass= in locator)
  -idfile   File with internal gene ids to be deleted
  -help     Summary of options

=head1 EXAMPLES

  perl delete_genes.pl -dbhost ecs2b -dbuser ensadmin -dbpass **** \
    -dbname rat_Jun03_mk2 -idfile genes_to_delete

  or

  perl delete_genes.pl -config_dbname GENEBUILD_DB -idfile genes_to_delete

=cut

use strict;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::Utilities;
use Bio::EnsEMBL::Utils::Exception;

my $host;
my $port;
my $dbname;
my $user;
my $pass;
my $idfile;
my $config_dbname;


GetOptions( 'dbhost:s'        => \$host,
            'dbport:n'        => \$port,
            'dbname:s'        => \$dbname,
            'dbuser:s'        => \$user,
            'dbpass:s'        => \$pass,
            'idfile:s'        => \$idfile,
            'config_dbname:s' => \$config_dbname, );


my $db;

if ($config_dbname) {
  $db = get_db_adaptor_by_string($config_dbname);
} elsif ( $dbname && $host ) {
  $db =
    new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $host,
                                        -user   => $user,
                                        -port   => $port,
                                        -dbname => $dbname,
                                        -pass   => $pass, );
} else {
  throw(   "Need to pass either -dbhost $host and -dbname $dbname or "
         . "-config_dbname $config_dbname for the script to work" );
}

my $gene_adaptor = $db->get_GeneAdaptor;

open(INFILE, "<$idfile") or die ("Can't read $idfile $! \n");

while (<INFILE>) {
  chomp;
  my $gene_id = $_;

  eval{
    my $gene = $gene_adaptor->fetch_by_dbID($gene_id);

    # it seems that some xrefs might not be deleted when using this method
    # Coud it be because the gene is lazy-loaded?

    $gene_adaptor->remove($gene);
    print STDERR "Deleted $gene_id\n";
  };
  if($@){
    print "Couldn't remove gene $gene_id ($@)\n";
  }
}
close(INFILE);
