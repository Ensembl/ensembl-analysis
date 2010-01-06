#!/usr/local/ensembl/bin/perl -w

=head1 NAME

  delete_genes.pl

=head1 SYNOPSIS
 
  delete_genes.pl
  deletes genes from given database whose ids are passed in through STDIN

=head1 DESCRIPTION


=head1 OPTIONS

    -dbhost      host name for database (gets put as host= in locator)
    -dbport      For RDBs, what port to connect to (port= in locator)
    -dbname    For RDBs, what name to connect to (dbname= in locator)
    -dbuser    For RDBs, what username to connect as (dbuser= in locator)
    -dbpass    For RDBs, what password to use (dbpass= in locator)
    -help      summary of options


=head2 EXAMPLES

./delete_genes.pl -dbhost ecs2b -dbuser ensadmin -dbpass **** -dbname rat_Jun03_mk2 genes_to_delete

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
my $config_dbname;

&GetOptions( 
            'dbhost:s'      => \$host,
            'dbport:n'      => \$port,
            'dbname:s'    => \$dbname,
            'dbuser:s'    => \$user,
            'dbpass:s'      => \$pass,
            'config_dbname:s' => \$config_dbname,
           );


my $db;

if($config_dbname){
  $db = get_db_adaptor_by_string($config_dbname);
}elsif($dbname && $host){
  $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                            -host   => $host,
                                            -user   => $user,
                                            -port   => $port,
                                            -dbname => $dbname,
                                            -pass => $pass,
                                           );
}else{
  throw("Need to pass either -dbhost $host and -dbname $dbname or ".
        "-config_dbname $config_dbname for the script to work");
}

my $gene_adaptor = $db->get_GeneAdaptor;

while(<>){
  chomp;
  my $gene_id= $_;
  
  #my $sth = $db->prepare("delete from gene where gene_id = $gene_id");
  #$sth->execute;
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

