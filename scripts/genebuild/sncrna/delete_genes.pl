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

=head1 NAME

  delete_genes.pl

=head1 DESCRIPTION

  Given a list of gene_ids or stable_ids, deletes the genes from the 
  specified database. If config_dbname is provided, the script reads
  database details from the Bio::EnsEMBL::Analysis::Config::GeneBuild::Databases
  configuration.

=head1 OPTIONS

  -dbhost     host name for database (gets put as host= in locator)
  -dbport     For RDBs, what port to connect to (port= in locator)
  -dbname     For RDBs, what name to connect to (dbname= in locator)
  -dbuser     For RDBs, what username to connect as (dbuser= in locator)
  -dbpass     For RDBs, what password to use (dbpass= in locator)
  -idfile     File with internal gene ids or stable IDs to be deleted
  -stable_id  A boolean flag to indicate that the file specified in
              -idfile contains gene stable IDs.
  -all        A boolean flag to delete all genes in the db, no file needed
=head1 EXAMPLES

  perl delete_genes.pl -dbhost my_host -dbuser ensadmin -dbpass **** \
    -dbname rat_Jun03_mk2 -idfile genes_to_delete.dbIDs

  or

  perl delete_genes.pl -dbhost my_host -dbuser ensadmin -dbpass **** \
    -dbname some_database -idfile my_ENSG_IDs.txt -stable_id

  or

  perl delete_genes.pl -config_dbname GENEBUILD_DB -idfile genes_to_delete

=cut

use warnings ;
use strict;
use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::Utilities;
use Bio::EnsEMBL::Utils::Exception;

my $host;
my $port=3306;
my $dbname;
my $user;
my $pass;
my $idfile;
my $stable_id = 0;
my $all = 0;
my $config_dbname;


GetOptions( 'dbhost|host|h:s'        => \$host,
            'dbport|port|P:n'        => \$port,
            'dbname|db|D:s'        => \$dbname,
            'dbuser|user|u:s'        => \$user,
            'dbpass|pass|p:s'        => \$pass,
            'idfile:s'        => \$idfile,
            'stable_id!'      => \$stable_id,
            'all!'      => \$all,
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

if($idfile) {
open(INFILE, "<$idfile") or die ("Can't read $idfile $! \n");

  while (<INFILE>) {
    chomp;
    my $gene_id = $_;

    eval{
      my $gene;
      if ($stable_id) {
        $gene = $gene_adaptor->fetch_by_stable_id($gene_id);
      } else {
        $gene = $gene_adaptor->fetch_by_dbID($gene_id);
      }

      # it seems that some xrefs might not be deleted when using this method
      # Could it be because the gene is lazy-loaded?
      $gene_adaptor->remove($gene);
      print STDERR "Deleted $gene_id\n";
    };

    if($@){
      print "Couldn't remove gene $gene_id ($@)\n";
    }
  }
  close(INFILE);
}

elsif($all) {
  my $genes = $gene_adaptor->fetch_all();
  foreach my $gene (@{$genes}) {
    eval{
      $gene_adaptor->remove($gene);
      print STDERR "Deleted gene ".$gene->dbID();
    };
    if($@){
      print "Couldn't remove gene ".$gene->dbID()." (".$@.")\n";
    }
  }
}

else {
  die "No delete option was selected, either use -idfile or -all";
}
