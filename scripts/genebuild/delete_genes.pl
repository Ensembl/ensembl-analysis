#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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
  -update_stable_event Add the relevant rows into the mapping_session and
              stable_id_event tables. It assumes that the name of the database
              follows the Ensembl convention
  -use_last_session This is usefull to delete a gene after a stable id mapping
              has been run without having to rerun the mapping again.
              Not recommended...
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
my $update_stable_event = 0;
my $use_last_session = 0;
my $config_dbname;


GetOptions( 'dbhost|host|h:s'        => \$host,
            'dbport|port|P:n'        => \$port,
            'dbname|db|D:s'        => \$dbname,
            'dbuser|user|u:s'        => \$user,
            'dbpass|pass|p:s'        => \$pass,
            'idfile:s'        => \$idfile,
            'stable_id!'      => \$stable_id,
            'all!'      => \$all,
            'update_stable_event!' => \$update_stable_event,
            'use_last_session!' => \$use_last_session,
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
my $session_id;
if ($update_stable_event) {
  $session_id = get_mapping_session_id($db, $use_last_session);
}

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
      if ($update_stable_event) {
        update_stable_event_data($gene, $session_id);
      }
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


=head2 get_mapping_session_id

 Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor
 Arg [2]    : Boolean
 Description: Fetch the last session id from mapping_session.
              If Arg[2] is true, simply use the latest session.
              It will be created if it doesn't exist.
              Otherwise create a new mapping session
 Returntype : Int
 Exceptions : Throws if the session id cannot be retrieved

=cut

sub get_mapping_session_id {
  my ($db, $use_last_session) = @_;

  my $session_id;
  my $dbc = $db->dbc;
  my $select_sth = $dbc->prepare('SELECT * FROM mapping_session ORDER BY mapping_session_id DESC LIMIT 1');
  my @dbname = split('_', $dbc->dbname);
  $dbname[-2]--; # we work on the assumption that it is a official Ensembl DB
  my $last_dbname = join('_', @dbname);
  if ($use_last_session) {
    $select_sth->execute;
    foreach my $row (@{$select_sth->fetchall_arrayref}) {
      $session_id = $row->[0];
    }
  }
  if (!$session_id) {
    my $cs_version = $db->get_CoordSystemAdaptor->get_default_version;
    my $insert_sth = $dbc->prepare(
      'INSERT INTO mapping_session'.
      ' (old_db_name, new_db_name, old_release, new_release, old_assembly, new_assembly, created)'.
      ' VALUES (?, ?, ?, ?, ?, ?, NOW())'
    );
    $insert_sth->bind_param(1, join('_', @dbname));
    $insert_sth->bind_param(2, $dbc->dbname);
    $insert_sth->bind_param(3, $dbname[-2]);
    $insert_sth->bind_param(4, $dbname[-2]+1);
    $insert_sth->bind_param(5, $cs_version);
    $insert_sth->bind_param(6, $cs_version);
    $insert_sth->execute;
    $select_sth->execute;
    foreach my $row (@{$select_sth->fetchall_arrayref}) {
      $session_id = $row->[0];
    }
  }
  throw("'Could not retrieve a session id: '$session_id'") unless ($session_id and $session_id > 0);
  return $session_id;
}


=head2 update_stable_event_data

 Arg [1]    : Bio::EnsEMBL::Gene, it must have an adaptor
 Arg [2]    : Int Session id
 Description: Add the necesary rows to the stable_id_event
              table for all translations, transcripts of the
              gene and the gene itself.
 Returntype : None
 Exceptions : None

=cut

sub update_stable_event_data {
  my ($gene, $session_id) = @_;

  my $set_sth = $gene->adaptor->db->dbc->prepare(
    'INSERT INTO stable_id_event'.
    ' (old_stable_id, old_version, mapping_session_id, type)'.
    " VALUES (?, ?, $session_id, ?)"
  );
  foreach my $transcript (@{$gene->get_all_Transcripts}) {
    if ($transcript->translation) {
      $set_sth->bind_param(1, $transcript->translation->stable_id);
      $set_sth->bind_param(2, $transcript->translation->version);
      $set_sth->bind_param(3, 'translation');
      $set_sth->execute;
    }
    $set_sth->bind_param(1, $transcript->stable_id);
    $set_sth->bind_param(2, $transcript->version);
    $set_sth->bind_param(3, 'transcript');
    $set_sth->execute;
  }
  $set_sth->bind_param(1, $gene->stable_id);
  $set_sth->bind_param(2, $gene->version);
  $set_sth->bind_param(3, 'gene');
  $set_sth->execute;
}
