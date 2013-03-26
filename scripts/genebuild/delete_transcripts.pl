#!/usr/bin/env perl
# $Source: /tmp/ENSCOPY-ENSEMBL-ANALYSIS/scripts/genebuild/delete_transcripts.pl,v $
# $Revision: 1.6 $
=head1 NAME

  delete_transcripts.pl

=head1 DESCRIPTION

  deletes transcripts from given database whose ids are passed in through STDIN as a file

=head1 OPTIONS

=head1 

=head2 DB connection 

    -dbhost         host name for database (gets put as host= in locator)
    -dbport         For RDBs, what port to connect to (port= in locator)
    -dbname         For RDBs, what name to connect to (dbname= in locator)
    -dbuser         For RDBs, what username to connect as (dbuser= in locator)
    -dbpass         For RDBs, what password to use (dbpass= in locator)

=head2 DB connection (alternative method)

    -config_dbname  the alias for the DB you want to connect to, as defined in 
                    Bio::EnsEMBL::Analysis::Config::Databases config file.
                    Use this option or the dbhost/dbport/dbname etc options above,
                    but not both.

=head2 Other options

    -stable_id      A boolean flag to indicate that the IDs passed in are stable IDs

=head1 EXAMPLES

./delete_transcripts.pl -dbhost ecs2b -dbuser ensadmin -dbpass **** -dbname rat_Jun03_mk2 transcripts_to_delete.ls

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
my $stable_id;
my $config_dbname;

GetOptions( 'dbhost|host|h:s'        => \$host,
            'dbport|port|P:n'        => \$port,
            'dbname|db|D:s'        => \$dbname,
            'dbuser|user|u:s'        => \$user,
            'dbpass|pass|p:s'        => \$pass,
            'stable_id!'      => \$stable_id,
            'config_dbname:s' => \$config_dbname );


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

my $transcript_adaptor = $db->get_TranscriptAdaptor;

while (<>) {
  chomp;
  my $transcript_id = $_;

  eval {
    my $transcript;
    if ($stable_id) {
      $transcript = $transcript_adaptor->fetch_by_stable_id($transcript_id);
    } else {
      $transcript = $transcript_adaptor->fetch_by_dbID($transcript_id);
    }
    $transcript_adaptor->remove($transcript);
    print STDERR "Deleted $transcript_id\n";
  };
  if ($@) {
    print "Couldn't remove transcript $transcript_id ($@)\n";
  }
}

