#!/usr/bin/env perl
# $Source: /tmp/ENSCOPY-ENSEMBL-ANALYSIS/scripts/genebuild/delete_transcripts.pl,v $
# $Revision: 1.8 $

=head1 NAME

  delete_transcripts.pl

=head1 DESCRIPTION

  Deletes transcripts from given database whose ids are passed in
  through STDIN.

=head1 OPTIONS

=head1 

=head2 DB connection 

  -dbhost   Host name for database.
  -dbport   What port to connect to.
  -dbname   What name to connect to.
  -dbuser   What username to connect as.
  -dbpass   What password to use.

=head2 DB connection (alternative method)

  --config_dbname   The alias for the database you
                    want to connect to, as defined in
                    Bio::EnsEMBL::Analysis::Config::Databases
                    configuration file.  Use this option or the
                    dbhost/dbport/dbname etc. options above, but not
                    both.

=head2 Other options

  --stable_id   A boolean flag to indicate that the IDs passed in are
                stable IDs.

=head1 EXAMPLES

  ./delete_transcripts.pl --dbhost=genebuild2 \
    --dbuser=ensadmin --dbpass=**** \
    --dbname=rat_Jun03_mk2 transcripts_to_delete.list

=cut

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::Utilities;
use Bio::EnsEMBL::Utils::Exception;

my $host;
my $port = 3306;
my $dbname;
my $user;
my $pass;
my $use_stable_ids;
my $config_dbname;

GetOptions( 'dbhost|host|h:s' => \$host,
            'dbport|port|P:n' => \$port,
            'dbname|db|D:s'   => \$dbname,
            'dbuser|user|u:s' => \$user,
            'dbpass|pass|p:s' => \$pass,
            'stable_id!'      => \$use_stable_ids,
            'config_dbname:s' => \$config_dbname ) ||
  die('Command option parsing error');

my $db;

if (defined($config_dbname)) {
  $db = get_db_adaptor_by_string($config_dbname);
}
elsif ( defined($dbname) && define($host) ) {
  $db =
    new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $host,
                                        -user   => $user,
                                        -port   => $port,
                                        -dbname => $dbname,
                                        -pass   => $pass, );
}
else {
  throw( "Need to pass either -dbhost $host and -dbname $dbname or " .
         "-config_dbname $config_dbname for the script to work" );
}

my $ta = $db->get_TranscriptAdaptor();
my $ga = $db->get_GeneAdaptor();

while ( my $transcript_id = <> ) {
  chomp($transcript_id);

  my $transcript;

  if ( defined($use_stable_ids) ) {
    $transcript = $ta->fetch_by_stable_id($transcript_id);
  }
  else {
    $transcript = $ta->fetch_by_dbID($transcript_id);
  }

  if ( !defined($transcript) ) {
    warn("Can't fetch transcript $transcript_id\n");
    next;
  }

  eval {
    $ta->remove($transcript);
    print("Deleted $transcript_id\n");
  };

  if ($@) {
    warn("Couldn't remove transcript $transcript_id ($@)\n");
    next;
  }

  # See if gene is now empty or split.

  my $gene = $ga->fetch_by_transcript_id( $transcript->dbID() );

  my @transcripts = @{ $gene->get_all_Transcripts() };

  if ( !@transcripts ) {
    warn( sprintf( "Gene %s (id = %d) is now empty, delete it!\n",
                   $gene->stable_id(), $gene->dbID() ) );
    next;
  }

} ## end while ( my $transcript_id...)
