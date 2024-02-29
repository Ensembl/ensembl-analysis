#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

  delete_transcripts.pl

=head1 DESCRIPTION

Deletes transcripts from given database whose IDs are passed in through
standard input or in a file named on the command line.

=head1 OPTIONS

=head1 

=head2 Database connection options

  -dbhost   Host name for database.
  -dbport   What port to connect to.
  -dbname   What name to connect to.
  -dbuser   What username to connect as.
  -dbpass   What password to use.

=head2 Database connection options (alternative method)

  --config_dbname   The alias for the database you
                    want to connect to, as defined in
                    Bio::EnsEMBL::Analysis::Config::Databases
                    configuration file.  Use this option or the
                    dbhost/dbport/dbname etc. options above, but not
                    both.

=head2 Other options

  --stable_id   A boolean flag to indicate that the IDs passed in are
                transcript stable IDs.

=head1 EXAMPLES

  ./delete_transcripts.pl --dbhost=genebuild2 \
    --dbuser=ensadmin --dbpass=**** \
    --dbname=rat_Jun03_mk2 transcripts_to_delete.list

=cut

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::Utilities;
use Bio::EnsEMBL::Utils::Exception;

$Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor::SILENCE_CACHE_WARNINGS = 1;

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
  die('Command line parsing error');

my $db;

if ( defined($config_dbname) ) {
  $db = get_db_adaptor_by_string($config_dbname);
}
elsif ( defined($dbname) && defined($host) ) {
  $db =
    new Bio::EnsEMBL::DBSQL::DBAdaptor( -host     => $host,
                                        -user     => $user,
                                        -port     => $port,
                                        -dbname   => $dbname,
                                        -pass     => $pass,
                                        -no_cache => 1 );
}
else {
  throw( "Need to pass either --dbhost, --dbuser, and --dbname, " .
         "or -config_dbname" );
}

my $ta = $db->get_TranscriptAdaptor();
my $ga = $db->get_GeneAdaptor();

my %gene_ids;

while ( my $transcript_id = <> ) {
  chomp($transcript_id);

  if ( !defined($transcript_id) || $transcript_id eq '' ) {
    last;
  }

  my $transcript;

  if ( defined($use_stable_ids) ) {
    $transcript = $ta->fetch_by_stable_id($transcript_id);
  }
  else {
    $transcript = $ta->fetch_by_dbID($transcript_id);
  }

  if ( !defined($transcript) ) {
    print("Can't fetch transcript $transcript_id\n");
    next;
  }

  my $gene_id =
    $ga->fetch_by_transcript_id( $transcript->dbID() )->dbID();

  my $has_transcript_stable_id = defined( $transcript->stable_id() );
  my $copy_of_dbID = $transcript->dbID();

  eval {
    # We want to update the gene in case the boundaries have changed
    $ta->remove($transcript, 1);
    if ($has_transcript_stable_id) {
      printf( "Deleted transcript %s (id = %d)\n",
              $transcript->stable_id(), $copy_of_dbID );
    }
    else {
      printf( "Deleted transcript (id = %d)\n", $copy_of_dbID );
    }
  };

  if ($@) {
    if ($has_transcript_stable_id) {
      printf( "Could not remove transcript %s (id = %d): %s\n",
              $transcript->stable_id(),
              $transcript->dbID(), $@ );
    }
    else {
      printf( "Could not remove transcript (id = %d): %s\n",
              $transcript->dbID(), $@ );
    }
    next;
  }

  $gene_ids{$gene_id}{ $copy_of_dbID } = 1;
} ## end while ( my $transcript_id...)

# Now get all those genes again and see if they are empty or split.

foreach my $gene_id ( keys(%gene_ids) ) {
  my $gene = $ga->fetch_by_dbID($gene_id);

  # See if gene is now empty.

  my $has_gene_stable_id = defined( $gene->stable_id() );

  my @transcripts = @{ $gene->get_all_Transcripts() };

  if ( !@transcripts ) {
    if ($has_gene_stable_id) {
      printf( "Gene %s (id = %d) is now empty, delete it!\n",
              $gene->stable_id(), $gene->dbID() );
    }
    else {
      printf( "Gene (id = %d) is now empty, delete it!\n",
              $gene->dbID() );
    }
    next;
  }

  # See if gene is now split.

  my $max_end;
  my @cluster;

  foreach my $t ( sort { $a->start() <=> $b->start() } @transcripts ) {
    if ( exists( $gene_ids{$gene_id}{ $t->dbID() } ) ) {
      next;
    }

    if ( !defined($max_end) ) {
      $max_end = $t->end();
    }
    elsif ( $t->start() > $max_end ) {
      # There's a gap between the ending of the last transcript and the
      # beginning of this, which means that the transcripts in @cluster
      # should now be removed from $gene and put into a new gene object.

      if ($has_gene_stable_id) {
        printf( "WARNING:\tFrom gene %s (id = %d),\n" .
                  "\t\ta new gene should be created\n" .
                  "\t\twith the following %d trancript(s):\n",
                $gene->stable_id(), $gene->dbID(), scalar(@cluster) );
        print(
          map {
            sprintf( "\t%s (id = %d)\n", $_->stable_id(), $_->dbID() )
          } @cluster );
      }
      else {
        printf( "WARNING:\tFrom gene (id = %d),\n" .
                  "\t\ta new gene should be created\n" .
                  "\t\twith the following %d trancript(s):\n",
                $gene->dbID(), scalar(@cluster) );
        print( map { sprintf( "\t(id = %d)\n", $_->dbID() ) }
               @cluster );

      }
      @cluster = ();
    } ## end elsif ( $t->start() > $max_end) [ if ( !defined($max_end...))]

    if ( $max_end < $t->end() ) {
      $max_end = $t->end();
    }

    push( @cluster, $t );
  } ## end foreach my $t ( sort { $a->start...})

} ## end foreach my $gene_id ( keys(...))
