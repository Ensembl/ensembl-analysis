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

  check_split_genes.pl

=head1 DESCRIPTION

Checks if any gene in the given database is split. It deletes the split genes if parameter 'delete' is specified.

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
  --delete   Delete the split genes.

=head1 EXAMPLES

  ./check_split_genes.pl --dbhost=genebuild2 \
    --dbuser=ensro --dbpass=**** \
    --dbname=rat_Jun03_mk2

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
my $config_dbname;
my $delete = 0;

GetOptions( 'dbhost|host|h:s' => \$host,
            'dbport|port|P:n' => \$port,
            'dbname|db|D:s'   => \$dbname,
            'dbuser|user|u:s' => \$user,
            'dbpass|pass|p:s' => \$pass,
            'config_dbname:s' => \$config_dbname,
            'delete!'          => \$delete) ||
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

my $ga = $db->get_GeneAdaptor();

my %gene_ids;

# Get all the genes and see if they are split.

foreach my $gene_id (@{$ga->list_dbIDs()}) {
  my $gene = $ga->fetch_by_dbID($gene_id);
  my $has_gene_stable_id = defined( $gene->stable_id() );
  my $max_end;
  my @cluster;

  my @transcripts = @{ $gene->get_all_Transcripts() };

  foreach my $t (sort {$a->start() <=> $b->start()} @transcripts) {
    if (exists($gene_ids{$gene_id}{ $t->dbID()})) {
      next;
    }

    if (!defined($max_end)) {
      $max_end = $t->end();
    }
    elsif ($t->start() > $max_end) {
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
      
      # delete the split gene if 'delete' option
      if ($delete) {
      	my $geneid = $gene->dbID();
      	$ga->remove($gene);
      	print "Deleted $geneid\n";
      	if ($@) {
          print "Couldn't remove gene $gene_id ($@)\n";
        }
      }
      
      @cluster = ();
    } ## end elsif ( $t->start() > $max_end) [ if ( !defined($max_end...))]

    if ($max_end < $t->end()) {
      $max_end = $t->end();
    }

    push(@cluster,$t);
  } ## end foreach my $t ( sort { $a->start...})

} ## end foreach my $gene_id ( keys(...))
