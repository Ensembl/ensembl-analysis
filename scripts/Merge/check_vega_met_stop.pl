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

# $Id: check_vega_met_stop.pl,v 1.9 2013-03-05 13:47:18 fm2 Exp $

=pod

  This is a modified version of Julio's script with the same name.
  It takes one, potentially two, databases and a given list of
  transcript_ids and checks if the transcripts have met to stop
  integrity.

=cut

use strict;
use warnings;

use Carp;

use Getopt::Long qw(:config no_ignore_case);

use Bio::EnsEMBL::DBSQL::TranscriptAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

sub usage {
  my $indent = ' ' x ( length($0) );

  print <<USAGE_END;
USAGE:

  $0 --host [--port] --dbname \\
  $indent [--user] [--pass] \\
  $indent [--dna_host] [--dna_port] [--dna_dbname] \\
  $indent [list_of_dbIDs.in] >transcript_report.out

  $0 --help

EXAMPLE:

  $0 --host=genebuild2 --dbname=ak4_zebrafish_vega_fixed \\
  $indent --dna_dbname=ak4_zebrafish_ensembl \\
    <havana_coding_transcript_ids.txt |
  awk '\$3 == \$4 { print \$2, "is ok" }
       \$3 != \$4 { print \$2, "needs checking" }'

OPTIONS:

  --host    [required]  The database server host name.
  --port    [optional]  The database server port (default: 3306).
  --dbname  [required]  The database name.
  --user    [optional]  The database user (default: ensro).
  --pass    [optional]  Tho database user password (no default).

  --dna_host    [optional]  The DNA database server host name.
  --dna_port    [optional]  The DNA database server port.
  --dna_dbname  [optional]  The DNA database name.

  --help    Displays this helpful text.

  The defaults for the DNA database options are whatever values are
  specified for the equivalent options for the non-DNA database options.
  It is assumed that the user name and password may be reused for the
  DNA database.

INPUT:

  On the command line, the script expects the name of a file containing
  a list of internal IDs of protein coding transcripts.  If no file
  name is given, the script will fetch all transcripts with biotype
  "protein_coding" from the given database instead.

OUTPUT:

  By default, the script outputs the internal IDs of transcripts
  together with its stable ID and two codes.

  The first code has the following meanings:

    - ok        It has proper start and end.
    - start     Its start codon is partial.
    - stop      Its stop codon is partial.
    - both      Its start and stop codon is partial.

  The second code has the following meanings:

    - ok        Havana has not added either the cds_start_NF nor the
                cds_end_NF attributes to the transcript.

    - start     Havana has added the cds_start_NF attribute to the
                transcript.

    - stop      Havana has added the cds_end_NF attribute to the
                transcript.

    - both      Havana has added both the cds_start_NF and the
                cds_end_NF attributes to the transcript.

USAGE_END
} ## end sub usage

my $host;
my $port = 3306;
my $user = 'ensro';
my $pass = undef;
my $dbname;

my $dnadbname;
my $dnahost;
my $dnaport;

if ( !GetOptions( 'host|dbhost|h=s'       => \$host,
                  'port|dbport|P=s'       => \$port,
                  'user|dbuser|u=s'       => \$user,
                  'pass|dbpass|p=s'       => \$pass,
                  'dbname|db|D=s'     => \$dbname,
                  'dna_host=s'   => \$dnahost,
                  'dna_port=s'   => \$dnaport,
                  'dna_dbname=s' => \$dnadbname,
                  'help!'        => sub { usage(); exit(0); } ) ||
     !( defined($host) && defined($dbname) ) )
{
  usage();
  croak("Failed while parsing command line arguments");
}

my $db =
  new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $host,
                                      -port   => $port,
                                      -user   => $user,
                                      -pass   => $pass,
                                      -dbname => $dbname );

my $coord_system_adaptor = $db->get_CoordSystemAdaptor();

if ( defined($dnadbname) ) {
  my $dnadb =
    new Bio::EnsEMBL::DBSQL::DBAdaptor( -host => $dnahost || $host,
                                        -port => $dnaport || $port,
                                        -user => $user,
                                        -pass => $pass,
                                        -dbname => $dnadbname );

  my $dna_coord_system_adaptor = $dnadb->get_CoordSystemAdaptor();
  $db->dnadb($dnadb);
}

my $transcript_adaptor = $db->get_TranscriptAdaptor();

my $use_stdin = 1;
my $things    = [];

while ( defined( my $line = <> ) ) {
  if ( defined($line) ) {
    chomp($line);
    push( @{$things}, $line );
  }
}
if ( !@{$things} ) {
  $things = $transcript_adaptor->fetch_all_by_biotype('protein_coding');
  $use_stdin = 0;
}

while (1) {
  my $is_done = 0;

  my $dbID;
  my $transcript;

  if ($use_stdin) {
    $dbID = shift( @{$things} );
    if ( !defined($dbID) ) {
      $is_done = 1;
    }
    else {
      $transcript = $transcript_adaptor->fetch_by_dbID($dbID);
    }
  }
  else {
    $transcript = shift( @{$things} );
    if ( !defined($transcript) ) {
      $is_done = 1;
    }
    else {
      $dbID = $transcript->dbID();
    }
  }

  if ($is_done) { last }

  my $cds       = $transcript->translateable_seq();
  my $cds_start = substr( $cds, 0, 3 );
  my $cds_end   = substr( $cds, -3, 3 );

  my $start_partial = 0;
  my $stop_partial  = 0;

  my $partial_code;
  my $nf_code;

  if ( $cds_start ne "ATG" ) {
    $start_partial = 1;
  }
  if ( $cds_end ne "TGA" && $cds_end ne "TAA" && $cds_end ne "TAG" ) {
    $stop_partial = 1;
  }

  if ( !$start_partial && !$stop_partial ) {
    $partial_code = 'ok';
  }
  elsif ( $start_partial && !$stop_partial ) {
    $partial_code = 'start';
  }
  elsif ( !$start_partial && $stop_partial ) {
    $partial_code = 'stop';
  }
  else {
    $partial_code = 'both';
  }

  my $start_nf_attr = $transcript->get_all_Attributes('cds_start_NF');
  my $start_nf =
    ( @{$start_nf_attr} && $start_nf_attr->[0]->value() != 0 );

  my $stop_nf_attr = $transcript->get_all_Attributes('cds_end_NF');
  my $stop_nf =
    ( @{$stop_nf_attr} && $stop_nf_attr->[0]->value() != 0 );

  if ( !$start_nf && !$stop_nf ) {
    $nf_code = 'ok';
  }
  elsif ( $start_nf && !$stop_nf ) {
    $nf_code = 'start';
  }
  elsif ( !$start_nf && $stop_nf ) {
    $nf_code = 'stop';
  }
  else {
    $nf_code = 'both';
  }

  printf( "%d\t%s\t%s\t%s\n",
          $dbID, $transcript->stable_id(),
          $partial_code, $nf_code );

} ## end while (1)
