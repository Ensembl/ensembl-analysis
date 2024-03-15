#!/usr/bin/env perl


# Copyright [2017-2024] EMBL-European Bioinformatics Institute
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

  copy_transcripts.pl

=head1 DESCRIPTION

Based off the copy_genes script, this script will take a list of transcript dbIDs and copy them
to an output db as individual single transcript genes. This is useful for something like copying
canonical transcripts to an output database

The source database (from which you copy genes FROM) must contain DNA,
or else the script will throw.  If your source database does not contain
DNA, then provide details of a DNA database using the options --dnauser,
--dnahost, --dnadbname etc.

=head1 OPTIONS

  --logic               This option will change the analysis of the
                        genes and transcripts being written to the output database to an
                        analysis with the specified logic_name.

  --biotype             This option will change the biotype of the
                        genes and transcripts being written to the output database to the
                        biotype specified.

  --source              This option will change the source of the
                        genes and transcripts being written to the output database to the
                        source specified.

  --file                Read gene dbIDs out of a supplied file and only
                        copy these genes.

=head1 EXAMPLE

  perl copy_transcripts.pl \
    --sourcehost=srchost --sourceuser=ensro \
    --sourcedbname=core_db \
    --sourcedbport=XXX
    --targethost=trghost --targetuser=user --targetpass=XXX \
    --targetdbname=canonical \
    --targetdbport=XXX \
    --file=gene_ids_to_copy

=cut

use strict;
use warnings;
use feature 'say';

use Carp;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( throw warning verbose);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(empty_Transcript);

my $sourcehost   = '';
my $sourceuser   = 'ensro';
my $sourcedbname = '';
my $sourcepass   = undef;
my $sourceport;

my $outhost   = '';
my $outuser   = '';
my $outpass   = '';
my $outdbname = undef;
my $outport;

my $dnahost   = '';
my $dnauser   = 'ensro';
my $dnapass   = '';
my $dnadbname = undef;
my $dnaport;

my $in_config_name;
my $out_config_name;

my $logic;
my $biotype;
my $source;
my $infile;

GetOptions( 'inhost|sourcehost:s'                  => \$sourcehost,
            'inuser|sourceuser:s'                  => \$sourceuser,
            'indbname|sourcedbname:s'              => \$sourcedbname,
            'inport|sourceport:n'                  => \$sourceport,
            'inpass|sourcepass:s'                  => \$sourcepass,
            'outhost|targethost:s'                 => \$outhost,
            'outuser|targetuser:s'                 => \$outuser,
            'outpass|targetpass:s'                 => \$outpass,
            'outdbname|targetdbname:s'             => \$outdbname,
            'outport|targetport:n'                 => \$outport,
            'dnahost:s'                            => \$dnahost,
            'dnauser:s'                            => \$dnauser,
            'dnadbname:s'                          => \$dnadbname,
            'dnaport:n'                            => \$dnaport,
            'logic:s'                              => \$logic,
            'biotype:s'                            => \$biotype,
            'source:s'                             => \$source,
            'file:s'             => \$infile ) ||
  throw("Error while parsing command line options");



my $attach_dna_db = 1 ;


my $sourcedb;
my $dnadb;


$sourcedb =
    new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $sourcehost,
                                        -user   => $sourceuser,
                                        -pass   => $sourcepass,
                                        -port   => $sourceport,
                                        -dbname => $sourcedbname );
  if ($attach_dna_db) {
    if ( !defined($dnadbname) ) {
      my $dna_query = q(SELECT count(1) FROM dna);
      my $dna_sth   = $sourcedb->dbc()->prepare($dna_query);
      $dna_sth->execute();

      my $dna_count = $dna_sth->fetchrow();

      if ( !$dna_count ) {
        croak( "\nYour source database does not contain DNA. " .
               "Please provide a database with genomic sequences.\n" );
      }
    }
    else {
      $dnadb =
        new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $dnahost,
                                            -user   => $dnauser,
                                            -port   => $dnaport,
                                            -dbname => $dnadbname );

      $sourcedb->dnadb($dnadb);
    }
  }


my $outdb;
  $outdb =
    new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $outhost,
                                        -user   => $outuser,
                                        -pass   => $outpass,
                                        -port   => $outport,
                                        -dbname => $outdbname );

my $ta = $sourcedb->get_TranscriptAdaptor();

my @transcripts;
my @copy_transcripts;


my @transcript_ids;
if ( defined($infile) ) {
  open( IN, $infile ) or die("Can not open '$infile' for reading: $!");

  while ( my $line = <IN> ) {
    chomp($line);
    push( @transcript_ids, $line );
  }

  close(IN);
}
else {
  @transcript_ids = @{ $ta->list_dbIDs() };
}

printf( "Got %d transcript IDs\n", scalar(@transcript_ids) );

my $outga = $outdb->get_GeneAdaptor();
my $analysis;
if ( defined($logic) ) {
  $analysis = $outdb->get_AnalysisAdaptor->fetch_by_logic_name($logic);
  if ( !defined($analysis) ) {
    $analysis = Bio::EnsEMBL::Analysis->new( -logic_name => $logic, );
  }
}

my $genes_stored = 0;

foreach  my $transcript_id (@transcript_ids) {
  my $transcript = $ta->fetch_by_dbID($transcript_id);
  $transcript->load();
  if (defined($logic)) {
    $transcript->analysis($analysis);
  }

  if (defined($biotype)) {
    $transcript->biotype($biotype);
  }

  if (defined($source)) {
    $transcript->source($source);
  }

  empty_Transcript($transcript);
  my $gene = Bio::EnsEMBL::Gene->new();
  $gene->analysis($transcript->analysis);
  $gene->biotype($transcript->biotype);
  $gene->stable_id($transcript->stable_id);
  $gene->source($transcript->source);
  $gene->add_Transcript($transcript);
  $outga->store($gene);
  $genes_stored++;
} # foreach  my $transcript_id

say "Copy transcripts complete. Stored ".$genes_stored++." transcripts into new genes";

exit;
