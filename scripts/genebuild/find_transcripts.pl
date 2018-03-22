#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

  find_transcripts.pl

=head1 DESCRIPTION

This script takes database options and a file of transcript stable ids
from the source database and finds equivalent transcripts based on the exon structure in the target database.

The source database must contain DNA,
or else the script will throw.  If your source database does not contain
DNA, then provide details of a DNA database using the options --dnauser,
--dnahost, --dnadbname etc.

=head1 OPTIONS

  --all                 This option will find all transcripts from
                        $sourcedbname in the target db.

  --file                Read gene dbIDs out of a supplied file and only
                        find these genes.

=head1 EXAMPLE

  perl find_transcripts.pl \
    --sourcehost=srchost --sourceuser=ensro \
    --sourcedbname=core_db \
    --targethost=trghost --targetuser=user --targetpass=XXX \
    --targetdbname=ens_db \
    --file=transcript_sids \
    [--nobiotype]

=cut

use strict;
use warnings;

use Carp;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( throw warning verbose);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils;
use Bio::EnsEMBL::Analysis::Tools::Utilities;

my $sourcehost   = '';
my $sourceuser   = 'ensro';
my $sourcedbname = '';
my $sourcepass   = undef;
my $sourceport   = 3306;

my $outhost   = '';
my $outuser   = '';
my $outpass   = '';
my $outdbname = undef;
my $outport   = 3306;

my $dnahost   = '';
my $dnauser   = 'ensro';
my $dnapass   = '';
my $dnadbname = undef;
my $dnaport   = 3306;

my $infile;
my $all;

my $nobiotype = 0; # include biotype as part of the exon key by default

verbose('EXCEPTION');

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
            'nobiotype!'                           => \$nobiotype,
            'all!'                                 => \$all,
            'file:s'                               => \$infile ) ||
  throw("Error while parsing command line options");

my $attach_dna_db = 1 ;

if ( $all && defined($infile) ) {
  throw("Specify either --all or --file, but not both");
}
elsif ( !$all && !defined($infile) ) {
  throw( "Specify either --all " .
         "(to find all genes from $sourcedbname) " .
         "or --file (to copy a list of transcript ids)" );
}

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
    $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $dnahost,
                                            -user   => $dnauser,
                                            -port   => $dnaport,
                                            -dbname => $dnadbname );
    $sourcedb->dnadb($dnadb);
  }
}

my $outdb = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $outhost,
                                        -user   => $outuser,
                                        -pass   => $outpass,
                                        -port   => $outport,
                                        -dbname => $outdbname );

my $ta = $sourcedb->get_TranscriptAdaptor();

my @source_transcript_ids;
my @target_transcript_ids;

if ( defined($infile) ) {
  open( IN, $infile ) or die("Can not open '$infile' for reading: $!");

  while ( my $line = <IN> ) {
    chomp($line);
    push(@source_transcript_ids,$line);
  }
  close(IN);
}
else {
  @source_transcript_ids = @{ $ta->list_stable_ids() };
}

my $outta = $outdb->get_TranscriptAdaptor();
@target_transcript_ids = @{ $outta->list_dbIDs() };

printf( "Got %d source transcript IDs\n", scalar(@source_transcript_ids) );
printf( "Got %d target transcript IDs\n", scalar(@target_transcript_ids) );

my %target_transcripts_ek = ();
foreach my $t (@{$outta->fetch_all_by_dbID_list([@target_transcript_ids])}) {
  $target_transcripts_ek{get_transcript_exon_key($t,$nobiotype)} = $t->dbID();
}

foreach my $t (@{$ta->fetch_all_by_stable_id_list([@source_transcript_ids])}) {
  my $source_t_ek = get_transcript_exon_key($t,$nobiotype);
  
  if ($target_transcripts_ek{$source_t_ek}) {
  	print($t->dbID()." ".$target_transcripts_ek{$source_t_ek}."\n");
  } else {
  	print("Couldn't find transcript ".$t->stable_id()." in target database.\n");
  }
}

sub get_transcript_exon_key {
  my $transcript = shift;
  my $nobiotype = shift;
  
  my $biotype_str = $transcript->biotype;
  if ($nobiotype) {
  	$biotype_str = '';
  }
  
  my $string = $transcript->slice->seq_region_name.":".$biotype_str.":".$transcript->seq_region_start.":".$transcript->seq_region_end.":".$transcript->seq_region_strand.":".@{$transcript->get_all_translateable_Exons()}.":";

  my $exons = sort_by_start_end_pos($transcript->get_all_Exons());
  foreach my $exon (@{$exons}) {
    $string .= ":".$exon->seq_region_start.":".$exon->seq_region_end;
  }

  return $string;
}

sub sort_by_start_end_pos {
  my ($unsorted) = @_;

  my @sorted = sort { if ($a->seq_region_start < $b->seq_region_start) {
        return -1;
    } elsif ($a->seq_region_start == $b->seq_region_start) {
      if ($a->seq_region_end < $b->seq_region_end) {
        return-1;
      } elsif ($a->seq_region_end == $b->seq_region_end) {
        return 0;
      } elsif ($a->seq_region_end > $b->seq_region_end) {
        return 1;
      }
        return 0;
    } elsif ($a->seq_region_start > $b->seq_region_start) {
        return 1;
    }
  } @$unsorted;

  return \@sorted;
}

1;
