# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

#!/usr/bin/env perl

# Script to check the integrity of some or all of the genes in an Ensembl 
# database 

# any questions please send to http://lists.ensembl.org/mailman/listinfo/dev

=head1 NAME

check_GBGenes.pl

=head1 SYNOPSIS

perl check_GBGenes.pl -chromosome 1 -coordsystem chromosome -transcripts

This will run the transcript checks on chromosome 1 on the final database and
reference database specified in Bio::EnsEMBL::Pipeline::Config::GeneBuild::Databases

Note: If submitting to LSF you will need to specify memory requiremensts of >100MB
      For instance, 1GB suffices all human chromosomes:

      bsub -R"select[mem>1000] rusage[mem=1000]" -M1000000 

=head1 DESCRIPTION

This script runs a series of tests on genes in the database it is pointed
too. These tests are described in more detail in the various perl
modules in this directory

=head1 OPTIONS

  -dbhost database host
  -dbuser database user
  -dbpass database password
  -dbport database port
  -dbname database name
  -dnahost host for dna database
  -dnaport port for dna database
  -dnadbname name for dna database

  These settings by default are taken from:
  Bio::EnsEMBL::Analysis::Config::Databases

  -ignorewarnings flag to specify whether to ignore the warnings in the 
                  code

  This is take by default from:
  Bio::EnsEMBL::Pipeline::Config::GeneBuild::GeneBuilder


  -chromosome name of seq region to get genes from
  -coordsystem name of the coordinate system the seq region belongs to

  these are both obligatory options

  -chrstart the start coordinate of the piece of seq region to fetch
  -chrend the end coordinate of the piece of seq region to fetch

  without these the whole seq region if fetched

  -duplicates check for duplicate exons (off by default)
  -transcripts check the transcripts (on by default)

  -help print the perl docs

=head1 EXAMPLES

perl check_GBGenes.pl -chromosome 1 -coordsystem chromosome -transcripts

runs the standard checks on chromosomes using the databases defined in
configuration


perl check_GBGenes.pl -chromosome 3 -coordsystem chromosome -duplicates
-host yourhost -user youruser -port 3306 -dbname yourdb -dnahost yourhost
-dnaport 3306 -dnadbname yourdnadb

runs the standard checks and checks for duplicate exons on the database
specified on the commandline

=head1 NOTES

Note this script uses several settings from the
Bio::EnsEMBL::Analysis::Config::GeneBuild::GeneBuilder

Also this scripts uses several modules which can be found in this
directory so you need to make sure you run the script in this
directory or alternatively put this directory in your PERL5LIB

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use TranscriptChecker;
use ContigGenesChecker;
use GeneChecker;
use buildchecks::ScriptUtils;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::BuildChecks qw (
                                                               MAXSHORTINTRONLEN
                                                               MINSHORTINTRONLEN
                                                               MINLONGINTRONLEN
                                                               MAX_EXONSTRANSCRIPT
                                                               MAXSHORTEXONLEN
                                                               MINSHORTEXONLEN
                                                               MINLONGEXONLEN 
                                                               MAXTRANSCRIPTS
                                                               MINTRANSLATIONLEN
                                                               IGNOREWARNINGS
                                                               MAXGENELEN
                                                              ); ;
use Bio::EnsEMBL::Analysis::Tools::Utilities;


$| = 1;

my $host;
my $dbname;
my $user;
my $pass;
my $port='3306';
my $gene_database_name;
my $path;

# default path comes out of database

my $maxshortintronlen  =  50;
if (defined($MAXSHORTINTRONLEN)) { $maxshortintronlen  = $MAXSHORTINTRONLEN; }

my $minshortintronlen  =  3;
if (defined($MINSHORTINTRONLEN)) { $minshortintronlen  = $MINSHORTINTRONLEN; }

my $minlongintronlen   =  100000;
if (defined($MINLONGINTRONLEN )) { $minlongintronlen   = $MINLONGINTRONLEN; }

my $maxexonstranscript =  150;
if (defined($MAX_EXONSTRANSCRIPT)) { $maxexonstranscript = $MAX_EXONSTRANSCRIPT; }

my $maxshortexonlen    =  10;
if (defined($MAXSHORTEXONLEN)) { $maxshortexonlen    = $MAXSHORTEXONLEN; }

my $minshortexonlen    =  3;
if (defined($MINSHORTEXONLEN)) { $minshortexonlen    = $MINSHORTEXONLEN; }

my $minlongexonlen     =  50000;
if (defined($MINLONGEXONLEN )) { $minlongexonlen     = $MINLONGEXONLEN; }

my $maxtranscripts     =  10; 
if (defined($MAXTRANSCRIPTS)) { $maxtranscripts     = $MAXTRANSCRIPTS; }

my $mintranslationlen  =  10; 
if (defined($MINTRANSLATIONLEN)) { $mintranslationlen  = $MINTRANSLATIONLEN; }

my $ignorewarnings     =  0; 
if (defined($IGNOREWARNINGS)) { $ignorewarnings     = $IGNOREWARNINGS; }

my $maxgenelen     =  2_000_000; 
if (defined($MAXGENELEN)) { $maxgenelen     = $MAXGENELEN; }

my @chromosomes;

my $specstart = 1;
my $specend   = undef;

my $dnadbname = "";
my $dnahost   = "";
my $dnaport   = '3306';

my $exon_dup_check = 0;

my $check_transcripts = 1;
my $schema = 20;
my $coordsystem = 'chromosome';
my $help;
my @genetypes;

GetOptions(
            'dbhost|host|h:s'    => \$host,
            'dbuser|user|u:s'    => \$user,
            'dbpass|pass|p:s'    => \$pass,
            'dbport|port|P:n'    => \$port,
            'dbname|db|D:s'         => \$dbname,
            'dnahost:s'        => \$dnahost,
            'dnaport:n'        => \$dnaport,
            'dnadbname:s'      => \$dnadbname,
            'gene_database:s'  => \$gene_database_name,
            'path:s'           => \$path, 
            'ignorewarnings!'  => \$ignorewarnings,
            'chromosomes:s'    => \@chromosomes,
            'coordsystem|cs_name:s'    => \$coordsystem,
            'chrstart:n'       => \$specstart,
            'chrend:n'         => \$specend,
            'schema:n'         => \$schema,
            'genetypes:s'      => \@genetypes,
            'duplicates!'      => \$exon_dup_check,
            'transcripts!'     => \$check_transcripts,
            'help!'            => \$help,
           ) or perldocs("Failed to get options");

if ((!defined($host) || !defined($dbname))){
  if(!defined($gene_database_name)) {
    die "ERROR: Must at least set host (-host), dbname (-dbname) or -gene_database\n" .
      "       to allow the connection to be fetched from\n".
        "     Bio::EnsEMBL::Analysis::Config::Databases.pm \n";
  }
}

if($help){
  perldocs();
}

if (scalar(@chromosomes)) {
  @chromosomes = split(/,/,join(',',@chromosomes));
}

if (scalar(@genetypes)) {
  @genetypes = split(/,/,join(',',@genetypes));
}

my $db;

if($gene_database_name){
  $db = get_db_adaptor_by_string($gene_database_name);
}else{
  $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $host,
                                           -user => $user,
                                           -port => $port,
                                           -dbname => $dbname,
                                           -pass => $pass);
  if ($dnadbname ne "") {
    if ($dnahost eq "") {
      $dnahost = $host;
    }
    my $dnadbase = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $dnahost,
                                                      -user   => $user,
                                                      -port   => $dnaport,
                                                      -dbname => $dnadbname,
                                                      -pass   => $pass,
                                                     );
    $db->dnadb($dnadbase);
  }
}

my $sa = $db->get_SliceAdaptor();
my $csa = $db->get_CoordSystemAdaptor();

if (!defined($path)) {
  $path = $csa->fetch_by_name($coordsystem)->version;
}


#Not practical to do any other way
if ($exon_dup_check) {
  print "Performing exon duplicate check for ALL exons\n";
  find_duplicate_exons($db);
  print "Done duplicate check\n";
}

my $chrhash;

if ($schema == 20) {
  $chrhash = get_chrlengths_v20($db, $path,$coordsystem);
} else {
  $chrhash = get_chrlengths_v19($db, $path);
}

filter_to_chr_list(\@chromosomes,$chrhash,$db->dbc->dbname);

#print "Start $specstart End $specend\n";


my @failed_transcripts;

my $total_transcripts_with_errors = 0;
my $total_genes_with_errors = 0;
my $total_genes = 0;
my $total_transcripts = 0;

# Begin testing genes
foreach my $chr (@{sort_chr_names($chrhash)}) {

  my $chrstart = $specstart;
  my $chrend = (defined ($specend) && $specend < $chrhash->{$chr}) ? $specend :
               $chrhash->{$chr};

  my $slice;
  my $slicename;
  if ($schema == 20) {
    $slicename = "$coordsystem:$path:$chr:$chrstart:$chrend:1";
    print "Slice = $slicename\n";
    $slice = $sa->fetch_by_name($slicename);
  } else {
    $slicename = "$chr:$chrstart-$chrend";
    print "Slice = $slicename\n";
    $slice = $sa->fetch_by_chr_start_end($chr,$chrstart,$chrend);
  }
  
  my $genes;

  if (!scalar(@genetypes)) {
    $genes = $slice->get_all_Genes;
  } else {

    foreach my $type (@genetypes) {
      push @$genes, @{$slice->get_all_Genes_by_type($type)};
    }
  }

  my $cgc = new ContigGenesChecker(
                                     -slice          => $slice,
                                     -genes          => $genes,
                                     -ignorewarnings => $ignorewarnings,
                                     -adaptor        => $db, 
                                     );
  $cgc->check;

  if ($cgc->has_Errors()) {
    print "--------------------------------------\n";
    print "Slice with errors: $slicename\n";
    $cgc->output;
  }

  GENE: foreach my $gene (@$genes) {
    $total_genes++;
  
    my $gc = new GeneChecker(-gene          => $gene,
                             -maxtransgene  => $maxtranscripts,
                             -maxgenelen    => $maxgenelen,
                             -ignorewarnings=> $ignorewarnings,
                             -adaptor       => $db, 
                             -slice         => $slice,
                            );
    $gc->check;
    my $nwitherror = 0;

    if ($gc->has_Errors()) {
      $gc->output;

      $total_genes_with_errors++;
      $nwitherror=1; 
    }
  

    if ($check_transcripts) {
      my @trans = @{$gene->get_all_Transcripts()};

      TRANSCRIPT: foreach my $transcript (@trans) {
        $total_transcripts++;
        my $tc = new TranscriptChecker(
                                       -transcript         => $transcript,
                                       -minshortintronlen  => $minshortintronlen,
                                       -maxshortintronlen  => $maxshortintronlen,
                                       -minlongintronlen   => $minlongintronlen,
                                       -minshortexonlen    => $minshortexonlen,
                                       -maxshortexonlen    => $maxshortexonlen,
                                       -minlongexonlen     => $minlongexonlen,
                                       -mintranslationlen  => $mintranslationlen,
                                       -maxexonstranscript => $maxexonstranscript,
                                       -ignorewarnings     => $ignorewarnings,
                                       -adaptor            => $db, 
                                       -slice              => $slice,
                                       );
        $tc->check;
  
        if ($tc->has_Errors()) {
          $total_transcripts_with_errors++;
          if (!$nwitherror) {
            $gc->output;
            $total_genes_with_errors++;
          }
          $tc->output;
          # Don't store for now!!!  push @failed_transcripts,$tc;
          $nwitherror++;
        }
      }
    }
  }
}

print "Summary:\n";
print "Number of genes checked           = $total_genes\n";
print "Number of transcripts checked     = $total_transcripts\n\n";
print "Number of transcripts with errors = $total_transcripts_with_errors\n";
print "Number of genes with errors       = $total_genes_with_errors\n\n";

#End of main

sub find_duplicate_exons {
  my $db = shift;
  
  if (!$db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
    die "find_duplicate_exons should be passed a Bio::EnsEMBL::DBSQL::DBAdaptor\n";
  }

  my $q = qq( SELECT e1.exon_id, e2.exon_id 
              FROM exon e1, exon e2 
              WHERE e1.exon_id<e2.exon_id AND e1.seq_region_start=e2.seq_region_start AND 
                    e1.seq_region_end=e2.seq_region_end AND e1.seq_region_id=e2.seq_region_id AND 
                    e1.seq_region_strand=e2.seq_region_strand AND e1.phase=e2.phase AND
                    e1.end_phase=e2.end_phase
              ORDER BY e1.exon_id 
            ); 
  my $sth = $db->dbc->prepare($q) || $db->throw("can't prepare: $q");
  my $res = $sth->execute || $db->throw("can't execute: $q");
 
  my $ndup = 0;
  while( my ($exon1_id, $exon2_id) = $sth->fetchrow_array) {
    print "ERROR: Exon duplicate pair: $exon1_id and $exon2_id\n"; 
    $ndup++;
  }
  print "Total number of duplicate pairs = $ndup\n";
}

sub perldocs{
  my ($msg) = @_;
  print $msg."\n" if($msg);
  exec('perldoc', $0);
  exit(0);
}
