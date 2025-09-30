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

=pod

  This script compares the loaded patches in a db to the ones
  downloaded from the original source. The regions that are compared
  include the head or tail of the patch.

  The script creates subdirectories using the downloaded sequences from the
  web in addition to creating a fasta file from the dumped sequences from the
  database. These are all temporary dirs and files and can be deleted after
  the test is done.

=cut

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;

$| = 1;

my $host;
my $user;
my $pass;
my $port;
my $dbname;
my $coord_system     = 'chromosome';
my $alt_scaf_file    = ''; # the alt_scaffold_placement.txt file from GRC. Describes how scaffolds fit into chromosomes.
my $write_dir        = ''; # a directory that various subdirectories and files will be written to while this script runs
my $downloaded_fasta = ''; # the downloaded scaffold/patch fasta from GRC
my $exploded_dir;          # subdir that will hold the exploded GRC file (i.e. a file for each patch)
my $subseq_dir;            # as exploded but with seqs trimmed to only the toplevel sections of the scaffolds
my $revcomp_dir;           # the revcomp seq if the scaffold has orientation -1 relative to the chromosome level
my $dna_dbname;
my $dna_host;
my $dna_port;
my $dna_user;

GetOptions( 'dbhost:s'           => \$host,
            'dbuser:s'           => \$user,
            'dbname:s'           => \$dbname,
            'dbpass:s'           => \$pass,
            'dbport:n'           => \$port,
            'downloaded_fasta:s' => \$downloaded_fasta,
            'write_dir:s'        => \$write_dir,
            'alt_scaf_file:s'    => \$alt_scaf_file,
            'dnadbname:s'        => \$dna_dbname,
            'dnahost:s'          => \$dna_host,
            'dnaport:s'          => \$dna_port,
            'dnauser:s'          => \$dna_user );

my $db =
  new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $host,
                                      -user   => $user,
                                      -pass   => $pass,
                                      -port   => $port,
                                      -dbname => $dbname );

my $dnadb;
if ($dna_dbname) {
  $dnadb =
    new Bio::EnsEMBL::DBSQL::DBAdaptor(-dbname => $dna_dbname,
                                       -host   => $dna_host,
                                       -port   => $dna_port,
                                       -user   => $dna_user);
  $db->dnadb($dnadb);
}

if ( !-d $write_dir ) {
  throw("Directory to write seq files to doesn't exist\n");
} else {
  $subseq_dir   = $write_dir . "/subseq";
  $exploded_dir = $write_dir . "/exploded";
  $revcomp_dir  = $write_dir . "/revcomp";

  if ( -d $revcomp_dir or -d $exploded_dir ) {
    throw( "Revcomp or exploded sub dir already exists in write dir. "
         . "Remove and restart.\n" );
  } else {
    mkdir $subseq_dir;
    mkdir $revcomp_dir;
    mkdir $exploded_dir;
  }
}

if ( !-d $subseq_dir ) {
  throw("Directory to hold subseqs doesn't exist\n");
}
if ( !-d $revcomp_dir ) {
  throw("Directory to hold revcomp seqs doesn't exist\n");
}
if ( !-d $exploded_dir ) {
  throw("Directory to hold exploded seqs doesn't exist\n");
}


my $acc_header = "downloaded_fasta_acc_header";
open( ACC_HEADER, ">" . $write_dir . "/" . $acc_header ) || die "Could not open file $acc_header";

open( DOWNLOAD_SEQ, "<" . $downloaded_fasta ) || die( "Could not open file $downloaded_fasta");
while(my $line = <DOWNLOAD_SEQ>){
  $line =~ s/\>\w+\|\w+\|\w+\|(\w+\.\w+).*/\>$1/;
  # should end up with acc beginning GL as header and
  # file name (GL339449.1 for example) matches entry in
  # alt_scaffold_placement.txt
  print ACC_HEADER $line;
}
close(DOWNLOAD_SEQ);
close(ACC_HEADER);

system("fastaexplode -f ".$write_dir."/".$acc_header." -d $exploded_dir");
my $seqIOout = Bio::SeqIO->new(
                             -file => '>' . $write_dir . '/all_patches_db.fa',
                             -format => 'fasta' );

my $sa = $db->get_SliceAdaptor;

open( ALT_SCAF, "<" . $alt_scaf_file ) || die("Could not open file $alt_scaf_file");
REGION: while (<ALT_SCAF>) {
  next REGION if m/^#/;
  chomp;

  my @arr          = split(/\t/, $_);
  my $region_name  = "CHR_".$arr[2];
  my $acc          = $arr[3];
  my $strand       = $arr[8];
  my $alt_start    = $arr[9];
  my $alt_stop     = $arr[10];
  my $parent_start = $arr[11];
  my $parent_stop  = $arr[12];
  my $head         = $arr[13];
  my $tail         = $arr[14];

  if ( $strand eq "+" ) {
    $strand = 1;
  } elsif ( $strand eq "-" ) {
    $strand = -1;
  } else {
    print "WARNING: Strand $strand info is wrong in alt scaf file "
        . "for $region_name, defaulting to 1\n";
    $strand = 1;
  }

  # get this region from db
  my $alt_chrom_slice = $sa->fetch_by_region('chromosome',  $region_name);
  my $asm_exc_adaptor = $db->get_AssemblyExceptionFeatureAdaptor();
  my @exc_feat = @{$asm_exc_adaptor->fetch_all_by_Slice($alt_chrom_slice)};
  my $ref_slice = $exc_feat[0]->alternate_slice();
  @exc_feat = @{$asm_exc_adaptor->fetch_all_by_Slice($ref_slice)};
  my $slice;
  foreach my $exc (@exc_feat){
    if($exc->alternate_slice->seq_region_name eq $region_name){
      $slice = $exc->alternate_slice();
    }  
  }
  print $slice->name . " is db slice name\n";
  $seqIOout->write_seq($slice);


  # revcomp dowloaded if necessary
  if ( $strand == -1 ) {
    system( "fastarevcomp -f " . $exploded_dir . "/" . $acc . ".fa  > " . $revcomp_dir . "/" . $acc . ".fa" );
    system( "cat " . $revcomp_dir . "/" . $acc . ".fa >> " . $write_dir . "/all_patches_grc.fa" );
  } else {
    system( "cat " . $exploded_dir . "/" . $acc . ".fa >> " . $write_dir . "/all_patches_grc.fa" );
  }

} ## end while (<ALT_SCAF>)
close(ALT_SCAF);

system( "fastadiff -1 " . $write_dir . "/all_patches_db.fa -2 " . $write_dir .  "/all_patches_grc.fa -a -c FALSE >> " . $write_dir . "/differences.out" );

print "NOTE: differences should have been written to "
    . $write_dir . "/differences.out. "
    . "Check the file for the fastadiff output.\n";
