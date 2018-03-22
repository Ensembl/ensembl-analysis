#!/usr/bin/env perl

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

use strict;
use Bio::EnsEMBL::Utils::Exception qw(stack_trace);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(fully_load_Gene);

my $dbname;
my $port = 3306;
my $host;
my $fdbuser;
my $fdbpass;
my $final_dbname;
my $final_port = 3306;
my $final_host;
my $write;
my $dna_dbname;
my $dna_port = 3306;
my $dna_host;
my $biotype;
my @chromosomes;
my $coordsystem ='toplevel';
my $gene_count;
my $start;
my $end;
my $all;
my $chrlist;
$| = 1;

&GetOptions(
  'fdbuser=s'  => \$fdbuser,
  'fdbpass=s'  => \$fdbpass,
  'write!'     => \$write,
  'fdbname=s'  => \$final_dbname,
  'fdbhost=s'  => \$final_host,
  'fdbport=s'  => \$final_port,
  'dnaname=s'  => \$dna_dbname,
  'dnahost=s'  => \$dna_host,
  'dnaport=s'  => \$dna_port,
  'dbname=s'   => \$dbname,
  'dbhost=s'   => \$host,
  'dbport=s'   => \$port,
  'biotype=s'  => \$biotype,  
  'chromosomes:s' => \@chromosomes,
  'chromosome_list:s' => \$chrlist,
  'coordsystem:s' => \$coordsystem,
  'all'        => \$all,
);

die ( "transfer_genes.pl 
  'write!'     => $write,
  'fdbpass=s'  => $fdbpass,
  'fdbuser=s'  => $fdbuser,
  'fdbname=s'  => $final_dbname,
  'fdbhost=s'  => $final_host,
  'fdbport=s'  => $final_port,
  'dnaname=s'  => $dna_dbname,
  'dnahost=s'  => $dna_host,
  'dnaport=s'  => $dna_port,
  'dbname=s'   => $dbname,
  'dbhost=s'   => $host,
  'dbport=s'   => $port,
  'biotype=s'  => $biotype,  
  'chromosomes:s' => @chromosomes,
  'chromosome_list:s' => $chrlist,
  'coordsystem:s' => $coordsystem,
  'all'        => (transfer all the genes in one go ),"
)
unless ( $fdbuser && $fdbpass && $final_port && $final_host && $final_dbname && $port && $host && $dbname );

print "Using data in $dbname\@$host:$port\n" ; 

my $sdb = new Bio::EnsEMBL::DBSQL::DBAdaptor
(
  -host   => $host,
  -user   => 'ensro',
  -port   => $port,
  -dbname => $dbname,
);


my $final_db = new Bio::EnsEMBL::DBSQL::DBAdaptor
(
  -host   => $final_host,
  -user   => $fdbuser,
  -port   => $final_port,
  -dbname => $final_dbname,
  -pass   => $fdbpass,
);
die ("Cannot find databases ") unless $final_db  && $sdb;
if ($dna_dbname){
  my $dna_db = new Bio::EnsEMBL::DBSQL::DBAdaptor
  (
    -host   => $dna_host,
    -user   => 'ensro',
    -port   => $dna_port,
    -dbname => $dna_dbname
  );
  $sdb->dnadb($dna_db);
}
if (-e $chrlist){
  open (LIST,$chrlist) or die ("Cannot open $chrlist\n");
  while (<LIST>){
    chomp;
    push @chromosomes, $_;
  }
}

my $final_ga = $final_db->get_GeneAdaptor;
my $ssa = $sdb->get_SliceAdaptor;
my $sga = $sdb->get_GeneAdaptor;
foreach my $chr (@{$ssa->fetch_all($coordsystem)}){
  my @genes;
  my @predicted_genes;
  if ($all){
    my $slice = $ssa->fetch_by_name($chr->name); 
    my @slicegenes;
    if ($biotype){
      @slicegenes = @{$slice->get_all_Genes_by_type($biotype,undef,1)};
    } else {
      @slicegenes = @{$slice->get_all_Genes(undef,undef,1)};
    }
    print STDERR $slice->name . " " . scalar(@slicegenes) . " genes\n";
    foreach my $gene (@slicegenes){
      if ($biotype){
        next unless $gene->biotype eq $biotype;
      }
      push @genes, fully_load_Gene($gene);
    }
  } else {   
    my $slice;
    foreach my $c ( @chromosomes ){
      $slice = $ssa->fetch_by_name($chr->name) if $c eq $chr->seq_region_name;
    }
    next unless $slice;
    my $out_slice;
    # fetch and write them by slice and lazy load them up front while they are still in memeory
    my @slicegenes;
    if ($biotype){
      @slicegenes = @{$slice->get_all_Genes_by_type($biotype,undef,1)};
    } else {
      @slicegenes = @{$slice->get_all_Genes(undef,undef,1)};
    }
    print STDERR $slice->name . " " . scalar(@slicegenes) . " genes\n";
    foreach my $gene (@slicegenes){
      push @genes, fully_load_Gene($gene);
    }
  }

  print "Starting with " . scalar(@genes) . " genes to copy \n";
  GENE:  foreach my $gene (@genes) {
    # lazyload it anyway - just in case
    if ($biotype){
      next unless $gene->biotype eq $biotype;
    }
    $gene->get_all_DBEntries;
    foreach my $trans (@{$gene->get_all_Transcripts}) {
      $trans->get_all_supporting_features;
      $trans->get_all_DBEntries;
      $trans->get_all_Attributes;    
      foreach my $exon (@{$trans->get_all_Exons}) {
        $exon->get_all_supporting_features;
        if($exon->start == 0){
          print STDOUT "Exon starts with zero - skipping\n";
          next GENE;
        }
      }
      $trans->translation;
    }
    push @predicted_genes, $gene;
  }
  if ($write) {
    foreach my $gene (@predicted_genes) {
      # store gene
      eval {
        $final_ga->store($gene);
      };
      if ( $@ ) {
        warn("UNABLE TO WRITE GENE:\n$@");
      }
#      print STDERR  "wrote gene $gene_count\tID " . $gene->dbID . "\n" ;
      $gene_count++;
    }
  } else {
    print "Not writing - write protect on!\n" ;
  }
  print "$gene_count genes copied\n";
}




