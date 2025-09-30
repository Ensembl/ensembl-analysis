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

#!/usr/bin/env perl

# This script calculates the repeat coverage of a genome which is repeatmasked by a specific analysis.
#
# output is in format : 
# 
#  Number bases masked = 267592
#  Number genes overlapped = 0
#  Total bases = 173499994
#  Total masked = 23746780
#  Total genes overlapped = 0
#
# To run on an assembly with coord_system.version COD_PRE and only toplevel regions and for the analysis repeatmasker_fish: 
#
# perl repeat_coverage.pl $dbcon -repeattypes repeatmasker_fish -path COD_PRE -coord_system toplevel 
# 

use strict;
use warnings;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Getopt::Long qw(:config no_ignore_case);


my $host   = '';
my $user   = 'ensro';
my $pass = '';
my $dbname = '';
my $port   = 3306;

my $dnahost   = '';
my $dnauser   = 'ensro';
my $dnapass = '';
my $dnadbname = undef;
my $dnaport = 3306;

my $genehost   = '';
my $geneuser   = 'ensro';
my $genepass = '';
my $genedbname = undef;
my $geneport   = 3306;

my @chromosomes;
my @genetypes = ('ensembl');
my @defrepeattypes = ('RepeatMask');
my @repeattypes;

my $coord_system= undef; # 'toplevel'
my $repeat_coord_system= undef; # 'contig';
my $path= '';

my $include_non_reference; # boolean flag
my $include_duplicates;


$| = 1;

GetOptions(
  'host|dbhost|h:s'       => \$host,
  'user|dbuser|u:s'       => \$user,
  'pass|dbpass|p:s'       => \$pass,
  'dbname|db|D:s'         => \$dbname,
  'path|cs_version:s'     => \$path,
  'port|dbport|P:n'       => \$port,
  'dnahost:s'             => \$dnahost,
  'dnauser:s'             => \$dnauser,
  'dnapass:s'             => \$dnapass,
  'dnadbname:s'           => \$dnadbname,
  'dnaport:n'             => \$dnaport,
  'genehost:s'            => \$genehost,
  'geneuser:s'            => \$geneuser,
  'genedbname:s'          => \$genedbname,
  'geneport:n'            => \$geneport,
  'chromosomes:s'         => \@chromosomes,
  'genetypes:s'           => \@genetypes,
  'repeattypes:s'         => \@repeattypes,
  'coord_system|cs_name:s'=> \$coord_system,          
  'rep_coord_system:s'    => \$repeat_coord_system,
  'include_non_reference' => \$include_non_reference, # boolean flag
  'include_duplicates'    => \$include_duplicates, # boolean flag
);

if(!$coord_system){
  throw("Specify -coord_system i.e. toplevel, chromosome, etc.");
}

if (scalar(@chromosomes)) {
  @chromosomes = split(/,/,join(',',@chromosomes));
}

if (scalar(@genetypes)) {
  @genetypes = split(/,/,join(',',@genetypes));
} 

if (scalar(@repeattypes)) {
  @repeattypes = split(/,/,join(',',@repeattypes));
}
else{
  @repeattypes = @defrepeattypes;
}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host   => $host,
  -user   => $user,
  -pass   => $pass,
  -port   => $port,
  -dbname => $dbname
);

my $genedb;
if ($genedbname) {
  print "Initing db for $genedbname\n";
  $genedb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host   => $genehost,
    -user   => $geneuser,
    -pass   => $genepass,
    -port   => $geneport,
    -dbname => $genedbname
  );
  print "Gene db = " . $genedb . "\n";
}
                   
my $dnadb;
if ($dnadbname) {
  $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host   => $dnahost,
    -user   => $dnauser,
    -pass   => $dnapass,
    -port   => $dnaport,
    -dbname => $dnadbname
  );

  $db->dnadb($dnadb);
}

my $sa     = $db->get_SliceAdaptor();

my $genesa;
if ($genedb) {
  $genesa = $genedb->get_SliceAdaptor();
}


my $seq_regions = $sa->fetch_all($coord_system, $path, $include_non_reference, $include_duplicates);

$| = 1;

my $total_masked = 0;
my $total_genes = 0;
my $total_overgene = 0;
my $total_bases = 0;
CHR:
foreach my $chr (@{$seq_regions}) {

  print STDERR "Slice " . $chr->seq_region_name . " from 1 to " . $chr->length. "\n";
  my $slicename = $chr->name();

  my $slice;
  my $geneslice;
  $slice = $sa->fetch_by_name($slicename);
  if ($genedb) {
    $geneslice = $genesa->fetch_by_name($slicename);
  }

  my @genes;
  if ($genedb) {
    print "Fetching genes\n";
  
    my %genes_hash;
    foreach my $genetype (@genetypes) {
      $genes_hash{$genetype} = $geneslice->get_all_Genes_by_type($genetype);
      print "Got " . scalar(@{$genes_hash{$genetype}}) . " $genetype genes\n";
      push @genes,@{$genes_hash{$genetype}};
    }
    print "Done fetching genes (fetched " . scalar(@genes) .")\n";
    $total_genes += scalar(@genes);
  } else {
    print "NO genedb specified so no genes fetched\n";
  }

  print "Fetching repeats\n";

  my @repeats;
  foreach my $repeattype (@repeattypes) {

    if ($repeat_coord_system && ($repeat_coord_system ne $coord_system)) {
      print STDERR "For repeattype $repeattype, repeat coord_system is $repeat_coord_system and gene coord_system is $coord_system\n";
      if ($repeat_coord_system ne 'contig') {
         # just a safety measure to make sure we are moving low->high
         # rather than high->low. 
         throw("Repeat_coord_system must be 'contig'");
      }

      # when $coord_system defined as 'chromosome' then
      # slice is a chromosome made up of contigs.
      my $contig_projection = $slice->project($repeat_coord_system);
      foreach my $segment (@$contig_projection) {
        my $contig = $segment->to_Slice();
        print $slice->seq_region_name(), ':', $segment->from_start(), '-',
              $segment->from_end(), ' -> ',
              $contig->seq_region_name(), ':', $contig->start(), '-',$contig->end(),
              ' ', $contig->strand(), "\n";
  
        # for human v55 the repeats were all stored on contig level so
        # need to get them up to toplevel
        foreach my $repeat (@{$contig->get_all_RepeatFeatures($repeattype)}) {
          my $transformed = $repeat->transform($coord_system);
          if (!defined $transformed) {
            warning("Transform of RepeatFeature start ".$repeat->start.
                    " end ".$repeat->end." on ".$contig->name." to $coord_system not possible"); 
          } else {
            push @repeats,$transformed;
          }
        }
  
      } # segment / contig
    } else {
      push @repeats,@{$slice->get_all_RepeatFeatures($repeattype)};
    }
  }
  print "Done fetching repeats (fetched " . scalar(@repeats) .")\n";

  @repeats = map { $_->[1] } sort { $a->[0] <=> $b->[0] } map { [$_->start, $_] } @repeats;

  my @repeat_blocks;

  my $curblock = undef;
  REPLOOP: foreach my $repeat (@repeats) {
    if ($repeat->start <= 0) { $repeat->start(1); }
    if (defined($curblock) && $curblock->end >= $repeat->start) {
      #print "Adding to block with " . $repeat->start . " to " . $repeat->end . "\n";
      #print "          Block was " . $curblock->start . " to " . $curblock->end . "\n";
      if ($repeat->end > $curblock->end) { $curblock->end($repeat->end); }
    } else {
      #print "Starting new block with " . $repeat->start . " to " . $repeat->end . "\n";
      #$curblock = Bio::EnsEMBL::SeqFeature->new(-START => $repeat->start, -END => $repeat->end);
      $curblock = Bio::EnsEMBL::Feature->new(-START => $repeat->start, -END => $repeat->end);
      push @repeat_blocks,$curblock;
    }
  }

  @repeat_blocks = map { $_->[1] } sort { $a->[0] <=> $b->[0] } map { [$_->start, $_] } @repeat_blocks;

  @genes = sort {$a->start <=> $b->start} @genes;
  my $ngene = scalar(@genes);


  my $novergene = 0;

  my $nmasked = 0;
  foreach my $block (@repeat_blocks) {
    $nmasked += $block->length;

   GENE: 
    for (my $ind=0; $ind < $ngene; $ind++) {
      my $gene = $genes[$ind];
      if (!defined($gene)) {
        next GENE;
      }
      if ($block->end >= $gene->start &&
        $block->start <= $gene->end) {
        foreach my $exon (@{$gene->get_all_Exons}) {
      #      print "\nprint block = $block\n";
      #      print "\nprint block seqname = ".$block->seqname."\n";
          if ($exon->overlaps($block)) {
            print "Overlap for gene " . get_gene_id($gene) . "\n";
            foreach my $trans (@{$gene->get_all_Transcripts}) {
              foreach my $tsf (@{$trans->get_all_supporting_features}) {
                print " Support " . $tsf->analysis->logic_name . " " . $tsf->hseqname . "\n";
              }
            }
            $novergene++;
            $genes[$ind] = undef;
            next GENE;
          }  
        }
      } elsif ($gene->start > $block->end) {
        last;
      } elsif ($gene->end < $block->start) {
        #print "Removing gene " . $genes[$ind]->stable_id . "\n";
        $genes[$ind] = undef;
      }
    }
  }
  print "Number bases masked = $nmasked\n";
  print "Number genes overlapped = $novergene\n";
 
  $total_overgene += $novergene;
  $total_masked += $nmasked;
  $total_bases += $slice->length;
}




my $ratio = ($total_masked / $total_bases) * 100  ;
print "Total bases = $total_bases\n";
print "Total masked = $total_masked\t";
print " ( $ratio % masked) \n" ;
print "Total genes = $total_genes\n";
print "Total genes overlapped = $total_overgene\n";
print "Done\n";





sub get_gene_id {
  my ($gene) = @_;
  return $gene->stable_id . " (" . $gene->dbID . ")" if ($gene->stable_id);
  return $gene->dbID;
}

sub get_transcript_id {
  my ($transcript) = @_;
  return $transcript->stable_id if ($transcript->stable_id);
  return $transcript->dbID;
}

sub get_exon_id {
  my ($exon) = @_;
  return $exon->stable_id if ($exon->stable_id);
  return $exon->dbID;
}


