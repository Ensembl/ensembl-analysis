#!/usr/bin/env perl

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

use warnings ;
use vars qw(%Config);
use Bio::EnsEMBL::Analysis::Config::Databases;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use strict;
use Getopt::Long;


my $blast_db;
my $write;
my $chromosome;
my @thresholds ;
my $file;
my $t = '100,80,50';

my $usage = "perl label_models_by_coverage.pl
-db        $blast_db, Database containing models with BLAST hits (use the key from Databases.pm)
-write     $write, Write the results
-thresholds Comma separated array of thresholds to separate and label the models.
           For example the default is 100,80,50 and this will separate the models into:
           80-100
           50-80
           0-50
-sql       $file, output file to write the sql to.
The script will  cluster adjacent models that have the same blast hit
and sum the coverage across them - they will be labeled using the total hcoverage
and given the biotype 'fragment'
The script will output the sql to update the models
";

$| = 1;

&GetOptions(
	    'db:s'       => \$blast_db,
  	    'write!'       => \$write,
	    'chr:s'        => \$chromosome,
	    'thresholds:s' => \$t,
	    'sql:s'        => \$file,
   );

@thresholds = split(",",$t) if $t;
foreach my $t ( @thresholds ) {
  throw("Threshold has to be a number between 0 and 100\n") unless ( $t =~ /\d+/ && $t >= 0 && $t <= 100);
}

die($usage) unless ($blast_db && $file);
open(SQL,">$file") or throw("Cannot open file for writing $file\n");
# get database hash
my %database_hash;
my %databases = %{$DATABASES};
for my $category (keys %databases ) {
  if(!$databases{$category}{-host} || !$databases{$category}{-dbname}){
    print STDERR "Skipping category ".$category." no dbinfo ".
      "defined\n";
    next;
  }
  print STDERR "\n$category: $databases{$category}{-host} $databases{$category}{-dbname} :\n--------------------------------------\n" ;
  my %constructor_args = %{$databases{$category}};
  my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor
    (
     %constructor_args,
    );
  $database_hash{$category} = $dba;
}
# test that dbs all exist
my $blast = $database_hash{$blast_db};
throw("Db $blast_db not found in Databases.pm\n") unless $blast;

my $sa = $blast->get_SliceAdaptor;
my @new_genes;
my $count = 0 ;
CHR: foreach my $chr ( @{$sa->fetch_all('toplevel')}){
  if ( $chromosome ) {
    next CHR unless $chr->seq_region_name eq $chromosome;
    print STDERR  "Running on " . $chr->name ."\n";
  }
  my @genes = @{$chr->get_all_Genes(undef, undef, 1)};
  # do 1 strand at a time
  my @fwd;
  my @rev;
  foreach my $gene ( @genes ) {
    push @fwd, $gene if $gene->strand == 1;
    push @rev, $gene if $gene->strand == -1;
  }
  group_models(\@fwd);
  group_models(\@rev);
}

sub group_models {
  my ( $genes ) = @_;
  my $last_hit;
  my @groups;
  foreach my $gene (  sort { $a->start <=> $b->start } @$genes ) {
    # check 
    my $hit = $gene->get_all_Transcripts->[0]->get_all_supporting_features->[0];
    if ( $hit ) {
      print STDERR  "HIT " . $hit->hseqname ."\n";
      if ( $last_hit && $hit->hseqname ne $last_hit->hseqname ) {
	# label the group
	label_singletons($groups[0]) if scalar(@groups) == 1;
	label_fragments(\@groups) if scalar(@groups) > 1;
	# empty the array
	@groups = [];
	pop @groups;
      }
    } else {
      print STDERR  "No hit\n";
      if ( scalar(@groups) > 0 ){
	label_singletons($groups[0]) if scalar(@groups) == 1;
	label_fragments(\@groups) if scalar(@groups) > 1;
	# empty the array
	@groups = [];
	pop @groups;
      }
    }
    push @groups, $gene if $hit;
    $last_hit = $hit;
  }
  # end of chromosome - check for any left in the array
  if ( scalar(@groups) > 0 ){
    label_singletons($groups[0]) if scalar(@groups) == 1;
    label_fragments(\@groups) if scalar(@groups) > 1;
  }
}

sub label_singletons {
  my ($gene,$hcov,$qcov) = @_;
  print STDERR "Singleton\n";
  my $biotype = "fragment";
  my $hit = $gene->get_all_Transcripts->[0]->get_all_supporting_features->[0];
  unless ( $hcov && $qcov ) {
    $hcov = $hit->hcoverage;
    $qcov = $hit->percent_id;
    $biotype = $gene->biotype;
  }
  print STDERR $hit->hseqname ." HCOV QCOV $hcov $qcov\n";
  
  # what are our thresholds?
  @thresholds = sort {$b <=> $a} @thresholds;
  for ( my $i = 0 ; $i < $#thresholds ; $i++ ) {
    if ( ($hcov <= $thresholds[$i] or $qcov <= $thresholds[$i]) && ($hcov >= $thresholds[$i+1] && $qcov >= $thresholds[$i+1]) ){
      print SQL "UPDATE gene SET biotype = '" . $thresholds[$i+1] ."-". $thresholds[$i].":" . $biotype . "' WHERE gene_id = " . $gene->dbID .";\n";
      return;
    }
  }
  print SQL "UPDATE gene SET biotype = '0-". $thresholds[-1].":" . $biotype . "' WHERE gene_id = " . $gene->dbID .";\n";
}

sub label_fragments {
  my ($genes) = @_;
  print STDERR "Group\n";
  my  ( $hcov, $qcov ) ;
  my $total_hcoverage = 0;
  my @range;
  my $singletons = 0 ;
  foreach my $gene ( @$genes ) {
    my $hit = $gene->get_all_Transcripts->[0]->get_all_supporting_features->[0];
    $hcov = $hit->hcoverage;
    $qcov = $hit->percent_id;  
    print STDERR  $hit->hseqname ." HCOV QCOV $hcov $qcov\n";
    # make sure the hsps do not overlap
    push @range,[$hit->hstart,$hit->hend];
    $total_hcoverage += $hcov;
  }
  print STDERR "Total hcoverage = $total_hcoverage\n";
  $total_hcoverage = 100 if   $total_hcoverage >  100;
  my $last_end = 0;
  foreach my $hsp ( sort { $a->[0] <=> $b->[0] } @range ) {
    print STDERR  $hsp->[0] . " " . $hsp->[1] ."\n";
    if ( $hsp->[0] < $last_end - 10  ) {
      print STDERR "HSPs are not contiguous, treating as singletons\n";
      $singletons = 1;
      last;
    }
    $last_end = $hsp->[1];
  }
  if ( $singletons ) {
    foreach my $gene ( @$genes ) {
      label_singletons($gene);
    }
  } else {  
    foreach my $gene ( @$genes ) {
      label_singletons($gene,$total_hcoverage, 100);
    }
  }
}
