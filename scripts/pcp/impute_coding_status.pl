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

# Script to impute proteing coding status into pcp database
# Reads in and validates rnasamaba and cpc2 output files and creates intersection
# Uses this to update biotype (impute proteincodign stayus effectively)

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case);

use Bio::EnsEMBL::DBSQL::DBAdaptor;

## Some defaults
my $coord_system = 'toplevel';
my $dbname = '';
my $user   = 'ensadmin';
my $host = $ENV{GBS5};
my $port = $ENV{GBP5};
my $pass = 'ensembl';
my ($cpc2_file, $rnasamba_file ) = ( '', '' );

my $options = GetOptions ("user|dbuser|u=s"     => \$user,
                          "host|dbhost|h=s"     => \$host,
                          "port|dbport|P=i"     => \$port,
                          "dbname|db|D=s"       => \$dbname,
                          "dbpass|pass|p=s"     => \$pass,
                          "cpc2=s"              => \$cpc2_file,
                          "rnaSamba|rnas=s"     => \$rnasamba_file,
                          "coords:sub"          => \$coord_system);


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $port,
  -user    => $user,
  -host    => $host,
  -dbname  => $dbname,
  -pass    => $pass);

my %cpc_results = parse_results($cpc2_file, 8);
my %rnasamba_results = parse_results($rnasamba_file, 2);
my @selected_genes = ();

if (%cpc_results ne %rnasamba_results) {
    die "Results files have different number of genes\n";
} else {
    my %compare = map { $_ => 1 } keys %cpc_results;
    for my $key (keys %rnasamba_results) {
        last unless exists $compare{$key};
        delete $compare{$key};
    }
    if (%compare) {
        die "Results files contain different genes\n";
    } else {
        ## SUCCESS
        print "Same genes compared, nothing missing\n";
    }
}

foreach my $key (keys %cpc_results) {
  ## Currently using AND to account for relative sensitivity  / specifity differences
  ## between the algorithms.
  if ($rnasamba_results{$key} eq 'coding' and $cpc_results{$key} eq 'coding') {
      push(@selected_genes, $key);
  }
}

my $slice_adaptor = $db->get_SliceAdaptor();
my $slices = $slice_adaptor->fetch_all('toplevel', undef, 1 );
my $gene_adaptor = $db->get_GeneAdaptor();

while ( my $slice = shift @{$slices} ) {
  my $genes = $slice->get_all_Genes();
  while ( my $gene = shift @{$genes} ) {
    my $check = $gene->dbID;
    if ( grep( /^$check$/, @selected_genes ) ) {
      my $current = $gene->biotype;
      $gene->biotype('pcp_protein_coding');
      $gene_adaptor->update($gene);
    }
  }
}

## Parses results from the classification program output
## Takes the filename and column where protein coding classifier stored
sub parse_results {
  my ($in_file, $column_num) = @_;
  my %results;

  open my $FILE, '<', $in_file or die $!;
  my $header = <$FILE>;
  ## header check maybe?

  while( <$FILE> ) {
    chomp( my @row = split'\t', $_ );
    my @temp = split':', $row[0];
    $results{$temp[-1]} = $row[$column_num];
  }
  return %results;
}
