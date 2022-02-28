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

=head2

  This script adds an alt allele group to the given database for each set of gene stable IDs in each line of the file 'file'.

=cut

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;

my $host;
my $dbname;
my $port;
my $user;
my $pass;
my $dnahost;
my $dnauser;
my $dnadbname;
my $dnaport;
my $store;
my $file;
my $path;

$| = 1;

&GetOptions(
  'host:s'      => \$host,
  'user:s'      => \$user,
  'pass:s'      => \$pass,
  'dbname:s'    => \$dbname,
  'path:s'      => \$path,
  'port:n'      => \$port,
  'dnahost:s'   => \$dnahost,
  'dnauser:s'   => \$dnauser,
  'dnadbname:s' => \$dnadbname,
  'dnaport:n'   => \$dnaport,
  'file:s'      => \$file,
  'store'       => \$store,
);

if (!defined($file)) {
  die "Must give file argument specifying file containing gene sets to link as alt alleles. An alt allele group per line with gene stable IDs separated by space.\n";
}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host   => $host,
  -user   => $user,
  -port   => $port,
  -pass   => $pass,
  -dbname => $dbname
);

my $dnadb;
if ($dnadbname) {
  $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host   => $dnahost,
    -user   => $dnauser,
    -port   => $dnaport,
    -dbname => $dnadbname
  );
  $db->dnadb($dnadb);
}

my $aaga = $db->get_AltAlleleGroupAdaptor();
my $ga = $db->get_geneAdaptor();

$| = 1;

open FP,"<$file" or die "Couldn't open file $file\n";

while (<FP>) {
  chomp;
  my (@geneids) = split();

  my @genes;
  foreach my $geneid (@geneids) {
    my $gene = undef; 
    eval {
      $gene = $ga->fetch_by_stable_id($geneid);
    };
    if ($@) {
      die "Failed to fetch $geneid at line $_";
    }
    if (!defined($gene)) {
      die "Failed to fetch $geneid at line $_";
    }
    push @genes,$gene;
  }
  
  if ($store) {
    print "Storing alt alleles for genes $_\n";
    
    my %flags = ('MANUALLY_ASSIGNED' => '1');
    
    my $first_gene = shift(@genes);
    my $new_group = Bio::EnsEMBL::AltAlleleGroup->new(-MEMBERS => [[$first_gene->dbID(),\%flags] ]);
    
    foreach my $gene (@genes) {
      $new_group->add_member($gene->dbID(),\%flags);
    }
    
    my $new_group_id = $aaga->store($new_group);
    print "Stored alt allele group ".$new_group_id."\n";
  } else {
    print "Not storing alt alleles for genes $_ because -store not set\n";
  }
}
