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

use warnings ;
use vars qw(%Config);
use Bio::EnsEMBL::Analysis::Config::Databases;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils ;
use strict;
use Getopt::Long;


my $core_db;
my $out_db;
my $write;
my $chromosome;


my $usage = "perl make_single_transcript_genes.pl
-core    $core_db, Database containing core gene set (use the key from Databases.pm)
-out     $out_db, Database to write the single transcript genes to  (use the key from Databases.pm)
-write   $write, Write the results";

$| = 1;

&GetOptions(
	    'core:s'     => \$core_db,
            'out:s'      => \$out_db,
	    'write!'     => \$write,
	    'chr:s'      => \$chromosome,
	   );

die($usage) unless ($core_db && $out_db);
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
my $core = $database_hash{$core_db};
my $out = $database_hash{$out_db};
throw("Db $core_db not found in Databases.pm\n") unless $core;
throw("Db $out_db not found in Databases.pm\n") unless $out;

my $sa = $core->get_SliceAdaptor;
my @new_genes;
my $count = 0 ;
CHR: foreach my $chr ( @{$sa->fetch_all('toplevel')}){
  if ( $chromosome ) {
    next CHR unless $chr->seq_region_name eq $chromosome;
    print "Running on " . $chr->name ."\n";
  }
  my @genes = @{$chr->get_all_Genes(undef, undef, 1)};
  foreach my $gene ( @genes ) {
    $count++;
  #  last CHR if $count > 100;
    foreach my $tran ( @{$gene->get_all_Transcripts} ){
      my $new_tran = clone_Transcript($tran);
      my ( $new_gene ) = @{convert_to_genes(($new_tran),$gene->analysis)};
      push @new_genes,$new_gene;
    }
  }
}
print "Had $count multi transcript genes now have " . scalar(@new_genes) ."\n";
# write the output


my $gene_adaptor = $out->get_GeneAdaptor;
GENE: foreach my $gene (@new_genes){
  unless ($write) {
    print "Not writing - write protect on\n";
    last GENE;
  }
  eval {
      $gene_adaptor->store($gene)
  };    
  if ($@){
    warning("Unable to store gene!!\n$@");
  }
}


