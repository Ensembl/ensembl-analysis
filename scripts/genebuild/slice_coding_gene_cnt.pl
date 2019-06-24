#!/usr/bin/env perl
#
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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
# limitations under the License

#This script checks slices for a given db to flag up cases where slices > 5mb have no protein coding gene

#

use strict;
use warnings;

use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use feature 'say';
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

my $dbname = '';
my $host = '';
my $user = '';
my $port = '';
my $pass = '';
my $driver = '';

GetOptions('dbname:s' => \$dbname,
           'host:s'  => \$host,
           'user:s' => \$user,
           'port:s' => \$port,
           'pass:s' => \$pass,
           'driver:s' => \$driver,
           );
my $slice_cnt = 0; my $gene_cnt = 0; my $slice_with_gene = 0; my $slice_no_gene = 0; my $size = 0;

my $db_adaptor = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -dbname => $dbname,
  -host   => $host,
  -port   => $port,
  -user   => $user,
  -pass   => $pass,
  -driver => $driver,
);
my $slice_adaptor = $db_adaptor->get_SliceAdaptor();
foreach my $slice ( @{ $slice_adaptor->fetch_all('toplevel') } ){
	
	my $gene_cnt = 0;
	$slice_cnt++;
		
	#retrieving gene from slice 
	foreach my $gene ( @{ $slice->get_all_Genes } ){
		if ($gene->biotype eq 'protein_coding'){
			#counting all protein coding genes
			$gene_cnt++;
		}
		
	}
	if ($slice->length >= 5000000) {#check that slice is bigger than 5mb
		
		if ($gene_cnt < 1){#check if slice contains protein coding genes
			 throw("slice has no protein coding gene");
			#say "slice " . $slice->name . " has no protein coding gene";
			$slice_no_gene++;
		}
		else{#slice contains protein coding genes
			#say "slice " . $slice->name . " has $gene_cnt protein coding genes";
			$slice_with_gene++;
		}
	}
	else{#count number of slice in database
		$size++;
	}
	
}
#print stats found
say "Total slice = $slice_cnt"; say "Slice with genes = $slice_with_gene"; say "slice with no gene = $slice_no_gene";
say "slice less 5mb = $size";
