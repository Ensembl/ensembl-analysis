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
# limitations under the License.

#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $fasta_file;
my $map_file;
my $edited_fasta_file;

&GetOptions(
            'fasta_file:s' => \$fasta_file,
            'map_file:s' => \$map_file,
            'edited_fasta_file:s' => \$edited_fasta_file
            );

my %map_ids;

open (EFF,">".$edited_fasta_file) || die "Could not open edited_fasta_file for writing\n";
open (MAP, "$map_file") or die "Can't open ".$map_file."\n";

while(<MAP>){
  chomp;
  my @values = split(/\t/,$_);
 
  if ($values[1] eq "EMBL-CDS"){
   $map_ids{$values[2]} = $values[0];
  }
}

close MAP;

open (FASTA, "$fasta_file") or die "Can't open ".$fasta_file."\n";

while (<FASTA>){
  chomp;

	if ($_ =~/^>/){
	 my @accessions = split(/\s+/,$_);
	  if($map_ids{$accessions[1]}){
      print EFF $accessions[0]." ".$map_ids{$accessions[1]}."\n";
	  }else{
		  print EFF $_."\n";
    }
	}else{
	  print EFF $_."\n";
	}
}
close EFF;
close FASTA;
