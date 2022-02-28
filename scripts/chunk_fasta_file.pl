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

#This is a script which will run through the fasta file you provide and
#chunk into into fasta files each containing the number of entries you
#specify. 
#the usage is
# chunk_fasta_file.pl fasta_file output_dir chunk_size

use warnings ;
use strict;


my $filename = shift;
my $output_dir = shift;
my $chunk_size = shift;

if(!$filename || !$output_dir || !$chunk_size){
  print "usage chunk_fasta_file.pl fasta_file output_dir chunk_size";
  exit;
}

if($filename eq '-h' || $filename eq '-help'){
  print "usage chunk_fasta_file.pl fasta_file output_dir chunk_size";
  exit;
}

&chunk_pepfile($filename, $output_dir, $chunk_size);

sub chunk_pepfile {
  my ($pepfile, $scratchdir, $size) = @_;
  
  #Chunk the peptide file
  open (PEPFILE, "$pepfile") or die "couldn't open $pepfile $!";
  my $count = 0;
  my $chunk = 1;
  #print STDERR "chunking peptide file\n";
  
  
  $/ = "\>";
  #print "have opened ".$pep_file."\n";
  while(<PEPFILE>){
    #print $_."\n";
    if ($_ ne "\>") {
      if ($count == 0) {
        open (CHUNK,">".$scratchdir."/".$pepfile."_chunk.$chunk") or die "couldn't open ".$scratchdir."/".$pepfile."_chunk.$chunk";
        #print "have opened ".$scratchdir."/chunks/chunk.$chunk\n";
      }
      
      $_ =~ s/\>$//;  
      
      print CHUNK ">$_";
      $count++;
      if ($count == $size) {
        $count = 0;
        $chunk++;
      }
    }
  }
  $/ = "\n";
}
