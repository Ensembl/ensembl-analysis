#!/usr/bin/env perl

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

use warnings ;
use vars qw(%Config);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use strict;
use Getopt::Long;

my $fastq_dir;
my $size;
my $meta_file;
my $fc ;
my @files;
my $chunks;
my $update;
my $new_meta;
my %meta_hash;
my $check;
my $header;

my $usage = "perl split_fastq.pl
-dir        $fastq_dir, Directory containg fastq files
-chunks_dir $chunks, Directory to put the fastq chunks in
-size       $size, number of fastq to have in each file
-meta       $meta_file, file containing the metadata about the RNASeq - this will get expanded to include the chunked fastq,
-file_col   $fc, Column in the meta data file that contains the fastq file names,
-update     $update, Once the chunking is complete run the script again with - update to check the chunks and update the meta file,
-check      $check, Count the reads in the chunks to check they are the same as the original file, (optional)
";
$| = 1;

&GetOptions(
  'dir:s'        => \$fastq_dir,
  'chunks_dir:s' => \$chunks,
  'size:s'       => \$size,
  'meta:s'       => \$meta_file,
  'file_col:s'   => \$fc,
  'update!'      => \$update,
  'check!'       => \$check,
);

die($usage) unless ($fastq_dir && $size && $meta_file

);
open (META,"$meta_file") or die("Cannot open meta file $meta_file\n");
my $line = 0;
while ( <META> ) {
  chomp;
  $line++;
  if ( $line == 1 ) {
    $header = "$_";
    next;
  }
  my @cells = split(/\t/,$_);
  push @files , $cells[$fc-1];
  $meta_hash{$cells[$fc-1]} = \@cells;
}

unless ( $update ) {
  print "Have these fastq:\n";
  foreach my $file (@files) {
    print "$file\n";
  }
  print "splitting into chunks containing $size reads\n";
  print "Please run the bsubs to split the files\n";
  foreach my $file (@files ) {
    throw("File $file not found\n")
    unless -e $fastq_dir."/".$file;
    print "bsub -o chunk_$file.%J.out -e chunk_$file.%J.err \"split $fastq_dir/$file -l " . $size*4 ." $chunks/$file-\"\n";  
  }
  print "Once the bsubs have finished run the script again with the -update flag to check that the chunking has worked and to update the meta file with the new fastq files\n";
  exit;
}
# check the chunks are correct
foreach my $file (@files ) {
  my $total_reads = 0;
  my $chunked_reads = 0;
  my @cells = @{$meta_hash{$file}};
  # count how many reads we have in the fasta file to start with
  my $command = "wc -l $fastq_dir/$file";
  print "File $file " ;
  $total_reads = count_reads($fastq_dir."/" . $file) if $check;
  $command = "ls $chunks/$file-*";
  open  ( my $fh,"$command 2>&1 |" ) || 
  throw("Error counting chunks");
  print "$file has the following chunks:\n";
  while (<$fh>){
    chomp;
    print "$_\n";
    $chunked_reads += count_reads($_) if $check;
    #write the new meta file
    for ( my $i =0 ; $i < scalar(@cells) ; $i++ ) {
      if ( $i == $fc -1  ) {
        my @path = split("\/",$_);
        my $file = pop(@path);
        print STDERR "FILE: $file\n";
        $new_meta .= "$file\t";
        next;
      }
      $new_meta .= $cells[$i] ."\t";
    }
    $new_meta .= "\n";
  }
  if ( $check ) {
    if ( $total_reads == $chunked_reads ) {
      print "all reads accounted for $total_reads in file and $chunked_reads in chunks\n";
    } else {
      throw("Reads in file and chunks do not add up $total_reads in $file and $chunked_reads in chunks\n");
    }
  }
}
print "Writing new meta file $meta_file.chunks \n";
open (METACHUNK,">$meta_file.chunks") or die ( "Cannot open new meta file for writing $meta_file.chunks\n");
print METACHUNK "$header\n";
print METACHUNK "$new_meta";
close METACHUNK;

exit;

sub count_reads {
  my ($file) = @_;
  my $total_reads = 0;
  my $command = "wc -l $file";
  eval  {
    open  ( my $fh,"$command 2>&1 |" ) || 
    throw("Error counting reads");
    while (<$fh>){
      chomp;
      if ( $_ =~ /(\d+) \S+/ ) {
        $total_reads = $1 / 4;
        print  " has $total_reads reads";
      }
    }; if($@){
      throw("Error processing alignment \n$@\n");
    }
  };
  print "\n";
  return $total_reads;
}

