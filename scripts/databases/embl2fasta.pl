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

use warnings;
use strict;
use feature 'say';


my $infile = $ARGV[0];
my $outfile = $infile;
$outfile =~ s/\.dat/\.fasta/;

my $state = 0;
my $current_accession = "";
my $current_version = "";
my $current_seq = "";
open(IN,$infile);
open(OUT,">$outfile");
while(<IN>) {
  my $line = $_;
  if($state == 0 && $line =~ /^AC +([a-zA-Z\d]+)\;/) {
    $current_accession = $1;
    $state = 1;
  }

  elsif($state == 1 && $line =~ /^DT.+sequence version (\d+)/) {
    $current_version = $1;
    $state = 2;
  }

  elsif($state == 1 && $line =~ /^SQ +SEQUENCE +/) {
    die "Failed to parse version for ".$current_accession."\nExiting";
  }

  elsif($state == 2 && $line =~ /^SQ +SEQUENCE +/) {
    $state = 3;
  }

  elsif($state == 3 && !($line =~ /^\/\//)) {
    $current_seq .= $line;
  }

  elsif($state == 3 && $line =~ /^\/\//) {
    $current_seq =~ s/ //g;
    print OUT '>'.$current_accession.".".$current_version."\n".$current_seq;
    $state = 0;
    $current_accession = '';
    $current_version = '';
    $current_seq = '';
  }

}
close OUT;
close IN;

exit;
