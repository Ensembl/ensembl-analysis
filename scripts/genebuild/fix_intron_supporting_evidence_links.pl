#!/usr/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016] EMBL-European Bioinformatics Institute
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

use strict;
use warnings;

use Getopt::Long;

my $daf_outfile;
my $daf_infile;
my $sf_infile;
my $sf_outfile;


#85      29      25790   37912   38072   1       AADN04000776*1:37912:38072:1:canon      10.000  DEPTH   1

GetOptions(
            'outdaf=s'       => \$daf_outfile,
            'outsf=s'        => \$sf_outfile,
            'indaf=s'        => \$daf_infile,
            'insf=s'         => \$sf_infile,
           );

if (!$daf_outfile || !$sf_outfile ) {
  die("Need 2 outfiles: -outdaf and -outsf");
}
if (!$daf_infile || !$sf_infile ) {
  die("Need 2 infiles: -indaf and -insf");
}

print STDERR "DAFIN\t$daf_infile\nSFIN\t$sf_infile\nDAFOUT\t$daf_outfile\nSFOUT\t$sf_outfile\n";

my %to_replace;
my %to_keep;
my @new_lines;
my $count = 1;
my %seen;
open(ISE, $daf_infile) || die("Could not open $daf_infile");
while (<ISE>) {
  my @line = split("\t", $_);
  if (exists $seen{$line[6]}) {
    if ($seen{$line[6]} >= $line[7]) {
      $to_replace{$line[0]} = $to_keep{$line[6]};
#      print STDERR $line[0], ' -> ', $to_replace{$line[0]}, "\n";
    }
    else {
#      print STDERR $line[0], ' -> ';
      $line[0] = $to_replace{$line[0]} = $to_keep{$line[6]};
      $new_lines[$to_keep{$line[6]}-1] = \@line;
#      print STDERR $to_keep{$line[6]}, "\n";
    }
  }
  else {
#    print STDERR $line[0], ' -> ';
    $to_keep{$line[6]} = $line[0] = $to_replace{$line[0]} = $count++;
#    print STDERR $to_keep{$line[6]}, "\n";
    push(@new_lines, \@line);
  }
  $seen{$line[6]} = $line[7];
}
close(ISE) || die("Could not close file $daf_infile");
local $, = "\t";
open(WISE, '>'.$daf_outfile) || die("Could not open $daf_outfile");
foreach my $line (@new_lines) {
  print WISE @$line;
}
close(WISE) || die("Could not close $daf_outfile");
open(TISE, $sf_infile) || die("Could not open $sf_infile");
open(WTISE, '>'.$sf_outfile) || die("Could not open $sf_outfile");
local $\ = "\n";
while(<TISE>) {
  my @line = split;
  $line[1] = $to_replace{$line[1]};
  print WTISE @line;
}
close(TISE) || die("Could not close $sf_infile");
close(WTISE) || die("Could not close $sf_outfile");
