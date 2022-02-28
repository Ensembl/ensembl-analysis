#!/usr/bin/env perl

# Copyright [2019-2022] EMBL-European Bioinformatics Institute
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
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw( throw warning verbose);

my $input_file;
my $output_file;
my $conversion_type;
my $remove_masking = 0;
GetOptions( 'input_file:s'      => \$input_file,
            'output_file:s'     => \$output_file,
            'conversion_type:s' => \$conversion_type,
            'remove_masking!'   => \$remove_masking);

unless($input_file && $output_file && $conversion_type) {
  throw("You must specify both an input file, an output file and a conversion type");
}

unless(-e $input_file) {
  throw("The input file specified does not exist. Input file: ".$input_file);
}

if($conversion_type eq "slice_name_to_seq_region_name") {
  slice_name_to_seq_region_name($input_file,$output_file,$remove_masking);
} else {
  throw("The conversion type you selected is not supported. Conversion type selected: ".$conversion_type);
}

exit;

sub slice_name_to_seq_region_name {
  my ($input_file,$output_file,$remove_masking) = @_;

  open(IN,$input_file);
  unless(open(OUT,">".$output_file)) {
    throw("Could not open output file for writing. Output file: ".$output_file);
  }

  while(<IN>) {
    my $line = $_;
    if($line =~ /^>/) {
      unless($line =~ /[^\:]+\:[^\:]+\:([^\:]+)\:/) {
        throw("Failed to parse the header line. Expected to find a seq region name after the second colon in header. Header used: ".$line);
      }
      my $header = ">".$1;
      say OUT $header;
    } else {
      if($remove_masking) {
        $line = uc($line);
      }
      print OUT $line;
    }
  }
  close OUT;
  close IN;
}
