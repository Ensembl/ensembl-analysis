#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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
use strict;
use Getopt::Long;
use Bio::SeqIO;

my ($input_pep_file,
    $output_chunk_dir,
    $chunk_size);

GetOptions('pepfile=s'   => \$input_pep_file,
           'chunkdir=s' => \$output_chunk_dir,
           'chunksize=s' => \$chunk_size);


die "You must supply a valid input peptide file\n"
    if not defined $input_pep_file or not -e $input_pep_file;
die "You must supply a valid output chunk directory\n"
    if not defined $output_chunk_dir or not -d $output_chunk_dir;

if (not defined $chunk_size or $chunk_size < 0) {
  warn "No/invalid chunk size given; defaulting to 20";
  $chunk_size = 20;
}

my $seqio = Bio::SeqIO->new(-format => 'fasta',
                            -file   => $input_pep_file);

my $count = 0;
my $chunk_num = 1;
my $outseqio;

while (my $seq = $seqio->next_seq) {
  if (not defined $outseqio) {
    $outseqio = Bio::SeqIO->new(-format => 'fasta',
                                -file   => ">$output_chunk_dir/chunk." . $chunk_num++);
  }

  $outseqio->write_seq($seq);
  $count++;

  if ($count >= $chunk_size) {
    $outseqio->close;
    $outseqio = undef;
    $count = 0;
  }
}
