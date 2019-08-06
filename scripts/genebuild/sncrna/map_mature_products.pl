#!/usr/env perl
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

use strict;
use warnings;

use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;

my $min_size = 18;
my $max_size = 25;

my ($mature_fa, $precursor_fa, $working_dir) = @ARGV;

GetOptions( 'maturefa|m:s' => \$mature_fa,
            'precursors|p:s' => \$precursor_fa,
            'output_dir|o:s' => \$working_dir);


# filter and collapse mature sequences
my %mature;
my $seqio = Bio::SeqIO->new(-file => $mature_fa, '-format' => 'Fasta');

while(my $seq = $seqio->next_seq){
 if( $seq->length ge $min_size && $seq->length le $max_size){
   $mature{$seq->seq} = $seq->primary_id;
  }
}

my $fn = $working_dir . "/mature_mirnas.fa";
open(FH, '>', $fn) or die "Could not write to $fn";

foreach my $seq (keys %mature) {
  print FH ">" .
    $mature{$seq} . "\n" .
    $seq . "\n";
}
close(FH);

# build bowtie index
my $command = "bowtie-build $precursor_fa $working_dir/precursors";

print STDERR "Building bowtie index for pre-miRNAs >>> $command\n";

if (system($command)){
  print("Error building bowtie index\nError code: $?\n");
}

# map mature products to precursors
my $options = " -p16 -v 3 -M 20 -k 5 --norc --best -f --sam --sam-nohead ";

$command = "bowtie $options $working_dir/precursors $working_dir/mature_mirnas.fa $working_dir/mature_mirnas_unfiltered.sam";

print STDERR "Mapping miRBase mature miRNAs to pre-cursors >>> $command\n";

if (system($command)){
    print STDERR "Error mapping mature miRNAs w/ Bowtie\nError code: $?\n";
}

# filter alignments

$fn = "$working_dir/mature_mirnas_unfiltered.sam";
my $fo = "$working_dir/mature_mirnas.sam";

open(FH, '<', $fn) or die "Could not open $fn";
open(FO, '>', $fo) or die "Could not open $fo for writing";
while(<FH>){
  my @contents = split /\t/, $_;

  if($contents[1] ne '4'){
    print FO $_;
  }
}

close($fn);
close($fo);

