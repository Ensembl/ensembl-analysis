#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

use strict;
use warnings;

use Getopt::Long;
use Bio::SeqIO;

my $infile;
my $outfile;
my $isofile;
my $use_version = 0;
my $use_description = 0;

&GetOptions (
            'i|infile=s'  => \$infile,
            'v|isofile=s' => \$isofile,
            'o|outfile=s' => \$outfile,
            'version!'    => \$use_version,
            'desc!'       => \$use_description,
        );

my %accessions;
open(RF, $infile) || die("Could not open $infile\n");
while (<RF>) {
    if (/^>(\w+)\.(\d+)/) {
        $accessions{$1} = $2;
    }
}
close(RF) || die("Could not close $infile\n");


my $sequences = Bio::SeqIO->new(-format => 'fasta', -file => $isofile);
my $writer = Bio::SeqIO->new(-format => 'fasta', -file => '>'.$outfile);
$writer->preferred_id_type('accession');
$writer->preferred_id_type('accession.version') if ($use_version);
while (my $seq = $sequences->next_seq()) {
    my ($accession, $isoform_id) = $seq->id =~ /[sptr]{2}\|(\w+)(-\d+)/;
    if (exists $accessions{$accession}) {
        $seq->accession_number($accession.$isoform_id);
        $seq->version($accessions{$accession}) if ($use_version);
        $seq->desc('') unless ($use_description);
        $writer->write_seq($seq);
    }
}
