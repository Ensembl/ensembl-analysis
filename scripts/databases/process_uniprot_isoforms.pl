#!/usr/env perl

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
