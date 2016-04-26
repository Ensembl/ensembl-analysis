#!/usr/env perl

use strict;
use warnings;

use Getopt::Long;

use IO::Zlib qw(:gzip_external 1);
use IO::File;

my $fastqfile;
my $outfile;
my $outdir;
my $min_polyt = 10;
my $n_as_t = 1;
my $compression_level = 6;

&GetOptions (
            'i|infile=s' => \$fastqfile,
            'o|outfile=s' => \$outfile,
            'd|outdir=s' => \$outdir,
            'm|min_polyt=i' => \$min_polyt,
            'c|convert_N!' => \$n_as_t,
            'C|compression=i' => \$compression_level,
        );


if ($outdir) {
    if (-d $outdir) {
        $outfile = "$outdir/$outfile";
    }
    else {
        die("$outdir does not exists!\n");
    }
}
die("$outfile already exists!\n") if (-e $outfile);

my $use_zlib = $fastqfile =~ /gz$/;
my $unchanged = 0;
my $in_fh = $use_zlib ? IO::Zlib->new($fastqfile, 'rb') : IO::File->new($fastqfile, 'r');
my $out_fh = $use_zlib ? IO::Zlib->new($outfile, "wb$compression_level") : IO::File->new($outfile, 'w');
while (my $id = <$in_fh>) {
    my $seq = <$in_fh>;
    my $tmpseq = $seq;
    my $plusline = <$in_fh>;
    my $quality = <$in_fh>;

    die("Problem parsing $fastqfile\n") unless ($seq and $plusline and $quality);
    if ($n_as_t) {
        $tmpseq =~ s/.{2}[TN]{$min_polyt}[TN]*//;
    }
    else {
        $tmpseq =~ s/.{2}T{$min_polyt}T*//;
    }
    $unchanged += $seq eq $tmpseq;
    print $out_fh $id, $tmpseq, $plusline, substr($quality, (length($seq)-length($tmpseq)));
}
print STDOUT $unchanged, "\n";

$out_fh->close || die("Could not close $outfile\n");
$in_fh->close || die("Could not close $fastqfile\n");
