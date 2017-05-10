#!/usr/bin/env perl

# This script creates bigwig files from all bam files in a given directory
# parameters can be stated on the command line or they can be filled in below
# usage:
# $ perl bam2bw_converter.pl
# or
# $ perl bam2bw_converter.pl -qdir /lustre/scratch110/ensembl/dm15/molly/PoeFor/output/rnaseq/bamtoBW2/ -outdir /lustre/scratch110/ensembl/dm15/molly/PoeFor/output/rnaseq/bamtoBW2/
# or you can read the parameters from a file and execute using awk if you're converting files from different directories:
# $ awk '{system ("perl bam2bw_converter.pl -qdir "$1 "-outdir "$2)}' filename.txt 

use strict;
use warnings;

use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

# specify the query and output directories
my $qdir = '/nfs/ensnfs-dev/staging/SPECIES/ASSEMBLY_NAME/rnaseq/';
my $outdir = '/nfs/ensnfs-dev/staging/SPECIES/ASSEMBLY_NAME/rnaseq/';

my $bedtools = '/nfs/software/ensembl/RHEL7/linuxbrew/bin/bedtools';
my $bed2bigwig = '/nfs/software/ensembl/RHEL7/linuxbrew/bin/bedGraphToBigWig';
# set memory requirements for tissue and merged files
my $mergebG_mem = 8000;
my $tisbG_mem = 3000;

&GetOptions (
  'qdir=s'          => \$qdir,
  'outdir=s'        => \$outdir,
  'bedtools=s'      => \$bedtools,
  'bed2bigwig=s'    => \$bed2bigwig,
  'merge_mem=n'     => \$mergebG_mem,
  'tissue_mem=n'    => \$tisbG_mem,
);

# make sure the directories stated above end with a backslash
unless ($qdir=~/\/$/){
  $qdir = $qdir . '/';
}
unless ($outdir=~/\/$/){
  $outdir = $outdir . '/';
}

die ("$qdir does not exist!\n") unless (-d $qdir);
die ("$outdir does not exist!\n") unless (-d $outdir);
# Now use bedtools and bedGraphToBigWig to convert the bam files
opendir (DIR, $qdir);
my @dir = readdir DIR;
foreach (@dir) {
  if ($_ =~ /bam$/) {
    my $mem_bG;
    my $mem_bw;
    my $job_name = $_;
    $job_name =~ s/\.bam//;

    # get all the chromosome lengths
    # system should return 0 if the command was successful
    if (system ("samtools view -H ". $qdir. $_ . " | grep '\@SQ'  | sed 's/SN://' | sed 's/LN://' | awk -v OFS='\t' '{ print \$2, \$3}' > $outdir/chrom_lengths.txt")) {
        die("samtools for getting the chromosome length died: $?\n");
    }

    # increase the memory for merge files
    if ($_ =~ /merge/) {
      $mem_bG = $mergebG_mem;
    }
    else {
      $mem_bG = $tisbG_mem;
    }
    # produce the bedGraph file, which will then be converted to BigWig
    if (system ('bsub -oo ' . $outdir . $job_name . '_genomecov.out -M ' . $mem_bG . ' -R "select[mem>' . $mem_bG . '] rusage[mem=' . $mem_bG . ']" -J ' . $job_name . ' "'.$bedtools.' genomecov -ibam ' . $qdir . $_ . ' -bg -split > ' . $outdir . $job_name . '.bedGraph"')) {
        die("Could not submit job for running bedtools\n");
    }

    # from the results of bedtools create a bigwig file
    if (system ('bsub -oo ' . $outdir . $job_name . '_bedgraph2bw.out -w "done(' . $job_name . ')" -M 3000 -R "select[mem>3000] rusage[mem=3000]" "'.$bed2bigwig.' ' . $outdir . $job_name . '.bedGraph ' . $outdir . 'chrom_lengths.txt ' . $outdir . $job_name . '.bam.bw"')) {
        die("Could not submit job for running Bed2BigWig\n");
    }
  }
}

exit;



