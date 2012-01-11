#!/usr/local/ensembl/bin/perl -w

# This is a script which will create the delimited text file containing
# the summary data describing the RNASeq reads and tissues. This text file is 
# required for running the RNASeq pipeline. The output format is just an example as
# it could vary depending on how your data looks.
# Please, have a look at ensembl-analysis/scripts/RNASeq/README for further information.

# The usage is:
# create_tissues_metadata_file.pl bam_files_path fastq_files_path description read_length paired output_file

# 'desc' will be the string for the DS (description) field
# 'bam_files_path', path to the directory where all of your bam files are located.
#       -Note that all of the bam files in this directory will be processed.
# 'fastq_files_path', path to the directory where all of your fastq files are located.
# 	-Note that fastq files will not be accessed. This parameter is only used for completing the metadata fields.
# 'read_length', assumed the same length for every read
# 'paired' can be either 0 (single end) or 1 (paired end).
#       -Note that this script only works correctly with paired end reads.
# 'output_file', path and filename where the text file contents will be written
#
# Example:
# perl $HOME/src/ensembl-analysis/scripts/RNASeq/create_tissues_metadata_file.pl $SCR/RNASeq/data/ $SCR/RNASeq/data/fastq/all/ "Canis Lupus familiaris" 101 1 $SCR/RNASeq/data/all_dog_tissues.txt


use strict;

my $bam_files_path = shift;
my $fastq_files_path = shift;
my $desc = shift;
my $read_length = shift;
my $paired = shift;
my $output_file = shift;

my $SAMTOOLS_PATH = "/software/solexa/bin/samtools";

if (!$bam_files_path | !$fastq_files_path || !$desc || !$read_length || !$paired || !$output_file){
  print "usage: create_tissues_metadata_file.pl bam_files_path fastq_files_path description read_length paired output_file";
  exit;
}

if ($bam_files_path eq '-h' || $bam_files_path eq '-help'){
  print "usage: create_tissues_metadata_file.pl bam_files_path fastq_files_path description read_length paired output_file";
  exit;
}

# open output file in write mode
open(OUTPUT_FILE, ">", $output_file) or die "cannot open ".$output_file;

# write headers in output file
# Please note that if your data looks in a different way, you'll may have to edit the setup_rnaseq_pipeline_config.pm file as well as the headers order below
print OUTPUT_FILE "ID\tPL\tPU\tLB\tST\tSM\tCN\tDS\tFILE\tLENGTH\tPAIRED\n";

# execute 'ls' command to get the bam filenames in the $bam_files_path directory given as parameter
my @bamfilenames;
open(LS, "ls $bam_files_path |");
while (<LS>) {
  chomp;
  if ($_ =~ m/.bam$/) {
    # only accepting the filenames whose end is equal to ".bam"
    push @bamfilenames, $_;
  }
}
close(LS);

# run samtools view for each bam file and parse their output to write the tab separated fields to the output file in the order found. As described above, in this case it is as follows:
# ID PL PU LB ST SM CN DS FILE LENGTH PAIRED
foreach my $filename (@bamfilenames) {
  open(SAMVIEW,$SAMTOOLS_PATH.' view -H '.$bam_files_path.$filename.' | grep @RG |');

  while (<SAMVIEW>) {
    my $platform_unit = "";
    my $current_line = "";

    chomp;
    my @bamfields = split(/\t/);
    
    foreach my $field (@bamfields) {
      (my $label,my $value) = split(':',$field);
      
      if ($label && $value) {
        $current_line .= "$value\t";

        if ($label eq "PU") {
          $platform_unit = $value;
        }
      }
    }

    # write 2 lines to the output file corresponding to both paired ends
    my $fastqfilename = $fastq_files_path.$platform_unit."_1.fastq";
    print OUTPUT_FILE "$current_line$desc\t$fastqfilename\t$read_length\t$paired\n";

    $fastqfilename = $fastq_files_path.$platform_unit."_2.fastq";
    print OUTPUT_FILE "$current_line$desc\t$fastqfilename\t$read_length\t$paired\n";
    
    $current_line = "";
  }
  close(SAMVIEW);
}

