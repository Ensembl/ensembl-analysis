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

use Test::More;
use File::Spec::Functions;

require LWP::UserAgent;

use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Analysis::Runnable::Sam2Bam;


#Generating the data needed
note('Preparing the data');
my $accession = 'NBDQ01000126.1';
my $rest_url = "http://rest.ensembl.org/sequence/region/fox/$accession?content-type=text/x-fasta";
my @samstring = (
'SRR5576245.15342458/1	3	NBDQ01000126.1	874344	0	104M162N11M10S	NBDQ01000126.1	1	115	CCCTTTTCTATATATTTTTAAGGTATTCATTTATCTGTTTCTAATTGGGGCACTTAATGTAGTCTTCTACCTGAAGGGTGAAACACTGCTATTTCTAGAAAAAAAAAAAAAAAAATATAACCCTT	*	RG:Z:SRR5576245_6',
'SRR5576245.16070957/1	3	NBDQ01000126.1	874346	0	103M161N16M6S	NBDQ01000126.1	1	119	CTTTTCTATATATTTTTAAGGTATTCATTTATCTGTTTCTAATTGGGGCACTTAATGTAGTCTTCTACCTGAAGGGTGAAACACTGCTATTTCTAGAAAAAAAAAAAAAAAAAAGAATACTTATA	*	RG:Z:SRR5576245_6',
'SRR5576247.27408061/1	3	NBDQ01000126.1	874336	0	105M165N15M5S	NBDQ01000126.1	1	120	AGCAGTTACCCTTTTCTATATATTTTTAAGGTATTCATTTATCTGTTTCTAATTGGGGCTCTTAATGTAGTCTTCTACCTGAAGGGTGTATCACTGCTATTTCTAGTAAAAAAAAAAAAAAAACA	*	RG:Z:SRR5576247_10',
'SRR5612849.8249029/2	3	NBDQ01000126.1	874336	0	112M162N13M	NBDQ01000126.1	1	125	AGCAGTTACCCTTTTCTATATATTTTTAAGGTATTCATTTATCTGCTTCTAATTGGGGCACTTAATGTAGTCTTCTACCTGAAGGGTGAAACACTGCTATTTCTAGAAAAAAAAAAAAAAAACGG	*	RG:Z:SRR5612849_3',
'SRR5612851.15542662/1	3	NBDQ01000126.1	874346	0	102M162N16M7S	NBDQ01000126.1	1	118	CTTTTCTATATATTTTTAAGGTATTCATTTATCTGTTTCTAATTGGGGCACTTAATGTAGTCTTCTACCTGAAGGGTGAAACACTGCTATTTCTAGAAAAAAAAAAAAAAAACCGAAATCCATTT	*	RG:Z:SRR5612851_6',
'SRR5612855.13813522/2	3	NBDQ01000126.1	874339	0	108M163N16M1S	NBDQ01000126.1	1	124	AGTTACCCTTTTCTATATATTTTTAAGGTATTCATTTATCTGTTTCTAATTGGGTCACTTAATGTAGTCTTCTACCTGAAGGGTGAAACACTGCTATTTCTAGAAAAAAAAAAAAAAAATGAGAC	*	RG:Z:SRR5612855_5'
);
my $root_dir = catdir('modules', 't');
my @samfiles = (catfile($root_dir, 'file1.sam'), catfile($root_dir, 'file2.sam'));
open(WH, '>'.$samfiles[0]) || die ('Could not open '.$samfiles[0]);
foreach my $line ((reverse @samstring)[0..3]) {
  print WH $line, "\n";
}
print WH '@EOF';
close(WH) || die('Could not close '.$samfiles[0]);

open(WH, '>'.$samfiles[1]) || die ('Could not open '.$samfiles[1]);
foreach my $line ((reverse @samstring)[4..5]) {
  print WH $line, "\n";
}
print WH '@EOF';
close(WH) || die('Could not close '.$samfiles[1]);

my $bamfile = catfile($root_dir, 'introns');
my $samtools = 'samtools';
my $headerfile = catfile($root_dir, 'merged_header.h');
my $genomefile = catfile($root_dir, 'fox_genome.fa');
my $ua = LWP::UserAgent->new();
$ua->env_proxy;
my $response = $ua->get($rest_url, ':content_file' => $genomefile.'.tmp');
die("Failed to fetch the genome sequence") unless ($response->is_success);
open(RH, $genomefile.'.tmp') || die("Could not open $genomefile.tmp");
open(WH, '>'.$genomefile) || die("Could not open $genomefile");
while (<RH>) {
  if (/^>/) {
    my $line = $_;
    $line =~ s/>.*/>$accession/;
    print WH $line;
  }
  else {
    print WH $_;
  }
}
close(RH) || die("Could not close $genomefile.tmp");
close(WH) || die("Could not close $genomefile");
unlink $genomefile.'.tmp';

my $bam_url = 'http://ftp.ensemblorg.ebi.ac.uk/pub/current_data_files/vulpes_vulpes/VulVul2.2/rnaseq/VulVul2.2.ENA.merged.1.bam';
my $command = "$samtools view -H $bam_url";
my @headers;
my $header_lcount = 0;
open(RH, $command.' | ') || die("Could not open command: '$command'");
open(WH, '>'.$headerfile) || die("Could not open header file $headerfile");
while (<RH>) {
  if (/^\@[RP]G/) {
    print WH $_;
    push(@headers, $_);
    ++$header_lcount;
  }
}
close(RH) || die("Could not close command: '$command'");
close(WH) || die("Could not close header file $headerfile");
# I need to substract one as samtools will add a line
--$header_lcount;

# Starting the tests
note('Starting the tests');
my $analysis = new_ok('Bio::EnsEMBL::Analysis', ['-logic_name', 'sam2bam_r_test'], 'Analysis');
my $runnable = new_ok('Bio::EnsEMBL::Analysis::Runnable::Sam2Bam', [
  '-analysis', $analysis,
  '-program', $samtools,
  '-samfiles', \@samfiles,
  '-genome', $genomefile,
  '-bamfile', $bamfile,
  '-header', $headerfile,
  ], 'Runnable'
  );

isa_ok($runnable->analysis, 'Bio::EnsEMBL::Analysis', 'Checking analysis');
ok(-e $runnable->program, 'Check program is executable');
like($runnable->program, qr/$samtools/, 'Check program is samtools');
cmp_ok(scalar(@{$runnable->samfiles}), '==', 2, 'Check samfiles');
is($runnable->genome, $genomefile, 'Check genome');
is($runnable->headerfile, $headerfile, 'Check headerfile');
is($runnable->bamfile, $bamfile, 'Check bamfile');
$runnable->run;
ok(-e $genomefile.'.fai', 'Genome index file exists');
ok(-e $bamfile.'.bam', 'BAM file exists');
ok(-e $bamfile.'.bam.bai', 'Index file exists');
ok(!-e $bamfile.'_unsorted.bam', 'Unsorted BAM deleted');
my $program = $runnable->program;
my $result_check = `$program view $bamfile.bam`;
is($result_check, join("\n", @samstring)."\n", 'Checking sort result');
$command = "$program view -H $bamfile.bam";
my @results;
open(RH, $command.' | ') || die("Could not open command: '$command'");
while (<RH>) {
  if (/^\@[RP]G/) {
    push(@results, $_);
  }
}
close(RH) || die("Could not close command: '$command'");
is_deeply(\@headers, [(@results)[0..$header_lcount]], 'Checking header of sort result');
done_testing();

#Cleaning
note('Cleaning');
foreach my $file ($genomefile, $genomefile.'.fai', @samfiles, $bamfile.'.bam', $bamfile.'.bam.bai', $bamfile.'.header', $headerfile) {
  unlink $file if (-e $file);
}
