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
my $accession = 'KQ759721.1';
my $rest_url = "http://rest.ensembl.org/sequence/region/chicken/$accession?content-type=text/x-fasta";
my @samstring = (
'ERR1298523.56710939/2	3	KQ759721.1	388	0	34M1330N66M	=	1	100	GGTGGGTTAGGGGTAGGGAAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGGTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGG	*	RG:Z:ERR1298523',
'ERR1298571.10022360/1	3	KQ759721.1	1460	0	26M3382N35M1N5M1N22M78N12M	=	1	100	GCTTAGGGTTAGGGTTAGGGTTACGGTAGGGTTAGGGTTAGGATTAGGGTTAGGGTTAGGGTAGGGTAGGGTTAGGGTTAGCGTTAGGGGTTAGGGTTAG	*	RG:Z:ERR1298571',
'ERR1298620.1954341/1	3	KQ759721.1	3741	0	6M1N25M465N16M629N8M1N35M1N10M	=	1	100	GTGTTAGGTTTAGGGTTAGGGTAGGGTTAGGGTTAGGGTTTGGGTAGGGTTAGGGTAGGGTTAGGGCTAGGGTTAGGGTTAGGGTTAGGGTAGGGTTAGG	*	RG:Z:ERR1298620',
'ERR1298619.3909349/1	3	KQ759721.1	3827	0	43M1069N49M1N8M	=	1	100	GGGTTAGGGTTAGGGTTAGGGTTAGGGTAGGGTTAGGGTTAGGTTAGCGTTAGGGTTAGGGTTAGGGTTTAGGGTTAGGGTTAGGGTTAGTGTAGGGTTA	*	RG:Z:ERR1298619',
'ERR1298538.12305296/2	3	KQ759721.1	3835	0	46M1087N32M1N22M	=	1	100	GTTAGGGTTAGGGTTAGGGTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTAGGTTTCAGGTTAGGGTTAGG	*	RG:Z:ERR1298538',
'ERR1298525.50972204/2	3	KQ759721.1	3857	0	29M860N61M10S	=	1	90	GGTTAGGGTTAGGGTTTGGGTTCGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGAGATCGGAAG	*	RG:Z:ERR1298525',
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
my $genomefile = catfile($root_dir, 'chicken_genome.fa');
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

my $bam_url = 'http://ftp.ensembl.org/pub/current_data_files/gallus_gallus/Gallus_gallus-5.0/rnaseq/Gallus_gallus-5.0.Roslin.merged.1.bam';
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
