# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

#!/usr/bin/env perl

use strict; 
use warnings;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

$| = 1; # disable buffering

my $infile;
my $blastfile;
my $threshold = 0.0001;
my $seqhash;
my $regex;

GetOptions( 
	    'infile:s'      => \$infile,
	    'blastfile:s'   => \$blastfile,
	    'threshold:s'   => \$threshold,
	    'regex:s'       => \$regex,
	   );

# usage
if(!defined $infile    ||
   !defined $blastfile   )
{
  print  "USAGE: fasta_blast_threshold.pl 
-infile     $infile, input fasta file
-blastfile  $blastfile , blast output file of the input vs a database
-threshold  $threshold , sequences scoring above and below the threshold ( top blast hit ) are written to separate files
-regex      $regex, Use a regex on the description to push sequences to the above file\n";
    exit(1);
}

my $seqin  = new Bio::SeqIO(-file   => "<$infile",
			    -format => 'fasta'
			   );
my $seqouta  = new Bio::SeqIO(-file   => ">$infile"."_keep",
			      -format => 'fasta'
			     );
my $seqoutb  = new Bio::SeqIO(-file   => ">$infile" ."_discard",
			      -format => 'fasta'
			     );
die ("cannot open files for writing\n")  unless ( $seqouta && $seqoutb ) ;

while ( my $seq = $seqin->next_seq ) {
  $seqhash->{$seq->display_id} = $seq;
}

my $seq;
my $score;
my $string;
my $go = 0;
open (BLAST,"$blastfile") or die ( "Cannot open $blastfile \n");
while (<BLAST>) 
{
  chomp;
  if ( $_ =~ /^EXIT CODE/ ) 
  {
    my $bioseq = $seqhash->{$seq};
    $bioseq->description($bioseq->description . " Blast hit $string") if $string;
    die("Cannot find sequnece $seq from hash\n")
      unless $bioseq;
    if ( $regex && $string =~ /$regex/ ) 
    {
      $score = 1;
    }
    # dont let though known repeats
    unless ( $seq =~ /Unknown/ ) 
    {
      $score = 1;
    }
    
    if ( $score && $score < $threshold ) 
    {
      $seqoutb->write_seq($bioseq);
    } 
    else 
    {
      $seqouta->write_seq($bioseq);
    }
    $seq = undef;
    $score = undef;
    $string = undef;
  }

  if ( $_ =~ /Query=  (\S+).+/)
  {
    $seq = $1;
    if ( $seqhash->{$1} ) 
    {}
  }

  if ( $_ =~ /Sequences producing High-scoring Segment Pairs/ )
  {
    $go = 1;
  }
  if ( $go )
  {
    my @array = split(/\s+/,$_);
    if (scalar(@array) > 3 &&  $array[scalar(@array)-1] =~ /\d+/ && 
        $array[scalar(@array)-2] =~ /\d+/ && 
        $array[scalar(@array)-3] =~ /\d+/) 
    {
      $string = $_;
      $go = 0;
      $score = $array[scalar(@array)-2];
    }
  }
}

