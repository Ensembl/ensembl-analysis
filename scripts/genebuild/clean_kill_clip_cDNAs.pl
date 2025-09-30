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


#!/usr/bin/env perl

=pod

=head1 NAME 

clean_kill_clip_cDNAs.pl

=head1 DESCRIPTION

=head1

This script performs all the necessary clean-ups and checks on a cDNA file to
prepare it for genomic alignment (e.g. using exonerate est2genome model). It
outputs a cleaned cDNA file and also keeps track of cDNAs discarded in the
process (how many of them, which ones, and why).

=head2 The basic functions performed are:


(1) Header cleaning, leaving behind the ">" fasta sign and the cDNA accession
    (including version number).

(2) Removal of cDNAs which are below a length threshold (default = 60bp)

(3) Removal of RefSeq predicted sequences (accession IDs starting with "XM"
    and "XR".

(4) Removal of duplicated cDNAs in the input file, ensuring there's only one
    entry per cDNA accession ID.

=head2 The optional functions performed are:


(1) Removal of cDNAs which are present in the Ensembl kill_list DB.

(2) Clipping of polyA tail of cDNAs which have passed all the above checks
    (i.e. not too short, not RefSeq predicted sequence etc)


=head1 OPTIONS

=head1

=head2 Mandatory options

        -cdnafile       The input cDNA file to be processed in FASTA format
        -outfile        The output file containing cleaned cDNAs, also in FASTA
                        format.

=head2 Optional

        -min_length     The minimum length threshold for cDNAs, set at 60bp
                        by default. The cDNA shouldn't be shorter than this
                        length before or after polyA clipping, or else it
                        will be discarded

        -clip           Boolean flag indicating whether you want polyA
                        clipping on the cDNAs. Set to "no" by default

        -use_kill_list  Boolean flag indicating whether you want to check
                        each input cDNA against the Ensembl kill_list DB.
                        if a cDNA is present in the kill_list DB, it will
                        be discarded. Set to "no" by default.

        -verbose        If you want wordy output to follow the script.
                        Switched off by default.


=head1 EXAMPLES

=head2 Clean headers, check cDNAs against kill_list, min_length = 60bp, don't clip polyA:

  perl clean_kill_clip_cDNAs.pl -cdnafile my_input_cDNAs.fa  -outfile cleaned_cDNAs.fa\
    -min_length 60 -use_kill_list  >& cleaning_cDNAs.log

=head2 Clean headers, check cDNAs against kill_list, clip polyA, min_length = 60bp:

  perl clean_kill_clip_cDNAs.pl -cdnafile my_input_cDNAs.fa  -outfile cleaned_cDNAs.fa\
    -min_length 60 -use_kill_list -clip >& cleaning_cDNAs.log

=cut


use strict;
use warnings;
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;
use Bio::EnsEMBL::Utils::PolyA;
use Bio::EnsEMBL::KillList::KillList;

$| = 1; # disable buffering

my $cdnafile;
my $seqoutfile;

my $min_length = 60;
my $clip       = 0;
my $kill       = 0;
my $verbose    = 0;

GetOptions( 'cdnafile=s'     => \$cdnafile,
            'outfile=s'      => \$seqoutfile,
            'min_length:n'   => \$min_length,
            'clip!'          => \$clip,
            'use_kill_list!' => \$kill,
            'verbose!'       => \$verbose, );

if (!$cdnafile || !$seqoutfile) {
  die print "You must provide both the input cDNA file using -cdnafile flag and ".
  "the output cDNA file using -outfile flag on the commandline.\n";
}

print STDOUT "Using minimum_length $min_length.\n" if ($verbose);

if ($kill) {
  print STDOUT "Will check cDNAs against kill_list.\n" if ($verbose);
}
if ($clip) {
  print STDOUT "Will clip polyA tails.\n" if ($verbose);
}

my $seqin = new Bio::SeqIO( -file   => "<$cdnafile",
                            -format => "Fasta", );

my $seqout = new Bio::SeqIO( -file   => ">$seqoutfile",
                             -format => "Fasta" );

my $total_cDNA = 0;                 # counter for total cDNAs processed

# More counters to keep track of cDNAs removed for various reasons

my $count_XM = 0;                   # predicted mDNAs
my $count_XR = 0;                   # predicted non-coding RNAs
my $count_too_short = 0;            # too short
my $count_too_short_after_clip = 0; # too short after polyA clipping
my $count_vanished_after_clip = 0;  # clipping causes the cDNA to disappear altogether
my $count_killed = 0;               # present in ensembl kill_list
my $written = 0;

my $kill_list_object = undef;
my %kill_list;

if ($kill) {
  $kill_list_object = Bio::EnsEMBL::KillList::KillList->new(-TYPE => "cDNA");
  %kill_list = %{ $kill_list_object->get_kill_list() };
}

my %hadhash;

SEQFETCH:
while ( my $cdna = $seqin->next_seq ) {
  $total_cDNA++;

  if ( $cdna->length < $min_length ) {
    if ($verbose) {
      print STDOUT $cdna->display_id . " is too short. Discarded.\n";
    }
    $count_too_short++;
    next;
  }

  my $display_id  = $cdna->display_id;
  my $description = $cdna->desc;

  # boring old fasta header
  if ( $display_id =~ /^(\S+\.\d+)$/ ) {
    $display_id = $1;
  # NCBI fasta headers
  } elsif ( $display_id =~ /gi\|\S+\|gb\|(\S+\.\d+)\|/ ) {
    $display_id = $1;
  } elsif ( $display_id =~ /gi\|\S+\|emb\|(\S+\.\d+)\|/ ) {
    $display_id = $1;
  } elsif ( $display_id =~ /gi\|\S+\|dbj\|(\S+\.\d+)\|/ ) {
    $display_id = $1;
  } elsif ( $display_id =~ /gi\|\S+\|sp\|(\S+)\|/ ) {
    $display_id = $1;
  } elsif ( $display_id =~ /gi\|\S+\|sp\|\|(\S+)/ ) {
    $display_id = $1;
  } elsif ( $display_id =~ /gi\|\S+\|ref\|(NM\S+\.\d+)\|/ ) {
    $display_id = $1;
  } elsif ( $display_id =~ /gi\|\S+\|ref\|(XM\S+\.\d+)\|/ ) {
    $count_XM++;
    next SEQFETCH;
  } elsif ( $display_id =~ /gi\|\S+\|ref\|(XR\S+\.\d+)\|/ ) {
    $count_XR++;
    next SEQFETCH;
  } elsif ( $display_id =~ /gi\|\S+\|tpe\|(\S+\.\d+)\|/ ) {
    $display_id = $1;
  } elsif ( $display_id =~ /gi\|\S+\|tpg\|(\S+\.\d+)\|/ ) {
    $display_id = $1;
  } elsif ( $display_id =~ /gi\|\S+\|prf\|\S*\|(\S+)/ ) {
    $display_id = $1;
  } elsif ( $display_id =~ /gi\|\S+\|pir\|\S*\|(\S+)/ ) {
    $display_id = $1;
  } elsif ( $display_id =~ /gi\|\S+\|pdb\|/ ) {
    next SEQFETCH;
  } elsif ( $display_id =~ /gi\|\S+\|(\S+)\|\S+\|/ ) {
    print STDERR "Unhandled subtype $1 for $display_id\n";
    next SEQFETCH;
  } elsif ($description) {
    my @labels = split /\s+/, $description;
    $display_id = $labels[0];
  } else {
    # Leave current display id
  }

  if ( exists( $hadhash{$display_id} ) ) {
    if ($verbose) {
      print STDOUT "Skipping $display_id as already had one with that name\n";
    }
    next SEQFETCH;
  }

  if ($kill) {
    my ($no_ver_id) = $display_id =~ /(\w+)\.\d/;
    if ( exists( $kill_list{$no_ver_id} ) ) {
      if ($verbose) {
        print STDOUT "$display_id is present in kill list DB, discarded.\n";
      }
      $count_killed++;
      next;
    }
  }

  $hadhash{$display_id} = 1;

  $cdna->display_id($display_id);
  $cdna->desc("");

  if ($@) {
    warn("can't parse sequence for [$description]:\n$@\n");
    next SEQFETCH;
  }

  my $polyA_clipper = Bio::EnsEMBL::Utils::PolyA->new();
  my $new_cdna;
  if ($clip) {
    $new_cdna = $polyA_clipper->clip($cdna);
    # where poly a clipping totally clips the entry the cdna is undef
    unless ( defined $new_cdna ) {
         print STDOUT "Clipping made a cdna vanish...\n";
         $count_vanished_after_clip++;
         next SEQFETCH;
    }
    if ( $new_cdna->length < $min_length ) {
      if ($verbose) {
        print STDOUT $new_cdna->display_id
            . " is under $min_length bp after clipping, discarded.\n";
      }
      $count_too_short_after_clip++;
      next SEQFETCH;
    }
  } else {
    $new_cdna = $cdna;
  }

  # write sequence
  $seqout->write_seq($new_cdna);
  $written ++
} ## end while ( my $cdna = $seqin...


print "\n";
printf STDOUT "Number of cDNAs processed: %50d\n", $total_cDNA;
printf STDOUT "Number of cDNAs written in output: %42d\n", $written;
printf STDOUT "Number of cDNAs shorter than $min_length bp: %41d\n", $count_too_short;
printf STDOUT "Number of XM_ predicted removed from the cDNA set: %26d\n", $count_XM;
printf STDOUT "Number of XR_ predicted non-coding RNAs removed from the cDNA set: %10d\n", $count_XR;
if ($clip) {
  printf STDOUT "Number of cDNAs shorter than $min_length bp "
              . "after polyA clipping: %20d\n", $count_too_short_after_clip;
  printf STDOUT "Number of cDNAs that vanished after clipping: %31d\n", $count_vanished_after_clip;
}

if ($kill) {
  printf STDOUT "Number of cDNAs found in the kill_list DB: %34d\n", $count_killed;
}
