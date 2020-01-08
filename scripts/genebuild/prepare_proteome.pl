#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

# call 
#
# perl ./genebuild/prepare_proteome.pl 
#  -proteome_file proteome_output_file.fa 
#  -file_info protein_input_file.fa='(\S+)' 
#
# use regex '/(\S+)/' if header format of your fasta is 
# >ABC234.1
# >ABC234.1
#

use strict;
use warnings;
use Bio::EnsEMBL::KillList::KillList;
use Bio::EnsEMBL::KillList::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Tools::Logger qw (logger_info logger_verbosity);
use Getopt::Long;

my $output_file;
my @file_info;
my $skip_xs = 5;
my $use_killlist = 1;
my $skip_xps = 1;
my $logger_verbosity = "NONE";
my $min_length = 15;

&GetOptions( "proteome_file|pmatch_output_file:s"    => \$output_file,
             "file_info:s@"       => \@file_info,
             "skip_xs:s"          => \$skip_xs,
             "use_killlist!"      => \$use_killlist,
             "min_length:i"       => \$min_length,
             "logger_verbosity:s" => \$logger_verbosity, );

logger_verbosity($logger_verbosity);
my $kill_list_object = Bio::EnsEMBL::KillList::KillList->new(-TYPE => 'PROTEIN');
my %kill_list = %{$kill_list_object->get_kill_list()};

my %files;
my %refseq;
my %uniprot;

foreach my $file_info (@file_info) {
  #print "Have file info ".$file_info."\n";
  my ( $file, $regex, $other ) = split /\=/, $file_info;
  #print "HAVE FILE ".$file." and regex ".$regex."\n";
  throw( "Need both a file " . $file . " and a regex " . $regex )
    if ( !$file || !$regex );
  $files{$file} = $regex;
  if ( $other && $other =~ /^refseq$/i ) {
    $refseq{$file} = 1;
  } elsif ( $other && $other =~ /^uniprot$/i ) {
    $uniprot{$file} = 1;
  }
}


if(-e $output_file){
  print "\nProtein file '" . $output_file . "' already exists, these ".
    "entries will be appended to the end of the file.\n";
  print "Do you want this? Answer y/n or d to delete the file: ";
  my $reply = get_input_arg();
  if($reply =~ m/n/i) {
    print "\nYou must delete or rename '" . $output_file . "' before you continue.\n";
    exit;
  } elsif($reply =~ m/d/i) {
    print "\nRunning rm -f '$output_file'\n";
    system("rm -f $output_file");
    print "File '$output_file' deleted.\n";
  } else {
    print "\nWill append new entries to '" . $output_file . "'.\n";
  }
}

my $output_io = Bio::SeqIO->new( -format => 'fasta',
                                 -file   => ">>" . $output_file, );

my %ids;
foreach my $file(keys(%files)) {
  my $refseq = $refseq{$file};
  my $uniprot = $uniprot{$file};
  my $regex = $files{$file};
  print "\nFile '" . $file . "' with regex: '" . $regex . "' was specified.\n\n";
  my $x_string;
  if($skip_xs){
    $x_string = "X" x $skip_xs;
  }

  my $io = Bio::SeqIO->new( -format => 'fasta',
                            -file   => $file, );

 SEQ:while(my $seq = $io->next_seq){  

    my $parseable_string = $seq->id; 
    $parseable_string.=" ".$seq->desc if $seq->desc() ;  

    # bit of a hack:
    my $sequence_version;
    if ($parseable_string =~ /SV=(\d+)$/) {
      $sequence_version = $1;
    }

    my ($id) = $parseable_string =~ /$regex/; 

    if($id =~ m/^1$/){
      print "got ".$id." from ".$parseable_string."\n";
      die;
    } 
    if(!$id){
      warn($regex." failed to parse an id out of ".
            $parseable_string);
        next SEQ;
    }
    if($x_string){
      if($seq->seq =~ /$x_string/){
        logger_info($id." contains too many X characters, SKIPPING");
        next SEQ;
      }
    }
    if (length($seq->seq) < $min_length) {
      logger_info($id." is shorter than $min_length");
      next SEQ;
    }
    if($refseq && !($id =~ /^NP/)){
      logger_info($id." isn't an NP so skipping");
      next SEQ;
    }
    if($uniprot && !($id =~ /^sp/ || $id =~ /^tr/)){
      logger_info($id." isn't a sp or tr ");
      throw("Don't know what to do with uniprot id $id");
    } elsif ($uniprot && $regex eq '(\S+)') {
      # need to fix the header because it will look like
      #sp|Q9I9D5|RX1_ASTFA
      #tr|O42292|O42292_ASTMX
      $id =~ /^(sp|tr)\|(\w+)\|\w+$/;
      $id = $2;
      if ($id !~ /\.\d+/ && defined $sequence_version) {
        $id .= ".".$sequence_version;
      }
      if ($id !~ /\w+/) {
        throw("uniprot id looks wrong: $id");
      }
    }
    if ($seq->seq =~/(B|Z|J)/) {
      warn("AMBIGUITY CODE FOUND!!! $id contains at least one residue of code $1. May interfere with Pmatch.");
    }
    my $no_version_id = $id;
    $no_version_id =~ s/\.\d+//;
    if($use_killlist && exists $kill_list{$no_version_id}){
      logger_info($id." on kill list as $no_version_id");
      next SEQ;
    }
    if(exists($ids{$id})){
      logger_info($id." has already been stored");
      next SEQ;
    }
    $ids{$id} = 1;
    $seq->desc("");
    $seq->id($id);
    my $seq_string = $seq->seq;
    $seq_string =~	s/U/X/g;
    $seq->seq($seq_string);
    $output_io->write_seq($seq);
  }
}

sub get_input_arg {
  #my $line;
  #print "Getting input arg\n";
  my $line = <>;
  return $line;
}
