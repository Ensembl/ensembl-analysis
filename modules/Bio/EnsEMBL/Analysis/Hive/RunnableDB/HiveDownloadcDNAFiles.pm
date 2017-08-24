#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadcDNAFiles;

use strict;
use warnings;
use feature 'say';

#use Bio::EnsEMBL::IO::Parser::Fasta;
#use Bio::Phylo::Parsers::Fasta;
use Bio::SeqIO;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;
  return 1;
}

sub run {
  my ($self) = shift;

  my $embl_hash = $self->param('embl_sequences');
  my $output_path = $embl_hash->{'output_path'};
  my $species = $embl_hash->{'species'};
  my $taxon_id;

  if ($species eq 'human') {
    $taxon_id = 9606;
  } elsif ($species eq 'mouse') {
    $taxon_id = 10090;
  }

  #$self->download_embl_seqs($species,$output_path);

  my $refseq_hash = $self->param('refseq_sequences');
  my $ftp_path = $refseq_hash->{'refseq_ftp'};
  #$self->download_refseq_seqs($ftp_path,$output_path);

  #$self->unzip($output_path);

  #$self->convert_embl_to_fasta($output_path,$taxon_id);

  #$self->concat_refseq($species,$output_path);

  return 1;
}

sub write_output {
  my $self = shift;
  return 1;
}

sub download_embl_seqs {
  my ($self,$species,$output_path) = @_;

  say "The cdnas will be downloaded from the ENA ftp site";

  # check if the output dir for contigs exists; otherwise, create it
  if (-e "$output_path") {
    say "Output path ".$output_path." found";
  } else {
    `mkdir -p $output_path`;
    if (-e "$output_path") {
      say "Output path $output_path not found.\n".$output_path." created successfully";
    } else {
      $self->throw("Cannot create output path for contigs ".$output_path);
    }
  }
  my @ftp_dirs = ("new/", "release/std/");
  my $ftp = "ftp://ftp.ebi.ac.uk/pub/databases/embl/";
  my @prefix = ("rel_htc_", "rel_std_", "cum_htc_", "cum_std_");

  my $abv;

  if ($species eq 'human') {
    $abv = 'hum';
  } elsif ($species eq 'mouse') {
    $abv = 'mus';
  }
  foreach my $dir (@ftp_dirs) {
    foreach my $pre (@prefix) {
      system("wget -nv $ftp$dir$pre$abv*dat.gz -P $output_path");
    }
  }
  say "Finished downloading EMBL cdna files";
}

sub unzip {
  my ($self,$output_path) = @_;
  say "Unzipping the compressed files...";
  $self->throw("gunzip operation failed. Please check your error log file.") if (system("gunzip -rf $output_path") == 1);
  say "Unzipping finished!";
}

sub convert_embl_to_fasta {
  my ($self,$dir,$taxon_id) = @_; 
  say "Converting EMBL files to Fasta...";
  opendir DIR, $dir or $self->throw("Could not open directory.");
  my @files= readdir DIR;
  closedir DIR;

  my $i = 0;

  foreach my $file (@files) {
    if ($file =~/.dat$/) {
      $file = $dir . $file;
      my $outfile = $file;
      $outfile =~ s/\.dat/\.fasta/;
      my $state = 0;
      my $current_accession = "";
      my $current_version = "";
      my $mol = "";
      my $class = "";
      my $info = "";
      my $current_seq = "";

      open(IN,$file) or die "Can't open '$file': $!";
      open(OUT,">$outfile");
      while(<IN>) {
        my $line = $_;
        if($state == 0 && $line =~ /^ID +([a-zA-Z\d]+)\;.*/) {
          my @line = split (/\s+/, $line);
          $current_accession = $line[1];
          $current_version = $line[3];
          $mol = $line[5];
          $class = $line[6];
          $state = 1;
        } elsif($state == 1 && $line =~ /^DE +([\S+a-zA-Z\d]+).*/) {
          $info = $line;
          $info =~ s/^DE +//;
          chomp $info;
          $state = 2;
        } elsif($state == 2 && ($line =~ /^DE +([\S+a-zA-Z\d]+).*/)) {
          $line =~ s/^DE +//;
          $info .= " ".$line;
          chomp $info;
        } elsif($state == 1 && $line =~ /^SQ +Sequence +/) {
          die "Failed to get description for ".$current_accession."\nExiting";
        } elsif($state == 2 && $line =~ /^SQ +Sequence +/) {
          $state = 3;
        } elsif($state == 3 && !($line =~ /^\/\//)) {
          $current_seq .= uc $line;
        } elsif($state == 3 && $line =~ /^\/\//) {
          $current_accession =~ s/;//g;
          $current_version =~ s/;//g;
          $mol =~ s/;//g;
          $class =~ s/;//g;
          $current_seq =~ s/ //g;
          $current_seq =~ s/\d+//g;
          if ($mol eq 'mRNA' && ($class eq 'STD' || $class eq 'HTC' || $class eq 'PAT')) {
            print OUT '>'.$current_accession.".".$current_version." ".$info."\n".$current_seq;
          }
          $state = 0;
          $current_accession = '';
          $current_version = '';
          $current_seq = '';
        }
      }
      close OUT;
      close IN;
    }
  }
  system("cat $dir/*fasta > $dir/embl_$taxon_id.fa");
  #system("rm *fasta");
  #system("rm *dat");
}

sub download_refseq_seqs {
  my ($self,$ftp,$output_path) = @_;

  say "The cdnas will be downloaded from the RefSeq ftp site";

  # check if the output dir for contigs exists; otherwise, create it
  if (-e "$output_path") {
    say "Output path ".$output_path." found";
  } else {
    `mkdir -p $output_path`;
    if (-e "$output_path") {
      say "Output path $output_path not found.\n".$output_path." created successfully";
    } else {
      $self->throw("Cannot create output path ".$output_path);
    }
  }
  system("wget $ftp*fna.gz -P $output_path");
  say "Finished downloading RefSeq cdna files";
}


sub concat_refseq {
  my ($self,$species,$output_path) = @_;
  opendir DIR, $output_path or $self->throw("Could not open directory.");
  my @files= grep /fna$/, readdir DIR;
  closedir DIR;

  my $binomial;
  my $outfile;
  if ($species eq 'human') {
    $binomial = 'Homo sapiens';
    $outfile = $output_path."/refseq_human.fa";
  } elsif ($species eq 'mouse') {
    $binomial = 'Mus musculus';
    $outfile = $output_path."/refseq_mouse.fa";
  }

  foreach my $file (@files) {
    if ($file =~/^vertebrate_mammalian/) {
      my $seq_file = new Bio::SeqIO( -file => $output_path."/".$file,
                                     -format => "Fasta",
                                   );

      while (my $seq = $seq_file->next_seq) {
        my $desc = $seq->desc;
        my $seqout = Bio::SeqIO->new( -file => ">>$outfile",
                                      -format => "Fasta",
                                    );
        if ($desc =~ /$binomial/) {
          $seqout->write_seq($seq);
        }
      }
    }
  }
}

#sub download_refseq_seqs {
#  my ($self,$refseq_file,$refseq_dir,$output_path) = @_;
#  my $cmd = 'scp -p '.$refseq_dir.'/'.$refseq_file.' '.$output_path.'/'. $refseq_file;
#  print "Copying RefSeq file containing cDNAs:\n", $cmd, "\n";

#  if (system($cmd) == 1) {
#    die "Failed to copy RefSeq file.\nExiting"; 
#  } else {
#    print "RefSeq file copied.\n";
#  }
#}

1;

