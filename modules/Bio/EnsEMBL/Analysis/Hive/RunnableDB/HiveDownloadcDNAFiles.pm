#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

use File::Spec::Functions qw(catfile);
use File::Path qw(make_path);

use Bio::SeqIO;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;

  my $output_path = $self->param_required('embl_sequences')->{output_path};
  if (!-d "$output_path") {
    make_path($output_path);
  }
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

  $self->download_embl_seqs($species,$output_path);

  my $refseq_hash = $self->param('refseq_sequences');
  my $ftp_path = $refseq_hash->{'refseq_ftp'};
  $self->download_refseq_seqs($ftp_path,$output_path);

  $self->unzip($output_path);

  $self->convert_embl_to_fasta($output_path,$taxon_id);

  $self->concat_refseq($species,$output_path);

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
  my @files= grep {$_ =~ /\.dat$/} readdir DIR;
  closedir DIR || $self->throw('Could not close directory');

  my $fasta_out = Bio::SeqIO->new(-format => 'fasta', -file => '>'.catfile($dir, "embl_$taxon_id.fa"));
  $fasta_out->preferred_id_type('accession.version');
  foreach my $file (@files) {
    my $fasta_in = Bio::SeqIO->new(-format => 'embl', -file => catfile($dir, $file));

    while(my $seq = $fasta_in->next_seq) {
      if ($seq->molecule eq 'mRNA') {
        $fasta_out->write_seq($seq);
      }
    }
  }
}

sub download_refseq_seqs {
  my ($self,$ftp,$output_path) = @_;

  say "The cdnas will be downloaded from the RefSeq ftp site";

  system("wget $ftp*.rna.fna.gz -P $output_path");
  say "Finished downloading RefSeq cdna files";
}


sub concat_refseq {
  my ($self,$species,$output_path) = @_;
  opendir DIR, $output_path or $self->throw('Could not open directory '.$output_path);
  my @files= grep /fna$/, readdir DIR;
  closedir(DIR) || $self->throw('Could not open directory '.$output_path);

  my $binomial;
  if ($species eq 'human') {
    $binomial = 'Homo sapiens';
  } elsif ($species eq 'mouse') {
    $binomial = 'Mus musculus';
  }

  my $seqout = Bio::SeqIO->new( -file => '>'.catfile($output_path, 'refseq_'.$species.'.fa'),
                                -format => "Fasta",
                              );
  foreach my $file (@files) {
    if ($file =~/^vertebrate_mammalian/) {
      my $seq_file = new Bio::SeqIO( -file => catfile($output_path, $file),
                                     -format => "Fasta",
                                   );

      while (my $seq = $seq_file->next_seq) {
        if ($seq->id =~ /^N[MR]_/ and $seq->desc =~ /$binomial/) {
          $seqout->write_seq($seq);
        }
      }
    }
  }
}


1;
