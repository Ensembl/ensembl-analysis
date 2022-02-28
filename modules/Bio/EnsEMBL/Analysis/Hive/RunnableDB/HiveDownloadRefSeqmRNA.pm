#!/usr/bin/env perl

# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadRefSeqmRNA;

use strict;
use warnings;
use feature 'say';
use Bio::EnsEMBL::Utils::PolyA;
use Bio::SeqIO;
#use Bio::EnsEMBL::IO::Parser::Fasta;
use POSIX qw(strftime);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;

  my $ftp_path = $self->param('ftp_path');
  my $output_path = $self->param('output_path');

  unless(-e $output_path) {
    system("mkdir -p ".$output_path);
  }

  my $cmd = "wget ".$ftp_path." -O ".$output_path."/refseq_cdnas.fa.gz";
  my $result = system($cmd);
  if($result) {
    $self->throw("Unable to download the RefSeq mRNA file. Commandline used:\n".$cmd);
  }

  $cmd = "gunzip -f ".$output_path."/refseq_cdnas.fa.gz";
  $result = system($cmd);
  if($result) {
    $self->throw("Unable to unzip the RefSeq mRNA file. Commandline used:\n".$cmd);
  }

  return 1;
}

sub run {
  my $self = shift;

  return 1;
}

sub write_output {
  my $self = shift;

  my $db_path = $self->param('output_path')."/refseq_cdnas.fa";
  my $min_length = $self->param("min_seq_length");
  unless($min_length) {
    $min_length = 60;
  }

#  my $parser = Bio::EnsEMBL::IO::Parser::Fasta->open($db_path);
  my $date = strftime "%m/%d/%Y", localtime;

  my $table_adaptor = $self->db->get_NakedTableAdaptor();
  $table_adaptor->table_name('refseq_sequences');

  my $seq_file = new Bio::SeqIO( -file   => $db_path,
                              -format => "Fasta", );

  while (my $seq = $seq_file->next_seq) {
    my $biotype = 'cdna';
    my $header = $seq->display_id;
    say "Header: ".$header;

    if($header =~ /^XR/ || $header =~ /^NR/) {
      say "Not storing RNA accession: ".$header;
      next;
    }

    unless($header =~ /^NM/ || $header =~ /^AM/) {
      $biotype .= "_predicted";
    }

    my $polyA_clipper = Bio::EnsEMBL::Utils::PolyA->new();

    my $clipped_seq = $polyA_clipper->clip($seq);
    unless (defined $clipped_seq) {
      say "Clipping made a cdna vanish...";
      next;
    }

    if ($clipped_seq->length < $min_length) {
      say "Sequence is below min length...";
      next;
    }

    my $output_hash = {};
    $output_hash->{'iid'} = [$header];

    my $db_row = [{ 'accession'  => $header,
                    'source_db'  => 'RefSeq',
                    'date'       => $date,
                    'seq'        => $clipped_seq->seq,
                    'biotype'    => $biotype,
                  }];
    $table_adaptor->store($db_row);

    $self->dataflow_output_id($output_hash,2);
  }
  return 1;
}

1;
