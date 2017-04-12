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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadcDNAs;

use strict;
use warnings;
use feature 'say';
#use Bio::EnsEMBL::IO::Parser::Fasta;
use Bio::SeqIO;
use Data::Dumper;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;

  return 1;
}

sub run {
  my $self = shift;

  return 1;
}

sub write_output {
  my $self = shift;

  my $db_path = $self->param('cdna_file');

  my $seq_file = new Bio::SeqIO( -file => $db_path,
                                 -format => "Fasta",
                               );

  my $header;
  my $seq;
  my $biotype = 'cdna';

  my $table_adaptor = $self->db->get_NakedTableAdaptor();
  $table_adaptor->table_name('cdna_sequences');

  while (my $seq = $seq_file->next_seq) {
    my $header = $seq->id;
    my $sequence = $seq->seq;
    $header =~ s/(\S+).*/$1/;

    my $output_hash = {};
    $output_hash->{'iid'} = $header;

    my $db_row = [{
      'accession'  => $header,
      'seq'        => $sequence,
      'biotype'    => $biotype,
    }];
    $table_adaptor->store($db_row);
    
    $self->dataflow_output_id($output_hash,2);
  }
  return 1;
}

1;
