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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadSeqsPipeDB;

use strict;
use warnings;
use feature 'say';
use Bio::EnsEMBL::IO::Parser::Fasta;
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

  my $db_path = $self->param('uniprot_db_path');
  my $index_path = $self->param('uniprot_index_path');

  my ($start_index,$end_index) = @{$self->param('uniprot_range')};

  my $parser = Bio::EnsEMBL::IO::Parser::Fasta->open($db_path);
  my $header;
  my $seq;

  my $table_adaptor = $self->db->get_NakedTableAdaptor();
  $table_adaptor->table_name('uniprot_sequences');

  my $index_count = 0;
  while($parser->next() && $index_count <= $end_index) {

    unless($index_count >= $start_index) {
      $index_count++;
      next;
    }

    $header = $parser->getHeader();
    $seq = $parser->getSequence();

    my $output_hash = {};
    $output_hash->{'iid'} = [$header];

    if($index_path) {
      my $cmd = "grep '^".$header."\:' $index_path";
      my $result = `$cmd`;
      chomp $result;

      unless($result) {
        $self->throw("You specified an index file to use but the accession wasn't found in it. Commandline used:\n".$cmd);
      }

      my @result_array = split(':',$result);
      my $database = $result_array[1];
      my $pe_level = $result_array[2];
      my $group = $result_array[3];
      my $biotype = $group."_".$database;
      if($biotype eq '_') {
        $self->throw("Found a malformed biotype based on parsing index. The accession in question was:\n".$header);
      }

      my $db_row = [{ 'accession'  => $header,
                      'source_db'  => $database,
                      'pe_level'   => $pe_level,
                      'biotype'    => $biotype,
                      'group_name' => $group,
                      'seq'        => $seq,
                   }];
      $table_adaptor->store($db_row);
    }

    $self->dataflow_output_id($output_hash,1);

    $index_count++;
  }
  return 1;
}

1;
