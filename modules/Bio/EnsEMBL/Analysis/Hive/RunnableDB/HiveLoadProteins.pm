#!/usr/bin/env perl

# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadProteins;

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

  my $db_path = $self->param('protein_file');
  my $protein_table_name = $self->param('protein_table_name');
  my $start_index = 0;

  my $parser = Bio::EnsEMBL::IO::Parser::Fasta->open($db_path);
  my $header;
  my $seq;


  my $table_adaptor = $self->db->get_NakedTableAdaptor();
  $table_adaptor->table_name($protein_table_name);

  while($parser->next()) {
    $header = $parser->getHeader();
    $seq = $parser->getSequence();

    my $output_hash = {};
    $output_hash->{'iid'} = [$header];

    my $db_row = [{
      'accession'  => $header,
      'seq'        => $seq,
    }];
    $table_adaptor->store($db_row);

    $self->dataflow_output_id($output_hash,1);
  }
  return 1;
}

1;
