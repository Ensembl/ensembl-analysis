#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

use parent ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults()},
    column_names => ['iid'],
    sequence_table_name => 'protein_sequences',
    load_biotype => 0,
  }
}

sub fetch_input {
  my $self = shift;

  my $parser = Bio::EnsEMBL::IO::Parser::Fasta->open($self->param_required('protein_file'));

  my $table_adaptor = $self->db->get_NakedTableAdaptor();
  $table_adaptor->table_name($self->param_required('sequence_table_name'));

  my @iids;
  while($parser->next()) {
    my ($accession) = $parser->getHeader =~ /^(\S+)/;
    my $db_row = [{
      'accession'  => $accession,
      'seq'        => $parser->getSequence,
    }];
    if ($self->param('load_biotype')) {
      if ($parser->getHeader =~ /\S+\s+(\S+)/) {
        $db_row->[0]->{biotype} = $1;
      }
      else {
        $self->warning('Could not find biotype for '.$accession);
      }
    }
    $table_adaptor->store($db_row);
    push(@iids, $accession);
  }
  $self->param('inputlist', \@iids);
}

1;
