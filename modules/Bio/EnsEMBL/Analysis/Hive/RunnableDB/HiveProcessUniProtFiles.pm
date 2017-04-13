#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2017] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProcessUniProtFiles;

use strict;
use warnings;
use feature 'say';
use Bio::EnsEMBL::IO::Parser::Fasta;
use Bio::EnsEMBL::KillList::KillList;

use parent ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    min_seq_length => 50,
    column_names => ['iid'],
  }
}


sub fetch_input {
  my $self = shift;

  my $file_path = $self->param_required('iid');
  $self->throw("The input id doesn't exist, offending path:\n$file_path")
    unless(-e $file_path);
  my $input_seq_count = 0;
  my $below_min_length_count = 0;
  my $killed_count = 0;

  my $kill_list_object = Bio::EnsEMBL::KillList::KillList->new(-TYPE => $self->param('killlist_type'),
                                                                -KILL_LIST_DB => $self->param('killlist_db'),
                                                                -FILTER_PARAMS => $self->param('KILL_LIST_FILTER'));
  my $kill_list = $kill_list_object->get_kill_list();


  my $min_seq_length = $self->param('min_seq_length');

  my $table_adaptor = $self->db->get_NakedTableAdaptor();
  $table_adaptor->table_name($self->param_required('sequence_table_name'));
  my ($group_name) = $file_path =~ /(\w+)\.fasta/;

  my $parser = Bio::EnsEMBL::IO::Parser::Fasta->open($file_path);
  my $header;
  my $seq;

  my @iids;
  while($parser->next()) {
    $input_seq_count++;
    $header = $parser->getHeader();
    $seq = $parser->getSequence();
    if ($header =~ /^([^\|]+)\|([^\|]+)\|.+ PE\=([1-5]) SV\=(\d+)/) {
      my ($source_db, $accession, $pe_level, $sequence_version) = ($1, $2, $3, $4);
      if(length($seq) < $min_seq_length) {
        $below_min_length_count++;
        next;
      } elsif(exists($kill_list->{$accession})) {
        say "Removing ".$accession." as it is present in kill list";
        $killed_count++;
        next;
      }

      my $versioned_accession = $accession.'.'.$sequence_version;
      my $db_row = [{ 'accession'  => $versioned_accession,
                      'source_db'  => $source_db,
                      'pe_level'   => $pe_level,
                      'biotype'    => $group_name.'_'.$source_db,
                      'group_name' => $group_name,
                      'seq'        => $seq,
                   }];
      $table_adaptor->store($db_row);
      push(@iids, $versioned_accession);
    }
    else {
      $self->throw("Matched a header but couldn't parse it fully. Expect to find sp/tr, accession, ".
                   "pe level and sequence version. Offending header:\n".$header);
    }
  }

  $self->warning("Total sequences: $input_seq_count\nSequences below $min_seq_length: $below_min_length_count\nKilled sequences: $killed_count\nStored sequences: ".scalar(@iids)."\n");
  $self->param('inputlist', \@iids);
}

1;
