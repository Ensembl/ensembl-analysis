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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProcessUniProtFiles;

use strict;
use warnings;
use feature 'say';
use Bio::EnsEMBL::IO::Parser::Fasta;
use Bio::EnsEMBL::KillList::HiveKillList;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;

  # I want this to take a set of input files in fasta format, each with a group name,
  # loop through each file, check seq against killlist, parse out the info from the header
  # and group and output to an index. I would also like to output a set of ranges to the
  # next analysis, basically for loading into the pipe db

  return 1;
}

sub run {
  my $self = shift;

  $self->process_fasta_files();

  return 1;
}

sub write_output {
  my $self = shift;

  my $output_ids = $self->output_ids();
  foreach my $output_id (@{$output_ids}) {
    my $output_hash = {};
    $output_hash->{'uniprot_range'} = $output_id;
    $self->dataflow_output_id($output_hash,1);
  }

  return 1;
}


sub process_fasta_files {
  my ($self) = @_;

  my $input_seq_count = 0;
  my $output_seq_count = 0;
  my $below_min_length_count = 0;
  my $killed_count = 0;

  my $output_batch_size = 200;
  my $min_seq_length = 50;
  my $index_file_name = 'uniprot_index';
  my $db_file_name = 'uniprot_db';

  my $input_files_hash = $self->param('iid');

  my $kill_list_object;
  my %kill_list;

  $kill_list_object = Bio::EnsEMBL::KillList::HiveKillList->new(-TYPE => $self->param('killlist_type'),
                                                                -KILL_LIST_DB => $self->param('killlist_db'),
                                                                -FILTER_PARAMS => $self->param('KILL_LIST_FILTER'));
  %kill_list = %{$kill_list_object->get_kill_list()};


  unless($input_files_hash) {
    $self->throw("No input id found, this module requires an input id to be passed in under 'iid'");
  }

  if($self->param('output_batch_size')) {
    $output_batch_size = $self->param('output_batch_size');
  } else {
    $self->warning("You have not provided an output id batch size using 'output_batch_size'. Will default to:\n".$output_batch_size);
  }

  if($self->param('min_seq_length')) {
    $min_seq_length = $self->param('min_seq_length');
  } else {
    $self->warning("You have not provided a minimum seq length using 'min_seq_length'. Will default to:\n".$min_seq_length);
  }

  if($self->param('db_file_name')) {
    $db_file_name = $self->param('db_file_name');
  } else {
    $self->warning("You have not provided an output db name using 'db_file_name'. Will default to:\n".$db_file_name);
  }

  if($self->param('index_file_name')) {
    $index_file_name = $self->param('index_file_name');
  } else {
    $self->warning("You have not provided an index file name using 'index_file_name'. Will default to:\n".$index_file_name);
  }

  my $index_file_path = $self->param('dest_dir').'/'.$index_file_name;
  my $db_file_path = $self->param('dest_dir').'/'.$db_file_name;

  unless(-e $self->param('dest_dir')) {
   `mkdir -p $self->param('dest_dir')`;
  }

  open(OUT_DB,">".$db_file_path);
  open(OUT_INDEX,">".$index_file_path);
  foreach my $group_name (keys(%$input_files_hash)) {
    my $file_path = $input_files_hash->{$group_name};
    unless(-e $file_path) {
      $self->throw("One or more of the files in the input id don't exist. Offending path:\n".$file_path);
    }

    my $parser = Bio::EnsEMBL::IO::Parser::Fasta->open($file_path);
    my $header;
    my $seq;

    while($parser->next()) {
      $input_seq_count++;
      $header = $parser->getHeader();
      $seq = $parser->getSequence();
      my ($accession,$sequence_version,$source_db,$pe_level) = @{$self->parse_header_info($header)};
      if(length($seq) < $min_seq_length) {
        $below_min_length_count++;
        next;
      } elsif(exists($kill_list{$accession})) {
        say "Removing ".$accession." as it is present in kill list";
        $killed_count++;
        next;
      }

      my $versioned_accession = $accession.'.'.$sequence_version;
      my $index_line = $versioned_accession.':'.$source_db.':'.$pe_level.':'.$group_name;
      my $record = ">".$versioned_accession."\n".$seq;
      say OUT_INDEX $index_line;
      say OUT_DB $record;
      $output_seq_count++;
    }
  }
  close OUT_DB;
  close OUT_INDEX;

  say "After cleaning:";
  say "Input sequence count: ".$input_seq_count;
  say "Output sequence count: ".$output_seq_count;
  say "Short sequence count: ".$below_min_length_count;
  say "Killed sequence count: ".$killed_count;

  my $output_ids = $self->build_output_ranges($output_seq_count,$output_batch_size);
  $self->output_ids($output_ids);

}

sub parse_header_info {
  my ($self,$header) = @_;

  unless($header =~ /^([^\|]+)\|([^\|]+)\|.+ PE\=([1-5]) SV\=(\d+)/) {
    $self->throw("Matched a header but couldn't parse it fully. Expect to find sp/tr, accession, ".
                 "pe level and sequence version. Offending header:\n".$header);
  }

  my $source_db = $1;
  my $accession = $2;
  my $pe_level = $3;
  my $sequence_version = $4;

  return([$accession,$sequence_version,$source_db,$pe_level]);
}


sub build_output_ranges {
  my ($self,$seq_count,$output_batch_size) = @_;

  my $output_ids = [];
  my $loop_count = int($seq_count / $output_batch_size);
  my $remainder = $seq_count % $output_batch_size;
  my $start_index = 0;
  my $end_index = 0;
  my $i=0;
  for($i=0; $i<$loop_count; $i++) {
    $start_index = $i * $output_batch_size;
    $end_index = $start_index + $output_batch_size - 1;
    push(@{$output_ids},[$start_index,$end_index]);
  }

  if($remainder) {
    push(@{$output_ids},[($end_index + 1),($seq_count - 1)]);
  }

  return($output_ids);
}

sub output_ids {
  my ($self,$val) = @_;

  if($val) {
    $self->param('_output_ids',$val);
  }

  return $self->param('_output_ids');
}

1;
