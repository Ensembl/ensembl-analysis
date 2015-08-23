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


use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;

  return 1;
}

sub run {
  my $self = shift;

  $self->index_fasta_files();

  return 1;
}

sub write_output {
  my $self = shift;

  my $input_files_hash = $self->param('iid');
  foreach my $group_name (keys(%$input_files_hash)) {
    my $file_path = $input_files_hash->{$group_name};
    my $output_path = $self->param('dest_dir');
    my $output_hash = {};
    $output_hash->{'iid'} = $file_path;
    $output_hash->{'chunk_dir_name'} = $group_name."_chunks";
    $self->dataflow_output_id($output_hash,1);
  }

  return 1;
}

sub index_fasta_files {
  my ($self) = @_;

  my $input_files_hash = $self->param('iid');
  my $index_file_path = $self->param('dest_dir').'/'.$self->param('index_file_name');
  unless(-e $self->param('dest_dir')) {
   `mkdir -p $self->param('dest_dir')`;
  }

  open(OUT, ">".$index_file_path);
  foreach my $group_name (keys(%$input_files_hash)) {
    my $file_path = $input_files_hash->{$group_name};
    unless(-e $file_path) {
     $self->throw("One or more of the files in the input id don't exist. Offending path:\n".$file_path);
    }

    open(IN,$file_path);
    while(<IN>) {
      my $line = $_;
      if($line =~ /^\>/) {
        unless($line =~ /^\>([^\|]+)\|([^\|]+)\|.+ PE\=([1-5]) SV\=(\d+)/) {
          $self->throw("Matched a header but couldn't parse it fully. Expect to find sp/tr, accession, ".
                       "pe level and sequence version. Offending header:\n".$line);
        }

        my $database = $1;
        my $accession = $2;
        my $pe_level = $3;
        my $sequence_version = $4;
        say OUT $accession.":".$database.":".$pe_level.":".$sequence_version.":".$group_name;
      }
    }
    close IN;
  }
  close OUT;
}

1;
