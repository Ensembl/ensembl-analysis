#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

use File::Spec::Functions qw(splitpath);
use Bio::EnsEMBL::IO::Parser::Fasta;
use Bio::EnsEMBL::KillList::KillList;

use Bio::SeqIO;
use Bio::Seq;

use parent ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    min_seq_length => 50,
    column_names => ['iid'],
    output_file => undef, # by default we do not write the file it's need for targetted at the moment
    skip_Xs => undef, # We only remove data with too many X in targetted, the value should be 5
    delete_file => 1, # Delete the output_file if it exists
  }
}


sub fetch_input {
  my $self = shift;

  my $files = $self->param_required('iid');
  my $input_seq_count = 0;
  my $below_min_length_count = 0;
  my $too_many_X = 0;
  my $killed_count = 0;

  my $kill_list;
  if ($self->param_is_defined('killlist_db')) {
    my $kill_list_object = Bio::EnsEMBL::KillList::KillList->new(-TYPE => $self->param_required('killlist_type'),
                                                                  -KILL_LIST_DB => $self->param('killlist_db'),
                                                                  -FILTER_PARAMS => $self->param('KILL_LIST_FILTER'));
    $kill_list = $kill_list_object->get_kill_list();
  }


  my $min_seq_length = $self->param('min_seq_length');

  my $table_adaptor = $self->db->get_NakedTableAdaptor();
  $table_adaptor->table_name($self->param_required('sequence_table_name'));
  if (ref($files) ne 'ARRAY') {
    $files = [$files];
  }
  my @iids;
  my @seqs;
  my $skip_X = $self->param('skip_Xs');
  my $output_file = $self->param('output_file');
  if ($output_file) {
    unlink $output_file if ($self->param('delete_file') and -e $output_file);
    $self->param('iid', $output_file);
    $output_file = Bio::SeqIO->new(-format => 'fasta', -file => ">$output_file");
  }
  foreach my $file_path (@$files) {
    $self->throw("The input id doesn't exist, offending path:\n$file_path")
      unless(-e $file_path);
    my (undef, undef, $group_name) = splitpath($file_path);
    $group_name =~ s/\.\w+$//;

    my $parser = Bio::EnsEMBL::IO::Parser::Fasta->open($file_path);

    while($parser->next()) {
      $input_seq_count++;
      my $seq = $parser->getSequence();
      if(length($seq) < $min_seq_length) {
        $below_min_length_count++;
      }
      elsif ($skip_X and $seq =~ /X{$skip_X}/) {
        $too_many_X++;
      }
      else {
        my $header = $parser->getHeader();
        if ($seq =~ /(B|J|Z)/) {
          $self->warning("You have a $1 for $header, this may cause problems");
        }
        my ($source_db, $accession, $pe_level, $sequence_version, $versioned_accession, $last_pe);
        if ($header =~ /^([^\|]+)\|([^\|]+)\|.+ PE\=([1-5]) SV\=(\d+)/) {
          ($source_db, $accession, $pe_level, $sequence_version) = ($1, $2, $3, $4);
          $versioned_accession = $accession.'.'.$sequence_version;
          $last_pe = $pe_level;
        }
        elsif ($header =~ /^([sptr]{2})\|([^\|]+-\d+)\|/) {
# If you ask for isoforms, the first sequence has the PE level
          ($source_db, $accession) = ($1, $2);
          $versioned_accession = $accession;
          $pe_level = $last_pe;
        }
        elsif ($header =~ /^(\w+\.\d+)/) {
          $versioned_accession = $1;
          next if ($versioned_accession =~ /^XP_/);
          if ($versioned_accession =~ /^NP_/) {
            $source_db = 'refseq';
          }
          else {
            $source_db = 'uniprot';
          }
        }
        else {
          $self->throw("Matched a header but couldn't parse it fully. Expect to find sp/tr, accession, ".
                       "pe level and sequence version. Offending header:\n".$header);
        }
        if(exists($kill_list->{$accession})) {
          say "Removing ".$accession." as it is present in kill list";
          $killed_count++;
        }
        else {
          my $db_row = [{ 'accession'  => $versioned_accession,
                          'source_db'  => $source_db,
                          'pe_level'   => $pe_level,
                          'biotype'    => $group_name.'_'.$source_db,
                          'group_name' => $group_name,
                          'seq'        => $seq,
                       }];
          $table_adaptor->store($db_row);
          push(@iids, $versioned_accession);
          if ($output_file) {
            $output_file->write_seq(Bio::Seq->new(-id => $versioned_accession, -seq => $seq));
          }
        }
      }
    }
  }

  $self->warning("Total sequences: $input_seq_count\nSequences below $min_seq_length: $below_min_length_count\nKilled sequences: $killed_count\nStored sequences: ".scalar(@iids)."\n");
  $self->param('inputlist', \@iids);
}

1;
