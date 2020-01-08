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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCleanFastaHeaders;

use strict;
use warnings;
use feature 'say';

use Bio::Seq;
use Bio::SeqIO;
use Bio::EnsEMBL::KillList::KillList;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;
  return 1;
}

sub run {
  my $self = shift;
  $self->clean_headers;
  return 1;
}

sub write_output {
  my $self = shift;

  my $file_path = $self->output_file_path;

  my $output_hash = {};
  $output_hash->{'iid'} = $file_path;

  # If this param exists pass it along. This is for chunking analyses further down, if there
  # are multiple input ids flowing into the chunking analysis then it needs the specific chunk
  # dir name for each file from the indexing step
  if($self->param('chunk_dir_name')) {
    $output_hash->{'chunk_dir_name'} = $self->param('chunk_dir_name');
  }
  $self->dataflow_output_id($output_hash,1);

  return 1;
}

sub clean_headers {
  my ($self) = @_;

  my $input_file;
  my $output_file;
  my $min_seq_length = 10;
  my $skip_broken_headers = 0;
  my $header_source;
  my $use_killlist = 0;
  my $kill_list_object;
  my %kill_list;

  # Either take the input file path from the parameters hash if it exists, otherwise parse it from
  # the input id. Throw is you can't find either file
  if($self->param_is_defined('input_file_path')) {
    $input_file = $self->param('input_file_path');
    unless(-e $input_file) {
      $self->throw("You have specified a path to the protein file in the pipeline config, but this file does not exist. Note that ".
                   "specifying a path overwrites any input id that might be present, so remove the path from the config/analysis_base ".
                   "table if you want to use the input id as the path");
    }
  } else {
    $input_file = $self->input_id;
    unless(-e $input_file) {
      $self->throw("You have not specified an input_file_path variable in your hive config, so the input id was used instead. When the ".
                   "input id was parsed the resulting file path does not exist.\nInput id (unparsed):\n".$self->input_id."\nParsed input id:\n".
                   $input_file);
    }
  }

  # Open output file, either from parameters hash first or failing that by putting _clean on the input file path
  if($self->param_is_defined('output_file_path')) {
    say "You have defined an output file path to write to";
    $output_file = $self->param('output_file_path');
  } else {
    $output_file = $input_file."_clean";
  }

  # Set this for write_output
  $self->output_file_path($output_file);

  # Set the min seq length if defined, otherwise use the default (10)
  if($self->param_is_defined('min_seq_length')) {
    $min_seq_length = $self->param('min_seq_length');
  }

  # Unless a source is defined throw
  unless($self->param_is_defined('header_source')) {
    $self->throw("You have not defined the header_source in your config file (e.g. header_source => 'uniprot')");
  }

  # Set the flag, unless the string is 'no'. Note that if it's 0 it will get set through this but will
  # be considered false later and therefore not used
  if($self->param_is_defined('skip_broken_headers') && ($self->param('skip_broken_headers') eq "1" ||
     $self->param('skip_broken_headers') eq "yes" || $self->param('skip_broken_headers') eq "YES")) {
    $skip_broken_headers = 1;
  }

  # Use the killlist or not
  if($self->param_is_defined('use_killlist') && ($self->param('use_killlist') eq "1" ||
     $self->param('use_killlist') eq "yes" || $self->param('use_killlist') eq "YES")) {
    say "Using the killlist";
    unless($self->param_is_defined('killlist_type')) {
      $self->throw("You have selected to use the killlist but haven't defined a killlist_type in your pipeline config, e.g ".
                   " 'killlist_type' => 'protein'");
    }

    unless($self->param_is_defined('killlist_db')) {
      $self->throw("You have selected to use the killlist but haven't defined a killlist_db in your pipeline config, e.g ".
                   "'killlist_db' => $self->o('killlist_db')");
    }

    unless($self->param_is_defined('KILL_LIST_FILTER')) {
      say "You have selected to use the killlist but haven't defined a KILL_LIST_FILTER hash in your pipeline config, ".
          "the KillList module will look for a default hash for your molecule type";
    }

    say "Killlist molecule type set to:\n".$self->param('killlist_type');
    $use_killlist = 1;
    $kill_list_object = Bio::EnsEMBL::KillList::KillList->new(-TYPE => $self->param('killlist_type'),
                                                                  -KILL_LIST_DB => $self->param('killlist_db'),
                                                                  -FILTER_PARAMS => $self->param('KILL_LIST_FILTER'),
                                                                 );
    %kill_list = %{ $kill_list_object->get_kill_list() };
  } else {
    say "Not using the killlist";
  }

  say "Reading from input file:\n".$input_file;
  say "Will write cleaned seqs to:\n".$output_file;
  say "Min seq length set to:\n".$min_seq_length;
  say "Skip broken headers set to:\n".$skip_broken_headers;
  say "Source set to:\n".$self->param('header_source');

  my $seqin  = new Bio::SeqIO( -file   => "<$input_file",
			       -format => "Fasta",
                             );

  my $seqout = new Bio::SeqIO( -file   => ">$output_file",
                               -format => "Fasta"
                             );

  # Some counts for the stats at the end
  my $input_seq_count = 0;
  my $output_seq_count = 0;
  my $short_count = 0;
  my $skip_count = 0;
  my $duplicate_count = 0;
  my $killed_count = 0;

  # This hash is just used to track duplicate accessions
  my %uniprot_accession_hash = ();

  # Loop through the set of sequences
  while(my $prot = $seqin->next_seq) {
    $input_seq_count++;

    # If it's below the min seq length then skip it
    if($prot->length <= $min_seq_length) {
      say STDERR "Length < ".$min_seq_length.": rejecting ".$prot->display_id." with length ".$prot->length;
      $short_count++;
      next;
    }

    # Try and find the parse the accession
    my $display_id = $prot->display_id;
    my $uniprot_accession;
    $uniprot_accession = $self->match_against_source($display_id);

    # If nothing was parsed then check if the skip headers flag is in use. If it is then skip, otherwise throw
    unless($uniprot_accession) {
      if($skip_broken_headers) {
        say STDERR "Skipping the following header/sequence as accession can't be parsed and skip_broken_headers is in use:\n".
                   $display_id;
        $skip_count++;
        next;
      } else {
        $self->throw("Could not match a uniprot accession in the header. Header:\n".$display_id);
      }
    }

    # Check if the accession is already in the hash, if it is then skip
    if($uniprot_accession_hash{$uniprot_accession}) {
      say STDERR "Skipping accession as it has already been seen. Accession:\n".$uniprot_accession;
      $duplicate_count++;
      next;
    }

    # If it's in the killlist then remove
    if ($use_killlist) {
      if (exists( $kill_list{$uniprot_accession})) {
        say STDERR "$uniprot_accession is present in kill list DB, discarded.\n";
        $killed_count++;
        next;
      }
    }

    # At this point this is the first time the accession was seen so set it
    $uniprot_accession_hash{$uniprot_accession} = 1;

    # Write out the cleaned sequence
    $prot->display_id($uniprot_accession);
    $prot->desc("");
    $prot->seq(uc($prot->seq));
    $seqout->write_seq($prot);
    $output_seq_count++;
  }

  say "After cleaning:";
  say "Input sequence count: ".$input_seq_count;
  say "Output sequence count: ".$output_seq_count;
  say "Skipped sequence count: ".$skip_count;
  say "Short sequence count: ".$short_count;
  say "Duplicate sequence count: ".$duplicate_count;
  say "Killed sequence count: ".$killed_count;
}

sub match_against_source {
  my ($self,$display_id) = @_;

  my $uniprot_accession;

  # Set the source
  my $header_source = $self->param('header_source');

  if($header_source eq 'uniprot') {
    # if this is working on something that has come from the processing module, then the headers will already be fixed and versioned
    if($self->param('header_pre_cleaned')) {
      if($display_id =~ /^.+\.\d+/) {
        $uniprot_accession = $&;
      }
    } elsif($display_id =~ /^(sp|tr)\|([^\|]+)\|.+ SV\=(\d+)/) {
      $uniprot_accession = $2.'.'.$3;
    }
  }

  else {
    $self->throw("You have entered a source that is not supported. The code must be updated to ".
                 "deal with it. Source:\n".$header_source);
  }

  return($uniprot_accession);
}

sub output_file_path {
  my ($self,$path) = @_;

  if(defined($path)) {
    $self->param('_output_file_path',$path);
  }

  return($self->param('_output_file_path'));
}

1;
