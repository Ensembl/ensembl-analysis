#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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
use Data::Dumper;
use File::Spec::Functions qw(splitpath);
use Bio::EnsEMBL::IO::Parser::Fasta;
use Bio::EnsEMBL::KillList::KillList;

use parent ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    min_seq_length => 50,
    column_names => ['iid'],
    skip_Xs => undef, # We only remove data with too many X in targetted, the value should be 5
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Parse the protein file to load into the database and be able to align the sequence
              If the protein has been added to the killlist database, it will be removed from the
              set.
              When skip_Xs is set, the protein is removed from the set if it contains as many X in
              the sequence as the value given in the parameter.
              If the sequence has B, J or Z, a warning is printed as some software may have problem
              to use them.
              When the protein has a uniprot header, the pe level is set accordingly and the source
              is set to uniprot.
              When the protein has a RefSeq accession, [AN]P_, the pe level is set to 2 and the source
              is set to refseq. Predicted proteins [YX]P_ are removed from the set
              If the protein has a selenocysteine, the source is set to seleno
 Returntype : None
 Exceptions : Throws if the header cannot be parsed
              Throws if it failed to store the sequence

=cut

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
  my $write_file = 0;
  if ($self->param_is_defined('output_file') and $self->param('output_file')) {
    $write_file = 1;
    open(FH, '>'.$self->param('output_file')) || $self->throw('Could not open '.$self->param('output_file'));
  }
  foreach my $file_path (@$files) {
    unless(-e $file_path){
      $self->warning("The input id doesn't exist, offending path:\n$file_path");
      next;
    }
    my (undef, undef, $group_name) = splitpath($file_path);
    $group_name =~ s/\.\w+$//;

    my $parser = Bio::EnsEMBL::IO::Parser::Fasta->open($file_path);
    my $last_pe = 2;
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
        my ($source_db, $accession, $pe_level, $sequence_version, $versioned_accession);
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
          next if ($versioned_accession =~ /^[YX]P_/);
          ($accession) = $versioned_accession =~ /^(.+)\.\d+$/;
          if ($versioned_accession =~ /^[AN]P_/) {
            $source_db = 'refseq';
            $pe_level = 2;
          }
          else {
            $source_db = 'uniprot';
          }
        }
        else {
          $self->throw("Matched a header but couldn't parse it fully. Expect to find sp/tr, accession, ".
                       "pe level and sequence version. Offending header:\n".$header);
        }
        if ($seq =~ /U/) {
          $source_db = 'seleno';
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
          eval {
            $table_adaptor->store($db_row);
          };

          if($@) {
            my $except = $@;
            unless($except =~ /Duplicate entry/) {
              $self->throw("Issue strong the following accession: ".$versioned_accession."\n".$except);
	    }
          }
          print FH ">$versioned_accession\n$seq\n" if ($write_file);
          push(@iids, $versioned_accession);
        }
      }
    }
  }
  if ($write_file) {
    close(FH) || $self->throw('Could not close '.$self->param('output_file'));
  }

  $self->warning("Total sequences: $input_seq_count\nSequences below $min_seq_length: $below_min_length_count\nKilled sequences: $killed_count\nStored sequences: ".scalar(@iids)."\n");

  if (@iids) {
    $self->param('inputlist', \@iids);
  }
  else {
    if ($write_file and -e $self->param('output_file')) {
       unlink $self->param('output_file');
    }
    $self->input_job->autoflow(0);
    $self->complete_early('No sequences have been stored or written to file');
  }
}

1;
