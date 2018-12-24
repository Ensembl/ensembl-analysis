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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveParseCsvIntoTable;

use strict;
use warnings;

use parent ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');


=head2 param_defaults

 Arg [1]    : None
 Description: Defaults parameters for the module
               _sample_column_size => 50, #Trying to detect using DBI otherwise it should be set to the value of the 'sample_name' column to be effective
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    _sample_column_size => 50,
  }
}


=head2 write_output

 Arg [1]    : None
 Description: Store the csv information in a custom table in the hive database
              via 'csvfile_table' and flowing the information from 'sample_column'
              into branch 'fan_branch_code', usually 2
              If you specify -1 for the 'is_mate_1' column, it will use the value
              in 'pairing_regex' to determine which mate correspond to the filename
              processed
 Returntype : None
 Exceptions : Throws if the file has more than two rows with the same ID

=cut

sub write_output {
    my $self = shift;

    my %keyword_hash;
    my %id_check;
    my $table_adaptor = $self->db->get_NakedTableAdaptor;
    $table_adaptor->table_name($self->param('csvfile_table'));
    my $db_column_info = $table_adaptor->dbc->db_handle->column_info(undef, undef, $self->param('csvfile_table'), $self->param('sample_column'));
    my $sample_column_size = $self->param('_sample_column_size');
    if ($db_column_info) {
      my $info = $db_column_info->fetchall_arrayref();
      $sample_column_size = $info->[0]->[6];
      $self->say_with_header($info->[0]->[6]);
    }
    foreach my $input_id (@{$self->param('output_ids')}) {
        $input_id->{$self->param('sample_column')} =~ tr/ :\t/_/;
        if (length($input_id->{$self->param('sample_column')}) > $sample_column_size) {
          $self->throw('Sample '.$input_id->{$self->param('sample_column')}.' from '.$self->param('sample_column')." is bigger than $sample_column_size");
        }
        if (exists $id_check{$input_id->{ID}}) {
          ++$id_check{$input_id->{ID}};
          $self->throw("You should only have one or two file with the same ID") if ($id_check{$input_id->{ID}} > 2);
        }
        else {
          $id_check{$input_id->{ID}} = 1;
        }
        if (!exists $input_id->{is_mate_1} or $input_id->{'is_mate_1'} == -1) {
            my $regex = $self->param('pairing_regex');
            my ($pair) = $input_id->{filename} =~ /$regex/;
            if (($pair == 1) or ($input_id->{is_paired} == 0)) {
                $input_id->{is_mate_1} = 1;
            }
            else {
                $input_id->{is_mate_1} = 0;
            }
        }

#get read_length from the read_length table
	my $length_table_adaptor = $self->db->get_NakedTableAdaptor;
	$length_table_adaptor->table_name($self->param('read_length_table'));

	my $split_fastq = $input_id->{filename};
	$split_fastq =~ m/([A-Z0-9]+)_.*([0-9]\.fastq\.gz)/;
	my $fastq = $1."_".$2;

  if($1 eq ""){
    $split_fastq =~ m/([a-zA-Z0-9]+)_.*([0-9]\.fq\.gz)/;
    $fastq = $1."_".$2;
  }
	my $db_row = $length_table_adaptor->fetch_by_dbID($fastq);
	my $read_length = $db_row->{read_length};
  print("DEBUG:: $fastq \t $read_length \t $split_fastq \t $db_row->{DS}\n");
	$input_id->{read_length} = $read_length;
#finish get read_length from table

        foreach my $key (keys %$input_id) {
          $input_id->{$key} =~ tr /:\t/ /;
        }

        my $table_adaptor = $self->db->get_NakedTableAdaptor;
        $table_adaptor->table_name($self->param('csvfile_table'));
        $table_adaptor->store([$input_id]);
        $keyword_hash{$input_id->{$self->param('sample_column')}} = 1;
    }
    my @output_ids;
    foreach my $key (keys %keyword_hash) {
        push(@output_ids, { sample_name => $key });
    }
    $self->dataflow_output_id(\@output_ids, $self->param('fan_branch_code'));
}

1;
