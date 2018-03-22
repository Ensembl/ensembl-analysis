#!/usr/bin/env perl

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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
    foreach my $input_id (@{$self->param('output_ids')}) {
        $input_id->{$self->param('sample_column')} =~ tr/ :\t/_/;
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
        foreach my $key (keys %$input_id) {
          $input_id->{$key} =~ tr /:\t/ /;
        }
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
