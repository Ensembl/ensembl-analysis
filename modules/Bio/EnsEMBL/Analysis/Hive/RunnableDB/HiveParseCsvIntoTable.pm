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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveParseCsvIntoTable;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');

sub write_output {
    my $self = shift;

    my %keyword_hash;
    my $table_adaptor = $self->db->get_NakedTableAdaptor;
    $table_adaptor->table_name($self->param('csvfile_table'));
    foreach my $input_id (@{$self->param('output_ids')}) {
        $input_id->{$self->param('sample_column')} =~ s/ /_/g;
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
