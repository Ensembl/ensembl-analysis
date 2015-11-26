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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateRefineGenesJobs;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');

sub fetch_input {
    my $self = shift;

    my $table_adaptor = $self->db->get_NakedTableAdaptor;
    $table_adaptor->table_name($self->param('csvfile_table'));
    my @output_ids;
    my %tissue_hash;
    my $results = $table_adaptor->fetch_all();
    foreach my $result (@$results) {
        push(@{$tissue_hash{$result->{$self->param('tissue_name')}}}, $result->{$self->param('sample_id')});
    }
    foreach my $key (keys %tissue_hash) {
        push(@output_ids, [$self->param('iid'), [file => $self->param('wide_intron_bam_file'), groupname => $tissue_hash{$key}, depth => 0, mixed_bam => 0], $self->param('wide_species').'_'.$key.'_rnaseq', $self->param('wide_species').'_'.$key.'_introns', "best_$key", "single_$key", '', '']);
    }
    push(@output_ids, [$self->param('iid'), [file => $self->param('wide_intron_bam_file'), groupname => [], depth => 0, mixed_bam => 0], $self->param('wide_species').'_merged_rnaseq', $self->param('wide_species').'_merged_introns', "best", "single", '', '']);
    $self->param('inputlist', \@output_ids);
    $self->param('column_names', ['iid', 'intron_bam_files', 'logic_names', 'introns_logic_name', 'best_score', 'single_exon_model', 'other_isoforms', 'bad_models']);
}

1;
