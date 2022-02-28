#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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


=head2 param_defaults

 Arg [1]    : None
 Description: Returns the default parameters:
               _ln_gene_ext => 'gene',
               _ln_introns_ext => 'daf',
 Returntype : Hash ref
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    _ln_gene_ext => 'gene',
    _ln_introns_ext => 'daf',
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Creates input id based on a custom table 'csvfile_table' in the hive database
              It will generate the parameters for creating the models with the RefineSolexaGenes
              module. If you specify 'single_tissue' to 1 it will generate input ids for the merged
              set AND for each samples in 'sample_column'
              It stores the input ids in 'inputlist'
 Returntype : None
 Exceptions : None

=cut

sub fetch_input {
    my $self = shift;

    my @output_ids;
    my $gene_ext = $self->param('_ln_gene_ext');
    my $introns_ext = $self->param('_ln_introns_ext');
    if ($self->param('single_tissue')) {
        my $table_adaptor = $self->db->get_NakedTableAdaptor;
        $table_adaptor->table_name($self->param('csvfile_table'));
        my %tissue_hash;
        my $results = $table_adaptor->fetch_all();
        foreach my $result (@$results) {
            $tissue_hash{$result->{$self->param('sample_column')}}->{$result->{$self->param('sample_id_column')}} = 1;
        }
        foreach my $key (keys %tissue_hash) {
            push(@output_ids, [$self->param('iid'), [{file => $self->param('wide_intron_bam_file').'.bam', groupname => [keys %{$tissue_hash{$key}}], depth => 0, mixed_bam => 0}], $self->param('wide_species').'_'.$key.'_rnaseq_'.$gene_ext, $self->param('wide_species').'_'.$key.'_rnaseq_'.$introns_ext, "best_$key", "single_$key", '', '']);
        }
    }
    push(@output_ids, [$self->param('iid'), [{file => $self->param('wide_intron_bam_file').'.bam', groupname => [], depth => 0, mixed_bam => 0}], $self->param('wide_species').'_merged_rnaseq_'.$gene_ext, $self->param('wide_species').'_merged_rnaseq_'.$introns_ext, "best", "single", '', '']);
    $self->param('inputlist', \@output_ids);
    $self->param('column_names', ['iid', 'intron_bam_files', 'logic_name', 'introns_logic_name', 'best_score', 'single_exon_model', 'other_isoforms', 'bad_models']);
}

1;
