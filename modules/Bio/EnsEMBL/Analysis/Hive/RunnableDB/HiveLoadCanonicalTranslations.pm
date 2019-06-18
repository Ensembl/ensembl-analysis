# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

=head1 NAME 

HiveCopyGenes.pm

=head1 DESCRIPTION

This module copies takes a gene id or set of gene ids as input and loads their canonical translation into the protein table.

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadCanonicalTranslations;

use strict;
use warnings;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub fetch_input {
  my $self = shift;
  my $input_dba = $self->hrdb_get_dba($self->param('target_db'));
  $self->hrdb_set_con($input_dba,'target_db');
  return 1;
}

sub run {

  my $self = shift;

  my $input_dba = $self->hrdb_get_con('target_db');
  my $input_ga = $input_dba->get_GeneAdaptor();
  my $output_transcripts = [];
#  my $input_gene_ids = $self->param('iid');
#  foreach my $gene_id (@{$input_gene_ids}) {
  my $gene = $input_ga->fetch_by_dbID($self->param('iid'));
  my $canonical_transcript = $gene->canonical_transcript();
  unless($canonical_transcript) {
    $self->throw("A canonical transcript could not be found for gene with dbID: ".$self->param('iid'));
  }

  if($canonical_transcript->translate) {
    push(@{$output_transcripts},$canonical_transcript);
    $self->output_transcripts($output_transcripts);
  }

#    unless($canonical_transcript->translate) {
#      next;
#    }
#    push(@{$output_transcripts},$canonical_transcript);
#  }
#  $self->output_transcripts($output_transcripts);

  return 1;
}

sub write_output {
  my $self = shift;

  my $table_adaptor = $self->db->get_NakedTableAdaptor();
  $table_adaptor->table_name('protein_sequences');

  my $input_dba = $self->hrdb_get_con('target_db');
  my $db_name = $input_dba->dbname();
  my $output_transcripts = $self->output_transcripts();

  foreach my $transcript (@{$output_transcripts}) {
    my $accession = $transcript->stable_id();
    my $seq = $transcript->translate->seq();
    my $db_row = [{ 'accession'  => $accession,
                    'source'     => $db_name,
                    'seq'        => $seq,
                   }];
    $table_adaptor->store($db_row);

    my $output_hash = {};
    $output_hash->{'iid'} = [$accession];
    $self->dataflow_output_id($output_hash,1);
  }

  return 1;
}


sub output_transcripts {
  my ($self,$output_genes) = @_;

  if($output_genes) {
    $self->param('_output_genes',$output_genes);
  }

  return($self->param('_output_genes'));
}

1;
