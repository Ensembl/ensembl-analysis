=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveInternalStopFix

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveInternalStopFix;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable::InternalStopFix;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
    my $self = shift;

    my $dnadb = $self->get_database_by_name('dna_db');
    my $db = $self->hrdb_set_con($self->get_database_by_name('source_db', $dnadb));
    if ($self->param_is_defined('target_db')) {
      $self->hrdb_set_con($self->get_database_by_name('source_db', $dnadb), 'target_db');
    }
    else {
      $self->hrdb_set_con($db, 'target_db');
    }
    my $slice = $db->get_SliceAdaptor->fetch_by_name($self->input_id);
    my $genes = $slice->get_all_Genes($self->LOGIC_NAME, undef, 1, $self->SOURCE, $self->BIOTYPE);
    if (scalar(@$genes)) {
      $self->analysis($genes->[0]->analysis);
      my $runnable = Bio::EnsEMBL::Analysis::Runnable::InternalStopFix->new(
          -analysis => $self->analysis,
          -stop_codon_biotype => $self->STOP_CODON_BIOTYPE,
          -edited_biotype => $self->EDITED_BIOTYPE,
          -genes => $genes,
      );
      $self->runnable($runnable);
    }
    else {
      $self->complete_early('No genes to process');
    }
}

sub write_output {
    my $self = shift;

    my $gene_adaptor = $self->hrdb_get_con('target_db')->get_GeneAdaptor();
    my $update_gene_adaptor = $self->hrdb_get_con->get_GeneAdaptor();
    my $update_transcript_adaptor = $self->hrdb_get_con->get_TranscriptAdaptor();
    foreach my $gene (@{$self->output}) {
        if ($gene->dbID) {
            $update_gene_adaptor->update($gene);
            foreach my $t (@{$gene->get_all_Transcripts}) {
              $update_transcript_adaptor->update($t);
            }
        }
        else {
            $gene_adaptor->store($gene);
        }
    }
}

sub EDITED_BIOTYPE {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->param('edited_biotype', $val);
  }
  if ($self->param_is_defined('edited_biotype')) {
    return $self->param('edited_biotype');
  }
  else {
    return;
  }
}

sub STOP_CODON_BIOTYPE {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->param('stop_codon_biotype', $val);
  }
  if ($self->param_is_defined('stop_codon_biotype')) {
    return $self->param('stop_codon_biotype');
  }
  else {
    return;
  }
}

sub LOGIC_NAME {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->param('logic_name', $val);
  }
  if ($self->param_is_defined('logic_name')) {
    return $self->param('logic_name');
  }
  else {
    return;
  }
}

sub SOURCE {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->param('source', $val);
  }
  if ($self->param_is_defined('source')) {
    return $self->param('source');
  }
  else {
    return;
  }
}

sub BIOTYPE {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->param('biotype', $val);
  }
  if ($self->param_is_defined('biotype')) {
    return $self->param('biotype');
  }
  else {
    return;
  }
}

1;

