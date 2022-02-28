=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::CleanUTRs

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::CleanUTRs;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene clean_utrs);

use parent qw(Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB);


sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    min_size_utr_exon => 30,
    ratio_5prime_utr => .3,
    ratio_3prime_utr => .6,
    ratio_same_transcript => .02,
    ratio_max_allowed_difference => .05,
    ratio_expansion => 3,
    minimum_expanding_number_for_single_transcript => 2,
    ratio_transcript_fragment => 3,
    ratio_exon_expansion => 2,
    ratio_utrs => 2,
    store_rejected => 0,
    copy_biotypes_to_ignore => {
      low_coverage => 1,
      CRISPR => 1,
      broken_gene => 1
    }
  }
}



sub fetch_input {
  my ($self) = @_;

  my $db = $self->get_database_by_name('source_db');
  if ($self->param_is_defined('target_db')) {
    $self->hrdb_set_con($self->get_database_by_name('target_db'), 'target_db');
  }
  else {
    $self->hrdb_set_con($db, 'target_db');
  }
  if ($self->param_is_defined('dna_db')) {
    my $dna_db = $self->get_database_by_name('dna_db');
    $db->dnadb($dna_db);
  }

  my $slice = $self->fetch_sequence($self->input_id, $db);
# We store the genes directly in output as we will store any genes but the transcripts will be modified
  my @genes;
  my @protein_coding;
  my $unwanted = $self->param('copy_biotypes_to_ignore');
  foreach my $gene (@{$slice->get_all_Genes}) {
    if (@{$gene->get_all_Transcripts}) {
      if (!exists $unwanted->{$gene->biotype}) {
        $gene->load;
        push(@genes, $gene);
        if ($gene->biotype eq 'protein_coding') {
          push(@protein_coding, $gene);
        }
      }
    }
  }
  if (@genes) {
    $self->output(\@genes);
    $self->param('protein_coding_genes', \@protein_coding);
  }
  else {
    $self->complete_early('No genes found for '.$self->input_id);
  }
}

sub run {
  my ($self) = @_;

  my ($extra_genes, $rejected) = clean_utrs($self->param('protein_coding_genes'), $self->param('min_size_utr_exon'),
          $self->param('ratio_5prime_utr'), $self->param('ratio_3prime_utr'), $self->param('ratio_same_transcript'),
          $self->param('ratio_max_allowed_difference'), $self->param('ratio_expansion'), $self->param('minimum_expanding_number_for_single_transcript'),
          $self->param('ratio_transcript_fragment'), $self->param('ratio_exon_expansion'), $self->param('ratio_utrs'), $self->param('store_rejected'));

  $self->output($extra_genes);
  $self->output($rejected) if ($self->param('store_rejected'));
}

sub write_output {
  my ($self) = @_;

  my $gene_adaptor = $self->hrdb_get_con('target_db')->get_GeneAdaptor;
  my $suffix = $self->param('suffix');
  foreach my $gene (@{$self->output}) {
    if ($suffix) {
      my $biotype = $gene->biotype."_$suffix";
      $gene->biotype($biotype);
      foreach my $transcript (@{$gene->get_all_Transcripts}) {
        $transcript->biotype($biotype);
      }
    }
    empty_Gene($gene);
    $gene_adaptor->store($gene);
  }
}

1;
