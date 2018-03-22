=head1 LICENSE

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveHomologyRNASeqIntronsCheck

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveHomologyRNASeqIntronsCheck;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene fully_load_Gene);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults()},
    target_db => undef,
    update_genes => 0,
    full_support_suffix => 'top', # We will concat this value to the previous biotype
  }
}


sub fetch_input {
  my ($self) = @_;

  $self->create_analysis;
  my $dna_db = $self->get_database_by_name('dna_db');
  my $source_db = $self->get_database_by_name('source_db', $dna_db);
  if ($self->param('target_db')) {
    $self->hrdb_set_con($self->get_database_by_name('target_db', $dna_db), 'target_db');
  }
  else {
    $self->hrdb_set_con($source_db, 'target_db');
    $self->param('update_genes', 1); # If people didn't specify a target_db, we update the biotype, otherwise they have to specify a target_db
  }
  my $slice = $source_db->get_SliceAdaptor->fetch_by_name($self->input_id);
  my $genes;
  if ($self->param_is_defined('source_logic_name')) {
    if (ref($self->param('source_logic_name')) eq 'ARRAY') {
      foreach my $logic_name (@{$self->param('source_logic_name')}) {
        push(@$genes, @{$slice->get_all_Genes($logic_name, undef, 1)});
      }
    }
    else {
      $genes = $slice->get_all_Genes($self->param('source_logic_name'), undef, 1);
    }
  }
  else {
    $genes = $slice->get_all_Genes(undef, undef, 1);
  }
  $self->hrdb_set_con($self->get_database_by_name('intron_db', $dna_db), 'intron_db');
  print STDERR 'Fetched ', scalar(@$genes), "\n";
  $self->output($genes);
}


sub run {
  my ($self) = @_;

  my $intron_adaptor = $self->hrdb_get_con('intron_db')->get_DnaAlignFeatureAdaptor;
  my %good_introns;
  my $full_support_suffix = $self->param('full_support_suffix');
  foreach my $gene (@{$self->output}) {
    my %introns;
    foreach my $transcript (@{$gene->get_all_Transcripts}) {
      my $intron_count = 0;
      my $good_intron_count = 0;
      foreach my $intron (@{$transcript->get_all_Introns}) {
        ++$intron_count;
        my $hashkey = sprintf("%d:%d:%d", $intron->seq_region_start, $intron->seq_region_end, $intron->strand);
        if (exists $good_introns{$hashkey}) {
          ++$good_intron_count;
        }
        elsif (exists $introns{$hashkey}) {
          ++$good_intron_count;
          $good_introns{$hashkey} = $introns{$hashkey};
        }
        elsif (scalar(keys %introns) == 0) {
          foreach my $intron (@{$intron_adaptor->fetch_all_by_Slice($gene->feature_Slice)}) {
            $introns{sprintf("%d:%d:%d", $intron->seq_region_start, $intron->seq_region_end, $intron->strand*$gene->strand)} = $intron->score;
          }
          if (exists $introns{$hashkey}) {
            ++$good_intron_count;
            $good_introns{$hashkey} = $introns{$hashkey};
          }
        }
      }
      if ($intron_count == $good_intron_count) {
        print STDERR 'Gene ', $gene->display_id, " updated\n";
        $gene->biotype($transcript->biotype.'_'.$full_support_suffix);
        $transcript->biotype($transcript->biotype.'_'.$full_support_suffix);
      }
      else {
        print STDERR 'Gene ', $gene->display_id, ' not fully supported', "\n";
      }
    }
  }
}


sub write_output {
  my ($self) = @_;

  my $analysis = $self->analysis;
  my $gene_adaptor = $self->hrdb_get_con('target_db')->get_GeneAdaptor;
  my $update_genes = $self->param('update_genes');
  foreach my $gene (@{$self->output}) {
    if ($update_genes) {
      foreach my $transcript (@{$gene->get_all_Transcripts}) {
        $transcript->adaptor->update($transcript);
      }
      $gene->adaptor->update($gene);
    }
    else {
      fully_load_Gene($gene);
      empty_Gene($gene);
      $gene->analysis($analysis);
      $gene_adaptor->store($gene);
    }
  }
}

1;
