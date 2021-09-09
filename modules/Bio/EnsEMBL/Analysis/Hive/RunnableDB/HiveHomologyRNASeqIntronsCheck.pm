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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveHomologyRNASeqIntronsCheck

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveHomologyRNASeqIntronsCheck;

use strict;
use warnings;
use feature 'say';

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene fully_load_Gene);
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults()},
    target_db => undef,
    update_genes => 0,
    full_support_suffix => 'top', # We will concat this value to the previous biotype
    classify_weak=> 1,
    classify_low => 2,
    classify_medium => 4,
    classify_high => 6,
    classify_top => 8,
    slice_strand => 0,
  }
}


sub fetch_input {
  my ($self) = @_;

  $self->create_analysis;
  my $source_db = $self->get_database_by_name('source_db');
  my $intron_db = $self->get_database_by_name('intron_db');
  my $target_db;
  if ($self->param('target_db')) {
     $target_db = $self->get_database_by_name('target_db');
  }

  if($self->param('use_genome_flatfile')) {
    say "Ingoring dna table and using fasta file for sequence fetching";
    unless($self->param_required('genome_file') && -e $self->param('genome_file')) {
      $self->throw("You selected to use a flatfile to fetch the genome seq, but did not find the flatfile. Path provided:\n".$self->param('genome_file'));
    }
    setup_fasta(
                 -FASTA => $self->param_required('genome_file'),
               );
  } else {
    say "Attaching dna db";
    my $dna_dba = $self->get_database_by_name('dna_db');
    $source_db->dnadb($dna_dba);
    $intron_db->dnadb($dna_dba);
    if($target_db) {
      $target_db->dnadb($dna_dba);
    }
  }

  if($target_db) {
    $target_db->dbc->disconnect_when_inactive(1);
    $self->hrdb_set_con($target_db,'target_db');
  } else {
    $self->hrdb_set_con($source_db, 'target_db');
    $self->param('update_genes', 1); # If people didn't specify a target_db, we update the biotype, otherwise they have to specify a target_db
  }

  $self->hrdb_set_con($intron_db,'intron_db');

  my $slice = $source_db->get_SliceAdaptor->fetch_by_name($self->input_id);
  my $genes;
  if ($self->param_is_defined('source_logic_name')) {
    if (ref($self->param('source_logic_name')) eq 'ARRAY') {
      foreach my $logic_name (@{$self->param('source_logic_name')}) {
        push(@$genes, @{$slice->get_all_Genes($logic_name)});
      }
    }
    else {
      $genes = $slice->get_all_Genes($self->param('source_logic_name'));
    }
  }
  else {
    $genes = $slice->get_all_Genes();
  }

  # If we have the slice_strand param set then we want to filter the genes based on the strand
  if($self->param('slice_strand')) {
    my $initial_genes = $genes;
    $genes = [];
    foreach my $gene (@{$initial_genes}) {
      unless($self->param('slice_strand') == $gene->strand) {
        next;
      }
      push(@$genes,$gene);
    }
  }

  print STDERR 'Fetched ', scalar(@$genes), "\n";
  $self->say_with_header('Fetched '.scalar(@$genes));
  $self->input_genes($genes);
}


sub run {
  my ($self) = @_;

  my $intron_adaptor = $self->hrdb_get_con('intron_db')->get_DnaAlignFeatureAdaptor;
  my %good_introns;
  my $full_support_suffix = $self->param('full_support_suffix');

  foreach my $gene (@{$self->input_genes}) {
    my %introns;
    my $transcripts = $gene->get_all_Transcripts;
    if(scalar(@$transcripts) > 1) {
      $self->throw("The module is not currently designed to work with multi-transcript genes");
    }

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

      if($self->param('classify_by_count')) {
        my $weak_count   = $self->param_required('classify_weak');
        my $low_count    = $self->param_required('classify_low');
        my $medium_count = $self->param_required('classify_medium');
        my $high_count   = $self->param_required('classify_high');
        my $top_count    = $self->param_required('classify_top');
        my $intron_support_diff = $good_intron_count - ($intron_count - $good_intron_count);
        if($intron_support_diff >= $top_count) {
          $transcript->biotype('genblast_rnaseq_top');
          $gene->biotype('genblast_rnaseq_top');
          $self->output([$gene]);
        } elsif($intron_support_diff >= $high_count) {
          $transcript->biotype('genblast_rnaseq_high');
          $gene->biotype('genblast_rnaseq_high');
          $self->output([$gene]);
        } elsif($intron_support_diff >= $medium_count) {
          $transcript->biotype('genblast_rnaseq_medium');
          $gene->biotype('genblast_rnaseq_medium');
          $self->output([$gene]);
        } elsif($intron_support_diff >= $low_count) {
          $transcript->biotype('genblast_rnaseq_low');
          $gene->biotype('genblast_rnaseq_low');
          $self->output([$gene]);
        } elsif($intron_support_diff >= $weak_count) {
          $transcript->biotype('genblast_rnaseq_weak');
          $gene->biotype('genblast_rnaseq_weak');
          $self->output([$gene]);
        } else {
          say "Skipping gene as there is no support";
        }
      } else {
        if ($intron_count == $good_intron_count) {
          $self->say_with_header('Gene '.$gene->display_id.' updated');
          $gene->biotype($transcript->biotype.'_'.$full_support_suffix);
          $transcript->biotype($transcript->biotype.'_'.$full_support_suffix);
        }
        else {
          $self->say_with_header('Gene '.$gene->display_id.' not fully supported');
        }
        $self->output([$gene]);
      }

    } # End foreach my $transcript
  } # End foreach my $gene
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

sub input_genes {
  my ($self,$genes) = @_;

  if($genes) {
    $self->param('_input_genes',$genes);
  }

  return($self->param('_input_genes'));
}

1;
