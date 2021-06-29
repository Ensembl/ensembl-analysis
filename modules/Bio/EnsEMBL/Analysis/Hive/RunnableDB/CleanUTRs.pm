=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils qw(make_types_hash cluster_Genes);

use parent qw(Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB);


sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    min_size_utr_exon => 30,
    ratio_5prime_utr => .3,
    ratio_3prime_utr => .6,
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

  my $slice = $self->fetch_sequence($self->input_id, $db);
# We store the genes directly in output as we will store any genes but the transcripts will be modified
  my @genes;
  my @protein_coding;
  foreach my $gene (@{$slice->get_all_Genes}) {
    if (@{$gene->get_all_Transcripts}) {
      $gene->load;
      push(@genes, $gene);
      if ($gene->biotype eq 'protein_coding') {
        push(@protein_coding, $gene);
      }
    }
    else {
      $self->say_with_header($gene->display_id.' has no transcript');
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

  my $min_size_utr_exon = $self->param('min_size_utr_exon');
  my $ratio_5prime_utr = $self->param('ratio_5prime_utr');
  my $ratio_3prime_utr = $self->param('ratio_3prime_utr');
  my ($clusters, $unclustered) = cluster_Genes($self->param('protein_coding_genes'), make_types_hash($self->param('protein_coding_genes'), undef, 'set1'));
  foreach my $cluster (@$clusters) {
    my @overlapping_genes = sort {$a->start <=> $b->start || $a->end <=> $b->end} @{$cluster->get_Genes_by_Set('set1')};
    for (my $gene_index = 0; $gene_index <= $#overlapping_genes; $gene_index++) {
      my $gene = $overlapping_genes[$gene_index];
      $self->say_with_header("Working on $gene_index ".$gene->display_id);
      my $transcripts = $gene->get_all_Transcripts;
      for (my $next_gene_index = 0; $next_gene_index <= $#overlapping_genes; $next_gene_index++) {
        my $next_gene = $overlapping_genes[$next_gene_index];
        if ($gene_index != $next_gene_index) {
          $self->say_with_header("Checking $next_gene_index ".$next_gene->display_id);
          if ($gene->overlaps_local($next_gene)) {
            $self->say_with_header('Comparing '.$gene->display_id.' '.$next_gene->display_id);
            my $change_happened = 0;
            foreach my $transcript (@$transcripts) {
              if ($transcript->overlaps_local($next_gene)) {
                my $cds_start_genomic = $transcript->coding_region_start;
                my $cds_end_genomic = $transcript->coding_region_end;
                my $cds_start_index = 0;
                my $cds_end_index = 0;
                $self->say_with_header($transcript->display_id.' overlaps '.$next_gene->display_id);
                my %overlapping_exons;
                my %exons_to_delete;
                my @exons = sort {$a->start <=> $b->start} @{$transcript->get_all_Exons};
                my $count = 0;
                foreach my $exon (@exons) {
                  if ($cds_start_genomic >= $exon->start and $cds_start_genomic <= $exon->end) {
                    $cds_start_index = $count;
                  }
                  if ($cds_end_genomic >= $exon->start and $cds_end_genomic <= $exon->end) {
                    $cds_end_index = $count;
                  }
                  if ($exon->overlaps_local($next_gene)) {
                    $overlapping_exons{$exon->start.':'.$exon->end} = $exon;
                    $self->say_with_header('OVERLAP '.join(' ', $exon->display_id, $exon->start, $exon->end));
                  }
                  ++$count;
                }
                if (scalar(keys %overlapping_exons)) {
                  $self->say_with_header("$cds_start_genomic $cds_start_index $cds_end_genomic $cds_end_index");
                  my $new_utr_exon_start = 0;
                  my $new_utr_exon_end = 0;
                  foreach my $next_transcript (@{$next_gene->get_all_Transcripts}) {
                    $self->say_with_header('Comparing '.$transcript->display_id.' against '.$next_transcript->display_id);
                    foreach my $cds_exon (@{$next_transcript->get_all_CDS}) {
                      foreach my $utr_exon (values %overlapping_exons) {
                        if ($cds_exon->overlaps_local($utr_exon)) {
                          $self->say_with_header($cds_exon->start.':'.$cds_exon->end.' is overlapped by '.$utr_exon->display_id);
                          $exons_to_delete{$utr_exon->start.':'.$utr_exon->end} = $utr_exon;
                          if (!exists $overlapping_exons{$cds_exon->start.':'.$cds_exon->end}) {
                            if ($utr_exon->start <= $cds_start_genomic and $utr_exon->end >= $cds_start_genomic
                                and $utr_exon->start <= $next_transcript->coding_region_end and $utr_exon->end >= $next_transcript->coding_region_end) {
                              $new_utr_exon_start = $cds_exon->end;
                            }
                            if ($utr_exon->start <= $cds_end_genomic and $utr_exon->end >= $cds_end_genomic
                                and $utr_exon->start <= $next_transcript->coding_region_start and $utr_exon->end >= $next_transcript->coding_region_start) {
                              $new_utr_exon_end = $cds_exon->start;
                            }
                            $self->say_with_header($cds_exon->start.':'.$cds_exon->end." is not in hash $new_utr_exon_start $new_utr_exon_end");
                          }
                        }
                      }
                    }
                  }
                  if (scalar(keys %exons_to_delete)) {
                    $change_happened = 1;
                    my $translation;
                    if ($new_utr_exon_end and $transcript->strand == -1) {
                      $translation = $transcript->translation;
                    }
                    elsif ($new_utr_exon_start and $transcript->strand == 1) {
                      $translation = $transcript->translation;
                    }
                    $transcript->flush_Exons;
                    my $start_index = 0;
                    my $end_index = $#exons;
                    $self->say_with_header("$end_index exons to start with");
                    $self->say_with_header("Working on 5' end");
                    for (my $index = $cds_start_index; $index >= 0; $index--) {
                      $self->say_with_header("$cds_start_index $index");
                      if (exists $exons_to_delete{$exons[$index]->start.':'.$exons[$index]->end}) {
                        if ($exons_to_delete{$exons[$index]->start.':'.$exons[$index]->end}->start <= $cds_start_genomic and $exons_to_delete{$exons[$index]->start.':'.$exons[$index]->end}->end >= $cds_start_genomic) {
                          $start_index = $index;
                        }
                        last;
                      }
                      else {
                        $start_index = $index;
                      }
                    }
                    $self->say_with_header("Working on 3' end");
                    for (my $index = $cds_end_index; $index <= $#exons; $index++) {
                      $self->say_with_header("$cds_end_index $index");
                      if (exists $exons_to_delete{$exons[$index]->start.':'.$exons[$index]->end}) {
                        if ($exons_to_delete{$exons[$index]->start.':'.$exons[$index]->end}->start <= $cds_end_genomic and $exons_to_delete{$exons[$index]->start.':'.$exons[$index]->end}->end >= $cds_end_genomic) {
                          $end_index = $index;
                        }
                        last;
                      }
                      else {
                        $end_index = $index;
                      }
                    }
                    $self->say_with_header("Result: $start_index $end_index");
                    foreach my $exon (@exons[$start_index..$end_index]) {
                      if ($new_utr_exon_start >= $exon->start and $new_utr_exon_start <= $exon->end) {
                        if ($exon->strand == -1) {
                          $new_utr_exon_start = $cds_start_genomic-int(($cds_start_genomic-$exon->start)*$ratio_3prime_utr);
                        }
                        else {
                          $new_utr_exon_start = $cds_start_genomic-int(($cds_start_genomic-$exon->start)*$ratio_5prime_utr);
                        }
                        $self->say_with_header($new_utr_exon_start);
                        if ($cds_start_genomic-$new_utr_exon_start < $min_size_utr_exon) {
                          $new_utr_exon_start = $cds_start_genomic-$min_size_utr_exon;
                        }
                        $self->say_with_header($new_utr_exon_start);
                        if ($new_utr_exon_start < $exon->start) {
                          $new_utr_exon_start = $exon->start;
                        }
                        $self->say_with_header($new_utr_exon_start);
                        $self->say_with_header($translation || 'NULL');
                        $translation->start($cds_start_genomic-$new_utr_exon_start+1) if ($translation);
                        $exon->start($new_utr_exon_start);
                      }
                      if ($new_utr_exon_end >= $exon->start and $new_utr_exon_end <= $exon->end) {
                        if ($exon->strand == -1) {
                          $new_utr_exon_end = $cds_end_genomic+int(($exon->end-$cds_end_genomic)*$ratio_5prime_utr);
                        }
                        else {
                          $new_utr_exon_end = $cds_end_genomic+int(($exon->end-$cds_end_genomic)*$ratio_3prime_utr);
                        }
                        $self->say_with_header($new_utr_exon_end);
                        if ($new_utr_exon_end-$cds_end_genomic < $min_size_utr_exon) {
                          $new_utr_exon_end = $cds_end_genomic+$min_size_utr_exon;
                        }
                        $self->say_with_header($new_utr_exon_end);
                        if ($new_utr_exon_end > $exon->end) {
                          $new_utr_exon_end = $exon->end;
                        }
                        $self->say_with_header($new_utr_exon_end);
                        $self->say_with_header($translation || 'NULL');
                        $translation->start($new_utr_exon_end-$cds_end_genomic+1) if ($translation);
                        $exon->end($new_utr_exon_end);
                      }
                      $transcript->add_Exon($exon);
                    }
                    if ($translation and $translation->start_Exon == $translation->end_Exon) {
                      if ($transcript->strand == -1) {
                        $translation->end($translation->start_Exon->end-$cds_start_genomic+1);
                      }
                      else {
                        $translation->end($cds_end_genomic-$translation->start_Exon->start+1);
                      }
                    }
                  }
                }
              }
            }
            if ($change_happened) {
              my %hashes;
              my $transcripts = $gene->get_all_Transcripts;
              foreach my $transcript (@$transcripts) {
                my $id = '';
                foreach my $exon (@{$transcript->get_all_Exons}) {
                  $id .= join(':', $exon->start, $exon->end, $exon->phase, $exon->end_phase);
                }
                if ($transcript->translation) {
                  $id .= $transcript->translation->start.':'.$transcript->translation->end;
                }
                push(@{$hashes{$id}}, $transcript);
              }
              if (scalar(keys %hashes) != @$transcripts) {
                $gene->flush_Transcripts;
                foreach my $item (values %hashes) {
                  $gene->add_Transcript($item->[0]);
                }
              }
              else {
                $gene->recalculate_coordinates;
              }
              $self->throw($gene->display_id.' has no transcript') unless (@{$gene->get_all_Transcripts});
            }
          }
        }
      }
    }
  }
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
