=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRNASeqGenes

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRNASeqGenes;

use strict;
use warnings;

use Bio::DB::HTS;

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(convert_to_ucsc_name);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    intron_max_distance => 15000,
    max_exon_size => 4000,
  }
}



sub fetch_input {
  my( $self) = @_;

  my $input_id = $self->input_id;
# I want to be able to use either slices or gene stable ids
  my $dna_db = $self->get_database_by_name('dna_db');
  $dna_db->dbc->disconnect_when_inactive(0);
  my $genes_db = $self->get_database_by_name('input_db', $dna_db);
  $genes_db->dbc->disconnect_when_inactive(0);
  $self->hrdb_set_con($self->get_database_by_name('output_db', $dna_db), 'output_db');
  my @rough_genes;
  my $chr_slice;
  my $dna_slice;
  if ($self->is_slice_name($input_id)) {
    my $genes;
    my $slice = $self->fetch_sequence($input_id, $genes_db);
    my $real_slice_start = $slice->start;
    my $real_slice_end = $slice->end;
    $chr_slice = $dna_db->get_SliceAdaptor->fetch_by_region( 'toplevel', $slice->seq_region_name);

    if ($self->param('model_ln')) {
      $genes = $slice->get_all_Genes_by_type( undef,$self->param('model_ln') );
      print STDERR "Got " .  scalar(@$genes) . " genes with logic name " . $self->param('model_ln') ."\n";
    }
    else {
      $genes = $slice->get_all_Genes(undef, undef, 1);
      print STDERR "Got " .  scalar(@$genes) . "  genes  \n";
    }
    foreach my $gene (@$genes) {
# put them on the chromosome
      $gene = $gene->transfer($chr_slice);
# reject genes that are from a different slice that overlap our slice at the start or end
# say the models has to be > 10% on the slice
      my $os = $slice->start;
      $os = $gene->start if($gene->start > $slice->start);
      my $oe = $slice->end;
      $oe = $gene->end if($gene->end < $slice->end);
      if ($slice->start != $os or $gene->end != $oe) {
        my $overlap = $oe - $os +1;
        my $gc = int(($overlap / $gene->length) * 1000) / 10;
        my $sc =  int(($overlap / $slice->length) * 1000) /10;
        if ( $gc <= 10 && $sc <= 10) {
          print "Gene ", $gene->display_id, " has $gc% overlap with the slice\nSlice has $sc% overlap with the gene\n  Rejecting\n";
          next;
        }
        $real_slice_start = $gene->start < $real_slice_start ? $gene->start : $real_slice_start;
        $real_slice_end = $gene->end > $real_slice_end ? $gene->end : $real_slice_end;
      }
      $gene->load;
      push(@rough_genes,$gene);
    }
    $dna_slice = $slice->sub_Slice($real_slice_start, $real_slice_end);
    print STDERR "Got " . scalar(@rough_genes) . "  genes after filtering boundary overlaps  \n";
  }
  else {
    my $gene = $genes_db->get_GeneAdaptor->fetch_by_stable_id($input_id);
    $chr_slice = $dna_db->get_SliceAdaptor->fetch_by_region( 'toplevel', $gene->slice->seq_region_name);
    $dna_slice = $gene->slice->sub_Slice($gene->seq_region_start, $gene->seq_region_end);
    $gene->load;
    push(@rough_genes, $gene);
  }
  $self->chr_slice($chr_slice);
  if (scalar(@rough_genes)) {
    $self->create_analysis;
    $self->load_introns($dna_slice);
    $genes_db->dbc->disconnect_when_inactive(1) if ($self->param('disconnect_jobs'));
  }
  else {
    $self->input_job->autoflow(0);
    $self->complete_early('No genes to process');
  }
}

sub run {
  my ($self) = @_;

  my $intron_clusters = $self->param('intron_clusters');
  my $slice = $self->chr_slice;
  my $forward_intron_clusters = $intron_clusters->{1};
  if (@$forward_intron_clusters) {
    my @score_sorted_introns = sort {$b->{depth} <=> $a->{depth}} @$forward_intron_clusters;
  }
  my $reverse_intron_clusters = $intron_clusters->{-1};
  if (@$reverse_intron_clusters) {
    my @depth_sorted_introns = sort {$b->{depth} <=> $a->{depth}} @$reverse_intron_clusters;
    foreach my $intron (@depth_sorted_introns) {
      print __LINE__, ' ', $intron->{exon_5_start}, ' ', $intron->{start}, ' ', $intron->{end}, ' ', $intron->{exon_3_end}, ' ', $intron->{depth}, ' ', $intron->{is_canonical}, ' ', scalar(@{$intron->{extra_exon}}), "\n";
      foreach my $key ('left_canonical', 'left_non_canonical', 'right_canonical', 'right_non_canonical') {
        if (exists $intron->{$key}) {
          foreach my $link (@{$intron->{$key}}) {
            print __LINE__, '   ', $key, ' ', $link->{exon_5_start}, ' ', $link->{start}, ' ', $link->{end}, ' ', $link->{exon_3_end}, ' ', $link->{depth}, ' ', $link->{is_canonical}, ' ', scalar(@{$link->{extra_exon}}), ' ', ($key =~ /right/ ? ($link->{start}-$intron->{end}) : ($intron->{start}-$link->{end})), "\n";
          }
        }
      }
    }
    my $final_introns;
    my $current_introns;
    foreach my $intron (@depth_sorted_introns) {
      next if ($intron->{processed});
      going_5prime_reverse($intron, $current_introns, $final_introns);
    }
    print __LINE__, " RESULTS\n";
    my $transcript = Bio::EnsEMBL::Transcript->new;
    foreach my $transcript (@$final_introns) {
      print __LINE__, "\n";
      my @exons = (Bio::EnsEMBL::Exon->new(
        start => $transcript->[0]->{exon_5_start},
        end => $transcript->[0]->{start}-1,
        strand => -1,
        slice => $slice,
      ));
      foreach my $intron (@$transcript) {
        print __LINE__, ' ', $intron->{exon_5_start}, ' ', $intron->{start}, ' ', $intron->{end}, ' ', $intron->{exon_3_end}, ' ', $intron->{depth}, ' ', $intron->{is_canonical}, ' ', scalar(@{$intron->{extra_exon}}), "\n";
        $exons[-1]->end($intron->{start}-1);
        push(@exons, Bio::EnsEMBL::Exon->new(
        start => $intron->{end}+1,
        end => $intron->{exon_3_end},
        strand => -1,
        slice => $slice,
        ));
      }
    }
  }
}


sub write_output {
  my ($self) = @_;
}


sub going_5prime_reverse {
  my ($intron, $current_introns, $final_introns) = @_;

  if (!$intron->{processed}) {
    foreach my $key ('left_canonical', 'left_non_canonical', 'right_canonical', 'right_non_canonical') {
      if (exists $intron->{$key}) {
        @{$intron->{$key}} = sort {$b->{depth} <=> $a->{depth}} @{$intron->{$key}};
      }
    }
  }
  print __LINE__, ' ', 'PROCESSING ', $intron->{start}, ' ', $intron->{end}, ' ', $intron->{depth}, "\n";
  push (@$current_introns, $intron);
  if (exists $intron->{right_canonical}) {
    foreach my $link (@{$intron->{right_canonical}}) {
      going_5prime_reverse($link, $current_introns, $final_introns);
    }
  }
  else {
    print __LINE__, ' ', '  ADDING T NO RIGHT ', $intron->{start}, ' ', $intron->{end}, ' ', $intron->{depth}, "\n";
    add_model_to_final($current_introns, $final_introns);
  }
  pop(@$current_introns) if (@$current_introns);
  $intron->{processed} = 1;
}


sub add_model_to_final {
  my ($current_introns, $final_introns) = @_;

  if ($final_introns) {
    push(@$final_introns, []);
  }
  else {
    $final_introns = [];
  }
  foreach my $intron (@$current_introns) {
    push(@{$final_introns->[-1]}, $intron);
  }
}


sub load_introns {
  my ($self, $slice) = @_;

  my %intron_clusters = (-1 => [], 1 => []);
  my %read_groups;
  my $do_tissue = 0;
  my $total_read_count = 0;
  my $read_count = 0;
  my $dna_sequence = $slice->seq;
  my $slice_start = $slice->seq_region_start;
  my $region_to_fetch = $self->param('use_ucsc_naming') ? convert_to_ucsc_name($slice->seq_region_name, $slice) : $slice->seq_region_name;
  $region_to_fetch .= ':'.$slice_start.'-'.$slice->seq_region_end;
  foreach my $intron_files (@{$self->param('intron_bam_files')}) {
    my $bam = Bio::DB::HTS->new(
      -bam => $intron_files->{file},
      -expand_flags => 1,
    );
    $self->throw("Bam file " . $intron_files->{file} . "  not found \n") unless ($bam);
    if ($intron_files->{groupname} and scalar(@{$intron_files->{groupname}})) {
      print "Limiting to read groups ";
      foreach my $group (@{$intron_files->{groupname}}) {
        print " $group";
        $read_groups{$group} = 1;
      }
      print "\n";
      $do_tissue = 1;
    }
    my $_process_reads = sub {
      my $read = shift;
      ++$total_read_count;
      return if ($read->unmapped);
      return if ($do_tissue and !exists $read_groups{$read->get_tag_values('RG')});
      ++$read_count;
      my $start  = $read->start;
      my $end    = $read->end;
      my $strand = $read->strand;
      my $cigar_str = $read->cigar_str;
      print "$total_read_count $read_count $start $end $strand $cigar_str\n";
      if ($cigar_str =~ /N/) {
        my @introns;
        my $exon_start = $start;
        my $exon_end = $start;
        my $intron_end = $start;
        while ($cigar_str =~ /(\d*)([[:upper:]])/gc) {
          if ($2 eq 'M') {
            $exon_end += $1;
          }
          elsif ($2 eq 'N') {
            if ($1 < 11) {
              $exon_end += $1;
            }
            else {
              $intron_end = $exon_end+$1;
              if (@introns) {
                $introns[-1]->[3] = $exon_end;
                $introns[-1]->[4] = 1;
              }
              push(@introns, [$exon_end, $intron_end-1, $exon_start, $intron_end+2]);
              $exon_start = $intron_end+2;
            }
          }
          elsif ($2 eq 'D') {
            $exon_end += $1;
          }
#          I need to check that our alignment are correct to the specifications
          elsif ($2 eq 'S') {
            $exon_start += $1;
          }
        }
        $introns[-1]->[3] = $exon_end-1;
        INTRON: foreach my $intron (@introns) {
          foreach my $cluster (reverse @{$intron_clusters{$strand}}) {
            if ($cluster->{start} == $intron->[0] and $cluster->{end} == $intron->[1]) {
              ++$cluster->{depth};
              if ($cluster->{exon_5_start} > $intron->[2]) {
                $cluster->{exon_5_start} = $intron->[2];
              }
              if ($cluster->{exon_3_end} < $intron->[3]) {
                $cluster->{exon_3_end} = $intron->[3];
              }
              if (@introns > 1 and @$intron == 5) {
                push(@{$cluster->{extra_exon}}, [$intron->[2], $intron->[0]-1]);
              }
              next INTRON;
            }
            elsif ($cluster->{start} > $intron->[1]) {
              last;
            }
          }
          my %intron_cluster = (
            start => $intron->[0],
            end => $intron->[1],
            exon_5_start => $intron->[2],
            exon_3_end => $intron->[3],
            depth => 1,
            is_canonical => 0,
            extra_exon => [],
            processed => 0,
          );
          my $left_splice = substr($dna_sequence, $intron_cluster{start}-$slice_start, 2);
          my $right_splice = substr($dna_sequence, $intron_cluster{end}-$slice_start-1, 2);
          if (($strand == 1 and $left_splice eq 'GT' and $right_splice eq 'AG') or
              ($strand == -1 and $left_splice eq 'CT' and $right_splice eq 'AC')) {
            $intron_cluster{is_canonical} = 1;
          }
          elsif ( $left_splice eq 'NN' && $right_splice eq 'NN' ) {
            warn("Cannot find dna sequence for detecting non cannonical splices\n");
          }
          push(@{$intron_clusters{$strand}}, \%intron_cluster);
        }
      }
    };
    $bam->fetch($region_to_fetch, $_process_reads);
  }
#  my $dna_sequence = $slice->seq;
#  my $slice_start = $slice->seq_region_start;
#  foreach my $strand (-1, 1) {
#    foreach my $intron (@{$intron_clusters{$strand}}) {
#      my $left_splice = substr($dna_sequence, $intron->{start}-$slice_start, 2);
#      my $right_splice = substr($dna_sequence, $intron->{end}-$slice_start-1, 2);
#      print "$left_splice $right_splice\n";
#      if (($strand == 1 and $left_splice eq 'GT' and $right_splice eq 'AG') or
#          ($strand == -1 and $left_splice eq 'CT' and $right_splice eq 'AC')) {
#        $intron->{is_canonical} = 1
#      }
#      elsif ( $left_splice eq 'NN' && $right_splice eq 'NN' ) {
#        warn("Cannot find dna sequence for detecting non cannonical splices\n");
#      }
#      print __LINE__, ' ', $strand, ' ', $intron->{exon_5_start}, ' ', $intron->{start}, ' ', $intron->{end}, ' ', $intron->{exon_3_end}, ' ', $intron->{depth}, ' ', $intron->{is_canonical}, ' ', scalar(@{$intron->{extra_exon}}), "\n";
#    }
#  }
  my $last_end = $slice->seq_region_start;
  my $last_start = $slice->seq_region_end;
  my @forward_intron_clusters = ([]);
  my @reverse_intron_clusters = ([]);
  my $max_exon_size = $self->param('max_exon_size');
  my $intron_max_distance = $self->param('intron_max_distance');
  foreach my $intron (sort {$a->{end} <=> $b->{end} || $a->{start} <=> $b->{start}} reverse @{$intron_clusters{1}}) {
    if ($last_start-$intron->{end} < $intron_max_distance) {
      push(@{$forward_intron_clusters[-1]}, $intron);
    }
    else {
      push(@forward_intron_clusters, [$intron]);
    }
    $last_start = $intron->{start};
  }
  foreach my $cluster (@forward_intron_clusters) {
    if (@$cluster) {
      my $cluster_size = @$cluster;
      for (my $i = 0; $i < $cluster_size-1; $i++) {
        my $right_intron = $cluster->[$i];
        my $closest_start = $cluster->[-1]->{start};
        for (my $j = $i+1; $j < $cluster_size; $j++) {
          my $left_intron = $cluster->[$j];
          if ($closest_start < $left_intron->{start} and $right_intron->{start} > $left_intron->{end}) {
            if ($left_intron->{is_canonical}) {
              if (exists $right_intron->{left_canonical}) {
                push(@{$right_intron->{left_canonical}}, $left_intron);
              }
              else {
                $right_intron->{left_canonical} = [$left_intron];
              }
              $closest_start = $left_intron->{start} if ($closest_start < $left_intron->{start});
            }
            else {
              if (exists $right_intron->{left_non_canonical}) {
                push(@{$right_intron->{left_non_canonical}}, $left_intron);
              }
              else {
                $right_intron->{left_non_canonical} = [$left_intron];
              }
            }
            if ($right_intron->{is_canonical}) {
              if (exists $left_intron->{right_canonical}) {
                push(@{$left_intron->{right_canonical}}, $right_intron);
              }
              else {
                $left_intron->{right_canonical} = [$right_intron];
              }
            }
            else {
              if (exists $left_intron->{right_non_canonical}) {
                push(@{$left_intron->{right_non_canonical}}, $right_intron);
              }
              else {
                $left_intron->{right_non_canonical} = [$right_intron];
              }
            }
          }
          elsif ($closest_start < $left_intron->{end}) {
            last;
          }
        }
      }
    }
  }
  foreach my $intron (sort {$a->{start} <=> $b->{start} || $b->{end} <=> $a->{end}} @{$intron_clusters{-1}}) {
    if ($intron->{start}-$last_end < $intron_max_distance) {
      push(@{$reverse_intron_clusters[-1]}, $intron);
#      print __LINE__, ' ', $intron->{exon_5_start}, ' ', $intron->{start}, ' ', $intron->{end}, ' ', $intron->{exon_3_end}, ' ', $intron->{depth}, ' ', $intron->{is_canonical}, ' ', scalar(@{$intron->{extra_exon}}), ' ', $last_end, ' ', scalar(@reverse_intron_clusters), "\n";
    }
    else {
      push(@reverse_intron_clusters, [$intron]);
#      print __LINE__, ' ', $intron->{exon_5_start}, ' ', $intron->{start}, ' ', $intron->{end}, ' ', $intron->{exon_3_end}, ' ', $intron->{depth}, ' ', $intron->{is_canonical}, ' ', scalar(@{$intron->{extra_exon}}), ' ', $last_end, ' ', scalar(@reverse_intron_clusters), "\n";
    }
    $last_end = $intron->{end};
  }
  foreach my $cluster (@reverse_intron_clusters) {
    if (@$cluster) {
      print __LINE__, ' ', 'CLUSTER', "\n";
      for (my $i = 0; $i < @$cluster-1; $i++) {
        my $left_intron = $cluster->[$i];
        my $closest_end = $cluster->[-1]->{end};
        for (my $j = $i+1; $j < @$cluster; $j++) {
          my $right_intron = $cluster->[$j];
          if ($closest_end > $right_intron->{start} and $left_intron->{end} < $right_intron->{start} and $right_intron->{start}-$left_intron->{end} < $max_exon_size) {
            if ($right_intron->{is_canonical}) {
              if (exists $left_intron->{right_canonical}) {
                push(@{$left_intron->{right_canonical}}, $right_intron);
              }
              else {
                $left_intron->{right_canonical} = [$right_intron];
              }
              $closest_end = $right_intron->{end} if ($closest_end > $right_intron->{end});
            }
            else {
              if (exists $left_intron->{right_non_canonical}) {
                push(@{$left_intron->{right_non_canonical}}, $right_intron);
              }
              else {
                $left_intron->{right_non_canonical} = [$right_intron];
              }
            }
            if ($left_intron->{is_canonical}) {
              if (exists $right_intron->{left_canonical}) {
                push(@{$right_intron->{left_canonical}}, $left_intron);
              }
              else {
                $right_intron->{left_canonical} = [$left_intron];
              }
            }
            else {
              if (exists $right_intron->{left_non_canonical}) {
                push(@{$right_intron->{left_non_canonical}}, $left_intron);
              }
              else {
                $right_intron->{left_non_canonical} = [$left_intron];
              }
            }
          }
          elsif ($closest_end < $right_intron->{start}) {
            last;
          }
        }
      }
    }
  }
#  foreach my $strand (-1, 1) {
#    foreach my $intron (@{$intron_clusters{$strand}}) {
#      print __LINE__, ' ', $strand, ' ', $intron->{exon_5_start}, ' ', $intron->{start}, ' ', $intron->{end}, ' ', $intron->{exon_3_end}, ' ', $intron->{depth}, ' ', $intron->{is_canonical}, ' ', scalar(@{$intron->{extra_exon}}), "\n";
#      foreach my $key ('left_canonical', 'left_non_canonical', 'right_canonical', 'right_non_canonical') {
#        if (exists $intron->{$key}) {
#          foreach my $link (@{$intron->{$key}}) {
#            print __LINE__, '   ', $key, ' ', $strand, ' ', $link->{exon_5_start}, ' ', $link->{start}, ' ', $link->{end}, ' ', $link->{exon_3_end}, ' ', $link->{depth}, ' ', $link->{is_canonical}, ' ', scalar(@{$link->{extra_exon}}), "\n";
#          }
#        }
#      }
#    }
#  }
  $self->param('intron_clusters', \%intron_clusters);
}


sub chr_slice {
  my ($self, $slice) = @_;

  if ($slice) {
    $self->throw("$slice is not a Bio::EnsEMBL::Slice") unless $slice->isa('Bio::EnsEMBL::Slice');
    $self->param('seq_region_slice', $slice);
  }
  if ($self->param_is_defined('seq_region_slice')) {
    return $self->param('seq_region_slice');
  }
  else {
    return;
  }
}

1;
