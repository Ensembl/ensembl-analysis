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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRoughTranscipts

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRoughTranscripts;

use strict;
use warnings;

use Bio::DB::HTS;
use Bio::DB::HTS::Constants qw(FLAGS);

use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Exon;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    stranded_read => 0,
    _repeat_libs => ['dust', 'trf'],
    _repeats_allowed => 0.95,
    use_ucsc_naming => 0,
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: It will fetch the sequence based on the input_id, retrieve all alignments
               overlapping the region and create the "exon" blocks. Because the alignment
               is made on the whole genome, we need to provide the runnable with a full length
               slice. If your genome use UCSC style names (chr1,...), set 'use_ucsc_naming' to 1
 Returntype : None
 Exceptions : None

=cut

sub fetch_input {
  my ($self) = @_;

  my $reference_db = $self->get_database_by_name('dna_db');
  my $slice_adaptor = $reference_db->get_SliceAdaptor;

  my $slice = $self->fetch_sequence($self->input_id, $reference_db);
  if ($self->param('disconnect_jobs')) {
    $reference_db->dbc->disconnect_if_idle;
    $self->dbc->disconnect_if_idle;
  }
  if ($slice) {
    my $bam = Bio::DB::HTS->new(-bam => $self->param_required('alignment_bam_file'));
    my $bam_file = $bam->hts_file;
    my $bam_index = $bam->hts_index;
    if ($bam_file and $bam_index) {
      $self->create_analysis;
      my $target_db = $self->get_database_by_name('output_db', $reference_db);
      $self->hrdb_set_con($target_db, 'target_db');
      #Coordinates are based on the full sequence so we need to make sure we have the full sequence in query
      $self->query($target_db->get_SliceAdaptor->fetch_by_region(
        $slice->coord_system->name,
        $slice->seq_region_name,
        undef,
        undef,
        undef,
        $slice->coord_system->version,
      ));
      my $seq_region_name = $self->param('use_ucsc_naming') ? convert_to_ucsc_name($slice->seq_region_name, $slice) : $slice->seq_region_name;
      my ($tid, $zbased_start, $zbased_end) = $bam->header->parse_region($seq_region_name.':'.$slice->start.'-'.$slice->end);
      $self->param('region_to_fetch', [$tid, $zbased_start, $zbased_end]);
      $self->param('bam_file', $bam_file);
      $self->param('bam_index', $bam_index);
    }
    else {
      $self->throw('Could not open bam to access the alignments');
    }
  }
  else {
    $self->throw('Could not find sequence '.$self->param('iid'));
  }
}


sub run {
  my ($self) = @_;

  my $slice = $self->query;
  my $last_pairing = 0;
  my $analysis = $self->analysis;
  my $callback_data = {
    -1 => [],
    1 => [],
    cluster_count => 1,
    read_count => 0,
    total_read_count => 0,
    stranded_reads => $self->param('stranded_read'),
    last_pairing => \$last_pairing,
    analysis => $analysis,
    slice => $slice,
  };
  my $bam_file = $self->param('bam_file');
  my $bam_index = $self->param('bam_index');
  $self->say_with_header('Starting the loop');
  my $_process_reads = sub {
    my ($read, $data) = @_;
    ++$data->{total_read_count};
    # It seems we always get the unmmapped mate, so we need to remove it
    return if ($read->unmapped);
    ++$data->{read_count};
    my $query = $read->query;
    my $start  = $read->start;
    my $end    = $read->end;
    my $paired = $read->proper_pair;
    my $real_strand = -1; # We are making sure that we are -1 if not stranded, which should probably be changed...
    if ($data->{stranded_reads}) {
        $real_strand = $read->strand;
        $real_strand *= -1 if (($read->flag & FLAGS->{FIRST_MATE}) != 0);
    }
    my $exon_clusters = $data->{$real_strand};
    for (my $index = @$exon_clusters-1; $index > -1; $index--) {
      my $exon_cluster = $exon_clusters->[$index];
      if ($exon_cluster->start <= $end and $exon_cluster->end >= $start) {
        $exon_cluster->end($end) if ($exon_cluster->end < $end);
        $exon_cluster->score($exon_cluster->score+1);
        if ($paired and $read->isize < 0) {
          process_mate_pairs($read, $exon_cluster, $exon_clusters, $index, $data->{last_pairing});
        }
        return;
      }
      elsif ($exon_cluster->end < $start) {
        last;
      }
    }
    my $feat = Bio::EnsEMBL::DnaDnaAlignFeature->new(
      -start      => $start,
      -end        => $end,
      -strand     => $real_strand,
      -slice      => $data->{slice},
      -hstart     => $query->start,
      -hend       => $query->end,
      -hstrand    => 1,
      -hseqname   => $data->{cluster_count}++,
      -analysis   => $data->{analysis},
      -cigar_string => '1M', # IT will be set properly later
      -align_type => 'ensembl',
    );
    $feat->score(1);
    $feat->percent_id(100);
    $feat->hcoverage(0);
    push(@$exon_clusters, $feat);
    if ($paired and $read->isize < 0) {
      process_mate_pairs($read, $feat, $exon_clusters, scalar(@$exon_clusters)-1, $data->{last_pairing});
    }
  };
  $self->say_with_header(join(' ', @{$self->param('region_to_fetch')}));
  $bam_index->fetch($self->param('bam_file'), @{$self->param('region_to_fetch')}, $_process_reads, $callback_data);
  $self->say_with_header($callback_data->{total_read_count}.' reads in the region, '.$callback_data->{read_count}.' reads processed');
  my @genes;
  my $max_intron_length = $self->param('max_intron_length');
  if ($self->param('paired')) {
    foreach my $stranded_exons ($callback_data->{-1}, $callback_data->{1}) {
      if (@$stranded_exons) {
        my $clusters = [];
        my $cluster_index = -1;
        my %objects;
        foreach my $object (@$stranded_exons) {
          $objects{$object->hseqname} = 1;
          if (exists $object->{links}) {
            foreach my $link (@{$object->{links}}) {
            }
          }
          if ($object->hcoverage == 0) {
            if (exists $object->{links}) {
              $cluster_index = @$clusters;
              my $tmp_index = process_link($object, $clusters, $cluster_index, '');
              if ($tmp_index != $cluster_index) {
              }
            }
            else {
              $object->hcoverage(1);
              push(@$clusters, [$object]);
              $cluster_index = @$clusters;
            }
          }
        }
        foreach my $cluster (@$clusters) {
          if(@$cluster) {
            my $transcript = Bio::EnsEMBL::Transcript->new();
            foreach my $object (@$cluster) {
              my $exon = Bio::EnsEMBL::Exon->new(
                -slice => $object->slice,
                -start => $object->start,
                -end => $object->end,
                -strand => $object->strand,
                -analysis => $analysis,
              );
              $exon->phase(-1);
              $exon->end_phase(-1);
              $object->cigar_string($exon->length.'M');
              $exon->add_supporting_features($object);
              $transcript->add_Exon($exon);
            }
            if ($self->is_transcript_valid($transcript)) {
              $transcript->analysis($analysis);
              $transcript->biotype($analysis->logic_name);
              my $gene = Bio::EnsEMBL::Gene->new();
              $gene->add_Transcript($transcript);
              $gene->analysis($analysis);
              $gene->biotype($analysis->logic_name);
              push(@genes, $gene);
            }
            else {
              $self->say_with_header('Transcript without exons or invalid');
            }
          }
        }
      }
    }
  }
  else {
    foreach my $stranded_exons ($callback_data->{-1}, $callback_data->{1}) {
      if (@$stranded_exons) {
        my @transcripts = (Bio::EnsEMBL::Transcript->new());
        foreach my $object (@$stranded_exons) {
          my $exon = Bio::EnsEMBL::Exon->new(
            -slice => $object->slice,
            -start => $object->start,
            -end => $object->end,
            -strand => $object->strand,
            -analysis => $analysis,
          );
          $exon->phase(-1);
          $exon->end_phase(-1);
          $exon->add_supporting_features($object);
          if ($transcripts[-1]->end and $exon->start-$transcripts[-1]->end >= $max_intron_length) {
            if ($self->is_transcript_valid($transcripts[-1])) {
              push(@transcripts, Bio::EnsEMBL::Transcript->new());
            }
            else {
              $transcripts[-1]->flush_Exons;
            }
          }
          $transcripts[-1]->add_Exon($exon);
        }
        foreach my $transcript (@transcripts) {
          if ($transcript->start and $self->is_transcript_valid($transcript)) {
            my $gene = Bio::EnsEMBL::Gene->new();
            $gene->add_Transcript($transcript);
            $gene->analysis($analysis);
            push(@genes, $gene);
          }
          else {
            $self->say_with_header('Transcript without exons or invalid');
          }
        }
      }
    }
  }
  if (@genes) {
    $self->say_with_header(scalar(@genes).' rough transcripts');
    $self->output(\@genes);
  }
  else {
    $self->input_job->autoflow(0);
    $self->complete_early('No alignment found');
  }
}


sub process_mate_pairs {
  my ($read, $exon_cluster, $exon_clusters, $index, $last_pairing) = @_;

  my $mate_start = $read->mate_start;
  my $mate_end = $read->mate_end;
#  my $name = $read->name;
  if (!($mate_end >= $exon_cluster->start and $mate_start <= $exon_cluster->end)) {
    my $j = $index-1;
    if ($exon_clusters->[$$last_pairing]->start <= $mate_end and $exon_clusters->[$$last_pairing]->end >= $mate_start) {
      my $mate_cluster = $exon_clusters->[$$last_pairing];
      if (exists $mate_cluster->{links}) {
        foreach my $link (reverse @{$mate_cluster->{links}}) {
          if ($link == $exon_cluster) {
#            print STDERR ' LINK LAST ', $name, ' ', $exon_cluster->start, ' ', $exon_cluster->end, ' ', $exon_cluster->strand, ' ', $mate_start, ' ', $mate_end, ' ',$read->start, ' ', $read->end, ' ', $index, ' ', $j, ' ', $mate_cluster->hseqname, ' ', $exon_cluster->hseqname, "\n";
            return;
          }
        }
      }
      push(@{$mate_cluster->{links}}, $exon_cluster);
#      print STDERR ' NEW LINK LAST ', $name, ' ', $exon_cluster->start, ' ', $exon_cluster->end, ' ', $exon_cluster->strand, ' ', $mate_start, ' ', $mate_end, ' ',$read->start, ' ', $read->end, ' ', $index, ' ', $j, ' ', $mate_cluster->hseqname, ' ', $exon_cluster->hseqname, "\n";
      return;
    }
    elsif ($exon_clusters->[$j]->start <= $mate_end and $exon_clusters->[$j]->end >= $mate_start) {
      my $mate_cluster = $exon_clusters->[$j];
      $$last_pairing = $j;
      if (exists $mate_cluster->{links}) {
        foreach my $link (reverse @{$mate_cluster->{links}}) {
          if ($link == $exon_cluster) {
#            print STDERR ' LINK PREVIOUS ', $name, ' ', $exon_cluster->start, ' ', $exon_cluster->end, ' ', $exon_cluster->strand, ' ', $mate_start, ' ', $mate_end, ' ',$read->start, ' ', $read->end, ' ', $index, ' ', $j, ' ', $mate_cluster->hseqname, ' ', $exon_cluster->hseqname, "\n";
            return;
          }
        }
      }
      push(@{$mate_cluster->{links}}, $exon_cluster);
#      print STDERR ' NEW LINK PREVIOUS ', $name, ' ', $exon_cluster->start, ' ', $exon_cluster->end, ' ', $exon_cluster->strand, ' ', $mate_start, ' ', $mate_end, ' ',$read->start, ' ', $read->end, ' ', $index, ' ', $j, ' ', $mate_cluster->hseqname, ' ', $exon_cluster->hseqname, "\n";
      return;
    }
    else {
      my $low_boundary = 0;
      my $high_boundary = $j;
      my $start = $read->start;
      my $isize = $read->isize;
      $j = int(($mate_start/$exon_clusters->[$j]->start)*$high_boundary);
      print STDERR "START $j $start $low_boundary $high_boundary $mate_start $mate_end ", $isize, "\n";
      while (!($exon_clusters->[$j]->start <= $mate_end and $exon_clusters->[$j]->end >= $mate_start)) {
        print STDERR "  LOOP $j $start $low_boundary $high_boundary $mate_start $mate_end ", $isize, ' ', $exon_clusters->[$j]->start, ' ', $exon_clusters->[$j]->end, "\n";
        if ($exon_clusters->[$j]->start > $mate_end) {
          $high_boundary = $j;
          $j = $high_boundary-(int(($mate_start/$exon_clusters->[$j]->start)*($high_boundary-$low_boundary)) || 1);
        }
        else {
          $low_boundary = $j;
          $j = $low_boundary+(int(($mate_start/$exon_clusters->[$high_boundary]->start)*($high_boundary-$low_boundary)) || 1);
        }
        print STDERR "  NEW J $j\n";
      }
      my $mate_cluster = $exon_clusters->[$j];
      $$last_pairing = $j;
      if (exists $mate_cluster->{links}) {
        foreach my $link (reverse @{$mate_cluster->{links}}) {
          if ($link == $exon_cluster) {
#            print STDERR ' LINK ', $name, ' ', $exon_cluster->start, ' ', $exon_cluster->end, ' ', $exon_cluster->strand, ' ', $mate_start, ' ', $mate_end, ' ',$read->start, ' ', $read->end, ' ', $index, ' ', $j, ' ', $mate_cluster->hseqname, ' ', $exon_cluster->hseqname, "\n";
            return;
          }
        }
      }
      push(@{$mate_cluster->{links}}, $exon_cluster);
#      print STDERR ' NEW LINK ', $name, ' ', $exon_cluster->start, ' ', $exon_cluster->end, ' ', $exon_cluster->strand, ' ', $mate_start, ' ', $mate_end, ' ',$read->start, ' ', $read->end, ' ', $index, ' ', $j, ' ', $mate_cluster->hseqname, ' ', $exon_cluster->hseqname, "\n";
      return;
    }
#    print STDERR 'NO LINK ', $name, ' ', $exon_cluster->start, ' ', $exon_cluster->end, ' ', $exon_cluster->strand, ' ', $mate_start, ' ', $mate_end, ' ',$read->start, ' ', $read->end, ' ', $index, ' ', $exon_cluster->hseqname, ' ', $exon_cluster->hseqname, "\n";
  }
#  print STDERR 'OVERLAP ', $name, ' ', $exon_cluster->start, ' ', $exon_cluster->end, ' ', $exon_cluster->strand, ' ', $mate_start, ' ', $mate_end, ' ',$read->start, ' ', $read->end, ' ', $index, ' ', $exon_cluster->hseqname, ' ', $exon_cluster->hseqname, "\n";
}


sub process_link {
  my ($object, $cluster, $index, $filler) = @_;

  $object->hcoverage($object->hcoverage+1);
#  print STDERR $filler.'Processing ', $object->hseqname, ' ', $object->start, ' ', $object->end, ' ', $object->hcoverage, "\n";
  if (exists $object->{index}) {
    if ($index != $object->{index}) {
      my $old_index = $object->{index};
      while(my $elm = pop(@{$cluster->[$old_index]})) {
        if (grep {$elm == $_} @{$cluster->[$index]}) {
#          print STDERR $filler.'Already there ', $elm->hseqname, ' ', $elm->start, ' ', $elm->end, "\n";
        }
        else {
#          print STDERR $filler.'Moving ', $elm->hseqname, ' ', $elm->start, ' ', $elm->end, ' from ', $elm->{index}, " to $index\n";
          push(@{$cluster->[$index]}, $elm);
          $elm->{index} = $index;
        }
      }
    }
    else {
#      print STDERR $filler.'Already processed ', $object->hseqname, ' ', $object->start, ' ', $object->end, "\n";
    }
  }
  else {
    $object->{index} = $index;
    if (exists $object->{links}) {
      foreach my $link (@{$object->{links}}) {
        my $tmpindex = process_link($link, $cluster, $index, $filler.' ');
        if ($index != $tmpindex) {
#          print STDERR $filler.'Changed index for ', $object->hseqname, ' ', $object->start, ' ', $object->end, ' ', $index, ' to ', $tmpindex, ' because of ', $link->hseqname, ' ', $link->start, ' ', $link->end, "\n";
          $index = $tmpindex;
        }
      }
    }
#  push(@{$cluster->[$index]}, $object) unless (grep {$object == $_} @{$cluster->[$index]});
    if (!(grep {$object == $_} @{$cluster->[$index]})) {
#      print STDERR $filler.'Adding ', $object->hseqname, ' ', $object->start, ' ', $object->end, " TO $index\n";
      push(@{$cluster->[$index]}, $object);
    }
    else {
#      print STDERR $filler.'DOUBLE ', $object->hseqname, ' ', $object->start, ' ', $object->end, " AT $index\n";
    }
  }
  return $index;
}


sub is_transcript_valid {
  my ($self, $transcript) = @_;

  if ($transcript->length < $self->param('min_length')) {
    $self->say_with_header('Rejected on length: '.$transcript->length.' < '.$self->param('min_length'));
    return 0;
  }
  if (scalar(@{$transcript->get_all_Exons}) == 1) {
    if ($transcript->length < $self->param('min_single_exon_length')) {
      $self->say_with_header('Rejecting single exon transcript because of length '.$transcript->length);
      return 0;
    }
  }
  else {
# filter span on multiexon genes
    if (($transcript->end-$transcript->start+1)/$transcript->length < $self->param('min_span')) {
      if ($transcript->length < $self->param('min_single_exon_length')) {
        $self->say_with_header('Rejecting because of span '.(($transcript->end-$transcript->start+1)/$transcript->length));
        return 0;
      }
    }
  }
  return 1;
}


sub filter_results {
  my ($self, $results) = @_;

  $self->throw('You did not pass an array ref but a '.ref($results)) unless (ref($results) eq 'ARRAY');
  my @filtered_results;
  my $repeat_libs = $self->param('_repeat_libs');
  my $repeats_allowed = $self->param('_repeats_allowed');
  foreach my $gene (@$results) {
    if (@{$gene->get_all_Exons} == 1) {
      my $masked_slice = $gene->feature_Slice->get_repeatmasked_seq($repeat_libs);
      my $N_count = $masked_slice->seq =~ tr/N/N/;
      if ($N_count > $masked_slice->length*$repeats_allowed) {
        $self->warning('Too many simple (TRF/dust) repeats, not storing single exon gene');
      }
      else {
        push(@filtered_results, $gene);
      }
    }
    else {
      push(@filtered_results, $gene);
    }
  }
  return \@filtered_results;
}

=head2 write_output

 Arg [1]    : None
 Description: Write the proto transcripts in the database specified as 'output_db'.
 Returntype : None
 Exceptions : Throws if it cannot write all the genes in the database for any reason

=cut

sub write_output{
  my ($self) = @_;

  my $outdb = $self->hrdb_get_con('target_db');
  my $gene_adaptor = $outdb->get_GeneAdaptor;

  $gene_adaptor->dbc->disconnect_when_inactive(0);
  foreach my $gene ( @{$self->output} ) {
    $gene_adaptor->store($gene);
  }
}


1;
