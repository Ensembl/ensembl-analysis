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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBamAlign2Gene

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBamAlign2Gene;

use strict;
use warnings;

use Bio::DB::HTS;

use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::EvidenceUtils qw(create_feature);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my ($self) = @_;

  my $target_db = $self->get_database_by_name('target_db', $self->get_database_by_name('dna_db'));
  my $slice = $target_db->get_SliceAdaptor->fetch_by_name($self->param_required('iid'));
  if ($slice) {
    my $bam = Bio::DB::HTS->new(
      -bam => $self->param_required('bam_file'),
      -expand_flags => 1,
    );
    my $iterator = $bam->features(
      -iterator => 1,
      -seq_id => $slice->seq_region_name,
      -start => $slice->start,
      -end => $slice->end,
    );
    if ($iterator) {
      $self->create_analysis;
      $self->param('bam_iterator', $iterator);
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
    }
    else {
      $self->throw('Could not create iterator to access the alignments');
    }
  }
  else {
    $self->throw('Could not find sequence '.$self->param('iid'));
  }
}


sub run {
  my ($self) = @_;

  my $iterator = $self->param('bam_iterator');
  my @genes;
  my $slice = $self->query;
  while (my $seq = $iterator->next_seq) {
    my $query = $seq->query;
    my $isoseq_length = length($query->dna);
    print STDERR $seq->start, ' ', $seq->end, ' ', $seq->strand, ' ', $query->name, ' ', $query->start, ' ', $query->end, ' ', $query->length, ' ', $isoseq_length, ' ', $seq->cigar_str, ' ', $seq->aux, "\n";
    foreach my $tag ($seq->get_all_tags) {
      print STDERR $tag, ' ', $seq->get_tag_values($tag), "\n";
    }
    my ($score) = $seq->get_tag_values('AS');
    my ($h_strand) = $seq->get_tag_values('ts');
    $h_strand = $h_strand eq '+' ? 1 : -1;
    my $cigar_line = $seq->cigar_str;
    my @exons;
    my @sfs;
    my $seqname = $query->name;
    my $strand = $seq->strand;
    my $slice_position = $slice->start;
    my $exon = Bio::EnsEMBL::Exon->new(
      -start => $slice_position,
      -end => $seq->end,
      -strand => $strand,
    );
    push(@exons, $exon);
    my $exon_cigar;
    my $sf = [
      $query->start,
      $query->end,
      '',
    ];
    push(@sfs, $sf);
    my $seq_length = 0;
    while ($cigar_line =~ /(\d+)(\D)/gc)  {
      my ($length, $operation) = ($1, $2);
      if ($operation eq 'N') {
        $slice_position += $length;
        $exon = Bio::EnsEMBL::Exon->new(
          -start => $slice_position,
          -end => $slice_position+1,
          -strand => $strand,
        );
        push(@exons, $exon);
        $sf = [
          $seq_length,
          $seq_length+1,
          '',
        ];
        push(@sfs, $sf);
      }
      elsif ($operation eq 'M') {
        if ($strand == 1) {
          $sfs[-1]->[2] .= $length.$operation;
        }
        else {
          $sfs[-1]->[2] = $length.$operation.$sfs[-1]->[2];
        }
        $slice_position += $length;
        $seq_length += $length;
        $exon->end($slice_position-1);
        $sfs[-1]->[1] += $seq_length-1;
      }
      elsif ($operation eq 'D') {
        if ($strand == 1) {
          $sfs[-1]->[2] .= $length.$operation;
        }
        else {
          $sfs[-1]->[2] = $length.$operation.$sfs[-1]->[2];
        }
        $slice_position += $length;
        $exon->end($slice_position-1);
      }
      elsif ($operation eq 'I') {
        if ($strand == 1) {
          $sfs[-1]->[2] .= $length.$operation;
        }
        else {
          $sfs[-1]->[2] = $length.$operation.$sfs[-1]->[2];
        }
        $seq_length += $length;
        $sfs[-1]->[1] += $seq_length-1;
      }
      elsif ($operation eq 'S') {
        $seq_length += $length;
      }
      elsif ($operation eq 'H') {
        $seq_length += $length;
      }
      else {
        $self->throw("Not expected operation $length$operation for $seqname on ".$slice->seq_region_name);
      }
    }
#  my ($feature_string, $start, $end, $strand, $hstart, $hend, 
#      $hstrand, $percent_id, $p_value, $score, $hseqname,
#      $cigar_string, $analysis, $external_db_id, $hcoverage, $slice) = @_;
    print STDERR "$seq_length $isoseq_length\n";
    my $index = 0;
    my $transcript = Bio::EnsEMBL::Transcript->new();
    foreach my $exon (@exons) {
      print STDERR $exon->start, ' ', $exon->end, ' ', $exon->strand, "\n";
      my $sf = $sfs[$index++];
      my $hcoverage = (($sf->[1]-$sf->[0]+1)*100)/$seq_length;
      my $start = $sf->[0];
      my $end = $sf->[1];
      if ($exon->strand == -1) {
        $start = $seq_length-$sf->[0]+1;
        $start = $seq_length-$sf->[1]+1;
      }
      my $evidence = create_feature(
        'Bio::EnsEMBL::DnaDnaAlignFeature',
        $exon->start, $exon->end, $exon->strand,
        $start, $end, $h_strand,
        100, 0, $score,
        $seqname, $sf->[2], $self->analysis, undef, $hcoverage, $slice,
      );
      $exon->add_supporting_features($evidence);
      print STDERR '  ', $evidence->hseqname, ' ', $evidence->hstart, ' ', $evidence->hend, ' ', $evidence->hstrand, ' ', $evidence->cigar_string, "\n";
      $transcript->add_Exon($exon);
    }
    $transcript->analysis($self->analysis);
    $transcript->stable_id($seqname);
    $transcript->slice($slice);
    my $gene = Bio::EnsEMBL::Gene->new();
    $gene->slice($slice);
    $gene->add_Transcript($transcript);
    $gene->analysis($self->analysis);
    $gene->stable_id($seqname);
    push(@genes, $gene);
  }
  if (@genes) {
    $self->output(\@genes);
  }
  else {
    $self->input_job->autoflow(0);
    $self->complete_early('No alignment found');
  }
}


sub get_adaptor {
  my ($self) = @_;

  return $self->hrdb_get_con('target_db')->get_GeneAdaptor;
}

1;
