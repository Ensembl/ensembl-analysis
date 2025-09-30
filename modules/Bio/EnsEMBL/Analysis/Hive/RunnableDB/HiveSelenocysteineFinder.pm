=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2024] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSelenocysteineFinder

=head1 SYNOPSIS


=head1 DESCRIPTION

This module will run exonerate first to find where are the best 3 alignments.
Then it will run exonerate with the exhaustive parameter to refine the alignment.
Finally it will run genewise and compare the two alignment and the original protein
to find the best match. To allow genewise to create a good model we have to modify
the dna sequence to change the TGA for selenocysteine to a TGC for a cysteine

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSelenocysteineFinder;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::IO::FASTASerializer;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript;
use Bio::EnsEMBL::Analysis::Runnable::Genewise;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(attach_Analysis_to_Gene attach_Slice_to_Gene);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(attach_Slice_to_Transcript);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(write_seqfile);

use Bio::SeqIO;

use parent qw(Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDBSeqFiles);

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults()},
    padding => 200000,
# We allow 1 missmatch or 3% of mismatches between the original sequence and the predicted sequence
    missmatch_allowed => 3,
# For some sequences the alignment with the highest score is not the best alignment
    score_threshold => 2,
    minimum_identity => 95,
    coverage_threshold => 90,
    biotype => 'seleno',
  }
}


sub fetch_input {
  my ($self) = @_;

  my $target = $self->param_required('genome');
  my $exonerate = $self->param_required('exonerate');
  $self->param_required('genewise');
  my $dnadb = $self->get_database_by_name('dna_db');
  my $db = $self->get_database_by_name('target_db', $dnadb);
  $self->hrdb_set_con($db, 'target_db');
  my $analysis = $self->create_analysis;

  my $querys;
  my $iid_type = $self->param_required('iid_type');
  if ($iid_type eq 'db_seq') {
    my $iid = $self->input_id;
    if (!ref($iid)) {
      $iid = [$iid];
    }
    $querys = $self->get_query_seqs($iid);
  }
  elsif ($iid_type eq 'filename') {
    $querys = get_fasta_file_from_uniprot($self->input_id);
  }

  my $seqs = $self->prepare_sequences($querys);
  if ($seqs) {
    foreach my $seq (@$seqs) {
      $self->runnable(Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript->new(
        -program => $exonerate,
        -analysis => $analysis,
        -query_seqs => [$seq],
        -query_type => 'protein',
        -target_file => $target,
        -options => '--model protein2genome --bestn 3 --softmasktarget true',
      ));
    }
  }
  else {
    $self->input_job->autoflow(0);
    $self->complete_early('No selenoprotein to use');
  }
 if ($self->param('disconnect_jobs')) {
     $db->dbc->disconnect_when_inactive(1);
  }
}


sub run {
  my ($self) = @_;
  $self->dbc->disconnect_when_inactive(1) if ($self->param('disconnect_jobs'));
  my $slice_adaptor = $self->hrdb_get_con('target_db')->get_SliceAdaptor;
  my $minimum_identity = $self->param('minimum_identity');
  my $coverage_threshold = $self->param('coverage_threshold')/100;
  my $score_threshold = $self->param('score_threshold')/100;
  my $padding = $self->param('padding');
  my $exonerate = $self->param('exonerate');
  foreach my $runnable (@{$self->runnable}) {
    $runnable->run;
    my %alignments;
    foreach my $transcript (@{$runnable->output}) {
      my $tsf = $transcript->get_all_supporting_features->[0];
      next if ($tsf->hcoverage < $self->param('coverage_threshold'));
      if (exists $alignments{$tsf->seqname}) {
        if ($tsf->score > ($alignments{$tsf->seqname}->{score}+$tsf->{score}*$score_threshold)
            or ($tsf->percent_id > $minimum_identity and
              $tsf->percent_id > $alignments{$tsf->seqname}->get_all_supporting_features->[0]->{percent_id} and
		(($tsf->hend-$tsf->hstart) > $transcript->{query_length}*$coverage_threshold)) ) {
	    $transcript->slice($alignments{$tsf->seqname}->slice);
	    $alignments{$tsf->seqname} = $transcript;
        }
      }
      else {
        $alignments{$tsf->seqname} = $transcript;
        $transcript->slice($slice_adaptor->fetch_by_name($tsf->seqname));
      }
    }
    my @selenocysteine_transcripts;
    foreach my $transcript ( sort {$a->{score} > $b->{score}} values %alignments) {
      my $accession = $transcript->{accession};
      my ($seq) = grep {$_->id eq $accession} @{$runnable->query_seqs};
      print 'Working on ', $accession, ' ', $transcript->seq_region_start, ' ', $transcript->seq_region_end, ' ', $transcript->seq_region_strand, ' and ', $transcript->slice->name, "\n";
      $self->create_target_file('selenocysteine', 'fa');
      my $sub_slice = $self->param('target_file');
      my $fasta_serial = Bio::EnsEMBL::Utils::IO::FASTASerializer->new($sub_slice);
      
      my $transcript_slice_start_temp = $transcript->seq_region_start-$padding;
      $transcript_slice_start_temp = 1 if ($transcript_slice_start_temp<1); 
      my $transcript_slice_end_temp = $transcript->seq_region_end+$padding;
      my $transcript_seq_region_length = $transcript->slice->seq_region_length();
      $transcript_slice_end_temp = $transcript_seq_region_length if ($transcript_slice_end_temp>$transcript_seq_region_length);  
      
      my $gene_slice = $transcript->slice->sub_Slice($transcript_slice_start_temp, $transcript_slice_end_temp);
      $fasta_serial->print_Seq($gene_slice);

      my $exonerate_runnable = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript->new(
          -program => $runnable->program,
          -analysis => $runnable->analysis,
          -query_seqs => [$seq],
          -query_type => $runnable->query_type,
          -target_file => $sub_slice->filename,
          -options => '--bestn 1 --exhaustive --model protein2genome',
          );
      $exonerate_runnable->run;
      foreach my $exonerate_transcript (@{$exonerate_runnable->output}) {
        attach_Slice_to_Transcript($exonerate_transcript, $gene_slice);
      print 'Working on ', $accession, ' ', $exonerate_transcript->seq_region_start, ' ', $exonerate_transcript->seq_region_end, ' ', $exonerate_transcript->seq_region_strand, ' and ', $exonerate_transcript->slice->name, "\n";
        my $num_seqEdits = $self->get_SeqEdit($exonerate_transcript, $seq);
        if ($num_seqEdits) {
          my $genewise_padding = 2*$self->param('padding');
# For all our transcripts we try to find the genomic position of the start so we can change TGA to TGC
# Then we can use GeneWise to predict a model.

          my $genewise_slice_start_temp = $exonerate_transcript->seq_region_start-$genewise_padding;
          $genewise_slice_start_temp = 1 if ($genewise_slice_start_temp<1); 
          my $genewise_slice_end_temp = $exonerate_transcript->seq_region_end+$genewise_padding;
          my $gene_seq_region_length = $transcript->slice->seq_region_length();
          $genewise_slice_end_temp = $gene_seq_region_length if ($genewise_slice_end_temp>$gene_seq_region_length);  
          my $genewise_slice = $transcript->slice->sub_Slice($genewise_slice_start_temp, $genewise_slice_end_temp);
          my $genewise_sequence = $self->mutate_dna($exonerate_transcript, $genewise_slice);
          if ($genewise_sequence) {
# A Slice with its sequence manually modified can't use the SequenceAdaptor.
# So we need to create a copy of the Slice. We will use the copy for GeneWise
# then use the original one to get the translation.
            my $edited_slice = Bio::EnsEMBL::Slice->new(
                -start => $genewise_slice->start,
                -end => $genewise_slice->end,
                -strand => $genewise_slice->strand,
                -seq_region_name => $genewise_slice->seq_region_name,
                -coord_system => $genewise_slice->coord_system,
                );
            $edited_slice->{seq} = $$genewise_sequence;
            my $genewise_runnable = Bio::EnsEMBL::Analysis::Runnable::Genewise->new (
                -protein => $seq,
                -query => $edited_slice,
                -reverse => $exonerate_transcript->strand eq '1' ? 0 : 1,
                -analysis => $runnable->analysis,
                );
            $genewise_runnable->run;
            foreach my $genewise_transcript (@{$genewise_runnable->output}) {
              attach_Slice_to_Transcript($genewise_transcript, $genewise_slice);
      print 'Working on ', $accession, ' ', $genewise_transcript->seq_region_start, ' ', $genewise_transcript->seq_region_end, ' ', $genewise_transcript->seq_region_strand, ' and ', $genewise_transcript->slice->name, "\n";
              $self->get_SeqEdit_by_transcript($genewise_transcript);
              my $seleno_transcript = $self->get_best_transcript($genewise_transcript, $exonerate_transcript, $seq);
              push(@selenocysteine_transcripts, $seleno_transcript) if ($seleno_transcript);
            }
          }
        }
      }
    }
    $self->output(\@selenocysteine_transcripts);
  }
  $self->dbc->disconnect_when_inactive(0);
}


sub get_best_transcript {
  my ($self, $genewise_transcript, $exonerate_transcript, $seq) = @_;

  my $SEQ_DISPLAY_LENGTH = 70;
  my $is_equal = 1;
  my $accession = $seq->id;
  if ($genewise_transcript->seq_region_start == $exonerate_transcript->seq_region_start
      and $genewise_transcript->seq_region_end == $exonerate_transcript->seq_region_end
      and $genewise_transcript->seq_region_strand == $exonerate_transcript->seq_region_strand) {
    my $exonerate_idx = 0;
    my @exonerate_exons = sort {$a->start <=> $b->start} @{$exonerate_transcript->get_all_Exons};
    foreach my $exon (sort {$a->start <=> $b->start} @{$genewise_transcript->get_all_Exons}) {
      if (defined $exonerate_exons[$exonerate_idx]) {
        if ($exon->seq_region_start != $exonerate_exons[$exonerate_idx]->seq_region_start or $exon->seq_region_end != $exonerate_exons[$exonerate_idx]->seq_region_end or $exon->seq_region_strand != $exonerate_exons[$exonerate_idx]->seq_region_strand) {
          $is_equal = 0;
          $self->warning("Exon $exonerate_idx created by Exonerate and Genewise for $accession differs!");
        }
      }
      else {
        $is_equal = 0;
        $self->warning("Exon $exonerate_idx created by Exonerate and Genewise for $accession differs!");
      }
      $exonerate_idx++;
    }
    if (defined $exonerate_exons[$exonerate_idx]) {
      $is_equal = 0;
      $self->warning("Exon $exonerate_idx created by Exonerate and Genewise for $accession differs!");

    }
  }
  else {
    $is_equal = 0;
    $self->warning("Transcripts created by Exonerate and Genewise for $accession differs!");
  }
  my $transcript_ori = $seq->{original_seq};
  my $exo_ori = $self->is_same_translation($exonerate_transcript, $transcript_ori, 'Exonerate');
  my $transcript_exo = $exonerate_transcript->translation->seq;
  my $transcript_genewise = $genewise_transcript->translation->seq;
  if ($is_equal) {
    if ($exo_ori >= 0) {
      $self->warning("GREAT MATCH!! The predicted proteins match the original protein $accession");
      return $exonerate_transcript;
    }
    else {
      $self->warning("The predicted proteins does not match the original protein $accession");
    }
  }
  else {
    my $genewise_ori = $self->is_same_translation($genewise_transcript, $transcript_ori, 'GeneWise');
    if ($exo_ori >= 0 and abs($genewise_ori) >= $exo_ori) {
      $self->warning("GOOD MATCH!! The Exonerate protein match the original protein $accession");
      display_alignment('EXONERATE', $transcript_exo, 'ORIGIN', $transcript_ori, $SEQ_DISPLAY_LENGTH) if ($self->debug > 2);
      return $exonerate_transcript;
    }
    elsif ($genewise_ori >= 0 and abs($exo_ori) >= $genewise_ori) {
# For some reason genewise has FeaturePair fo the transcript_supporting_features.
# We might need to chage it but there are some warning in the GeneWise module...
# So let's change it to a DnaPepAlignFeature
      $genewise_transcript->flush_supporting_features;
      my @sfs;
      foreach my $exon (@{$genewise_transcript->get_all_Exons}) {
        push(@sfs, @{$exon->get_all_supporting_features});
      }
      my $new_tsf = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@sfs, -align_type => 'ensembl');
      $genewise_transcript->add_supporting_features($new_tsf);
      $self->warning("GOOD MATCH!! The GeneWise protein match the original protein $accession");
      display_alignment('GENEWISE', $transcript_genewise, 'ORIGIN', $transcript_ori, $SEQ_DISPLAY_LENGTH) if ($self->debug > 2);
      return $genewise_transcript;
    }
    else {
      display_alignment('EXONERATE', $transcript_exo, 'GENEWISE', $transcript_genewise, $SEQ_DISPLAY_LENGTH) if ($self->debug > 1);
      $self->warning("The predicted protein does not match the original protein $accession");
    }

  }
  return;
}

sub write_output {
  my ($self) = @_;

  my $gene_adaptor = $self->hrdb_get_con('target_db')->get_GeneAdaptor;
  $gene_adaptor->dbc->disconnect_when_inactive(0);
  my $analysis = $self->analysis;
  my $biotype = $self->param('biotype');
  foreach my $transcript (@{$self->output}) {
    $transcript->biotype($biotype);
    my $gene = Bio::EnsEMBL::Gene->new();
    $gene->add_Transcript($transcript);
    $gene->biotype($biotype);
    attach_Analysis_to_Gene($gene, $analysis);
    attach_Slice_to_Gene($gene, $transcript->slice);
    $gene_adaptor->store($gene);
  }
}

sub prepare_sequences {
  my ($self, $sequences) = @_;

  my @seqs;
  foreach my $seq (@$sequences) {
    next unless ($seq->seq =~ /U/);
    my $protein = $seq->seq;
    my @positions;
    while($protein =~ /U/gc) {
      push(@positions, pos($protein));
    }
    $seq->{positions} = \@positions;
    $seq->{original_seq} = $protein;
    $protein =~ s/U/C/g;
    $seq->seq($protein);
    push(@seqs, $seq);
  }
  if (@seqs) {
    return \@seqs;
  }
  else {
    return;
  }
}

sub get_fasta_file_from_uniprot {
    my ($proteinfile) = @_;

    my $fasta = Bio::SeqIO->new( -file => $proteinfile, -format => 'fasta');
    my @sequences;
    while( my $seq = $fasta->next_seq) {
        my ($db_name, $accession) = $seq->id =~ /(\w{2})\|([^|]+)/;
        $seq->id($accession);
        push(@sequences, $seq);
    }
    return \@sequences;
}


sub mutate_dna {
  my ($self, $transcript, $slice) = @_;

  my $selenocysteines = $transcript->translation->get_all_SeqEdits('_selenocysteine');
  return unless (scalar(@$selenocysteines));
  my $sequence_length = 0;
  my $index = 0;
  my $dna = $slice->seq;
  
  foreach my $exon (@{$transcript->get_all_Exons}) {
    last unless (defined $selenocysteines->[$index]);
    $sequence_length += $exon->length;
    while (defined $selenocysteines->[$index] and $sequence_length > $selenocysteines->[$index]->start*3) {
      my $position = ($selenocysteines->[$index]->start*3)-2;
      my $exon_seq = $exon->seq->seq;
      my $exon_position = $position-($sequence_length-$exon->length);
# The exon sequence is on the translateable strand so we can search for TGA.
# But the genomic sequence will always be on the 5'->3', so we need to do math
# to get the right position on the reverse strand
      if (substr($exon_seq, $exon_position-1, 3) eq 'TGA') {
        if ($exon->strand == 1) {
          substr($dna, ($exon->seq_region_start+$exon_position+2)-$slice->start-1, 1, 'C');
        }
        else {
          substr($dna, (($exon->seq_region_end-$exon_position+1)-$slice->start)-2, 1, 'G');
        }
      }
      else {
        $self->warning("Could not find the stop coding for position: ".$selenocysteines->[$index]->start."($exon_position) for the protein ".$transcript->{accession}."\n".substr($exon_seq, $exon_position-5, 13));
        return;
      }
      $index++;
    }
  }
  return \$dna;
}

sub is_same_translation {
    my ($self, $tseq1, $seq2, $label) = @_;

    my $seq1_start = $tseq1->get_all_supporting_features->[0]->hstart-1;
    my $seq1 = $tseq1->translation->seq;
    my $seq2_start = 0;
    if (ref($seq2) eq 'Bio::EnsEMBL::Transcript') {
        $seq2_start = $seq2->get_all_supporting_features->[0]->hstart-1;
        $seq2 = $seq2->translation->seq;
    }
    return 0 if ($seq1 eq $seq2);
    my $missmatch = $seq1_start;
    my $longest = length($seq2);
    if (length($seq1) > $longest) {
      $longest = length($seq1);
    }
    else {
      $missmatch += $longest-length($seq1);
    }
    my @seq1_split = split('', $seq1);
    my @seq2_split = split('', $seq2);
    for (my $dix = 0; $dix < $longest; $dix++) {
        last unless (defined $seq1_split[$dix+$seq2_start] and defined $seq2_split[$dix+$seq1_start]);
        if ($seq1_split[$dix+$seq2_start] ne $seq2_split[$dix+$seq1_start]) {
            if ($seq1_split[$dix+$seq2_start] eq '*' or $seq2_split[$dix+$seq1_start] eq '*') {
                if ($seq1_split[$dix+$seq2_start] eq 'C' or $seq2_split[$dix+$seq1_start] eq 'C') {
                    $self->warning('You have a stop codon and a cysteine at position '.($dix+1).'. you probably want to edit this model');
                }
                else {
                    if ($dix-1 > (length($seq2)*($self->param('coverage_threshold')/100))) {
                      $self->warning('One of your sequence has a stop codon at position '.($dix+1).' and could have been saved');
                    }
                    else {
                      $self->warning('One of your sequence has a stop codon at position '.($dix+1));
                    }
                }
                return -100;
            }
            $missmatch++;
        }
    }
    if ($missmatch == 1 or ($missmatch/@seq1_split) < ($self->param('missmatch_allowed')/100)) {
        $self->warning("only $missmatch missmatches, your model from $label is OK");
        return $missmatch;
    }
    else {
        $self->warning("Too many missmatches ($missmatch)");
        return -$missmatch;
    }
}

sub display_alignment {
    my ($label1, $seq1, $label2, $seq2, $length) = @_;

    my $max_length = length($seq1) > length($seq2) ? length($seq1) : length($seq2);
    for (my $i = 0; $i < $max_length; $i = $i+$length) {
        print sprintf("% 10s: ", ''), sprintf("% 10s", '|')x($length/10), "\n", sprintf("% 10s: ", $label1), substr($seq1, $i, $length), "\n", sprintf("% 10s: ", $label2), substr($seq2, $i, $length), "\n";
    }
}


sub get_SeqEdit {
    my ($self, $transcript, $seq) = @_;

    my @attributes;
    my $pidx = 0;
    my @positions;
    my $translation = $transcript->translation->seq;
    my $tsf = $transcript->get_all_supporting_features->[0];
    my $accession = $tsf->hseqname;
    my $q_start = $tsf->hstart;
    my $cigar_string = $tsf->cigar_string;
    while ($translation =~ /\*/g) {
        push(@positions, pos($translation));
    }
    my $real_positions = $seq->{positions};
    for (my $idx = 0; $idx < @$real_positions; $idx++) {
        my $rpos = $real_positions->[$idx];
        my $adjusted_position = $rpos-$q_start+1;
        my $pos = $positions[$pidx];
        while (defined $positions[$pidx] and $rpos != $positions[$pidx]) {
          my $guessed_pos = $positions[$pidx]+$q_start-1;
          if ($cigar_string =~ /D/) {
            my $calc_pos = $q_start;
            while ($cigar_string =~ /(\d*)([MD])/gc) {
              if ($2 eq 'M') {
                $calc_pos += $1/3;
              }
              elsif ($2 eq 'D') {
                $calc_pos -= $1/3;
                $adjusted_position -= $1/3;
              }
              last unless ($calc_pos < $pos);
            }
          }
          if ($pos == $adjusted_position) {
              $self->warning("Stop codon found at position $pos for $accession");
              last;
          }
          else {
              $self->warning("Could not find the stop codon at position $pos");
              $pidx++;
          }
        }
        last unless (defined $positions[$pidx]);
        if ($adjusted_position == $pos) {
            my $seqedit = Bio::EnsEMBL::SeqEdit->new(
                                   -CODE    => '_selenocysteine',
                                   -NAME    => 'Selenocysteine',
                                   -DESC    => 'Selenocysteine',
                                   -START   => $pos,#$total_residues_before_stop_codon+1,
                                   -END     => $pos,#$total_residues_before_stop_codon+1,
                                   -ALT_SEQ => 'U'
                                   );
            push(@attributes, $seqedit->get_Attribute);
            $self->warning("We found a selenocysteine for $accession $pos\n");
        }
        else {
            $self->warning("We're missing some selenocysteine for $accession:$rpos $adjusted_position $pos\n");
        }
        $pidx++;
    }
    if (scalar(@attributes)) {
      $transcript->translation->add_Attributes(@attributes);
      return scalar(@attributes);
    }
    else {
      return;
    }
}

sub get_SeqEdit_by_transcript {
    my ($self, $transcript) = @_;

    my @attributes;
    my $cdna = $transcript->translateable_seq;
    while ($cdna =~ /\G(\w\w\w)*?TGA/g) {
        my $pos = pos($cdna)/3;
        my $seqedit = Bio::EnsEMBL::SeqEdit->new(
                -CODE    => '_selenocysteine',
                -NAME    => 'Selenocysteine',
                -DESC    => 'Selenocysteine',
                -START   => $pos,#$total_residues_before_stop_codon+1,
                -END     => $pos,#$total_residues_before_stop_codon+1,
                -ALT_SEQ => 'U'
                );
        push(@attributes, $seqedit->get_Attribute);
        $self->warning('We found a selenocysteine for '. $transcript->get_all_supporting_features->[0]->hseqname. " $pos\n");
    }
    $transcript->translation->add_Attributes(@attributes);
    return scalar(@attributes);
}

1;
