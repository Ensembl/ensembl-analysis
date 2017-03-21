=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_verbosity);

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
    update => 0,
  }
}


sub fetch_input {
  my ($self) = @_;

  logger_verbosity('INFO');
  $self->param_required('genome');
  $self->param_required('exonerate');
  $self->param_required('genewise');
  my $dnadb = $self->get_database_by_name('dna_db');
  my $db = $self->get_database_by_name('target_db', $dnadb);
  $self->hrdb_set_con($db, 'target_db');
  if ($self->param_is_defined('source_db')) {
    $self->hrdb_set_con($self->get_database_by_name('source_db', $dnadb), 'source_db');
  }
# I want to get the proper version by asking UniProt
  $self->create_analysis;

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
  $self->param('querys', $self->prepare_sequence($querys));
  if (!scalar(@{$self->param('querys')})) {
    $self->complete_early('No selenoprotein to use');
    $self->input_job->autoflow(0);
  }
}


sub run {
  my ($self) = @_;

  my $slice_adaptor = $self->hrdb_get_con('target_db')->get_SliceAdaptor;
  foreach my $seq (@{$self->param('querys')}) {
    my $alignments = $self->exonerate_sequences(write_seqfile($seq));

    my $transcripts = $self->selenocysteine_finder($alignments, $seq, $slice_adaptor);
    my $genewise_sequences = $self->mutate_genomic_sequence($transcripts, $slice_adaptor);

    my $selenocysteine_transcripts = $self->genewise_sequences($genewise_sequences, $seq, $alignments);
    $self->output($selenocysteine_transcripts);
  }
}


sub write_output {
  my ($self) = @_;

  my $gene_adaptor = $self->hrdb_get_con('target_db')->get_GeneAdaptor;
  my $update = $self->param('update');
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

sub prepare_sequence {
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
  return \@seqs;
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

sub exonerate_sequences {
    my ($self, $seleno_fasta) = @_;

    my $basic_options = '--model protein2genome --bestn 3 --ryo "RESULT: %S %pi %ql %tl %V\n" --softmasktarget true --showalignment false --showsugar false';
#    my $basic_options ||= "--showsugar false --showvulgar false --showalignment false --ryo \"RESULT: %S %pi %ql %tl %g %V\\n\" ";
    my $exonerate_cmd = $self->param('exonerate').' '.$basic_options.' --target '.$self->param('genome')." --query $seleno_fasta";
    $exonerate_cmd =~ s/(showalignment\s+)false/$1true/ if ($self->debug > 2);
#    print STDERR $exonerate_cmd, "\n" if ($self->debug);
    my %alignments;
    my $minimum_identity = $self->param('minimum_identity');
    my $coverage_threshold = $self->param('coverage_threshold')/100;
    my $score_threshold = $self->param('score_threshold')/100;
    open(my $fh, "$exonerate_cmd |") || die("Could not execute exonerate: $exonerate_cmd\n");
    while(my $line = <$fh>) {
        print $line if ($self->debug > 2);
        next unless ($line =~ s/^RESULT: //);
        my ($accession, $q_start, $q_end, $q_orient, $t_seqname, $t_start, $t_end, $t_strand, $score, $perc_ident, $q_length, $t_length) = split(' ', $line);
        ($t_start, $t_end) = ($t_end, $t_start) if ($t_strand eq '-');
        # correct for exonerate half-open coords
        $q_start += 1;
        $t_start += 1;
        if (exists $alignments{$accession}) {
            if ($score > ($alignments{$accession}{score}+$score*$score_threshold) or ($perc_ident > $alignments{$accession}{identity} and $perc_ident > $minimum_identity and (($q_end-$q_start) > $q_length*$coverage_threshold)) ) {
                $alignments{$accession} = {
                    t_seqname => $t_seqname,
                    t_start => $t_start,
                    t_end => $t_end,
                    t_strand => $t_strand,
                    score => $score,
                    identity => $perc_ident,
                    };
            }
        }
        else {
            $alignments{$accession} = {
                t_seqname => $t_seqname,
                t_start => $t_start,
                t_end => $t_end,
                t_strand => $t_strand,
                score => $score,
                identity => $perc_ident,
                };
        }
    }
    close($fh) || die("Could not close exonerate filehandle $?");
    return \%alignments;
}


sub selenocysteine_finder {
    my ($self, $alignments, $sequence, $slice_adaptor) = @_;

    my @transcripts;
    my $padding = $self->param('padding');
    my $exonerate = $self->param('exonerate');
    my $analysis = $self->analysis;
    foreach my $accession (keys %$alignments) {
        my $slice = $slice_adaptor->fetch_by_name($alignments->{$accession}->{t_seqname});
        my $gene_slice = $slice->sub_Slice($alignments->{$accession}{t_start}-$padding, $alignments->{$accession}{t_end}+$padding);
        print 'Working on ', $accession, ' and ', $gene_slice->name, "\n";
        $self->create_target_file('selenocysteine', 'fa');
        my $sub_slice = $self->param('target_file');
        my $fasta_serial = Bio::EnsEMBL::Utils::IO::FASTASerializer->new($sub_slice);
        $fasta_serial->print_Seq($gene_slice);
        my $runnable = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript->new(
            -program => $exonerate,
            -analysis => $analysis,
            -query_seqs => [$sequence],
            -query_type => 'protein',
            -target_file => $sub_slice->filename,
            -options => '--bestn 1 --exhaustive --model protein2genome',
            );
        $runnable->verbose(1) if ($self->debug > 2);
        if ($self->debug > 1) {
            my $basic_options = $runnable->options;
            $basic_options =~ s/(showalignment\s+)false/$1true/;
            $runnable->options($basic_options);
        }
        $runnable->run;
        foreach my $transcript (@{$runnable->output}) {
            attach_Slice_to_Transcript($transcript, $gene_slice);
            my ($is_ok, $attributes) = $self->get_SeqEdit($transcript->translate->seq, $sequence, $accession, $transcript->get_all_supporting_features->[0]->hstart, $transcript->get_all_supporting_features->[0]->cigar_string);
            if ($is_ok) {
                $transcript->translation->add_Attributes(@$attributes);
                push(@transcripts, $transcript);
            }
        }
    }
    return \@transcripts;
}

sub mutate_genomic_sequence {
    my ($self, $transcripts, $slice_adaptor) = @_;

    my $padding = 2*$self->param('padding');
    my %genewise_sequences;
    # For all our transcripts we try to find the genomic position of the start so we can change TGA to TGC
    # Then we can use GeneWise to predict a model.
    foreach my $transcript (@$transcripts) {
        my $do_push = 1;
        my $genewise_slice = $slice_adaptor->fetch_by_region('toplevel', $transcript->seq_region_name, $transcript->seq_region_start-$padding, $transcript->seq_region_end+$padding);
        my $genewise_sequence = $genewise_slice->seq;
        my $accession = $transcript->get_all_supporting_features()->[0]->hseqname;
        my $selenocysteines = $transcript->translation->get_all_SeqEdits('_selenocysteine');
        my $sequence_length = 0;
        my $index = 0;
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
                        substr($genewise_sequence, ($exon->seq_region_start+$exon_position+2)-$genewise_slice->start-1, 1, 'C');
                    }
                    else {
                        substr($genewise_sequence, (($exon->seq_region_end-$exon_position+1)-$genewise_slice->start)-2, 1, 'G');
                    }
                }
                else {
                    $self->warning("Could not find the stop coding for position: ".$selenocysteines->[$index]->start."(".$exon_position.") for the protein $accession\n".substr($exon_seq, $exon_position-5, 13));
                    $do_push = 0;
                }
                $index++;
            }
        }
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
        $edited_slice->{seq} = $genewise_sequence;
        $genewise_sequences{$accession}{edited_sequence} = $edited_slice;
        $genewise_sequences{$accession}{sequence} = $genewise_slice;
        $genewise_sequences{$accession}{transcript} = $transcript if ($do_push);
    }
    return \%genewise_sequences;
}

sub genewise_sequences {
    my ($self, $genewise_sequences, $seq, $alignments) = @_;

    my @selenocysteine_transcripts;
    my $analysis = $self->param('analysis');
    my $SEQ_DISPLAY_LENGTH = 70;
    foreach my $accession (keys %$genewise_sequences) {
        my $edited_slice = $genewise_sequences->{$accession}{edited_sequence};
        my $slice = $genewise_sequences->{$accession}{sequence};
        my $runnable = Bio::EnsEMBL::Analysis::Runnable::Genewise->new (
            -protein => $seq,
            -query => $edited_slice,
            -reverse => $alignments->{$accession}{t_strand} eq '+' ? 0 : 1,
            -analysis => $analysis,
            );
        $runnable->verbose(1) if ($self->debug > 2);
        $runnable->run;
        my $exonerate_transcript = $genewise_sequences->{$accession}{transcript};
        foreach my $genewise_transcript (@{$runnable->output}) {
            my $is_equal = 1;
            attach_Slice_to_Transcript($genewise_transcript, $slice);
            my ($is_ok, $attributes) = $self->get_SeqEdit_by_transcript($genewise_transcript);
            if ($is_ok) {
                $genewise_transcript->translation->add_Attributes(@$attributes);
            }
            if ($genewise_transcript->seq_region_start == $exonerate_transcript->seq_region_start and $genewise_transcript->seq_region_end == $exonerate_transcript->seq_region_end and $genewise_transcript->seq_region_strand == $exonerate_transcript->seq_region_strand) {
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
                foreach my $t ($genewise_transcript, $exonerate_transcript) {
                    print $t->get_all_supporting_features->[0]->hseqname, ': ', $t->seq_region_start, ' - ', $t->seq_region_end, ' @ ', $t->seq_region_strand, "\n";
                    my $idx = 1;
                    foreach my $e (@{$t->get_all_Exons}) {
                        print "\tE$idx ", $e->seq_region_start, ' - ', $e->seq_region_end, "\n";
                        $idx++;
                    }
                }
            }
            my $transcript_exo = $exonerate_transcript->translation->seq;
            my $transcript_genewise = $genewise_transcript->translation->seq;
            my $transcript_ori = $seq->{original_seq};
            my $exo_ori = $self->is_same_translation($exonerate_transcript, $transcript_ori, 'Exonerate');
            if ($is_equal) {
                if ($exo_ori >= 0) {
                    push(@selenocysteine_transcripts, $exonerate_transcript);
                    $self->warning("GREAT MATCH!! The predicted proteins match the original protein $accession");
                }
                else {
                    $self->warning("The predicted proteins does not match the original protein $accession");
                }
            }
            else {
                my $genewise_ori = $self->is_same_translation($genewise_transcript, $transcript_ori, 'GeneWise');
                if ($exo_ori >= 0 and abs($genewise_ori) >= $exo_ori) {
                    push(@selenocysteine_transcripts, $exonerate_transcript);
                    $self->warning("GOOD MATCH!! The Exonerate protein match the original protein $accession");
                    display_alignment('EXONERATE', $transcript_exo, 'ORIGIN', $transcript_ori, $SEQ_DISPLAY_LENGTH) if ($self->debug > 2);
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
                    my $new_tsf = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@sfs);
                    $genewise_transcript->add_supporting_features($new_tsf); 
                    push(@selenocysteine_transcripts, $genewise_transcript);
                    $self->warning("GOOD MATCH!! The GeneWise protein match the original protein $accession");
                    display_alignment('GENEWISE', $transcript_genewise, 'ORIGIN', $transcript_ori, $SEQ_DISPLAY_LENGTH) if ($self->debug > 2);
                }
                else {
                    display_alignment('EXONERATE', $transcript_exo, 'GENEWISE', $transcript_genewise, $SEQ_DISPLAY_LENGTH) if ($self->debug > 1);
                    $self->warning("The predicted protein does not match the original protein $accession");
                }

            }
        }
    }
    return \@selenocysteine_transcripts;
}

sub is_same_translation {
    my ($self, $seq1, $seq2, $label) = @_;

    my $seq1_start = 0;
    my $seq2_start = 0;
    if (ref($seq1) eq 'Bio::EnsEMBL::Transcript') {
        $seq1_start = $seq1->get_all_supporting_features->[0]->hstart-1;
        $seq1 = $seq1->translation->seq;
    }
    if (ref($seq2) eq 'Bio::EnsEMBL::Transcript') {
        $seq2_start = $seq2->get_all_supporting_features->[0]->hstart-1;
        $seq2 = $seq2->translation->seq;
    }
    return 0 if ($seq1 eq $seq2);
    print "$seq1_start ## $seq2_start\n";
    my $missmatch = $seq1_start;
#    my $longest = length($seq1) >= length($seq2) ? length($seq1) : length($seq2);
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
            print STDOUT ($dix+$seq2_start), ': ', $seq1_split[$dix+$seq2_start], ' %% ', ($dix+$seq1_start), ': ', $seq2_split[$dix+$seq1_start], "\n";
            if ($seq1_split[$dix+$seq2_start] eq '*' or $seq2_split[$dix+$seq1_start] eq '*') {
                if ($seq1_split[$dix+$seq2_start] eq 'C' or $seq2_split[$dix+$seq1_start] eq 'C') {
                    $self->warning('You have a stop codon and a cysteine at position '.($dix+1).'. you probably want to edit this model');
                }
                else {
                    $self->warning('One of your sequence has a stop codon at position '.($dix+1));
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
    my ($self, $translation, $seq, $accession, $q_start, $cigar_string) = @_;

    my $is_ok = 1;
    my @attributes;
    my $pidx = 0;
    my @positions;
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
#              $pos = $guessed_pos;
              $self->warning('Stop codon found at position '.$pos." for $accession");
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
            $self->warning("We found a selenocysteine for $accession", $pos,"\n");
        }
        else {
            $self->warning("We're missing some selenocysteine for $accession:$rpos $adjusted_position $pos\n");
            $is_ok = 0;
        }
        $pidx++;
    }
    return $is_ok, \@attributes;
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
        $self->warning("We found a selenocysteine for ". $transcript->get_all_supporting_features->[0]->hseqname. " $pos\n");
    }
    return 1, \@attributes;
}

1;
