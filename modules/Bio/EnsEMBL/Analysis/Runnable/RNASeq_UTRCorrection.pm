=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::RNASeq_UTRCorrection

=head1 SYNOPSIS

  my $runnable =
    Bio::EnsEMBL::Analysis::Runnable::RNASeq_UTRCorrection->new(
    );

 $runnable->run;
 my @results = $runnable->output;

=head1 DESCRIPTION

Take RNASeq gene models as input to correct the UTR. It tries to find reads
with polyA tail to define the 3' UTR. In any other cases it search for a drop
in the read coverage. When there is no UTR it tries to create one.
Because we overpredict the length of the UTR, if a model already have UTR,
it can only be shorten but never extended. We need better data like 3' pulldown
data to be able to clearly define the UTR.

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Runnable::RNASeq_UTRCorrection;

use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(attach_Analysis_to_Gene);

use Bio::DB::Sam;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

sub new {
  my ( $class, @args ) = @_;
  my $self = $class->SUPER::new(@args);
  my ( $bam_file, $genes, $slice ) = rearrange( [qw(BAM_FILE GENES SLICE)],@args );
  $self->bam_file($bam_file);
  $self->genes($genes);
  $self->slice($slice);
  return $self;
}

sub run  {
    my ( $self) = @_;
    my $bam = Bio::DB::Bam->open($self->bam_file);
    my $header = $bam->header();
    my $slice = $self->slice;
    my ($tid, $start, $end) = $header->parse_region($slice->seq_region_name.':'.$slice->start.'-'.$slice->end);
    my $bam_index = Bio::DB::Bam->index_open($self->bam_file);
    $self->throw("Bam file " . $self->bam_file . "  not found \n") unless $bam_index;
    my $genes = $self->genes;
    my @new_genes;
    foreach my $gene (@$genes) {
        my $transcripts = $gene->get_all_Transcripts();
        $gene->flush_Transcripts;
        foreach my $transcript (@$transcripts) {
            # TODO: refactor this part of the code as I did it one strand at a time
            if ($transcript->strand == 1) {
                my $real_coding_start = $transcript->coding_region_start;
                my $real_coding_end = $transcript->coding_region_end;
                my $slice_shift = $slice->start-1;
                # The BAM file is in real coordinates so it's better to use them for everything...
                print STDOUT 'DEBUG: ', $gene->display_id, ' TRANSCRIPT: ', $transcript->display_id, ' SRS: ', $transcript->seq_region_start, ' SRE: ', $transcript->seq_region_end, ' % CRS: ', $real_coding_start, ' CRE: ', $real_coding_end, ' STRAND: ', $transcript->strand, "\n";
                my $possible_utr_start = $transcript->seq_region_start;
                my $possible_utr_end = $transcript->seq_region_end;
                if ($transcript->start == $real_coding_start) {
                    $possible_utr_start = $transcript->seq_region_start-200;
                }
                if ( $transcript->end == $real_coding_end) {
                    $possible_utr_end = $transcript->seq_region_end+200;
                }
                my $exons = $transcript->get_all_Exons;
                # flush_Exons removes strand, start, stop but it will be updated when we add the exons
                $transcript->flush_Exons;
                # First let's do the start
                print STDOUT '#### 5\' END  ####', "\n";
                my $read_count = 0;
                my $result = 0;
                my $coding_start = $exons->[0]->end < $real_coding_start ? $exons->[0]->seq_region_end : $real_coding_start+$slice_shift;
                my @callback_data = (\$result, \$read_count, $coding_start, 0, $exons->[0]->strand);
                $bam_index->fetch($bam, $tid, $coding_start, $coding_start+1, \&_process_reads, \@callback_data);
                print STDOUT 'DEBUG: CODING READ COUNT: ', $read_count, "\n";
                my $utr_start = 0;
                my $last_read_count = $read_count;
                my $max_read_count = $read_count;
                my $threshold = int($max_read_count/10) || 1;
                print STDOUT 'DEBUG: WORKING ON: ', $possible_utr_start, ' TO ', $coding_start, "\n";
                for (my $i = $coding_start-1; $i > $possible_utr_start; $i--) {
                    $read_count = 0;
                    $callback_data[2] = $i;
                    $bam_index->fetch($bam, $tid, $i, $i+1, \&_process_reads, \@callback_data);
                    if ($read_count < $threshold and $read_count > $last_read_count) {
                        last;
                    }
                    else {
                        $utr_start = $i unless ($read_count == $last_read_count);
                        $threshold = int($max_read_count/10) || 1 if ($read_count > $max_read_count);
                        $last_read_count = $read_count;
                    }
                }
                if ($utr_start) {
                    if ($utr_start < $exons->[0]->start and $exons->[0]->start != $real_coding_start) {
                        print STDOUT 'DEBUG: NO 5\' CHANGE (UTR EXTENSION)', "\n";
                    }
                    else {
                        $utr_start -= $slice_shift;
                        $exons->[0]->phase(-1);
                        $exons->[0]->start($utr_start);
                        print STDOUT 'DEBUG: NEW UTR 5\': ', $utr_start, "\n";
                    }
                }
                else {
                    print STDOUT 'DEBUG: NO 5\' CHANGE', "\n";
                }
                print STDOUT '#### 3\' END  ####', "\n";
                $read_count = 0;
                my $coding_end = $exons->[-1]->start > $real_coding_end ? $exons->[-1]->seq_region_start : $real_coding_end+$slice_shift;
                $callback_data[2] = $coding_end;
                $callback_data[3] = 1;
                $bam_index->fetch($bam, $tid, $coding_end-1, $coding_end, \&_process_reads, \@callback_data);
                print STDOUT 'DEBUG: CODING READ COUNT: ', $read_count, "\n";
                my $utr_end = 0;
                $last_read_count = $read_count;
                $max_read_count = $read_count;
                $threshold = int($max_read_count/10) || 1;
                print STDOUT 'DEBUG: WORKING ON: ', $coding_end, ' TO ', $possible_utr_end, ' ^ ', $threshold, "\n";
                for (my $i = $coding_end; $i < $possible_utr_end; $i++) {
                    $read_count = 0;
                    $callback_data[2] = $i;
                    $bam_index->fetch($bam, $tid, $i-1, $i, \&_process_reads, \@callback_data);
                    if ($result > 0) {
                        $utr_end = $result;
                        print STDOUT 'DEBUG: POLYA TAIL: ', $utr_end, "\n";
                        last;
                    }
                    if ($read_count < $threshold and $read_count > $last_read_count) {
                        last;
                    }
                    else {
                        $utr_end = $i-1 unless ($read_count == $last_read_count);
                        $threshold = int($max_read_count/10) || 1 if ($read_count > $max_read_count);
                        $last_read_count = $read_count;
                    }
                }
                if ($utr_end) {
                    $utr_end -= $slice_shift;
                    $exons->[-1]->end_phase(-1);
                    $exons->[-1]->end($utr_end);
                    print STDOUT 'DEBUG: NEW UTR 3\': ', $utr_end, "\n";
                }
                else {
                    print STDOUT 'DEBUG: NO 3\' CHANGE', "\n";
                }
                my $start_exon_index = 0;
                my $end_exon_index = 0;
                my $new_coding_start = 1;
                my $new_coding_end = 1;
                foreach my $exon (@$exons) {
                    $transcript->add_Exon($exon);
                    if ($exon->end < $real_coding_start) {
                        $start_exon_index++;
                        $end_exon_index++;
                    }
                    else {
                        if ($exon->start <= $real_coding_start) {
                            $new_coding_start = $real_coding_start-$exon->start+1;
                        }
                        if ($exon->end < $real_coding_end) {
                            $end_exon_index++;
                        }
                        elsif ($exon->start < $real_coding_end) {
                            $new_coding_end = $real_coding_end-$exon->start+1;
                        }
                    }
                }
                print STDOUT 'NEW TRANSCRIPT ',  $transcript->start, '(',$transcript->seq_region_start,') - ', $transcript->end, '(', $transcript->seq_region_end, ') ', $transcript->coding_region_start, ' TO ', $transcript->coding_region_end, "\n";
                print STDOUT 'DEBUG: EXON_START ', $start_exon_index, ' EXON_END ', $end_exon_index, "\n";
                    print STDOUT "DEBUG CREATE TRANSLATION ", $exons->[$start_exon_index]->start, ' ', $exons->[$start_exon_index]->end, ' ', $exons->[$start_exon_index]->strand, ' ', $exons->[$end_exon_index]->start, ' ', $exons->[$end_exon_index]->end, ' ', $exons->[$end_exon_index]->strand, "\n";
                    my $translation = Bio::EnsEMBL::Translation->new(
                            -start_exon => $exons->[$start_exon_index],
                            -end_exon => $exons->[$end_exon_index],
                            -seq_start => $new_coding_start,
                            -seq_end => $new_coding_end,
                            );
                    $transcript->translation($translation);
                print STDOUT 'NEW TRANSLATION ', $transcript->translation->start, '(',$transcript->translation->genomic_start,') - ', $transcript->translation->end, '(', $transcript->translation->genomic_end, ') ', "\n";
                $gene->add_Transcript($transcript);
                print STDOUT '#### DONE  ####', "\n";
            }
            else {
                my $real_coding_start = $transcript->coding_region_start;
                my $real_coding_end = $transcript->coding_region_end;
                my $slice_shift = $slice->start-1;
                # The BAM file is in real coordinates so it's better to use them for everything...
                print STDOUT 'DEBUG: ', $gene->display_id, ' TRANSCRIPT: ', $transcript->display_id, ' SRS: ', $transcript->seq_region_start, ' SRE: ', $transcript->seq_region_end, ' % CRS: ', $real_coding_start, ' CRE: ', $real_coding_end, ' STRAND: ', $transcript->strand, "\n";
                my $possible_utr_start = $transcript->seq_region_start;
                my $possible_utr_end = $transcript->seq_region_end;
                if ($transcript->start == $real_coding_start) {
                    $possible_utr_start = $transcript->seq_region_start-200;
                }
                if ( $transcript->end == $real_coding_end) {
                    $possible_utr_end = $transcript->seq_region_end+200;
                }
                my $exons = $transcript->get_all_Exons;
                # flush_Exons removes strand, start, stop but it will be updated when we add the exons
                $transcript->flush_Exons;
                # First let's do the start
                print STDOUT '#### 5\' END  ####', "\n";
                my $read_count = 0;
                my $result = 0;
                my $coding_start = $exons->[-1]->end < $real_coding_start ? $exons->[-1]->seq_region_end : $real_coding_start+$slice_shift;
                my @callback_data = (\$result, \$read_count, $coding_start, 1, $exons->[-1]->strand);
                $bam_index->fetch($bam, $tid, $coding_start-1, $coding_start, \&_process_reads, \@callback_data);
                print STDOUT 'DEBUG: CODING READ COUNT: ', $read_count, ' FOR ', $coding_start, "\n";
                my $utr_start = 0;
                my $last_read_count = $read_count;
                my $max_read_count = $read_count;
                my $threshold = int($max_read_count/10) || 1;
                print STDOUT 'DEBUG: WORKING ON: ', $possible_utr_start, ' TO ', $coding_start, "\n";
                for (my $i = $coding_start-1; $i > $possible_utr_start; $i--) {
                    $read_count = 0;
                    $callback_data[2] = $i;
                    $bam_index->fetch($bam, $tid, $i, $i+1, \&_process_reads, \@callback_data);
                    if ($result > 0) {
                        if ($result < $real_coding_start) {
                            $utr_start = $result;
                            print STDOUT 'DEBUG: POLYA TAIL: ', $utr_start, "\n";
                            last;
                        }
                        else {
                            $result = 0;
                            # I'm faking the read count as it is wrong for this iteration
                            $read_count = $last_read_count;
                            print STDOUT 'DEBUG: ARGH WRONG POLYA TAIL: ', "\n";
                        }
                    }
                    if ($read_count < $threshold and $read_count > $last_read_count) {
                        last;
                    }
                    else {
                        $utr_start = $i unless ($read_count == $last_read_count);
                        $threshold = int($max_read_count/10) || 1 if ($read_count > $max_read_count);
                        $last_read_count = $read_count;
                    }
                }
                if ($utr_start) {
                    $utr_start -= $slice_shift;
                    $exons->[-1]->phase(-1);
                    $exons->[-1]->start($utr_start);
                    print STDOUT 'DEBUG: NEW UTR 5\': ', $utr_start, "\n";
                }
                else {
                    print STDOUT 'DEBUG: NO 5\' CHANGE', "\n";
                }
                print STDOUT '#### 3\' END  ####', "\n";
                $read_count = 0;
                my $coding_end = $exons->[0]->start > $real_coding_end ? $exons->[0]->seq_region_start : $real_coding_end+$slice_shift;
                $callback_data[2] = $coding_end;
                $callback_data[3] = 0;
                $bam_index->fetch($bam, $tid, $coding_end-1, $coding_end, \&_process_reads, \@callback_data);
                print STDOUT 'DEBUG: CODING READ COUNT: ', $read_count, ' FOR ', $coding_end, "\n";
                my $utr_end = 0;
                $last_read_count = $read_count;
                $max_read_count = $read_count;
                $threshold = int($max_read_count/10) || 1;
                print STDOUT 'DEBUG: WORKING ON: ', $coding_end, ' TO ', $possible_utr_end, ' ^ ', $threshold, "\n";
                for (my $i = $coding_end; $i < $possible_utr_end; $i++) {
                    $read_count = 0;
                    $callback_data[2] = $i;
                    $bam_index->fetch($bam, $tid, $i-1, $i, \&_process_reads, \@callback_data);
                    if ($read_count < $threshold and $read_count > $last_read_count) {
                        last;
                    }
                    else {
                        $utr_end = $i-1 unless ($read_count == $last_read_count);
                        $threshold = int($max_read_count/10) || 1 if ($read_count > $max_read_count);
                        $last_read_count = $read_count;
                    }
                }
                if ($utr_end) {
                    if ($utr_end > $exons->[0]->end and $exons->[0]->end != $real_coding_end) {
                        print STDOUT 'DEBUG: NO 3\' CHANGE (UTR EXTENSION)', "\n";
                    }
                    else {
                        $utr_end -= $slice_shift;
                        $exons->[0]->end_phase(-1);
                        $exons->[0]->end($utr_end);
                        print STDOUT 'DEBUG: NEW UTR 3\': ', $utr_end, "\n";
                    }
                }
                else {
                    print STDOUT 'DEBUG: NO 3\' CHANGE', "\n";
                }
                my $start_exon_index = 0;
                my $end_exon_index = 0;
                my $new_coding_start = 1;
                my $new_coding_end = 1;
                foreach my $exon (@$exons) {
                    $transcript->add_Exon($exon);
                    if ($exon->start > $real_coding_end) {
                        $start_exon_index++;
                        $end_exon_index++;
                    }
                    else {
                        if ($exon->end >= $real_coding_end) {
                            $new_coding_start = $exon->end-$real_coding_end+1;
                        }
                        if ($exon->start > $real_coding_start) {
                            $end_exon_index++;
                        }
                        elsif ($exon->end >= $real_coding_start) {
                            $new_coding_end = $exon->end-$real_coding_start+1;
                        }
                    }
                }
                print STDOUT 'NEW TRANSCRIPT ',  $transcript->start, '(',$transcript->seq_region_start,') - ', $transcript->end, '(', $transcript->seq_region_end, ') ', $transcript->coding_region_start, ' TO ', $transcript->coding_region_end, "\n";
                print STDOUT 'DEBUG: EXON_START ', $start_exon_index, ' EXON_END ', $end_exon_index, "\n";
                    print STDOUT "DEBUG CREATE TRANSLATION ", $exons->[$start_exon_index]->start, ' ', $exons->[$start_exon_index]->end, ' ', $new_coding_start, ' ', $exons->[$end_exon_index]->start, ' ', $exons->[$end_exon_index]->end, ' ', $new_coding_end, "\n";
                    my $translation = Bio::EnsEMBL::Translation->new(
                            -start_exon => $exons->[$start_exon_index],
                            -end_exon => $exons->[$end_exon_index],
                            -seq_start => $new_coding_start,
                            -seq_end => $new_coding_end,
                            );
                    $transcript->translation($translation);
                print STDOUT 'NEW TRANSLATION ', $transcript->translation->start, '(',$transcript->translation->genomic_start,') - ', $transcript->translation->end, '(', $transcript->translation->genomic_end, ') ', "\n";
                $gene->add_Transcript($transcript);
                print STDOUT '#### DONE  ####', "\n";
            }
        }
        attach_Analysis_to_Gene($gene, $self->analysis);
        push(@new_genes, $gene);
    }
    foreach my $gene (@new_genes) {
        print STDOUT $gene->display_id, ' @ ', $gene->start, '(',$gene->seq_region_start,') - ', $gene->end, '(', $gene->seq_region_end, ')', "\n";
        foreach my $transcript (@{$gene->get_all_Transcripts}) {
            print STDOUT '  ', $transcript->display_id, ' @ ', $transcript->start, '(',$transcript->seq_region_start,') - ', $transcript->end, '(', $transcript->seq_region_end, ') ', $transcript->coding_region_start, ' == ', $transcript->coding_region_end, "\n";
            foreach my $exon (@{$transcript->get_all_Exons}) {
                print STDOUT '    ', $exon->display_id, ' @ ', $exon->start, '(',$exon->seq_region_start,') - ', $exon->end, '(', $exon->seq_region_end, ')', "\n";

            }
        }
    }
    $self->output(\@new_genes);
}

sub _process_reads {
    my ($read, $callback_data) = @_;

    return if ($read->aux_get('XS'));
    return if ($read->unmapped());
    my ($result, $read_count, $index, $is_cdna_3prime, $strand) = @$callback_data;
    $$read_count++;
    if ($is_cdna_3prime) {
        my $seq = $read->query->seq->seq;
        my $read_start = $read->start;
        return if ($index != $read_start);
        my $read_strand = $read->strand;
        my $read_name = $read->query->name;
        my $seq_len = int(length($seq)/10);
        return if ($read->aux_get('NM') < $seq_len);
        my $count = 0;
        my $position = 0;
        if ($strand == 1) {
            if ($read_strand == 1) {
                ($count) = $seq =~ /(A{$seq_len}A+).{0,3}$/;
                $position = $-[0];
            }
            else {
                ($count) = $seq =~ /(T{$seq_len}T+).{0,3}$/;
                $position = $-[0];
            }
            if (defined $position and (length($count) > ($seq_len*2))) {
                $$result = $read_start+$position-1;
                print STDOUT 'DEBUG: SEQ: ', $read_name, ' @ ', $read_start, ' % ', $seq, ' ^ ', $read_strand, "\n";
            }
        }
        else {
            if ($read_strand == 1) {
                ($count) = $seq =~ /^.{0,3}(A{$seq_len}A+)/;
                $position = $+[0];
            }
            else {
                ($count) = $seq =~ /^.{0,3}(T{$seq_len}T+)/;
                $position = $+[0];
            }
            if (defined $position and (length($count) > ($seq_len*2))) {
                $$result = $read_start+$position;
                print STDOUT 'DEBUG: SEQ: ', $read_name, ' @ ', $read_start, ' % ', $seq, ' ^ ', $read_strand, "\n";
            }
        }
    }
}



###########################################################
# containers

sub genes {
    my ($self, $genes) = @_;

    $self->{genes} = $genes if ($genes);
    return $self->{genes};
}

sub slice {
    my ($self, $slice) = @_;

    $self->{slice} = $slice if ($slice);
    return $self->{slice};
}

sub bam_file {
    my ($self, $bam_file) = @_;

    $self->{bam_file} = $bam_file if ($bam_file);
    return $self->{bam_file};
}

1;
