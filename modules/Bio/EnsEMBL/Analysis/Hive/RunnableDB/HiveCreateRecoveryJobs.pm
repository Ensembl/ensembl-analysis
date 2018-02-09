#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateRecoveryJobs;

use strict;
use warnings;

use Bio::Seq;
use Bio::SeqIO;
use base ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');

sub fetch_input {
    my $self = shift;

    my @reads;
    my $filename;
    if ($self->param_is_defined('iid')) {
        if (-e $self->param('iid')) {
            $filename = $self->param('iid');
        }
        else {
            $self->throw('Input id '.$self->param('iid').' does not exists or is not a file');
        }
    }
    else {
        $filename = $self->param('wide_output_dir').'/'.$self->param('read_id_tag').'_sorted.bam';
    }
    my $command = join(' ', $self->param('wide_samtools'), 'view', $self->param('bam_flags'), $filename, '| ');
    open (BH, $command) || die('Could not open bam file '.$filename);
    my $batch_size = $self->param('batch_size');
    my $file_index = 1;
    my @output_ids;
    while (<BH>) {
        my @line = split("\t", $_);
        next if ($line[9] =~ y/N/N/ > 3);
#        next if ($line[10] =~ y/!"#$%&'()*+/!"#$%&'()*+/ > length($line[10])/3);
        next if ($line[10] =~ y/@ABCDEFGHIJKLMNOPQRSTUVWXYZ/@ABCDEFGHIJKLMNOPQRSTUVWXYZ/ > length($line[10])/3);
        my $mate;
        if ($line[1] & (1 << 6)) {
            $mate = 1;
        }
        else {
            $mate = 2;
        }
        my $read_group = grep {/RG:Z:/} @line;
        $read_group =~ s/RG:Z://;
        chomp($read_group);
        my $bioseq = Bio::Seq->new(
                -seq        => $line[9],
                -display_id => $line[1].'/'.$mate,
                -desc       => $read_group,
                );
        push(@reads, $bioseq);
        if (scalar(@reads) == $batch_size) {
            my ($name) = $filename =~ /([^\/]+)$/;
            my $file = $self->param('wide_recovery_dir').'/'.$name.'_'.$file_index.'.fa';
            my $fh = Bio::SeqIO->new(-format => 'fasta', -file => ">$file");
            $file_index++;
            while(my $seq = pop(@reads)) {
                $fh->write_seq($seq);
            }
            push(@output_ids, [$file]);
        }
    }
#    my $batch_size = (scalar(@reads)%$self->param('batch_size')) ? int(scalar(@reads)/int(scalar(@reads)/$self->param('batch_size'))+1) : $self->param('batch_size');
#    my $batch = 0;
#    my $file_index = 1;
#    my $fh;
#    my @output_ids;
#    while (pop(@reads)) {
#        if ($batch < $batch_size) {
#            $batch++;
#        }
#        else {
#            my $file = $self->param('wide_output_dir').'/'.$filename.'_'.$file_index.'.fa';
#            push(@output_ids, [$file]);
#            $fh = Bio::SeqIO->new(-format => 'fasta', -file => ">$file");
#            $file_index++;
#        }
#    }
    $self->param('column_names', ['iid']);
    $self->param('inputlist', \@output_ids);
}

1;
