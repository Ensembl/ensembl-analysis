#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadcDNAs;

use strict;
use warnings;

use Bio::SeqIO;
use Bio::EnsEMBL::Analysis::Tools::PolyAClipping qw(clip_if_necessary);

use parent ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults()},
    sequence_biotype => 'cdna',
    column_names => ['iid'],
    sequence_table_name => 'cdna_sequences',
    iid_type  => 'db_seq',
    format => 'fasta',
  }
}


sub fetch_input {
  my $self = shift;

  my $process_polyA = 0;
  my $parser = Bio::SeqIO->new(-format => $self->param('format'), -file => $self->param_required('cdna_file'));
  if ($self->param_is_defined('process_polyA') and $self->param('process_polyA')) {
    $process_polyA = 1;
  }
  my $biotype = $self->param('sequence_biotype');

  my $adaptor;
  my $write_to_file = $self->param('iid_type') eq 'db_seq' ? 0 : 1;
  if ($write_to_file) {
    $adaptor = Bio::SeqIO->new(-format => 'fasta', -file => '>'.$self->param_required('output_file'));
  }
  else {
    $adaptor = $self->db->get_NakedTableAdaptor();
    $adaptor->table_name($self->param('sequence_table_name'));
  }

  my @iids;
  while(my $bioseq = $parser->next_seq) {
    my $header = $bioseq->id;
    if ($process_polyA) {
      ($bioseq, undef, undef) = clip_if_necessary($bioseq);
      if (!$bioseq) {
        $self->warning('Sequence full of polyA for '.$header);
        next;
      }
    }

    $header =~ s/^\w*\|\w*\|//;
    if ($write_to_file) {
      $bioseq->id($header);
      $adaptor->write_seq($bioseq);
    }
    else {
      my $db_row = [{
        'accession'  => $header,
        'seq'        => $bioseq->seq,
        'biotype'    => $biotype,
      }];
      $adaptor->store($db_row);
    }
    push(@iids, $header);
  }
  $self->param('inputlist', \@iids);
}

1;
