#!/usr/bin/env perl

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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadcDNAs;

use strict;
use warnings;

use Bio::EnsEMBL::IO::Parser::Fasta;
use Bio::EnsEMBL::Analysis::Tools::PolyAClipping qw(clip_if_necessary);

use parent ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults()},
    sequence_biotype => 'cdna',
    column_names => ['iid'],
    sequence_table_name => 'cdna_sequences',
  }
}


sub fetch_input {
  my $self = shift;

  my $process_polyA = 0;
  my $parser = Bio::EnsEMBL::IO::Parser::Fasta->open($self->param_required('cdna_file'));
  if ($self->param_is_defined('process_polyA') and $self->param('process_polyA')) {
    $process_polyA = 1;
  }
  my $header;
  my $seq;
  my $biotype = $self->param('sequence_biotype');

  my $table_adaptor = $self->db->get_NakedTableAdaptor();
  $table_adaptor->table_name($self->param('sequence_table_name'));

  my @iids;
  while($parser->next()) {
    $header = $parser->getHeader();
    $seq = $parser->getSequence();
    $header =~ s/(\S+).*/$1/;
    if ($process_polyA) {
      my $bioseq = Bio::Seq->new(-id => $header, -seq => $seq);
      ($bioseq, undef, undef) = clip_if_necessary($bioseq);
      $seq = $bioseq->seq;
    }

    my $db_row = [{
      'accession'  => $header,
      'seq'        => $seq,
      'biotype'    => $biotype,
    }];
    $table_adaptor->store($db_row);
    push(@iids, $header);
  }
  $self->param('inputlist', \@iids);
}

1;
