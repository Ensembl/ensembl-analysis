=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSequencesToFiles

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSequencesToFiles;

use strict;
use warnings;

use Bio::SeqIO;
use Bio::Seq;

use parent ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    format => 'fasta',
    chunk => 0,
    chunk_size => 10,
    column_names => ['filename'],
    inputlist => ['#filename#'],
  }
}


sub fetch_input {
  my ($self) = @_;

  my $sth = $self->db->dbc->prepare('SELECT accession, seq FROM '.$self->param_required('sequence_table_name'));
  $sth->execute();
  my $parser = Bio::SeqIO->new(-format => $self->param('format'), -file => '>'.$self->param_required('filename'));
  while (my $row = $sth->fetchrow_arrayref) {
    my $seq = Bio::Seq->new(-id => $row->[0], -seq => $row->[1]);
    $parser->write_seq($seq);
  }
}

1;
