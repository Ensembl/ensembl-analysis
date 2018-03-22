=head1 LICENSE

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadIMGT

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadIMGT;

use strict;
use warnings;

use Bio::EnsEMBL::IO::Parser::Fasta;

use parent ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    min_seq_length => 10,
    column_names => ['iid'],
  }
}


sub fetch_input {
  my $self = shift;

  my $files = $self->param_required('iid');
  my $input_seq_count = 0;
  my $below_min_length_count = 0;
  my $contains_stop = 0;

  my $min_seq_length = $self->param('min_seq_length');

  my $table_adaptor = $self->db->get_NakedTableAdaptor();
  $table_adaptor->table_name($self->param_required('sequence_table_name'));
  if (ref($files) ne 'ARRAY') {
    $files = [$files];
  }
  my @iids;
  my $sth = $table_adaptor->prepare('UPDATE '.$self->param('sequence_table_name').' SET seq = ? WHERE accession = ?');
  foreach my $file_path (@$files) {
    $self->throw("The input id doesn't exist, offending path:\n$file_path")
      unless(-e $file_path);
    my $parser = Bio::EnsEMBL::IO::Parser::Fasta->open($file_path);

    while($parser->next()) {
      ++$input_seq_count;
      my $seq = $parser->getSequence();
      my $header = $parser->getHeader();
      if ($header =~ /([^\|]+)\|([^\|]+)\|(\w+ [^\|_ ]+)/) {
        my $region_name = $2;
        my $accession = $1.'_'.$region_name;
        my $db_row = $table_adaptor->fetch_by_dbID($accession);
        if ($db_row) {
          $db_row->{seq} .= $seq;
          $sth->bind_param(1, $db_row->{seq});
          $sth->bind_param(2, $accession);
          $sth->execute();
        }
        else {
          my $species = lc($3);
          $species =~ s/ /_/;
          $db_row = [{ 'accession'  => $accession,
                      'source_db'  => 'imgt',
                      'pe_level'   => 1,
                      'biotype'    => lc(substr($region_name,0,2)).'_'.substr($region_name,3,1),
                      'group_name' => $species,
                      'seq'        => $seq,
                       }];
          $table_adaptor->store($db_row);
        }
        push(@iids, $accession);
      }
    }
  }

  $self->warning("Total sequences: $input_seq_count\nSequences below $min_seq_length: $below_min_length_count\nSequences with stops: $contains_stop\nStored sequences: ".scalar(@iids)."\n");
  $self->param('inputlist', \@iids);
}

1;
