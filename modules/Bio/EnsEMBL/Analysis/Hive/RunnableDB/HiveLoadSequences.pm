=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadSequences

=head1 SYNOPSIS


=head1 DESCRIPTION

Base module to load data into Hive custom tables. Modules should override
create_row_data which will be called before loading data into the table.
the accession field will be return in an arrayref on branch 2

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadSequences;

use strict;
use warnings;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

=head2 fetch_input

 Arg [1]    : None
 Description: Check that sequence_file and table_name exists. Load the module
              needed to parse the file
 Returntype :
 Exceptions :

=cut

sub fetch_input {
  my $self = shift;

  $self->param_required('sequence_file');
  $self->param_required('sequence_table_name');
  $self->require_module($self->get_module_name($self->param('filetype')));
  my $parser = $self->get_module_name($self->param('filetype'))->open($self->param('sequence_file'));
  $self->param('seq_parser', $parser);
}


=head2 get_module_name

 Arg [1]    : String $filetype, type of the file
 Description: Return the name of the ensembl-io module to use to parse the file
              Filetype supported are:
                fasta
                genbank
 Returntype : String module name
 Exceptions : Throws if Arg[1] is not a supported filetype

=cut

sub get_module_name {
  my ($self, $filetype) = @_;

  my $module_name = 'Bio::EnsEMBL::IO::Parser::';
  if ($filetype eq 'fasta') {
    $module_name .= 'Fasta';
  }
  elsif ($filetype eq 'genbank') {
    $module_name .= 'Genbank';
  }
  else {
    $self->throw('Filetype unknown, you may need to create a Parser for '.$filetype);
  }
  return $module_name
}


=head2 run

 Arg [1]    : None
 Description: We don't want to store the sequences before writing as the file could be big
              so everything is done in write_output
 Returntype :
 Exceptions :

=cut

sub run {
  return 1;
}

=head2 write_output

 Arg [1]    : None
 Description: Write the accession, sequence and other information in the table 'table_name'
              It dataflows the sequence accession taken from the accession field of the row
              returned by create_row_data on channel 2, 'iid'
 Returntype : None
 Exceptions :

=cut

sub write_output {
  my $self = shift;

  my $table_name = $self->param_required('sequence_table_name');
  my $parser = $self->param('seq_parser');

  my $table_adaptor = $self->db->get_NakedTableAdaptor();
  $table_adaptor->table_name($table_name);

  my $branch_to_flow_to = $self->param('_branch_to_flow_to');
  while($parser->next()) {
    my $row = $self->create_row_data($parser);
    $table_adaptor->store($row);

    $self->dataflow_output_id({iid => [$row->[0]->{accession}]}, $branch_to_flow_to);
  }
}


=head2 create_row_data

 Arg [1]    : Bio::EnsEMBL::IO::Parser $parser, an already open file
 Description: Inheriting class should  implement this method, it should return
              an hash reference in an array reference. The data will be stored
              by the write_output method in the 'table_name' table of the Hive
              pipeline database.
                return [{
                  'accession'  => $parser->get_accession,
                  'seq'        => $parser->get_sequence,
                }];
 Returntype : Array ref
 Exceptions : Throws if the method has not been implemented

=cut

sub create_row_data {
  my ($self, $parser) = @_;

  $self->throw('This method has to be implemented in '.ref($self));
}

1;
