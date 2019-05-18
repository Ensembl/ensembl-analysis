=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadmRNAs

=head1 SYNOPSIS


=head1 DESCRIPTION

Module to load mRNA into a customised table in the Hive database

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadmRNAs;

use strict;
use warnings;
use POSIX qw(strftime);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadSequences');

=head2 create_row_data

 Arg [1]    : Bio::EnsEMBL::IO::Parser object
 Description: It will truncate the header to have only the accession, then return
              the accession and the sequence to be stored in the table 'table_name'
              of the Hive pipeline database
 Returntype : Array ref
 Exceptions : None

=cut

sub create_row_data {
  my ($self, $parser) = @_;

  my ($accession) = $parser->getHeader =~ /^\s*(\S+)/;
  my $source = 'INSDC';
  if ($accession =~ /^NM/) {
    $source = 'RefSeq';
  }
  my $biotype = 'mRNA';
  my $date = strftime "%Y/%m/%d", localtime;
#  return [{accession => $accession, seq => $parser->getSequence, source => $source, biotype => $biotype, date => $date}];
  return [{accession => $accession, seq => $parser->getSequence, source_db => $source, biotype => $biotype, date => $date}];
}

1;

