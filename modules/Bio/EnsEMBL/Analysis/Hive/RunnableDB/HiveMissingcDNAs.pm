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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMissingcDNAs;

use strict;
use warnings;
use feature 'say';


use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');



sub fetch_input {
  my $self = shift;
  return 1;
}

sub run {
  my ($self) = shift;

  my $output_path = $self->param('dest_dir');
  my $qdb = $self->param('query_db');
  my $cdna_file = $self->param('cdna_file');

  $self->find_missing_cdnas($qdb,$output_path,$cdna_file);

  return 1;
}

sub write_output {
  my $self = shift;
  return 1;
}

sub find_missing_cdnas {
  my ($self,$cdna_db,$output_dir,$clipped_file) = @_;

  # Mysql queries involving temporary tables
  my $dbhost = $cdna_db->{'-host'};
  my $dbport = $cdna_db->{'-port'};
  my $dbuser = $cdna_db->{'-user'};
  my $dbname = $cdna_db->{'-dbname'};
  my $dbpass = $cdna_db->{'-pass'};

  my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                                 -host    => $dbhost,
                                                 -port    => $dbport,
                                                 -user    => $dbuser,
                                                 -dbname  => $dbname,
                                                 -pass    => $dbpass,
                                                 );

  my $sql = ("SELECT distinct hit_name FROM dna_align_feature");

  my $q1 = $db->dbc->prepare($sql) or croak("Sql error.$!");
  $q1->execute();

  # Make list of cdnas with hits in the database.
  my (%cdna_hits);
  while ( my $cdna = $q1->fetchrow_array ) {
    $cdna_hits{$cdna} = 1;
  }

  # Now go through clipped sequence file and extract those sequences
  # which do not have any hits in the database.
  open( OUT, ">" . $output_dir . "/missing_cdnas.fasta" )
    or croak("Can't open file missing_cdnas.fasta");

  local $/ = "\n>";

  open( IN, "<$clipped_file" ) or croak("Can't open file $clipped_file!\n$!");
  while (<IN>) {
    my $seq = $_;

    if ( $seq =~ /(\w+\.\d+)\n/ ) {
      if ( !exists $cdna_hits{$1} ) {
        $seq =~ s/>//g;
        print OUT ">$seq\n";
      }
    }
  }
  close IN;
  close OUT;
}

1;

