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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivecDNAManyHits;

use strict;
use warnings;
use feature 'say';


#use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');



sub fetch_input {
  my $self = shift;
  return 1;
}

sub run {
  my ($self) = shift;

  my $output_path = $self->param('dest_dir') . '/' . $self->param('file_dir');
  my $qdb = $self->param('query_db');

  $self->find_many_hits($qdb);

  return 1;
}

sub write_output {
  my $self = shift;
  return 1;
}

sub find_many_hits {
  my ($self,$cdna_db) = @_;

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

  # Make a table containing each transcript matching a cDNA
  my $sql1 =
    (   "CREATE temporary table tmp1 "
      . "SELECT hit_name, exon_transcript.transcript_id "
      . "FROM dna_align_feature, supporting_feature, exon_transcript "
      . "WHERE dna_align_feature_id = feature_id "
      . "AND supporting_feature.exon_id = exon_transcript.exon_id "
      . "GROUP by hit_name, exon_transcript.transcript_id" );

  my $q1 = $db->dbc()->prepare($sql1) or croak("Sql error 1.\n$!");
  $q1->execute();

  # Group these to find the number of hits per cDNA
  my $sql2 = (  "CREATE temporary table tmp2 SELECT hit_name, "
              . "count(*) as hits FROM tmp1 GROUP by hit_name" );

  my $q2 = $db->dbc()->prepare($sql2) or croak("Sql error 2.\n$!");
  $q2->execute();

  # Examine those which hit more than 20 places
  my $sql3 = ("SELECT * FROM tmp2 WHERE hits > 20 ORDER by hits desc");

  my $q3 = $db->dbc()->prepare($sql3) or croak("Sql error 3\n$!");
  $q3->execute();

  my $many_hits_flag = 0;
  while ( my ( $cdna, $hits ) = $q3->fetchrow_array ) {
    print "$cdna\t$hits\n";
    $many_hits_flag = 1;
  }

  if ($many_hits_flag) {
     $self->throw ("\nIt might be worth investigating these sequences to "
          . "see whether these are likely to be genuine hits.\n"
          . "If we don't want them in the database you "
          . "should add them to the kill list\n\n" );
  }
} ## end sub find_many_hits


1;

