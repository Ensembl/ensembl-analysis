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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadNewcDNAs;

use strict;
use warnings;
use feature 'say';
use Bio::SeqIO;
use Data::Dumper;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;

  return 1;
}

sub run {
  my $self = shift;

  return 1;
}

sub write_output {
  my $self = shift;

  my $db_path = $self->param('cdna_file');

  my $strategy = $self->param('strategy');

  my $seq_file = new Bio::SeqIO( -file => $db_path,
                                 -format => "Fasta",
                               );

  # create a hash that has the cdnas that have already been aligned previously as the key
  my %seen_cdna;

  my $old_cdna_db = $self->param('old_cdna_db');

  my $old_dbhost = $old_cdna_db->{'-host'};
  my $old_dbport = $old_cdna_db->{'-port'};
  my $old_dbuser = $old_cdna_db->{'-user'};
  my $old_dbname = $old_cdna_db->{'-dbname'};

  my $old_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                                 -host    => $old_dbhost,
                                                 -port    => $old_dbport,
                                                 -user    => $old_dbuser,
                                                 -dbname  => $old_dbname,
                                                 );

  my $new_cdna_db = $self->param('new_cdna_db');

  my $new_dbhost = $new_cdna_db->{'-host'};
  my $new_dbport = $new_cdna_db->{'-port'};
  my $new_dbuser = $new_cdna_db->{'-user'};
  my $new_dbname = $new_cdna_db->{'-dbname'};
  my $new_dbpass = $new_cdna_db->{'-pass'};

  my $new_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                                 -host    => $new_dbhost,
                                                 -port    => $new_dbport,
                                                 -user    => $new_dbuser,
                                                 -dbname  => $new_dbname,
                                                 -pass    => $new_dbpass,
                                                 );
  
  # only need to get the hit names if we're only doing an update and not complete
  if ($strategy eq 'update') {
    my $sql_1 = ("SELECT DISTINCT(hit_name) FROM dna_align_feature");
    my $query_1 = $old_db->dbc()->prepare($sql_1) or croak("Sql error\n$!");
    $query_1->execute();
    while ( my $cdna  = $query_1->fetchrow_array ) {
      $seen_cdna {$cdna} = 1;
    }
    my $sql_2 = ("SELECT DISTINCT(identifier) FROM unmapped_object");
    my $query_2 = $old_db->dbc()->prepare($sql_2) or croak("Sql error\n$!");
    $query_2->execute();
    while ( my $cdna = $query_2->fetchrow_array ) {
      $seen_cdna {$cdna} = 1;
    }
  }
  my $biotype = 'cdna';

  my $table_adaptor = $self->db->get_NakedTableAdaptor();
  $table_adaptor->table_name('cdna_sequences');

  my $gene_id_file = $self->param('retire_gene_file');
  open (GENEIDFILE, ">>$gene_id_file");

  while (my $seq = $seq_file->next_seq) {
    my $header = $seq->id;
    my $sequence = $seq->seq;
    $header =~ s/(\S+).*/$1/;

    if ($seen_cdna {$header}) {
      delete $seen_cdna {$header};  
    } else {
      my $output_hash = {};
      $output_hash->{'iid'} = $header;

      my $db_row = [{
        'accession'  => $header,
        'seq'        => $sequence,
        'biotype'    => $biotype,
      }];
      $table_adaptor->store($db_row);
    
      $self->dataflow_output_id($output_hash,2);
      delete $seen_cdna {$header};
    }
  }
  # now get a list of the retired genes
  if ($strategy eq 'update') {
    foreach my $key (keys %seen_cdna) { 
      my $get_gene_id_sql = ("SELECT DISTINCT(g.gene_id) from gene g, transcript t, transcript_supporting_feature tsf, dna_align_feature daf where g.gene_id = t.gene_id AND t.transcript_id = tsf.transcript_id AND feature_type = 'dna_align_feature' AND tsf.feature_id = daf.dna_align_feature_id AND daf.hit_name = '".$key."'");
      my $gene_id_query = $new_db->dbc()->prepare($get_gene_id_sql) or croak("Sql error\n$!");
      $gene_id_query->execute();

      while ( my $gene_id  = $gene_id_query->fetchrow_array ) {
        print GENEIDFILE $gene_id, "\n";  
      }
    }
    close GENEIDFILE;
  }
  return 1;
}

1;
