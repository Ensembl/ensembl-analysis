#!/usr/bin/env perl

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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveComparecDNAdbs;

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

  my $old_cdna_db = $self->param('old_cdna_db');
  my $new_cdna_db = $self->param('new_cdna_db');

  my $outfile = $self->param('output_file');

  $self->compare($old_cdna_db,$new_cdna_db,$outfile);

  return 1;
}

sub write_output {
  my $self = shift;
  return 1;
}

# Compare results to previous data as a health-check
# can also bsub a further function call for every chromosome
sub compare {
  my ($self,$old_cdna_db,$new_cdna_db,$outfile) = @_;

  my %chromosomes_old;
  my %chromosomes_new;
  my %hits_per_chrom_old;
  my %hits_per_chrom_new;

  my $hit_count_new = 0;
  my $hit_count_old = 0;

  my $old_dbhost = $old_cdna_db->{'-host'};
  my $old_dbport = $old_cdna_db->{'-port'};
  my $old_dbuser = $old_cdna_db->{'-user'};
  my $old_dbname = $old_cdna_db->{'-dbname'};

  my $new_dbhost = $new_cdna_db->{'-host'};
  my $new_dbport = $new_cdna_db->{'-port'};
  my $new_dbuser = $new_cdna_db->{'-user'};
  my $new_dbname = $new_cdna_db->{'-dbname'};
  my $new_dbpass = $new_cdna_db->{'-pass'};

  my $old_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                                 -host   => $old_dbhost,
                                                 -port   => $old_dbport,
                                                 -user   => $old_dbuser,
                                                 -dbname => $old_dbname,
                                                 );

  my $new_db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                                 -host    => $new_dbhost,
                                                 -port    => $new_dbport,
                                                 -user    => $new_dbuser,
                                                 -dbname  => $new_dbname,
                                                 -pass    => $new_dbpass,  
                                                 );

  my $query1 = "SELECT coord_system_id FROM coord_system WHERE name='chromosome' && attrib='default_version'";
  my $old_q1 = $old_db->dbc()->prepare($query1) or croak("Sql error!$!\n");;
  $old_q1->execute();
  my ($old_coord_system_id) = $old_q1->fetchrow_array;
  my $new_q1 = $new_db->dbc()->prepare($query1) or croak("Sql error!$!\n");;
  $new_q1->execute();
  my ($new_coord_system_id) = $new_q1->fetchrow_array;


  my $query2 = "SELECT seq_region_id, name FROM seq_region WHERE coord_system_id = "
               . $old_coord_system_id . " AND name not like '%NT%'";

  my $old_q2 = $old_db->dbc()->prepare($query2) or croak("Sql error!$!\n");;
  $old_q2->execute();
  while (my($seq_region_id,$name) = $old_q2->fetchrow_array ) {
    $chromosomes_old{$name} = $seq_region_id;
  }

  my $query3 = "SELECT seq_region_id, name FROM seq_region WHERE coord_system_id = "
               . $new_coord_system_id . " AND name not like '%NT%'";

  my $new_q3 = $new_db->dbc()->prepare($query3) or croak("Sql error!$!\n");;
  $new_q3->execute();
  while (my($seq_region_id,$name) = $new_q3->fetchrow_array ) {
    $chromosomes_new{$name} = $seq_region_id;
  }

  open (OUTFILE, ">$outfile");

  print OUTFILE "\nGetting hits per chromosome\n" . "\told\tnew\tdiff\n";


  # Check hits per chromosome
  my $query4 = "SELECT COUNT(distinct hit_name) FROM dna_align_feature daf, analysis a WHERE a.logic_name = 'cdna_update'"
            . " and a.analysis_id = daf.analysis_id and daf.seq_region_id = ?";

  my $old_q4 = $old_db->dbc->prepare($query4) or croak("Sql error!$!\n");

  my $new_q4 = $new_db->dbc->prepare($query4) or croak("Sql error!$!\n");

  my @sorted_chromosomes = sort keys %chromosomes_new;
    
  foreach my $chromosome (@sorted_chromosomes) {
    $old_q4->bind_param(1,$chromosomes_old{$chromosome});
    $new_q4->bind_param(1,$chromosomes_new{$chromosome});
    $old_q4->execute();
    $new_q4->execute();

    $hits_per_chrom_old{$chromosome} = $old_q4->fetchrow_array;
    $hits_per_chrom_new{$chromosome} = $new_q4->fetchrow_array;
    my $diff = $hits_per_chrom_new{$chromosome} - $hits_per_chrom_old{$chromosome};
    print OUTFILE "\n$chromosome:" . "\t"
          . $hits_per_chrom_old{$chromosome} . "\t"
          . $hits_per_chrom_new{$chromosome} . "\t"
          . $diff;
    $hit_count_old += $hits_per_chrom_old{$chromosome};
    $hit_count_new += $hits_per_chrom_new{$chromosome};
  }
  print OUTFILE "\n\nsum:" . "\t" . $hit_count_old . "\t" . $hit_count_new . "\n\n";

  close (OUTFILE);
} ## end sub compare

1;
