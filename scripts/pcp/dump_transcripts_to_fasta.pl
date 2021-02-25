#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

# Script to impute proteing coding status into pcp database
# Reads in and validates rnasamaba and cpc2 output files and creates intersection
# Uses this to update biotype (impute proteincodign stayus effectively)

use warnings;
use strict;

use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long qw(:config no_ignore_case);

my $coord_system = 'toplevel';
my $dna_dbname = 'homo_sapiens_core_101_38';
my $dna_user   = 'ensro';
my $dna_host   = 'mysql-ens-mirror-1.ebi.ac.uk';
my $dna_port   = '4240';
my $dna_pass;

my $options = GetOptions ("user|dbuser|u=s"	 => \$dna_user,
                          "host|dbhost|h=s"	 => \$dna_host,
                          "port|dbport|P=i"	 => \$dna_port,
                          "dbname|db|D=s"    => \$dna_dbname,
                          "dbpass|pass|p=s" => \$dna_pass);

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $dna_port,
  -user    => $dna_user,
  -host    => $dna_host,
  -dbname  => $dna_dbname,
  -pass    => $dna_pass);

my $slice_adaptor = $db->get_SliceAdaptor();
my $slices = $slice_adaptor->fetch_all($coord_system);
my $gene_adaptor = $db->get_GeneAdaptor();

while ( my $slice = shift @{$slices} ) {
  my $genes = $slice->get_all_Genes();
  while ( my $gene = shift @{$genes} ) {
    my @transcripts = @{$gene->get_all_Transcripts};
    while ( my $transcript = shift @transcripts ) {
      my $name =  ">". $transcript->seq_region_name . ":" . $transcript->seq_region_start . ":" . $transcr$
      print $name, "\n", $transcript->seq->seq, "\n";
    }
  }
}
