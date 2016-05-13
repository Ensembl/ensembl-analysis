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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadMitochondrion;

use strict;
use warnings;
use feature 'say';

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;

  unless ( $self->param('target_db') ) {
    $self->throw(
       "target_db not passed into parameters hash. The core db to load the assembly info " . "into must be passed in with write access" );
  }

  unless ( $self->param('enscode_dir') ) {
    $self->throw("enscode_dir not passed into parameters hash. You need to specify where your code checkout is");
  }

  unless ( $self->param('species_name') ) {
    $self->throw("species_name not passed into parameters hash. You need to specify what species you're working on");
  }

  unless ( $self->param('output_path') ) {
    $self->throw("output not passed into parameters hash. You need to specify where the output dir will be");
  }

  return 1;
}

sub run {
  my $self = shift;

  my $species_name = $self->param('species_name');
  my $output_path  = $self->param('output_path');
  my $mito_file    = $output_path . '/' . $species_name . "_mito_query.txt";

  my $cmd =
    'wget http://www.ncbi.nlm.nih.gov/nuccore/?term=' . $species_name . '%5BPORGN%5D+AND+biomol_genomic%5BPROP%5D+AND+' .
    'gene_in_mitochondrion%5BPROP%5D+AND+srcdb_refseq%5BPROP%5D+AND+complete%5BPROP%5D -O ' . $mito_file;

  my $result = system($cmd);

  unless ( $result == 0 ) {
    $self->throw( "Command to query the genbank nucleotide core for the species mito failed. Commandline used:\n" . $cmd );
  }

  $cmd    = "grep 'NCBI Reference Sequence:' " . $mito_file;
  $result = "";
  $result = `$cmd`;

  unless ( $result =~ /NC\_\d+\.\d+/ ) {
    $self->throw( "Failed to parse the NC accession for the mito out of the mito output file. Commandline used\n" . $cmd );
  }

  my $mito_accession = $&;

  say "Found mito accession for " . $species_name . ": " . $mito_accession;

  return 1;
} ## end sub run

sub write_output {
  my $self = shift;

  return 1;
}

1;
