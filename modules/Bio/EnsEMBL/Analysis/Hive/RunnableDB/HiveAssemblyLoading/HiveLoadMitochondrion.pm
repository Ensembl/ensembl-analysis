#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2017] EMBL-European Bioinformatics Institute
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
use File::Fetch;
use File::Spec::Functions;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;

  unless($self->param('target_db')) {
    $self->throw("target_db not passed into parameters hash. The core db to load the assembly info ".
                 "into must be passed in with write access");
  }

  unless($self->param('enscode_root_dir')) {
    $self->throw("enscode_root_dir not passed into parameters hash. You need to specify where your code checkout is");
  }

  unless($self->param('species_name')) {
    $self->throw("species_name not passed into parameters hash. You need to specify what species you're working on");
  }

  unless($self->param('output_path')) {
    $self->throw("output_path not passed into parameters hash. You need to specify where the output dir will be");
  }

  unless($self->param('mito_index_path')) {
    $self->throw("mito_index_path not passed into parameters hash. You need to specify where the index of mito accession is");
  }

  unless($self->param_is_defined('chromosomes_present')) {
    $self->throw("Need to pass in the chromosomes_present param to define the toplevel");
  }

  if ($self->param_is_defined('mt_accession')) {
    chdir($self->param('output_path'));
    my $fetcher = File::Fetch->new(uri => 'http://www.ncbi.nlm.nih.gov/nuccore/'.$self->param('mt_accession'));
    $self->param('mt_filename', $fetcher->fetch());
  }

  return 1;
}

sub run {
  my $self = shift;

  my $species_name = $self->param('species_name');
  my $output_path = $self->param('output_path');
  my $mito_index = $self->param('mito_index_path');
  my $enscode_dir = $self->param('enscode_root_dir');

  unless(-e $output_path) {
    my $return = system("mkdir -p ".$output_path);
    if($return) {
      $self->throw("Output dir didn't exist and failed to create it");
    }
  }

  my $mt_filename;
  my $mt_accession;
  if ($self->param_is_defined('mt_filename')) {
    $mt_filename = $self->param('mt_filename');
    $mt_accession = $self->param('mt_accession');
  }
  else {
    $self->throw("Could not locate the mito index on the path provided. Path:\n".$mito_index)
      unless (-e $mito_index);
    open(IN, $mito_index) || $self->throw("Could not open $mito_index");
    while(<IN>) {
      my ($accession,$index_species_name) = split('=',$_);
      chomp $index_species_name;
      if($index_species_name eq $species_name) {
        say "Found mito accession for ".$species_name.": ".$accession;
        $mt_filename = catfile($output_path, $accession.'.gb');
        $mt_accession = $accession;
        last;
      }
    }
    close IN;

    $self->throw("Could not fine species accession in the index file")
      unless ($mt_accession);
  }
  my $toplevel;
  my $chromosome_flag = "";
  my $scaffold_flag = "";
  my $contig_flag = " -contig ".$mt_accession;
  if($self->param('chromosomes_present')) {
    $toplevel = "chromosome";
    $chromosome_flag = " -chromosome  MT";
    $scaffold_flag = " -scaffold ".$mt_accession;
  } else {
    $toplevel = "scaffold";
    $scaffold_flag = " -scaffold MT";
  }

  my $target_db = $self->param('target_db');
  my $user = $target_db->{'-user'};
  my $pass = $target_db->{'-pass'};
  my $host = $target_db->{'-host'};
  my $port = $target_db->{'-port'};
  my $dbname = $target_db->{'-dbname'};

  my $cmd =  "perl ".$enscode_dir."/ensembl-pipeline/scripts/DataConversion/mitochondria/load_mitochondria.pl".
             " -dbhost ".$host.
             " -dbuser ".$user.
             " -dbport ".$port.
             " -dbpass ".$pass.
             " -dbname ".$dbname.
             $chromosome_flag.
             " -name MT".
             $scaffold_flag.
             $contig_flag.
             " -toplevel ".$toplevel.
             " -gene_type protein_coding".
             " -trna_type Mt_tRNA".
             " -rrna_type Mt_rRNA".
             " -codon_table 2".
             " -logic_name mt_genbank_import".
             " -source RefSeq".
             " -accession ".$mt_accession.
             " -genbank_file ".$mt_filename.
             " -non_interactive".
             " -download";

  my $return = system($cmd);
  if($return) {
    $self->throw("The load_mitochondria script returned a non-zero exit code. Commandline used:\n".$cmd);
  }

  return 1;
}

sub write_output {
  my $self = shift;

  return 1;
}

1;
