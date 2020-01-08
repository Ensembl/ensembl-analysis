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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadMitochondrion;

use strict;
use warnings;
use feature 'say';
use File::Fetch;
use File::Spec::Functions;
use File::Path qw(make_path);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;

  if($self->param_is_defined('skip_analysis') && $self->param('skip_analysis')) {
    $self->complete_early('The skip_analysis flag has been set to 1, skipping loading');
  }

  $self->param_required('target_db');
  $self->param_required('enscode_root_dir');
  $self->param_required('species_name');

  unless($self->param_is_defined('chromosomes_present')) {
    my $target_db = $self->get_database_by_name('target_db');
    my $cs_adaptor = $target_db->get_CoordSystemAdaptor;
    my $cs_rank1 = $cs_adaptor->fetch_by_rank(1);
    if ($cs_rank1->name eq 'primary_assembly' or $cs_rank1->name eq 'chromosome') {
      $self->param('chromosomes_present', $cs_rank1->name);
    }
  }

  if (!-e $self->param_required('output_path')) {
    make_path($self->param('output_path'));
  }

  if ($self->param_is_defined('mt_accession') and $self->param('mt_accession')) {
    chdir($self->param('output_path'));
    my $fetcher = File::Fetch->new(uri => 'http://www.ncbi.nlm.nih.gov/nuccore/'.$self->param('mt_accession'));
    $self->param('mt_filename', $fetcher->fetch());
  }
  else {
    $self->complete_early('No mitochondria for this species');
  }
}

sub run {
  my $self = shift;

  my $species_name = $self->param('species_name');
  my $output_path = $self->param('output_path');
  my $enscode_dir = $self->param('enscode_root_dir');

  my $mt_filename = $self->param('mt_filename');
  my $mt_accession = $self->param('mt_accession');
  my $toplevel = 'chromosome';
  my $chromosome_flag = "";
  my $scaffold_flag = "";
  if($self->param('chromosomes_present')) {
    $toplevel = 'primary_assembly' if ($self->param('chromosomes_present') eq 'primary_assembly');
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

  my $cmd =  'perl '.catfile($enscode_dir, 'ensembl-analysis', 'scripts', 'refseq', 'load_mitochondria.pl').
             " -dbhost ".$host.
             " -dbuser ".$user.
             " -dbport ".$port.
             " -dbpass ".$pass.
             " -dbname ".$dbname.
             $chromosome_flag.
             " -name MT".
             $scaffold_flag.
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
