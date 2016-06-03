#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDumpGenome;

use strict;
use warnings;
use feature 'say';

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;

  unless($self->param('target_db')) {
    $self->throw("target_db not passed into parameters hash. The core db to load the assembly info ".
                 "into must be passed in with write access. You need to pass in the connection hash with 'target_db'");
  }

  unless($self->param('species_name')) {
    $self->throw("species_name not passed into parameters hash. You need to specify what species you're working on with 'species_name'");
  }

  unless($self->param('output_path')) {
    $self->throw("Output path not passed into parameters hash. You need to specify where the output dir will be with 'output_path'");
  }

  unless($self->param('coord_system_name')) {
    $self->throw("Coord system name was not passed in using the 'coord_system_version' parameter");
  }

  unless($self->param('enscode_root_dir')) {
    $self->throw("enscode_dir not passed into parameters hash. You need to specify where your code checkout is");
  }

  return 1;
}

sub run {
  my $self = shift;

  my $species_name = $self->param('species_name');
  my $output_path = $self->param('output_path');
  my $output_file = $species_name."_softmasked_toplevel.fa";

  unless(-e $output_path) {
    `mkdir -p $output_path`;
    `lfs setstripe -c -1 $output_path`;
  }

  # The dumping script generally seems to ignore the output_dir flag for some reason, so change into the dir to be sure
  chdir($output_path);

  my $db_info = $self->param('target_db');
  my $target_host = $db_info->{'-host'};
  my $target_name = $db_info->{'-dbname'};
  my $coord_system_name = $self->param('coord_system_name');
  my $enscode_dir = $self->param('enscode_root_dir');
  my $repeat_logic_names = $self->param('repeat_logic_names');
  my $repeat_string = "";
  if(scalar(@{$repeat_logic_names})) {
    $repeat_string .= " -mask -softmask ";
    foreach my $repeat_logic_name (@{$repeat_logic_names}) {
      $repeat_string .= " -mask_repeat ".$repeat_logic_name;
    }
  }

  my $cmd = "perl ".$enscode_dir."/ensembl-analysis/scripts/sequence_dump.pl".
            " -dbuser ensro".
            " -dbport 3306".
            " -dbhost ".$target_host.
            " -dbname ".$target_name.
            " -coord_system_name ".$coord_system_name.
            " -toplevel".
            " -onefile".
            " -nonref".
            $repeat_string.
            " -output_dir ".$output_path.
            " -filename ".$output_file;

  say "Running command:\n".$cmd;

  my $result = system($cmd);

  unless($result == 0) {
    $self->throw("Command to dump the genome failed. Commandline used:\n".$cmd);
  }


  return 1;
}

sub write_output {
  my $self = shift;

  return 1;
}

1;
