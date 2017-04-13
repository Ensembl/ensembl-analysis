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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivecDNAUpdate::HivecDNAUpdate;

use strict;
use warnings;
use feature 'say';

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;

  unless($self->param('target_db')) {
    $self->throw("target_db flag not passed into parameters hash. The target db to store the results info ".
                 "into must be passed in with write access");
  }

  unless($self->param('pipeline_db')) {
    $self->throw("pipeline_db flag not passed into parameters hash. The pipeline db must be passed in with write access");
  }

  unless($self->param('working_dir')) {
    $self->throw("working_dir flag not passed into parameters hash. This should be the path to the working directory ".
                 "where cDNA sequences will be downloaded");
  }

  unless($self->param('enscode_root_dir')) {
    $self->throw("enscode_dir flag not passed into parameters hash. You need to specify where your code checkout is");
  }

  return 1;
}

sub run {
  my $self = shift;

  say "Preparing cDNA sequences";
  my $target_db = $self->param('target_db');
  my $path_to_files = $self->param('working_dir')."/".$self->param('species_name')."/".$self->param('primary_assembly_dir_name')."/AGP/";
  my $enscode_dir = $self->param('enscode_root_dir');

  $self->load_assembly($target_db,$path_to_files,$enscode_dir);

  say "Finished downloading contig files";
  return 1;
}

sub write_output {
  my $self = shift;

  return 1;
}

