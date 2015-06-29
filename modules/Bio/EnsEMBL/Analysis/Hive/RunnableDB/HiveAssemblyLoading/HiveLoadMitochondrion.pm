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

  unless($self->param('core_db')) {
    $self->throw("core_db flag not passed into parameters hash. The core db to load the assembly info ".
                 "into must be passed in with write access");
  }

  unless($self->param('enscode_dir')) {
    $self->throw("enscode_dir flag not passed into parameters hash. You need to specify where your code checkout is");
  }

  return 1;
}

sub run {
  my $self = shift;

  say "Loading meta information seq region synonyms into reference db\n";
  my $core_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(%{$self->param('core_db')});
  my $enscode_dir = $self->param('enscode_dir');
  my $output_path = $self->param('output_path');


  return 1;
}

sub write_output {
  my $self = shift;

  return 1;
}

1;
