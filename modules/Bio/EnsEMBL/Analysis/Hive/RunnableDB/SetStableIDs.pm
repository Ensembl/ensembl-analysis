#!/usr/bin/env perl

# Copyright [2017] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::SetStableIDs;

use strict;
use warnings;
use feature 'say';
use Data::Dumper;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;

  return 1;
}


sub run {
  my $self = shift;

  if($self->param('mapping_required')) {
    $self->run_stable_id_mapping();
  } else {
    $self->generate_stable_ids();
  }

  return 1;
}


sub write_output {
  my ($self) = @_;

  return 1;
}


sub run_stable_id_mapping {
  my ($self) = @_;


}

sub generate_stable_ids {
  my ($self) = @_;

  my $target_db = $self->param('target_db');
  my $start = $self->param_required('id_start');
  my $cmd = 'perl '.$self->param_required('enscode_root_dir').'/ensembl/misc-scripts/generate_stable_ids.pl'.
            ' -user '.$target_db->{'-user'}.
            ' -pass '.$target_db->{'-pass'}.
            ' -host '.$target_db->{'-host'}.
            ' -port '.$target_db->{'-port'}.
            ' -dbname '.$target_db->{'-dbname'}.
            ' -start '.$start.' > '.$self->param('output_path').'/stable_ids.sql';
  my $result = system($cmd);
  if($result) {
    $self->throw("The stable id generation script failed. Commandline used:\n".$cmd);
  }

  $cmd = 'mysql -u'.$target_db->{'-user'}.
         ' -p'.$target_db->{'-pass'}.
         ' -h'.$target_db->{'-host'}.
         ' -P'.$target_db->{'-port'}.
         ' -D'.$target_db->{'-dbname'}.
         ' < '.$self->param('output_path').'/stable_ids.sql';
  $result = system($cmd);
  if($result) {
    $self->throw("Loading stable ids failed. Commandline used:\n".$cmd);
  }

}

1;
