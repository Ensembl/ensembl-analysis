#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HivePopulateProductionTables;

use strict;
use warnings;
use feature 'say';

use File::Spec::Functions;
use File::Path qw(make_path);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;

  unless($self->param('target_db')) {
    $self->throw("target_db flag not passed into parameters hash. The target db to load the assembly info ".
                 "into must be passed in with write access");
  }

  unless($self->param('production_db')) {
    $self->throw("You have used the populate_production_tables flag but have not passed in the production db connection hash with ".
          "the production_db flag");
  }


  unless($self->param('enscode_root_dir')) {
    $self->throw("enscode_root_dir flag not passed into parameters hash. You need to specify where your code checkout is");
  }

  unless($self->param('output_path')) {
    $self->throw("You have not specified the path to the main output directory with the -output_path flag, ".
                 "this is needed to dump the backup tables into");
  }
  my $dump_path = catdir($self->param('output_path'), 'populate_script_dump');
  make_path($dump_path) unless (-d $dump_path);
  $self->param('dump_path', $dump_path);

  return 1;
}

sub run {
  my $self = shift;


  my $target_db = $self->param('target_db');
  my $production_db = $self->param('production_db');
  my $enscode_dir = $self->param('enscode_root_dir');


  say "Running populate script on target db...\n";
  $self->populate_production_db_tables($target_db,$production_db,$enscode_dir,$self->param('dump_path'));
  say "...finished running script on target db\n";
  return 1;
}

sub write_output {
  my $self = shift;

  return 1;
}


sub populate_production_db_tables {
  my ($self,$target_db,$production_db,$enscode_dir,$dump_path) = @_;

  my $target_dbname;
  my $target_host;
  my $target_port;
  my $target_user;
  my $target_pass;

  $target_dbname = $target_db->{'-dbname'};
  $target_host = $target_db->{'-host'};
  $target_port = $target_db->{'-port'};
  $target_user = $target_db->{'-user'};
  $target_pass = $target_db->{'-pass'};


  unless($target_dbname && $target_host && $target_user && $target_pass) {
    $self->throw("The connection info for the database could not be recovered from the parameters passed in. ".
          "Make sure you pass in user_w, pass_w and either a db string or a hash to target_db");
  }

  my $production_host = $production_db->{'-host'};
  my $production_user = $production_db->{'-user'};
  my $production_port = $production_db->{'-port'};
  my $production_dbname = $production_db->{'-dbname'};

  my $populate_script = catfile($self->param('enscode_root_dir'), 'ensembl-production', 'scripts', 'production_database', 'populate_production_db_tables.pl');
  my $cmd = "perl ".$populate_script.
            " -h ".$target_host.
            " -u ".$target_user.
            " -p ".$target_pass.
            " -d ".$target_dbname.
            " -mh ".$production_host.
            " -md ".$production_dbname.
            " -mu ".$production_user.
            " -dB ".
            " -dp ".$dump_path;

  if($target_port) {
    $cmd .= " -P ".$target_port;
  }

  if($production_port) {
    $cmd .= " -mP ".$production_port;
  }

  my $return = system($cmd);
  if($return) {
    $self->throw("The populate_production_db_tables script returned a non-zero exit value. Commandline used:\n".$cmd);
  }

}


1;
