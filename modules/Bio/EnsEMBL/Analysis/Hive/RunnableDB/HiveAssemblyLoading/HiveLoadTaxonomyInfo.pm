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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadTaxonomyInfo;

use strict;
use warnings;
use feature 'say';

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;

  unless($self->param('target_db')) {
    $self->throw("target_db flag not passed into parameters hash. The target db to load the assembly info ".
                 "into must be passed in with write access");
  }

  unless($self->param('production_db')) {
    throw("You have used the populate_production_tables flag but have not passed in the production db connection hash with ".
          "the production_db flag");
  }


  unless($self->param('enscode_root_dir')) {
    $self->throw("enscode_root_dir flag not passed into parameters hash. You need to specify where your code checkout is");
  }

  unless($self->param('output_path')) {
    $self->throw("You have not specified the path to the main output directory with the -output_path flag, ".
                 "this is needed to dump the backup tables into");
  }

  unless($self->param('taxonomy_db')) {
    throw("You have not passed in the taxonomy db connection hash with the taxonomy_db flag");
  }

  return 1;
}

sub run {
  my $self = shift;

  my $target_db = $self->param('target_db');
  my $taxonomy_db = $self->param('taxonomy_db');
  my $enscode_dir = $self->param('enscode_root_dir');
  my $taxon_id = $self->param('taxon_id');

  say "Running load_taxonomy script on target db...\n";
  $self->load_taxonomy($target_db,$taxonomy_db,$enscode_dir,$taxon_id);
  say "...finished running script on target db\n";
  return 1;
}

sub write_output {
  my $self = shift;

  return 1;
}

sub load_taxonomy {
  my ($self,$target_db,$taxonomy_db,$enscode_dir,$taxon_id) = @_;


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
    throw("The connection info for the database could not be recovered from the parameters passed in. ".
          "Make sure you pass in user_w, pass_w and either a db string or a hash to target_db");
  }

  unless($target_dbname && $target_host && $target_user && $target_pass) {
    throw("The connection info for the database could not be recovered from the parameters passed in. ".
          "Make sure you pass in user_w, pass_w and either a db string or a hash to target_db");
  }

  my $taxonomy_host = $taxonomy_db->{'-host'};
  my $taxonomy_user = $taxonomy_db->{'-user'};
  my $taxonomy_port = $taxonomy_db->{'-port'};
  my $taxonomy_dbname = $taxonomy_db->{'-dbname'};

  my $taxonomy_script = $self->param('enscode_root_dir')."/ensembl-pipeline/scripts/load_taxonomy.pl";
  my $cmd = "perl ".$taxonomy_script.
            " -lcdbhost ".$target_host.
            " -lcdbuser ".$target_user.
            " -lcdbpass ".$target_pass.
            " -lcdbname ".$target_dbname.
            " -lcdbport ".$target_port.
            " -taxondbhost ".$taxonomy_host.
            " -taxondbport ".$taxonomy_port.
            " -taxondbname ".$taxonomy_dbname.
            " -taxon_id ".$taxon_id;
            # Note script has no flag for -taxonuser, maybe it should be added. Hardcoded to ensro at the moment

  my $return = system($cmd);
  if($return) {
    throw("The load_taxonomy script returned a non-zero exit value. Commandline used:\n".$cmd);
  }

}


1;
