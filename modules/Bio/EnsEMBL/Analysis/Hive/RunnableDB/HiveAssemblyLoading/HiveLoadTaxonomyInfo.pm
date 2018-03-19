#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

#  unless($self->param('production_db')) {
#    throw("You have used the populate_production_tables flag but have not passed in the production db connection hash with ".
#          "the production_db flag");
#  }


  unless($self->param('enscode_root_dir')) {
    $self->throw("enscode_root_dir flag not passed into parameters hash. You need to specify where your code checkout is");
  }


#  unless($self->param('taxonomy_db')) {
#    throw("You have not passed in the taxonomy db connection hash with the taxonomy_db flag");
#  }

  return 1;
}

sub run {
  my $self = shift;

  my $target_db = $self->param('target_db');
#  my $taxonomy_db = $self->param('taxonomy_db');
  my $enscode_dir = $self->param('enscode_root_dir');
#  my $taxon_id = $self->param('taxon_id');

  say "Running load_taxonomy script on target db...\n";
  $self->load_taxonomy($target_db,$enscode_dir);
#  $self->load_taxonomy($target_db,$taxonomy_db,$enscode_dir,$taxon_id);
  say "...finished running script on target db\n";
  return 1;
}

sub write_output {
  my $self = shift;

  return 1;
}

sub load_taxonomy {
 my ($self,$target_db,$enscode_dir) = @_;
#  my ($self,$target_db,$taxonomy_db,$enscode_dir,$taxon_id) = @_;

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

#  my $taxonomy_host = $taxonomy_db->{'-host'};
#  my $taxonomy_user = $taxonomy_db->{'-user'};
#  my $taxonomy_port = $taxonomy_db->{'-port'};
#  my $taxonomy_dbname = $taxonomy_db->{'-dbname'};



  my $taxonomy_script = $self->param('enscode_root_dir')."/ensembl-production/scripts/production_database/populate_species_meta.pl";
  my $cmd = 'perl '.$taxonomy_script.
            ' -h '.$target_host.
            ' -u '.$target_user.
            ' -p '.$target_pass.
            ' -d '.$target_dbname.
            ' -P '.$target_port; # release not needed 
            # ' -release '.$ENV{ENSEMBL_RELEASE};

  if ($self->param_is_defined('production_db')) {
    $cmd .= ' -mh '.$self->param('production_db')->{-host}.
            ' -mu '.$self->param('production_db')->{-user}.
            ' -md '.$self->param('production_db')->{-dbname}.
            ' -mP '.$self->param('production_db')->{-port};
    $cmd .= ' -mp '.$self->param('production_db')->{-pass} if (defined $self->param('production_db')->{-pass});
  }
  if ($self->param_is_defined('taxonomy_db')) {
    $cmd .= ' -th '.$self->param('taxonomy_db')->{-host}.
            ' -tu '.$self->param('taxonomy_db')->{-user}.
            ' -td '.$self->param('taxonomy_db')->{-dbname}.
            ' -tP '.$self->param('taxonomy_db')->{-port};
    $cmd .= ' -tp '.$self->param('taxonomy_db')->{-pass} if (defined $self->param('taxonomy_db')->{-pass});
  }

  if(system($cmd)) {
    $self->throw("The load_taxonomy script returned a non-zero exit value. Commandline used:\n".$cmd);
  }

}


1;
