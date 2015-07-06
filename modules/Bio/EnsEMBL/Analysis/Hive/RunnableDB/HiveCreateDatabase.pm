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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase;

use strict;
use warnings;
use feature 'say';

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

use Data::Dumper;

sub param_defaults {
    return {

      # used by all create types
      create_type => '',
      source_db => '',
      target_db => '',

      # used by create_type = 'clone'
      script_path => '~/enscode/ensembl-personal/genebuilders/scripts/clone_database.ksh',

      # used by create_type = 'copy'
      db_dump_file => "/tmp/source_db_".time().".tmp",
      pass_r => '',
      pass_w => '',
      user_r => '',
      user_w => '',

    }
}

sub fetch_input {
  my $self = shift;
  return 1;
}

sub run {
  my $self = shift;

  unless(defined($self->param('create_type'))) {
    throw("You have not defined a create_type flag in you parameters hash! You must provide a create_type (e.g. 'clone')");
  }

  $self->create_db();

  if($self->param('populate_production_tables')) {
    $self->populate_production_db_tables();
  }

  if($self->param('load_taxonomy')) {
    $self->load_taxonomy();
  }
  return 1;
}

sub write_output {
  my $self = shift;

  return 1;
}

sub create_db {
  my $self = shift;
  my $create_type = $self->param('create_type');
  if($create_type eq 'clone') {
    $self->clone_db();
  } elsif ($create_type eq 'copy') {
    $self->copy_db();
  } elsif($create_type eq 'core_only') {
    $self->core_only_db();
  } else {
    throw("You have specified a create type of ".$create_type.", however this is not supported by the module");
  }

}

sub clone_db {
  my $self = shift;

  unless($self->param('source_db') && $self->param('target_db')) {
    throw("You have specified a create type of clone but you don't have both a source_db and target_db specified in your config");
  }

  unless($self->param('script_path')) {
    throw("You have specified a create type of clone but you don't have the script_path set in your config, e.g.:\n".
          "~/enscode/ensembl-personal/genebuilders/scripts/clone_database.ksh");
  }

  unless(-e $self->param('script_path')) {
    throw("The path to the clone script you have specified does not exist. Path:\n".$self->param('script_path'));
  }

  my $source_string;
  my $target_string;
  if(ref($self->param('source_db')) eq 'HASH') {
    $source_string = $self->convert_hash_to_db_string($self->param('source_db'));
  } else {
    $self->check_db_string($self->param('source_db'));
    $source_string = $self->param('source_db');
  }

  if(ref($self->param('target_db')) eq 'HASH') {
    $target_string = $self->convert_hash_to_db_string($self->param('target_db'));
  } else {
    $self->check_db_string($self->param('target_db'));
    $target_string = $self->param('target_db');
  }

  my $command = "ksh ".$self->param('script_path')." -s ".$source_string." -t ".$target_string;
  say "COMMAND: ".$command;

  my $exit_code = system($command);
  if($exit_code) {
    throw("The clone script exited with a non-zero exit code (did the target_db already exist?): ".$exit_code);
  }

}

sub copy_db {
  my $self = shift;

  if (not $self->param('source_db') or not $self->param('target_db')) {
    throw("You have specified a create type of copy but you don't have both a source_db and target_db specified in your config.");
  }

  if (not $self->param('user_r') or not $self->param('user_w') or not $self->param('pass_w')) {
    throw("You have specified a create type of copy but you haven't specified the user_r, user_w and pass_w.\n");
  }

  if (not $self->param('db_dump_file')) {
    throw("You have specified a create type of copy but you haven't specified a db_dump_file which will be used as a temporary file.");
  }

  my $source_string;
  my $target_string;
  if (ref($self->param('source_db')) eq 'HASH') {
    $source_string = $self->convert_hash_to_db_string($self->param('source_db'));
  } else {
    $self->check_db_string($self->param('source_db'));
    $source_string = $self->param('source_db');
  }

  if (ref($self->param('target_db')) eq 'HASH') {
    $target_string = $self->convert_hash_to_db_string($self->param('target_db'));
  } else {
    $self->check_db_string($self->param('target_db'));
    $target_string = $self->param('target_db');
  }

  my @source_string_at_split = split('@',$source_string);
  my $source_dbname = shift(@source_string_at_split);
  my @source_string_colon_split = split(':',shift(@source_string_at_split));
  my $source_host = shift(@source_string_colon_split);
  my $source_port = shift(@source_string_colon_split);

  my @target_string_at_split = split('@',$target_string);
  my $target_dbname = shift(@target_string_at_split);
  my @target_string_colon_split = split(':',shift(@target_string_at_split));
  my $target_host = shift(@target_string_colon_split);
  my $target_port = shift(@target_string_colon_split);

  dump_database($source_host,$source_port,$self->param('user_r'),$self->param('pass_r'),$source_dbname,$self->param('db_dump_file'));
  create_database($target_host,$target_port,$self->param('user_w'),$self->param('pass_w'),$target_dbname);
  load_database($target_host,$target_port,$self->param('user_w'),$self->param('pass_w'),$target_dbname,$self->param('db_dump_file'));
  remove_file($self->param('db_dump_file'));
}

sub core_only_db {
  my $self = shift;

  unless ($self->param('target_db')) {
    throw("You have specified a create type of core_only but you don't have a target_db specified in your config.");
  }

  unless($self->param('user_w') && $self->param('pass_w')) {
    throw("You have specified a create type of core_only but you haven't specified the user_w and pass_w.\n");
  }

  unless($self->param('enscode_root_dir')) {
    throw("You have specified a create type of core_only but you haven't specified the enscode_root_dir path\n");
  }

  my $table_file = $self->param('enscode_root_dir')."/ensembl/sql/table.sql";
  unless(-e $table_file) {
    throw("You have specified a create type of core_only but the path from enscode_root_dir to the table file is incorrect:\n".
          $self->param('enscode_root_dir')."/ensembl/sql/table.sql");
  }

  my $target_string;
  if (ref($self->param('target_db')) eq 'HASH') {
    $target_string = $self->convert_hash_to_db_string($self->param('target_db'));
  } else {
    $self->check_db_string($self->param('target_db'));
    $target_string = $self->param('target_db');
  }

  my @target_string_at_split = split('@',$target_string);
  my $target_dbname = shift(@target_string_at_split);
  my @target_string_colon_split = split(':',shift(@target_string_at_split));
  my $target_host = shift(@target_string_colon_split);
  my $target_port = shift(@target_string_colon_split);
  my $target_user = $self->param('user_w');
  my $target_pass = $self->param('pass_w');

  my $command;
  # Create the empty db
  if($target_port) {
    $command = "mysql -h".$target_host." -u".$target_user." -p".$target_pass." -P".$target_port." -e 'CREATE DATABASE ".$target_dbname."'";
  } else {
    $command = "mysql -h".$target_host." -u".$target_user." -p".$target_pass." -e 'CREATE DATABASE ".$target_dbname."'";
  }
  say "COMMAND: ".$command;

  my $exit_code = system($command);
  if($exit_code) {
    throw("The create database command exited with a non-zero exit code: ".$exit_code);
  }

  # Load core tables
  if($target_port) {
    $command = "mysql -h".$target_host." -u".$target_user." -p".$target_pass." -P".$target_port." -D".$target_dbname." < ".$table_file;
  } else {
    $command = "mysql -h".$target_host." -u".$target_user." -p".$target_pass." -D".$target_dbname." < ".$table_file;
  }
  $exit_code = system($command);
  if($exit_code) {
    throw("The load tables command exited with a non-zero exit code: ".$exit_code);
  }

}

sub convert_hash_to_db_string {
  my ($self,$connection_info) = @_;

  unless(defined($connection_info->{'-dbname'}) && defined($connection_info->{'-host'})) {
    throw("You have passed in a hash as your db info however the hash is missing either the dbname or host key,".
          " both are required if a hash is being passed in");
  }

  my $port = $connection_info->{'-port'};

  my $db_string;
  $db_string = $connection_info->{'-dbname'}.'@'.$connection_info->{'-host'};
  if(defined($port)) {
    $db_string .= ':'.$port;
  }

  $self->check_db_string($db_string);
  return($db_string);
}

sub check_db_string {
  my ($self,$db_string) = @_;

  unless($db_string =~ /[^\@]+\@[^\:]+\:\d+/ || $db_string =~ /[^\@]+\@[^\:]+/) {
    throw("Parsing check on the db string failed.\ndb string:\n".$db_string.
          "\nExpected format:\nmy_db_name\@myserver:port or my_db_name\@myserver");
  }
}

sub dump_database {

  my ($dbhost,$dbport,$dbuser,$dbpass,$dbname,$db_file) = @_;

  print "\nDumping database $dbname"."@"."$dbhost:$dbport...\n";

  my $command;
  if (!$dbpass) { # dbpass for read access can be optional
  	$command = "mysqldump --skip-opt -h$dbhost -P$dbport -u$dbuser $dbname > $db_file";
  } else {
  	$command = "mysqldump --skip-opt -h$dbhost -P$dbport -u$dbuser -p$dbpass $dbname > $db_file";
  }

  if (system($command)) {
    throw("The dump was not completed. Please, check that you have enough disk space in the output path $db_file as well as writing permission.");
  } else {
    print("The database dump was completed successfully into file $db_file\n");
  }
}

sub create_database {
  my ($dbhost,$dbport,$dbuser,$dbpass,$dbname) = @_;
  print "Creating database $dbname"."@"."$dbhost:$dbport...\n";
  if (system("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -e'CREATE DATABASE $dbname' ")) {
    throw("Couldn't create database  $dbname"."@"."$dbhost:$dbport. Please, check that it does not exist and you have write access to be able to perform this operation.");
  } else {
    print("Database $dbname"."@"."$dbhost:$dbport created successfully.\n");
  }
}

sub load_database {

  my ($dbhost,$dbport,$dbuser,$dbpass,$dbname,$db_file) = @_;

  print "\nLoading file $db_file into database $dbname"."@"."$dbhost:$dbport...\n";
  if (system("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname < $db_file")) {
    throw("The database loading process failed. Please, check that you have access to the file $db_file and the database you are trying to write to.");
  } else {
    print("\nThe database loading process was completed successfully from file $db_file into $dbname.\n");
  }
}

sub remove_file {
  my $db_file = shift;

  if (-e $db_file) {
  	print "Deleting file $db_file\n";
  	if (system("rm -f $db_file")) {
  	  throw("Couldn't delete file $db_file");
  	} else {
  	  print "File $db_file has been deleted.\n";
  	}
  }
}

sub populate_production_db_tables {
  my ($self) = @_;

  unless($self->param('enscode_root_dir')) {
    throw("You have used the populate_production_tables flag but have not passed in the location of the enscode dir with ".
          "the enscode_root_dir flag");
  }

  unless($self->param('output_path')) {
    throw("You have used the populate_production_tables flag but have not passed in the location of the output dir with ".
          "the output_path flag");
  }

  unless(-d $self->param('output_path')) {
    `mkdir -p $self->param('output_path')`;
  }

  my $target_dbname;
  my $target_host;
  my $target_port;
  my $target_user;
  my $target_pass;

  my $target_string;
  if (ref($self->param('target_db')) eq 'HASH') {
    my $db_hash = $self->param('target_db');
    $target_dbname = $db_hash->{'-dbname'};
    $target_host = $db_hash->{'-host'};
    $target_port = $db_hash->{'-port'};
    $target_user = $self->param('user_w');
    $target_pass = $self->param('pass_w');
  } else {
    $self->check_db_string($self->param('target_db'));
    $target_string = $self->param('target_db');
    my @target_string_at_split = split('@',$target_string);
    $target_dbname = shift(@target_string_at_split);
    my @target_string_colon_split = split(':',shift(@target_string_at_split));
    $target_host = shift(@target_string_colon_split);
    $target_port = shift(@target_string_colon_split);
    $target_user = $self->param('user_w');
    $target_pass = $self->param('pass_w');
  }

  unless($target_dbname && $target_host && $target_user && $target_pass) {
    throw("The connection info for the database could not be recovered from the parameters passed in. ".
          "Make sure you pass in user_w, pass_w and either a db string or a hash to target_db");
  }

  my $production_db = $self->param('production_db');
  unless($production_db) {
    throw("You have used the populate_production_tables flag but have not passed in the production db connection hash with ".
          "the production_db flag");
  }

  my $production_host = $production_db->{'-host'};
  my $production_user = $production_db->{'-user'};
  my $production_port = $production_db->{'-port'};
  my $production_dbname = $production_db->{'-dbname'};

  my $populate_script = $self->param('enscode_root_dir')."/ensembl-production/scripts/production_database/populate_production_db_tables.pl";
  my $cmd = "perl ".$populate_script.
            " -h ".$target_host.
            " -u ".$target_user.
            " -p ".$target_pass.
            " -d ".$target_dbname.
            " -mh ".$production_host.
            " -md ".$production_dbname.
            " -mu ".$production_user.
            " -dp ".$self->param('output_path');

  if($target_port) {
    $cmd .= " -P ".$target_port;
  }

  if($production_port) {
    $cmd .= " -mP ".$production_port;
  }

  my $return = system($cmd);
  if($return) {
    throw("The populate_production_db_tables script returned a non-zero exit value. Commandline used:\n".$cmd);
  }

}

sub load_taxonomy {
  my ($self) = @_;

  unless($self->param('enscode_root_dir')) {
    throw("You have used the populate_production_tables flag but have not passed in the location of the enscode dir with ".
          "the enscode_root_dir flag");
  }

  unless($self->param('taxon_id') || $self->param('taxon_name')) {
    throw("You have used the load_taxonomy parameter but have not passed in the taxon id or taxon name using the ".
          "taxon_id or taxon_name flag in the config");
  }

  my $taxon_id = $self->param('taxon_id');
  my $taxon_name = $self->param('taxon_name');
  my $target_dbname;
  my $target_host;
  my $target_port;
  my $target_user;
  my $target_pass;

  my $target_string;
  if (ref($self->param('target_db')) eq 'HASH') {
    my $db_hash = $self->param('target_db');
    $target_dbname = $db_hash->{'-dbname'};
    $target_host = $db_hash->{'-host'};
    $target_port = $db_hash->{'-port'};
    $target_user = $self->param('user_w');
    $target_pass = $self->param('pass_w');
  } else {
    $self->check_db_string($self->param('target_db'));
    $target_string = $self->param('target_db');
    my @target_string_at_split = split('@',$target_string);
    $target_dbname = shift(@target_string_at_split);
    my @target_string_colon_split = split(':',shift(@target_string_at_split));
    $target_host = shift(@target_string_colon_split);
    $target_port = shift(@target_string_colon_split);
    $target_user = $self->param('user_w');
    $target_pass = $self->param('pass_w');
  }

  unless($target_dbname && $target_host && $target_user && $target_pass) {
    throw("The connection info for the database could not be recovered from the parameters passed in. ".
          "Make sure you pass in user_w, pass_w and either a db string or a hash to target_db");
  }

  my $taxonomy_db = $self->param('taxonomy_db');
  unless($taxonomy_db) {
    throw("You have used the populate_production_tables flag but have not passed in the production db connection hash with ".
          "the production_db flag");
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
            " -taxondbname ".$taxonomy_dbname;
            # Note script has no flag for -taxonuser, maybe it should be added. Hardcoded to ensro at the moment

  if($taxon_id) {
    $cmd .= " -taxon_id ".$taxon_id;
  } else {
    $cmd .= " -name '".$taxon_name."'";
  }

  my $return = system($cmd);
  if($return) {
    throw("The load_taxonomy script returned a non-zero exit value. Commandline used:\n".$cmd);
  }

}

1;






