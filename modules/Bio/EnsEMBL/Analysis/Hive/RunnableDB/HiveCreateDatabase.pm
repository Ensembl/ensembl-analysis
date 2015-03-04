package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase;

use strict;
use warnings;
use feature 'say';

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

use Data::Dumper;

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

sub convert_hash_to_db_string {
  my ($self,$connection_info) = @_;

  unless(defined($connection_info->{'-dbname'}) && defined($connection_info->{'-host'})) {
    throw("You have passed in a hash as your db info however the hash is missing either the dbname or host key,".
          " both are required if a hash is being passed in");
  }

  my $port = $connection_info->{'port'};

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
          "Expected format:\nmy_db_name\@myserver:port or my_db_name\@myserver");
  }
}

1;
