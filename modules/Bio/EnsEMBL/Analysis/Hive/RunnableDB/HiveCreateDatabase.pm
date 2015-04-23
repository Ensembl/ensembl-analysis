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
          "\nExpected format:\nmy_db_name\@myserver:port or my_db_name\@myserver");
  }
}

sub dump_database {

  my ($dbhost,$dbport,$dbuser,$dbpass,$dbname,$db_file) = @_;

  print "\nDumping database $dbname"."@"."$dbhost:$dbport...\n";
  
  my $command;
  if (!$dbpass) { # dbpass for read access can be optional
  	$command = "mysqldump -h$dbhost -P$dbport -u$dbuser $dbname > $db_file";
  } else {
  	$command = "mysqldump -h$dbhost -P$dbport -u$dbuser -p$dbpass $dbname > $db_file";
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

sub DESTROY {
  my $self = shift;
  
  my $db_file = $self->param('db_dump_file');
  
  if (-e $db_file) {
  	print "Deleting temporary file $db_file\n";
  	if (system("rm -f $db_file")) {
  	  throw("Couldn't delete temporary file $db_file");
  	} else {
  	  print "Temporary file $db_file deleted\n";
  	}
  }	
}

1;






