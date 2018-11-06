=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase;

use strict;
use warnings;
use feature 'say';

use File::Basename;
use File::Path qw(make_path);
use File::Spec::Functions qw(catfile);

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name execute_with_wait);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

=head2 param_defaults

 Arg [1]    : None
 Description: Default parameters
               ignore_dna => 0, # if set to 1, the dna table won't be dumped
               force_drop => 0,
               _lock_tables => 'true', # can be true or false but should not be changed unless you have problems
               tables_to_reset => [
                 'analysis',
                 'gene',
                 'transcript',
                 'translation',
                 'exon',
                 'dna_align_feature',
                 'protein_align_feature',
                 'intron_supporting_evidence',
               ],
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    ignore_dna => 0, # if set to 1, the dna table won't be dumped
    _lock_tables => 'true', # can be true or false but should not be changed unless you have problems
    force_drop => 0,
    tables_to_reset => [
      'analysis',
      'gene',
      'transcript',
      'translation',
      'exon',
      'dna_align_feature',
      'protein_align_feature',
      'intron_supporting_evidence',
    ],
    vital_tables => vital_tables(),
    use_bash_pipefail => 1,
    use_bash_errexit => 1,
  }
}

=head2 fetch_input

 Arg [1]    : None
 Description: Check that 'create_type' has been specified. These are the type allowed:
                clone
                dna_db
                copy
                core_only
                backup DEPRECATED, please use Bio::EnsEMBL::Hive::RunnableDB::DatabaseDumper with "src_db_conn" and "output_file"
              If 'db_dump_file' is specified, creates the directory
              If 'output_path' is specified and does not exist, creates the directory
 Returntype : None
 Exceptions : Throws if 'create_type' is not specified

=cut

sub fetch_input {
  my $self = shift;

  if ($self->param_is_defined('skip_analysis') and $self->param('skip_analysis')) {
    $self->complete_early('You asked to sip the analysis');
  }
  my $create_type = $self->param_required('create_type');
  if ($self->param_is_defined('db_dump_file')) {
    my $dump_dir = dirname($self->param('db_dump_file'));
    make_path($dump_dir);
  }
  if ($self->param_is_defined('output_path')) {
    if (! -e $self->param('output_path')) {
      make_path($self->param('output_path'));
    }
  }
  my @commands;
  if($self->param('force_drop')) {
    push(@commands, $self->drop_database($self->param('target_db')));
  }

  if ($create_type eq 'clone' or $create_type eq 'dna_db') {
    push(@commands, $self->create_database($self->param_required('target_db')));
    push(@commands, @{$self->clone_database($self->param_required('source_db'), $self->param_required('target_db'))});
    my $cmd = $self->reset_autoincrement($self->param('target_db'));
    push(@commands, $cmd) if ($cmd);
    if($create_type eq 'dna_db') {
      if (not $self->param_is_defined('db_dump_file')) {
        $self->param('db_dump_file', create_file_name('source_db'));
      }
      push(@commands,
        $self->dump_database($self->param('source_db'), undef, 0, 0, ['dna', 'repeat_feature', 'repeat_consensus']).' | '.
        $self->load_database($self->param('target_db')));
    }
  }
  elsif ($create_type eq 'copy') {
    if (not $self->param_is_defined('db_dump_file')) {
      $self->param('db_dump_file', create_file_name('source_db'));
    }
    push(@commands, $self->create_database($self->param_required('target_db')));
    push(@commands,
      $self->dump_database($self->param('source_db'), undef, $self->param('ignore_dna')).' | '.
      $self->load_database($self->param('target_db')));
    my $cmd = $self->reset_autoincrement($self->param('target_db'));
    push(@commands, $cmd) if ($cmd);
  }
  elsif($create_type eq 'core_only') {
    my $table_file = catfile($self->param_required('enscode_root_dir'), 'ensembl', 'sql', 'table.sql');
    $self->throw("You have specified a create type of core_only but the path from enscode_root_dir to the table file is incorrect:\n$table_file")
      unless(-e $table_file);
    push(@commands, $self->create_database($self->param_required('target_db')));
    push(@commands, $self->load_database($self->param('target_db'), $table_file));
  }
  elsif($create_type eq 'backup') {
    $self->warning('Create type "backup" is deprecated, please use Bio::EnsEMBL::Hive::RunnableDB::DatabaseDumper instead with "src_db_conn" and "output_file"');
    my $dump_file = catfile($self->param_required('output_path'), $self->param_required('backup_name'));
    push(@commands, $self->dump_database($self->param_required('source_db'), $dump_file, $self->param('ignore_dna'), 1));
  }
  else {
    $self->throw("You have specified a create type of ".$create_type.", however this is not supported by the module");
  }
  $self->param('commands', \@commands);
}

=head2 run

 Arg [1]    : None
 Description: Run all the command generated in fetch_input
 Returntype : None
 Exceptions : Throws if one of the command fails

=cut

sub run {
  my $self = shift;

  my $command_options = {
    use_bash_pipefail => $self->param('use_bash_pipefail'),
    use_bash_errexit => $self->param('use_bash_errexit'),
  };
  foreach my $cmd (@{$self->param('commands')}) {
    $self->say_with_header($cmd);
    if ($self->run_system_command($cmd, $command_options)) {
      $self->throw("Could not execute '$cmd'");
    }
  }
}


=head2 write_output

 Arg [1]    : None
 Description: If 'db_dump_file' is defined, delete the file if it still exists
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my $self = shift;

  if ($self->param_is_defined('db_dump_file') and -e $self->param_is_defined('db_dump_file')) {
    unlink $self->param_is_defined('db_dump_file');
  }
}


=head2 clone_database

 Arg [1]    : Hashref source db connection details
 Arg [2]    : Hashref target db connection details
 Description: Creates the commands to:
               dump the schema from Arg[1] and load it into Arg[2]
               dump the data from the vital tables in 'vital_tables' from Arg[1] and load it into Arg[2]
 Returntype : Arrayref of String
 Exceptions : Throw if 'extra_data_tables' is defined but not an arrayref

=cut

sub clone_database {
  my ($self, $source_db, $target_db) = @_;

  my $do_lock = $self->param('_lock_tables');
  my $vital_tables = $self->param('vital_tables');
  if ($self->param_is_defined('extra_data_tables') and $self->param('extra_data_tables')) {
    $self->throw('You should give an array ref for your extra tables to dump') unless (ref($self->param('extra_data_tables')) eq 'ARRAY');
    push(@$vital_tables, @{$self->param('extra_data_tables')});
  }
  my ($source_dbhost, $source_dbport, $source_dbname, $source_dbuser, $source_dbpass) = $self->db_connection_details($source_db, 1);
  my ($target_dbhost, $target_dbport, $target_dbname, $target_dbuser, $target_dbpass) = $self->db_connection_details($target_db);
  my $base_source_cmd = "mysqldump --lock-tables=$do_lock -h $source_dbhost -P $source_dbport -u $source_dbuser";
  $base_source_cmd .= ' -p'.$source_dbpass if ($source_dbpass);
  my $base_target_cmd = "mysql -h $target_dbhost -P $target_dbport -u $target_dbuser -D $target_dbname";
  $base_target_cmd .= ' -p'.$target_dbpass if ($target_dbpass);
  my @commands = (
    "$base_source_cmd --no-data $source_dbname | $base_target_cmd",
    join(' ', $base_source_cmd, $source_dbname, @$vital_tables, '|',$base_target_cmd),
  );

  return \@commands;
}


=head2 dump_database

 Arg [1]    : Hashref connection details
 Arg [1]    : String filename
 Arg [1]    : Boolean true to ignore the dna table
 Arg [1]    : Boolean true to compress the file using gzip
 Arg [1]    : Arrayref a list of tables to dump
 Description: Create the command to dump a database. It will also check
              the value for max_allowed_packet to set it corrrectly if
              needed
 Returntype : String
 Exceptions : Throws if it cannot get the max_allowed_packet variable

=cut

sub dump_database {

  my ($self, $db, $db_file, $ignore_dna, $compress, $tables) = @_;

  my ($dbhost, $dbport, $dbname, $dbuser, $dbpass) = $self->db_connection_details($db, 1);

  my $command = "mysqldump -h$dbhost -P$dbport -u$dbuser ";
  if ($dbpass) { # dbpass for read access can be optional
    $command .= " -p$dbpass";
  }
  if ($ignore_dna) {
    $command .= " --ignore-table=".$dbname.".dna ";
  }
  else {
    # Check the max_allowed_packet before dumping tables
    my $max_allowed_packet;
    my $checkcmd = $command;
    $checkcmd =~ s/dump/admin/;
    open(RH, $checkcmd.' variables |') || $self->throw("Could not execute command: $checkcmd variables");
    while(<RH>) {
      if (/max_allowed_packet\D+(\d+)/) {
        $max_allowed_packet = $1;
        last;
      }
    }
    close(RH) || $self->throw("Could not close: $checkcmd variables");
    if ($max_allowed_packet) {
      $command .= ' --max_allowed_packet '.$max_allowed_packet;
    }
  }
  $command .= " $dbname";
  if ($tables and @$tables) {
    $command .= ' '.join(' ', @$tables);
  }
  if ($db_file) {
    if ($compress) {
      # If pipefail is not set your command can fail without telling you
      $command = "$command | gzip > $db_file.gz";
    }
    else {
      $command .= " > $db_file";
    }
  }
  return $command;
}


=head2 create_database

 Arg [1]    : Hashref connection details
 Description: Create the command to create a new database
 Returntype : String
 Exceptions : None

=cut

sub create_database {
  my ($self, $db) = @_;

  my ($dbhost, $dbport, $dbname, $dbuser, $dbpass) = $self->db_connection_details($db);
  return "mysql -h$dbhost -P$dbport -u$dbuser".($dbpass ? " -p$dbpass " : " ")."-e'CREATE DATABASE $dbname'";
}


=head2 load_database

 Arg [1]    : Hashref connection detailes
 Arg [2]    : String filename
 Description: Create the command to load data into the database
 Returntype : String
 Exceptions : None

=cut

sub load_database {
  my ($self, $db, $db_file) = @_;

  my $base_cmd = $self->get_client_cmd($self->db_connection_details($db));
  if ($db_file) {
    return "$base_cmd < $db_file";
  }
  else {
    return $base_cmd;
  }
}


=head2 drop_database

 Arg [1]    : Hashref connection details
 Description: Create the command to drop the database given in Arg[1]
 Returntype : String
 Exceptions : None

=cut

sub drop_database {
  my ($self, $db) = @_;

  my ($dbhost, $dbport, $dbname, $dbuser, $dbpass) = $self->db_connection_details($db);
  $self->warning("Dropping existing database if it exists $dbname"."@"."$dbhost:$dbport...");
  return "mysql -h$dbhost -P$dbport -u$dbuser".($dbpass ? " -p$dbpass " : " ")."-e 'DROP DATABASE IF EXISTS $dbname'";
}


=head2 remove_file

 Arg [1]    : String filename
 Description: Remove the dump file if it still exists
 Returntype : None
 Exceptions : None

=cut

sub remove_file {
  my ($self, $db_file) = @_;

  if (-e $db_file) {
    print "Deleting file $db_file\n";
    unlink $db_file;
  }
}


=head2 reset_autoincrement

 Arg [1]    : Hashref connection details
 Description: Command to eset the auto increment on a series of tables with 'tables_to_reset'.
              The default array ref should be enough
 Returntype : String
 Exceptions : None

=cut

sub reset_autoincrement {
  my ($self, $db) = @_;

  if (@{$self->param('tables_to_reset')}) {
    my $query;
    foreach my $table (@{$self->param('tables_to_reset')}) {
      $query .= "ALTER TABLE $table AUTO_INCREMENT = 1;";
    }
    my $base_cmd = $self->get_client_cmd($self->db_connection_details($db));
    return "$base_cmd -e '$query'";
  }
  else {
    $self->warning('You have no tables you want to reset AUTO_INCREMENT, strange...');
    return;
  }
}


=head2 db_connection_details

 Arg [1]    : Hashref
 Arg [2]    : Boolean true if read only user
 Description: Getter for the connection details of a database.
              If Arg[2] is true, dbuser will be 'user_r' if defined and dbpass
              will be undef unless 'pass_r' is defined.
              Otherwise dbuser is 'user_w' if defined and dbpass is 'pass_w' if
              defined.
 Returntype : Array dbhost, dbport, dbname, dbuser, dbpass, driver
 Exceptions : None

=cut

sub db_connection_details {
  my ($self, $db, $read_only) = @_;

  my $dbhost = $db->{'-host'};
  my $dbport = $db->{'-port'};
  my $dbuser = $db->{'-user'};
  my $dbpass = $db->{'-pass'};
  my $dbname = $db->{'-dbname'};
  if ($read_only) {
    if ($self->param_is_defined('user_r')) {
      $dbuser = $self->param('user_r');
      $dbpass = undef;
    }
    if ($self->param_is_defined('pass_r')) {
      $dbpass = $self->param('pass_r');
    }
  }
  else {
    if ($self->param_is_defined('user_w')) {
      $dbuser = $self->param('user_w');
    }
    if ($self->param_is_defined('pass_w')) {
      $dbpass = $self->param('pass_w');
    }
  }
  return $dbhost, $dbport, $dbname, $dbuser, $dbpass, ((exists $db->{'-driver'} && $db->{'-driver'}) ? $db->{'-driver'} : 'mysql');
}


=head2 get_client_cmd

 Arg [1]    : String host
 Arg [2]    : Int port
 Arg [3]    : String dbname
 Arg [4]    : String user
 Arg [5]    : String password
 Arg [6]    : String driver
 Description: Create the CLI command to be fed for create the real command
              The aim is to be able to get the correct command based on
              the driver
 Returntype : String
 Exceptions : None

=cut

sub get_client_cmd {
  my ($self, $dbhost, $dbport, $dbname, $dbuser, $dbpass, $driver) = @_;

  my $binary = 'mysql';
  if ($driver eq 'postgres') {
    return "psql -w -h $dbhost -p $dbport -U $dbuser -D $dbname"
  }
  else {
    return "mysql -h$dbhost -P$dbport -u$dbuser".($dbpass ? " -p$dbpass " : " ")."-D $dbname";
  }
}


=head2 vital_tables

 Arg [1]    : None
 Description: Getter of the tables that need to be synchronised between databases
              It can be used as a standalone method
              At the moment the tables are:
                analysis
                analysis_description
                assembly
                assembly_exception
                attrib_type
                coord_system
                external_db
                meta
                misc_set
                seq_region
                seq_region_attrib
                seq_region_synonym
                karyotype
                mapping_set
                seq_region_mapping
 Returntype : Arrayref
 Exceptions : None

=cut

sub vital_tables {
  return [
    'analysis',
    'analysis_description',
    'assembly',
    'assembly_exception',
    'attrib_type',
    'coord_system',
    'external_db',
    'meta',
    'misc_set',
    'seq_region',
    'seq_region_attrib',
    'seq_region_synonym',
    'karyotype',
    'mapping_set',
    'seq_region_mapping',
  ];
}

1;
