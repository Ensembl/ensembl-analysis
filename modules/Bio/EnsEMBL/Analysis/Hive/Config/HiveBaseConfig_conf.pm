=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 DESCRIPTION

This base config should be used by all pipeline created for the ensembl-annotation pipelines

=head1 METHODS

default_options: returns the default options from HiveGeneric_conf and it adds pipeline_db,
  dna_db and use_tokens. The inheriting class needs to specify; pipe_db_name, pipe_db_server,
  port, user, password, reference_db_name, reference_db_server, user_r, dna_db_name, dna_db_server

lsf_resource_builder: returns the parameters string for LSF meadow_type

=cut

package Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf;

use strict;
use warnings;
use feature 'say';

use File::Spec::Functions qw(catdir catfile);
use parent ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

use Bio::EnsEMBL::ApiVersion qw/software_version/;

=head2 default_options

 Arg [1]    : None
 Description: It returns a hashref containing the default options for HiveGeneric_conf
                use_tokens => 0,
                drop_databases => 0, # ONLY USE THIS PARAMETER ON THE COMMAND LINE
                databases_to_delete => [], # example: ['blast_db', 'refine_db', 'rough_db'],
                password_r => undef,

                ensembl_release => $ENV{ENSEMBL_RELEASE},
                genebuilder_id => $ENV{GENEBUILDER_ID} || 0,
                email => $ENV{HIVE_EMAIL},
                enscode_root_dir => $ENV{ENSCODE},
                software_base_path => $ENV{LINUXBREW_HOME},
                binary_base => catdir($self->o('software_base_path'), 'bin'),

                # Usefull if you want to use one server for all your databases, not great but ok
                data_db_server => $self->o('host'),
                data_db_host => $self->o('data_db_server'), # data_dbs_server will be deprecated
                data_db_port => $self->o('port'),
                data_db_user => $self->o('user'),
                data_db_password => $self->o('password'),
                data_db_driver => $self->o('hive_driver'),
                data_db_pass => $self->o('data_db_password'), # data_dbs_password will be deprecated

                dna_db_server => $self->o('host'),
                dna_db_host => $self->o('dna_db_server'), # dna_db_server will be deprecated
                dna_db_name => undef,
                dna_db_port => $self->o('port'),
                dna_db_user => $self->o('user_r'),
                dna_db_password => $self->o('password_r'),
                dna_db_pass => $self->o('dna_db_password'), # dna_db_password will be deprecated
                dna_db_driver => $self->o('hive_driver'),
                dna_db_host => $self->o('dna_db_host'),

                pipe_dbname => $self->o('dbowner').'_'.$self->o('pipeline_name').'_pipe',
                pipe_db_name => $self->o('pipe_dbname'), # pipe_dbname will be deprecated
                pipe_db_server => $self->o('host'),
                pipe_db_host => $self->o('pipe_db_server'), # pipe_db_server will be deprecated
                pipe_db_port => $self->o('port'),
                pipe_db_user => $self->o('user'),
                pipe_db_password => $self->o('password'),
                pipe_db_pass => $self->o('pipe_db_password'), # pipe_db_password will be deprecated
                pipe_db_driver => $self->o('hive_driver'),
              and two DB connection hash: pipeline_db and dna_db
 Returntype : Hashref
 Exceptions : None

=cut

sub default_options {
    my ($self) = @_;
    return {
        # inherit other stuff from the base class
        %{ $self->SUPER::default_options() },

        use_tokens => 0,
        drop_databases => 0, # ONLY USE THIS PARAMETER ON THE COMMAND LINE
        databases_to_delete => [], # example: ['blast_db', 'refine_db', 'rough_db'],
        password_r => undef,

        ensembl_release => $ENV{ENSEMBL_RELEASE},
        genebuilder_id => $ENV{GENEBUILDER_ID} || 0,
        email_address => $ENV{HIVE_EMAIL},
        enscode_root_dir => $ENV{ENSCODE},
        software_base_path => $ENV{LINUXBREW_HOME},
        binary_base => catdir($self->o('software_base_path'), 'bin'),

        guihive_host => 'http://guihive.ebi.ac.uk',
        guihive_port => 8080,

        data_db_server => $self->o('host'),
        data_db_host => $self->o('data_db_server'),
        data_db_port => $self->o('port'),
        data_db_user => $self->o('user'),
        data_db_password => $self->o('password'),
        data_db_driver => $self->o('hive_driver'),
        data_db_pass => $self->o('data_db_password'),

        dna_db_server => $self->o('host'),
        dna_db_host => $self->o('dna_db_server'),
        dna_db_name => undef,
        dna_db_port => $self->o('port'),
        dna_db_user => $self->o('user_r'),
        dna_db_password => $self->o('password_r'),
        dna_db_pass => $self->o('dna_db_password'),
        dna_db_driver => $self->o('hive_driver'),

        pipe_dbname => $self->o('dbowner').'_'.$self->o('pipeline_name').'_pipe',
        pipe_db_name => $self->o('pipe_dbname'),
        pipe_db_server => $self->o('host'),
        pipe_db_host => $self->o('pipe_db_server'),
        pipe_db_port => $self->o('port'),
        pipe_db_user => $self->o('user'),
        pipe_db_password => $self->o('password'),
        pipe_db_pass => $self->o('pipe_db_password'),
        pipe_db_driver => $self->o('hive_driver'),

        killlist_db_name => 'gb_kill_list',

        'pipeline_db' => {
            -dbname => $self->o('pipe_db_name'),
            -host   => $self->o('pipe_db_host'),
            -port   => $self->o('pipe_db_port'),
            -user   => $self->o('pipe_db_user'),
            -pass   => $self->o('pipe_db_pass'),
            -driver => $self->o('pipe_db_driver'),
        },

        'dna_db' => {
            -dbname => $self->o('dna_db_name'),
            -host   => $self->o('dna_db_host'),
            -port   => $self->o('dna_db_port'),
            -user   => $self->o('dna_db_user'),
            -pass   => $self->o('dna_db_pass'),
            -driver => $self->o('dna_db_driver'),
        },
    };
}


=head2 pipeline_create_commands

 Arg [1]    : None
 Description: If drop_databases is set to 1, it will delete all databases
              in databases_to_delete. -driver has to be set in the database
              you want to delete
 Returntype : Arrayref
 Exceptions : None

=cut

sub pipeline_create_commands {
  my ($self) = @_;

  my $drop_commands = $self->SUPER::pipeline_create_commands;
# To be able to have drop_databases as a commandline option we have to trick Hive so
# we store parameters and doing a string comparison on a number otherwise it fails...
  my $drop_databases = $self->o('drop_databases');
  my $databases_to_delete = $self->o('databases_to_delete');
  if ($drop_databases eq 1) {
    foreach my $database (@{$self->o('databases_to_delete')}) {
      push(@{$drop_commands}, $self->db_cmd('DROP DATABASE IF EXISTS', $self->dbconn_2_url($database, 1)));
    }
  }
  return $drop_commands;
}


=head2 hive_data_table

 Arg [1]    : String type, the type of data table you need: protein, cdna, refseq
 Arg [2]    : String name, the name of the table
 Description: Creates a table for protein or cdnas which can be used later in the pipeline
 Returntype : String mysql_query
 Exceptions : None

=cut

sub hive_data_table {
  my ($self, $type, $table_name) = @_;

  my %table_types = (
      protein => 'CREATE TABLE '.$table_name.' ('.
                  'accession varchar(50) NOT NULL,'.
                  'source_db varchar(50) NOT NULL,'.
                  'pe_level varchar(50) NOT NULL,'.
                  'biotype varchar(255) NOT NULL,'.
                  'group_name varchar(255) NOT NULL,'.
                  'seq text NOT NULL,'.
                  'PRIMARY KEY (accession))',
      refseq =>  'CREATE TABLE '.$table_name.' ('.
                  'accession varchar(50) NOT NULL,'.
                  'source_db varchar(50) NOT NULL,'.
                  'biotype varchar(25) NOT NULL,'.
                  'date varchar(50) NOT NULL,'.
                  'seq text NOT NULL,'.
                  'PRIMARY KEY (accession))',
      cdna =>    'CREATE TABLE '.$table_name.' ('.
                  'accession VARCHAR(50) NOT NULL,'.
                  'source VARCHAR(50) NOT NULL,'.
                  'date DATE DEFAULT 0,'.
                  'db_version VARCHAR(50) NOT NULL,'.
                  'species VARCHAR(50) NOT NULL,'.
                  'seq TEXT NOT NULL,'.
                  'protein_accession VARCHAR(50),'.
                  'PRIMARY KEY (accession))',
  );

  return $self->db_cmd($table_types{$type});
}


=head2 lsf_resource_builder

 Arg [1]    : String $queue, name of the queue to submit to, default is 'normal'
 Arg [2]    : Integer $mem, memory required in MB
 Arg [3]    : Arrayref String, list of server you will connect to during the analysis
 Arg [4]    : Arrayref Integer, list of tokens value for each server, default is 10
 Arg [5]    : Integer $num_threads, number of cores, use only if you ask for multiple cores
 Arg [6]    : String $extra_requirements, any other parameters you want to give to LSF option -R
 Arg [7]    : Arrayref String, any parameters related to your file system if you need to use -R"select[gpfs]"
 Example    : '1GB' => { LSF => $self->lsf_resource_builder('normal', 1000, [$self->default_options->{'pipe_db_server'}])},
              '3GB_multithread' => { LSF => $self->lsf_resource_builder('long', 3000, [$self->default_options->{'pipe_db_server'}], undef, 3)},
 Description: It will return the LSF requirement parameters you require based on the queue, the memory, the database servers, the number
              of CPUs. It uses options -q, -n, -M and -R. If you need any other other options you will need to add it to the returned string.
              If you want multiple cores from multiple CPUs, you will have to edit the returned string.
              If a server appears more than once, it will use the highest token value
              The command below will create a resource string for a job in the long queue, asking for 3GB of memory and it will reserve 10 token
              for server1, 20 for server2 and 10 for server3.
              $self->lsf_resource_builder('long', 3000, ['server1', 'server2', 'server3'], [, 20])
 Returntype : String
 Exceptions : None


=cut

sub lsf_resource_builder {
    my ($self, $queue, $memory, $servers, $tokens, $threads, $extra_requirements, $paths) = @_;

    my $lsf_requirement = '-q '.($queue || 'production-rh7');
    my @lsf_rusage;
    my @lsf_select;
    $extra_requirements = '' unless (defined $extra_requirements);
    if (defined $memory) {
        $lsf_requirement .= ' -M '.$memory;
        push(@lsf_rusage, 'mem='.$memory);
        push(@lsf_select, 'mem>'.$memory);
    }
    if (($self->o('use_tokens') eq 1) and defined $servers) {
        my $i = 0;
        my %seen;
        foreach my $server (@$servers) {
            if (! exists $seen{$server}) {
                my ($server_id) = $server =~ /(\d+)$/;
                push(@lsf_rusage, 'myens_build'.$server_id.'tok='.($tokens->[$i] || 10));
                $seen{$server} = $i;
            }
            else {
                my ($token) = $lsf_rusage[$seen{$server}] =~ /tok=(\d+)/;
                if (defined $token and ($tokens->[$i] || 10) > $token) {
                    $token = $tokens->[$i];
                    $lsf_rusage[$seen{$server}] =~ s/tok=\d+/tok=$token/;
                }
            }
            $i++;
        }
    }
    if (defined $threads) {
        $lsf_requirement .= ' -n '.$threads;
        $extra_requirements .= ' span[hosts=1]';
    }
    if (defined $paths) {
      foreach my $path (@$paths) {
        if ($path =~ '/gpfs/') {
          push(@lsf_select, 'gpfs');
        }
      }
    }
    return $lsf_requirement.' -R"select['.join(', ', @lsf_select).'] rusage['.join(', ', @lsf_rusage).'] '.$extra_requirements.'"';
}


=head2 slurm_resource_builder

 Arg [1]    : String $queue, name of the queue to submit to, default is 'standard'
 Arg [2]    : Integer $mem, memory required in MB
 Arg [3]    : String $time, allocated time for the job, format is 'mm', 'hh:mm:ss', 'dd-hh'
 Arg [4]    : Integer $num_threads, number of cores, use only if you ask for multiple cores
 Arg [5]    : String $extra_requirements, any other parameters you want to give to SLURM
 Example    : '1GB' => { SLURM => $self->slurm_resource_builder('standard', 1000, 5)},
              '3GB_multithread' => { SLURM => $self->slurm_resource_builder('long', 3000, '12:00:00', 3)},
 Description: It will return the SLURM requirement parameters you require based on the queue, the memory, time limit and the number
              of CPUs. If you need any other other options you can add it with Arg[5].
 Returntype : String
 Exceptions : None


=cut

sub slurm_resource_builder {
    my ($self, $queue, $memory, $time, $threads, $extra_requirements) = @_;

    my $slurm_requirement = '--partition='.($queue || 'standard')
                          .' --mem='.($memory || 1000)
                          .' --time='.($time || '1:00:00');
    if ($threads) {
        $slurm_requirement .= " --cpus-per-task=$threads";
    }
    return $slurm_requirement.($extra_requirements || '');
}


=head2 create_database_hash

 Arg [1]    : String host
 Arg [2]    : Int port
 Arg [3]    : String user name
 Arg [4]    : String password
 Arg [5]    : String dbname
 Arg [6]    : String driver
 Description: Create a database hash based on parameters given
              or using get_server_port_lists to randomly assign
              the location of the database
 Returntype : Hashref db connection details
 Exceptions : None

=cut

sub create_database_hash {
  my ($self, $host, $port, $user, $password, $dbname, $driver) = @_;

  my ($server_list, $port_list) = $self->get_server_port_lists;
  my $random = int(rand(@$server_list));
  my %dbconn = (
   -host => $host || $server_list->[$random],
   -port => $port || $port_list->[$random],
   -user => $user || $self->o('user'),
   -pass => $password || $self->o('password'),
   -dbname => $dbname || lc($self->o('dbowner').'_'.$self->o('pipeline_name').'_pipe'),
   -driver => $driver || $self->o('hive_driver'),
 );
 $self->warning("Your dbname has upper case character, it might cause problems, ".$dbconn{'-dbname'})
   if ($dbconn{'-dbname'} =~ /[[:upper:]]/);
  return \%dbconn;
}


=head2 get_server_port_lists

 Arg [1]    : None
 Description: Create 2 listes, one for the hosts and another for the ports
              They have to be synchronised to be able to fetch the port of
              a server in the port list with the index of the server in the
              hosts list
 Returntype : Array of two array ref
 Exceptions : None

=cut

sub get_server_port_lists {
  my ($self) = @_;

  my @server_list;
  my @port_list;
  foreach my $i (1..7) {
    push(@server_list, 'mysql-ens-genebuild-prod-'.$i);
    push(@port_list, 4526+$i);
  }
  return \@server_list, \@port_list;
}


=head2 get_meta_db_information

 Arg [1]    : Hashref, DB connection details
 Arg [2]    : String, database name
 Description: Creates the connection details for a Hive pipeline database, the url and
              the HTML link tag (<a>) to be able to easily connect to the meta database
              from the job panel in GuiHive.
              Arg[1] and Arg[2] are mutually exclusive
              If Arg[1] is not provided, it will create the connection details based on 'user',
              'password', Arg[2] and randomly choose a server.
 Returntype : Arrayref, composed of 3 elements: an arrayref for the database conenction,
              a string representing the url and a string containing HTML code
 Exceptions : None

=cut

sub get_meta_db_information {
  my ($self, $db_conn, $dbname) = @_;

  my $url;
  my $guiurl;
  if ($self->_is_second_pass) {
    if (!$db_conn) {
      $db_conn = $self->create_database_hash(undef, undef, $self->o('user'), $self->o('password'), $dbname);
    }
    $self->root->{_tmp_db_conn} = $db_conn;
    $url = $self->dbconn_2_url('_tmp_db_conn', 1);
    $guiurl = sprintf("<a target='_blank' href='%s:%s/?driver=%s&username=%s&passwd=%s&host=%s&port=%d&dbname=%s'>guihive</a>", $self->o('guihive_host'), $self->o('guihive_port'), $db_conn->{'-driver'}, $db_conn->{'-user'}, $db_conn->{'-pass'}, $db_conn->{'-host'}, $db_conn->{'-port'}, $db_conn->{'-dbname'});
  }
  return $db_conn, $url, $guiurl;
}


=head2 _is_second_pass

 Arg [1]    : String (optional), a parameter name from the default_options hash
              Default parameter is 'pipe_db_name'
 Description: Check if all substitutions have been done and that $self->o() can
              be used safely
 Returntype : Boolean, true if substitutions have been done
 Exceptions : None

=cut

sub _is_second_pass {
  my ($self, $param_name) = @_;

  $param_name ||= 'pipe_db_name';

  return $self->o($param_name) !~ /:subst/;
}

1;
