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

=cut

package Hive_exonerate_protein_file_conf;

use strict;
use warnings;
use feature 'say';

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

use Bio::EnsEMBL::ApiVersion qw/software_version/;

sub default_options {
    my ($self) = @_;
    return {
        # inherit other stuff from the base class
        %{ $self->SUPER::default_options() },

############################################################################
#
# MOSTLY CHANGE
#
############################################################################

        'pipeline_name'              => '',
        'pipe_db_name'                => '',
        'pipe_db_server'             => '',

        'reference_db_name'           => '', # This should be a core db with DNA
        'reference_db_server'        => '',

        'exonerate_output_db_name'    => '', # This will be created for you automatically
        'exonerate_output_db_server' => '',

        'user_w'                     => '',
        'password'                   => '',
        'user_r'                     => 'ensro',
        'port'                       => 3306,

        'exonerate_pid'              => '50', # Cut-off for percent id
        'exonerate_cov'              => '50', # Cut-off for coverage

        'output_path'                => '',
        'protein_file'               => '',
        'genome_file'                => '',

        'repeat_masking_logic_names' => [], # e.g ['repeatmask_repbase_baboon','dust']


############################################################################
#
# MOSTLY CONSTANT
#
############################################################################

        'create_type'                          => 'clone',
        'num_tokens'                           => '10',
        'user'                                 => 'ensro',
        'driver'                               => 'mysql',
        'protein_table_name'                   => 'protein_sequences',

        'exonerate_path'                       => '/software/ensembl/genebuild/usrlocalensemblbin/exonerate-0.9.0',
        'exonerate_calculate_coverage_and_pid' => '1',

        'default_mem'                => '900',
        'exonerate_mem'              => '2900',
        'exonerate_retry_mem'        => '5900',

        'pipeline_db' => {
            -host   => $self->o('pipe_db_server'),
            -port   => $self->o('port'),
            -user   => $self->o('user_w'),
            -pass   => $self->o('password'),
            -dbname => $self->o('pipe_db_name'),
            -driver => $self->o('driver'),
        },

        'reference_db' => {
                            -dbname => $self->o('reference_db_name'),
                            -host   => $self->o('reference_db_server'),
                            -port   => $self->o('port'),
                            -user   => $self->o('user_r'),

                          },

        'dna_db' => {
                      -dbname => $self->o('reference_db_name'),
                      -host   => $self->o('reference_db_server'),
                      -port   => $self->o('port'),
                      -user   => $self->o('user_r'),
                    },

        'exonerate_output_db' => {
                                   -dbname => $self->o('exonerate_output_db_name'),
                                   -host   => $self->o('exonerate_output_db_server'),
                                   -port   => $self->o('port'),
                                   -user   => $self->o('user_w'),
                                   -pass   => $self->o('password'),
                                 },

    };
}

sub pipeline_create_commands {
    my ($self) = @_;
    return [
    # inheriting database and hive tables' creation
      @{$self->SUPER::pipeline_create_commands},

      $self->db_cmd('CREATE TABLE '.$self->o('protein_table_name').' ('.
                    'accession varchar(50) NOT NULL,'.
                    'source varchar(50) NOT NULL,'.
                    'date varchar(50) NOT NULL,'.
                    'db_version varchar(50) NOT NULL,'.
                    'species varchar(50) NOT NULL,'.
                    'seq text NOT NULL,'.
                    'PRIMARY KEY (accession))'),
    ];
}

sub hive_meta_table {
    my ($self) = @_;
    return {
            %{$self->SUPER::hive_meta_table},       # here we inherit anything from the base class

        'hive_use_param_stack'  => 1,           # switch on the new param_stack mechanism
    };
}


## See diagram for pipeline structure
sub pipeline_analyses {
    my ($self) = @_;

    return [

      {
        -logic_name => 'create_exonerate_output_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('reference_db'),
                         target_db => $self->o('exonerate_output_db'),
                         create_type => $self->o('create_type'),
                       },
        -rc_name    => 'default',
        -input_ids => [{}],
      },

      {
        -logic_name => 'load_protein_file',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadProteins',
        -parameters => {
                         protein_table_name => $self->o('protein_table_name'),
                       },
        -flow_into => {
                        1 => ['exonerate'],
                      },
        -input_ids => [{'protein_file' => $self->o('protein_file')},
                      ],
        -rc_name    => 'default',
      },

     {
        -logic_name => 'exonerate',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
        -rc_name    => 'exonerate',
        -wait_for   => ['create_exonerate_output_db'],
        -parameters => {
                         iid_type => 'db_seq',
                         query_table_name => 'protein_sequences',
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('exonerate_output_db'),
                         logic_name => 'exonerate',
                         module     => 'HiveExonerate2Genes',
                         config_settings => $self->get_config_settings('exonerate_protein','exonerate'),
                         query_seq_dir => $self->o('output_path'),
                         calculate_coverage_and_pid => $self->o('exonerate_calculate_coverage_and_pid'),
                      },
        -flow_into => {
                        -1 => ['exonerate_retry'],
                        -2 => ['failed_exonerate_proteins'],
                      },
        -batch_size => 100,
        -failed_job_tolerance => 5,
     },

     {
        -logic_name => 'exonerate_retry',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
        -rc_name    => 'exonerate',
        -parameters => {
                         iid_type => 'db_seq',
                         query_table_name => 'protein_sequences',
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('exonerate_output_db'),
                         logic_name => 'exonerate',
                         module     => 'HiveExonerate2Genes',
                         config_settings => $self->get_config_settings('exonerate_protein','exonerate'),
                         query_seq_dir => $self->o('output_path'),
                         calculate_coverage_and_pid => $self->o('exonerate_calculate_coverage_and_pid'),
                      },
        -flow_into => {
                        -2 => ['failed_exonerate_proteins'],
                      },
        -batch_size => 100,
        -failed_job_tolerance => 100,
        -can_be_empty => 1,
      },


      {
        -logic_name => 'failed_exonerate_proteins',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {},
        -rc_name     => 'default',
        -can_be_empty  => 1,
        -failed_job_tolerance => 100,
      },

    ];
}

sub pipeline_wide_parameters {
    my ($self) = @_;

    return {
        %{ $self->SUPER::pipeline_wide_parameters() },  # inherit other stuff from the base class
    };
}

# override the default method, to force an automatic loading of the registry in all workers
#sub beekeeper_extra_cmdline_options {
#    my $self = shift;
#    return "-reg_conf ".$self->o("registry");
#}

sub resource_classes {
    my $self = shift;


  my $pipe_db_server = $self->default_options()->{'pipe_db_server'};
  my $dna_db_server = $self->default_options()->{'reference_db_server'};
  my $exonerate_output_db_server = $self->default_options()->{'exonerate_output_db_server'};

  my $pipe_db_server_number;
  my $dna_db_server_number;
  my $exonerate_output_db_server_number;

  my $default_mem = $self->default_options()->{'default_mem'};
  my $exonerate_mem = $self->default_options()->{'exonerate_mem'};
  my $exonerate_retry_mem = $self->default_options()->{'exonerate_retry_mem'};


  my $num_tokens = $self->default_options()->{'num_tokens'};

    unless($pipe_db_server =~ /(\d+)$/) {
    die "Failed to parse the server number out of the pipeline db server name. This is needed for setting tokens\n".
        "pipe_db_server: ".$pipe_db_server;
  }

  $pipe_db_server_number = $1;

    unless($dna_db_server =~ /(\d+)$/) {
    die "Failed to parse the server number out of the pipeline db server name. This is needed for setting tokens\n".
        "dna_db_server: ".$dna_db_server;
  }

  $dna_db_server_number = $1;

  unless($exonerate_output_db_server =~ /(\d+)$/) {
    die "Failed to parse the server number out of the exonerate db server name. This is needed for setting tokens\n".
        "exonerate_output_db_server: ".$exonerate_output_db_server;
  }

  $exonerate_output_db_server_number = $1;

  unless($num_tokens) {
    die "num_tokens is uninitialised or zero. num_tokens needs to be present in default_options and not zero\n".
        "num_tokens: ".$num_tokens;
  }

    return {
      'default' => { LSF => '-q normal -M'.$default_mem.' -R"select[mem>'.$default_mem.'] '.
                            'rusage[mem='.$default_mem.','.
                            'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

      'exonerate' => { LSF => '-q normal -W 120 -M'.$exonerate_mem.' -R"select[mem>'.$exonerate_mem.'] '.
                              'rusage[mem='.$exonerate_mem.','.
                              'myens_build'.$exonerate_output_db_server_number.'tok='.$num_tokens.','.
                              'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
                              'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

      'exonerate_retry' => { LSF => '-q normal -W 120 -M'.$exonerate_retry_mem.' -R"select[mem>'.$exonerate_retry_mem.'] '.
                                    'rusage[mem='.$exonerate_retry_mem.','.
                                    'myens_build'.$exonerate_output_db_server_number.'tok='.$num_tokens.','.
                                    'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
                                    'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },
    }
}

sub get_config_settings {

   # Shift in the group name (a hash that has a collection of logic name hashes and a default hash)
   # Shift in the logic name of the specific analysis
   my $self = shift;
   my $config_group = shift;
   my $config_logic_name = shift;

   # And additional hash keys will be stored in here
   my @additional_configs = @_;

   # Return a ref to the master hash for the group using the group name
   my $config_group_hash = $self->master_config_settings($config_group);
   unless(defined($config_group_hash)) {
     die "You have asked for a group name in master_config_settings that doesn't exist. Group name:\n".$config_group;
   }
   # Final hash is the hash reference that gets returned. It is important to note that the keys added have
   # priority based on the call to this subroutine, with priority from left to right. Keys assigned to
   # $config_logic_name will have most priority, then keys in any additional hashes, then keys from the
   # default hash. A default hash key will never override a $config_logic_name key
   my $final_hash;

   # Add keys from the logic name hash
   my $config_logic_name_hash = $config_group_hash->{$config_logic_name};
   unless(defined($config_logic_name_hash)) {
     die "You have asked for a logic name hash that doesn't exist in the group you specified.\n".
         "Group name:\n".$config_group."\nLogic name:\n".$config_logic_name;
   }

   $final_hash = $self->add_keys($config_logic_name_hash,$final_hash);

   # Add keys from any additional hashes passed in, keys that are already present will not be overriden
   foreach my $additional_hash (@additional_configs) {
     my $config_additional_hash = $config_group_hash->{$additional_hash};
     $final_hash = $self->add_keys($config_additional_hash,$final_hash);
   }

   # Default is always loaded and has the lowest key value priority
   my $config_default_hash = $config_group_hash->{'Default'};
   $final_hash = $self->add_keys($config_default_hash,$final_hash);

   return($final_hash);
}

sub add_keys {
  my ($self,$hash_to_add,$final_hash) = @_;

  foreach my $key (keys(%$hash_to_add)) {
    unless(exists($final_hash->{$key})) {
      $final_hash->{$key} = $hash_to_add->{$key};
    }
  }

  return($final_hash);
}

sub master_config_settings {

  my ($self,$config_group) = @_;
  my $master_config_settings = {

  exonerate_protein => {
    Default => {
                 IIDREGEXP           => '(\d+):(\d+)',
                 OPTIONS             => '--model protein2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE --bestn 1',
                 COVERAGE_BY_ALIGNED => 0,
                 QUERYTYPE           => 'protein',
                 GENOMICSEQS         => $self->o('genome_file'),
                 PROGRAM             => $self->o('exonerate_path'),
                 SOFT_MASKED_REPEATS => $self->o('repeat_masking_logic_names'),
               },

    exonerate => {
                   FILTER                        => {
                     OBJECT                      => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
                     PARAMETERS                  => {
                       -coverage                 => $self->o('exonerate_cov'),
                       -percent_id               => $self->o('exonerate_pid'),
                       -best_in_genome           => 1,
                       -reject_processed_pseudos => 1,
                     },
                   },
                 },

    },
  };

  return($master_config_settings->{$config_group});

}

1;
