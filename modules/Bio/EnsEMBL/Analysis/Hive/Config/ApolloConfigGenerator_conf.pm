package ApolloConfigGenerator_conf;

use strict;
use warnings;
use File::Spec::Functions;

use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_analysis_settings);
use base ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

sub default_options {
  my ($self) = @_;
  return {
    # inherit other stuff from the base class
    %{ $self->SUPER::default_options() },

# The database file should be a list of core dbs, with the core db server shortcut, followed the output server shortcut
# You can see the shortcuts in the actuall apollo script. The format of the file is:
## core_db_name@core_server_shortcut:output_db_shortcut
# prim96_pongo_abelii_core_96@gb6:gb5
# prim96_prolemur_simus_core_96@gb6:gb5
# prim96_theropithecus_gelada_core_96@gb6:gb5

'database_file_path'          => '', # The text file to generate the config from
'apollo_config_script'        => '', # The apollo config generator script
'output_dir'                  => '', # output dir for the config file
'sleep_timer'                 => 1800, # The time in seconds to sleep between generating new config

'user_r'                      => '',
'user_w'                      => '',
'password'                    => '',

'pipe_db_host'                => '',
'port'                        => '',

# Mostly constant stuff
'pipeline_name'               => 'apollo_config',
'default_mem'                 => '900',

'user'                        => $self->o('user_w'),
'dna_db_name'                 => '', # Leave blank, just needs to be present

 } # end return
} # end default_options


sub pipeline_create_commands {
  my ($self) = @_;
  return [
    @{$self->SUPER::pipeline_create_commands},
  ];
} # end pipeline_create_commands



sub pipeline_analyses {
  my ($self) = @_;
  return [

    {
      -logic_name => 'create_apollo_config',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
                      cmd => 'perl '.$self->o('apollo_config_script').
                             ' -output_dir '.$self->o('output_dir').
                             ' -database_file_path '.$self->o('database_file_path'),
                     },
      -rc_name    => 'default',
      -max_retry_count => 0,
      -flow_into  => {
                       '1' => ['sleep_cycle'],
                     },
      -input_ids  => [{}],
    },

    {
      -logic_name => 'sleep_cycle',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
                       cmd => 'sleep '.$self->o('sleep_timer'),
                     },
      -rc_name    => 'default',
      -max_retry_count => 0,
      -flow_into  => {
                       '1' => ['create_apollo_config'],
                     },
    },

  ]
} # end pipeline analyses

sub resource_classes {
  my $self = shift;

  return {
    'default'       => { 'LSF' => $self->lsf_resource_builder('production-rh7',$self->default_options->{'default_mem'}) },
  }
} # end resource_classes

1;
