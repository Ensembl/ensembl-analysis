=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

package GenomePreparation_conf;

use strict;
use warnings;
use feature 'say';
use Data::Dumper;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

use Bio::EnsEMBL::ApiVersion qw/software_version/;

sub default_options {
  my ($self) = @_;
  return {
    # inherit other stuff from the base class
    %{ $self->SUPER::default_options() },

'assembly_data_directory' => '',
'assembly_name'           => '',
'ensembl_api_dir'         => '',

# Repeatmasker parameters
'repeatmasker_path'       => '',
'repeatmasker_library'    => '',
'repeatmasker_engine'     => '',

# Genscan parameters
'genscan_path'            => '',
'genscan_matrix_path'     => '',

# Database settings
'pipe_dbname'             => '',
'pipeline_name'           => '',
'pipe_db_server'          => '',
'core_dbname'             => '',
'core_db_server'          => '',

'user'                    => '',
'user_r'                  => '',
'user_w'                  => '',
'pass_w'                  => '',
'port'                    => ,

'driver' => 'mysql',
'create_type' => 'clone',
'create_db_script_path' => $self->o('ensembl_api_dir').'/ensembl-analysis/scripts/clone_database.ksh',

# DB connection info
'pipeline_db' => {
                    # connection parameters
                   -host   => $self->o('pipe_db_server'),
                   -port   => $self->o('port'),
                   -user   => $self->o('user_w'),
                   -pass   => $self->o('pass_w'),
                   -dbname => $self->o('pipe_dbname'),
                   -driver => $self->o('driver'),
                 },

'core_db' => {
               -dbname => $self->o('core_dbname'),
               -host   => $self->o('core_db_server'),
               -port   => $self->o('port'),
               -user   => $self->o('user_w'),
               -pass   => $self->o('pass_w'),
               -driver => $self->o('driver'),
             },

'dna_db' => {
              -dbname => $self->o('core_dbname'),
              -host   => $self->o('core_db_server'),
              -port   => $self->o('port'),
              -user   => $self->o('user_r'),
            },

    };
}

sub pipeline_create_commands {
    my ($self) = @_;
    return [
    # inheriting database and hive tables' creation
      @{$self->SUPER::pipeline_create_commands},
    ];
}


## See diagram for pipeline structure
sub pipeline_analyses {
    my ($self) = @_;

    return [

            {
              -logic_name => 'create_emtpy_core_db',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
              -parameters => {
                               'target_db'        => $self->o('core_db'),
                               'user_w'           => $self->o('user_w'),
                               'pass_w'           => $self->o('pass_w'),
                               'enscode_root_dir' => $self->o('ensembl_api_dir'),
                               'create_type'      => 'core_only',
                             },
              -rc_name    => 'local',
              -flow_into => {
                              '1' => ['load_assembly_info'],
                            },
              -input_ids  => [ {} ],
            },

            {
              -logic_name => 'load_assembly_info',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadAssembly',
              -parameters => {
                               'simple_load' => 1,
                               'target_db' => $self->o('core_db'),
                               'assembly_data_directory' => $self->o('assembly_data_directory'),
                               'assembly_name' => $self->o('assembly_name'),
                               'enscode_root_dir' => $self->o('ensembl_api_dir'),
                             },
              -flow_into => {
                              '1' => ['create_200kb_slice_ids'],
                            },

              -rc_name    => 'local',

            },


            {
              # Create 200kb slices, each species flow into this independantly
              -logic_name => 'create_200kb_slice_ids',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
              -parameters => {
                               target_db => $self->o('core_db'),
                               coord_system_name => 'toplevel',
                               iid_type => 'slice',
                               slice_size => 200000,
                               include_non_reference => 0,
                               top_level => 1,
                             },
              -flow_into => {
                              '2' => ['run_repeatmasker'],
                            },

            },


            {
              -logic_name => 'run_repeatmasker',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveRepeatMasker',
              -parameters => {
                               target_db => $self->o('core_db'),
                               logic_name => 'repeatmasker_repbase_'.$self->o('repeatmasker_library'),
                               module => 'HiveRepeatMasker',
                               repeatmasker_path => $self->o('repeatmasker_path'),
                               commandline_params => '-nolow -lib "'.$self->o('repeatmasker_library').'" -engine "'.$self->o('repeatmasker_engine').'"',
                             },
               -rc_name    => 'local',
               -flow_into => {
                               1 => ['run_genscan'],
                             },

            },

            {
              # Run genscan, uses 200kb slices from repeatmasker. Flows into create_prediction_transcript_ids which
              # then takes these 1mb slices and converts them into individual prediction transcript input ids based
              # on the dbID of each feature generate by this analysis
              -logic_name => 'run_genscan',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveGenscan',
              -parameters => {
                               target_db => $self->o('core_db'),
                               logic_name => 'genscan',
                               module => 'HiveGenscan',
                               genscan_path => $self->o('genscan_path'),
                               genscan_matrix_path => $self->o('genscan_matrix_path'),
                               repeat_masking_logic_names => ['repeatmasker_'.$self->o('repeatmasker_library')],
                       },
              -rc_name    => 'local',
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
    return {
      'local' => {'LOCAL' => ''},
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
   my $config_default_hash = $config_group_hash->{'default'};
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
  };

  return($master_config_settings->{$config_group});

}

1;
