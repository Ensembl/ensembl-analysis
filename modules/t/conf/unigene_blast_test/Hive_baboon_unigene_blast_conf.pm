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

package Hive_baboon_unigene_blast_conf;

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

        'pipeline_name' => 'BaboonBlastTest',
        'pipe_db_server' => 'genebuild13',
        'reference_db_server' => 'genebuild11',
        'user' => 'ensadmin',
        'password' => 'ensembl',
        'port' => 3306,
        'driver' => 'mysql',

        'pipe_dbname' => 'fm2_baboon_unigene_blast_test_pipe',
        'pipeline_db' => {
            # connection parameters
            -host   => $self->o('pipe_db_server'),
            -port   => $self->o('port'),
            -user   => $self->o('user'),
            -pass   => $self->o('password'),
            -dbname => $self->o('pipe_dbname'),
            -driver => $self->o('driver'),
        },

        'reference_db' => {
                            -dbname => 'fm2_papio_anubis_ref',
                            -host   => $self->o('reference_db_server'),
                            -port   => $self->o('port'),
                            -user   => $self->o('user'),
                            -pass   => $self->o('password'),
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
        -logic_name => 'run_unigene_blast_genscan',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlastGenscanDNA',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         blast_exe_path => '/software/ensembl/genebuild/usrlocalensemblbin/wutblastn',
                         blast_program => 'wutblastn',
                         commandline_params => '-cpus 3 -hitdist 40',
                         repeat_masking_logic_names => ['run_repeatmasker_wublast'],
                         prediction_transcript_logic_names => ['run_genscan_wu'],
                         logic_name => 'run_unigene_blast_genscan',
                         module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlastGenscanDNA',
                         config_settings => $self->get_config_settings('HiveBlast','HiveBlastGenscanDNA'),
                      },
        -rc_name    => 'blast',
        -input_ids => [{"iid" => "chromosome:Panu_2.0:4:140000:200000:1"}],
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
      'default' => { LSF => '-q normal -M900 -R"select[mem>900] rusage[mem=900,myens_build11tok=10,myens_build12tok=10,myens_build13tok=10]"' },
      'blast' => { LSF => '-q normal -M2900 -n 3 -R "select[mem>2900] rusage[mem=2900,myens_build11tok=10,myens_build12tok=10,myens_build13tok=10] span[hosts=1]"' },
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
  HiveBlast => {
    HiveBlastGenscanDNA => {
                       BLAST_PARAMS => {
                         -unknown_error_string => 'FAILED',
                         -type => 'wu',
                       },
                       BLAST_PARSER => 'Bio::EnsEMBL::Analysis::Tools::FilterBPlite',
                       PARSER_PARAMS => {
                        -regex => '\/ug\=([\w\.]+)',
                        -query_type => 'pep',
                        -database_type => 'dna',
                        -threshold_type => 'PVALUE',
                        -threshold => 0.001,
                       },
                       BLAST_FILTER => 'Bio::EnsEMBL::Analysis::Tools::FeatureFilter',
                       FILTER_PARAMS => {
                        -prune => 1,
                       },
               },

    },
  };

  return($master_config_settings->{$config_group});

}

1;
