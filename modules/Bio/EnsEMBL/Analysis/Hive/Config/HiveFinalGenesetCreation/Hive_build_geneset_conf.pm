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

package Hive_build_geneset_conf;

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

#############################################################
#                                                           #
# THINGS TO FILL IN                                         #
#                                                           #
#############################################################

        'pipeline_name' => '',

        'user_r'                    => '',
        'user_w'                    => '',
        'password'                  => '',
        'port'                      => '',

        'pipe_db_name'               => '',
        'reference_db_name'          => '',
        'dna_db_name'                => '',
        'genblast_db_name'           => '',
        'genewise_db_name'           => '',
        'exonerate_db_name'          => '',
        'rnaseq_db_name'             => '',
        'projection_db_name'         => '',
        'initial_cluster_db_name'    => '',
        'final_geneset_db_name'      => '',

        'pipe_db_server'            => '',
        'reference_db_server'       => '',
        'dna_db_server'             => '',
        'genblast_db_server'        => '',
        'genewise_db_server'        => '',
        'exonerate_db_server'       => '',
        'rnaseq_db_server'          => '',
        'projection_db_server'      => '',
        'initial_cluster_db_server' => '',
        'final_geneset_db_server'   => '',

        'single_exon_support_penalty' => 2,

        'input_gene_dbs' => [$self->o('genblast_db'),
                             $self->o('genewise_db'),
                             $self->o('exonerate_db'),
                             $self->o('projection_db'),
                             $self->o('rnaseq_db'),
                            ],


        'allowed_input_sets' => {'rnaseq_blast'          => {'rnaseq_80_100' => 1,
                                                             'rnaseq_50_80'  => 1,
                                                             'rnaseq_0_50'   => 1,
                                                            },
                                 'genblast_human'        => 1,
                                 'genblast_primates'     => 1,
                                 'genblast_mammals'      => 1,
                                 'genblast_vert'         => 1,
                                 'genblast_primates_345' => 1,
                                },

        'logic_name_weights' => {'rnaseq_blast'          => {'rnaseq_80_100' => 1,
                                                             'rnaseq_50_80'  => 5,
                                                             'rnaseq_0_50'   => 7,
                                                            },
                                 'genblast_human'        => 2,
                                 'genblast_primates'     => 3,
                                 'genblast_mammals'      => 4,
                                 'genblast_vert'         => 5,
                                 'genblast_primates_345' => 8,
                                },


#############################################################
#                                                           #
# THINGS THAT MOSTLY DON'T CHANGE                           #
#                                                           #
#############################################################


        'create_type' => 'clone',
        'driver' => 'mysql',

        'pipeline_db' => {
            -dbname => $self->o('pipe_db_name'),
            -host   => $self->o('pipe_db_server'),
            -port   => $self->o('port'),
            -user   => $self->o('user_w'),
            -pass   => $self->o('password'),
            -driver => $self->o('driver'),
        },

        'reference_db' => {
                            -dbname => $self->o('reference_db_name'),
                            -host   => $self->o('reference_db_server'),
                            -port   => $self->o('port'),
                            -user   => $self->o('user_r'),
                          },

        'dna_db' => {
                      -dbname => $self->o('dna_db_name'),
                      -host   => $self->o('dna_db_server'),
                      -port   => $self->o('port'),
                      -user   => $self->o('user_r'),
                    },

        'genblast_db' => {
                           -dbname => $self->o('genblast_db_name'),
                           -host   => $self->o('genblast_db_server'),
                           -port   => $self->o('port'),
                           -user   => $self->o('user_r'),
                         },

        'genewise_db' => {
                           -dbname => $self->o('genewise_db_name'),
                           -host   => $self->o('genewise_db_server'),
                           -port   => $self->o('port'),
                           -user   => $self->o('user_r'),
                        },

        'exonerate_db' => {
                           -dbname => $self->o('exonerate_db_name'),
                           -host   => $self->o('exonerate_db_server'),
                           -port   => $self->o('port'),
                           -user   => $self->o('user_r'),
                         },

        'projection_db' => {
                            -dbname => $self->o('projection_db_name'),
                            -host   => $self->o('projection_db_server'),
                            -port   => $self->o('port'),
                            -user   => $self->o('user_r'),
                          },

        'rnaseq_db' => {
                            -dbname => $self->o('rnaseq_db_name'),
                            -host   => $self->o('rnaseq_db_server'),
                            -port   => $self->o('port'),
                            -user   => $self->o('user_r'),
                          },

        'initial_cluster_db' => {
                                 -dbname => $self->o('initial_cluster_db_name'),
                                 -host   => $self->o('initial_cluster_db_server'),
                                 -port   => $self->o('port'),
                                 -user   => $self->o('user_w'),
                                 -pass   => $self->o('password'),
                               },

         'final_geneset_db' => {
                                -dbname => $self->o('final_geneset_db_name'),
                            -host   => $self->o('final_geneset_db_server'),
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
    ];
}


## See diagram for pipeline structure
sub pipeline_analyses {
    my ($self) = @_;

    return [

      {
        -logic_name => 'create_initial_cluster_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('reference_db'),
                         target_db => $self->o('initial_cluster_db'),
                         create_type => $self->o('create_type'),
                       },
        -rc_name    => 'default',
        -input_ids => [{}],
      },

      {
        -logic_name => 'create_final_geneset_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('reference_db'),
                         target_db => $self->o('final_geneset_db'),
                         create_type => $self->o('create_type'),
                       },
        -rc_name    => 'default',
        -input_ids => [{}],
      },

     {
        -logic_name => 'create_toplevel_slices',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('dna_db'),
                         coord_system_name => 'toplevel',
                         slice => 1,
                         include_non_reference => 0,
                         top_level => 1,
                       },
        -flow_into => {
                       1 => ['cluster_input_genes'],
                      },
        -rc_name    => 'default',
        -input_ids => [{}],
      },

      {
        -logic_name => 'cluster_input_genes',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClusterSourceGenes',
        -parameters => {
                         logic_name => 'cluster_input_genes',
                         dna_db => $self->o('dna_db'),
                         output_db => $self->o('initial_cluster_db'),
                         input_gene_dbs => $self->o('input_gene_dbs'),
                         allowed_input_sets => $self->o('allowed_input_sets'),
                       },
        -rc_name    => 'cluster_input_genes',
        -wait_for => ['create_initial_cluster_db','create_final_geneset_db'],
        -flow_into => {
                        1 => ['fan_toplevel_slices'],
                      },
      },

      {
        -logic_name => 'fan_toplevel_slices',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('dna_db'),
                         split_slice => 1,
                         slice_size => 50000000,
                        },

        -rc_name    => 'default',
        -flow_into => {
                        1 => ['finalise_geneset'],
                      },
      },

      {
        -logic_name => 'finalise_geneset',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProcessGeneClusters',
        -parameters => {
                         logic_name => 'finalise_geneset',
                         cluster_db => $self->o('initial_cluster_db'),
                         processed_cluster_db => $self->o('final_geneset_db'),
                         dna_db => $self->o('reference_db'),
                         iid_type => 'slice',
                         logic_name_weights => $self->o('logic_name_weights'),
                         single_exon_support_penalty => $self->o('single_exon_support_penalty'),
                       },
        -rc_name    => 'finalise_geneset',
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
      'default' => { LSF => '-q normal -M1900 -R"select[mem>1900] rusage[mem=1900,myens_build2tok=10,myens_build3tok=10,myens_build13tok=10]"' },
      'cluster_input_genes' => { LSF => '-q normal -M1900 -R"select[mem>1900] rusage[mem=1900,myens_build2tok=10,myens_build3tok=10,myens_build4tok=10]"' },
      'finalise_geneset' => { LSF => '-q normal -M1900 -R"select[mem>1900] rusage[mem=1900,myens_build2tok=10,myens_build5tok=10,myens_build6tok=10]"' },
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
  my $master_config_settings = {};

  return($master_config_settings->{$config_group});

}

1;
