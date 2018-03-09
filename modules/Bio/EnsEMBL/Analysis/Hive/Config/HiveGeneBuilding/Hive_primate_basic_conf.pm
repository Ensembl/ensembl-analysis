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

package Hive_primate_basic_conf;

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

##########################################################################
#                                                                        #
# CHANGE STUFF HERE                                                      #
#                                                                        #
##########################################################################
        'pipeline_name'              => '',

        'user_r'                     => '',
        'user_w'                     => '',
        'password'                   => '',
        'port'                       => '',

        'pipe_db_name'                => '',
        'reference_db_name'           => '',
        'dna_db_name'                 => '',
        'genblast_output_db_name'     => '',
        'genewise_output_db_name'     => '',
        'exonerate_output_db_name'    => '',

        'pipe_db_server'             => '',
        'reference_db_server'        => '',
        'dna_db_server'              => '',
        'genblast_output_db_server'  => '',
        'genewise_output_db_server'  => '',
        'exonerate_output_db_server' => '',
        'killlist_db_server'         => '',
        'output_path'                => '',
        'genome_file'                => '',

        'repeat_masking_logic_names' => [],



##########################################################################
#                                                                        #
# MOSTLY STAYS CONSTANT, MOSTLY                                          #
#                                                                        #
##########################################################################

        'seqs_per_chunk'             => 10,
        'fastasplit_random_path'     => '/software/ensembl/bin/fastasplit_random',

        'human_taxon_id'             => '9606',
        'primates_taxon_id'          => '9443',
        'mammals_taxon_id'           => '40674',
        'vert_taxon_id'              => '7742',

        'uniprot_index_name'         => 'uniprot_index',
        'uniprot_db_name'            => 'uniprot_db',
        'uniprot_query_dir_name'     => 'uniprot_temp',
        'uniprot_genblast_batch_size' => 10,
        'uniprot_table_name'         => 'uniprot_sequences',

        'genblast_path'              => 'genblast',
        'genblast_eval'              => '1e-20',
        'genblast_cov'               => '0.5',
        'genblast_rechunk_dir_name'  => 'genblast_fail_chunks',

        'exonerate_path'             => '/software/ensembl/genebuild/usrlocalensemblbin/exonerate-0.9.0',
        'exonerate_pid'              => '50',
        'exonerate_cov'              => '50',
        'exonerate_rechunk_dir_name' => 'exonerate_fail_chunks',

        'default_mem'                => '900',
        'genblast_mem'               => '1900',
        'genblast_retry_mem'         => '4900',
        'exonerate_mem'              => '3900',
        'exonerate_retry_mem'        => '5900',


        'create_type' => 'clone',
        'driver' => 'mysql',
        'num_tokens' => 10,
        'user' => 'ensro',

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

        'genblast_output_db' => {
                           -dbname => $self->o('genblast_output_db_name'),
                           -host   => $self->o('genblast_output_db_server'),
                           -port   => $self->o('port'),
                           -user   => $self->o('user_w'),
                           -pass   => $self->o('password'),
                         },

        'genewise_output_db' => {
                           -dbname => $self->o('genewise_output_db_name'),
                           -host   => $self->o('genewise_output_db_server'),
                           -port   => $self->o('port'),
                           -user   => $self->o('user_w'),
                           -pass   => $self->o('password'),
                        },

        'exonerate_output_db' => {
                           -dbname => $self->o('exonerate_output_db_name'),
                           -host   => $self->o('exonerate_output_db_server'),
                           -port   => $self->o('port'),
                           -user   => $self->o('user_w'),
                           -pass   => $self->o('password'),
                         },


        'killlist_db' => {
                           -dbname    => $self->o('killlist_db_name'),
                           -host      => $self->o('killlist_db_server'),
                           -port      => $self->o('port'),
                           -user      => $self->o('user_r'),
                         },
    };
}



sub pipeline_create_commands {
    my ($self) = @_;
    return [
    # inheriting database and hive tables' creation
      @{$self->SUPER::pipeline_create_commands},

      $self->db_cmd('CREATE TABLE '.$self->o('uniprot_table_name').' ('.
                    'accession varchar(50) NOT NULL,'.
                    'source_db varchar(50) NOT NULL,'.
                    'pe_level varchar(50) NOT NULL,'.
                    'biotype varchar(255) NOT NULL,'.
                    'group_name varchar(255) NOT NULL,'.
                    'seq text NOT NULL,'.
                    'PRIMARY KEY (accession))'),
    ];
}


## See diagram for pipeline structure
sub pipeline_analyses {
    my ($self) = @_;

    return [

      {
        -logic_name => 'create_genblast_output_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('reference_db'),
                         target_db => $self->o('genblast_output_db'),
                         create_type => $self->o('create_type'),
                       },
        -rc_name    => 'default',
        -input_ids => [{}],
      },

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
        -logic_name => 'download_uniprot_files',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadUniProtFiles',
        -parameters => {
                         multi_query_download => {
                           human_pe12 => {
                             file_name => 'human_pe12.fasta',
                             taxon_id  => $self->o('human_taxon_id'),
                             dest_dir  => $self->o('output_path'),
                             compress  => 1,
                             pe_level  => [1,2],},

                           primates_pe12 => {
                             file_name  => 'primate_pe12.fasta',
                             taxon_id   => $self->o('primates_taxon_id'),
                             exclude_id => $self->o('human_taxon_id'),
                             dest_dir   => $self->o('output_path'),
                             compress   => 1,
                             pe_level   => [1,2],},

                           primates_pe345 => {
                             file_name => 'primate_pe345.fasta',
                             taxon_id  => $self->o('primates_taxon_id'),
                             dest_dir  => $self->o('output_path'),
                             compress  => 1,
                             pe_level  => [3,4,5],
                           },

                           mammals_pe12 => {
                             file_name  => 'mammal_pe12.fasta',
                             taxon_id   => $self->o('mammals_taxon_id'),
                             exclude_id => $self->o('primates_taxon_id'),
                             dest_dir   => $self->o('output_path'),
                             compress   => 1,
                             pe_level   => [1,2],
                           },

                           vert_pe12 => {
                             file_name  => 'vert_pe12.fasta',
                             taxon_id   => $self->o('vert_taxon_id'),
                             exclude_id => $self->o('mammals_taxon_id'),
                             dest_dir   => $self->o('output_path'),
                             compress   => 1,
                             pe_level   => [1,2],
                           },
                         },
                       },

        -rc_name          => 'default',
        -input_ids  => [ {} ],
        -flow_into => {
                        1 => ['process_uniprot_files'],
                      },
      },

      {
        -logic_name => 'process_uniprot_files',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProcessUniProtFiles',
        -parameters => {
                         uniprot_db_name => $self->o('uniprot_db_name'),
                         uniprot_index_name => $self->o('uniprot_index_name'),
                         dest_dir   => $self->o('output_path'),
                         killlist_type => 'protein',
                         killlist_db => $self->o('killlist_db'),
                      },
        -rc_name => 'default',
        -flow_into => {
                        1 => ['load_uniprot_seqs'],
                      },
      },

      {
        -logic_name => 'load_uniprot_seqs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadSeqsPipeDB',
        -parameters => {
                         uniprot_db_path => $self->o('output_path').'/'.$self->o('uniprot_db_name'),
                         uniprot_index_path => $self->o('output_path').'/'.$self->o('uniprot_index_name'),
                       },
        -rc_name          => 'default',
        -wait_for     => ['process_uniprot_files'],
        -flow_into => {
                        1 => ['exonerate'],
                      },
      },


      {
        -logic_name => 'generate_genblast_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         uniprot_accession => 1,
                         uniprot_batch_size => $self->o('uniprot_genblast_batch_size'),
                         uniprot_table_name => $self->o('uniprot_table_name'),
                       },
        -rc_name      => 'default',
        -wait_for     => ['load_uniprot_seqs'],
        -input_ids  => [ {} ],
        -flow_into => {
                        1 => ['genblast'],
                      },
      },

      {
        -logic_name => 'genblast',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
        -parameters => {
                         iid_type => 'db_seq',
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('genblast_output_db'),
                         logic_name => 'genblast',
                         module => 'HiveGenblast',
                         genblast_path => $self->o('genblast_path'),
                         genblast_db_path => $self->o('genome_file'),
                         commandline_params => ' -P wublast -gff -e '.$self->o('genblast_eval').' -c '.$self->o('genblast_cov').' ',
                         query_seq_dir => $self->o('output_path').'/'.$self->o('uniprot_query_dir_name'),
                       },
        -rc_name    => 'genblast',
        -wait_for => ['create_genblast_output_db'],
        -flow_into => {
                        -1 => ['genblast_himem'],
                      },
        -failed_job_tolerance => 50,
      },

      {
        -logic_name => 'genblast_himem',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
        -parameters => {
                         iid_type => 'db_seq',
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('genblast_output_db'),
                         logic_name => 'genblast',
                         module => 'HiveGenblast',
                         genblast_path => $self->o('genblast_path'),
                         genblast_db_path => $self->o('genome_file'),
                         commandline_params => ' -P wublast -gff -e '.$self->o('genblast_eval').' -c '.$self->o('genblast_cov').' ',
                         query_seq_dir => $self->o('output_path').'/'.$self->o('uniprot_query_dir_name'),
                       },
        -rc_name          => 'genblast_retry',
        -can_be_empty  => 1,
        -failed_job_tolerance => 100,
      },

      {
        -logic_name => 'exonerate',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
        -rc_name          => 'exonerate',
        -wait_for => ['create_exonerate_output_db'],
        -parameters => {
                         iid_type => 'db_seq',
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('exonerate_output_db'),
                         logic_name => 'exonerate',
                         module     => 'HiveExonerate2Genes',
                         config_settings => $self->get_config_settings('exonerate_protein','exonerate'),
                         query_seq_dir => $self->o('output_path').'/'.$self->o('uniprot_query_dir_name'),
                      },
        -flow_into => {
                        -1 => ['exonerate_himem'],
                      },
        -failed_job_tolerance => 50,
      },


      {
        -logic_name => 'exonerate_himem',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
        -parameters => {
                         iid_type => 'db_seq',
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('exonerate_output_db'),
                         logic_name => 'exonerate',
                         module     => 'HiveExonerate2Genes',
                         config_settings => $self->get_config_settings('exonerate_protein','exonerate'),
                         query_seq_dir => $self->o('output_path').'/'.$self->o('uniprot_query_dir_name'),
                       },
        -rc_name          => 'exonerate_retry',
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

  # Note that this builds resource requests off some of the variables at the top. Works off the idea
  # that the references all get put on one server and the pipe db is on another
  my $pipe_db_server = $self->default_options()->{'pipe_db_server'};
  my $dna_db_server = $self->default_options()->{'dna_db_server'};
  my $genblast_output_db_server = $self->default_options()->{'genblast_output_db_server'};
  my $exonerate_output_db_server = $self->default_options()->{'exonerate_output_db_server'};
  my $killlist_db_server = $self->default_options()->{'killlist_db_server'};

  my $default_mem = $self->default_options()->{'default_mem'};
  my $genblast_mem = $self->default_options()->{'genblast_mem'};
  my $genblast_retry_mem = $self->default_options()->{'genblast_retry_mem'};
  my $exonerate_mem = $self->default_options()->{'exonerate_mem'};
  my $exonerate_retry_mem = $self->default_options()->{'exonerate_retry_mem'};

  my $pipe_db_server_number;
  my $dna_db_server_number;
  my $genblast_output_db_server_number;
  my $exonerate_output_db_server_number;
  my $killlist_db_server_number;


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

  unless($genblast_output_db_server =~ /(\d+)$/) {
    die "Failed to parse the server number out of the pipeline db server name. This is needed for setting tokens\n".
        "genblast_output_db_server: ".$genblast_output_db_server;
  }

  $genblast_output_db_server_number = $1;

    unless($exonerate_output_db_server =~ /(\d+)$/) {
    die "Failed to parse the server number out of the pipeline db server name. This is needed for setting tokens\n".
        "exonerate_output_db_server: ".$exonerate_output_db_server;
  }

  $exonerate_output_db_server_number = $1;

    unless($killlist_db_server=~ /(\d+)$/) {
    die "Failed to parse the server number out of the pipeline db server name. This is needed for setting tokens\n".
        "killlist_db_server: ".$killlist_db_server;
  }

  $killlist_db_server_number = $1;

  unless($num_tokens) {
    die "num_tokens is uninitialised or zero. num_tokens needs to be present in default_options and not zero\n".
        "num_tokens: ".$num_tokens;
  }

    return {
      'default' => { LSF => '-q normal -M'.$default_mem.' -R"select[mem>'.$default_mem.'] '.
                            'rusage[mem='.$default_mem.','.
                            'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

      'genblast' => { LSF => '-q normal -M'.$genblast_mem.' -R"select[mem>'.$genblast_mem.'] '.
                             'rusage[mem='.$genblast_mem.','.
                             'myens_build'.$genblast_output_db_server_number.'tok='.$num_tokens.','.
                             'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
                             'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

      'genblast_retry' => { LSF => '-q normal -M'.$genblast_retry_mem.' -R"select[mem>'.$genblast_retry_mem.'] '.
                                   'rusage[mem='.$genblast_retry_mem.','.
                                   'myens_build'.$genblast_output_db_server_number.'tok='.$num_tokens.','.
                                   'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
                                   'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

      'exonerate' => { LSF => '-q normal -M'.$exonerate_mem.' -R"select[mem>'.$exonerate_mem.'] '.
                              'rusage[mem='.$exonerate_mem.','.
                              'myens_build'.$exonerate_output_db_server_number.'tok='.$num_tokens.','.
                              'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
                              'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

      'exonerate_retry' => { LSF => '-q normal -M'.$exonerate_retry_mem.' -R"select[mem>'.$exonerate_retry_mem.'] '.
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

    killlist_protein => {
                         KILLLISTDB          => $self->o('killlist_db'),
                         USE_KILL_LIST       => 1,
                         KILL_TYPE           => 'protein',
                         KILL_LIST_FILTER    => {
                                                  -only_mol_type        => 'protein',
                                                  -user_id              => undef,
                                                  -from_source_species  => undef,
                                                  -before_date          => undef,
                                                  -having_status        => undef,
                                                  -reasons              => [],
                                                  -for_analyses         => [],
                                                  -for_species          => [],
                                                  -for_external_db_ids  => [],
                                                },
                       },
    },
  };

  return($master_config_settings->{$config_group});

}

1;
