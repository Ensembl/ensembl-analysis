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

# this draft will just do the first part of the pipeline - downlaod the fasta files, remove kill list obj and then clean and clip


package Hive_cDNA_update_conf;

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

    'pipeline_name'              => 'human_cdna_hive',

    'reference_dbname'           => 'homo_sapiens_core_82_38',
    'reference_db_server'        => 'ens-livemirror',

    'species'                    => 'human',
    'user_r'                     => 'ensro',
    'user_w'                     => 'ensadmin',
    'password'                   => 'ensembl',
    'port'                       => '3306',

    'pipe_dbname'                => 'dm15_human_cdna_hive',
    'pipe_db_server'             => 'genebuild11',

    'killlist_db_server'         => 'genebuild6',
    'output_path'                => '/lustre/scratch109/ensembl/dm15/hive_cdna/',
    'gss_file'                   => '/nfs/users/nfs_d/dm15/cvs_checkout_head/ensembl-personal/genebuilders/cDNA_update/gss_acc.txt',
    'polyA_script'               => '~/enscode/ensembl-pipeline/scripts/EST/new_polyA_clipping.pl',

    'refseq_path'                => '/data/blastdb/Ensembl/RefSeq_2015_06/',
    'refseq_file'                => 'hs.fna',

    'embl_file'                  => 'embl_9606.fa',

    'output_dbname'              => 'dm15_human_cdnatest',
    'output_db_server'           => 'genebuild13',

    'clone_db_script_path'       => '/nfs/users/nfs_d/dm15/cvs_checkout_head/ensembl-personal/genebuilders/scripts/clone_database.ksh',

    'genomic_seq'                => '/data/blastdb/Ensembl/human/GRCh38/genome/softmasked_dusted/toplevel.with_nonref_and_GRCh38_p3.no_duplicate.softmasked_dusted.fa',

##########################################################################
#                                                                        #
# MOSTLY STAYS CONSTANT, MOSTLY                                          #
#                                                                        #
##########################################################################

    'seqs_per_chunk'             => 10,
    'fastasplit_random_path'     => '/software/ensembl/bin/fastasplit_random',
    'killlist_dbname'            => 'gb_kill_list',

    'cdna_file_name'             => 'cdna_update',

    'human_taxon_id'             => '9606',
    'mouse_taxon_id'             => '10090',

    #'cdna_query_dir_name'        => 'cdna_temp',

    #'cdna_rechunk_dir_name'      => 'cdna_fail_chunks',

    'cdna_table_name'            => 'cdna_sequences',
    'cdna_batch_size'            => '10',

    'default_mem'                => '900',

    #'create_type' => 'clone',
    'driver' => 'mysql',
    'num_tokens' => 10,
    'user' => 'ensro',

    'create_type' => 'clone',

    'pipeline_db' => {
      -dbname => $self->o('pipe_dbname'),
      -host   => $self->o('pipe_db_server'),
      -port   => $self->o('port'),
      -user   => $self->o('user_w'),
      -pass   => $self->o('password'),
      -driver => $self->o('driver'),
    },
 
    'reference_db' => {
      -dbname => $self->o('reference_dbname'),
      -host   => $self->o('reference_db_server'),
      -port   => $self->o('port'),
      -user   => $self->o('user_r'),
    },

    'cdna_output_db' => {
      -dbname => $self->o('output_dbname'),
      -host   => $self->o('output_db_server'),
      -port   => $self->o('port'),
      -user   => $self->o('user_w'),
      -pass   => $self->o('password'),
      -driver => $self->o('driver'),
    },

    'killlist_db' => {
      -dbname    => $self->o('killlist_dbname'),
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

    $self->db_cmd('CREATE TABLE '.$self->o('cdna_table_name').' ('.
      'accession varchar(50) NOT NULL,'.
      'seq text NOT NULL,'.
      'PRIMARY KEY (accession))'),
  ];
}

## See diagram for pipeline structure
sub pipeline_analyses {
  my ($self) = @_;

  return [
    {
      -logic_name => 'create_output_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db => $self->o('reference_db'),
        target_db => $self->o('cdna_output_db'),
        create_type => $self->o('create_type'),
        script_path => $self->o('clone_db_script_path'),
      },
      -rc_name    => 'default',
      -input_ids => [{}],
    },
    {
      -logic_name => 'download_cdnas',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadcDNAFiles',
      -parameters => {
        embl_sequences => {
          output_path => $self->o('output_path'),
          output_file => $self->o('cdna_file_name'),
          species => $self->o('species'),
          taxon_id => $self->o('human_taxon_id'),
        },
        refseq_sequences => {
          refseq_path => $self->o('refseq_path'),
          refseq_file => $self->o('refseq_file'),
        },
      },  
      -rc_name   => 'default',
      -input_ids => [ {} ],
      -flow_into => {
        1 => ['prepare_cdnas'],
      },
    },
    {
      -logic_name => 'prepare_cdnas',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivePreparecDNAs',
      -parameters => {
        prepare_seqs => {
          dest_dir   => $self->o('output_path'),
          embl_file => $self->o('embl_file'),
          refseq_file => $self->o('refseq_file'),
          gss_file => $self->o('gss_file'),
          killlist_type => 'cdna_update',
          killlist_db => $self->o('killlist_db'),
          polyA_script => $self->o('polyA_script'),
          cdna_file => $self->o('cdna_file_name'),
          species => $self->o('species'),
        },
      },
      -rc_name => 'default',
      -flow_into => {
        1 => ['load_cdnas'],
      },
    },
    {
      -logic_name => 'load_cdnas',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadcDNAs',
      -parameters => {
        cdna_file => $self->o('output_path').'/'.$self->o('cdna_file_name').'.clipped',
        species => $self->o('species'),
      },
      -rc_name => 'default',
      -wait_for => ['prepare_cdnas'],
#      -flow_into => {
#        1 => ['generate_jobs'],
#      },
    },
    {
      -logic_name => 'generate_jobs',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        cdna_accession => 1,
        cdna_batch_size => $self->o('cdna_batch_size'),
        cdna_table_name => $self->o('cdna_table_name'),
      },
      -rc_name => 'default',
      -wait_for => ['load_cdnas'],
      -input_ids => [ {} ],
      -flow_into => {
        1 => ['exonerate'],
      },
    },
    {
      -logic_name => 'exonerate',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
#      -rc_name          => 'exonerate',
      -wait_for => ['create_output_db'],
      -parameters => {
        iid_type => 'db_seq',
##      dna_db => $self->o('dna_db'),
#        target_db => $self->o('output_db'),
#        logic_name => 'exonerate',
        module     => 'HiveExonerate2Genes',
#        config_settings => $self->get_config_settings('exonerate_protein','exonerate'),
#      query_seq_dir => $self->o('output_path').'/'.$self->o('uniprot_query_dir_name'),
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
#      dna_db => $self->o('dna_db'),
#        target_db => $self->o('output_db'),
#        logic_name => 'exonerate',
        module     => 'HiveExonerate2Genes',
#        config_settings => $self->get_config_settings('exonerate_protein','exonerate'),
#      query_seq_dir => $self->o('output_path').'/'.$self->o('uniprot_query_dir_name'),
      },
#      -rc_name          => 'exonerate_retry',
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

#sub resource_classes {
#  my $self = shift;
  # Note that this builds resource requests off some of the variables at the top. Works off the idea
  # that the references all get put on one server and the pipe db is on another
#  my $pipe_db_server = $self->default_options()->{'pipe_db_server'};
#  my $killlist_db_server = $self->default_options()->{'killlist_db_server'};

#  my $default_mem = $self->default_options()->{'default_mem'};

#  my $pipe_db_server_number;
#  my $killlist_db_server_number;

#  my $num_tokens = $self->default_options()->{'num_tokens'};

#  unless($pipe_db_server =~ /(\d+)$/) {
#    die "Failed to parse the server number out of the pipeline db server name. This is needed for setting tokens\n".
#        "pipe_db_server: ".$pipe_db_server;
#  }

#  $pipe_db_server_number = $1;

#  unless($killlist_db_server=~ /(\d+)$/) {
#    die "Failed to parse the server number out of the pipeline db server name. This is needed for setting tokens\n".
#        "killlist_db_server: ".$killlist_db_server;
#  }

#  $killlist_db_server_number = $1;

#  unless($num_tokens) {
#    die "num_tokens is uninitialised or zero. num_tokens needs to be present in default_options and not zero\n".
#        "num_tokens: ".$num_tokens;
#  }

#  return {
#    'default' => { LSF => '-q normal -M'.$default_mem.' -R"select[mem>'.$default_mem.'] '.
#    'rusage[mem='.$default_mem.','.
#    'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

#      'exonerate_retry' => { LSF => '-q normal -M'.$exonerate_retry_mem.' -R"select[mem>'.$exonerate_retry_mem.'] '.
#                                    'rusage[mem='.$exonerate_retry_mem.','.
#                                    'myens_build'.$exonerate_output_db_server_number.'tok='.$num_tokens.','.
#                                    'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
#                                    'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },
#  }
#}

sub resource_classes {
  my ($self) = @_;

  return {
    'default' => {'LSF' => "-q normal -M 1000 -R 'select[mem>=1000] rusage[mem=1000]' " },
  };
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
      killlist_cdna => {
        KILLLISTDB          => $self->o('killlist_db'),
        USE_KILL_LIST       => 1,
        KILL_TYPE           => 'cdna',
        KILL_LIST_FILTER    => {
          -only_mol_type        => 'cdna',
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
