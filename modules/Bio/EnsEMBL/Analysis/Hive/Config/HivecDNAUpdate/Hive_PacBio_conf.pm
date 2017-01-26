=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016] EMBL-European Bioinformatics Institute

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


package Hive_PacBio_conf;

use strict;
use warnings;
use feature 'say';

use Bio::EnsEMBL::Analysis::Tools::Utilities qw (get_analysis_settings);
use parent ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

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
    'pipe_dbname'               => $ENV{USER}.'_'.$self->o('pipeline_name').'_hive',
    'pipe_db_server'             => 'genebuild11',

    'dna_dbname'                => $ENV{USER}.'_human_83_copy',
    'dna_db_server'              => 'genebuild12',

    'killlist_db_server'         => 'genebuild6',

    'exonerate_output_db_name'   => $ENV{USER}.'_'.$self->o('pipeline_name').'_exonerate',
    'exonerate_output_db_server' => 'genebuild13',

    'output_path'                => '',

    'exonerate_batch_size'       => '50',

    'ensembl_repo_root'          => $ENV{ENSCODE},
    'clone_db_script_path'       => $self->o('ensembl_repo_root').'/ensembl-analysis/scripts/clone_database.ksh',

    'genome_file'                => '/data/blastdb/Ensembl/human/GRCh38/genome/softmasked_dusted/toplevel.with_nonref_and_GRCh38_p5.no_duplicate.softmasked_dusted.fa',

    'repeat_masking_logic_names' => ['repeatmask_repbase_human'],


##########################################################################
#                                                                        #
# MOSTLY STAYS CONSTANT, MOSTLY                                          #
#                                                                        #
##########################################################################

    'fastasplit_random_path'     => '/software/ensembl/bin/fastasplit_random',
    'killlist_db_name'           => 'gb_kill_list',

    'cdna_file_name'             => 'cdna_update',

    'user_r'                     => 'ensro',
    'user'                     => 'ensadmin',
    'password'                   => '',
    'port'                       => '3306',

    'cdna_query_dir_name'        => 'cdna_temp',

    'many_hits_dir'              => 'many_hits',

    'cdna_table_name'            => 'cdna_sequences',
    'cdna_batch_size'            => '1',

    'default_mem'                => '900',
    'exonerate_mem'              => '3900',
    'exonerate_retry_mem'        => '5900',

    'exonerate_path'             => '/software/ensembl/genebuild/usrlocalensemblbin/exonerate-0.9.0',
    'exonerate_pid'              => '97',
    'exonerate_cov'              => '90',

    'driver'                     => 'mysql',
    'num_tokens'                 => 10,

    'create_type'                => 'clone',
    production_dbname => 'ensembl_production',
    production_server => 'ens-staging1',
    production_ro_user => $self->o('user_r'),

    'exonerate_output_db' => {
      -dbname => $self->o('exonerate_output_db_name'),
      -host => $self->o('exonerate_output_db_server'),
      -port => $self->o('port'),
      -user => $self->o('user'),
      -pass => $self->o('password'),
    },

    'killlist_db' => {
      -dbname => $self->o('killlist_db_name'),
      -host => $self->o('killlist_db_server'),
      -port => $self->o('port'),
      -user => $self->o('user_r'),
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
      'biotype varchar(50) NOT NULL,'.
      'PRIMARY KEY (accession))'),
  ];
}

## See diagram for pipeline structure
sub pipeline_analyses {
  my ($self) = @_;

  return [
    {
      # need to make sure the database actually copies as if it doesn't the job does appear to complete according to eHive
      -logic_name => 'create_output_db',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db => $self->o('dna_db'),
        target_db => $self->o('exonerate_output_db'),
        create_type => $self->o('create_type'),
        script_path => $self->o('clone_db_script_path'),
      },
      -rc_name => 'default',
      -input_ids => [{
        cdna_file => $self->o('output_path').'/'.$self->o('cdna_file_name'),
      }],
      -max_retry_count => 0,
      -flow_into => {
        1 => ['populate_production'],
      }
    },
    {
      -logic_name => 'populate_production',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('ensembl_repo_root').'/ensembl-production/scripts/production_database/populate_production_db_tables.pl'.
               ' -dp '.$self->o('output_path').
               ' -d '.$self->o('exonerate_output_db','-dbname').
               ' -h '.$self->o('exonerate_output_db','-host').
               ' -u '.$self->o('exonerate_output_db','-user').
               ' -p '.$self->o('exonerate_output_db','-pass').
               ' -md '.$self->o('production_dbname').' -mh '.$self->o('production_server').' -mu '.$self->o('production_ro_user').
               ' -t external_db -t attrib_type -t misc_set -t unmapped_reason'
      },
      -max_retry_count => 0,
      -flow_into => {
        1 => ['load_cdnas'],
      }
    },
    {
      -logic_name => 'load_cdnas',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadcDNAs',
      -rc_name => 'default',
      -max_retry_count => 0,
      -flow_into => {
        1 => ['generate_jobs'],
      }
    },
    {
      -logic_name => 'generate_jobs',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        cdna_accession => 1,
        cdna_batch_size => $self->o('cdna_batch_size'),
        cdna_table_name => $self->o('cdna_table_name'),
      },
      -rc_name => 'default',
      -flow_into => {
        2 => ['exonerate'],
      },
    },
    {
      -logic_name => 'exonerate',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes_cdna',
      -parameters => {
        iid_type => 'db_seq',
        dna_db => $self->o('dna_db'),
        target_db => $self->o('exonerate_output_db'),
        logic_name => 'cdna_update',
        module => 'HiveExonerate2Genes_cdna',
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate')},
        query_seq_dir => $self->o('output_path').'/'.$self->o('cdna_query_dir_name'),
        stdout_file => $self->o('output_path').'/exonerate_1.out',
        GENOMICSEQS         => $self->o('genome_file'),
        PROGRAM             => $self->o('exonerate_path'),
        SOFT_MASKED_REPEATS => $self->o('repeat_masking_logic_names'),
      },
      -flow_into => {
        2 => [ 'exonerate_second_run' ],
        -1 => ['exonerate_himem'],
      },
      -rc_name => 'exonerate',
      -failed_job_tolerance => 50,
      -batch_size => $self->o('exonerate_batch_size'),
    },
    {
      -logic_name => 'exonerate_himem',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes_cdna',
      -parameters => {
        iid_type => 'db_seq',
        dna_db => $self->o('dna_db'),
        target_db => $self->o('exonerate_output_db'),
        logic_name => 'cdna_update',
        module => 'HiveExonerate2Genes_cdna',
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate')},
        query_seq_dir => $self->o('output_path').'/'.$self->o('cdna_query_dir_name'),
        stdout_file => $self->o('output_path').'/exonerate_himem.out',
        GENOMICSEQS         => $self->o('genome_file'),
        PROGRAM             => $self->o('exonerate_path'),
        SOFT_MASKED_REPEATS => $self->o('repeat_masking_logic_names'),
      },
      -flow_into => {
        2 => ['exonerate_second_run'],
      },
      -rc_name => 'exonerate_himem',
      -can_be_empty => 1,
      -failed_job_tolerance => 100,
    },
    {
      -logic_name => 'exonerate_second_run',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes_cdna',
      -parameters => {
        iid_type => 'db_seq',
        dna_db => $self->o('dna_db'),
        target_db => $self->o('exonerate_output_db'),
        logic_name => 'cdna_update',
        module => 'HiveExonerate2Genes_cdna',
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_2')},
        query_seq_dir => $self->o('output_path').'/'.$self->o('cdna_query_dir_name'),
        stdout_file => $self->o('output_path').'/exonerate_2.out',
        GENOMICSEQS         => $self->o('genome_file'),
        PROGRAM             => $self->o('exonerate_path'),
        SOFT_MASKED_REPEATS => $self->o('repeat_masking_logic_names'),
      },
      -rc_name => 'exonerate_2',
      -can_be_empty => 1,
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

sub resource_classes {
  my $self = shift;

  return {
    'default' => { LSF => $self->lsf_resource_builder('normal', $self->default_options->{'default_mem'}, [$self->default_options->{'pipe_db_server'}], [$self->default_options->{'num_tokens'}]) },
    'exonerate' => { LSF => $self->lsf_resource_builder('normal', $self->default_options->{'exonerate_mem'}, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'exonerate_output_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]) },
    'exonerate_himem' => { LSF => $self->lsf_resource_builder('normal', $self->default_options->{'exonerate_retry_mem'}, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'exonerate_output_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]) },
    'exonerate_2' => { LSF => $self->lsf_resource_builder('normal', $self->default_options->{'exonerate_retry_mem'}, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'exonerate_output_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]) },
  }
}

1;
