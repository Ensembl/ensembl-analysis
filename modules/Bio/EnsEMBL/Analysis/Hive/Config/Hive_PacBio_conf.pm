=head1 LICENSE

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

# this draft will just do the first part of the pipeline - downlaod the fasta files, remove kill list obj and then clean and clip


package Hive_PacBio_conf;

use strict;
use warnings;

use File::Spec::Functions qw(catfile);

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

    'pipeline_name' => '',
    'pipe_db_name'   => $self->o('dbowner').'_'.$self->o('pipeline_name').'_hive',

    'host'     => '',
    'user_r'   => '',
    'user'     => '',
    'password' => '',
    'port'     => '',

    'dna_db_name'    => '',
    'dna_db_server' => '',
    'dna_db_port'   => '',

    'killlist_db_server' => 'mysql-ens-genebuild-prod-6',
    'killlist_db_port'   => '4532',

    'exonerate_db_name'   => $self->o('dbowner').'_'.$self->o('pipeline_name').'_exonerate',
    'exonerate_db_server' => '',
    'exonerate_db_port'   => '4527',

    'exonerate_batch_size' => '50',
    'exonerate_path'       => catfile($self->o('software_base_path'), 'opt', 'exonerate09', 'bin', 'exonerate'), #You may need to specify the full path to the exonerate binary version 0.9.0
    'genome_file'          => '',

    'repeat_masking_logic_names' => [],
    'isoseq_dir' => '',

##########################################################################
#                                                                        #
# MOSTLY STAYS CONSTANT, MOSTLY                                          #
#                                                                        #
##########################################################################

    'isoseq_table_name'            => 'cdna_sequences',
    'batch_size'            => '1',

    'create_type'                => 'clone',

    'exonerate_db' => {
      -dbname => $self->o('exonerate_db_name'),
      -host => $self->o('exonerate_db_server'),
      -port => $self->o('exonerate_db_port'),
      -user => $self->o('user'),
      -pass => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'killlist_db' => {
      -dbname => $self->o('killlist_db_name'),
      -host   => $self->o('killlist_db_server'),
      -port   => $self->o('killlist_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -driver => $self->o('hive_driver'),
    },

    databases_to_delete => ['exonerate_db'],
  };
}

sub pipeline_create_commands {
  my ($self) = @_;
  return [
  # inheriting database and hive tables' creation
    @{$self->SUPER::pipeline_create_commands},
    $self->hive_data_table('refseq', $self->o('isoseq_table_name')),
#    $self->db_cmd('CREATE TABLE '.$self->o('isoseq_table_name').' ('.
#      'accession varchar(40) NOT NULL,'.
#      'seq text NOT NULL,'.
#      'biotype varchar(40) NOT NULL,'.
#      'PRIMARY KEY (accession))'),
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
        target_db => $self->o('exonerate_db'),
        create_type => $self->o('create_type'),
      },
      -rc_name => 'default',
      -input_ids => [{cdna_file => $self->o('cdna_file')}],
      -max_retry_count => 0,
      -flow_into => {
        1 => ['load_isoseqs'],
      }
    },
    {
      -logic_name => 'load_isoseqs',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadcDNAs',
      -rc_name => 'default',
      -parameters => {
        process_polyA => 1,
      },
      -max_retry_count => 0,
      -flow_into => {
        1 => ['generate_jobs'],
      }
    },
    {
      -logic_name => 'generate_jobs',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        iid_type => 'sequence_accession',
        batch_size => $self->o('batch_size'),
        sequence_table_name => $self->o('isoseq_table_name'),
      },
      -rc_name => 'exonerate_himem',
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
        target_db => $self->o('exonerate_db'),
        logic_name => 'isoseq',
        exonerate_cdna_cov => 90,
        exonerate_cdna_pid => 95,
        exonerate_bestn => 1,
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_cov_per_bestn_sub')},
      },
      -flow_into => {
        2 => [ 'exonerate_second_run' ],
        -1 => ['exonerate_himem'],
      },
      -rc_name => 'exonerate',
      -batch_size => $self->o('exonerate_batch_size'),
    },
    {
      -logic_name => 'exonerate_himem',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes_cdna',
      -parameters => {
        iid_type => 'db_seq',
        dna_db => $self->o('dna_db'),
        target_db => $self->o('exonerate_db'),
        logic_name => 'isoseq',
        exonerate_cdna_cov => 90,
        exonerate_cdna_pid => 95,
        exonerate_bestn => 1,
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_cov_per_bestn_sub')},
      },
      -flow_into => {
        2 => ['exonerate_second_run'],
      },
      -rc_name => 'exonerate_himem',
      -can_be_empty => 1,
    },
    {
      -logic_name => 'exonerate_second_run',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes_cdna',
      -parameters => {
        iid_type => 'db_seq',
        dna_db => $self->o('dna_db'),
        target_db => $self->o('exonerate_db'),
        logic_name => 'isoseq',
        exonerate_cdna_cov => 90,
        exonerate_cdna_pid => 95,
        exonerate_bestn => 1,
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_cov_per_bestn_loose_sub')},
      },
      -rc_name => 'exonerate_2',
      -can_be_empty => 1,
    },
  ];
}

sub pipeline_wide_parameters {
  my ($self) = @_;
  return {
    %{ $self->SUPER::pipeline_wide_parameters() },  # inherit other stuff from the base class
    GENOMICSEQS         => $self->o('genome_file'),
    PROGRAM             => $self->o('exonerate_path'),
    SOFT_MASKED_REPEATS => $self->o('repeat_masking_logic_names'),
  };
}

sub resource_classes {
  my $self = shift;

  return {
    'default' => { LSF => $self->lsf_resource_builder('production-rh7', 900, [$self->default_options->{'pipe_db_server'}]) },
    'exonerate' => { LSF => $self->lsf_resource_builder('production-rh7', 3900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'exonerate_output_db_server'}, $self->default_options->{'dna_db_server'}]) },
    'exonerate_himem' => { LSF => $self->lsf_resource_builder('production-rh7', 5900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'exonerate_output_db_server'}, $self->default_options->{'dna_db_server'}]) },
    'exonerate_2' => { LSF => $self->lsf_resource_builder('production-rh7', 5900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'exonerate_output_db_server'}, $self->default_options->{'dna_db_server'}]) },
  }
}

1;
