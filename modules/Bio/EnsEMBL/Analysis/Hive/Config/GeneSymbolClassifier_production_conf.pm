=head1 LICENSE

See the NOTICE file distributed with this work for additional information
regarding copyright ownership.

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

=pod

=head1 NAME

Bio::EnsEMBL::Production::Pipeline::PipeConfig::GeneSymbolClassifier_conf

=head1 DESCRIPTION

eHive pipeline to assign gene symbols to the protein coding gene sequences in an Ensembl core database using a neural network gene symbol classifier.

=cut

package Bio::EnsEMBL::Production::Pipeline::PipeConfig::GeneSymbolClassifier_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Production::Pipeline::PipeConfig::Base_conf');

use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Hive::Version 2.5;
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub default_options {
  my ($self) = @_;

  # inherited from the Base pipeline config
  my $gsc_data_directory = $self->o('pipeline_dir');

  my $core_db_name = $self->o('core_db_name');

  my $pipeline_name => $self-o('pipeline_name') || 'GeneSymbolClassifier_production';
  my $pipe_db_name = $ENV{USER}.$self->o('pipeline_name');

  my $gsc_scripts_directory = $ENV{ENSCODE}.'/ensembl-genes/pipelines/gene_symbol_classifier';

  my $protein_sequences_fasta_path = "${gsc_data_directory}/${core_db_name}_protein_sequences.fa";
  my $gene_symbols_csv_path = "${gsc_data_directory}/${core_db_name}_protein_sequences_symbols.csv";
  my $filtered_assignments_csv_path = "${gsc_data_directory}/${core_db_name}_protein_sequences_symbols_filtered.csv";

  return {
    %{$self->SUPER::default_options},

    species      => [],
    division     => [],
    run_all      => 0,
    antispecies  => [],
    meta_filters => {},

    history_file => undef,

    pipeline_name => $self->o('pipeline_name'),

    pipeline_db => {
      -driver => 'mysql',
      -host   => $self->o('host'),
      -port   => $self->o('port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -dbname => $self->o('pipe_db_name'),
    },

    user_r => 'ensro',

    loading_threshold => $self->o('loading_threshold'),

    gsc_scripts_directory => $gsc_scripts_directory,
    gsc_data_directory => $gsc_data_directory,
    protein_sequences_fasta_path => $protein_sequences_fasta_path,
    gene_symbols_csv_path => $gene_symbols_csv_path,
    filtered_assignments_csv_path => $filtered_assignments_csv_path,
  };
}

# Implicit parameter propagation throughout the pipeline.
sub hive_meta_table {
  my ($self) = @_;

  return {
    %{$self->SUPER::hive_meta_table},
    'hive_use_param_stack'  => 1,
  };
}

sub pipeline_create_commands {
    my ($self) = @_;

    my $gsc_data_directory = $self->o('gsc_data_directory');

    return [
        @{ $self->SUPER::pipeline_create_commands },

        "mkdir --parents --verbose $gsc_data_directory",
    ];
}

sub pipeline_analyses {
  my ($self) = @_;

  return [
    {
      -logic_name      => 'InitialisePipeline',
      -module          => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -input_ids       => [ {} ],
      -max_retry_count => 1,
      -flow_into       => {
        '1->A' => ['SpeciesFactory_All'],
        'A->1' => ['Notify'],
      },
      -rc_name         => 'normal',
    },

    {
      -logic_name      => 'SpeciesFactory_All',
      -module          => 'Bio::EnsEMBL::Production::Pipeline::Common::SpeciesFactory',
      -parameters      => {
        species      => $self->o('species'),
        division     => $self->o('division'),
        run_all      => $self->o('run_all'),
        antispecies  => $self->o('antispecies'),
        meta_filters => $self->o('meta_filters'),
      },
      -max_retry_count => 1,
      -flow_into       => {
        '2' => ['dump_protein_sequences'],
      },
      -rc_name         => 'default',
    },

    {
      # input: Ensembl core db
      # output: FASTA file with protein coding genes canonical sequences
      -logic_name => 'dump_protein_sequences',
      -comment    => 'Retrieve the protein coding gene sequences from a Ensembl core database and store them as a FASTA file. The analysis is auto-seeded with a job for the target core database.',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        'cmd' => 'perl '.$self->o('gsc_scripts_directory').'/dump_protein_sequences.pl --db_host '.$self->o('core_db_server_host').' --db_port '.$self->o('core_db_server_port').' --db_name '.$self->o('core_db_name').' --output_file '.$self->o('protein_sequences_fasta_path'),
      },
      -rc_name    => 'default',
      -flow_into  => {
        1 => 'assign_gene_symbols',
      },
    },

    {
      # input: FASTA file with protein coding genes canonical sequences
      # output: CSV file with symbol assignments and metadata
      -logic_name => 'assign_gene_symbols',
      -comment    => 'Use a gene symbol classifier neural network to assign gene symbols to protein sequences in the FASTA file and save the assignments to a CSV file.',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        'cmd' => 'singularity run --bind '.$self->o('classifier_directory').':/app/checkpoints --bind '.$self->o('gsc_data_directory').':/app/data '.$self->o('singularity_image').' --checkpoint /app/checkpoints/'.$self->o('classifier_filename').' --sequences_fasta /app/data/'.$self->o('core_db_name').'_protein_sequences.fa --scientific_name "'.$self->o('scientific_name').'"',
      },
      -rc_name    => 'default',
      -flow_into  => {
        1 => 'filter_assignments',
      },
    },

    {
      # input: CSV file with symbol assignments and metadata
      # output: CSV file with assignments to be loaded to the Ensembl core db
      -logic_name => 'filter_assignments',
      -comment    => 'Filter assignments using the threshold probability and save them to a separate CSV file.',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        'cmd' => 'singularity run --bind '.$self->o('gsc_data_directory').':/app/data '.$self->o('ehive_singularity_image').' --symbol_assignments '.$self->o('gene_symbols_csv_path').' --threshold '.$self->o('loading_threshold'),
      },
      -rc_name    => 'default',
      -flow_into  => {
        1 => 'load_gene_symbols',
      },
    },

    {
      # input: CSV file with assignments to be loaded to the Ensembl core db
      # output: gene symbols loaded to the Ensembl core db
      -logic_name => 'load_gene_symbols',
      -comment    => 'Read gene symbols assignments from a CSV file and load them to the Ensembl core database.',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        'cmd' => 'perl '.$self->o('gsc_scripts_directory').'/load_gene_symbols.pl --db_host '.$self->o('core_db_server_host').' --db_port '.$self->o('core_db_server_port').' --db_name '.$self->o('core_db_name').' --username '.$self->o('user').' --password '.$self->o('password').' --symbol_assignments '.$self->o('filtered_assignments_csv_path'),
      },
      -rc_name    => 'default',
      -flow_into  => {
        1 => 'RunDataChecks',
      },
    },

    {
      -logic_name      => 'RunDataChecks',
      -module          => 'Bio::EnsEMBL::DataCheck::Pipeline::RunDataChecks',
      -parameters      => {
        datacheck_names => [],
        history_file    => $self->o('history_file'),
        failures_fatal  => 1,
      },
      -max_retry_count => 1,
      -hive_capacity   => 50,
      -batch_size      => 10,
      -rc_name         => 'default',
    },

    {
      -logic_name => 'Notify',
      -module     => 'Bio::EnsEMBL::Production::Pipeline::Production::EmailSummaryCore',
      -parameters => {
        email   => $self->o('email'),
        subject => $self->o('pipeline_name').' has finished',
      },
      -rc_name    => 'normal',
    },

  ];
}

sub resource_classes {
  my ($self) = @_;

  return {
    %{$self->SUPER::resource_classes},
    'default' => { 'LSF' => '-q production'},
    'normal'  => { 'LSF' => '-q production -M 500 -R "rusage[mem=500]"'},
  }
}

1;
