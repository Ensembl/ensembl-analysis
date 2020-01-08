=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

ExonerateIMGT_conf

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package ExonerateIMGT_conf;

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
    'pipeline_name' => '',

    'user'     => '',
    'password' => '',
    'port'     => 4533,
    'host'     => '',
    'user_r'   => '',

    'genome_file'          => '',
    'fasta_file'           => '',
    'repeatmasker_library' => '',

    'exonerate_max_intron'                 => '5000', # Max intron size, default should be 200000
    'exonerate_pid'                        => '50', # Cut-off for percent id
    'exonerate_cov'                        => '50', # Cut-off for coverage
    'exonerate_path'                       => catfile($self->o('software_base_path'), 'opt', 'exonerate09', 'bin', 'exonerate'),
    'exonerate_calculate_coverage_and_pid' => 1,

    'imgt_table_name'     => 'protein_sequences',
    'sequence_batch_size' => 1,

    'imgt_db_name'    => $self->o('dbowner').'_'.$self->o('pipeline_name').'_imgt',
    'imgt_db_server' => '',
    'imgt_db_port'   => 4529,

    'dna_db_name'    => '',
    'dna_db_server' => '',
    'dna_db_port'   => '',

    'imgt_db' => {
      -dbname => $self->o('imgt_db_name'),
      -host   => $self->o('imgt_db_server'),
      -port   => $self->o('imgt_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },
    databases_to_delete => ['imgt_db'],
  };
}

sub pipeline_create_commands {
    my ($self) = @_;
    return [
    # inheriting database and hive tables' creation
      @{$self->SUPER::pipeline_create_commands},
      $self->hive_data_table('protein', $self->o('imgt_table_name')),
    ];
}
sub pipeline_analyses {
  my ($self) = @_;

  my %commandline_params = (
      'ncbi' => '-num_threads 3 -window_size 40',
      'wu' => '-cpus 3 -hitdist 40',
      'legacy_ncbi' => '-a 3 -A 40',
      );

  return [
    {
      -logic_name => 'create_output_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db => $self->o('dna_db'),
        target_db => $self->o('imgt_db'),
        create_type => 'clone',
      },
      -rc_name    => 'default',
      -flow_into => {
        1 => ['load_seqs'],
      },
      -input_ids  => [{iid => $self->o('fasta_file')}],
    },

    {
      -logic_name => 'load_seqs',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadIMGT',
      -parameters => {
        sequence_table_name => $self->o('imgt_table_name'),
      },
      -rc_name => 'default',
      -flow_into => {
        1 => ['generate_jobs'],
      },
    },

    {
      -logic_name => 'generate_jobs',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        iid_type => 'sequence_accession',
        batch_size => $self->o('sequence_batch_size'),
        sequence_table_name => $self->o('imgt_table_name'),
      },
      -rc_name      => 'default',
      -flow_into => {
        2 => ['exonerate'],
      },
    },

    {
      -logic_name => 'exonerate',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
      -rc_name    => 'align4GB',
      -parameters => {
        exonerate_path => $self->o('exonerate_path'),
        iid_type => 'db_seq',
        sequence_table_name => $self->o('imgt_table_name'),
        dna_db => $self->o('dna_db'),
        target_db => $self->o('imgt_db'),
        logic_name => 'exonerate',
        module     => 'HiveExonerate2Genes',
        calculate_coverage_and_pid => $self->o('exonerate_calculate_coverage_and_pid'),
        exonerate_cdna_pid => $self->o('exonerate_pid'),
        exonerate_cdna_cov => $self->o('exonerate_cov'),
        exonerate_max_intron => $self->o('exonerate_max_intron'),
        genome_file => $self->o('genome_file'),
        exonerate_path => $self->o('exonerate_path'),
        repeat_libraries => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
        exonerate_bestn => 1,
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','protein_cov_per_bestn_maxintron_sub')},
      },
    },
  ];
}


sub resource_classes {
  my $self = shift;

  return {
    'default' => { 'LSF' => $self->lsf_resource_builder('production-rh7', 900) },
    'align4GB' => { 'LSF' => $self->lsf_resource_builder('production-rh7', 4000) },
  }
}
1;
