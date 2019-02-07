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

package Bio::EnsEMBL::Analysis::Hive::Config::ExonerateIMGT_conf;

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
    pipeline_name => '',

    user     => '',
    password => '',
    port     => 4533,
    host     => '',
    user_r   => '',
    species_name => 'salmo_salar',
    release_number => '96',

    genome_file          => '/hps/nobackup2/production/ensembl/genebuild/production/fish/salmo_salar/GCA_000233375.4/genome_dumps/salmo_salar_softmasked_toplevel.fa',
    output_dir           => '/hps/nobackup2/production/ensembl/genebuild/production/fish/salmo_salar/GCA_000233375.4/imgt',
    repeatmasker_library => 'teleost',

    exonerate_path                       => catfile($self->o('software_base_path'), 'opt', 'exonerate09', 'bin', 'exonerate'),
    exonerate_calculate_coverage_and_pid => 1,

    imgt_table_name     => 'protein_sequences',

    imgt_db_name    => $self->o('dbowner').'_'.$self->o('species_name').'_imgt_'.$self->o('release_number'),
    imgt_db_server => 'mysql-ens-genebuild-prod-3',
    imgt_db_port   => 4529,

    dna_db_name    => 'fish1_salmo_salar_core_95',
    dna_db_server => 'mysql-ens-genebuild-prod-2',
    dna_db_port   => 4528,

    imgt_db => {
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
        1 => ['fetch_data_file'],
      },
      -input_ids  => [{}],
    },
    {
      -logic_name => 'fetch_data_file',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'wget -P #output_dir# "ftp://ftp.ebi.ac.uk/pub/databases/imgt/LIGM-DB/imgt.dat.Z"',
        output_dir => $self->o('output_dir'),
      },
      -rc_name => 'default',
      -flow_into => {
        1 => ['unzip_data_file'],
      },
    },
    {
      -logic_name => 'unzip_data_file',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'gunzip #output_dir#/imgt.dat.Z',
        output_dir => $self->o('output_dir'),
      },
      -rc_name => 'default',
      -flow_into => {
        1 => ['load_seqs'],
      },
    },

    {
      -logic_name => 'load_seqs',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadIMGT',
      -parameters => {
        sequence_table_name => $self->o('imgt_table_name'),
        iid => catfile($self->o('output_dir'), 'imgt.dat'),
      },
      -rc_name => 'default',
      -flow_into => {
        1 => [
          'generate_ig_c_jobs',
          'generate_ig_d_jobs',
          'generate_ig_j_jobs',
          'generate_ig_v_jobs',
          'generate_full_gene_jobs',
        ],
      },
    },

    {
      -logic_name => 'generate_ig_c_jobs',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        inputquery => 'SELECT accession FROM #sequence_table_name# WHERE biotype IN ("IG_#ig_type#_gene", "TR_#ig_type#_gene", "#ig_type#_gene")',
        sequence_table_name => $self->o('imgt_table_name'),
        column_names => ['iid'],
        ig_type => 'C',
      },
      -rc_name      => 'default',
      -flow_into => {
        2 => ['exonerate_ig_c'],
      },
    },

    {
      -logic_name => 'exonerate_ig_c',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
      -rc_name    => '4GB',
      -parameters => {
        exonerate_path => $self->o('exonerate_path'),
        iid_type => 'db_seq',
        sequence_table_name => $self->o('imgt_table_name'),
        dna_db => $self->o('dna_db'),
        target_db => $self->o('imgt_db'),
        calculate_coverage_and_pid => $self->o('exonerate_calculate_coverage_and_pid'),
        genome_file => $self->o('genome_file'),
        exonerate_path => $self->o('exonerate_path'),
        repeat_libraries => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
        exonerate_bestn => 5,
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','c_segment')},
      },
    },

    {
      -logic_name => 'generate_ig_d_jobs',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        inputquery => 'SELECT accession FROM #sequence_table_name# WHERE biotype IN ("IG_#ig_type#_gene", "TR_#ig_type#_gene", "#ig_type#_gene")',
        sequence_table_name => $self->o('imgt_table_name'),
        column_names => ['iid'],
        ig_type => 'D',
      },
      -rc_name      => 'default',
      -flow_into => {
        2 => ['exonerate_ig_d'],
      },
    },

    {
      -logic_name => 'exonerate_ig_d',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
      -rc_name    => '4GB',
      -parameters => {
        exonerate_path => $self->o('exonerate_path'),
        iid_type => 'db_seq',
        sequence_table_name => $self->o('imgt_table_name'),
        dna_db => $self->o('dna_db'),
        target_db => $self->o('imgt_db'),
        calculate_coverage_and_pid => $self->o('exonerate_calculate_coverage_and_pid'),
        genome_file => $self->o('genome_file'),
        exonerate_path => $self->o('exonerate_path'),
        repeat_libraries => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
        exonerate_bestn => 5,
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','d_segment')},
      },
    },

    {
      -logic_name => 'generate_ig_j_jobs',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        inputquery => 'SELECT accession FROM #sequence_table_name# WHERE biotype IN ("IG_#ig_type#_gene", "TR_#ig_type#_gene", "#ig_type#_gene")',
        sequence_table_name => $self->o('imgt_table_name'),
        column_names => ['iid'],
        ig_type => 'J',
      },
      -rc_name      => 'default',
      -flow_into => {
        2 => ['exonerate_ig_j'],
      },
    },

    {
      -logic_name => 'exonerate_ig_j',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
      -rc_name    => '4GB',
      -parameters => {
        exonerate_path => $self->o('exonerate_path'),
        iid_type => 'db_seq',
        sequence_table_name => $self->o('imgt_table_name'),
        dna_db => $self->o('dna_db'),
        target_db => $self->o('imgt_db'),
        calculate_coverage_and_pid => $self->o('exonerate_calculate_coverage_and_pid'),
        genome_file => $self->o('genome_file'),
        exonerate_path => $self->o('exonerate_path'),
        repeat_libraries => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
        exonerate_bestn => 5,
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','j_segment')},
      },
    },

    {
      -logic_name => 'generate_ig_v_jobs',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        inputquery => 'SELECT accession FROM #sequence_table_name# WHERE biotype IN ("IG_#ig_type#_gene", "TR_#ig_type#_gene", "#ig_type#_gene")',
        sequence_table_name => $self->o('imgt_table_name'),
        column_names => ['iid'],
        ig_type => 'V',
      },
      -rc_name      => 'default',
      -flow_into => {
        2 => ['exonerate_ig_v'],
      },
    },

    {
      -logic_name => 'exonerate_ig_v',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
      -rc_name    => '4GB',
      -parameters => {
        exonerate_path => $self->o('exonerate_path'),
        iid_type => 'db_seq',
        sequence_table_name => $self->o('imgt_table_name'),
        dna_db => $self->o('dna_db'),
        target_db => $self->o('imgt_db'),
        calculate_coverage_and_pid => $self->o('exonerate_calculate_coverage_and_pid'),
        genome_file => $self->o('genome_file'),
        exonerate_path => $self->o('exonerate_path'),
        repeat_libraries => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
        exonerate_bestn => 5,
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','v_segment')},
      },
    },

    {
      -logic_name => 'generate_full_gene_jobs',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        inputquery => 'SELECT accession FROM #sequence_table_name# WHERE biotype IN ("IG_#ig_type#_gene", "#ig_type#_gene")',
        sequence_table_name => $self->o('imgt_table_name'),
        column_names => ['iid'],
        ig_type => 'full',
      },
      -rc_name      => 'default',
      -flow_into => {
        2 => ['exonerate_ig_full'],
      },
    },

    {
      -logic_name => 'exonerate_ig_full',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
      -rc_name    => '4GB',
      -parameters => {
        exonerate_path => $self->o('exonerate_path'),
        iid_type => 'db_seq',
        sequence_table_name => $self->o('imgt_table_name'),
        dna_db => $self->o('dna_db'),
        target_db => $self->o('imgt_db'),
        calculate_coverage_and_pid => $self->o('exonerate_calculate_coverage_and_pid'),
        genome_file => $self->o('genome_file'),
        exonerate_path => $self->o('exonerate_path'),
        repeat_libraries => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
        exonerate_bestn => 5,
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','c_segment')},
      },
    },
  ];
}


sub resource_classes {
  my $self = shift;

  return {
    'default' => { 'LSF' => $self->lsf_resource_builder('production-rh7', 900) },
    '4GB' => { 'LSF' => $self->lsf_resource_builder('production-rh7', 4000) },
  }
}

1;
