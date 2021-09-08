
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

=cut

package Bio::EnsEMBL::Analysis::Hive::Config::GenblastHomology;

use strict;
use warnings;
use File::Spec::Functions;

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_analysis_settings);
use base ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

sub default_options {
  my ($self) = @_;
  return {

    # inherit other stuff from the base class
    %{ $self->SUPER::default_options() },

######################################################
    #
    # Variable settings- You change these!!!
    #
######################################################
########################
    # Misc setup info
########################
    'dbowner' => '' || $ENV{EHIVE_USER} || $ENV{USER},
    'pipeline_name' => '' || $self->o('production_name') . '_' . $self->o('ensembl_release'),
    'user_r'           => '',                                  # read only db user
    'user'             => '',                                  # write db user
    'password'         => '',                                  # password for write db user
    'server_set'       => '',                                  # What server set to user, e.g. set1
    'pipe_db_host'     => '',                                  # host for pipe db
    'databases_host'   => '',                                  # host for general output dbs
    'dna_db_host'      => '',                                  # host for dna db
    'pipe_db_port'     => '',                                  # port for pipeline host
    'databases_port'   => '',                                  # port for general output db host
    'dna_db_port'      => '',                                  # port for dna db host
    'release_number'   => '' || $self->o('ensembl_release'),
    'species_name'     => '',                                  # e.g. mus_musculus
    'production_name'  => '',                                  # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
    'taxon_id'         => '',                                  # should be in the assembly report file
    'uniprot_set'      => '',                                  # e.g. mammals_basic, check UniProtCladeDownloadStatic.pm module in hive config dir for suitable set,
    'output_path'      => '',                                  # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses

########################
    # Pipe and ref db info
########################

# The following might not be known in advance, since the come from other pipelines
# These values can be replaced in the analysis_base table if they're not known yet
# If they are not needed (i.e. no projection or rnaseq) then leave them as is

    'pipe_db_name' => $self->o('dbowner') . '_' . $self->o('production_name') . '_pipe_' . $self->o('release_number'),
    'dna_db_name'  => $self->o('dbowner') . '_' . $self->o('production_name') . '_core_' . $self->o('release_number'),

    'genblast_db_host'   => $self->o('databases_host'),
    'genblast_db_port'   => $self->o('databases_port'),

    'genblast_nr_db_host'   => $self->o('databases_host'),
    'genblast_nr_db_port'   => $self->o('databases_port'),

    'killlist_db_host'   => $self->o('databases_host'),
    'killlist_db_port'   => $self->o('databases_port'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    ensembl_release => $ENV{ENSEMBL_RELEASE},    # this is the current release version on staging to be able to get the correct database

    databases_to_delete => [ 'genblast_db', 'genblast_nr_db' ],

######################################################
    #
    # Mostly constant settings
    #
######################################################

    genome_dumps => catdir( $self->o('output_path'), 'genome_dumps' ),

# This one is used by most analyses that run against a genome flatfile like exonerate, genblast etc. Has slice name style headers. Is softmasked
    softmasked_genome_file => catfile( $self->o('genome_dumps'), $self->o('species_name') . '_softmasked_toplevel.fa' ),

    use_genome_flatfile           => '1',# This will read sequence where possible from a dumped flatfile instead of the core db
    genome_dumps                  => catdir($self->o('output_path'), 'genome_dumps'),
    faidx_genome_file             => catfile($self->o('genome_dumps'), $self->o('species_name').'_toplevel.fa'),

    'min_toplevel_slice_length' => 250,

    'homology_models_path' => catdir( $self->o('output_path'), 'homology_models' ),

    # This is used for "messaging" other sub pipeline
    transcript_selection_url => undef,

########################
    # Executable paths
########################

    'blast_type'                  => 'ncbi',                                                                                # It can be 'ncbi', 'wu', or 'legacy_ncbi'
    'exonerate_path'              => catfile( $self->o('linuxbrew_home_path'), 'opt', 'exonerate09', 'bin', 'exonerate' ),
    'uniprot_genblast_batch_size' => 15,
    'uniprot_table_name'          => 'uniprot_sequences',
    genewise_path   => catfile( $self->o('binary_base'), 'genewise' ),
    'genblast_path' => catfile( $self->o('binary_base'), 'genblast' ),
    'genblast_eval' => $self->o('blast_type') eq 'wu' ? '1e-20' : '1e-1',
    'genblast_cov'  => '0.5',
    'genblast_pid'  => '30',
    'genblast_max_rank'           => '5',
    'genblast_flag_small_introns' => 1,
    'genblast_flag_subpar_models' => 0,

########################
    # db info
########################

    'genblast_db' => {
      -dbname => $self->o('dbowner') . '_' . $self->o('production_name') . '_genblast_' . $self->o('release_number'),
      -host   => $self->o('genblast_db_host'),
      -port   => $self->o('genblast_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'genblast_nr_db' => {
      -dbname => $self->o('dbowner') . '_' . $self->o('production_name') . '_genblast_nr_' . $self->o('release_number'),
      -host   => $self->o('genblast_nr_db_host'),
      -port   => $self->o('genblast_nr_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'killlist_db' => {
      -dbname => $self->o('killlist_db_name'),
      -host   => $self->o('killlist_db_host'),
      -port   => $self->o('killlist_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -driver => $self->o('hive_driver'),
    },

  };
}

sub pipeline_create_commands {
  my ($self) = @_;

  return [

    # inheriting database and hive tables' creation
    @{ $self->SUPER::pipeline_create_commands },

    $self->hive_data_table( 'protein', $self->o('uniprot_table_name') ),

  ];
}


sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
    use_genome_flatfile => $self->o('use_genome_flatfile'),
    genome_file => $self->o('faidx_genome_file'),
  }
}


## See diagram for pipeline structure
sub pipeline_analyses {
  my ($self) = @_;
  my %genblast_params = (
    wu_genome   => '-P wublast -gff -e #blast_eval# -c #blast_cov#',
    ncbi_genome => '-P blast -gff -e #blast_eval# -c #blast_cov# -W 3 -softmask -scodon 50 -i 30 -x 10 -n 30 -d 200000 -g T',
  );
  my %commandline_params = (
    'ncbi'        => '-num_threads 3 -window_size 40',
    'wu'          => '-cpus 3 -hitdist 40',
    'legacy_ncbi' => '-a 3 -A 40',
  );
  return [
######################################################################################
    #
    # Protein models (genblast and genewise)
    #
######################################################################################
    {
      -logic_name => 'create_genblast_output_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('dna_db'),
        target_db   => $self->o('genblast_db'),
        create_type => 'clone',
      },
      -rc_name   => 'default',
      -input_ids  => [{}],
      -flow_into => {
        '1->A' => ['download_uniprot_files'],
        'A->1' => ['classify_genblast_models'],
      },
    },

    {
      -logic_name => 'download_uniprot_files',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadUniProtFiles',
      -parameters => {
        multi_query_download => get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::UniProtCladeDownloadStatic', $self->o('uniprot_set') ),
        taxon_id             => $self->o('taxon_id'),
        output_path          => $self->o('homology_models_path'),
      },
      -rc_name   => 'default',
      -flow_into => {
        '2->A' => ['process_uniprot_files'],
        'A->1' => ['generate_genblast_jobs'],
      },
    },

    {
      -logic_name => 'process_uniprot_files',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProcessUniProtFiles',
      -parameters => {
        killlist_type       => 'protein',
        killlist_db         => $self->o('killlist_db'),
        sequence_table_name => $self->o('uniprot_table_name'),
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'generate_genblast_jobs',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        iid_type            => 'sequence_accession',
        batch_size          => $self->o('uniprot_genblast_batch_size'),
        sequence_table_name => $self->o('uniprot_table_name'),
      },
      -rc_name   => '3GB',
      -flow_into => {
        2 => ['genblast'],
        1 => ['create_seleno_homology_jobs'],
      },
    },

    {
      -logic_name => 'create_seleno_homology_jobs',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        inputquery   => 'SELECT accession FROM ' . $self->o('uniprot_table_name') . ' WHERE source_db = "seleno"',
        column_names => ['iid'],
      },
      -rc_name   => 'default',
      -flow_into => {
        2 => ['process_homology_selenocysteine'],
      },
    },

    {
      -logic_name => 'process_homology_selenocysteine',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSelenocysteineFinder',
      -parameters => {
        target_db           => $self->o('genblast_db'),
        dna_db              => $self->o('dna_db'),
        genome              => $self->o('softmasked_genome_file'),
        exonerate           => $self->o('exonerate_path'),
        genewise            => $self->o('genewise_path'),
        iid_type            => 'db_seq',
        sequence_table_name => $self->o('uniprot_table_name'),
        biotype             => 'seleno_other',
        missmatch_allowed   => 10,
      },
      -rc_name => '3GB',
    },

    {
      -logic_name => 'genblast',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
      -parameters => {
        iid_type            => 'db_seq',
        dna_db              => $self->o('dna_db'),
        target_db           => $self->o('genblast_db'),
        logic_name          => 'genblast',
        module              => 'HiveGenblast',
        genblast_path       => $self->o('genblast_path'),
        genblast_db_path    => $self->o('softmasked_genome_file'),
        commandline_params  => $genblast_params{ $self->o('blast_type') . '_genome' },
        sequence_table_name => $self->o('uniprot_table_name'),
        max_rank            => $self->o('genblast_max_rank'),
        genblast_pid        => $self->o('genblast_pid'),
        timer               => '2h',
        blast_eval          => $self->o('genblast_eval'),
        blast_cov           => $self->o('genblast_cov'),
        flag_small_introns  => $self->o('genblast_flag_small_introns'),
        flag_subpar_models  => $self->o('genblast_flag_subpar_models'),
      },
      -rc_name   => 'genblast',
      -flow_into => {
        -1 => ['split_genblast_jobs'],
        -2 => ['split_genblast_jobs'],
        -3 => ['split_genblast_jobs'],
      },
      -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
    },

    {
      -logic_name => 'split_genblast_jobs',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        iid_type   => 'rechunk',
        batch_size => 1,
      },
      -rc_name      => 'default',
      -can_be_empty => 1,
      -flow_into    => {
        2 => ['genblast_retry'],
      },
    },

    {
      -logic_name => 'genblast_retry',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
      -parameters => {
        iid_type            => 'db_seq',
        dna_db              => $self->o('dna_db'),
        target_db           => $self->o('genblast_db'),
        logic_name          => 'genblast',
        module              => 'HiveGenblast',
        genblast_path       => $self->o('genblast_path'),
        genblast_db_path    => $self->o('softmasked_genome_file'),
        commandline_params  => $genblast_params{ $self->o('blast_type') . '_genome' },
        sequence_table_name => $self->o('uniprot_table_name'),
        max_rank            => $self->o('genblast_max_rank'),
        genblast_pid        => $self->o('genblast_pid'),
        timer               => '1h',
        blast_eval          => $self->o('genblast_eval'),
        blast_cov           => $self->o('genblast_cov'),
        flag_small_introns  => $self->o('genblast_flag_small_introns'),
        flag_subpar_models  => $self->o('genblast_flag_subpar_models'),
      },
      -rc_name              => 'genblast_retry',
      -can_be_empty         => 1,
      -failed_job_tolerance => 100,
      -flow_into            => {
        -1 => ['failed_genblast_proteins'],
        -2 => ['failed_genblast_proteins'],
        -3 => ['failed_genblast_proteins'],
      },
      -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
    },

    {
      -logic_name => 'failed_genblast_proteins',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -parameters => {
      },
      -rc_name              => 'default',
      -can_be_empty         => 1,
      -failed_job_tolerance => 100,
    },

    {
      -logic_name => 'classify_genblast_models',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport',
      -parameters => {
        classification_type => 'standard',
        update_gene_biotype => 1,
        target_db           => $self->o('genblast_db'),
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['genblast_sanity_checks'],
      },
    },

    {
      -logic_name => 'genblast_sanity_checks',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
      -parameters => {
        target_db                  => $self->o('genblast_db'),
        sanity_check_type          => 'gene_db_checks',
        min_allowed_feature_counts => get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
          'gene_db_checks' )->{ $self->o('uniprot_set') }->{'genblast'},
      },
      -rc_name   => '4GB',
      -flow_into => {
        1 => ['create_genblast_nr_db'],
      },
    },

    {
      -logic_name => 'create_genblast_nr_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('genblast_db'),
        target_db   => $self->o('genblast_nr_db'),
        create_type => 'copy',
      },
      -rc_name   => 'default',
      -flow_into => {
        '1' => ['create_genblast_nr_slices'],
      },
    },

    {
      -logic_name => 'create_genblast_nr_slices',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db             => $self->o('dna_db'),
        coord_system_name     => 'toplevel',
        iid_type              => 'slice',
        slice_size            => 20000000,
        include_non_reference => 0,
        top_level             => 1,
        min_slice_length      => $self->o('min_toplevel_slice_length'),
        batch_slice_ids       => 1,
        batch_target_size     => 20000000,
      },
      -rc_name   => '2GB',
      -flow_into => {
        '2->A' => ['remove_redundant_genblast_genes'],
        'A->1'  => ['notification_pipeline_is_done'],
      },
    },

    {
      -logic_name => 'remove_redundant_genblast_genes',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::RemoveRedundantGenes',
      -parameters => {
        target_db   => $self->o('genblast_nr_db'),
        target_type => 'biotype_priority',
        layers      => get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::LayerAnnotationStatic', $self->o('uniprot_set'), undef, 'ARRAY' ),
      },
      -rc_name => '5GB',
    },

    {
      -logic_name => 'notification_pipeline_is_done',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::MessagePipeline',
      -parameters => {
        messages   => [
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'create_toplevel_slices',
          param => 'feature_dbs',
          data => $self->o('genblast_nr_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'split_slices_on_intergenic',
          param => 'input_gene_dbs',
          data => $self->o('genblast_nr_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'layer_annotation',
          param => 'SOURCEDB_REFS',
          data => $self->o('genblast_nr_db'),
          update => 1,
        }],
        tweak_script => catfile($self->o('enscode_root_dir'), 'ensembl-hive', 'scripts', 'tweak_pipeline.pl'),
      },
      -rc_name    => 'default',
    },

  ];
}

sub resource_classes {
  my $self = shift;
  return {
    '2GB'     => { LSF => $self->lsf_resource_builder( 'production', 2000 ) },
    '3GB'     => { LSF => $self->lsf_resource_builder( 'production', 3000 ) },
    '4GB'     => { LSF => $self->lsf_resource_builder( 'production', 4000 ) },
    '5GB'     => { LSF => $self->lsf_resource_builder( 'production', 5000 ) },
    'default' => { LSF => $self->lsf_resource_builder( 'production', 900 ) },
    'genblast'       => { LSF => $self->lsf_resource_builder( 'production', 3900 ) },
    'genblast_retry' => { LSF => $self->lsf_resource_builder( 'production', 4900 ) },
    }
}

1;
