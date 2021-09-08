
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Analysis::Hive::Config::LastZ;

use strict;
use warnings;
use File::Spec::Functions;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
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
    'user_r'         => '',    # read only db user
    'user'           => '',    # write db user
    'password'       => '',    # password for write db user
    'pipe_db_host'   => '',    # host for pipe db
    'dna_db_host'    => '',    # host for dna db
    'pipe_db_port'   => '',    # port for pipeline host
    'dna_db_port'    => '',    # port for dna db host

    'release_number' => '' || $self->o('ensembl_release'),
    'species_name'    => '',                                                          # e.g. mus_musculus
    'production_name' => '',                                                          # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
    'output_path'     => '',                                                          # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses
    'registry_file'   => '' || catfile( $self->o('output_path'), "Databases.pm" ),    # Path to databse registry for LastaZ and Production sync

########################
    # Pipe and ref db info
########################

    'projection_source_db_name'         => '',                                        # This is generally a pre-existing db, like the current human/mouse core for example
    'projection_source_db_host'         => 'mysql-ens-mirror-1',
    'projection_source_db_port'         => '4240',
    'projection_source_production_name' => '',

    'compara_db_name'   => 'leanne_ensembl_compara_95',
    'compara_db_host'   => 'mysql-ens-genebuild-prod-5',
    'compara_db_port'   => 4531,

    'pipe_db_name' => $self->o('dbowner') . '_' . $self->o('production_name') . '_pipe_' . $self->o('release_number'),
    'dna_db_name'  => $self->o('dbowner') . '_' . $self->o('production_name') . '_core_' . $self->o('release_number'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'ensembl_release' => $ENV{ENSEMBL_RELEASE},    # this is the current release version on staging to be able to get the correct database

########################
    # Extra db settings
########################

    'compara_master'                    => 'compara_master',
    'compara_conf_file'                 => '',
    'compara_innodb_schema'             => 1,
    'compara_genome_db_update_path'     => catfile( $self->o('enscode_root_dir'), '/ensembl-compara/scripts/pipeline/update_genome.pl' ),
    'compara_mlss_script_path'          => catfile( $self->o('enscode_root_dir'), '/ensembl-compara/scripts/pipeline/create_mlss.pl' ),
    'compara_mlss_reg_conf_path'        => catfile( $self->o('enscode_root_dir'), '/ensembl-compara/scripts/pipeline/production_reg_ensembl_conf.pl' ),
    'compara_populate_new_database_exe' => catfile( $self->o('enscode_root_dir'), 'ensembl-compara/scripts/pipeline/populate_new_database.pl' ),
    'compara_only_cellular_component'   => undef,
    'compara_dump_dir'                  => catdir( $self->o('output_path'), 'lastz' ),

    'mlss_id_list'       => undef,
    'compara_collection' => '',

    'compara_ref_species'     => $self->o('projection_source_production_name'),
    'compara_non_ref_species' => $self->o('production_name'),
    'only_cellular_component' => undef,                                           # Do we load *all* the dnafrags or only the ones from a specific cellular-component ?
    'mix_cellular_components' => 0,                                               # Do we try to allow the nuclear genome vs MT, etc ?
    'dump_min_nib_size'       => 11500000,
    'dump_min_chunk_size'     => 1000000,
    'dump_min_chunkset_size'  => 1000000,
    'quick'                   => 1,
    'default_chunks'          => {
      'reference' => {
        'homo_sapiens' => {
          'chunk_size'            => 30000000,
          'overlap'               => 0,
          'include_non_reference' => -1,                                          #1  => include non_reference regions (eg human assembly patches)
                                                                                  #0  => do not include non_reference regions
                                                                                  #-1 => auto-detect (only include non_reference regions if the non-reference species is high-coverage
                                                                                  #ie has chromosomes since these analyses are the only ones we keep up-to-date with the patches-pipeline)
          'masking_options'       => '{default_soft_masking => 1}',

# if you have a specific selection of repeat elements for the masking
#'masking_options_file' => $self->check_file_in_ensembl('ensembl-compara/scripts/pipeline/human36.spec'),
        },

        #non human example
        'default' => {
          'chunk_size'      => 10000000,
          'overlap'         => 0,
          'masking_options' => '{default_soft_masking => 1}'
        },
      },
      'non_reference' => {
        'chunk_size'      => 10100000,
        'group_set_size'  => 10100000,
        'overlap'         => 100000,
        'masking_options' => '{default_soft_masking => 1}'
      },
    },

    'compara_window_size'             => 10000,
    'filter_duplicates_rc_name'       => '2GB_lastz',
    'filter_duplicates_himem_rc_name' => '8GB_lastz',

    #
    #Default pair_aligner
    #
    'pair_aligner_method_link' => [ 1001, 'LASTZ_RAW' ],
    'pair_aligner_logic_name'  => 'LastZ',
    'pair_aligner_module'      => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::LastZ',

    'pair_aligner_options' => {
      default => 'T=1 L=3000 H=2200 O=400 E=30 --ambiguous=iupac',                                                                                                          # ensembl genomes settings
      7742    => 'T=1 K=3000 L=3000 H=2200 O=400 E=30 --ambiguous=iupac',                                                                                                   # vertebrates - i.e. ensembl-specific
      9526    => 'T=1 K=5000 L=5000 H=3000 M=10 O=400 E=30 Q=' . $self->check_file_in_ensembl('ensembl-compara/scripts/pipeline/primate.matrix') . ' --ambiguous=iupac',    # primates
      33554   => 'T=1 K=5000 L=5000 H=3000 M=10 O=400 E=30 --ambiguous=iupac',                                                                                              # carnivora
      3913    => 'T=1 L=3000 H=2200 O=400 E=30 --ambiguous=iupac --matchcount=1000',
      4070    => 'T=1 L=3000 H=2200 O=400 E=30 --ambiguous=iupac --matchcount=1000',
    },

    #
    #Default chain
    #
    'chain_input_method_link'  => [ 1001, 'LASTZ_RAW' ],
    'chain_output_method_link' => [ 1002, 'LASTZ_CHAIN' ],

    #linear_gap=>medium for more closely related species, 'loose' for more distant
    'linear_gap' => 'medium',
    'chain_parameters' => { 'max_gap' => '50', 'linear_gap' => $self->o('linear_gap'), 'faToNib' => $self->o('faToNib_exe'), 'lavToAxt' => $self->o('lavToAxt_exe'), 'axtChain' => $self->o('axtChain_exe'), 'max_blocks_for_chaining' => 100000 },

    #
    #Default patch_alignments
    #
    'patch_alignments' => 0,    #set to 1 to align the patches of a species to many other species

    #
    #Default net
    #
    'net_input_method_link'  => [ 1002, 'LASTZ_CHAIN' ],
    'net_output_method_link' => [ 16,   'LASTZ_NET' ],
    'net_ref_species' => $self->o('compara_ref_species'),                                 #default to ref_species
    'net_parameters'  => { 'max_gap' => '50', 'chainNet' => $self->o('chainNet_exe') },
    'bidirectional'   => 0,

    #
    #Default healthcheck
    #
    'previous_db'               => 'compara_prev',
    'prev_release'              => 0,                # 0 is the default and it means "take current release number and subtract 1"
    'max_percent_diff'          => 20,
    'max_percent_diff_patches'  => 99.99,
    'do_pairwise_gabs'          => 1,
    'do_compare_to_previous_db' => 0,

    'compara_bed_dir'     => $self->o('compara_dump_dir') . '/bed_dir',
    'compara_feature_dir' => $self->o('compara_dump_dir') . '/feature_dumps',

    #
    #Default pairaligner config
    #
    'skip_pairaligner_stats' => 1,    #skip this module if set to 1

    'pair_aligner_method_link' => [ 1001, 'LASTZ_RAW' ],
    'pair_aligner_logic_name'  => 'LastZ',
    'pair_aligner_module'      => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::LastZ',
    'chain_input_method_link'  => [ 1001, 'LASTZ_RAW' ],
    'chain_output_method_link' => [ 1002, 'LASTZ_CHAIN' ],
    'linear_gap'               => 'medium',
    'net_input_method_link'    => [ 1002, 'LASTZ_CHAIN' ],
    'net_output_method_link'   => [ 16,   'LASTZ_NET' ],

    # Capacities
    'pair_aligner_analysis_capacity'  => 700,
    'pair_aligner_batch_size'         => 40,
    'chain_hive_capacity'             => 200,
    'chain_batch_size'                => 10,
    'net_hive_capacity'               => 300,
    'net_batch_size'                  => 10,
    'filter_duplicates_hive_capacity' => 200,
    'filter_duplicates_batch_size'    => 10,

    # LastZ is used to align the genomes
    opt_dir                      => catdir($self->o('linuxbrew_home_path'), 'opt'),
    pair_aligner_exe             => catfile( $self->o('opt_dir'), 'lastz', 'bin', 'lastz' ),
    axtChain_exe                 => catfile( $self->o('opt_dir'), 'kent', 'bin', 'axtChain' ),
    chainNet_exe                 => catfile( $self->o('opt_dir'), 'kent', 'bin', 'chainNet' ),
    faToNib_exe                  => catfile( $self->o('opt_dir'), 'kent', 'bin', 'faToNib' ),
    lavToAxt_exe                 => catfile( $self->o('opt_dir'), 'kent', 'bin', 'lavToAxt' ),
    compara_scripts              => catdir($self->o('enscode_root_dir'), 'ensembl-compara', 'scripts'),
    compare_beds_exe             => catfile( $self->o('enscode_root_dir'), 'pipeline', 'compare_beds.pl' ),
    create_pair_aligner_page_exe => catfile( $self->o('enscode_root_dir'), 'report', 'create_pair_aligner_page.pl' ),
    dump_features_exe            => catfile( $self->o('enscode_root_dir'), 'dumps', 'DumpMultiAlign.pl' ),

########################
    # db info
########################
    'projection_source_db' => {
      -dbname => $self->o('projection_source_db_name'),
      -host   => $self->o('projection_source_db_host'),
      -port   => $self->o('projection_source_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -driver => $self->o('hive_driver'),
    },

    'compara_db' => {
      -dbname => $self->o('compara_db_name'),
      -host   => $self->o('compara_db_host'),
      -port   => $self->o('compara_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },
  };
}

sub pipeline_create_commands {
  my ($self) = @_;

  my $second_pass = $self->_is_second_pass;
  return $self->SUPER::pipeline_create_commands if $self->can('no_compara_schema');
  my $pipeline_url = $self->pipeline_url();
  my $parsed_url   = $second_pass && Bio::EnsEMBL::Hive::Utils::URL::parse($pipeline_url);
  my $driver       = $second_pass ? $parsed_url->{'driver'} : '';

  return [

    # inheriting database and hive tables' creation
    @{ $self->SUPER::pipeline_create_commands },

    'mkdir -p ' . $self->o('compara_dump_dir'),
    'mkdir -p ' . $self->o('compara_bed_dir'),

# Compara 'release' tables will be turned from MyISAM into InnoDB on the fly by default:
    ( $self->o('compara_innodb_schema') ? "sed 's/ENGINE=MyISAM/ENGINE=InnoDB/g' " : 'cat ' )
      . $self->check_file_in_ensembl('ensembl-compara/sql/table.sql') . ' | ' . $self->db_cmd(),

# Compara 'pipeline' tables are already InnoDB, but can be turned to MyISAM if needed:
    ( $self->o('compara_innodb_schema') ? 'cat ' : "sed 's/ENGINE=InnoDB/ENGINE=MyISAM/g' " )
      . $self->check_file_in_ensembl('ensembl-compara/sql/pipeline-tables.sql') . ' | ' . $self->db_cmd(),

    # MySQL specific procedures
    $driver eq 'mysql' ? ( $self->db_cmd() . ' < ' . $self->check_file_in_ensembl( 'ensembl-compara/sql/procedures.' . $driver ) ) : (),
  ];
}

sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{ $self->SUPER::pipeline_wide_parameters },
    }
}

sub pipeline_analyses {
  my ($self) = @_;

  return [
    {
      -logic_name => 'setup_lastz',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::LastZSetup',
      -parameters => {
        compara_genome_db_update_path => $self->o('compara_genome_db_update_path'),
        compara_mlss_script_path      => $self->o('compara_mlss_script_path'),
        compara_mlss_reg_conf_path    => $self->o('compara_mlss_reg_conf_path'),
        compara_db                    => $self->o('compara_db'),
        projection_source_db          => $self->o('projection_source_db'),
        target_db                     => $self->o('dna_db'),
        pipeline_db                   => $self->o('pipeline_db'),
        output_path                   => $self->o('output_path'),
        compara_db_url                => 'mysql://' . $self->o('user') . ':' . $self->o('password') . '@' . $self->o('compara_db_host') . ':' . $self->o('compara_db_port') . '/' . $self->o('compara_db_name'),
        registry_file                 => $self->o('registry_file'),
      },
      -rc_name   => '2GB_lastz',
      -input_ids  => [{}],
      -flow_into => {
        1 => { 'get_species_list' => { 'mlss_id' => '#mlss_id#', 'reg_conf' => '#reg_conf#' } },
      },
    },

    {
      -logic_name => 'get_species_list',
      -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::ParsePairAlignerConf',
      -parameters => {
        'master_db'        => $self->o('compara_master'),
        'get_species_list' => 1,
        'core_dbs'         => undef,
      },
      -rc_name   => '2GB_lastz',
      -flow_into => {
        1 => { 'populate_new_database' => { 'mlss_id' => '#mlss_id#', 'reg_conf' => '#reg_conf#' } },
      },
    },

    {
      -logic_name => 'populate_new_database',
      -module     => 'Bio::EnsEMBL::Compara::RunnableDB::GenomicAlignBlock::PopulateNewDatabase',
      -parameters => {
        'program'                 => $self->o('compara_populate_new_database_exe'),
        'mlss_id_list'            => undef,
        'collection'              => $self->o('compara_collection'),
        'master_db'               => $self->o('compara_master'),
        'only_cellular_component' => $self->o('only_cellular_component'),
      },
      -rc_name   => '2GB_lastz',
      -flow_into => {
        1 => { 'add_method_link_species_link_tag' => { 'mlss_id' => '#mlss_id#', 'reg_conf' => '#reg_conf#' } },
      },
    },

    {
      -logic_name => 'add_method_link_species_link_tag',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('pipeline_db'),
        sql     => [
          'INSERT INTO method_link_species_set_tag VALUES (#mlss_id#,"reference_species","' . $self->o('compara_ref_species') . '")',
        ],
      },
      -rc_name   => '2GB_lastz',
      -flow_into => {
        1 => { 'parse_pair_aligner_conf' => { 'mlss_id' => '#mlss_id#', 'reg_conf' => '#reg_conf#' } }
      },
    },

    {
      -logic_name => 'parse_pair_aligner_conf',
      -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::ParsePairAlignerConf',
      -parameters => {
        'conf_file'                 => $self->o('compara_conf_file'),
        'ref_species'               => $self->o('compara_ref_species'),
        'non_ref_species'           => $self->o('compara_non_ref_species'),
        'dump_dir'                  => $self->o('compara_dump_dir'),
        'default_chunks'            => $self->o('default_chunks'),
        'default_pair_aligner'      => $self->o('pair_aligner_method_link'),
        'default_parameters'        => $self->o('pair_aligner_options'),
        'default_chain_output'      => $self->o('chain_output_method_link'),
        'default_net_output'        => $self->o('net_output_method_link'),
        'default_chain_input'       => $self->o('chain_input_method_link'),
        'default_net_input'         => $self->o('net_input_method_link'),
        'net_ref_species'           => $self->o('net_ref_species'),
        'mlss_id_list'              => $self->o('mlss_id_list'),
        'collection'                => $self->o('compara_collection'),
        'master_db'                 => $self->o('compara_master'),
        'do_pairwise_gabs'          => $self->o('do_pairwise_gabs'),             #healthcheck options
        'do_compare_to_previous_db' => $self->o('do_compare_to_previous_db'),    #healthcheck options
        'bidirectional'             => $self->o('bidirectional'),
      },
      -rc_name   => '2GB_lastz',
      -flow_into => {
        1  => ['create_pair_aligner_jobs'],
        2  => ['chunk_and_group_dna'],
        3  => ['create_filter_duplicates_jobs'],
        4  => ['no_chunk_and_group_dna'],
        5  => ['create_alignment_chains_jobs'],
        6  => ['create_alignment_nets_jobs'],
        10 => ['create_filter_duplicates_net_jobs'],
        7  => ['pairaligner_stats'],
        8  => ['healthcheck'],
      },
    },

    {
      -logic_name => 'chunk_and_group_dna',
      -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::ChunkAndGroupDna',
      -parameters => {
        'only_cellular_component' => $self->o('only_cellular_component'),
        'mix_cellular_components' => $self->o('mix_cellular_components'),
      },
      -rc_name   => '2GB_lastz',
      -flow_into => {
        2 => ['store_sequence'],
      },
    },

    {
      -logic_name    => 'store_sequence',
      -hive_capacity => 100,
      -module        => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::StoreSequence',
      -parameters    => {
        'dump_min_chunkset_size' => $self->o('dump_min_chunkset_size'),
        'dump_min_chunk_size'    => $self->o('dump_min_chunk_size'),
      },
      -rc_name   => '2GB_lastz',
      -flow_into => {
        -1 => ['store_sequence_again'],
      },
    },

    {
      -logic_name    => 'store_sequence_again',
      -hive_capacity => 50,
      -module        => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::StoreSequence',
      -parameters    => {
        'dump_min_chunkset_size' => $self->o('dump_min_chunkset_size'),
        'dump_min_chunk_size'    => $self->o('dump_min_chunk_size'),
      },
      -can_be_empty => 1,
      -rc_name      => '4GB_lastz',
    },

    {
      -logic_name => 'create_pair_aligner_jobs',                                                #factory
      -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::CreatePairAlignerJobs',
      -parameters => {
        'mix_cellular_components' => $self->o('mix_cellular_components'),
      },
      -hive_capacity => 10,
      -wait_for      => [ 'store_sequence', 'store_sequence_again', 'chunk_and_group_dna' ],
      -rc_name       => '2GB_lastz',
      -flow_into     => {
        1 => ['check_no_partial_gabs'],
        2 => [ $self->o('pair_aligner_logic_name') ],
      },
    },

    {
      -logic_name        => $self->o('pair_aligner_logic_name'),
      -module            => $self->o('pair_aligner_module'),
      -analysis_capacity => $self->o('pair_aligner_analysis_capacity'),
      -batch_size        => $self->o('pair_aligner_batch_size'),
      -parameters        => {
        'pair_aligner_exe' => $self->o('pair_aligner_exe'),
      },
      -wait_for  => ['create_pair_aligner_jobs'],
      -rc_name   => '2GB_lastz',
      -flow_into => {
        -1 => [ $self->o('pair_aligner_logic_name') . '_himem1' ],    # MEMLIMIT
      },
    },

    {
      -logic_name        => $self->o('pair_aligner_logic_name') . "_himem1",
      -module            => $self->o('pair_aligner_module'),
      -analysis_capacity => $self->o('pair_aligner_analysis_capacity'),
      -parameters        => {
        'pair_aligner_exe' => $self->o('pair_aligner_exe'),
      },
      -wait_for     => ['create_pair_aligner_jobs'],
      -batch_size   => $self->o('pair_aligner_batch_size'),
      -can_be_empty => 1,
      -rc_name      => '8GB_lastz',
    },

    {
      -logic_name => 'check_no_partial_gabs',
      -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::SqlHealthChecks',
      -parameters => {
        'mode' => 'gab_inconsistencies',
      },
      -wait_for => [ $self->o('pair_aligner_logic_name'), $self->o('pair_aligner_logic_name') . "_himem1" ],
      -flow_into => {
        1 => ['update_max_alignment_length_before_FD'],
      },
    },

    {
      -logic_name => 'update_max_alignment_length_before_FD',
      -module     => 'Bio::EnsEMBL::Compara::RunnableDB::GenomicAlignBlock::UpdateMaxAlignmentLength',
      -parameters => {
        'quick' => $self->o('quick'),
      },
      -rc_name   => '2GB_lastz',
      -flow_into => {
        1 => ['update_max_alignment_length_after_FD'],
      },
    },

    {
      -logic_name => 'create_filter_duplicates_jobs',                                                #factory
      -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::CreateFilterDuplicatesJobs',
      -parameters => {},
      -wait_for   => ['update_max_alignment_length_before_FD'],
      -rc_name    => '2GB_lastz',
      -flow_into  => {
        2 => { 'filter_duplicates' => INPUT_PLUS() },
      },
    },

    {
      -logic_name => 'filter_duplicates',
      -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::FilterDuplicates',
      -parameters => {
        'window_size' => $self->o('compara_window_size')
      },
      -hive_capacity => $self->o('filter_duplicates_hive_capacity'),
      -batch_size    => $self->o('filter_duplicates_batch_size'),
      -can_be_empty  => 1,
      -rc_name       => '2GB_lastz',
      -flow_into     => {
        -1 => ['filter_duplicates_himem'],    # MEMLIMIT
      },
    },

    {
      -logic_name => 'filter_duplicates_himem',
      -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::FilterDuplicates',
      -parameters => {
        'window_size' => $self->o('compara_window_size')
      },
      -hive_capacity => $self->o('filter_duplicates_hive_capacity'),
      -batch_size    => $self->o('filter_duplicates_batch_size'),
      -can_be_empty  => 1,
      -rc_name       => '8GB_lastz',
    },

    {
      -logic_name => 'update_max_alignment_length_after_FD',
      -module     => 'Bio::EnsEMBL::Compara::RunnableDB::GenomicAlignBlock::UpdateMaxAlignmentLength',
      -parameters => {
        'quick' => $self->o('quick'),
      },
      -wait_for => [ 'create_filter_duplicates_jobs', 'filter_duplicates', 'filter_duplicates_himem' ],
      -rc_name => '2GB_lastz',
    },

    {
      -logic_name => 'no_chunk_and_group_dna',
      -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::ChunkAndGroupDna',
      -parameters => {
        'only_cellular_component' => $self->o('only_cellular_component'),
        'mix_cellular_components' => $self->o('mix_cellular_components'),
      },
      -wait_for  => ['update_max_alignment_length_after_FD'],
      -rc_name   => '2GB_lastz',
      -flow_into => {
        2 => ['dump_large_nib_for_chains'],
      },
    },

    {
      -logic_name => 'dump_large_nib_for_chains',
      -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::DumpDnaCollection',
      -parameters => {
        'faToNib_exe'       => $self->o('faToNib_exe'),
        'dump_min_nib_size' => $self->o('dump_min_nib_size'),
        'overwrite'         => 1,
      },
      -hive_capacity => 10,
      -rc_name       => '2GB_lastz',
      -flow_into     => {
        -1 => ['dump_large_nib_for_chains_himem'],    # MEMLIMIT
      },
    },

    {
      -logic_name => 'dump_large_nib_for_chains_himem',
      -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::DumpDnaCollection',
      -parameters => {
        'faToNib_exe'       => $self->o('faToNib_exe'),
        'dump_min_nib_size' => $self->o('dump_min_nib_size'),
        'overwrite'         => 1,
      },
      -hive_capacity => 10,
      -can_be_empty  => 1,
      -rc_name       => '8GB_lastz',
    },

    {
      -logic_name => 'create_alignment_chains_jobs',
      -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::CreateAlignmentChainsJobs',
      -parameters => {},
      -wait_for   => [ 'no_chunk_and_group_dna', 'dump_large_nib_for_chains', 'dump_large_nib_for_chains_himem' ],
      -rc_name    => '2GB_lastz',
      -flow_into  => {
        1 => ['remove_inconsistencies_after_chain'],
        2 => ['alignment_chains'],
      },
    },

    {
      -logic_name      => 'alignment_chains',
      -hive_capacity   => $self->o('chain_hive_capacity'),
      -batch_size      => $self->o('chain_batch_size'),
      -module          => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::AlignmentChains',
      -parameters      => $self->o('chain_parameters'),
      -max_retry_count => 4,
      -wait_for        => ['create_alignment_chains_jobs'],
      -rc_name         => '2GB_lastz',
      -flow_into       => {
        -1 => ['alignment_chains_himem'],
      },
    },

    {
      -logic_name      => 'alignment_chains_himem',
      -hive_capacity   => $self->o('chain_hive_capacity'),
      -batch_size      => 1,
      -module          => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::AlignmentChains',
      -parameters      => $self->o('chain_parameters'),
      -can_be_empty    => 1,
      -max_retry_count => 4,
      -rc_name         => '8GB_lastz',
      -wait_for        => ['alignment_chains'],
      -flow_into       => {
        -1 => ['alignment_chains_super_himem'],
      },
      -can_be_empty => 1,
    },

    {
      -logic_name      => 'alignment_chains_super_himem',
      -hive_capacity   => $self->o('chain_hive_capacity'),
      -batch_size      => 1,
      -module          => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::AlignmentChains',
      -parameters      => $self->o('chain_parameters'),
      -can_be_empty    => 1,
      -max_retry_count => 4,
      -rc_name         => '15GB_lastz',
      -wait_for        => ['alignment_chains'],
      -can_be_empty    => 1,
    },

    {
      -logic_name => 'remove_inconsistencies_after_chain',
      -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::RemoveAlignmentDataInconsistencies',
      -wait_for   => [ 'alignment_chains', 'alignment_chains_himem', 'alignment_chains_super_himem' ],
      -rc_name    => '2GB_lastz',
      -flow_into  => {
        1 => ['update_max_alignment_length_after_chain'],
      },
    },

    {
      -logic_name => 'update_max_alignment_length_after_chain',
      -module     => 'Bio::EnsEMBL::Compara::RunnableDB::GenomicAlignBlock::UpdateMaxAlignmentLength',
      -parameters => {
        'quick' => $self->o('quick'),
      },
      -rc_name => '2GB_lastz',
    },

    {
      -logic_name => 'create_alignment_nets_jobs',
      -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::CreateAlignmentNetsJobs',
      -parameters => {},
      -wait_for   => [ 'update_max_alignment_length_after_chain', 'remove_inconsistencies_after_chain' ],
      -rc_name    => '2GB_lastz',
      -flow_into  => {
        1 => ['remove_inconsistencies_after_net'],
        2 => ['alignment_nets'],
      },
    },

    {
      -logic_name    => 'alignment_nets',
      -hive_capacity => $self->o('net_hive_capacity'),
      -batch_size    => $self->o('net_batch_size'),
      -module        => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::AlignmentNets',
      -parameters    => $self->o('net_parameters'),
      -rc_name       => '2GB_lastz',
      -flow_into     => {
        -1 => ['alignment_nets_himem'],
      },
    },

    {
      -logic_name    => 'alignment_nets_himem',
      -hive_capacity => $self->o('net_hive_capacity'),
      -batch_size    => $self->o('net_batch_size'),
      -module        => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::AlignmentNets',
      -parameters    => $self->o('net_parameters'),
      -can_be_empty  => 1,
      -rc_name       => '4GB_lastz',
      -flow_into     => {
        -1 => ['alignment_nets_hugemem'],
      },
    },

    {
      -logic_name    => 'alignment_nets_hugemem',
      -hive_capacity => $self->o('net_hive_capacity'),
      -batch_size    => $self->o('net_batch_size'),
      -module        => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::AlignmentNets',
      -parameters    => $self->o('net_parameters'),
      -can_be_empty  => 1,
      -rc_name       => '8GB_lastz',
    },

    {
      -logic_name => 'remove_inconsistencies_after_net',
      -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::RemoveAlignmentDataInconsistencies',
      -wait_for   => [ 'alignment_nets', 'alignment_nets_himem', 'alignment_nets_hugemem', 'create_alignment_nets_jobs' ],    # Needed because of bi-directional netting: 2 jobs in create_alignment_nets_jobs can result in 1 job here
      -rc_name    => '2GB_lastz',
      -flow_into  => {
        1 => ['update_max_alignment_length_after_net'],
      },
    },

    {
      -logic_name   => 'create_filter_duplicates_net_jobs',
      -module       => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::CreateFilterDuplicatesJobs',
      -parameters   => {},
      -wait_for     => ['remove_inconsistencies_after_net'],
      -can_be_empty => 1,
      -rc_name      => '2GB_lastz',
      -flow_into    => {
        2 => { 'filter_duplicates_net' => INPUT_PLUS() },
      },
    },

    {
      -logic_name => 'filter_duplicates_net',
      -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::FilterDuplicates',
      -parameters => {
        'window_size'           => $self->o('compara_window_size'),
        'filter_duplicates_net' => 1,
      },
      -hive_capacity => $self->o('filter_duplicates_hive_capacity'),
      -batch_size    => $self->o('filter_duplicates_batch_size'),
      -can_be_empty  => 1,
      -rc_name       => '2GB_lastz',
      -flow_into     => {
        -1 => ['filter_duplicates_net_himem'],    # MEMLIMIT
      },
    },

    {
      -logic_name => 'filter_duplicates_net_himem',
      -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::FilterDuplicates',
      -parameters => {
        'window_size'           => $self->o('compara_window_size'),
        'filter_duplicates_net' => 1,
      },
      -hive_capacity => $self->o('filter_duplicates_hive_capacity'),
      -batch_size    => $self->o('filter_duplicates_batch_size'),
      -can_be_empty  => 1,
      -rc_name       => '8GB_lastz',
    },

    {
      -logic_name => 'update_max_alignment_length_after_net',
      -module     => 'Bio::EnsEMBL::Compara::RunnableDB::GenomicAlignBlock::UpdateMaxAlignmentLength',
      -rc_name    => '2GB_lastz',
      -wait_for   => [ 'create_filter_duplicates_net_jobs', 'filter_duplicates_net', 'filter_duplicates_net_himem' ],
      -flow_into  => ['set_internal_ids_collection'],
    },

    {
      -logic_name => 'set_internal_ids_collection',
      -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::SetInternalIdsCollection',
      -parameters => {
        'skip' => $self->o('patch_alignments'),
      },
      -flow_into => {
        2 => ['set_internal_ids_slow'],
      },
      -analysis_capacity => 1,
    },

    {
      -logic_name        => 'set_internal_ids_slow',
      -module            => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::SetInternalIdsSlow',
      -analysis_capacity => 1,
      -rc_name           => '8GB_lastz',
    },

    {
      -logic_name => 'healthcheck',
      -module     => 'Bio::EnsEMBL::Compara::RunnableDB::HealthCheck',
      -parameters => {
        'previous_db'      => $self->o('previous_db'),
        'ensembl_release'  => $self->o('ensembl_release'),
        'prev_release'     => $self->o('prev_release'),
        'max_percent_diff' => $self->o('patch_alignments') ? $self->o('max_percent_diff_patches') : $self->o('max_percent_diff'),
      },
      -wait_for => ['set_internal_ids_collection'],
      -rc_name  => '2GB_lastz',
    },

    {
      -logic_name => 'pairaligner_stats',
      -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::PairAlignerStats',
      -parameters => {
        'skip'                     => $self->o('skip_pairaligner_stats'),
        'dump_features'            => $self->o('dump_features_exe'),
        'compare_beds'             => $self->o('compare_beds_exe'),
        'create_pair_aligner_page' => $self->o('create_pair_aligner_page_exe'),
        'bed_dir'                  => $self->o('compara_bed_dir'),
        'ensembl_release'          => $self->o('ensembl_release'),
        'output_dir'               => $self->o('compara_feature_dir'),
      },
      -wait_for => ['healthcheck'],
      -rc_name  => '2GB_lastz',
    },
  ];
}

sub resource_classes {
  my $self = shift;

  return {
    '2GB_lastz'  => { LSF => [ $self->lsf_resource_builder( 'production', 2000 ), undef, '-reg_conf ' . $self->default_options->{registry_file} ] },
    '4GB_lastz'  => { LSF => [ $self->lsf_resource_builder( 'production', 4000 ), undef, '-reg_conf ' . $self->default_options->{registry_file} ] },
    '8GB_lastz'  => { LSF => [ $self->lsf_resource_builder( 'production', 8000 ), undef, '-reg_conf ' . $self->default_options->{registry_file} ] },
    '15GB_lastz' => { LSF => [ $self->lsf_resource_builder( 'production', 15000 ), undef, '-reg_conf ' . $self->default_options->{registry_file} ] },
    'default' => { LSF => $self->lsf_resource_builder( 'production', 900 ) },
    }
}

sub check_file_in_ensembl {
  my ( $self, $file_path ) = @_;
  push @{ $self->{'_ensembl_file_paths'} }, $file_path;
  return $self->o('enscode_root_dir') . '/' . $file_path;
}

1;

