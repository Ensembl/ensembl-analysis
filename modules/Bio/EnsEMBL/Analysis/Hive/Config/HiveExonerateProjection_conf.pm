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

package HiveExonerateProjection_conf;

use strict;
use warnings;
use File::Spec::Functions;

use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_analysis_settings);
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

    'species_name'              => '', # scientific name of the target species to project genes to, e.g. mus_musculus
    'production_name'           => '' || $self->o('species_name'),
    'production_name_modifier'  => '', # Do not set unless working with non-reference strains, breeds etc. Must include _ in modifier, e.g. _hni for medaka strain HNI
    'ensembl_release'           => $ENV{ENSEMBL_RELEASE}, # this is the current release version on staging to be able to get the correct database

    'dbowner'                   => '' || $ENV{EHIVE_USER} || $ENV{USER},
    'pipeline_name'             => '' || $self->o('production_name').$self->o('production_name_modifier').'_'.$self->o('ensembl_release'),
    
    'user_r'                    => '', # read only db user
    'user'                      => '', # write db user
    'password'                  => '', # password for write db user
    
    'server_set'                => '', # What server set to user, e.g. set1
    'pipe_db_server'            => '', # host for pipe db
    'databases_server'          => '', # host for general output dbs
    'dna_db_server'             => '', # host for dna db
    'pipe_db_port'              => '', # port for pipeline host
    'databases_port'            => '', # port for general output db host
    'dna_db_port'               => '', # port for dna db host
    'release_number'            => '' || $self->o('ensembl_release'),
    'output_path'               => '/path/to/output_path/', # output dir
    'compara_registry_file'     => $self->o('output_path').'/Databases.pm',

    'projection_source_db_name'    => '', # This is generally a pre-existing db, like the current human/mouse core for example
    'projection_source_db_server'  => '',
    'projection_source_db_port'    => '',
    'projection_source_production_name' => '', # e.g. 'homo_sapiens',

    'dna_db_name'                   => '',# dna db name of the target db, $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_core_'.$self->o('release_number'),

######################################################
##
## Mostly constant settings
##
#######################################################

    'skip_projection'           => '0',
    'skip_lastz'                => '0',

    'compara_db_name'    => 'leanne_ensembl_compara_95', # e.g. XXX_ensembl_compara_95
    'compara_db_server'  => 'mysql-ens-genebuild-prod-5',
    'compara_db_port'    => 4531,

    'projection_lastz_db_name'     => $self->o('pipe_db_name'),
    'projection_lastz_db_server'   => $self->o('pipe_db_server'),
    'projection_lastz_db_port'     => $self->o('pipe_db_port'),

    'pipe_db_name'                  => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_pipe_'.$self->o('release_number'),

    'reference_db_name'            => $self->o('dna_db_name'),
    'reference_db_server'          => $self->o('dna_db_server'),
    'reference_db_port'            => $self->o('dna_db_port'),

    'projection_db_server'  => $self->o('databases_server'),
    'projection_db_port'    => $self->o('databases_port'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'ensembl_release'              => $ENV{ENSEMBL_RELEASE}, # this is the current release version on staging to be able to get the correct database
    'production_db_server'         => 'mysql-ens-meta-prod-1',
    'production_db_port'           => '4483',

    ensembl_analysis_script           => catdir($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts'),

########################
# Extra db settings
########################

    'num_tokens' => 10,
    mysql_dump_options => '--max_allowed_packet=1000MB',

########################
# Executable paths
########################

    'minimap2_genome_index'  => $self->o('faidx_genome_file').'.mmi',
    'minimap2_path'          => '/hps/nobackup2/production/ensembl/fergal/coding/long_read_aligners/new_mm2/minimap2/minimap2',
    'paftools_path'          => '/hps/nobackup2/production/ensembl/fergal/coding/long_read_aligners/new_mm2/minimap2/misc/paftools.js',
    'minimap2_batch_size'    => '5000',

    'blast_type' => 'ncbi', # It can be 'ncbi', 'wu', or 'legacy_ncbi'
    'exonerate_path'         => catfile($self->o('software_base_path'), 'opt', 'exonerate09', 'bin', 'exonerate'),
    indicate_path  => catfile($self->o('binary_base'), 'indicate'),
    pmatch_path  => catfile($self->o('binary_base'), 'pmatch'),
    exonerate_annotation => catfile($self->o('binary_base'), 'exonerate'),
    samtools_path => catfile($self->o('binary_base'), 'samtools'), #You may need to specify the full path to the samtools binary
    picard_lib_jar => catfile($self->o('software_base_path'), 'Cellar', 'picard-tools', '2.6.0', 'libexec', 'picard.jar'), #You need to specify the full path to the picard library
    'cesar_path' => catdir($self->o('software_base_path'),'opt','cesar','bin'),

# Max internal stops for projected transcripts
    'projection_exonerate_padding'          => 5000,

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# No option below this mark should be modified
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#

    'flag_potential_pseudogenes_script' => catfile($self->o('ensembl_analysis_script'),'genebuild','flag_potential_pseudogenes.pl'),

##################################
# Memory settings for the analyses
##################################
    'default_mem'          => '900',
    'genblast_mem'         => '1900',
    'genblast_retry_mem'   => '4900',
    'genewise_mem'         => '3900',
    'genewise_retry_mem'   => '5900',
    'refseq_mem'           => '9900',
    'projection_mem'       => '1900',
    'layer_annotation_mem' => '3900',
    'genebuilder_mem'      => '1900',

########################
# LastZ
########################

    'method_link_type'     => 'LASTZ_NET',
    'compara_databases_conf_filename' => 'Databases.pm',

    'compara_master'             => 'compara_master',
    'compara_conf_file'             => '',
    'compara_innodb_schema'         => 1,
    'compara_genome_db_update_path' => catfile($self->o('enscode_root_dir'),'/ensembl-compara/scripts/pipeline/update_genome.pl'),
    'compara_mlss_script_path'      => catfile($self->o('enscode_root_dir'),'/ensembl-compara/scripts/pipeline/create_mlss.pl'),
    'compara_mlss_reg_conf_path'    => catfile($self->o('enscode_root_dir'),'/ensembl-compara/scripts/pipeline/production_reg_ensembl_conf.pl'),
    'compara_populate_new_database_exe' => catfile($self->o('enscode_root_dir'),'ensembl-compara/scripts/pipeline/populate_new_database.pl'),
    'compara_only_cellular_component' => undef,
    'compara_dump_dir'              => catdir($self->o('output_path'),'lastz'),

    'mlss_id_list' => undef,
    'compara_collection' => '',

    'compara_ref_species'       => $self->o('projection_source_production_name'),
    'compara_non_ref_species'   => $self->o('production_name'),
    'only_cellular_component'   => undef,   # Do we load *all* the dnafrags or only the ones from a specific cellular-component ?
    'mix_cellular_components'   => 0,       # Do we try to allow the nuclear genome vs MT, etc ?
    'dump_min_nib_size'         => 11500000,
    'dump_min_chunk_size'       => 1000000,
    'dump_min_chunkset_size'    => 1000000,
    'quick' => 1,
    'default_chunks' => {
      'reference'   => {
        'homo_sapiens' => {
          'chunk_size' => 30000000,
          'overlap'    => 0,
          'include_non_reference' => -1, #1  => include non_reference regions (eg human assembly patches)
                                         #0  => do not include non_reference regions
                                         #-1 => auto-detect (only include non_reference regions if the non-reference species is high-coverage
                                         #ie has chromosomes since these analyses are the only ones we keep up-to-date with the patches-pipeline)
          'masking_options' => '{default_soft_masking => 1}',
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

    'compara_window_size' => 10000,
    'filter_duplicates_rc_name' => '2GB_lastz',
    'filter_duplicates_himem_rc_name' => '8GB_lastz',

   #
    #Default pair_aligner
    #
    'pair_aligner_method_link' => [1001, 'LASTZ_RAW'],
    'pair_aligner_logic_name' => 'LastZ',
    'pair_aligner_module' => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::LastZ',

    'pair_aligner_options' => {
       default => 'T=1 L=3000 H=2200 O=400 E=30 --ambiguous=iupac', # ensembl genomes settings
       7742    => 'T=1 K=3000 L=3000 H=2200 O=400 E=30 --ambiguous=iupac', # vertebrates - i.e. ensembl-specific
       9526    => 'T=1 K=5000 L=5000 H=3000 M=10 O=400 E=30 Q=' . $self->check_file_in_ensembl('ensembl-compara/scripts/pipeline/primate.matrix').' --ambiguous=iupac', # primates
       33554   => 'T=1 K=5000 L=5000 H=3000 M=10 O=400 E=30 --ambiguous=iupac', # carnivora
       3913    => 'T=1 L=3000 H=2200 O=400 E=30 --ambiguous=iupac --matchcount=1000',
       4070    => 'T=1 L=3000 H=2200 O=400 E=30 --ambiguous=iupac --matchcount=1000',
    },

    #
    #Default chain
    #
    'chain_input_method_link' => [1001, 'LASTZ_RAW'],
    'chain_output_method_link' => [1002, 'LASTZ_CHAIN'],

    #linear_gap=>medium for more closely related species, 'loose' for more distant
    'linear_gap' => 'medium',

    'chain_parameters' => {'max_gap'=>'50','linear_gap'=> $self->o('linear_gap'), 'faToNib' => $self->o('faToNib_exe'), 'lavToAxt'=> $self->o('lavToAxt_exe'), 'axtChain'=>$self->o('axtChain_exe'), 'max_blocks_for_chaining' => 100000},

    #
    #Default patch_alignments
    #
    'patch_alignments' => 0,  #set to 1 to align the patches of a species to many other species

    #
    #Default net
    #
    'net_input_method_link' => [1002, 'LASTZ_CHAIN'],
    'net_output_method_link' => [16, 'LASTZ_NET'],
    'net_ref_species' => $self->o('compara_ref_species'),  #default to ref_species
    'net_parameters' => {'max_gap'=>'50', 'chainNet'=>$self->o('chainNet_exe')},
    'bidirectional' => 0,

    #
    #Default healthcheck
    #
    'previous_db' => 'compara_prev',
    'prev_release' => 0,   # 0 is the default and it means "take current release number and subtract 1"
    'max_percent_diff' => 20,
    'max_percent_diff_patches' => 99.99,
    'do_pairwise_gabs' => 1,
    'do_compare_to_previous_db' => 0,

    'compara_bed_dir' => $self->o('compara_dump_dir').'/bed_dir',
    'compara_feature_dir' => $self->o('compara_dump_dir').'/feature_dumps',

    #
    #Default pairaligner config
    #
    'skip_pairaligner_stats' => 1, #skip this module if set to 1

    'pair_aligner_method_link' => [1001, 'LASTZ_RAW'],
    'pair_aligner_logic_name' => 'LastZ',
    'pair_aligner_module' => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::LastZ',
    'chain_input_method_link' => [1001, 'LASTZ_RAW'],
    'chain_output_method_link' => [1002, 'LASTZ_CHAIN'],
    'linear_gap' => 'medium',
    'net_input_method_link' => [1002, 'LASTZ_CHAIN'],
    'net_output_method_link' => [16, 'LASTZ_NET'],

    # Capacities
    'pair_aligner_analysis_capacity' => 700,
    'pair_aligner_batch_size' => 40,
    'chain_hive_capacity' => 200,
    'chain_batch_size' => 10,
    'net_hive_capacity' => 300,
    'net_batch_size' => 10,
    'filter_duplicates_hive_capacity' => 200,
    'filter_duplicates_batch_size' => 10,

    # LastZ is used to align the genomes
    'pair_aligner_exe'  => $self->o('lastz_exe'),
    'cellar_dir'                        => '/nfs/software/ensembl/RHEL7-JUL2017-core2/linuxbrew/Cellar/',
    'lastz_exe'                         => catfile($self->o('cellar_dir'),'lastz/1.04.00/bin/lastz'),
    'axtChain_exe'                      => catfile($self->o('cellar_dir'),'kent/v335_1/bin/axtChain'),
    'chainNet_exe'                      => catfile($self->o('cellar_dir'),'kent/v335_1/bin/chainNet'),
    'faToNib_exe'                       => catfile($self->o('cellar_dir'),'kent/v335_1/bin/faToNib'),
    'lavToAxt_exe'                      => catfile($self->o('cellar_dir'),'kent/v335_1/bin/lavToAxt'),
    'compare_beds_exe'                  => catfile($self->o('enscode_root_dir'),'ensembl-compara/scripts/pipeline/compare_beds.pl'),
    'create_pair_aligner_page_exe'      => catfile($self->o('enscode_root_dir'),'ensembl-compara/scripts/report/create_pair_aligner_page.pl'),
    'dump_features_exe'                 => catfile($self->o('enscode_root_dir'),'ensembl-compara/scripts/dumps/DumpMultiAlign.pl'),


########################
# db info
########################
    'reference_db' => {
      -dbname => $self->o('reference_db_name'),
      -host   => $self->o('reference_db_server'),
      -port   => $self->o('reference_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

   'compara_db' => {
      -dbname => $self->o('compara_db_name'),
      -host   => $self->o('compara_db_server'),
      -port   => $self->o('compara_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'projection_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_proj_'.$self->o('release_number'),
      -host   => $self->o('projection_db_server'),
      -port   => $self->o('projection_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'selected_projection_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_sel_proj_'.$self->o('release_number'),
      -host   => $self->o('projection_db_server'),
      -port   => $self->o('projection_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'projection_source_db' => {
      -dbname => $self->o('projection_source_db_name'),
      -host   => $self->o('projection_source_db_server'),
      -port   => $self->o('projection_source_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -driver => $self->o('hive_driver'),
    },

    'projection_lastz_db' => {
      -dbname => $self->o('projection_lastz_db_name'),
      -host   => $self->o('projection_lastz_db_server'),
      -port   => $self->o('projection_lastz_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -driver => $self->o('hive_driver'),
    },

  };
}

sub pipeline_create_commands {
    my ($self) = @_;

################
# LastZ
################

    my $second_pass     = exists $self->{'_is_second_pass'};
    $self->{'_is_second_pass'} = $second_pass;
    return $self->SUPER::pipeline_create_commands if $self->can('no_compara_schema');
    my $pipeline_url    = $self->pipeline_url();
    my $parsed_url      = $second_pass && Bio::EnsEMBL::Hive::Utils::URL::parse( $pipeline_url );
    my $driver          = $second_pass ? $parsed_url->{'driver'} : '';

################
# /LastZ
################

    return [
    # inheriting database and hive tables' creation
      @{$self->SUPER::pipeline_create_commands},

#################
# LastZ
#################

     'mkdir -p '.$self->o('compara_dump_dir'),
     'mkdir -p '.$self->o('compara_bed_dir'),
      # Compara 'release' tables will be turned from MyISAM into InnoDB on the fly by default:
      ($self->o('compara_innodb_schema') ? "sed 's/ENGINE=MyISAM/ENGINE=InnoDB/g' " : 'cat ')
      . $self->check_file_in_ensembl('ensembl-compara/sql/table.sql').' | '.$self->db_cmd(),

      # Compara 'pipeline' tables are already InnoDB, but can be turned to MyISAM if needed:
      ($self->o('compara_innodb_schema') ? 'cat ' : "sed 's/ENGINE=InnoDB/ENGINE=MyISAM/g' ")
      . $self->check_file_in_ensembl('ensembl-compara/sql/pipeline-tables.sql').' | '.$self->db_cmd(),

      # MySQL specific procedures
      $driver eq 'mysql' ? ($self->db_cmd().' < '.$self->check_file_in_ensembl('ensembl-compara/sql/procedures.'.$driver)) : (),

#################
# /LastZ
#################


    ];
}


sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
    skip_lastz => $self->o('skip_lastz'),
    wide_ensembl_release => $self->o('ensembl_release'),
  }
}

## See diagram for pipeline structure
sub pipeline_analyses {
    my ($self) = @_;

    my %genblast_params = (
      wu    => '-P wublast -gff -e #blast_eval# -c #blast_cov#',
      ncbi  => '-P blast -gff -e #blast_eval# -c #blast_cov# -W 3 -softmask -scodon 50 -i 30 -x 10 -n 30 -d 200000 -g T',
      wu_genome    => '-P wublast -gff -e #blast_eval# -c #blast_cov#',
      ncbi_genome  => '-P blast -gff -e #blast_eval# -c #blast_cov# -W 3 -softmask -scodon 50 -i 30 -x 10 -n 30 -d 200000 -g T',
      wu_projection    => '-P wublast -gff -e #blast_eval# -c #blast_cov# -n 100 -x 5 ',
      ncbi_projection  => '-P blast -gff -e #blast_eval# -c #blast_cov# -W 3 -scodon 50 -i 30 -x 10 -n 30 -d 200000 -g T',
      );
    my %commandline_params = (
      'ncbi' => '-num_threads 3 -window_size 40',
      'wu' => '-cpus 3 -hitdist 40',
      'legacy_ncbi' => '-a 3 -A 40',
      );
    my $header_line = create_header_line($self->default_options->{'file_columns'});

    return [

########################################################################
#
# Projection analyses
#
########################################################################
     
      {
        -logic_name => 'fan_projection',
        -input_ids  => [
                        {
          },
        ],
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'if [ "#skip_projection#" -ne "0" ]; then exit 42; else exit 0;fi',
          return_codes_2_branches => {'42' => 2},
        },
        -rc_name    => 'default',
        -flow_into  => {
           '1->A' => ['fan_lastz'],
           'A->1' => ['create_projection_db'],
        },
      },

      {
        -logic_name => 'create_projection_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('projection_db'),
                         create_type => 'clone',
                         #script_path => $self->o('clone_db_script_path'),
                         #user_r => $self->o('user_r'),
                         #user_w => $self->o('user_w'),
                         #pass_w => $self->o('password'),
                       },
        -rc_name    => 'default',
        -flow_into  => {
          '1->A' => ['exonerate_create_projection_input_ids'],#['cesar_create_projection_input_ids','wga_create_projection_input_ids'],
          'A->1' => ['classify_projected_genes'],
        },
      },

      {
        -logic_name => 'exonerate_create_projection_input_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         #skip_analysis => $self->o('skip_projection'),
                         target_db => $self->o('projection_source_db'),
                         iid_type => 'feature_id',
                         feature_type => 'transcript',
                         feature_restriction => 'projection',
                         biotypes => {
                                       'protein_coding' => 1,
                                       'Mt_tRNA' => 1,
                                       'Mt_rRNA' => 1,
                                       'misc_RNA' => 1,
                                       'snoRNA' => 1,
                                       'miRNA' => 1,
                                       'protein_coding' => 1,
                                       'Mt_tRNA' => 1,
                                       'Mt_rRNA' => 1,
                                       'misc_RNA' => 1,
                                       'snoRNA' => 1,
                                       'miRNA' => 1,
                                       'snRNA' => 1,
                                       'pseudogene' => 1,
                                       'scaRNA' => 1,
                                       'rRNA' => 1,
                                       'processed_pseudogene' => 1,
                                       'ribozyme' => 1,
                                       'sRNA' => 1,
                                       'lincRNA' => 1,
                                       'IG_LV_gene' => 1,
                                       'processed_transcript' => 1,
                                       'TEC' => 1,
                                       'transcribed_processed_pseudogene' => 1,
                                       'unprocessed_pseudogene' => 1,
                                       'sense_intronic' => 1,
                                       'antisense' => 1,
                                       'bidirectional_promoter_lncRNA' => 1,
                                       'transcribed_unprocessed_pseudogene' => 1,
                                       'sense_overlapping' => 1,
                                       'unitary_pseudogene' => 1,
                                       'transcribed_unitary_pseudogene' => 1,
                                       'TR_J_gene' => 1,
                                       'TR_V_gene' => 1,
                                       'IG_C_gene' => 1,
                                       'TR_C_gene' => 1,
                                       'IG_C_pseudogene' => 1,
                                       'IG_J_gene' => 1,
                                       'polymorphic_pseudogene' => 1,
                                       'IG_D_gene' => 1,
                                       '3prime_overlapping_ncRNA' => 1,
                                       'IG_D_pseudogene' => 1,
                                       'IG_pseudogene' => 1,
                                       'IG_V_pseudogene' => 1,
                                       'IG_V_gene' => 1,
                                       'TR_V_pseudogene' => 1,
                                       'TR_D_gene' => 1,
                                       'TR_J_pseudogene' => 1,
                                       'scRNA' => 1,
                                       'translated_unprocessed_pseudogene' => 1,
                                       'macro_lncRNA' => 1,
                                    }, 
                         batch_size => 25,
        },
        -flow_into => {
                        2 => ['exonerate_project_transcripts'],
        },
        -rc_name    => 'default',
      },

      {
        -logic_name => 'exonerate_project_transcripts',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProjectionExonerate',
        -parameters => {
                         'logic_name' => 'exonerate_projection',
                         'source_dna_db' => $self->default_options()->{'projection_source_db'},
                         'target_dna_db' => $self->o('dna_db'),
                         'source_db' => $self->o('projection_source_db'),
                         'target_db' => $self->o('projection_db'),
                         'compara_db' => $self->o('compara_db'),
                         'method_link_type' => $self->o('method_link_type'),
                         'exon_region_padding' => $self->o('projection_exonerate_padding'),
                         'exonerate_path' => '/nfs/software/ensembl/RHEL7/linuxbrew/bin/exonerate',
                         'exonerate_coverage' => 80,
                         #'exonerate_percent_id' => 80,
                         'exonerate_percent_id' => 60,
                         #'calculate_coverage_and_pid' => 1,
                         'generate_annotation_file' => 1,
                         #%{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_projection_coding')},
                         #%{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_cdna')},
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','cdna2genome')},
                       },
        -rc_name    => 'default',
        -hive_capacity => 900,
        #-input_ids => $self->o('test_input_id'),
        -flow_into => {
                        -3 => ['failed_coding_jobs'],
                      },
      },
      
      {
        -logic_name => 'failed_coding_jobs',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name          => 'default',
        -can_be_empty  => 1,
      },

      {
        -logic_name => 'classify_projected_genes',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport',
        -parameters => {
          skip_analysis => $self->o('skip_projection'),
          classification_type => 'standard',
          update_gene_biotype => 1,
          target_db => $self->o('projection_db'),
        },
        -rc_name    => '2GB',
        -flow_into => {
          1 => ['flag_problematic_projections'],
#          When the realign part is fix, the line above needs to be deleted and the lines below uncommented
#          All analysis below create_projection_realign_db need to be uncommented too
#          '1->A' => ['fix_unaligned_protein_hit_names'],
#          'A->1' => ['create_projection_realign_db'],
        },
      },

      {
        -logic_name => 'flag_problematic_projections',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('flag_potential_pseudogenes_script').
                                ' -host '.$self->o('projection_db','-host').
                                ' -port '.$self->o('projection_db','-port').
                                ' -user_w '.$self->o('projection_db','-user').
                                ' -pass '.$self->o('projection_db','-pass').
                                ' -dbname '.$self->o('projection_db','-dbname').
                                ' -dna_host '.$self->o('dna_db','-host').
                                ' -dna_port '.$self->o('dna_db','-port').
                                ' -user_r '.$self->o('dna_db','-user').
                                ' -dna_dbname '.$self->o('dna_db','-dbname'),
                       },
        -rc_name => '2GB',
        -flow_into  => {
          1 => ['fix_projection_db_issues'],
        },
      },

      {
        # This will fix issues when proteins that were too long for MUSCLE alignment didn't have proper ENST accessions
        # Will also update the transcript table to match the gene table once the pseudo/canon genes are tagged in the
        # previous analysis
        -logic_name => 'fix_projection_db_issues',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn    => $self->o('projection_db'),
          sql => [
            'UPDATE gene JOIN transcript USING(gene_id) SET transcript.biotype=gene.biotype',
            'UPDATE protein_align_feature JOIN transcript_supporting_feature ON feature_id = protein_align_feature_id'.
              ' JOIN transcript USING(transcript_id) SET hit_name = stable_id',
            'UPDATE protein_align_feature JOIN supporting_feature ON feature_id = protein_align_feature_id'.
              ' JOIN exon_transcript USING(exon_id) JOIN transcript USING(transcript_id) SET hit_name = stable_id',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
          1 => ['projection_sanity_checks'],
        },
      },

      {
        -logic_name => 'projection_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
          target_db => $self->o('projection_db'),
          sanity_check_type => 'gene_db_checks',
          min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
            'gene_db_checks')->{$self->o('uniprot_set')}->{'projection_coding'},
        },
        -rc_name    => '4GB',
        #-flow_into => {
        #                1 => ['select_projected_genes'],
        #              },
      },

######################################################################################
#
# LastZ pipeline
#
######################################################################################

      {
        -logic_name => 'fan_lastz',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'if [ "#skip_lastz#" -ne "0" ]; then exit 42; else exit 0;fi',
                         return_codes_2_branches => {'42' => 2},
                       },
        -rc_name => '2GB_lastz',
        -flow_into  => {
          '1' => ['insert_projection_source_assembly_into_compara_db'],
        },
      },

      {
        -logic_name => 'insert_projection_source_assembly_into_compara_db',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => { 
                         cmd => 'perl '.$self->o('compara_genome_db_update_path').
                                ' --reg_conf '.$self->o('output_path').'/'.$self->o('compara_databases_conf_filename').
                                ' --compara mysql://'.$self->o('user').':'.$self->o('password').'@'.$self->o('compara_db_server').':'.$self->o('compara_db_port').'/'.$self->o('compara_db_name').
                                ' --species "'.$self->o('projection_source_production_name').'" --force'
                       },
        -flow_into  => {
          '1' => ['setup_lastz'],
        },
      },

	    {
        -logic_name => 'setup_lastz',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::LastZSetup',
        -parameters => {
                         compara_genome_db_update_path => $self->o('compara_genome_db_update_path'),
                         compara_mlss_script_path => $self->o('compara_mlss_script_path'),
                         compara_mlss_reg_conf_path => $self->o('compara_mlss_reg_conf_path'),
                         compara_db => $self->o('compara_db'),
                         projection_source_db => $self->o('projection_source_db'),
                         target_db => $self->o('dna_db'),
                         pipeline_db => $self->o('pipeline_db'),
                         output_path => $self->o('output_path'),
                         compara_db_url => 'mysql://'.$self->o('user').':'.$self->o('password').'@'.$self->o('compara_db_server').':'.$self->o('compara_db_port').'/'.$self->o('compara_db_name'),
                       },

        -rc_name => '2GB_lastz',
             -flow_into      => {
                             1 => {'get_species_list' => {'mlss_id' => '#mlss_id#', 'reg_conf' => '#reg_conf#'}},
                           },
       },


	    {
         -logic_name    => 'get_species_list',
         -module        => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::ParsePairAlignerConf',
         -parameters    => {
                            'master_db' => $self->o('compara_master'),
                            'get_species_list' => 1,
                            'core_dbs' => undef,
                          },

         -rc_name => '2GB_lastz',
         -flow_into      => {
                              1 => {'populate_new_database' => {'mlss_id' => '#mlss_id#', 'reg_conf' => '#reg_conf#'}},
                            },
       },


	    {
         -logic_name => 'populate_new_database',
         -module     => 'Bio::EnsEMBL::Compara::RunnableDB::GenomicAlignBlock::PopulateNewDatabase',
         -parameters    => {
                             'program'        => $self->o('compara_populate_new_database_exe'),
                             'mlss_id_list'   => undef,
                             'collection'     => $self->o('compara_collection'),
                             'master_db'      => $self->o('compara_master'),
                             'only_cellular_component' => $self->o('only_cellular_component'),
                           },
         -flow_into => {
           1 => {'add_method_link_species_link_tag' => {'mlss_id' => '#mlss_id#', 'reg_conf' => '#reg_conf#'}},
         },
         -rc_name => '2GB_lastz',
       },


	    {
        -logic_name => 'add_method_link_species_link_tag',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('pipeline_db'),
          sql => [
            'INSERT INTO method_link_species_set_tag VALUES (#mlss_id#,"reference_species","'.$self->o('compara_ref_species').'")',
          ],
        },
        -rc_name    => '2GB_lastz',
        -flow_into => {
          1 => {'parse_pair_aligner_conf' => {'mlss_id' => '#mlss_id#', 'reg_conf' => '#reg_conf#'}}
        },
      },


	    {
         -logic_name    => 'parse_pair_aligner_conf',
         -module        => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::ParsePairAlignerConf',
         -parameters    => {
                             'conf_file' => $self->o('compara_conf_file'),
                             'ref_species' => $self->o('compara_ref_species'),
                             'non_ref_species' => $self->o('compara_non_ref_species'),
                             'dump_dir' => $self->o('compara_dump_dir'),
                             'default_chunks' => $self->o('default_chunks'),
                             'default_pair_aligner' => $self->o('pair_aligner_method_link'),
                             'default_parameters' => $self->o('pair_aligner_options'),
                             'default_chain_output' => $self->o('chain_output_method_link'),
                             'default_net_output' => $self->o('net_output_method_link'),
                             'default_chain_input' => $self->o('chain_input_method_link'),
                             'default_net_input' => $self->o('net_input_method_link'),
                             'net_ref_species' => $self->o('net_ref_species'),
                             'mlss_id_list' => $self->o('mlss_id_list'),
                             'collection' => $self->o('compara_collection'),
                             'master_db' => $self->o('compara_master'),
                             'do_pairwise_gabs' => $self->o('do_pairwise_gabs'), #healthcheck options
                             'do_compare_to_previous_db' => $self->o('do_compare_to_previous_db'), #healthcheck options
                             'bidirectional' => $self->o('bidirectional'),
                           },
         -flow_into => {
                         1 => [ 'create_pair_aligner_jobs'],
                         2 => [ 'chunk_and_group_dna' ],
                         3 => [ 'create_filter_duplicates_jobs' ],
                         4 => [ 'no_chunk_and_group_dna' ],
                         5 => [ 'create_alignment_chains_jobs' ],
                         6 => [ 'create_alignment_nets_jobs' ],
                        10 => [ 'create_filter_duplicates_net_jobs' ],
                         7 => [ 'pairaligner_stats' ],
                         8 => [ 'healthcheck' ],
                       },
         -rc_name => '2GB_lastz',
       },


	    {
         -logic_name => 'chunk_and_group_dna',
         -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::ChunkAndGroupDna',
         -parameters => {
                          'only_cellular_component' => $self->o('only_cellular_component'),
                          'mix_cellular_components' => $self->o('mix_cellular_components'),
                        },
         -flow_into => {
           2 => [ 'store_sequence' ],
         },
         -rc_name => '2GB_lastz',
       },


	    {
         -logic_name => 'store_sequence',
         -hive_capacity => 100,
         -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::StoreSequence',
         -parameters => {
           'dump_min_chunkset_size' => $self->o('dump_min_chunkset_size'),
           'dump_min_chunk_size' => $self->o('dump_min_chunk_size'),
         },
         -flow_into => {
           -1 => [ 'store_sequence_again' ],
         },
         -rc_name => '2GB_lastz',
       },


	    {
         -logic_name => 'store_sequence_again',
         -hive_capacity => 50,
         -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::StoreSequence',
         -parameters => {
           'dump_min_chunkset_size' => $self->o('dump_min_chunkset_size'),
           'dump_min_chunk_size' => $self->o('dump_min_chunk_size'),
          },
          -can_be_empty  => 1,
          -rc_name => '4GB_lastz',
       },


	    {
         -logic_name => 'create_pair_aligner_jobs',  #factory
         -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::CreatePairAlignerJobs',
         -parameters => {
                          'mix_cellular_components' => $self->o('mix_cellular_components'),
                        },
         -hive_capacity => 10,
         -wait_for => [ 'store_sequence', 'store_sequence_again', 'chunk_and_group_dna'  ],
         -flow_into => {
                         1 => [ 'check_no_partial_gabs' ],
                         2 => [ $self->o('pair_aligner_logic_name')  ],
                       },
         -rc_name => '2GB_lastz',
       },


	    {
         -logic_name => $self->o('pair_aligner_logic_name'),
         -module     => $self->o('pair_aligner_module'),
         -analysis_capacity => $self->o('pair_aligner_analysis_capacity'),
         -batch_size => $self->o('pair_aligner_batch_size'),
         -parameters => {
                          'pair_aligner_exe' => $self->o('pair_aligner_exe'),
                        },
         -wait_for  => [ 'create_pair_aligner_jobs'  ],
         -flow_into => {
                         -1 => [ $self->o('pair_aligner_logic_name') . '_himem1' ],  # MEMLIMIT
                       },
         -rc_name => '2GB_lastz',
       },


	    {
         -logic_name => $self->o('pair_aligner_logic_name') . "_himem1",
         -module     => $self->o('pair_aligner_module'),
         -analysis_capacity => $self->o('pair_aligner_analysis_capacity'),
         -parameters => {
                          'pair_aligner_exe' => $self->o('pair_aligner_exe'),
                        },
         -wait_for   => [ 'create_pair_aligner_jobs'  ],
         -batch_size => $self->o('pair_aligner_batch_size'),
         -can_be_empty  => 1,
         -rc_name => '8GB_lastz',
       },


	    {
         -logic_name => 'check_no_partial_gabs',
         -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::SqlHealthChecks',
         -parameters => {
           'mode'          => 'gab_inconsistencies',
         },
         -wait_for =>  [ $self->o('pair_aligner_logic_name'), $self->o('pair_aligner_logic_name') . "_himem1" ],
         -flow_into => {
           1 => [ 'update_max_alignment_length_before_FD' ],
         },
       },


	    {
         -logic_name => 'update_max_alignment_length_before_FD',
         -module     => 'Bio::EnsEMBL::Compara::RunnableDB::GenomicAlignBlock::UpdateMaxAlignmentLength',
         -parameters => {
           'quick' => $self->o('quick'),
         },
         -flow_into => {
           1 => [ 'update_max_alignment_length_after_FD' ],
         },
         -rc_name => '2GB_lastz',
       },


	    {
         -logic_name => 'create_filter_duplicates_jobs', #factory
         -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::CreateFilterDuplicatesJobs',
         -parameters => { },
         -wait_for =>  [ 'update_max_alignment_length_before_FD' ],
         -flow_into => {
           2 => { 'filter_duplicates' => INPUT_PLUS() },
         },
         -rc_name => '2GB_lastz',
       },


	    {
         -logic_name   => 'filter_duplicates',
         -module        => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::FilterDuplicates',
         -parameters    => {
           'window_size' => $self->o('compara_window_size')
         },
         -hive_capacity => $self->o('filter_duplicates_hive_capacity'),
         -batch_size    => $self->o('filter_duplicates_batch_size'),
         -can_be_empty  => 1,
         -flow_into => {
                         -1 => [ 'filter_duplicates_himem' ], # MEMLIMIT
         },
         -rc_name => $self->o('filter_duplicates_rc_name'),
       },


	    {
        -logic_name   => 'filter_duplicates_himem',
        -module        => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::FilterDuplicates',
        -parameters    => {
                            'window_size' => $self->o('compara_window_size')
         },
         -hive_capacity => $self->o('filter_duplicates_hive_capacity'),
         -batch_size    => $self->o('filter_duplicates_batch_size'),
         -can_be_empty  => 1,
         -rc_name => $self->o('filter_duplicates_himem_rc_name'),
      },


	    {
        -logic_name => 'update_max_alignment_length_after_FD',
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::GenomicAlignBlock::UpdateMaxAlignmentLength',
        -parameters => {
          'quick' => $self->o('quick'),
        },
        -wait_for =>  [ 'create_filter_duplicates_jobs', 'filter_duplicates', 'filter_duplicates_himem' ],
        -rc_name => '2GB_lastz',
      },


	    {
        -logic_name => 'no_chunk_and_group_dna',
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::ChunkAndGroupDna',
        -parameters => {
          'only_cellular_component' => $self->o('only_cellular_component'),
          'mix_cellular_components' => $self->o('mix_cellular_components'),
        },
        -flow_into => {
          2 => [ 'dump_large_nib_for_chains' ],
        },
        -wait_for  => ['update_max_alignment_length_after_FD' ],
        -rc_name => '2GB_lastz',
      },


	    {
        -logic_name => 'dump_large_nib_for_chains',
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::DumpDnaCollection',
        -parameters => {
          'faToNib_exe' => $self->o('faToNib_exe'),
          'dump_min_nib_size' => $self->o('dump_min_nib_size'),
          'overwrite'=>1,
        },
        -hive_capacity => 10,
        -flow_into => {
          -1 => [ 'dump_large_nib_for_chains_himem' ],  # MEMLIMIT
        },
        -rc_name => '2GB_lastz',
      },


	    {
        -logic_name => 'dump_large_nib_for_chains_himem',
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::DumpDnaCollection',
        -parameters => {
                         'faToNib_exe' => $self->o('faToNib_exe'),
                         'dump_min_nib_size' => $self->o('dump_min_nib_size'),
                         'overwrite'=>1,
        },
        -hive_capacity => 10,
        -can_be_empty  => 1,
        -rc_name => '8GB_lastz',
      },



	    {
        -logic_name => 'create_alignment_chains_jobs',
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::CreateAlignmentChainsJobs',
        -parameters => { },
        -flow_into => {
                        1 => [ 'remove_inconsistencies_after_chain' ],
                        2 => [ 'alignment_chains' ],
        },
        -wait_for => [ 'no_chunk_and_group_dna', 'dump_large_nib_for_chains', 'dump_large_nib_for_chains_himem' ],
        -rc_name => '2GB_lastz',
      },


	    {
        -logic_name => 'alignment_chains',
        -hive_capacity => $self->o('chain_hive_capacity'),
        -batch_size => $self->o('chain_batch_size'),
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::AlignmentChains',
        -parameters => $self->o('chain_parameters'),
        -max_retry_count => 10,
        -flow_into => {
          -1 => [ 'alignment_chains_himem' ],  # MEMLIMIT
        },
        -wait_for   => [ 'create_alignment_chains_jobs' ],
        -rc_name => '2GB_lastz',
      },


	    {
        -logic_name => 'alignment_chains_himem',
        -hive_capacity => $self->o('chain_hive_capacity'),
        -batch_size => 1,
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::AlignmentChains',
        -parameters => $self->o('chain_parameters'),
        -can_be_empty  => 1,
        -max_retry_count => 10,
        -rc_name => '8GB_lastz',
        -flow_into => {
          -1 => [ 'alignment_chains_super_himem' ],  # MEMLIMIT
         },
         -wait_for => ['alignment_chains'],
         -can_be_empty  => 1,
      },


	    {
        -logic_name => 'alignment_chains_super_himem',
        -hive_capacity => $self->o('chain_hive_capacity'),
        -batch_size => 1,
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::AlignmentChains',
        -parameters => $self->o('chain_parameters'),
        -can_be_empty  => 1,
        -max_retry_count => 10,
        -rc_name => '15GB_lastz',
        -wait_for => ['alignment_chains'],
         -can_be_empty  => 1,
      },

	    {
        -logic_name => 'remove_inconsistencies_after_chain',
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::RemoveAlignmentDataInconsistencies',
        -flow_into => {
          1 => [ 'update_max_alignment_length_after_chain' ],
        },
        -wait_for =>  [ 'alignment_chains', 'alignment_chains_himem', 'alignment_chains_super_himem' ],
        -rc_name => '2GB_lastz',
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
        -parameters => { },
        -flow_into => {
          1 => [ 'remove_inconsistencies_after_net' ],
          2 => [ 'alignment_nets' ],
        },
        -wait_for => [ 'update_max_alignment_length_after_chain', 'remove_inconsistencies_after_chain' ],
        -rc_name => '2GB_lastz',
      },


	    {
        -logic_name => 'alignment_nets',
        -hive_capacity => $self->o('net_hive_capacity'),
        -batch_size => $self->o('net_batch_size'),
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::AlignmentNets',
        -parameters => $self->o('net_parameters'),
        -flow_into => {
          -1 => [ 'alignment_nets_himem' ],  # MEMLIMIT
         },
         -rc_name => '2GB_lastz',
      },


	    {
        -logic_name => 'alignment_nets_himem',
        -hive_capacity => $self->o('net_hive_capacity'),
        -batch_size => $self->o('net_batch_size'),
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::AlignmentNets',
        -parameters => $self->o('net_parameters'),
        -can_be_empty => 1,
        -flow_into => {
          -1 => [ 'alignment_nets_hugemem' ],  # MEMLIMIT
        },
        -rc_name => '4GB_lastz',
      },


	    {
        -logic_name     => 'alignment_nets_hugemem',
        -hive_capacity  => $self->o('net_hive_capacity'),
        -batch_size     => $self->o('net_batch_size'),
        -module         => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::AlignmentNets',
        -parameters     => $self->o('net_parameters'),
        -can_be_empty   => 1,
        -rc_name        => '8GB_lastz',
      },


	    {
        -logic_name => 'remove_inconsistencies_after_net',
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::RemoveAlignmentDataInconsistencies',
        -flow_into => {
          1 => [ 'update_max_alignment_length_after_net' ],
        },
        -wait_for =>  [ 'alignment_nets', 'alignment_nets_himem', 'alignment_nets_hugemem', 'create_alignment_nets_jobs' ],    # Needed because of bi-directional netting: 2 jobs in create_alignment_nets_jobs can result in 1 job here
        -rc_name => '2GB_lastz',
      },


	    {
        -logic_name => 'create_filter_duplicates_net_jobs', #factory
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::CreateFilterDuplicatesJobs',
        -parameters => { },
        -wait_for =>  [ 'remove_inconsistencies_after_net' ],
        -flow_into => {
          2 => { 'filter_duplicates_net' => INPUT_PLUS() },
        },
        -can_be_empty  => 1,
        -rc_name => '2GB_lastz',
      },


	    {
        -logic_name   => 'filter_duplicates_net',
        -module        => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::FilterDuplicates',
        -parameters    => {
          'window_size' => $self->o('compara_window_size'),
          'filter_duplicates_net' => 1,
        },
        -hive_capacity => $self->o('filter_duplicates_hive_capacity'),
        -batch_size    => $self->o('filter_duplicates_batch_size'),
        -flow_into => {
          -1 => [ 'filter_duplicates_net_himem' ], # MEMLIMIT
        },
        -can_be_empty  => 1,
        -rc_name => $self->o('filter_duplicates_rc_name'),
      },


	    {
        -logic_name   => 'filter_duplicates_net_himem',
        -module        => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::FilterDuplicates',
        -parameters    => {
          'window_size' => $self->o('compara_window_size'),
          'filter_duplicates_net' => 1,
        },
       -hive_capacity => $self->o('filter_duplicates_hive_capacity'),
       -batch_size    => $self->o('filter_duplicates_batch_size'),
       -can_be_empty  => 1,
       -rc_name => $self->o('filter_duplicates_himem_rc_name'),
      },


	    {
        -logic_name => 'update_max_alignment_length_after_net',
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::GenomicAlignBlock::UpdateMaxAlignmentLength',
        -rc_name => '2GB_lastz',
        -wait_for =>  [ 'create_filter_duplicates_net_jobs', 'filter_duplicates_net', 'filter_duplicates_net_himem' ],
        -flow_into => [ 'set_internal_ids_collection' ],
      },

	    {
        -logic_name => 'set_internal_ids_collection',
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::SetInternalIdsCollection',
        -parameters => {
          'skip' => $self->o('patch_alignments'),
         },
         -flow_into => {
           2 => [ 'set_internal_ids_slow' ],
          },
         -analysis_capacity => 1,
      },


	    {
        -logic_name => 'set_internal_ids_slow',
        -module     => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::SetInternalIdsSlow',
        -analysis_capacity => 1,
        -rc_name => '8GB_lastz',
      },

	    {
        -logic_name => 'healthcheck',
        -module => 'Bio::EnsEMBL::Compara::RunnableDB::HealthCheck',
        -parameters => {
          'previous_db' => $self->o('previous_db'),
          'ensembl_release' => $self->o('ensembl_release'),
          'prev_release' => $self->o('prev_release'),
          'max_percent_diff' => $self->o('patch_alignments') ? $self->o('max_percent_diff_patches') : $self->o('max_percent_diff'),
        },
        -wait_for => [ 'set_internal_ids_collection' ],
        -rc_name => '2GB_lastz',
      },


	    {
        -logic_name => 'pairaligner_stats',
        -module => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::PairAlignerStats',
        -parameters => {
          'skip' => $self->o('skip_pairaligner_stats'),
          'dump_features' => $self->o('dump_features_exe'),
          'compare_beds' => $self->o('compare_beds_exe'),
          'create_pair_aligner_page' => $self->o('create_pair_aligner_page_exe'),
          'bed_dir' => $self->o('compara_bed_dir'),
          'ensembl_release' => $self->o('ensembl_release'),
          'output_dir' => $self->o('compara_feature_dir'),
        },
        -wait_for =>  [ 'healthcheck' ],
        -rc_name => '2GB_lastz',
      },

    ];
}


sub resource_classes {
  my $self = shift;

  return {
    '1GB' => { LSF => $self->lsf_resource_builder('production-rh74', 1000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '2GB_lastz' => { LSF => [$self->lsf_resource_builder('production-rh74', 2000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]), '-reg_conf '.$self->default_options->{compara_registry_file}]},
    '2GB' => { LSF => $self->lsf_resource_builder('production-rh74', 2000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '3GB' => { LSF => $self->lsf_resource_builder('production-rh74', 3000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '4GB_lastz' => { LSF => [$self->lsf_resource_builder('production-rh74', 4000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]), '-reg_conf '.$self->default_options->{compara_registry_file}]}, 
    '4GB' => { LSF => $self->lsf_resource_builder('production-rh74', 4000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '5GB' => { LSF => $self->lsf_resource_builder('production-rh74', 5000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '6GB' => { LSF => $self->lsf_resource_builder('production-rh74', 6000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '6GB_registry' => { LSF => [$self->lsf_resource_builder('production-rh74', 6000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]), '-reg_conf '.$self->default_options->{registry_file}]},
    '7GB' => { LSF => $self->lsf_resource_builder('production-rh74', 7000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '8GB_lastz' => { LSF => [$self->lsf_resource_builder('production-rh74', 8000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]), '-reg_conf '.$self->default_options->{compara_registry_file}]},
    '8GB' => { LSF => $self->lsf_resource_builder('production-rh74', 8000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '9GB' => { LSF => $self->lsf_resource_builder('production-rh74', 9000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '10GB' => { LSF => $self->lsf_resource_builder('production-rh74', 10000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '15GB_lastz' => { LSF => [$self->lsf_resource_builder('production-rh74', 15000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]), '-reg_conf '.$self->default_options->{compara_registry_file}]},
    '15GB' => { LSF => $self->lsf_resource_builder('production-rh74', 15000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '20GB' => { LSF => $self->lsf_resource_builder('production-rh74', 20000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '25GB' => { LSF => $self->lsf_resource_builder('production-rh74', 25000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '30GB' => { LSF => $self->lsf_resource_builder('production-rh74', 30000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '35GB' => { LSF => $self->lsf_resource_builder('production-rh74', 35000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '40GB' => { LSF => $self->lsf_resource_builder('production-rh74', 40000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '50GB' => { LSF => $self->lsf_resource_builder('production-rh74', 50000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '75GB' => { LSF => $self->lsf_resource_builder('production-rh74', 75000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '80GB' => { LSF => $self->lsf_resource_builder('production-rh74', 80000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '100GB' => { LSF => $self->lsf_resource_builder('production-rh74', 100000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'default' => { LSF => $self->lsf_resource_builder('production-rh74', 900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'repeatmasker' => { LSF => $self->lsf_resource_builder('production-rh74', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'repeatmasker_rebatch' => { LSF => $self->lsf_resource_builder('production-rh74', 5900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'simple_features' => { LSF => $self->lsf_resource_builder('production-rh74', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'genscan' => { LSF => $self->lsf_resource_builder('production-rh74', 3900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'genscan_short' => { LSF => $self->lsf_resource_builder('production-rh74', 5900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'blast' => { LSF => $self->lsf_resource_builder('production-rh74', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], undef, 3)},
    'blast10GB' => { LSF => $self->lsf_resource_builder('production-rh74', 10000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], undef, 3)},
    'blast_retry' => { LSF => $self->lsf_resource_builder('production-rh74', 5900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], undef, 3)},
    'project_transcripts' => { LSF => $self->lsf_resource_builder('production-rh74', 7200, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'projection_db_server'}, $self->default_options->{'projection_lastz_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'refseq_import' => { LSF => $self->lsf_resource_builder('production-rh74', 9900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'refseq_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
  }
}

sub hive_capacity_classes {
  my $self = shift;

  return {
           'hc_low'    => 200,
           'hc_medium' => 500,
           'hc_high'   => 1000,
         };
}


sub check_file_in_ensembl {
  my ($self, $file_path) = @_;
  push @{$self->{'_ensembl_file_paths'}}, $file_path;
  return $self->o('enscode_root_dir').'/'.$file_path;
}

1;
