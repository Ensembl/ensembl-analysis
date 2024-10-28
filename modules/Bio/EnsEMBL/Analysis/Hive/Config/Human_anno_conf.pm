=head1 LICENSE

Copyright [2021] EMBL-European Bioinformatics Institute

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

package Human_anno_conf;

use strict;
use warnings;
use File::Spec::Functions;

use Bio::EnsEMBL::ApiVersion qw/software_version/;
use base ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

sub default_options {
  my ($self) = @_;
  return {
    # inherit other stuff from the base class
    %{ $self->SUPER::default_options() },
    'reference_fasta'           => '',
    'dbowner'                   => '' || $ENV{EHIVE_USER} || $ENV{USER},
    'user'                      => '', # write db user
    'password'                  => '', # password for write db user
    'base_output_dir'           => '', # Where to write the files to
    'registry_file'             => catfile($self->o('base_output_dir'), 'Databases.pm'), # This needs to be a standard registry with the production/meta data/taxonomy db adaptors in there
    'pipeline_name'             => '' || $self->o('production_name').$self->o('production_name_modifier').'_'.$self->o('release_number'),
    'production_name'           => 'homo_sapiens' || $self->o('species_name'), # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
    'release_number'            => '',
    'ref_db_server'             => '', # host for dna db
    'ref_db_port'               => '',
    'ref_db_name'               => '',
    'user_r'                    => '', # read only db user
    'current_genebuild'            => 0,
    'assembly_accession'           => '', #the pipeline is initialed via standalone job  # Versioned GCA assembly accession, e.g. GCA_001857705.1
    'num_threads' => 20,
    'pipe_db_server'            => $ENV{GBS7}, # host for pipe db
    'dna_db_server'             => $ENV{GBS6}, # host for dna db
    'pipe_db_port'              => $ENV{GBP7}, # port for pipeline host
    'dna_db_port'               => $ENV{GBP6}, # port for dna db host
    'registry_db_server'        => $ENV{GBS1}, # host for registry db
    'registry_db_port'          => $ENV{GBP1}, # port for registry db
    'registry_db_name'          => 'gb_assembly_registry', # name for registry db
    'species_name'              => '', # e.g. mus_musculus
    'use_genome_flatfile'       => '1',# This will read sequence where possible from a dumped flatfile instead of the core db
    'production_name_modifier'  => '', # Do not set unless working with non-reference strains, breeds etc. Must include _ in modifier, e.g. _hni for medaka strain HNI
    species_division            => 'EnsemblVertebrates',
    strain_type                 => 'haplotype',
    initial_release_date        => '',
    source_assembly_name        => 'GRCh38',

########################
# Cut-offs
########################

    # cut-offs for problematic transcripts and anchors in LowDivergenceProjection
    coverage_cutoff                  => 98,
    perc_id_cutoff                   => 99,
    extended_length_variation_cutoff => 0.1,
    anchor_coverage_cutoff           => 0.95,
    anchor_perc_id_cutoff            => 0.99,

########################
# Pipe and ref db info
########################

    'provider_name'                => 'Ensembl',
    'provider_url'                 => 'www.ensembl.org',

    'pipe_db_name'                  => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_pipe_'.$self->o('release_number'),
    'dna_db_name'                   => $self->o('dbowner').'_'.$self->o('production_name').$self->o('production_name_modifier').'_core_'.$self->o('release_number'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'production_db_server'         => 'mysql-ens-meta-prod-1',
    'production_db_port'           => '4483',

    databases_to_delete => [],

######################################################
#
# Mostly constant settings
#
######################################################

    ensembl_analysis_script           => catdir($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts'),
    ensembl_analysis_script_genebuild => catdir($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'genebuild'),
    mapping_stats_script              => catfile($self->o('ensembl_analysis_script'), 'genebuild', 'calculate_remapping_stats.pl'),

    xy_scanner_path => catfile($self->o('ensembl_analysis_script'), 'pangenome', 'xy_scanner.py'),
    xy_scanner_data_path => '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/hprc/',
    x_marker_fasta_path => catfile($self->o('xy_scanner_data_path'), 'x_markers.fa'),
    y_marker_fasta_path => catfile($self->o('xy_scanner_data_path'), 'y_markers.fa'),

    ensembl_misc_script        => catdir($self->o('enscode_root_dir'), 'ensembl', 'misc-scripts'),
    repeat_types_script        => catfile($self->o('ensembl_misc_script'), 'repeats', 'repeat-types.pl'),
    meta_coord_script          => catfile($self->o('ensembl_misc_script'), 'meta_coord', 'update_meta_coord.pl'),
    meta_levels_script         => catfile($self->o('ensembl_misc_script'), 'meta_levels.pl'),
    frameshift_attrib_script   => catfile($self->o('ensembl_misc_script'), 'frameshift_transcript_attribs.pl'),
    select_canonical_script    => catfile($self->o('ensembl_misc_script'),'canonical_transcripts', 'select_canonical_transcripts.pl'),


########################
# Executable paths
########################
    minimap2_path => catfile($self->o('binary_base'), 'minimap2'),
    paftools_path => catfile($self->o('binary_base'), 'paftools.js'),
    samtools_path => catfile($self->o('binary_base'), 'samtools'), #You may need to specify the full path to the samtools binary

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# No option below this mark should be modified
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
########################################################
# URLs for retrieving the INSDC contigs and RefSeq files
########################################################
    'ncbi_base_ftp'           => 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all',
    'insdc_base_ftp'          => $self->o('ncbi_base_ftp').'/#expr(substr(#assembly_accession#, 0, 3))expr#/#expr(substr(#assembly_accession#, 4, 3))expr#/#expr(substr(#assembly_accession#, 7, 3))expr#/#expr(substr(#assembly_accession#, 10, 3))expr#/#assembly_accession#_#assembly_name#',
    'assembly_ftp_path'       => $self->o('insdc_base_ftp'),

########################
# db info
########################
    'core_db' => {
      -dbname => $self->o('dna_db_name'),
      -host   => $self->o('dna_db_server'),
      -port   => $self->o('dna_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'reference_db' => {
      -dbname => $self->o('ref_db_name'),
      -host   => $self->o('ref_db_server'),
      -port   => $self->o('ref_db_port'),
      -user   => $self->o('user_r'),
      -driver => $self->o('hive_driver'),
    },


    'production_db' => {
      -host   => $self->o('production_db_server'),
      -port   => $self->o('production_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -dbname => 'ensembl_production',
      -driver => $self->o('hive_driver'),
    },

    'taxonomy_db' => {
      -host   => $self->o('production_db_server'),
      -port   => $self->o('production_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -dbname => 'ncbi_taxonomy',
      -driver => $self->o('hive_driver'),
    },

    'registry_db' => {
      -host   => $self->o('registry_db_server'),
      -port   => $self->o('registry_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -dbname => $self->o('registry_db_name'),
      -driver => $self->o('hive_driver'),
    },

  };
}


sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
    wide_reference_fasta => $self->o('reference_fasta'),
    use_genome_flatfile => $self->o('use_genome_flatfile'),
  }
}

## See diagram for pipeline structure
sub pipeline_analyses {
    my ($self) = @_;

    return [


###############################################################################
#
# ASSEMBLY LOADING ANALYSES
#
###############################################################################
# 1) Process GCA - works out settings, flows them down the pipeline -> this should be seeded by another analysis later
# 2) Standard create core, populate tables, download data etc
# 3) Either run gbiab or setup gbiab
# 4) Finalise steps


      {
        # Creates a reference db for each species
        -logic_name => 'process_gca',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::ProcessGCA',
        -parameters => {
                         'num_threads'      => $self->o('num_threads'),
                         'dbowner'                     => $self->o('dbowner'),
                         'core_db'          => $self->o('core_db'),
                         'ensembl_release'  => $self->o('release_number'),
                         'base_output_dir'  => $self->o('base_output_dir'),
                         'registry_db'      => $self->o('registry_db'),
                         'registry_file'      => $self->o('registry_file'),
                         'current_genebuild'           => $self->o('current_genebuild'),
                         'assembly_accession'     =>$self->o('assembly_accession'),
                       },
        -rc_name    => '4GB',

        -flow_into  => {1 => ['create_core_db'],
          #                         1 => {'create_core_db' => {
          #                 assembly_accession => '#assembly_accession#',
          #                 assembly_name => '#assembly_name#',
          #                 core_db => '#core_db#',
          #                 output_path => '#output_path#',
          #                 stable_id_prefix => '#stable_id_prefix#',
          #                 species_url => '#species_url#',
          #                 species_name => '#species_name#',
          #                 species_display_name => '#species_display_name#',
          #                 production_name => '#production_name#',
          #                 toplevel_genome_file => '#toplevel_genome_file#',
          #                           reheadered_toplevel_genome_file => '#reheadered_toplevel_genome_file#',
          #                 core_dbname => '#core_dbname#',
          #                 stable_id_start => '#stable_id_start#',
          #                 }},
                       },
        -analysis_capacity => 1,
        -input_ids  => [],
      },


      {
        # Creates a reference db for each species
        -logic_name => 'create_core_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         'target_db'        => '#core_db#',
                         'enscode_root_dir' => $self->o('enscode_root_dir'),
                         'create_type'      => 'core_only',
                       },
        -rc_name    => '4GB',

        -flow_into  => {
                         1 => ['populate_production_tables'],
                       },
      },


      {
        # Load production tables into each reference
        -logic_name => 'populate_production_tables',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HivePopulateProductionTables',
        -parameters => {
                         'target_db'        => '#core_db#',
                         'output_path'      => '#output_path#',
                         'enscode_root_dir' => $self->o('enscode_root_dir'),
                         'production_db'    => $self->o('production_db'),
                       },
        -rc_name    => '4GB',

        -flow_into  => {
                         1 => ['process_assembly_info'],
                       },
      },


      {
        # Download the files and dir structure from the NCBI ftp site. Uses the link to a species in the ftp_link_file
        -logic_name => 'process_assembly_info',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProcessAssemblyReport',
        -parameters => {
          full_ftp_path => $self->o('assembly_ftp_path'),
          output_path   => '#output_path#',
          target_db     => '#core_db#',
        },
        -rc_name    => '4GB',
      	-max_retry_count => 0,
        -flow_into  => {
          1 => ['load_meta_info'],
        },
        -analysis_capacity => 20,
      },


      {
        # Load some meta info and seq_region_synonyms
        -logic_name => 'load_meta_info',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => '#core_db#',
          sql => [
            'INSERT IGNORE INTO meta (species_id, meta_key, meta_value) VALUES '.
              '(1, "annotation.provider_name", "'.$self->o('provider_name').'"),'.
              '(1, "annotation.provider_url", "'.$self->o('provider_url').'"),'.
              '(1, "assembly.coverage_depth", "high"),'.
              '(1, "assembly.provider_name", NULL),'.
              '(1, "assembly.provider_url", NULL),'.
              '(1, "assembly.ucsc_alias", NULL),'.
              '(1, "species.stable_id_prefix", "#stable_id_prefix#"),'.
              '(1, "species.url", "#species_url#"),'.
              '(1, "species.display_name", "#species_display_name#"),'.
              '(1, "species.division", "'.$self->o('species_division').'"),'.
              '(1, "species.strain", REPLACE("#assembly_name#", "-", "_")),'.
              '(1, "repeat.analysis", "repeatdetector"),' .
              '(1, "repeat.analysis", "dust"),' .
              '(1, "repeat.analysis", "trf"),' .
              '(1, "repeat.analysis", "repeatmask_repbase_human"),' .
              '(1, "species.production_name", "#production_name#"),'.
              '(1, "strain.type", "'.$self->o('strain_type').'"),'.
              '(1, "genebuild.initial_release_date", "'.$self->o('initial_release_date').'"),'.
              '(1, "genebuild.last_geneset_update", "'.$self->o('initial_release_date').'"),'.
              '(1, "genebuild.projection_source_db", "'.$self->o('ref_db_name').'"),'.
              '(1, "genebuild.id", '.$self->o('genebuilder_id').'),'.
              '(1, "genebuild.method", "projection_build"),'.
              '(1, "genebuild.method_display", "Projection from '.$self->o('source_assembly_name').'")',
          ],
        },
        -max_retry_count => 0,
        -rc_name    => '4GB',
        -flow_into  => {
          1 => ['load_taxonomy_info'],
        },
      },

      {
        -logic_name => 'load_taxonomy_info',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadTaxonomyInfo',
        -parameters => {
                         'target_db'        => '#core_db#',
                         'taxonomy_db'      => $self->o('taxonomy_db'),
                       },
        -rc_name    => '4GB',

        -flow_into  => {
                          1 => ['dump_toplevel_file'],#['load_windowmasker_repeats'],# 'fan_refseq_import'],
                       },
      },


      {
        -logic_name => 'dump_toplevel_file',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDumpGenome',
        -parameters => {
                         'coord_system_name'    => 'toplevel',
                         'target_db'            => '#core_db#',
                         'output_path'          => '#output_path#',
                         'enscode_root_dir'     => $self->o('enscode_root_dir'),
                         'species_name'         => '#species_name#',
                         'repeat_logic_names'   => [], # This is emtpy as we just use masking present in downloaded file
                       },
        -flow_into => {
          1 => ['reheader_toplevel_file'],
        },
        -analysis_capacity => 20,
        -rc_name    => '4GB',
      },


      {
        -logic_name => 'reheader_toplevel_file',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.catfile($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'genebuild', 'convert_genome_dump.pl').
                                      ' -conversion_type slice_name_to_seq_region_name'.
                                      ' -input_file #toplevel_genome_file#'.
                                      ' -output_file #reheadered_toplevel_genome_file#',
                       },
        -rc_name => '4GB',
        -flow_into => {
          1 => ['create_faidx'],
        },
     },


     {
        -logic_name => 'create_faidx',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -rc_name => '4GB',
        -parameters => {
          cmd => 'if [ ! -e "'.'#reheadered_toplevel_genome_file#'.'.fai" ]; then '.$self->o('samtools_path').' faidx '.'#reheadered_toplevel_genome_file#'.';fi',
        },

       -flow_into  => {
          1 => ['create_minimap2_index'],
        },

      },


      {
        -logic_name => 'create_minimap2_index',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'if [ ! -e "'.'#reheadered_toplevel_genome_file#.mmi'.'" ]; then '.$self->o('minimap2_path').
                 ' -d '.'#reheadered_toplevel_genome_file#.mmi'.' '.'#reheadered_toplevel_genome_file#'.';fi',
        },
        -flow_into  => {
          1 => ['check_index_not_empty'],
        },
        -rc_name => '20GB',
      },


      {
        -logic_name => 'check_index_not_empty',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'if [ -s "'.'#reheadered_toplevel_genome_file#.mmi'.'" ]; then exit 0; else exit 42;fi',
                         return_codes_2_branches => {'42' => 2},
        },
        -flow_into  => {
          1 => ['run_anno_repeats'],
        },
        -rc_name => '9GB',
      },

      {
        -logic_name => 'run_anno_repeats',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'python ' . catfile( $self->o('enscode_root_dir'), 'ensembl-anno', 'ensembl_anno.py' ) . ' #anno_repeats_commandline#',
      },
        -rc_name         => 'anno',
        -max_retry_count => 0,
        -flow_into       => {
        1 => ['xy_scanner']
      },
      },

      {
        -logic_name => 'xy_scanner',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::XYScanner',
        -parameters => {
                         xy_scanner_path     => $self->o('xy_scanner_path'),
                         target_genome_index => '#reheadered_toplevel_genome_file#'.'.mmi',
                         x_marker_fasta_path => $self->o('x_marker_fasta_path'),
                         y_marker_fasta_path => $self->o('y_marker_fasta_path'),
                         output_dir          => '#output_path#',
                       },
        -rc_name    => '9GB',
        -flow_into  => {
          '1->A' => ['create_remap_jobs'],
          'A->1' => ['map_remaining_genes'],
        },
      },


      {
        -logic_name => 'create_remap_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
          target_db           => $self->o('reference_db'),
          iid_type            => 'feature_id',
          feature_type        => 'gene',
          feature_restriction => 'readthrough',
          batch_size          => 100,
          id_output_file_path => '#output_path#/gene_ids_to_map.txt',
        },
        -rc_name    => '15GB',
        -flow_into => {
          2 => {'project_gene_batches' => {'xy_scanner' => '#xy_scanner#', 'core_db' => '#core_db#','genome_index' => '#reheadered_toplevel_genome_file#'.'.mmi','iid' => '#iid#'}},
        },
      },


      {
        -logic_name => 'project_gene_batches',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::LowDivergenceProjection',
        -parameters => {
                         genome_index   => '#genome_index#',
                         source_dna_db => $self->o('reference_db'),
                         source_gene_db => $self->o('reference_db'),
                         target_dna_db => '#core_db#',
                         target_gene_db => '#core_db#',
                         paftools_path => $self->o('paftools_path'),
                         minimap2_path => $self->o('minimap2_path'),
                         source_dna_fasta => '#wide_reference_fasta#',
                         disconnect_jobs => 1,
                         coverage_cutoff => $self->o('coverage_cutoff'),
                         perc_id_cutoff => $self->o('perc_id_cutoff'),
                         extended_length_variation_cutoff => $self->o('extended_length_variation_cutoff'),
                         anchor_coverage_cutoff => $self->o('anchor_coverage_cutoff'),
                         anchor_perc_id_cutoff => $self->o('anchor_perc_id_cutoff'),
                       },
        -rc_name    => '15GB',
        -hive_capacity => 700,
        -max_retry_count => 2,
      },


      {
        -logic_name => 'map_remaining_genes',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::Minimap2Remap',
        -parameters => {
                         genome_index   => '#reheadered_toplevel_genome_file#'.'.mmi',
                         logic_name => 'minimap2remap',
                         source_dna_fasta => '#wide_reference_fasta#',
                         source_dna_db => $self->o('reference_db'),
                         source_gene_db => $self->o('reference_db'),
                         target_dna_db => '#core_db#',
                         target_gene_db => '#core_db#',
                         paftools_path => $self->o('paftools_path'),
                         minimap2_path => $self->o('minimap2_path'),
                         input_id_file => '#output_path#/gene_ids_to_map.txt',
                       },
        -rc_name    => '35GB',
        -max_retry_count => 0,
        -hive_capacity => 700,
        -flow_into  => {
          1 => ['create_paralogue_jobs'],
	      },
      },


      {
        -logic_name => 'create_paralogue_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
          target_db           => '#core_db#',
          iid_type            => 'feature_id',
          feature_type        => 'gene',
          batch_size          => 1000,
        },
        -rc_name    => '15GB',
        -flow_into => {
          '2->A' => {'find_paralogues' => {'core_db' => '#core_db#','genome_file' => '#reheadered_toplevel_genome_file#','genome_index' => '#reheadered_toplevel_genome_file#'.'.mmi','iid' => '#iid#'}},
          'A->1' => ['collapse_paralogues'],
        },
      },


      {
        -logic_name => 'find_paralogues',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::FindRecentParalogues',
        -parameters => {
                         genome_index    => '#genome_index#',
                         genome_file     => '#reheadered_toplevel_genome_file#',
                         target_dna_db   => '#core_db#',
                         target_db  => '#core_db#',
                         paftools_path => $self->o('paftools_path'),
                         minimap2_path => $self->o('minimap2_path'),
                       },
        -rc_name    => '15GB',
        -hive_capacity => 700,
        -max_retry_count => 0,
      },


      {
        -logic_name => 'collapse_paralogues',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::CollapseParalogues',
        -parameters => {
                         genome_file     => '#reheadered_toplevel_genome_file#',
                         target_dna_db   => '#core_db#',
                         target_db  => '#core_db#',
                       },
        -rc_name    => '15GB',
        -max_retry_count => 0,
        -flow_into  => {
          1 => ['finalise_geneset'],
        },
      },


      {
        -logic_name => 'finalise_geneset',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::FinaliseRemappedGeneset',
        -parameters => {
                         genome_file     => '#reheadered_toplevel_genome_file#',
                         source_gene_db => $self->o('reference_db'),
                         target_dna_db   => '#core_db#',
                         target_db  => '#core_db#',
                       },
        -rc_name    => '15GB',
        -max_retry_count => 0,
        -flow_into  => {
          1 => ['final_cleaning'],
        },
      },

      {
        -logic_name => 'final_cleaning',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => '#core_db#',
          sql => [
            'DELETE exon FROM exon LEFT JOIN exon_transcript ON exon.exon_id = exon_transcript.exon_id WHERE exon_transcript.exon_id IS NULL',
            'TRUNCATE supporting_feature',
            'TRUNCATE transcript_supporting_feature',
            'TRUNCATE dna_align_feature',
            'UPDATE analysis SET logic_name="ensembl" WHERE logic_name="minimap2remap"',
            'UPDATE gene SET analysis_id = (SELECT analysis_id FROM analysis WHERE logic_name = "ensembl")'.
            ' WHERE analysis_id IN'.
            ' (SELECT analysis_id FROM analysis WHERE logic_name IN ("find_paralogues","collapse_paralogues"))',
            'DELETE FROM analysis WHERE logic_name IN ("find_paralogues","collapse_paralogues")',
            'DELETE FROM ad USING analysis_description ad LEFT JOIN analysis a ON ad.analysis_id = a.analysis_id WHERE a.analysis_id IS NULL',
            'UPDATE transcript JOIN gene USING(gene_id) SET transcript.analysis_id = gene.analysis_id',
            'UPDATE repeat_feature SET repeat_start = 1 WHERE repeat_start < 1',
            'UPDATE repeat_feature SET repeat_end = 1 WHERE repeat_end < 1',
            'UPDATE repeat_feature JOIN seq_region USING(seq_region_id) SET repeat_end = length WHERE repeat_end > length',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
                       1 => ['set_meta_coords'],
        },

      },



      {
        -logic_name => 'set_meta_coords',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('meta_coord_script').
                                ' -user '.$self->o('user').
                                ' -pass '.$self->o('password').
                                ' -host '.$self->o('core_db','-host').
                                ' -port '.$self->o('core_db','-port').
                                ' -dbpattern '.'#core_dbname#'
                       },
        -rc_name => 'default',
        -flow_into => {
                        1 => ['set_meta_levels'],
                      },
      },


      {
        -logic_name => 'set_meta_levels',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('meta_levels_script').
                                ' -user '.$self->o('user').
                                ' -pass '.$self->o('password').
                                ' -host '.$self->o('core_db','-host').
                                ' -port '.$self->o('core_db','-port').
                                ' -dbname '.'#core_dbname#'
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['set_frameshift_introns'] },
      },


      {
        -logic_name => 'set_frameshift_introns',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('frameshift_attrib_script').
                                ' -user '.$self->o('user').
                                ' -pass '.$self->o('password').
                                ' -host '.$self->o('core_db','-host').
                                ' -port '.$self->o('core_db','-port').
                                ' -dbpattern '.'#core_dbname#'
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['set_canonical_transcripts'] },
      },


      {
        -logic_name => 'set_canonical_transcripts',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('select_canonical_script').
                                ' -dbuser '.$self->o('user').
                                ' -dbpass '.$self->o('password').
                                ' -dbhost '.$self->o('core_db','-host').
                                ' -dbport '.$self->o('core_db','-port').
                                ' -dbname '.'#core_dbname#'.
                                ' -coord toplevel -write'
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['null_columns'] },
      },

      {
        -logic_name => 'null_columns',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => '#core_db#',
          sql => [
            'UPDATE gene SET stable_id = NULL',
            'UPDATE transcript SET stable_id = NULL',
            'UPDATE translation SET stable_id = NULL',
            'UPDATE exon SET stable_id = NULL',
            'UPDATE protein_align_feature set external_db_id = NULL',
            'UPDATE dna_align_feature set external_db_id = NULL',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['run_stable_ids'],
                      },
      },


      {
        -logic_name => 'run_stable_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::SetStableIDs',
        -parameters => {
                         enscode_root_dir => $self->o('enscode_root_dir'),
                         mapping_required => 0,
                         target_db => '#core_db#',
                         id_start => '#stable_id_prefix#'.'#stable_id_start#',
                         output_path => '#output_path#',
                       },
        -rc_name    => 'default',
        -flow_into => {
                         1 => ['generate_mapping_stats'],
                      },
      },


      {
        -logic_name => 'generate_mapping_stats',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('mapping_stats_script').
                                ' -xy_scanner '.'#xy_scanner#'.
                                ' -query_user '.$self->o('user_r').
                                ' -query_host '.$self->o('core_db','-host').
                                ' -query_port '.$self->o('core_db','-port').
                                ' -query_dbname '.'#core_dbname#'.
                                ' -reference_user '.$self->o('user_r').
                                ' -reference_host '.$self->o('ref_db_server').
                                ' -reference_port '.$self->o('ref_db_port').
                                ' -reference_dbname '.$self->o('ref_db_name').
                                ' -output_dir '.'#output_path#'.
                                ' -output_file_prefix '.'#assembly_accession#'."_mapping_stats"
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['add_placeholder_sample_location'] },
      },


      {
        -logic_name => 'add_placeholder_sample_location',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAddPlaceholderLocation',
        -parameters => {
                        input_db => '#core_db#',
                       },
        -rc_name    => 'default',
        -flow_into => {
                       1 => ['populate_analysis_descriptions'],
        },
      },

      {
        -logic_name => 'populate_analysis_descriptions',
        -module     => 'Bio::EnsEMBL::Production::Pipeline::ProductionDBSync::PopulateAnalysisDescription',
        -parameters => {
                        species => '#production_name#',
                        group => 'core',
                       },
        -rc_name    => '4GB_registry',
        -flow_into => {
                       1 => ['run_data_checks'],
        },
      },

      {
        -logic_name      => 'run_data_checks',
        -module          => 'Bio::EnsEMBL::DataCheck::Pipeline::RunDataChecks',
        -parameters      => {
          datacheck_groups => ['core'],
          failures_fatal  => 1,
          output_file     => catfile('#output_path#', '#production_name#_dc.log'),
          history_file     => catfile('#output_path#', '#production_name#_dc.json'),
          registry_file   => $self->o('registry_file'),
          species         => '#production_name#',
        },
        -max_retry_count => 1,
        -hive_capacity   => 50,
        -batch_size      => 10,
        -rc_name         => '4GB_registry',
        -flow_into => {
                       1 => {'create_dump_dir' => {core_dbname => '#core_dbname#', production_name => '#production_name#', output_path => '#output_path#'}},
        },
      },

      {
        -logic_name => 'create_dump_dir',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'mkdir -p #output_dir#',
          output_dir => catdir('#output_path#', 'GFF3'),
        },
        -rc_name => '4GB_registry',
        -flow_into => {
          1 => ['dump_gff'],
        },
      },

      {
        -logic_name  => 'dump_gff',
        -module      => 'Bio::EnsEMBL::Production::Pipeline::FileDump::Geneset_GFF3',
        -parameters  => {
          output_dir => catdir('#output_path#', 'GFF3'),
          species    => '#production_name#',
          custom_filenames => {
            genes => catfile('#output_dir#', '#core_dbname#.gff3'),
          }
        },
        -max_retry_count => 1,
        -hive_capacity   => 50,
        -rc_name         => '4GB_registry',
      },

    ];
}


sub resource_classes {
  my $self = shift;
  return {
    'default' => {
      LSF => $self->lsf_resource_builder('production', 3000),
      SLURM => $self->slurm_resource_builder(3000,'7-00:00:00'),
    },
    '4GB_registry' => {
      LSF => [$self->lsf_resource_builder('production', 4000), '-reg_conf '.$self->o('registry_file')],
      SLURM => [ $self->slurm_resource_builder( 4000, '1-00:00:00', undef), ' -reg_conf ' . $self->o('registry_file')]
    },
    '4GB' => {
      LSF => $self->lsf_resource_builder('production', 4000),
      SLURM => $self->slurm_resource_builder(4000,'7-00:00:00'),
    },
    '5GB' => {
      LSF => $self->lsf_resource_builder('production', 5000),
      SLURM => $self->slurm_resource_builder(5000,'7-00:00:00'),
    },
    '9GB' => {
      LSF => $self->lsf_resource_builder('production', 9000),
      SLURM => $self->slurm_resource_builder(9000,'7-00:00:00'),
    },
    '12GB' => {
      LSF => $self->lsf_resource_builder('production', 12000),
      SLURM => $self->slurm_resource_builder(12000,'7-00:00:00'),
    },
    '15GB' => {
      LSF => $self->lsf_resource_builder('production', 15000),
      SLURM => $self->slurm_resource_builder(15000,'7-00:00:00'),
    },
    '20GB' => {
      LSF => $self->lsf_resource_builder('production', 20000),
      SLURM => $self->slurm_resource_builder(20000,'7-00:00:00'),
    },
    '35GB' => {
      LSF => $self->lsf_resource_builder('production', 35000),
      SLURM => $self->slurm_resource_builder(35000,'7-00:00:00'),
    },
    'anno'             => {
     LSF => $self->lsf_resource_builder( 'production', 50000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ], $self->default_options->{'num_threads'} ),
     SLURM =>  $self->slurm_resource_builder(50000, '7-00:00:00', $self->default_options->{'num_threads'} ),
     }
  };
}

1;
