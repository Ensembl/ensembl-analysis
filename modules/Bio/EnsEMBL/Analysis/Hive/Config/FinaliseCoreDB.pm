
=Head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2024] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Analysis::Hive::Config::FinaliseCoreDB;

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
    'user_r'             => '',                                                                                                # read only db user
    'user'               => '',                                                                                                # write db user
    'password'           => '',                                                                                                # password for write db user
    'pipe_db_host'       => '',                                                                                                # host for pipe db
    'databases_host'     => '',                                                                                                # host for general output dbs
    'dna_db_host'        => '',                                                                                                # host for dna db
    'pipe_db_port'       => '',                                                                                                # port for pipeline host
    'databases_port'     => '',                                                                                                # port for general output db host
    'dna_db_port'        => '',                                                                                                # port for dna db host
    'registry_host'      => '',                                                                                                # host for registry db
    'registry_port'      => '',                                                                                                # port for registry db
    'registry_db'        => '',                                                                                                # registry DB (active annotation tracking registry)
    'species_name'       => '',                                                                                                # e.g. mus_musculus
    'production_name'    => '',                                                                                                # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
    dbname_accession          => '', # This is the assembly accession without [._] and all lower case, i.e gca001857705v1
    'release_number'     => '' || $self->o('ensembl_release'),
    'uniprot_set'        => '',                                                                                                # e.g. mammals_basic, check UniProtCladeDownloadStatic.pm module in hive config dir for suitable set,
    'sanity_set'         => '',                                                                                                # sanity checks
    'output_path'        => '',# Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses
    'assembly_name'      => '',                                                                                                # Name (as it appears in the assembly report file)
    'assembly_accession' => '',                                                                                                # Versioned GCA assembly accession, e.g. GCA_001857705.1
    'taxon_id'  => '', # should be in the assembly report file
    'stable_id_prefix'   => '',                                                                                                # e.g. ENSPTR. When running a new annotation look up prefix in the assembly registry db
    'stable_id_start'    => '0',                                                                                               # When mapping is not required this is usually set to 0
    'skip_projection'    => '0',                                                                                               # Will skip projection process if 1
    'skip_rnaseq'        => '0',                                                                                               # Will skip rnaseq analyses if 1
    'mapping_required'   => '0',                                                                                               # If set to 1 this will run stable_id mapping sometime in the future. At the moment it does nothing
    'mapping_db'         => '',                                                                                                # Tied to mapping_required being set to 1, we should have a mapping db defined in this case, leave undef for now
    'uniprot_version'    => 'uniprot_2021_04',                                                                                 # What UniProt data dir to use for various analyses
    'protein_entry_loc'  => catfile( $self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), 'entry_loc' ),    # Used by genscan blasts and optimise daf/paf. Don't change unless you know what you're doing
    'registry_file'      => catfile($self->o('output_path'), 'Databases.pm'), # Path to databse registry for LastaZ and Production sync
    'gst_dir'            => catfile($self->o('output_path'), 'gst'),
    'busco_singularity_image'      => '/hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/busco_v5.8.2_cv1.sif',
    'busco_download_path'          => '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/data/busco_data/data_odb12/',
    'busco_genome_dir'   => catdir($self->o('output_path'), 'busco_genome'),
    'busco_protein_dir'   => catdir($self->o('output_path'), 'busco_protein'),
    'busco_datasets_file'  => catfile( $self->o('enscode_root_dir'), 'ensembl-genes', 'src', 'python', 'ensembl', 'genes', 'metrics', 'busco_lineage.json'),
    'busco_threshold' => 90, # If the busco score is above this threshold, the pre-release files will be produced
    'busco_lower_threshold' => 80, # If the busco score is above this threshod and the difference less than 'busco_difference_threshold', the pre-release files will be produced
    'busco_difference_threshold' => 10, # If the difference between the gene and protein busco score is less than this value, the pre-release files will be produced as long as the busco score is above 'busco_lower_threshold'
     'ftp_release'            => catfile($self->o('output_path'), 'ftp_release'),

########################
# Pipe and ref db info
########################

    'projection_source_db_name'         => '',                                                                                 # This is generally a pre-existing db, like the current human/mouse core for example

    'pipe_db_name' => $self->o('dbowner') . '_' . $self->o('dbname_accession') . '_pipe_' . $self->o('release_number'),
    'dna_db_name'  => $self->o('dbowner') . '_' . $self->o('dbname_accession') . '_core_' . $self->o('release_number'),

    'reference_db_name'   => $self->o('dna_db_name'),
    'reference_db_host'   => $self->o('dna_db_host'),
    'reference_db_port'   => $self->o('dna_db_port'),

    'final_geneset_db_host'   => $self->o('databases_host'),
    'final_geneset_db_port'   => $self->o('databases_port'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'ensembl_release'      => $ENV{ENSEMBL_RELEASE},     # this is the current release version on staging to be able to get the correct database
    'production_db_host'   => 'mysql-ens-meta-prod-1',
    'production_db_port'   => '4483',
    'production_ftp_dir' => '/nfs/ftp/public/databases/ensembl/pre-release',

##############################
### Set rnaseq file paths
################################
'rnaseq_dir' => catdir( $self->o('output_path'), 'rnaseq' ),
'merge_dir'  => catdir( $self->o('rnaseq_dir'),  'merge' ),
'output_dir' => catdir( $self->o('rnaseq_dir'),  'output' ),
#################################


########################
# BLAST db paths
########################
    'base_blast_db_path' => $ENV{BLASTDB_DIR},

    ensembl_analysis_script  => catdir( $self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts' ),
    load_optimise_script     => catfile( $self->o('ensembl_analysis_script'), 'genebuild', 'load_external_db_ids_and_optimize_af.pl' ),
    ensembl_misc_script      => catdir( $self->o('enscode_root_dir'), 'ensembl', 'misc-scripts' ),
    meta_coord_script        => catfile( $self->o('ensembl_misc_script'), 'meta_coord', 'update_meta_coord.pl' ),
    meta_levels_script       => catfile( $self->o('ensembl_misc_script'),     'meta_levels.pl' ),
    frameshift_attrib_script => catfile( $self->o('ensembl_misc_script'),     'frameshift_transcript_attribs.pl' ),
    select_canonical_script  => catfile( $self->o('ensembl_misc_script'),     'canonical_transcripts', 'select_canonical_transcripts.pl' ),
    assembly_name_script     => catfile( $self->o('ensembl_analysis_script'), 'update_assembly_name.pl' ),
    core_metadata_script     => catdir( $self->o('enscode_root_dir'), 'ensembl-genes', 'src', 'python', 'ensembl', 'genes', 'metadata', 'core_meta_data.py'),
    core_stats_script        => catdir( $self->o('enscode_root_dir'), 'ensembl-genes', 'src', 'perl', 'ensembl', 'genes', 'generate_species_homepage_stats.pl'),	
    ensembl_gst_script       => catdir( $self->o('enscode_root_dir'), 'ensembl-genes', 'pipelines' , 'gene_symbol_classifier'),   
    gst_dump_proteins_script => catfile( $self->o('ensembl_gst_script'), 'dump_protein_sequences.pl' ),
    gst_load_symbols_script  => catfile( $self->o('ensembl_gst_script'), 'load_gene_symbols.pl' ),
    busco_lineage_selector_script => catfile( $self->o('enscode_root_dir'), 'ensembl-genes', 'src', 'python', 'ensembl', 'genes', 'metrics', 'busco_lineage_selector.py'),
    busco_metakeys_patch_script   => catfile( $self->o('enscode_root_dir'), 'ensembl-genes', 'src', 'python', 'ensembl', 'genes', 'metrics', 'busco_metakeys_patch.py'),
    registry_status_update_python_script => catfile( $self->o('enscode_root_dir'), 'ensembl-genes', 'src', 'python', 'ensembl', 'genes', 'info_from_registry', 'update_assembly_registry.py' ),
    write_metrics_to_registry_script     => catfile( $self->o('enscode_root_dir'), 'ensembl-genes', 'src', 'python', 'ensembl', 'genes', 'info_from_registry', 'write_metrics_to_registry.py' ),
    pre_release_ftp_script            => catdir( $self->o('enscode_root_dir'), 'ensembl-genes', 'src', 'python', 'ensembl', 'genes', 'automation', 'pre_release_ftp.py' ),
 
# Genes biotypes to ignore from the final db when copying to core
    copy_biotypes_to_ignore => {
      'low_coverage' => 1,
      'CRISPR'       => 1,
      broken_gene    => 1, # This is a biotype to quickly skip bad genes
    },

########################
# Extra db settings
########################
    mysql_dump_options => '--max_allowed_packet=1000MB',
    desired_slice_length => 10000000,
    store_rejected => 0,

########################
# db info
########################
    'reference_db' => {
      -dbname => $self->o('reference_db_name'),
      -host   => $self->o('reference_db_host'),
      -port   => $self->o('reference_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'final_geneset_db' => {
      -dbname => $self->o('dbowner') . '_' . $self->o('dbname_accession') . '_final_' . $self->o('release_number'),
      -host   => $self->o('final_geneset_db_host'),
      -port   => $self->o('final_geneset_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'production_db' => {
      -host   => $self->o('production_db_host'),
      -port   => $self->o('production_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -dbname => 'ensembl_production',
      -driver => $self->o('hive_driver'),
    },

########################
# FTP Dump 
########################
    ## gff3 & gtf parameter
    'abinitio'               => 1,
    'gene'                   => 1,

    ## gtf parameters, e! specific
    'gtftogenepred_exe'      => 'gtfToGenePred',
    'genepredcheck_exe'      => 'genePredCheck',

    ## gff3 parameters
    'gt_exe'                 => 'gt',
    'gff3_tidy'              => $self->o('gt_exe') . ' gff3 -tidy -sort -retainids -fixregionboundaries -force',
    'gff3_validate'          => $self->o('gt_exe') . ' gff3validator',

    'feature_type'           => [ 'Gene', 'Transcript', 'SimpleFeature' ], #'RepeatFeature'
    'per_chromosome'         => 1,
    'include_scaffold'       => 1,
    'logic_name'             => [],
    'db_type'                => 'core',
    'out_file_stem'          => undef,
    'xrefs'                  => 0,
########################
########################

  };
}

## See diagram for pipeline structure
sub pipeline_analyses {
  my ($self) = @_;

  return [

    {
      -logic_name => 'create_toplevel_slices',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db             => $self->o('final_geneset_db'),
        iid_type              => 'slice',
        coord_system_name     => 'toplevel',
        include_non_reference => 0,
        top_level             => 1,
      },
      -input_ids  => [{}],
      -flow_into => {
        '2->A' => ['split_slices_on_intergenic'],
        'A->1' => ['add_missing_analyses'],
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'split_slices_on_intergenic',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db             => $self->o('final_geneset_db'),
        iid_type              => 'intergenic_slice',
        coord_system_name     => 'toplevel',
        include_non_reference => 0,
        top_level             => 1,
        desired_slice_length  => $self->o('desired_slice_length'),
      },
      -batch_size => 300,
      -flow_into => {
        2 => ['clean_utr'],
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'clean_utr',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::CleanUTRs',
      -parameters => {
        source_db => $self->o('final_geneset_db'),
        target_db => $self->o('reference_db'),
        dna_db    => $self->o('dna_db'),
        store_rejected => $self->o('store_rejected'),
        copy_biotypes_to_ignore => $self->o('copy_biotypes_to_ignore'),
      },
      -rc_name => '4GB',
      -max_retry_count => 1,
      -hive_capacity => $self->o('hc_normal'),
    },

    {
      -logic_name => 'add_missing_analyses',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('reference_db'),
        sql     => [
          'INSERT IGNORE INTO analysis (created, logic_name) VALUES (NOW(), "rnaseq_intron_support")',
          'INSERT IGNORE INTO analysis (created, logic_name) VALUES (NOW(), "cdna_alignment_core")',
          'INSERT IGNORE INTO analysis (created, logic_name, db) VALUES (NOW(), "other_protein", "uniprot")',
          'INSERT IGNORE INTO analysis (created, logic_name, db) VALUES (NOW(), "projected_transcript", "'.$self->o('projection_source_db_name').'")',
          'INSERT IGNORE INTO analysis (logic_name, db_version, db_file, program_file, module) VALUES ("rfamcmsearch", "14.0", "'.$self->o('output_path').'", "cmsearch", "HiveCMSearch")',
        ],
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['update_biotypes_and_analyses'],
      },
    },

    {
      -logic_name => 'update_biotypes_and_analyses',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('reference_db'),
        sql     => [
          'UPDATE analysis SET logic_name = "cdna2genome" WHERE logic_name = "best_targetted"',
          'UPDATE gene SET biotype = "protein_coding" WHERE biotype = "ensembl"',
          'UPDATE gene SET biotype = "vault_RNA" WHERE biotype IN ("Vault_RNA","vaultRNA")',
          'UPDATE gene SET analysis_id = (SELECT analysis_id FROM analysis WHERE logic_name = "ensembl")'.
            ' WHERE analysis_id NOT IN (SELECT analysis_id FROM analysis WHERE logic_name IN ("ncrna", "mt_genbank_import"))',
          'UPDATE transcript JOIN gene USING(gene_id) SET transcript.biotype = gene.biotype',
          'UPDATE transcript JOIN gene USING(gene_id) SET transcript.analysis_id = gene.analysis_id',
          'UPDATE protein_align_feature SET analysis_id =' .
            '(SELECT analysis_id FROM analysis WHERE logic_name = "projected_transcript") WHERE analysis_id IN ' .
            '(SELECT analysis_id FROM analysis WHERE logic_name IN ("project_transcripts","cesar"))',
          'UPDATE protein_align_feature SET analysis_id =' .
            '(SELECT analysis_id FROM analysis WHERE logic_name = "other_protein") WHERE analysis_id NOT IN' .
            '(SELECT analysis_id FROM analysis WHERE logic_name IN ("uniprot","projected_transcript"))',
          'UPDATE dna_align_feature SET analysis_id =' .
            '(SELECT analysis_id FROM analysis WHERE logic_name = "projected_transcript") WHERE analysis_id IN' .
            '(SELECT analysis_id FROM analysis WHERE logic_name IN ("project_lincrna","project_pseudogene"))',
          'UPDATE dna_align_feature SET analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name = "rfamcmsearch")'.
            ' WHERE analysis_id = (SELECT analysis_id FROM analysis WHERE logic_name = "ncrna")',
          'UPDATE dna_align_feature SET analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name = "cdna_alignment_core")'.
            ' WHERE analysis_id IN (SELECT analysis_id FROM analysis WHERE logic_name IN ("exonerate","cdna2genome","best_targetted"))',
	    'UPDATE intron_supporting_evidence SET analysis_id = (SELECT analysis_id FROM analysis WHERE logic_name = "rnaseq_intron_support")',
	    'UPDATE repeat_feature SET repeat_start = 1 WHERE repeat_start < 1',
	    'UPDATE repeat_feature SET repeat_end = 1 WHERE repeat_end < 1',
	    'UPDATE repeat_feature JOIN seq_region USING(seq_region_id) SET seq_region_end = length WHERE seq_region_end > length',
        ],
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['null_columns'],
      },
    },

    {
      -logic_name => 'null_columns',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('reference_db'),
        sql     => [
          'UPDATE gene SET stable_id = NULL',
          'UPDATE transcript SET stable_id = NULL',
          'UPDATE translation SET stable_id = NULL',
          'UPDATE exon SET stable_id = NULL',
          'UPDATE protein_align_feature SET external_db_id = NULL',
          'UPDATE dna_align_feature SET external_db_id = NULL',
        ],
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['clean_unused_analyses'],
      },
    },

    {
      -logic_name => 'clean_unused_analyses',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('reference_db'),
        sql     => [
          'DELETE FROM analysis WHERE logic_name IN'.
            ' ("spliced_elsewhere","pseudogenes","genblast","genblast_not_best","project_pseudogene","blast_long_read",'.
            ' "project_lincrna","project_transcripts","ig_tr_collapse", "exonerate", "cdna2genome", "best_targetted",'.
            ' "filter_lncrnas", "blast", "process_homology_selenocysteine")',
          'DELETE from analysis where logic_name like "%\_rnaseq\_%"',
        ],
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['xref_cleaning'],
      },
    },

    {
      -logic_name => 'xref_cleaning',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('reference_db'),
        sql     => [
          'TRUNCATE associated_xref',
          'TRUNCATE dependent_xref',
          'TRUNCATE identity_xref',
          'TRUNCATE object_xref',
          'TRUNCATE ontology_xref',
          'TRUNCATE xref',
          'UPDATE gene SET display_xref_id = NULL',
          'UPDATE transcript SET display_xref_id = NULL',
        ],
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['set_meta_coords'],
      },
    },

    {
      -logic_name => 'set_meta_coords',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . $self->o('meta_coord_script') .
          ' -user ' . $self->o('user') .
          ' -pass ' . $self->o('password') .
          ' -host ' . $self->o( 'reference_db', '-host' ) .
          ' -port ' . $self->o( 'reference_db', '-port' ) .
          ' -dbpattern ' . $self->o( 'reference_db', '-dbname' )
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['set_meta_levels'],
      },
    },

    {
      -logic_name => 'set_meta_levels',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . $self->o('meta_levels_script') .
          ' -user ' . $self->o('user') .
          ' -pass ' . $self->o('password') .
          ' -host ' . $self->o( 'reference_db', '-host' ) .
          ' -port ' . $self->o( 'reference_db', '-port' ) .
          ' -dbname ' . $self->o( 'reference_db', '-dbname' )
      },
      -rc_name => 'default',
      -flow_into => { 1 => ['set_frameshift_introns'] },
    },

    {
      -logic_name => 'set_frameshift_introns',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . $self->o('frameshift_attrib_script') .
          ' -user ' . $self->o('user') .
          ' -pass ' . $self->o('password') .
          ' -host ' . $self->o( 'reference_db', '-host' ) .
          ' -port ' . $self->o( 'reference_db', '-port' ) .
          ' -dbpattern ' . $self->o( 'reference_db', '-dbname' )
      },
      -rc_name => '10GB',
      -flow_into => { 1 => ['set_canonical_transcripts'] },
    },

    {
      -logic_name => 'set_canonical_transcripts',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . $self->o('select_canonical_script') .
          ' -dbuser ' . $self->o('user') .
          ' -dbpass ' . $self->o('password') .
          ' -dbhost ' . $self->o( 'reference_db', '-host' ) .
          ' -dbport ' . $self->o( 'reference_db', '-port' ) .
          ' -dbname ' . $self->o( 'reference_db', '-dbname' ) .
          ' -coord toplevel -write'
      },
      -rc_name => '10GB',
      -flow_into => { 1 => ['run_stable_ids'] },
    },

    {
      -logic_name => 'run_stable_ids',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::SetStableIDs',
      -parameters => {
        enscode_root_dir => $self->o('enscode_root_dir'),
        mapping_required => $self->o('mapping_required'),
        target_db        => $self->o('reference_db'),
        mapping_db       => $self->o('mapping_db'),
        id_start         => $self->o('stable_id_prefix') . $self->o('stable_id_start'),
        output_path      => $self->o('output_path'),
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['backup_core_db_pre_optimise'],
      },
    },

    {
      #Creates a reference db for each species
      -logic_name => 'backup_core_db_pre_optimise',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::DatabaseDumper',
      -parameters => {
        src_db_conn  => $self->o('dna_db'),
        output_file  => catfile( $self->o('output_path'), 'core_post_stable_idsbak.sql.gz' ),
        dump_options => $self->o('mysql_dump_options'),
      },
      -rc_name => 'default',
      -flow_into => { 1 => ['load_external_db_ids_and_optimise_af_tables'] },
    },

    {
      -logic_name => 'load_external_db_ids_and_optimise_af_tables',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . $self->o('load_optimise_script') .
          ' -output_path ' . catdir( $self->o('output_path'), 'optimise' ) .
          ' -uniprot_filename ' . $self->o('protein_entry_loc') .
          ' -dbuser ' . $self->o('user') .
          ' -dbpass ' . $self->o('password') .
          ' -dbport ' . $self->o( 'reference_db', '-port' ) .
          ' -dbhost ' . $self->o( 'reference_db', '-host' ) .
          ' -dbname ' . $self->o( 'reference_db', '-dbname' ) .
          ' -prod_dbuser ' . $self->o('user_r') .
          ' -prod_dbhost ' . $self->o( 'production_db', '-host' ) .
          ' -prod_dbname ' . $self->o( 'production_db', '-dbname' ) .
          ' -prod_dbport ' . $self->o( 'production_db', '-port' ) .
          ' -ise -core'
      },
      -max_retry_count => 0,
      -rc_name         => '8GB',
      -flow_into       => {
        1 => ['drop_backup_tables_job'],
      },
    },


    {
      -logic_name => 'drop_backup_tables_job',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        inputquery   => 'SHOW TABLES LIKE "%bak%"',
        column_names => ['table'],
        db_conn      => $self->o('reference_db'),
      },
      -rc_name   => 'default',
      -flow_into => {
        '2->A' => ['drop_backup_tables'],
        'A->1' => ['final_meta_updates'],
      },
    },

    {
      -logic_name => 'drop_backup_tables',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        sql     => 'DROP TABLE #table#',
        db_conn => $self->o('reference_db'),
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'final_meta_updates',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('reference_db'),
        sql     => [
          'INSERT INTO meta (species_id, meta_key, meta_value) VALUES ' .
            '(1, "genebuild.last_geneset_update", (SELECT CONCAT((EXTRACT(YEAR FROM now())),"-",(LPAD(EXTRACT(MONTH FROM now()),2,"0"))))),'.
            '(1, "genebuild.method", "full_genebuild"),'.
            '(1, "genebuild.method_display", "Ensembl Genebuild"),'.
            '(1, "species.annotation_source", "ensembl")',
        ],
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['gst_dump_protein_sequences'],
      },
    },

    {
      -logic_name => 'gst_dump_protein_sequences',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
	  cmd => 'perl ' . $self->o('gst_dump_proteins_script') . ' --group core --species ' . $self->o('production_name') . ' --registry ' . $self->o('registry_file') . ' --output_file ' . $self->o('gst_dir') . '/' . $self->o('production_name') . '_protein_sequences.fa',
      },
	  -rc_name => 'default_registry',
	  -flow_into       => { 1 => ['gst_assign_gene_symbols'], },

    },

    {
      -logic_name => 'gst_assign_gene_symbols',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
	  cmd => 'singularity run --bind /hps/software/users/ensembl/genebuild/gene_symbol_classifier/data:/app/checkpoints --bind ' . $self->o('gst_dir') . ':/app/data /hps/software/users/ensembl/genebuild/gene_symbol_classifier/singularity/gene_symbol_classifier_0.12.1.sif --checkpoint /app/checkpoints/mlp_10_min_frequency_2022-01-29_03.08.32.ckpt --sequences_fasta /app/data/' . $self->o('production_name') .  '_protein_sequences.fa --scientific_name ' . $self->o('species_name'),
      },
	  -rc_name => 'default_registry',
	  -flow_into       => { 1 => ['gst_filter_assignments'], },
    },

    {
     -logic_name => 'gst_filter_assignments',
     -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
     -parameters => {
	 cmd => 'singularity run --bind  ' . $self->o('gst_dir') . ':/app/data /hps/software/users/ensembl/genebuild/gene_symbol_classifier/singularity/gene_symbol_classifier_filter_0.3.0.sif --symbol_assignments /app/data/' . $self->o('production_name') . '_protein_sequences_symbols.csv --threshold 0.7',
     },
	 -rc_name => 'default_registry',
	 -flow_into       => { 1 => ['gst_load_gene_symbols'], },	 
    },

    {
     -logic_name => 'gst_load_gene_symbols',
     -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
     -parameters => {
	 cmd => 'perl ' . $self->o('gst_load_symbols_script') . ' --species ' . $self->o('production_name') . ' --group core --registry ' . $self->o('registry_file') . ' --symbol_assignments ' .  $self->o('gst_dir') . '/' . $self->o('production_name') . '_protein_sequences_symbols_filtered.csv --primary_ids_file /hps/software/users/ensembl/genebuild/gene_symbol_classifier/data/display_name_dbprimary_acc_105.dat --program_version 0.12.1',
     },
	 -rc_name => 'default_registry',
	 -flow_into       => { 1 => ['add_placeholder_sample_location'], },
    },
      
    {
      -logic_name => 'add_placeholder_sample_location',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAddPlaceholderLocation',
      -parameters => {
        input_db => $self->o('reference_db'),
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['populate_analysis_descriptions'],
      },
    },

    {
      -logic_name => 'populate_analysis_descriptions',
      -module     => 'Bio::EnsEMBL::Production::Pipeline::ProductionDBSync::PopulateAnalysisDescription',
      -parameters => {
        species => $self->o('production_name'),
        group   => 'core',
      },
      -rc_name   => 'default_registry',
      -flow_into => {
        1 => ['run_meta_updates'],
      },
    },

    {
      -logic_name => 'run_meta_updates',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
          cmd => 'python ' . $self->o('core_metadata_script')  .  ' -o ' . $self->o('output_path') . ' -d '  . $self->o('reference_db_name') . ' -s ' . $self->o('reference_db_host') . ' -p ' .$self->o
('reference_db_port') . ' --team genebuild',
      },
      -rc_name => '1GB',
      -flow_into       => { 1 => ['load_meta_updates'], },
    },

    {
      -logic_name => 'load_meta_updates',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
          cmd => '/hps/software/users/ensembl/ensw/mysql-cmds/ensembl/ensadmin/' . $self->o('reference_db_host') . ' ' . $self->o('reference_db_name') . ' <' . $self->o('output_path') . '/' .$self->o('reference_db_name') . '.sql',
      },
      -rc_name => 'default_registry',
      -flow_into       => { 1 => ['run_core_stats'], },
    },

    {
      -logic_name => 'run_core_stats',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
          cmd => 'perl ' . $self->o('core_stats_script')  .  ' -dbname '  . $self->o('reference_db_name') . ' -host ' .  $self->o('reference_db_host') . ' -port ' .$self->o('reference_db_port') . ' -production_name ' . $self->o('production_name') . ' -output_dir ' . $self->o('output_path'),
      },
      -rc_name => '5GB',
      -flow_into       => { 1 => ['load_core_stats'], },
    },

    {
      -logic_name => 'load_core_stats',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
          cmd => '/hps/software/users/ensembl/ensw/mysql-cmds/ensembl/ensadmin/' . $self->o('reference_db_host') . ' ' . $self->o('reference_db_name') . ' <' . $self->o('output_path') . '/stats_'   .$self->o('reference_db_name') . '.sql',
      },
      -rc_name => 'default_registry',
      -flow_into       => { 1 => ['core_gene_set_sanity_checks'], },
    },
      
    {
      -logic_name => 'core_gene_set_sanity_checks',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
      -parameters => {
        skip_rnaseq                => $self->o('skip_rnaseq'),
        skip_projection            => $self->o('skip_projection'),
        target_db                  => $self->o('reference_db'),
        sanity_check_type          => 'gene_db_checks',
        min_allowed_feature_counts => get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
          'gene_db_checks' )->{ $self->o('sanity_set') }->{'core'},
      },
      -rc_name   => '4GB',
      -flow_into => {
        1 => ['core_healthchecks'],
      },
    },

    {
      -logic_name => 'core_healthchecks',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveHealthcheck',
      -parameters => {
        input_db         => $self->o('reference_db'),
        species          => $self->o('species_name'),
        group            => 'core_handover',
        enscode_root_dir => $self->o('enscode_root_dir'),
      },
      -rc_name         => 'default',
      -max_retry_count => 0,
      -flow_into       => {
        1 => ['core_assembly_name_update'], },
    },

    {
      -logic_name => 'core_assembly_name_update',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . $self->o('assembly_name_script') .
          ' -user ' . $self->o('user') .
          ' -pass ' . $self->o('password') .
          ' -host ' . $self->o( 'reference_db', '-host' ) .
          ' -port ' . $self->o( 'reference_db', '-port' ) .
          ' -dbname ' . $self->o( 'reference_db', '-dbname' ) .
          ' -driver ' . $self->o('hive_driver') .
          ' -assembly_accession ' . $self->o('assembly_accession') .
          ' -assembly_name ' . $self->o('assembly_name') .
          ' -registry_host ' . $self->o('registry_host') .
          ' -registry_port ' . $self->o('registry_port') .
          ' -registry_db ' . $self->o('registry_db'),
      },
	  -rc_name => 'default',
	  -flow_into       => {
             1 => ['pepstats'], },
    },

{
          -logic_name => 'pepstats',
          -module     => 'Bio::EnsEMBL::Production::Pipeline::Production::PepStatsBatch',
          -parameters => {
              dbtype => 'core',
              species => $self->o('production_name'),
              pepstats_binary => 'pepstats',
              reg_conf => $self->o('registry_file'),
              tmpdir => $self->o('output_path'),
          },
          -max_retry_count => 1,
          -hive_capacity   => 50,
          -rc_name => '50GB_registry',
          -flow_into       => {
                  1 => ['delete_short_reads'], },
      },
{
      -logic_name => 'delete_short_reads',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'rm -r '. $self->o('merge_dir') . '/* ' .$self->o('output_dir') . '/* ',
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['fetch_taxon_id'], },
 },


{
  -logic_name => 'fetch_taxon_id',
-module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
-parameters => {
	db_conn => $self->o('reference_db'),
	input_query => 'SELECT TRIM(meta_value) FROM meta WHERE meta_key="species.taxonomy_id";',
	prepend => ['-NB', '-q'],
	output_file => catfile($self->o('output_path'), 'taxon.txt'),
},
  -rc_name   => 'default',
  -flow_into => {
    1 => ['select_busco_lineage']
  },
},


{
  -logic_name => 'select_busco_lineage',
  -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRunExternalCmd',
  -parameters => {
    indir    => $self->o('output_path'),
    infile   => 'taxon.txt',
    cmd      => 'python ' . $self->o('busco_lineage_selector_script') .
                ' -d ' . $self->o('busco_datasets_file') .
                ' -t ' .'$'.'(cat '.catfile($self->o('output_path'), 'taxon.txt'). ')',
},
  -rc_name => 'default',
  -flow_into => {
    1 => { 'run_busco_genome' => { lineage => '#stdout#' } },
  },
},

     {
       -logic_name => 'run_busco_genome',
       -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
       -parameters => {
	 cmd => 'cd '. $self->o('output_path') . ';' .
	 'singularity exec ' . $self->o('busco_singularity_image') . ' busco -f '. '-i '. $self->o('output_path') . 'genome_dumps/' . $self->o('species_name').'_toplevel.fa '.'-m genome -l #lineage# -c 30 ' . '-o busco_genome ' . '--offline --download_path '. $self->o('busco_download_path') . ';' 
      },
           -rc_name => '120GB_30cpus',
           -flow_into => {
              1 => {'run_busco_protein' => {lineage => '#lineage#'} },
      },
   },

    {
       -logic_name => 'run_busco_protein',
       -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
       -parameters => {
         cmd => 'cd '. $self->o('output_path') . ';' .
	 'singularity exec ' . $self->o('busco_singularity_image') . ' busco -f '. '-i '. $self->o('gst_dir') . '/' . $self->o('production_name') . '_protein_sequences.fa ' . '-m protein -l #lineage# -c 30 ' . '-o busco_protein '. '--offline --download_path '. $self->o('busco_download_path') . ';' 
      },
           -rc_name => '80GB_30cpus',
           -flow_into => {
              1 => ['load_genome_busco_into_core'],
      },
   },

{
        -logic_name => 'load_genome_busco_into_core',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
          cmd => 'python ' . $self->o('busco_metakeys_patch_script') .
              ' -db '. $self->o('reference_db_name') . ' ' .  
			  ' -host ' .  $self->o('reference_db_host') . ' ' .
              ' -port ' . $self->o('reference_db_port') . ' ' .
              ' -user ' . $self->o('user') . 
              ' -password ' . $self->o('password') . 
              ' -file ' . $self->o('busco_genome_dir') . '/' . '*.txt' .
              ' -output_dir ' .  $self->o('busco_genome_dir') . ' -run_query true',
          },
         -rc_name         => 'default',
         -flow_into       => { 1 => ['load_protein_busco_into_core'], },
    }, 

{
        -logic_name => 'load_protein_busco_into_core',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
          cmd => 'python ' . $self->o('busco_metakeys_patch_script') .
              ' -db '. $self->o('reference_db_name') . ' ' .  
			  ' -host ' .  $self->o('reference_db_host') . ' ' .
              ' -port ' . $self->o('reference_db_port') . ' ' .
              ' -user ' . $self->o('user') . 
              ' -password ' . $self->o('password') . 
              ' -file ' . $self->o('busco_protein_dir') . '/' . '*.txt' .
              ' -output_dir ' .  $self->o('busco_protein_dir') . ' -run_query true',
          },
         -rc_name         => 'default',
         -flow_into       => { 1 => ['populate_registry_metrics'], },
    }, 
  {
      -logic_name => 'populate_registry_metrics',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
          cmd => 'python ' . $self->o('write_metrics_to_registry_script') .
                 ' --registry_host ' . $self->o('registry_host') .
                 ' --registry_port ' . $self->o('registry_port') .
                 ' --registry_user ' . $self->o('user') .
                 ' --registry_password ' . $self->o('password') .
                 ' --registry_db ' . $self->o('registry_db') .
                 ' --core_host ' . $self->o('reference_db', '-host') .
                 ' --core_port ' . $self->o('reference_db', '-port') .
                 ' --core_user ' . $self->o('user_r') .
                 ' --core_password ' . $self->o('password_r') .
                 ' --core_db ' . $self->o('reference_db', '-dbname') .
                 ' --genebuilder $USER' .
                 ' --assembly ' . $self->o('assembly_accession'),
      },
      -rc_name => '1GB',
      -flow_into => {
          1 => ['check_busco_score'],
      },
  },
{
        -logic_name => 'check_busco_score',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
            cmd => 'if python ' .  catfile( $self->o('enscode_root_dir'), 'ensembl-genes','src','python','ensembl','genes','metrics', 'check_busco_score.py' ) .
            ' --genome ' . $self->o('busco_genome_dir')  . '/' .  $self->o('reference_db_name').'*.json ' .
            ' --protein '. $self->o('busco_protein_dir') . '/' .  $self->o('reference_db_name').'*.json'  .
            ' --min_range_protein_score "' . $self->o('busco_lower_threshold') . '"' .
            ' --max_range_protein_score "' . $self->o('busco_threshold') . '"' .
            ' --diff_prot_gen_mode "' . $self->o('busco_difference_threshold') . '"' .
            '; then exit 0; else exit 42; fi',
          return_codes_2_branches => { '42' => 2 },
        },
        -rc_name => 'default',
        -flow_into  => {
            1 => 'update_registry_completed',
            2 => 'update_registry_as_check',
        }
    },

  {
      -logic_name => 'update_registry_completed',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
          cmd => 'python ' . $self->o('registry_status_update_python_script') .
              ' --host ' . $self->o('registry_host') .
              ' --port ' . $self->o('registry_port') .
              ' --user ' . $self->o('user') .
              ' --password ' . $self->o('password') .
              ' --database ' . $self->o('registry_db') .
              ' --assembly ' . $self->o('assembly_accession') .
              ' --status completed' .
              ' --genebuilder $USER' .
              ' --annotation_source ensembl' .
              ' --annotation_method full_genebuild',
      },
      -rc_name => 'default',
      -flow_into => {
          1 => ['backbone_job_pipeline'],
      },
  },


  {
      -logic_name => 'update_registry_as_check',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
          cmd => 'python ' . $self->o('registry_status_update_python_script') .
              ' --host ' . $self->o('registry_host') .
              ' --port ' . $self->o('registry_port') .
              ' --user ' . $self->o('user') .
              ' --password ' . $self->o('password') .
              ' --database ' . $self->o('registry_db') .
              ' --assembly ' . $self->o('assembly_accession') .
              ' --status check_busco' .
              ' --genebuilder $USER' .
              ' --annotation_source ensembl' .
              ' --annotation_method full_genebuild',
      },
      -rc_name => 'default',
  },

{
    -logic_name     => 'backbone_job_pipeline',
    -module         => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
    -hive_capacity  => -1,
    -flow_into      => {
      '1->A'  => ['gff3','gtf', 'gzip_softmasked_fasta', 'prepare_twobit'],
      'A->1'  => ['prepare_pre_release_ftp'],
    }
  },

{
    -logic_name      => 'prepare_pre_release_ftp',
    -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -parameters => {
          cmd => 'python ' . $self->o('pre_release_ftp_script').
	  ' -p '  . $self->o('output_path') .
	  ' -g '  . $self->o('assembly_accession') .
	  ' -s '  . $self->o('species_name'),
      },
    -hive_capacity => 10,
    -flow_into => {
      '1' => ['rsync_ftp_release'],
    },
},

{
    -logic_name => 'gzip_softmasked_fasta',
    -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -parameters => {
        cmd => 'gzip -c ' .
            catfile( $self->o('output_path'), 'red_output', 'mask_output', $self->o('species_name') . '_toplevel.msk')
            . ' > ' .
            catfile( $self->o('output_path'), $self->o('species_name') . '_softmasked_toplevel.fa.gz' ),
    },
    -hive_capacity => 10,
    -rc_name       => '2GB',
},

{
    -logic_name => 'prepare_twobit',
    -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -parameters => {
        cmd => 'faToTwoBit '
            . catfile(
                $self->o('output_path'),
                'red_output',
                'mask_output',
                $self->o('species_name') . '_toplevel.msk'
            )
            . ' '
            . catfile(
                $self->o('output_path'),
                $self->o('species_name') . '_softmasked_toplevel.2bit'
            ),
    },
    -hive_capacity => 10,
    -rc_name       => '2GB',
},

  ### GTF
  { -logic_name      => 'gtf',
      -module        => 'Bio::EnsEMBL::Production::Pipeline::GTF::DumpFile',
      -parameters    => {
          gtf_to_genepred => $self->o('gtftogenepred_exe'),
          gene_pred_check => $self->o('genepredcheck_exe'),
          abinitio        => $self->o('abinitio'),
          gene            => $self->o('gene'),
          species             => $self->o('production_name'),
          group               => 'core',
          base_path           => $self->o('output_path'),
          release             => $self->o('ensembl_release'),
      },
      -hive_capacity => 50,
      -rc_name       => '2GB_registry',
      -flow_into     => { '-1' => 'gtf_32GB'},
  },

  { -logic_name      => 'gtf_32GB',
      -module        => 'Bio::EnsEMBL::Production::Pipeline::GTF::DumpFile',
      -parameters    => {
          gtf_to_genepred => $self->o('gtftogenepred_exe'),
          gene_pred_check => $self->o('genepredcheck_exe'),
          abinitio        => $self->o('abinitio'),
          gene            => $self->o('gene'),
          species             => $self->o('production_name'),
          group               => 'core',
          base_path           => $self->o('output_path'),
          release             => $self->o('ensembl_release'),
      },
      -hive_capacity => 50,
      -rc_name       => '32GB_registry',
  },


  ### GFF3
  { -logic_name      => 'gff3',
      -module        => 'Bio::EnsEMBL::Production::Pipeline::GFF3::DumpFile',
      -parameters    => {
          feature_type        => $self->o('feature_type'),
          per_chromosome      => $self->o('per_chromosome'),
          include_scaffold    => $self->o('include_scaffold'),
          logic_name          => $self->o('logic_name'),
          db_type             => $self->o('db_type'),
          abinitio            => $self->o('abinitio'),
          gene                => $self->o('gene'),
          out_file_stem       => $self->o('out_file_stem'),
          xrefs               => $self->o('xrefs'),
          species             => $self->o('production_name'),
          group               => 'core',
          base_path           => $self->o('output_path'),
          release             => $self->o('ensembl_release'),
      },
      -hive_capacity  => 50,
      -rc_name        => '2GB_registry',
      -flow_into      => { '-1' => 'gff3_32GB', '1' => 'tidy_gff3', },
  },

  { -logic_name      => 'gff3_32GB',
      -module        => 'Bio::EnsEMBL::Production::Pipeline::GFF3::DumpFile',
      -parameters    => {
          feature_type        => $self->o('feature_type'),
          per_chromosome      => $self->o('per_chromosome'),
          include_scaffold    => $self->o('include_scaffold'),
          logic_name          => $self->o('logic_name'),
          db_type             => $self->o('db_type'),
          abinitio            => $self->o('abinitio'),
          gene                => $self->o('gene'),
          out_file_stem       => $self->o('out_file_stem'),
          xrefs               => $self->o('xrefs'),
          species             => $self->o('production_name'),
          group               => 'core',
          base_path           => $self->o('output_path'),
          release             => $self->o('ensembl_release'),
      },
      -hive_capacity => 50,
      -rc_name       => '32GB_registry',
      -flow_into     => { '1' => 'tidy_gff3', },
  },

  ### GFF3:post-processing
  { -logic_name      => 'tidy_gff3',
      -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters    => { cmd => $self->o('gff3_tidy') . ' -gzip -o #out_file#.sorted.gz #out_file#', },
      -hive_capacity => 10,
      -batch_size    => 10,
      -flow_into     => 'validate_gff3',
      -rc_name       => '2GB',
  },
  {
      -logic_name    => 'validate_gff3',
      -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters    => { cmd => $self->o('gff3_validate') . ' #out_file#', },
      -hive_capacity => 10,
      -batch_size    => 10,
      -rc_name       => '2GB',
  },

  {
    -logic_name => 'rsync_ftp_release',
    -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -parameters => {
        cmd => 'sudo -u genebuild rsync -ahvW ' . $self->o('ftp_release'). '/' . ' '. $self->o('production_ftp_dir')
            },
    -rc_name => 'datamover',
    -flow_into => {
        1 => ['update_registry_pre_release'],
    },
  },

  {
      -logic_name => 'update_registry_pre_release',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
          cmd => 'python ' . $self->o('registry_status_update_python_script') .
              ' --host ' . $self->o('registry_host') .
              ' --port ' . $self->o('registry_port') .
              ' --user ' . $self->o('user') .
              ' --password ' . $self->o('password') .
              ' --database ' . $self->o('registry_db') .
              ' --assembly ' . $self->o('assembly_accession') .
              ' --status pre_released' .
              ' --genebuilder $USER' .
              ' --annotation_source ensembl' .
              ' --annotation_method full_genebuild',
      },
      -rc_name => '1GB',
      -flow_into => { 1 => ['create_target_db_gb1'], },
  },

{
      -logic_name => 'create_target_db_gb1',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
          'target_db' => {
              -dbname => $self->o('production_name') . '_core_' . $self->o('release_number') . '_1',
              -host   => $ENV{GBS1},
              -port   => $ENV{GBP1},
              -user   => $self->o('user'),
              -pass   => $self->o('password'),
              -driver => $self->o('hive_driver'),
          },
          'create_type' => 'core_only',
	  'enscode_root_dir' => $self->o('enscode_root_dir'),
      },
      -rc_name => 'default',
      -flow_into => {
          1 => ['copy_core_db_to_gb1'],
      },
  },

  {
      -logic_name => 'copy_core_db_to_gb1',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::DatabaseDumper',
      -parameters => {
          'src_db_conn' => $self->o('dna_db'),
          'output_db' => {
              -dbname => $self->o('production_name') . '_core_' .  $self->o('release_number') . '_1',
              -host   => $ENV{GBS1},
              -port   => $ENV{GBP1},
              -user   => $self->o('user'),
              -pass   => $self->o('password'),
              -driver => $self->o('hive_driver'),
          },
          'exclude_ehive' => 1,
	  'dump_options' => $self->o('mysql_dump_options'),
      },
      -rc_name => '10GB',
  },

  ];
}

sub resource_classes {
  my $self = shift;

  return {
    #inherit other stuff from the base class
     %{ $self->SUPER::resource_classes() },
    }
}

1;
