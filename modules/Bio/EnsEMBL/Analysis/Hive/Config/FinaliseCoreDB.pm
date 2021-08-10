
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
    'pipe_db_server'     => '',                                                                                                # host for pipe db
    'databases_server'   => '',                                                                                                # host for general output dbs
    'dna_db_server'      => '',                                                                                                # host for dna db
    'pipe_db_port'       => '',                                                                                                # port for pipeline host
    'databases_port'     => '',                                                                                                # port for general output db host
    'dna_db_port'        => '',                                                                                                # port for dna db host
    'registry_host'      => '',                                                                                                # host for registry db
    'registry_port'      => '',                                                                                                # port for registry db
    'registry_db'        => '',                                                                                                # name for registry db
    'species_name'       => '',                                                                                                # e.g. mus_musculus
    'production_name'    => '',                                                                                                # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
    'release_number'     => '' || $self->o('ensembl_release'),
    'uniprot_set'        => '',                                                                                                # e.g. mammals_basic, check UniProtCladeDownloadStatic.pm module in hive config dir for suitable set,
    'output_path'        => '',                                                                                                # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses
    'assembly_name'      => '',                                                                                                # Name (as it appears in the assembly report file)
    'assembly_accession' => '',                                                                                                # Versioned GCA assembly accession, e.g. GCA_001857705.1
    'stable_id_prefix'   => '',                                                                                                # e.g. ENSPTR. When running a new annotation look up prefix in the assembly registry db
    'stable_id_start'    => '0',                                                                                               # When mapping is not required this is usually set to 0
    'skip_projection'    => '0',                                                                                               # Will skip projection process if 1
    'skip_rnaseq'        => '0',                                                                                               # Will skip rnaseq analyses if 1
    'mapping_required'   => '0',                                                                                               # If set to 1 this will run stable_id mapping sometime in the future. At the moment it does nothing
    'mapping_db'         => '',                                                                                                # Tied to mapping_required being set to 1, we should have a mapping db defined in this case, leave undef for now
    'uniprot_version'    => 'uniprot_2019_04',                                                                                 # What UniProt data dir to use for various analyses
    'protein_entry_loc'  => catfile( $self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), 'entry_loc' ),    # Used by genscan blasts and optimise daf/paf. Don't change unless you know what you're doing

########################
# Pipe and ref db info
########################

    'projection_source_db_name'         => '',                                                                                 # This is generally a pre-existing db, like the current human/mouse core for example

    'pipe_db_name' => $self->o('dbowner') . '_' . $self->o('production_name') . '_pipe_' . $self->o('release_number'),
    'dna_db_name'  => $self->o('dbowner') . '_' . $self->o('production_name') . '_core_' . $self->o('release_number'),

    'reference_db_name'   => $self->o('dna_db_name'),
    'reference_db_server' => $self->o('dna_db_server'),
    'reference_db_port'   => $self->o('dna_db_port'),

    'final_geneset_db_server' => $self->o('databases_server'),
    'final_geneset_db_port'   => $self->o('databases_port'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'ensembl_release'      => $ENV{ENSEMBL_RELEASE},     # this is the current release version on staging to be able to get the correct database
    'production_db_server' => 'mysql-ens-meta-prod-1',
    'production_db_port'   => '4483',

########################
# BLAST db paths
########################
    'base_blast_db_path' => $ENV{BLASTDB_DIR},

    ensembl_analysis_script => catdir( $self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts' ),
    load_optimise_script => catfile( $self->o('ensembl_analysis_script'), 'genebuild', 'load_external_db_ids_and_optimize_af.pl' ),
    ensembl_misc_script => catdir( $self->o('enscode_root_dir'), 'ensembl', 'misc-scripts' ),
    meta_coord_script => catfile( $self->o('ensembl_misc_script'), 'meta_coord', 'update_meta_coord.pl' ),
    meta_levels_script       => catfile( $self->o('ensembl_misc_script'),     'meta_levels.pl' ),
    frameshift_attrib_script => catfile( $self->o('ensembl_misc_script'),     'frameshift_transcript_attribs.pl' ),
    select_canonical_script  => catfile( $self->o('ensembl_misc_script'),     'canonical_transcripts', 'select_canonical_transcripts.pl' ),
    assembly_name_script     => catfile( $self->o('ensembl_analysis_script'), 'update_assembly_name.pl' ),

    # Genes biotypes to ignore from the final db when copying to core
    copy_biotypes_to_ignore => {
      'low_coverage' => 1,
      'CRISPR'       => 1,
    },

########################
# Extra db settings
########################
    mysql_dump_options => '--max_allowed_packet=1000MB',

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

    'final_geneset_db' => {
      -dbname => $self->o('dbowner') . '_' . $self->o('production_name') . '_final_' . $self->o('release_number'),
      -host   => $self->o('final_geneset_db_server'),
      -port   => $self->o('final_geneset_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
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

  };
}

## See diagram for pipeline structure
sub pipeline_analyses {
  my ($self) = @_;

  return [

    {
      -logic_name => 'create_gene_ids_to_copy',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db    => $self->o('final_geneset_db'),
        iid_type     => 'feature_id',
        feature_type => 'gene',
        batch_size   => 100,
      },
      -input_ids  => [{}],
      -flow_into => {
        '2->A' => ['copy_genes_to_core'],
        'A->1' => ['update_biotypes_and_analyses'],
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'copy_genes_to_core',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCopyGenes',
      -parameters => {
        copy_genes_directly     => 1,
        source_db               => $self->o('final_geneset_db'),
        dna_db                  => $self->o('dna_db'),
        target_db               => $self->o('reference_db'),
        copy_biotypes_to_ignore => $self->o('copy_biotypes_to_ignore'),
      },
      -rc_name => '1GB',
    },

    {
      -logic_name => 'update_biotypes_and_analyses',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('reference_db'),
        sql     => [
          'UPDATE gene SET biotype = "protein_coding" WHERE biotype = "ensembl"',
          'UPDATE gene SET biotype = "vault_RNA" WHERE biotype in ("Vault_RNA","vaultRNA")',
          'UPDATE analysis set logic_name="cdna2genome" where logic_name="best_targetted"',
          'UPDATE gene SET analysis_id = (SELECT analysis_id FROM analysis WHERE logic_name = "ensembl")' .
            ' WHERE analysis_id IN' .
            ' (SELECT analysis_id FROM analysis' .
            ' WHERE logic_name IN ("spliced_elsewhere","pseudogenes","genblast","genblast_not_best","project_transcripts","project_pseudogene","project_lincrna","cesar"))',
          'UPDATE transcript JOIN gene USING(gene_id) SET transcript.biotype = gene.biotype',
          'UPDATE transcript JOIN gene USING(gene_id) SET transcript.analysis_id = gene.analysis_id',
          'INSERT IGNORE into analysis (created,logic_name,db) VALUES (NOW(),"other_protein","uniprot")',
          'INSERT IGNORE into analysis (created,logic_name,db) VALUES (NOW(),"projected_transcript","' . $self->o('projection_source_db_name') . '")',
          'UPDATE protein_align_feature SET analysis_id =' .
            '(SELECT analysis_id FROM analysis WHERE logic_name = "projected_transcript") WHERE analysis_id IN ' .
            '(SELECT analysis_id FROM analysis WHERE logic_name IN ("project_transcripts","cesar"))',
          'UPDATE protein_align_feature SET analysis_id =' .
            '(SELECT analysis_id FROM analysis WHERE logic_name = "other_protein") WHERE analysis_id NOT IN' .
            '(SELECT analysis_id FROM analysis WHERE logic_name IN ("uniprot","projected_transcript"))',
          'UPDATE dna_align_feature SET analysis_id =' .
            '(SELECT analysis_id FROM analysis WHERE logic_name = "projected_transcript") WHERE analysis_id IN' .
            '(SELECT analysis_id FROM analysis WHERE logic_name IN ("project_lincrna","project_pseudogene"))',
          'UPDATE repeat_feature SET repeat_start = 1 WHERE repeat_start < 1',
          'UPDATE repeat_feature SET repeat_end = 1 WHERE repeat_end < 1',
          'UPDATE gene SET analysis_id=(select analysis_id from analysis where logic_name="ensembl") WHERE analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name="filter_lncrnas")',
          'UPDATE transcript SET analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name="ensembl") WHERE analysis_id=(SELECT analysis_id from analysis WHERE logic_name="filter_lncrnas")',
          'DELETE FROM analysis WHERE logic_name="filter_lncrnas"',
          'UPDATE gene SET display_xref_id=NULL',
          'UPDATE transcript SET display_xref_id=NULL',
          'INSERT IGNORE into analysis (logic_name, db_version, db_file, program_file, module) values ("rfamcmsearch", "14.0", "' . $self->o('output_path') . '", "/nfs/software/ensembl/RHEL7-JUL2017-core2/linuxbrew/bin/cmsearch", "HiveCMSearch")',
          'INSERT IGNORE into analysis (created,logic_name) values (now(),"cdna_alignment_core")',
          'UPDATE dna_align_feature set analysis_id=(select analysis_id from analysis where logic_name="rfamcmsearch") where analysis_id=(select analysis_id from analysis where logic_name="ncrna")',
          'UPDATE dna_align_feature set analysis_id=(select analysis_id from analysis where logic_name="cdna_alignment_core") where analysis_id in (select analysis_id from analysis where logic_name in ("exonerate","cdna2genome","best_targetted"))',
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
      -flow_into => { 1 => ['null_columns'] },
    },

    {
      -logic_name => 'null_columns',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('reference_db'),
        sql     => [
          'UPDATE transcript set stable_id = NULL',
          'UPDATE translation set stable_id = NULL',
          'UPDATE exon set stable_id = NULL',
          'UPDATE protein_align_feature set external_db_id = NULL',
          'UPDATE dna_align_feature set external_db_id = NULL',
        ],
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['run_stable_ids'],
      },
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
        1 => ['update_ISE'],
      },
    },

    {
      -logic_name => 'update_ISE',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('reference_db'),
        sql     => [
          'INSERT IGNORE into analysis (created,logic_name) values (now(),"rnaseq_intron_support")',
          'UPDATE intron_supporting_evidence SET analysis_id=(select analysis_id from analysis where logic_name="rnaseq_intron_support")',
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
          'DELETE FROM analysis WHERE logic_name IN' .
            ' ("spliced_elsewhere","pseudogenes","genblast","genblast_not_best","project_pseudogene","blast_long_read",' .
            ' "project_lincrna","project_transcripts","ig_tr_collapse", "exonerate", "cdna2genome", "best_targetted")',
          'DELETE FROM ad USING analysis_description ad LEFT JOIN analysis a ON ad.analysis_id = a.analysis_id' .
            ' WHERE a.analysis_id IS NULL',
          'DELETE from analysis where logic_name like "%_rnaseq_ise%"',
        ],
      },
      -rc_name   => 'default',
      -flow_into => {
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
            '(1, "genebuild.last_geneset_update", (SELECT CONCAT((EXTRACT(YEAR FROM now())),"-",(LPAD(EXTRACT(MONTH FROM now()),2,"0")))))'
        ],
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['final_cleaning'],
      },
    },

    {
      -logic_name => 'final_cleaning',
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
          'DELETE exon FROM exon LEFT JOIN exon_transcript ON exon.exon_id = exon_transcript.exon_id WHERE exon_transcript.exon_id IS NULL',
          'DELETE supporting_feature FROM supporting_feature LEFT JOIN exon ON supporting_feature.exon_id = exon.exon_id WHERE exon.exon_id IS NULL',
          'DELETE supporting_feature FROM supporting_feature LEFT JOIN dna_align_feature ON feature_id = dna_align_feature_id WHERE feature_type="dna_align_feature" AND dna_align_feature_id IS NULL',
          'DELETE supporting_feature FROM supporting_feature LEFT JOIN protein_align_feature ON feature_id = protein_align_feature_id WHERE feature_type="protein_align_feature" AND protein_align_feature_id IS NULL',
          'DELETE transcript_supporting_feature FROM transcript_supporting_feature LEFT JOIN dna_align_feature ON feature_id = dna_align_feature_id WHERE feature_type="dna_align_feature" AND dna_align_feature_id IS NULL',
          'DELETE transcript_supporting_feature FROM transcript_supporting_feature LEFT JOIN protein_align_feature ON feature_id = protein_align_feature_id WHERE feature_type="protein_align_feature" AND protein_align_feature_id IS NULL',
        ],
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['add_placeholder_sample_location'],
      },
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
        1 => ['core_gene_set_sanity_checks'],
      },
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
          'gene_db_checks' )->{ $self->o('uniprot_set') }->{'core'},
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
    },

  ];
}

sub resource_classes {
  my $self = shift;

  return {
    'default' => { LSF => $self->lsf_resource_builder( 'production-rh74', 900, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    'default_registry' => { LSF => [ $self->lsf_resource_builder( 'production-rh74', 900, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ), '-reg_conf ' . $self->default_options->{registry_file} ] },
    '1GB'  => { LSF => $self->lsf_resource_builder( 'production-rh74', 1000,  [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    '4GB'  => { LSF => $self->lsf_resource_builder( 'production-rh74', 4000,  [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    '8GB'  => { LSF => $self->lsf_resource_builder( 'production-rh74', 8000,  [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    '10GB' => { LSF => $self->lsf_resource_builder( 'production-rh74', 10000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    }
}

1;
