=head1 LICENSE

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


package Bio::EnsEMBL::Analysis::Hive::Config::production_rnaseq_db;

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
    'dbowner'                   => '' || $ENV{EHIVE_USER} || $ENV{USER},
    'pipeline_name'             => '' || $self->o('production_name').'_'.$self->o('ensembl_release'),
    'user_r'                    => '', # read only db user
    'user'                      => '', # write db user
    'password'                  => '', # password for write db user
    'pipe_db_host'              => '', # host for pipe db
    'pipe_db_port'              => '', # port for pipeline host
    'databases_host'            => '', # host for general output dbs
    'databases_port'            => '', # port for general output db host
    'dna_db_host'               => '', # host for dna db
    'dna_db_port'               => '', # port for dna db host
    'registry_host'             => '', # host for registry db
    'registry_port'             => '', # port for registry db
    'registry_db'               => '', # name for registry db

    'release_number'            => '' || $self->o('ensembl_release'),
    'species_name'              => '', # e.g. mus_musculus
    'production_name'           => '', # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
    annotation_source           => 'ensembl', # Used to define the structure of the ftp, default value ensembl
    dbname_accession          => '', # This is the assembly accession without [._] and all lower case, i.e gca001857705v1
    'taxon_id'                  => '', # should be in the assembly report file
    'output_path'               => '', # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses
    'assembly_name'             => '', # Name (as it appears in the assembly report file)
    'assembly_accession'        => '', # Versioned GCA assembly accession, e.g. GCA_001857705.1

    'mapping_required'          => '0', # If set to 1 this will run stable_id mapping sometime in the future. At the moment it does nothing
    'mapping_db'                => '', # Tied to mapping_required being set to 1, we should have a mapping db defined in this case, leave undef for now

    'uniprot_version'           => 'uniprot_2021_04', # What UniProt data dir to use for various analyses
    'sanity_set'                => '',
    # Keys for custom loading, only set/modify if that's what you're doing
    'protein_entry_loc'            => catfile($self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), 'entry_loc'), # Used by genscan blasts and optimise daf/paf. Don't change unless you know what you're doing

    ########################
    # Pipe and ref db info
    ########################
    'rnaseq_db_host'               => $self->o('databases_host'),
    'rnaseq_db_port'               => $self->o('databases_port'),
    'rnaseq_db_name'               => $self->o('dbowner').'_'.$self->o('dbname_accession').'_rnaseq_'.$self->o('release_number'),

    'rnaseq_refine_db_host'         => $self->o('databases_host'),
    'rnaseq_refine_db_port'         => $self->o('databases_port'),
    'rnaseq_refine_db_name'         => $self->o('dbowner').'_'.$self->o('dbname_accession').'_refine_'.$self->o('release_number'),

    'rnaseq_blast_db_host'         => $self->o('databases_host'),
    'rnaseq_blast_db_port'         => $self->o('databases_port'),
    'rnaseq_blast_db_name'         => $self->o('dbowner').'_'.$self->o('dbname_accession').'_scallop_blast_'.$self->o('release_number'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'ensembl_release'              => $ENV{ENSEMBL_RELEASE}, # this is the current release version on staging to be able to get the correct database
    'production_db_host'           => 'mysql-ens-meta-prod-1',
    'production_db_port'           => '4483',

    databases_to_delete => ['rnaseq_db'],

    ########################
    # BLAST db paths
    ########################
    'base_blast_db_path'        => $ENV{BLASTDB_DIR},

    ######################################################
    #
    # Mostly constant settings
    #
    ######################################################
    ensembl_analysis_script           => catdir($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts'),
    load_optimise_script              => catfile($self->o('ensembl_analysis_script'), 'genebuild', 'load_external_db_ids_and_optimize_af.pl'),
    delete_genes_script     => catfile($self->o('ensembl_analysis_script'), 'genebuild', 'delete_genes.pl'),
    add_analyses_script => catfile( $self->o('ensembl_analysis_script'), 'genebuild', 'add_rnaseq_analysis_descriptions.pl' ),

    ensembl_misc_script        => catdir($self->o('enscode_root_dir'), 'ensembl', 'misc-scripts'),
    meta_coord_script          => catfile($self->o('ensembl_misc_script'), 'meta_coord', 'update_meta_coord.pl'),
    meta_levels_script         => catfile($self->o('ensembl_misc_script'), 'meta_levels.pl'),
    assembly_name_script       => catfile($self->o('ensembl_analysis_script'), 'update_assembly_name.pl'),

    rnaseq_daf_introns_file => catfile($self->o('output_dir'), 'rnaseq_daf_introns.dat'),

    gene_ids_per_file => 600,

    ########################
    # Executable paths
    ########################
    deeptools_bamcoverage_path => catfile($self->o('software_base_path'), 'pyenv', 'versions', 'genebuild', 'bin', 'bamCoverage'),

    use_threads => 4,

    # RNA-seq pipeline stuff
    'rnaseq_dir'    => catdir($self->o('output_path'), 'rnaseq'),
    output_dir    => catdir($self->o('rnaseq_dir'),'output'),
    merge_dir     => catdir($self->o('rnaseq_dir'),'merge'),

    delete_genes_dir => catdir($self->o('rnaseq_dir'),'delete_merge_genes'),
    delete_genes_prefix => catfile($self->o('delete_genes_dir'), 'genes_to_delete.'),
    optimise_dir => catdir($self->o('rnaseq_dir'),'optimise_rnaseq'),
    production_ftp_dir => '/nfs/production/flicek/ensembl/production/ensemblftp/rapid-release/',
    species_list => '/nfs/production/flicek/ensembl/genebuild/do_not_delete/main_species.csv',

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # No option below this mark should be modified
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ########################
    # db info
    ########################
    'rnaseq_db' => {
      -dbname => $self->o('rnaseq_db_name'),
      -host   => $self->o('rnaseq_db_host'),
      -port   => $self->o('rnaseq_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'rnaseq_blast_db' => {
      -dbname => $self->o('rnaseq_blast_db_name'),
      -host   => $self->o('rnaseq_blast_db_host'),
      -port   => $self->o('rnaseq_blast_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'rnaseq_refine_db' => {
      -dbname => $self->o('rnaseq_refine_db_name'),
      -host   => $self->o('rnaseq_refine_db_host'),
      -port   => $self->o('rnaseq_refine_db_port'),
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

  };
}


sub pipeline_create_commands {
  my ($self) = @_;

  return [
    # inheriting database and hive tables' creation
    @{$self->SUPER::pipeline_create_commands},

    'mkdir -p '.$self->o('delete_genes_dir'),
    'mkdir -p '.$self->o('optimise_dir'),
  ];
}


sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
    production_ftp_dir => $self->o('production_ftp_dir'),
  }
}

## See diagram for pipeline structure
sub pipeline_analyses {
  my ($self) = @_;

  return [
    {
      -logic_name => 'create_rnaseq_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db => $self->o('rnaseq_blast_db'),
        target_db => $self->o('rnaseq_db'),
        create_type => 'copy',
      },
      -max_retry_count => 0,
      -rc_name    => 'default',
      -input_ids  => [{}],
      -flow_into => {
        '1' => ['prepare_rnaseq_meta_data'],
      },
    },

    {
      -logic_name => 'prepare_rnaseq_meta_data',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('rnaseq_db'),
        sql => [
          'TRUNCATE dna_align_feature',
          'DELETE FROM transcript_supporting_feature WHERE feature_type = "dna_align_feature"',
          'DELETE FROM supporting_feature WHERE feature_type = "dna_align_feature"',
          'DELETE FROM analysis WHERE logic_name NOT LIKE "%rnaseq%"',
          'DELETE FROM analysis WHERE logic_name LIKE "%merged%"',
          'INSERT INTO analysis (logic_name, module, created, db_version) VALUES ("other_protein", "HiveBlastRNAseq", NOW(), "#uniprot_version#")',
          'INSERT INTO meta (species_id, meta_key, meta_value) VALUES (1, "species.annotation_source", "ensembl")',
          'UPDATE protein_align_feature paf, analysis a SET paf.analysis_id = a.analysis_id WHERE a.logic_name = "other_protein"',
          'DELETE FROM meta WHERE meta_key LIKE "assembly.web_accession%"',
          'DELETE FROM meta WHERE meta_key LIKE "removed_evidence_flag\.%"',
          'DELETE FROM meta WHERE meta_key LIKE "marker\.%"',
          'DELETE FROM meta WHERE meta_key IN ("genebuild.method","genebuild.projection_source_db","genebuild.start_date","repeat.analysis", "genebuild.method_display")',
          'UPDATE gene JOIN transcript USING(gene_id) SET canonical_transcript_id = transcript_id',
          'UPDATE transcript JOIN translation USING(transcript_id) SET canonical_translation_id = translation_id',
          'UPDATE intron_supporting_evidence ise, analysis a1, analysis a2 SET ise.analysis_id = a2.analysis_id'.
            ' WHERE ise.analysis_id = a1.analysis_id AND a2.logic_name = REPLACE(a1.logic_name, "rnaseq_gene", "rnaseq_ise")',
          'UPDATE gene SET source = "ensembl", biotype = "protein_coding", stable_id = NULL',
          'UPDATE transcript SET source = "ensembl", biotype = "protein_coding", stable_id = NULL',
          'UPDATE translation SET stable_id = NULL',
          'UPDATE exon SET stable_id = NULL',
        ],
        uniprot_version => $self->o('uniprot_version'),
      },
      -rc_name    => 'default',

      -flow_into => {
        '1' => ['generate_rnaseq_stable_ids'],
      },
    },

    {
      -logic_name => 'generate_rnaseq_stable_ids',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::SetStableIDs',
      -parameters => {
        enscode_root_dir => $self->o('enscode_root_dir'),
        mapping_required => 0,
        target_db => $self->o('rnaseq_db'),
        id_start => 'RNASEQ01',
        output_path => $self->o('output_path'),
        _stable_id_file => 'rnaseq_stable_ids.sql',
      },
      -rc_name    => '2GB',
      -flow_into => {
        1 => ['dump_daf_introns'],
      },
    },

    {
      -logic_name => 'dump_daf_introns',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
      -parameters => {
        db_conn => $self->o('rnaseq_refine_db'),
        input_query => 'SELECT daf.* FROM dna_align_feature daf, analysis a'.
          ' WHERE daf.analysis_id = a.analysis_id AND a.logic_name != "rough_transcripts" AND a.logic_name NOT LIKE "%\_merged_rnaseq_daf"',
        command_out => q(sort -nk2 -nk3 -nk4 | sed 's/NULL/\\N/g;s/^[0-9]\+/\\N/' > #daf_file#),
        daf_file => $self->o('rnaseq_daf_introns_file'),
        prepend => ['-NB', '-q'],
      },
      -rc_name => 'default',
      -flow_into => {
        1 => ['load_daf_introns'],
      },
    },

    {
      -logic_name => 'load_daf_introns',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
      -parameters => {
        db_conn => $self->o('rnaseq_db'),
        input_query => 'LOAD DATA LOCAL INFILE "#daf_file#" INTO TABLE dna_align_feature',
        daf_file => $self->o('rnaseq_daf_introns_file'),
      },
      -rc_name => 'default',
      -flow_into => {
        1 => ['update_dna_align_feature_table'],
      },
    },

    {
      -logic_name => 'update_dna_align_feature_table',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('rnaseq_db'),
        sql => [
          'UPDATE dna_align_feature SET external_db_id = NULL',
          'UPDATE dna_align_feature SET align_type = "ensembl"',
        ],
      },
      -rc_name    => 'default',
      -flow_into => {
        1 => ['set_rnaseq_meta_coords'],
      },
    },

    {
      -logic_name => 'set_rnaseq_meta_coords',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('meta_coord_script').
          ' -user '.$self->o('user').
          ' -pass '.$self->o('password').
          ' -host '.$self->o('rnaseq_db','-host').
          ' -port '.$self->o('rnaseq_db','-port').
          ' -dbpattern '.$self->o('rnaseq_db','-dbname')
      },
      -rc_name => 'default',
      -flow_into => {
        1 => ['set_rnaseq_meta_levels'],
      },
    },

    {
      -logic_name => 'set_rnaseq_meta_levels',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('meta_levels_script').
          ' -user '.$self->o('user').
          ' -pass '.$self->o('password').
          ' -host '.$self->o('rnaseq_db','-host').
          ' -port '.$self->o('rnaseq_db','-port').
          ' -dbname '.$self->o('rnaseq_db','-dbname')
      },
      -rc_name => 'default',
      -flow_into => { 1 => ['optimise_rnaseq'] },
    },

    {
      -logic_name => 'optimise_rnaseq',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('load_optimise_script').
          ' -output_path '.$self->o('optimise_dir').
          ' -uniprot_filename '.$self->o('protein_entry_loc').
          ' -dbuser '.$self->o('user').
          ' -dbpass '.$self->o('password').
          ' -dbport '.$self->o('rnaseq_db','-port').
          ' -dbhost '.$self->o('rnaseq_db','-host').
          ' -dbname '.$self->o('rnaseq_db','-dbname').
          ' -prod_dbuser '.$self->o('user_r').
          ' -prod_dbhost '.$self->o('production_db','-host').
          ' -prod_dbname '.$self->o('production_db','-dbname').
          ' -prod_dbport '.$self->o('production_db','-port').
          ' -nodaf -ise'
      },
      -max_retry_count => 0,
      -rc_name => '8GB',
      -flow_into => {
        1 => ['rnaseq_gene_sanity_checks'],
      },
    },

    {
      -logic_name => 'rnaseq_gene_sanity_checks',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
      -parameters => {
        target_db => $self->o('rnaseq_db'),
        sanity_check_type => 'gene_db_checks',
        min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
          'gene_db_checks')->{ $self->o('sanity_set') }->{'rnaseq_final'},
      },

      -rc_name    => '8GB',
      -flow_into => {
        1 => ['create_bam_file_job'],
      },
    },

    {
      -logic_name => 'create_bam_file_job',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        inputcmd => 'cd #working_dir#; ls *.bam',
        column_names => ['bam_file'],
        working_dir => $self->o('merge_dir'),
      },
      -rc_name => 'default',
      -flow_into  => {
        '2->A' => ['bam2bigwig'],
        'A->1' => ['concat_md5_sum'],
      },
    },

    {
      -logic_name => 'bam2bigwig',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        deeptools_bamcoverage => $self->o('deeptools_bamcoverage_path'),
        num_cpus => $self->o('use_threads'),
        cmd => 'TMPDIR=#working_dir# ; #deeptools_bamcoverage# --numberOfProcessors #num_cpus# --binSize 1 -b '.catfile('#working_dir#','#bam_file#').' -o '.catfile('#working_dir#', "#expr(join('', split('.bam',#bam_file#)))expr#").'.bw',
        working_dir => $self->o('merge_dir'),
      },
      -rc_name => '10GB_multithread',
      -flow_into  => {
        1 => ['md5_sum'],
      },
    },

    {
      -logic_name => 'md5_sum',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'cd #working_dir#;' .
        'md5sum ' . "#expr(join('', split('.bam',#bam_file#)))expr#" . '.bw #bam_file# #bam_file#.csi > #bam_file#.md5',
        working_dir => $self->o('merge_dir'),
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'concat_md5_sum',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'cat '.catfile('#working_dir#', '*md5').' | grep ".bw" > '.catfile('#working_dir#', 'md5sum.txt.1'),
        working_dir => $self->o('merge_dir'),
      },
      -rc_name => 'default',
      -flow_into  => {
        1 => ['clean_concat_md5_sum'],
      },
    },

    {
      -logic_name => 'clean_concat_md5_sum',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'rm '.catfile('#working_dir#', '*md5'),
        working_dir => $self->o('merge_dir'),
      },
      -rc_name => 'default',
      -flow_into  => {
        1 => ['create_readme'],
      },
    },

    {
      -logic_name => 'create_readme',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'cd #working_dir#;
        FILES=($(ls *.bw));
        printf "#free_text#" | sed "s/NUM/$((${#FILES[*]}))/g;s/ \([a-z]\)\([a-z]\+_\)/ \U\1\E\2/;s/_/ /g" > README.1
        IFS=$\'\n\';
        echo "${FILES[*]}" >> README.1',
        working_dir => $self->o('merge_dir'),
        species_name  => $self->o('species_name'),
        free_text => 'Note\n------\n\n'.
          'The RNASeq data for #species_name# consists of NUM individual sample(s).\n\n'.
          'The BigWig file (.bw) contains the coverage information.\n\n'.
          'Use the md5sum.txt.1 file to check the integrity of the downloaded files.\n\n'.
          'Files\n-----\n',
      },
      -rc_name => 'default',
      -flow_into  => {
        1 => ['rnaseq_healthchecks'],
      },
    },

    {
      -logic_name => 'rnaseq_healthchecks',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveHealthcheck',
      -parameters => {
        input_db => $self->o('rnaseq_db'),
        species  => $self->o('species_name'),
        group    => 'rnaseq_handover',
      },
      -max_retry_count => 0,

      -rc_name    => '4GB',
      -flow_into  => {
        1 => ['add_analyses_descriptions'],
      },
    },

    {
      -logic_name => 'add_analyses_descriptions',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('add_analyses_script').' -host #host# -port #port# -dbname #dbname#',
        dbname => $self->o('rnaseq_db_name'),
        host   => $self->o('rnaseq_db_host'),
        port   => $self->o('rnaseq_db_port'),
      },
      -rc_name => 'default',
      -flow_into  => {
        1 => ['rnaseq_assembly_name_update'],
      },
    },

    {
      -logic_name => 'rnaseq_assembly_name_update',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('assembly_name_script').
          ' -user '.$self->o('user').
          ' -pass '.$self->o('password').
          ' -host '.$self->o('rnaseq_db','-host').
          ' -port '.$self->o('rnaseq_db','-port').
          ' -dbname '.$self->o('rnaseq_db','-dbname').
          ' -driver '.$self->o('hive_driver').
          ' -assembly_accession '.$self->o('assembly_accession').
          ' -assembly_name '.$self->o('assembly_name').
          ' -registry_host '.$self->o('registry_host').
          ' -registry_port '.$self->o('registry_port').
          ' -registry_db '.$self->o('registry_db').
          ' -working_dir '.$self->o('merge_dir'),
      },

      -rc_name => 'default',
      -flow_into  => {
        1 => ['create_species_ftp_dir'],
      },
    },

    {
      -logic_name => 'create_species_ftp_dir',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
	      cmd => 'sudo -u genebuild mkdir -p ' . catdir('#production_ftp_dir#', 'species', ucfirst($self->o('species_name')), $self->o('assembly_accession'), $self->o('annotation_source'), 'rnaseq'),
                     },
      -rc_name    => '2GB',
      -flow_into => {
	1 => ['copy_bigwig_files'],
      }
    },

   {
    -logic_name => 'copy_bigwig_files',
    -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -parameters => {
      cmd => 'sudo -u genebuild rsync -ahvW '.$self->o('merge_dir').'/*.bw '.catdir('#production_ftp_dir#', 'species', ucfirst($self->o('species_name')), $self->o('assembly_accession'), $self->o('annotation_source'), 'rnaseq').' && rsync -avhc ' . $self->o('merge_dir') . '/*.bw '.catdir('#production_ftp_dir#', 'species', ucfirst($self->o('species_name')), $self->o('assembly_accession'), $self->o('annotation_source'), 'rnaseq'),
     },
     -rc_name    => '2GB',
     -flow_into => {
	1 => ['copy_readme_md5sum'],
      },
   },

  {
    -logic_name => 'copy_readme_md5sum',
    -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -parameters => {
      cmd => 'sudo -u genebuild rsync -ahvW '.$self->o('merge_dir').'/*.1 '.catdir('#production_ftp_dir#', 'species', ucfirst($self->o('species_name')), $self->o('assembly_accession'), $self->o('annotation_source'), 'rnaseq').' && rsync -avhc ' . $self->o('merge_dir') . '/*.1 '.catdir('#production_ftp_dir#', 'species', ucfirst($self->o('species_name')), $self->o('assembly_accession'), $self->o('annotation_source'), 'rnaseq'),
     },
    -rc_name    => '2GB',
    -flow_into => {
      1 => ['set_dir_permission'],
    },
  },

   {
     -logic_name => 'set_dir_permission',
     -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
     -parameters => {
       cmd => "sudo -u genebuild find " .catdir('#production_ftp_dir#', 'species', ucfirst($self->o('species_name'))). " -user genebuild -exec chmod g+w {} \\;"},
     -rc_name    => '2GB',
     -flow_into => {
        1 => ['check_main'],
     },
   },

   {
      -logic_name => 'check_main',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
         cmd => 'if [[ $(#search_query#) ]]; then exit 0; else exit 42;fi',
         return_codes_2_branches => {'42' => 2},
	 #note that this query will search for the first part of the binomial species name, i.e. the genus, but is safest for now
	 search_query => 'search_name=`echo "${'.$self->o('species_name').'^}" | cut -d\'_\' -f1`; grep -w $search_name '.$self->o('species_list'),
      },
    -rc_name => '2GB',
    -flow_into => {
         1 => ['copy_to_main_ftp'],
         2 => ['delete_data_files'],
     }
   },

   {
     -logic_name => 'delete_data_files',
     -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
     -parameters => {
        cmd => 'rm -r '. $self->o('merge_dir') . '/* ' .$self->o('output_dir') . '/* ',
     },
   },

   { #for now we won't do anything with the file, need to get SOP for moving to the main FTP space
      -logic_name => 'copy_to_main_ftp',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name => 'default',
   }

  ];
}


sub resource_classes {
  my $self = shift;

  return {
     #inherit other stuff from the base class
     %{ $self->SUPER::resource_classes() },
    '10GB_multithread' => { LSF => $self->lsf_resource_builder('production', 10000, undef, undef, $self->default_options->{'use_threads'}),
                            SLURM =>  $self->slurm_resource_builder(10000, '7-00:00:00',$self->default_options->{'use_threads'}),},
  }
}


1;
