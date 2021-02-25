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


package Bio::EnsEMBL::Analysis::Hive::Config::production_rnaseq_db_conf;

use strict;
use warnings;
use File::Spec::Functions;

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
    'pipe_db_server'            => '', # host for pipe db
    'pipe_db_port'              => '', # port for pipeline host
    'databases_server'          => '', # host for general output dbs
    'databases_port'            => '', # port for general output db host
    'dna_db_server'             => '', # host for dna db
    'dna_db_port'               => '', # port for dna db host
    'registry_host'             => '', # host for registry db
    'registry_port'             => '', # port for registry db
    'registry_db'               => '', # name for registry db

    'release_number'            => '' || $self->o('ensembl_release'),
    'species_name'              => '', # e.g. mus_musculus
    'production_name'           => '', # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
    'taxon_id'                  => '', # should be in the assembly report file
    'output_path'               => '', # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses
    'assembly_name'             => '', # Name (as it appears in the assembly report file)
    'assembly_accession'        => '', # Versioned GCA assembly accession, e.g. GCA_001857705.1

    'mapping_required'          => '0', # If set to 1 this will run stable_id mapping sometime in the future. At the moment it does nothing
    'mapping_db'                => '', # Tied to mapping_required being set to 1, we should have a mapping db defined in this case, leave undef for now

    'uniprot_version'           => 'uniprot_2019_04', # What UniProt data dir to use for various analyses

    # Keys for custom loading, only set/modify if that's what you're doing
    'protein_entry_loc'            => catfile($self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), 'entry_loc'), # Used by genscan blasts and optimise daf/paf. Don't change unless you know what you're doing

    ########################
    # Pipe and ref db info
    ########################
    'rnaseq_db_server'             => $self->o('databases_server'),
    'rnaseq_db_port'               => $self->o('databases_port'),

    'rnaseq_refine_db_server'       => $self->o('databases_server'),
    'rnaseq_refine_db_port'         => $self->o('databases_port'),

    'rnaseq_blast_db_server'       => $self->o('databases_server'),
    'rnaseq_blast_db_port'         => $self->o('databases_port'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'ensembl_release'              => $ENV{ENSEMBL_RELEASE}, # this is the current release version on staging to be able to get the correct database
    'production_db_server'         => 'mysql-ens-meta-prod-1',
    'production_db_port'           => '4483',

    databases_to_delete => ['rnaseq_db'],

    ######################################################
    #
    # Mostly constant settings
    #
    ######################################################
    ensembl_analysis_script           => catdir($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts'),
    load_optimise_script              => catfile($self->o('ensembl_analysis_script'), 'genebuild', 'load_external_db_ids_and_optimize_af.pl'),

    ensembl_misc_script        => catdir($self->o('enscode_root_dir'), 'ensembl', 'misc-scripts'),
    meta_coord_script          => catfile($self->o('ensembl_misc_script'), 'meta_coord', 'update_meta_coord.pl'),
    meta_levels_script         => catfile($self->o('ensembl_misc_script'), 'meta_levels.pl'),
    assembly_name_script       => catfile($self->o('ensembl_analysis_script'), 'update_assembly_name.pl'),

    rnaseq_daf_introns_file => catfile($self->o('output_dir'), 'rnaseq_daf_introns.dat'),

    ########################
    # Extra db settings
    ########################
    'num_tokens' => 10,

    ########################
    # Executable paths
    ########################
    deeptools_bamcoverage_path => '/nfs/software/ensembl/RHEL7-JUL2017-core2/pyenv/versions/genebuild/bin/bamCoverage',

    'output_dir'    => catdir($self->o('rnaseq_dir'),'output'),
    'merge_dir'     => catdir($self->o('rnaseq_dir'),'merge'),

    use_threads => 3,
    rnaseq_merge_threads => 12,
    rnaseq_merge_type => 'samtools',
    read_min_paired => 50,
    read_min_mapped => 50,
    other_isoforms => 'other', # If you don't want isoforms, set this to undef
    maxintron => 200000,

    # Please assign some or all columns from the summary file to the
    # some or all of the following categories.  Multiple values can be
    # separted with commas. ID, SM, DS, CN, is_paired, filename, read_length, is_13plus,
    # is_mate_1 are required. If pairing_regex can work for you, set is_mate_1 to -1.
    # You can use any other tag specified in the SAM specification:
    # http://samtools.github.io/hts-specs/SAMv1.pdf

    ####################################################################
    # This is just an example based on the file snippet shown below.  It
    # will vary depending on how your data looks.
    ####################################################################
    file_columns      => ['SM', 'ID', 'is_paired', 'filename', 'is_mate_1', 'read_length', 'is_13plus', 'CN', 'PL', 'DS'],
    long_read_columns => ['sample','filename'],

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # No option below this mark should be modified
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ########################
    # db info
    ########################
    'rnaseq_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_rnaseq_'.$self->o('release_number'),
      -host   => $self->o('rnaseq_db_server'),
      -port   => $self->o('rnaseq_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'rnaseq_blast_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_rnaseq_blast_'.$self->o('release_number'),
      -host   => $self->o('rnaseq_blast_db_server'),
      -port   => $self->o('rnaseq_blast_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'rnaseq_refine_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_refine_'.$self->o('release_number'),
      -host   => $self->o('rnaseq_refine_db_server'),
      -port   => $self->o('rnaseq_refine_db_port'),
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


sub pipeline_create_commands {
  my ($self) = @_;

  return [
    # inheriting database and hive tables' creation
    @{$self->SUPER::pipeline_create_commands},

    'mkdir -p '.$self->o('rnaseq_dir'),
  ];
}


sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
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
      -rc_name    => 'default',

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
          'INSERT INTO analysis (logic_name, module, created, db_version) VALUES ("other_protein", "HiveBlastRNAseq", NOW(), "#uniprot_version#")',
          'UPDATE protein_align_feature paf, analysis a SET paf.analysis_id = a.analysis_id WHERE a.logic_name = "other_protein"',
          'DELETE FROM meta WHERE meta_key LIKE "assembly.web_accession%"',
          'DELETE FROM meta WHERE meta_key LIKE "removed_evidence_flag\.%"',
          'DELETE FROM meta WHERE meta_key LIKE "marker\.%"',
          'DELETE FROM meta WHERE meta_key IN ("genebuild.method","genebuild.projection_source_db","genebuild.start_date","repeat.analysis")',
          'UPDATE gene JOIN transcript USING(gene_id) SET canonical_transcript_id = transcript_id',
          'UPDATE transcript JOIN translation USING(transcript_id) SET canonical_translation_id = translation_id',
          'UPDATE intron_supporting_evidence ise, analysis a1, analysis a2 SET ise.analysis_id = a2.analysis_id WHERE ise.analysis_id = a1.analysis_id AND a2.logic_name = REPLACE(a1.logic_name, "rnaseq_gene", "rnaseq_ise")',
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
        id_start => 'RNASEQ',
        output_path => $self->o('output_path'),
        _stable_id_file => 'rnaseq_stable_ids.sql',
      },
      -rc_name    => '2GB',
      -flow_into => {
        1 => ['populate_production_tables_rnaseq'],
      },
    },

    {
      -logic_name => 'populate_production_tables_rnaseq',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HivePopulateProductionTables',
      -parameters => {
        'target_db'        => $self->o('rnaseq_db'),
        'output_path'      => $self->o('output_path'),
        'enscode_root_dir' => $self->o('enscode_root_dir'),
        'production_db'    => $self->o('production_db'),
      },
      -rc_name    => 'default',
      -flow_into  => {
        1 => ['dump_daf_introns'],
      },
    },

    {
      -logic_name => 'dump_daf_introns',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
      -parameters => {
        db_conn => $self->o('rnaseq_refine_db'),
        input_query => 'SELECT daf.* FROM dna_align_feature daf, analysis a WHERE daf.analysis_id = a.analysis_id AND a.logic_name != "rough_transcripts"',
        command_out => q(sort -nk2 -nk3 -nk4 | sed 's/NULL/\\N/g;s/^[0-9]\+/\\N/' | awk -F \t '{$15="NULL"; print $0}' > #daf_file#),
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
          ' -output_path '.catfile($self->o('rnaseq_dir'), 'optimise_rnaseq').
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
        1 => ['null_rnaseq_columns'],
      },
    },

    {
      -logic_name => 'null_rnaseq_columns',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('rnaseq_db'),
        sql => [
          'UPDATE dna_align_feature SET external_db_id = NULL',
        ],
      },
      -rc_name    => 'default',
      -flow_into => {
        1 => ['rnaseq_set_align_type'],
      },
    },

    {
      -logic_name => 'rnaseq_set_align_type',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('rnaseq_db'),
        sql => [
          'UPDATE dna_align_feature SET align_type = "ensembl"',
        ],
      },
      -rc_name    => 'default',
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
          'gene_db_checks')->{'rnaseq_final'},
      },

      -rc_name    => '4GB',
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
        cmd => 'TMPDIR=#working_dir# ; #deeptools_bamcoverage# --numberOfProcessors 4 --binSize 1 -b '.catfile('#working_dir#','#bam_file#').' -o '.catfile('#working_dir#','#bam_file#').'.bw',
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
        cmd => 'cd #working_dir#;md5sum #bam_file#.bw #bam_file# #bam_file#.bai > #bam_file#.md5',
        working_dir => $self->o('merge_dir'),
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'concat_md5_sum',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'cat '.catfile('#working_dir#', '*md5').' > '.catfile('#working_dir#', 'md5sum.txt.1'),
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
        FILES=($(ls *.bam));
        if [[ $((${#FILES[*]}-1)) > 0 ]]
        then printf "#free_text_multi#" | sed "s/NUM/$((${#FILES[*]}-1))/g;s/ \([a-z]\)\([a-z]\+_\)/ \U\1\E\2/;s/_/ /g" > README.1
        else printf "#free_text_single#" | sed "s/ \([a-z]\)\([a-z]\+_\)/ \U\1\E\2/;s/_/ /g" >> README.1
        fi
        IFS=$\'\n\';
        echo "${FILES[*]}" >> README.1',
        working_dir => $self->o('merge_dir'),
        species_name  => $self->o('species_name'),
        free_text_single => 'Note\n------\n\n'.
          'The RNASeq data for #species_name# consists of 1 individual sample.\n\n'.
          'The bam file has an index file (.bai) and a BigWig file (.bw) which contains the coverage information.\n\n'.
          'Use the md5sum.txt file to check the integrity of the downloaded files.\n\n'.
          'Files\n-----\n',
        free_text_multi => 'Note\n------\n\n'.
          'The RNASeq data for #species_name# consists of NUM individual samples and one merged set containing all NUM samples.\n\n'.
          'All files have an index file (.bai) and a BigWig file (.bw) which contains the coverage information.\n\n'.
          'Use the md5sum.txt file to check the integrity of the downloaded files.\n\n'.
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
    },
  ];
}


sub resource_classes {
  my $self = shift;

  return {
    '2GB' => { LSF => $self->lsf_resource_builder('production-rh74', 2000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '4GB' => { LSF => $self->lsf_resource_builder('production-rh74', 4000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '8GB' => { LSF => $self->lsf_resource_builder('production-rh74', 8000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '10GB_multithread' => { LSF => $self->lsf_resource_builder('production-rh74', 10000, [$self->default_options->{'pipe_db_server'}], undef, ($self->default_options->{'use_threads'}+1))},
    'default' => { LSF => $self->lsf_resource_builder('production-rh74', 900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
  }
}


1;
