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

package Bio::EnsEMBL::Analysis::Hive::Config::OtherFeatureDb;

use strict;
use warnings;
use File::Spec::Functions;

use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_analysis_settings);
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf qw(INPUT_PLUS);
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
    'user_r'                    => '',                                  # read only db user
    'user'                      => '',                                  # write db user
    'password'                  => '',                                  # password for write db user
    'pipe_db_host'              => '',                                  # host for pipe db
    'databases_host'            => '',                                  # host for general output dbs
    'dna_db_host'               => '',                                  # host for dna db
    'pipe_db_port'              => '',                                  # port for pipeline host
    'databases_port'            => '',                                  # port for general output db host
    'dna_db_port'               => '',                                  # port for dna db host
    'registry_host'             => '',                                  # host for registry db
    'registry_port'             => '',                                  # port for registry db
    'registry_db'               => '',                                  # name for registry db
    'release_number'            => '' || $self->o('ensembl_release'),
    'species_name'              => '',                                  # e.g. mus_musculus
    'production_name'           => '',                                  # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
    'output_path'               => '',                                  # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses
    'assembly_name'             => '',                                  # Name (as it appears in the assembly report file)
    'assembly_accession'        => '',                                  # Versioned GCA assembly accession, e.g. GCA_001857705.1
    'uniprot_version'           => 'uniprot_2019_04',                   # What UniProt data dir to use for various analyses

    # Keys for custom loading, only set/modify if that's what you're doing
    'protein_entry_loc' => catfile( $self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), 'entry_loc' ),    # Used by genscan blasts and optimise daf/paf. Don't change unless you know what you're doing

########################
    # Pipe and ref db info
########################

    'pipe_db_name' => $self->o('dbowner') . '_' . $self->o('production_name') . '_pipe_' . $self->o('release_number'),
    'dna_db_name'  => $self->o('dbowner') . '_' . $self->o('production_name') . '_core_' . $self->o('release_number'),

    otherfeatures_db_name => $self->o('dbowner').'_'.$self->o('production_name').'_otherfeatures_'.$self->o('release_number'),
    'otherfeatures_db_host'   => $self->o('databases_host'),
    'otherfeatures_db_port'   => $self->o('databases_port'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'ensembl_release'      => $ENV{ENSEMBL_RELEASE},     # this is the current release version on staging to be able to get the correct database
    'production_db_host'   => 'mysql-ens-meta-prod-1',
    'production_db_port'   => '4483',

    cdna_prefix => 'cdna',
    long_read_prefix => 'lrinitial',
    refseq_prefix => 'refseq',
    databases_to_check => [
      [$self->o('cdna_prefix'), $self->o('cdna_logic_name'), 'cDNA'],
      [$self->o('refseq_prefix'), undef, undef],
      [$self->o('long_read_prefix'), undef, 'cDNA'],
    ],
    databases_to_delete => ['otherfeatures_db'],    #, 'projection_realign_db'

########################
    # BLAST db paths
########################
    'base_blast_db_path' => $ENV{BLASTDB_DIR},

######################################################
    #
    # Mostly constant settings
    #
######################################################

    genome_dumps => catdir( $self->o('output_path'), 'genome_dumps' ),
    # This one is used in replacement of the dna table in the core db, so where analyses override slice->seq. Has simple headers with just the seq_region name. Also used by bwa in the RNA-seq analyses. Not masked
    faidx_genome_file => catfile( $self->o('genome_dumps'), $self->o('species_name') . '_toplevel.fa' ),
    use_genome_flatfile => 1,

    min_threshold_genes => 10,
    cdna_logic_name => 'cdna_alignment',
    ensembl_analysis_script => catdir( $self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts' ),
    load_optimise_script => catfile( $self->o('ensembl_analysis_script'), 'genebuild', 'load_external_db_ids_and_optimize_af.pl' ),
    add_analyses_script => catfile( $self->o('ensembl_analysis_script'), 'genebuild', 'add_rnaseq_analysis_descriptions.pl' ),

    ensembl_misc_script => catdir( $self->o('enscode_root_dir'), 'ensembl', 'misc-scripts' ),
    meta_coord_script => catfile( $self->o('ensembl_misc_script'), 'meta_coord', 'update_meta_coord.pl' ),
    meta_levels_script       => catfile( $self->o('ensembl_misc_script'),     'meta_levels.pl' ),
    frameshift_attrib_script => catfile( $self->o('ensembl_misc_script'),     'frameshift_transcript_attribs.pl' ),
    select_canonical_script  => catfile( $self->o('ensembl_misc_script'),     'canonical_transcripts', 'select_canonical_transcripts.pl' ),
    assembly_name_script     => catfile( $self->o('ensembl_analysis_script'), 'update_assembly_name.pl' ),


########################
    # db info
########################
    'production_db' => {
      -host   => $self->o('production_db_host'),
      -port   => $self->o('production_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -dbname => 'ensembl_production',
      -driver => $self->o('hive_driver'),
    },

    'otherfeatures_db' => {
      -dbname => $self->o('otherfeatures_db_name'),
      -host   => $self->o('otherfeatures_db_host'),
      -port   => $self->o('otherfeatures_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

  };
}

sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
    genome_file         => $self->o('faidx_genome_file'),
    use_genome_flatfile => $self->o('use_genome_flatfile'),
  }
}


## See diagram for pipeline structure
sub pipeline_analyses {
  my ($self) = @_;

  return [

    {
      -logic_name => 'create_otherfeatures_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('dna_db'),
        target_db   => $self->o('otherfeatures_db'),
        create_type => 'clone',
      },
      -input_ids  => [{}],
      -rc_name   => 'default',
      -flow_into => {
        1 => ['prepare_otherfeatures_db'],
      },
    },

    {
      -logic_name => 'prepare_otherfeatures_db',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('otherfeatures_db'),
        sql     => [
          'TRUNCATE analysis',
          'DELETE FROM meta WHERE meta_key LIKE "assembly.web_accession%"',
          'DELETE FROM meta WHERE meta_key LIKE "sample.%"',
          'DELETE FROM meta WHERE meta_key IN'.
            ' ("repeat.analysis", "genebuild.last_geneset_update","genebuild.method","genebuild.projection_source_db","genebuild.start_date")',
          'INSERT INTO meta (species_id,meta_key,meta_value) VALUES (1,"genebuild.last_otherfeatures_update",NOW())',
        ],
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['create_databases_jobs'],
      },
    },
    {
      -logic_name => 'create_databases_jobs',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        inputlist    => $self->o('databases_to_check'),
        column_names => ['prefix', 'logic_name', 'biotype'],
      },
      -rc_name   => 'default',
      -flow_into => {
        '2->A' => ['fan_database_exists'],
        'A->1' => ['set_otherfeatures_meta_coords'],
      },
    },

    {
      -logic_name => 'fan_database_exists',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        databases_host => $self->o('databases_host'),
        databases_port => $self->o('databases_port'),
        read_only_user => $self->o('user_r'),
        dbowner => $self->o('dbowner'),
        release_number => $self->o('release_number'),
        production_name => $self->o('production_name'),
        cmd => q(if [[ -n `mysql -h #databases_host# -P #databases_port# -u #read_only_user# -NB -e "SHOW DATABASES LIKE '#dbowner#_#production_name#_#prefix#_#release_number#'"` ]]; then exit 0; else exit 42;fi),

        return_codes_2_branches => {'42' => 2},
      },
      -batch_size => 5,
      -rc_name => 'default',
      -flow_into  => {
        1 => ['fan_has_enough_genes'],
      },
    },

    {
      -logic_name => 'fan_has_enough_genes',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        databases_host => $self->o('databases_host'),
        databases_port => $self->o('databases_port'),
        read_only_user => $self->o('user_r'),
        dbowner => $self->o('dbowner'),
        release_number => $self->o('release_number'),
        production_name => $self->o('production_name'),
        min_threshold => $self->o('min_threshold_genes'),,
        cmd => q(if [[ `mysql -h #databases_host# -P #databases_port# -u #read_only_user# #dbowner#_#production_name#_#prefix#_#release_number# -NB -e "SELECT COUNT(*) FROM gene"` -gt #min_threshold# ]]; then exit 0; else exit 42;fi),

        return_codes_2_branches => {'42' => 2},
      },
      -batch_size => 5,
      -rc_name => 'default',
      -flow_into  => {
        1 => ['create_gene_ids_to_copy'],
      },
    },

    {
      -logic_name => 'create_gene_ids_to_copy',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db => {
          -dbname => $self->o('dbowner') . '_' . $self->o('production_name') . '_#prefix#_' . $self->o('release_number'),
          -host   => $self->o('databases_host'),
          -port   => $self->o('databases_port'),
          -user   => $self->o('user'),
          -pass   => $self->o('password'),
          -driver => $self->o('hive_driver'),
        },
        iid_type     => 'feature_id',
        feature_type => 'gene',
        batch_size   => 500,
      },
      -flow_into => {
        '2->A' => {'copy_genes_to_otherfeatures' => INPUT_PLUS()},
        'A->1' => ['dummy_semaphore_fan'],
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'copy_genes_to_otherfeatures',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCopyGenes',
      -parameters => {
        copy_genes_directly => 1,
        logic_name => '#logic_name#',
        biotype => '#biotype#',
        source_db           => {
          -dbname => $self->o('dbowner') . '_' . $self->o('production_name') . '_#prefix#_' . $self->o('release_number'),
          -host   => $self->o('databases_host'),
          -port   => $self->o('databases_port'),
          -user   => $self->o('user'),
          -pass   => $self->o('password'),
          -driver => $self->o('hive_driver'),
        },
        dna_db              => $self->o('dna_db'),
        target_db           => $self->o('otherfeatures_db'),
      },
      -rc_name => 'default',
      -batch_size => 30,
      -hive_capacity => 200,
    },

    {
      -logic_name => 'dummy_semaphore_fan',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -rc_name   => 'default',
      -flow_into  => {
        1 => [
          'fan_long_read_analyses',
          'fan_cdna_analyses'
        ],
      },
    },

    {
      -logic_name => 'fan_long_read_analyses',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'if [[ "#prefix#" = "'.$self->o('long_read_prefix').'" ]]; then exit 0; else exit 42;fi',
        return_codes_2_branches => {'42' => 2},
      },
      -rc_name => 'default',
      -flow_into  => {
        1 => ['update_long_read_data'],
      },
    },

    {
      -logic_name => 'update_long_read_data',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('otherfeatures_db'),
        stable_id_pattern => 'SRR',
        sql     => [
          # The version is set to 1 as default
          'UPDATE gene JOIN analysis USING(analysis_id) SET version = 1 WHERE logic_name LIKE "%_isoseq"',
          'UPDATE transcript JOIN analysis USING(analysis_id) SET version = 1 WHERE logic_name LIKE "%_isoseq"',
          'UPDATE gene g, transcript t, analysis a SET g.canonical_transcript_id = t.transcript_id'.
            ' WHERE g.gene_id = t.gene_id AND g.analysis_id = a.analysis_id AND a.logic_name LIKE "%_isoseq"',
        ],
      },
      -rc_name   => 'default',
      -flow_into  => {
        1 => ['add_long_read_analyses'],
      },
    },

    {
      -logic_name => 'add_long_read_analyses',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('add_analyses_script').' -host #host# -port #port# -dbname #dbname#',
        dbname => $self->o('dbowner') . '_' . $self->o('production_name') . '_#prefix#_' . $self->o('release_number'),
        host   => $self->o('databases_host'),
        port   => $self->o('databases_port'),
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'fan_cdna_analyses',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'if [[ "#prefix#" = "'.$self->o('cdna_prefix').'" ]]; then exit 0; else exit 42;fi',
        return_codes_2_branches => {'42' => 2},
      },
      -rc_name => 'default',
      -flow_into  => {
        1 => ['update_cdna_data'],
      },
    },

    {
      -logic_name => 'update_cdna_data',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('otherfeatures_db'),
        logic_name => $self->o('cdna_logic_name'),
        sql     => [
          'UPDATE gene g, transcript t, analysis a SET g.canonical_transcript_id = t.transcript_id'.
            ' WHERE g.gene_id = t.gene_id AND g.analysis_id = a.analysis_id AND a.logic_name = "#logic_name#"',
          'UPDATE transcript JOIN transcript_supporting_feature USING(transcript_id)' .
            ' JOIN dna_align_feature ON feature_id = dna_align_feature_id'.
            ' JOIN analysis ON analysis.analysis_id = transcript.analysis_id'.
            ' SET stable_id = SUBSTRING_INDEX(hit_name, ".", 1), version = SUBSTRING_INDEX(hit_name, ".", -1)'.
            ' WHERE feature_type = "dna_align_feature" AND logic_name = "#logic_name#"',
          'UPDATE gene JOIN transcript ON transcript_id = canonical_transcript_id' .
            ' JOIN analysis ON analysis.analysis_id = transcript.analysis_id'.
            ' SET gene.stable_id = transcript.stable_id, gene.version = transcript.version WHERE logic_name = "#logic_name#"',
        ],
      },
      -rc_name   => 'default',
    },

    {
      -logic_name => 'set_otherfeatures_meta_coords',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . $self->o('meta_coord_script') .
          ' -user ' . $self->o('user') .
          ' -pass ' . $self->o('password') .
          ' -host ' . $self->o( 'otherfeatures_db', '-host' ) .
          ' -port ' . $self->o( 'otherfeatures_db', '-port' ) .
          ' -dbpattern ' . $self->o( 'otherfeatures_db', '-dbname' )
      },
      -rc_name   => '4GB',
      -flow_into => {
        1 => ['set_otherfeatures_meta_levels'],
      },
    },

    {
      -logic_name => 'set_otherfeatures_meta_levels',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . $self->o('meta_levels_script') .
          ' -user ' . $self->o('user') .
          ' -pass ' . $self->o('password') .
          ' -host ' . $self->o( 'otherfeatures_db', '-host' ) .
          ' -port ' . $self->o( 'otherfeatures_db', '-port' ) .
          ' -dbname ' . $self->o( 'otherfeatures_db', '-dbname' )
      },
      -rc_name => '4GB',
      -flow_into => { 1 => ['null_otherfeatures_columns'] },
    },

    {
      -logic_name => 'null_otherfeatures_columns',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('otherfeatures_db'),
        sql     => [
          'UPDATE dna_align_feature SET external_db_id = NULL',
          'UPDATE protein_align_feature SET external_db_id = NULL',
        ],
      },
      -rc_name   => '4GB',
      -flow_into => {
        1 => ['load_external_db_ids_and_optimise_otherfeatures'],
      },
    },

    {
      -logic_name => 'load_external_db_ids_and_optimise_otherfeatures',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . $self->o('load_optimise_script') .
          ' -output_path ' . catdir( $self->o('output_path'), 'optimise_otherfeatures' ) .
          ' -uniprot_filename ' . $self->o('protein_entry_loc') .
          ' -dbuser ' . $self->o('user') .
          ' -dbpass ' . $self->o('password') .
          ' -dbport ' . $self->o( 'otherfeatures_db', '-port' ) .
          ' -dbhost ' . $self->o( 'otherfeatures_db', '-host' ) .
          ' -dbname ' . $self->o( 'otherfeatures_db', '-dbname' ) .
          ' -prod_dbuser ' . $self->o('user_r') .
          ' -prod_dbhost ' . $self->o( 'production_db', '-host' ) .
          ' -prod_dbname ' . $self->o( 'production_db', '-dbname' ) .
          ' -prod_dbport ' . $self->o( 'production_db', '-port' ) .
          ' -verbose'
      },
      -max_retry_count => 0,
      -rc_name         => '4GB',
      -flow_into       => {
        1 => ['otherfeatures_sanity_checks'],
      },
    },

    {
      -logic_name => 'otherfeatures_sanity_checks',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
      -parameters => {
        target_db                  => $self->o('otherfeatures_db'),
        sanity_check_type          => 'gene_db_checks',
        min_allowed_feature_counts => get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
          'gene_db_checks' )->{'otherfeatures'},
      },
      -rc_name   => '4GB',
      -flow_into => {
        1 => ['otherfeatures_healthchecks'],
      },
    },

    {
      -logic_name => 'otherfeatures_healthchecks',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveHealthcheck',
      -parameters => {
        input_db => $self->o('otherfeatures_db'),
        species  => $self->o('species_name'),
        group    => 'otherfeatures_handover',
      },
      -max_retry_count => 0,
      -rc_name   => '4GB',
      -flow_into => {
        1 => ['otherfeatures_assembly_name_update'],
      },
    },

    {
      -logic_name => 'otherfeatures_assembly_name_update',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . $self->o('assembly_name_script') .
          ' -user ' . $self->o('user') .
          ' -pass ' . $self->o('password') .
          ' -host ' . $self->o( 'otherfeatures_db', '-host' ) .
          ' -port ' . $self->o( 'otherfeatures_db', '-port' ) .
          ' -dbname ' . $self->o( 'otherfeatures_db', '-dbname' ) .
          ' -driver ' . $self->o('hive_driver') .
          ' -assembly_accession ' . $self->o('assembly_accession') .
          ' -assembly_name ' . $self->o('assembly_name') .
          ' -registry_host ' . $self->o('registry_host') .
          ' -registry_port ' . $self->o('registry_port') .
          ' -registry_db ' . $self->o('registry_db'),
      },
      -rc_name => '4GB',
    },

  ];
}

sub resource_classes {
  my $self = shift;
  return {
    '4GB'     => { LSF => $self->lsf_resource_builder( 'production', 4000 ) },
    'default' => { LSF => $self->lsf_resource_builder( 'production', 900 ) },
    }
}

1;
