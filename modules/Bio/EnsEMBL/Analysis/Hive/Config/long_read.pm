=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Analysis::Hive::Config::long_read;

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
    ########################
    # Misc setup info
    ########################
    'dbowner' => '' || $ENV{EHIVE_USER} || $ENV{USER},
    'pipeline_name' => '' || $self->o('production_name') . '_' . $self->o('ensembl_release'),
    'user_r'           => '',    # read only db user
    'user'             => '',    # write db user
    'password'         => '',    # password for write db user
    'pipe_db_host'     => '',    # host for pipe db
    'databases_host'   => '',    # host for general output dbs
    'dna_db_host'      => '',    # host for dna db
    'databases_port'   => '',    # port for general output db host

    'long_read_summary_file'       => '' || catfile( $self->o('long_read_dir'), $self->o('species_name') . '_long_read.csv' ),        # csv file for minimap2, should have 2 columns tab separated cols: sample_name\tfile_name
    'long_read_summary_file_genus' => '' || catfile( $self->o('long_read_dir'), $self->o('species_name') . '_long_read_gen.csv' ),    # csv file for minimap2, should have 2 columns tab separated cols: sample_name\tfile_name
    'long_read_fastq_dir' => '' || catdir( $self->o('long_read_dir'), 'input' ),

    'release_number' => '' || $self->o('ensembl_release'),
    'species_name'    => '',                                                                                                          # e.g. mus_musculus
    'production_name' => '',                                                                                                          # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
    dbname_accession          => '', # This is the assembly accession without [._] and all lower case, i.e gca001857705v1
    'output_path'     => '',                                                                                                          # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses
    'uniprot_version' => 'uniprot_2021_04',                                                                                           # What UniProt data dir to use for various analyses

    # Keys for custom loading, only set/modify if that's what you're doing
    'protein_blast_db' => '' || catfile( $self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), ($self->o('is_non_vert') eq '1') ? 'PE12' : 'PE12_vertebrata' ),    # Blast database for comparing the final models to.
    'protein_blast_index' => '' || catdir( $self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), ($self->o('is_non_vert') eq '1') ? 'PE12_index' : 'PE12_vertebrata_index' ),    # Indicate Index for the blast database.

    ########################
    # Pipe and ref db info
    ########################
    'long_read_initial_db_host'   => $self->o('databases_host'),
    'long_read_initial_db_port'   => $self->o('databases_port'),

    'long_read_collapse_db_host'   => $self->o('databases_host'),
    'long_read_collapse_db_port'   => $self->o('databases_port'),

    'long_read_final_db_host'   => $self->o('databases_host'),
    'long_read_final_db_port'   => $self->o('databases_port'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'ensembl_release' => $ENV{ENSEMBL_RELEASE},    # this is the current release version on staging to be able to get the correct database

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

    # This is used for "messaging" other sub pipeline
    transcript_selection_url => undef,


    ########################
    # Executable paths
    ########################
    'minimap2_genome_index' => $self->o('faidx_genome_file') . '.mmi',
    'minimap2_path'         => catfile( $self->o('binary_base'), 'minimap2' ),
    'paftools_path'         => catfile( $self->o('binary_base'), 'paftools.js' ),
    'minimap2_batch_size'   => '5000',

    'blast_type' => 'ncbi',    # It can be 'ncbi', 'wu', or 'legacy_ncbi'

    'uniprot_blast_exe_path' => catfile( $self->o('binary_base'), 'blastp' ),

    samtools_path => catfile( $self->o('binary_base'), 'samtools' ),    #You may need to specify the full path to the samtools binary

    'long_read_dir'       => catdir( $self->o('output_path'),   'long_read' ),
    'long_read_fastq_dir' => catdir( $self->o('long_read_dir'), 'input' ),

    use_threads => 3,

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
    long_read_columns => [ 'sample', 'filename', 'description', 'fastq_file', 'fastq_md5' ],
    download_method => 'ftp',
    databases_to_delete => ['long_read_initial_db', 'long_read_collapse_db', 'long_read_final_db'],

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # No option below this mark should be modified
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ########################
    # db info
    ########################
    long_read_initial_db => {
      -dbname => $self->o('dbowner') . '_' . $self->o('dbname_accession') . '_lrinitial_' . $self->o('release_number'),
      -host   => $self->o('long_read_initial_db_host'),
      -port   => $self->o('long_read_initial_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    long_read_collapse_db => {
      -dbname => $self->o('dbowner') . '_' . $self->o('dbname_accession') . '_lrcollapse_' . $self->o('release_number'),
      -host   => $self->o('long_read_collapse_db_host'),
      -port   => $self->o('long_read_collapse_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    long_read_final_db => {
      -dbname => $self->o('dbowner') . '_' . $self->o('dbname_accession') . '_lrfinal_' . $self->o('release_number'),
      -host   => $self->o('long_read_final_db_host'),
      -port   => $self->o('long_read_final_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

  };
}

sub pipeline_create_commands {
  my ($self) = @_;

  return [

    # inheriting database and hive tables' creation
    @{ $self->SUPER::pipeline_create_commands },

    'mkdir -p ' . $self->o('long_read_fastq_dir'),
  ];
}

sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{ $self->SUPER::pipeline_wide_parameters },
    genome_file         => $self->o('faidx_genome_file'),
    use_genome_flatfile => $self->o('use_genome_flatfile'),
    is_non_vert         => $self->o('is_non_vert'),
  }
}

## See diagram for pipeline structure
sub pipeline_analyses {
  my ($self) = @_;

  ########################################################################
  #
  # Minimap2 long read analyses
  #
  ########################################################################

  return [
    {
      -logic_name => 'create_fastq_dir',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name    => 'default',
      -parameters => {
        cmd => 'if [ ! -e "' . $self->o('long_read_fastq_dir') . '" ]; then mkdir -p ' . $self->o('long_read_fastq_dir') . ';fi',
      },
      -input_ids  => [{}],
      -flow_into => {
        '1' => ['create_initial_db'],
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'create_initial_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('dna_db'),
        target_db   => $self->o('long_read_initial_db'),
        create_type => 'clone',
      },
      -rc_name   => 'default',
      -flow_into => {
        '1' => ['create_minimap2_index'],
      },
    },

    {
      -logic_name => 'create_minimap2_index',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'if [ ! -e "' . $self->o('minimap2_genome_index') . '" ]; then ' . $self->o('minimap2_path') .
          ' -d ' . $self->o('minimap2_genome_index') . ' ' . $self->o('faidx_genome_file') . ';fi',
      },
      -flow_into => {
        1 => ['check_index_not_empty'],
      },
      -rc_name => '20GB',
    },

    {
      -logic_name => 'check_index_not_empty',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'if [ -s "' . $self->o('minimap2_genome_index') . '" ]; then exit 0; else exit 42;fi',
        return_codes_2_branches => { '42' => 2 },
      },
      -flow_into => {
        1 => ['create_fastq_download_jobs'],
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'create_fastq_download_jobs',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        inputfile    => $self->o('long_read_summary_file'),
        column_names => $self->o('long_read_columns'),
        delimiter    => '\t',
      },
      -flow_into => {
        '2->A' => { 'download_fastq' => { 'url' => '#fastq_file#' } },
        'A->1' => ['create_collapse_db'],
      },
    },

    {
      -logic_name => 'download_fastq',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadData',
      -parameters => {
        samtools => $self->o('samtools_path'),
        output_dir => $self->o('long_read_fastq_dir'),
        download_method => $self->o('download_method'),
        uncompress => 1,
        create_faidx  => 1,
      },
      -rc_name           => 'default',
      -analysis_capacity => 50,
      -flow_into         => {
        2 => { 'generate_minimap2_jobs' => { 'fastq_file' => '#filename#' } },
      },
    },

    {
      -logic_name => 'generate_minimap2_jobs',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        iid_type   => 'fastq_range',
        batch_size => $self->o('minimap2_batch_size'),
      },
      -rc_name   => 'default',
      -flow_into => {
        2 => { 'minimap2' => { 'input_file' => '#fastq_file#', 'iid' => '#iid#' } },
      },
    },

    {
      -logic_name => 'minimap2',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::Minimap2',
      -parameters => {
        genome_file                  => $self->o('faidx_genome_file'),
        long_read_summary_file       => $self->o('long_read_summary_file'),
        long_read_summary_file_genus => $self->o('long_read_summary_file_genus'),
        minimap2_genome_index        => $self->o('minimap2_genome_index'),
        minimap2_path                => $self->o('minimap2_path'),
        paftools_path                => $self->o('paftools_path'),
        target_db                    => $self->o('long_read_initial_db'),
        logic_name                   => 'minimap2',
        module                       => 'Minimap2',
      },
      -rc_name   => '15GB',
      -hive_capacity => $self->o('hc_normal'),
      -flow_into => {
        -1 => { 'minimap2_himem' => { 'input_file' => '#input_file#', 'iid' => '#iid#' } },
      },
    },

    {
      -logic_name => 'minimap2_himem',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::Minimap2',
      -parameters => {
        genome_file                  => $self->o('faidx_genome_file'),
        long_read_summary_file       => $self->o('long_read_summary_file'),
        long_read_summary_file_genus => $self->o('long_read_summary_file_genus'),
        minimap2_genome_index        => $self->o('minimap2_genome_index'),
        minimap2_path                => $self->o('minimap2_path'),
        paftools_path                => $self->o('paftools_path'),
        target_db                    => $self->o('long_read_initial_db'),
        logic_name                   => 'minimap2',
        module                       => 'Minimap2',
      },
      -hive_capacity => $self->o('hc_normal'),
      -rc_name => '25GB',
    },

    {
      -logic_name => 'create_collapse_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('long_read_initial_db'),
        target_db   => $self->o('long_read_collapse_db'),
        create_type => 'clone',
      },
      -rc_name         => 'default',
      -max_retry_count => 0,
      -flow_into       => {
        1 => ['create_final_db'],
        }
    },

    {
      -logic_name => 'create_final_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('long_read_initial_db'),
        target_db   => $self->o('long_read_final_db'),
        create_type => 'clone',
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['generate_collapse_jobs'],
      },
    },

    {
      -logic_name => 'generate_collapse_jobs',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db          => $self->o('long_read_collapse_db'),
        feature_dbs        => [ $self->o('long_read_initial_db') ],
        coord_system_name  => 'toplevel',
        iid_type           => 'stranded_slice',
        feature_constraint => 1,
        feature_type       => 'gene',
        top_level          => 1,
      },
      -rc_name         => 'default',
      -max_retry_count => 1,
      -flow_into       => {
        '2->A' => ['split_slices_on_intergenic'],
        'A->1' => ['classify_models'],
      },
    },

    {
      -logic_name => 'split_slices_on_intergenic',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveFindIntergenicRegions',
      -parameters => {
        dna_db         => $self->o('dna_db'),
        input_gene_dbs => [ $self->o('long_read_initial_db') ],
        iid_type       => 'slice',
        use_strand     => 1,
      },
      -batch_size => 100,
      -hive_capacity => $self->o('hc_normal'),
      -rc_name    => '5GB',
      -flow_into  => {
        2 => { 'collapse_transcripts' => { 'slice_strand' => '#slice_strand#', 'iid' => '#iid#' } },
      },
    },

    {
      -logic_name => 'collapse_transcripts',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveTranscriptCoalescer',
      -parameters => {
        target_db             => $self->o('long_read_collapse_db'),
        dna_db                => $self->o('dna_db'),
        source_dbs            => [ $self->o('long_read_initial_db') ],
        biotypes              => [ "isoseq", "cdna" ],
        reduce_large_clusters => 1,
      },
      -rc_name   => '5GB',
      -flow_into => {
        1 => ['blast'],
        -1 => { 'collapse_transcripts_20GB' => { 'slice_strand' => '#slice_strand#', 'iid' => '#iid#' } },
      },
      -batch_size        => 100,
      -hive_capacity => $self->o('hc_normal'),
    },

    {
      -logic_name => 'blast',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastRNASeqPep',
      -parameters => {
        input_db       => $self->o('long_read_collapse_db'),
        output_db      => $self->o('long_read_final_db'),
        source_db      => $self->o('long_read_collapse_db'),
        target_db      => $self->o('long_read_final_db'),
        dna_db         => $self->o('dna_db'),
        indicate_index => $self->o('protein_blast_index'),
        uniprot_index  => [ $self->o('protein_blast_db') ],
        blast_program  => $self->o('uniprot_blast_exe_path'),
        %{ get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic', 'BlastGenscanPep', { BLAST_PARAMS => { -type => $self->o('blast_type') } } ) },
        commandline_params => $self->o('blast_type') eq 'wu' ? '-cpus=' . $self->o('use_threads') . ' -hitdist=40' : '-num_threads ' . $self->o('use_threads') . ' -window_size 40 -seg no',
      },
      -rc_name   => 'blast',
      -hive_capacity => $self->o('hc_normal'),
      -flow_into => {
        -1 => ['blast_10G'],
      },
    },

    {
      -logic_name => 'blast_10G',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastRNASeqPep',
      -parameters => {
        input_db       => $self->o('long_read_collapse_db'),
        output_db      => $self->o('long_read_final_db'),
        dna_db         => $self->o('dna_db'),
        source_db      => $self->o('long_read_collapse_db'),
        target_db      => $self->o('long_read_final_db'),
        indicate_index => $self->o('protein_blast_index'),
        uniprot_index  => [ $self->o('protein_blast_db') ],
        blast_program  => $self->o('uniprot_blast_exe_path'),
        %{ get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic', 'BlastGenscanPep', { BLAST_PARAMS => { -type => $self->o('blast_type') } } ) },
        commandline_params => $self->o('blast_type') eq 'wu' ? '-cpus=' . $self->o('use_threads') . ' -hitdist=40' : '-num_threads ' . $self->o('use_threads') . ' -window_size 40 -seg no',
      },
      -hive_capacity => $self->o('hc_normal'),
      -rc_name => 'blast10GB',
    },

    {
      -logic_name => 'collapse_transcripts_20GB',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveTranscriptCoalescer',
      -parameters => {
        target_db             => $self->o('long_read_collapse_db'),
        dna_db                => $self->o('dna_db'),
        source_dbs            => [ $self->o('long_read_initial_db') ],
        biotypes              => [ "isoseq", "cdna" ],
        reduce_large_clusters => 1,
      },
      -rc_name   => '20GB',
      -flow_into => {
        1  => { 'blast' => { 'slice_strand' => '#slice_strand#', 'iid' => '#iid#' } },
        -1 => { 'failed_collapse' => { 'slice_strand' => '#slice_strand#', 'iid' => '#iid#' } },
      },
      -batch_size        => 10,
      -hive_capacity => $self->o('hc_normal'),
    },

    {
      -logic_name => 'failed_collapse',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveTranscriptCoalescer',
      -parameters => {
        target_db  => $self->o('long_read_collapse_db'),
        dna_db     => $self->o('dna_db'),
        source_dbs => [ $self->o('long_read_initial_db') ],
        biotypes   => [ "isoseq", "cdna" ],
        copy_only  => 1,
      },
      -hive_capacity => $self->o('hc_normal'),
      -rc_name   => '10GB',
      -flow_into => {
        1 => { 'blast' => { 'slice_strand' => '#slice_strand#', 'iid' => '#iid#' } },
      },
    },

    {
      -logic_name => 'classify_models',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport',
      -parameters => {
        classification_type => 'long_read',
        update_gene_biotype => 1,
        target_db           => $self->o('long_read_final_db'),
      },
      -rc_name => 'default',
      -flow_into  => {
        '1'  => ['delete_long_reads'],
      },
    },

   {
     -logic_name => 'delete_long_reads',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'rm '.$self->o('long_read_fastq_dir').'/*',
      },
      -rc_name => 'default',
      -flow_into  => {
	'1'  => ['notification_pipeline_is_done'],
      }
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
          data => $self->o('long_read_final_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'split_slices_on_intergenic',
          param => 'input_gene_dbs',
          data => $self->o('long_read_final_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'layer_annotation',
          param => 'SOURCEDB_REFS',
          data => $self->o('long_read_final_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'run_utr_addition',
          param => 'donor_dbs',
          data => $self->o('long_read_final_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'run_utr_addition_10GB',
          param => 'donor_dbs',
          data => $self->o('long_read_final_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'run_utr_addition_30GB',
          param => 'donor_dbs',
          data => $self->o('long_read_final_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'utr_memory_failover',
          param => 'donor_dbs',
          data => $self->o('long_read_final_db'),
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
    'default' => { LSF => $self->lsf_resource_builder( 'production', 900 ) },
    '5GB'     => { LSF => $self->lsf_resource_builder( 'production', 5000 ) },
    '10GB'    => { LSF => $self->lsf_resource_builder( 'production', 10000 ) },
    '15GB'    => { LSF => $self->lsf_resource_builder( 'production', 15000 ) },
    '20GB'    => { LSF => $self->lsf_resource_builder( 'production', 20000 ) },
    '25GB'    => { LSF => $self->lsf_resource_builder( 'production', 25000 ) },
    'blast'     => { LSF => $self->lsf_resource_builder( 'production', 2900, undef, 3 ) },
    'blast10GB' => { LSF => $self->lsf_resource_builder( 'production', 10000, undef, undef, 3 ) },
    }
}

1;
