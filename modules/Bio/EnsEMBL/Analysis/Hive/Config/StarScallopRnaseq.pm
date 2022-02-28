
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

package Bio::EnsEMBL::Analysis::Hive::Config::StarScallopRnaseq;

use strict;
use warnings;
use File::Spec::Functions;

use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_analysis_settings);
use base ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

# this is required for eHive's WHEN ELSE
use  Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use parent ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

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
    'user_r'                    => '',                                                                                # read only db user
    'user'                      => '',                                                                                # write db user
    'password'                  => '',                                                                                # password for write db user
    'server_set'                => '',                                                                                # What server set to user, e.g. set1
    'pipe_db_host'              => '',                                                                                # host for pipe db
    'databases_host'            => '',                                                                                # host for general output dbs
    'dna_db_host'               => '',                                                                                # host for dna db
    'pipe_db_port'              => '',                                                                                # port for pipeline host
    'databases_port'            => '',                                                                                # port for general output db host
    'dna_db_port'               => '',                                                                                # port for dna db host
    'rnaseq_summary_file'       => '' || catfile( $self->o('rnaseq_dir'), $self->o('species_name') . '.csv' ),        # Set this if you have a pre-existing cvs file with the expected columns
    'rnaseq_summary_file_genus' => '' || catfile( $self->o('rnaseq_dir'), $self->o('species_name') . '_gen.csv' ),    # Set this if you have a pre-existing genus level cvs file with the expected columns
    rnaseq_study_accession => '',
    'rnasamba_tsv_file'         => '' || catfile( $self->o('pcp_dir'), $self->o('pcp_db_name') . '_RNAsamba.tsv' ),
    'cpc2_fasta_file'           => '' || catfile( $self->o('pcp_dir'), $self->o('pcp_db_name') . '.fasta' ),
    'cpc2_file'                 => '' || catfile( $self->o('pcp_dir'), $self->o('pcp_db_name') . '_cpc2' ),
    'cpc2_txt_file'             => '' || catfile( $self->o('pcp_dir'), $self->o('pcp_db_name') . '_cpc2.txt' ),
    'release_number' => '' || $self->o('ensembl_release'),
    'is_non_vert'    => '' || $self->o('is_non_vert'),
    'protein_blast_db_file' => 'PE12_vertebrata', # use PE12 for non-vertebrates'. Note there must also be a PE12_index file in the same directory.
    'protein_entry_loc_file' => 'entry_loc',
    'species_name'        => '',                                                                                      # e.g. mus_musculus
    'production_name'     => '',                                                                                      # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
    'taxon_id'            => '',                                                                                      # should be in the assembly report file
    'genus_taxon_id'      => $self->o('taxon_id'),
    'uniprot_set'         => '',                                                                                      # e.g. mammals_basic, check UniProtCladeDownloadStatic.pm module in hive config dir for suitable set,
    'output_path'         => '',                                                                                      # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses
    'assembly_name'       => '',                                                                                      # Name (as it appears in the assembly report file)
    'uniprot_version'     => 'uniprot_2021_04',                                                                       # What UniProt data dir to use for various analyses
    'paired_end_only'     => '1',                                                                                     # Will only use paired-end rnaseq data if 1

    # Keys for custom loading, only set/modify if that's what you're doing
    'protein_blast_db' => '' || catfile( $self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), $self->o('protein_blast_db_file') ), # Blast database for comparing the final models to.
    'protein_blast_index' => '' || catdir( $self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), $self->o('protein_blast_db_file').'_index' ), # Indicate Index for the blast database.
    'protein_entry_loc' => catfile( $self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), $self->o('protein_entry_loc_file') ), # Used by genscan blasts and optimise daf/paf. Don't change unless you know what you're doing

########################
# Pipe and ref db info
########################

    'pipe_db_name' => $self->o('dbowner') . '_' . $self->o('production_name') . '_scallop_pipe_' . $self->o('release_number'),
    'dna_db_name'  => $self->o('dbowner') . '_' . $self->o('production_name') . '_core_' . $self->o('release_number'),

    'reference_db_name'   => $self->o('dna_db_name'),
    'reference_db_host'   => $self->o('dna_db_host'),
    'reference_db_port'   => $self->o('dna_db_port'),

    'rnaseq_for_layer_db_host'   => $self->o('databases_host'),
    'rnaseq_for_layer_db_port'   => $self->o('databases_port'),
    'rnaseq_for_layer_db_name'   => $self->o('dbowner') . '_' . $self->o('production_name') . '_star_rs_layer_' . $self->o('release_number'),

    'rnaseq_for_layer_nr_db_host'   => $self->o('databases_host'),
    'rnaseq_for_layer_nr_db_port'   => $self->o('databases_port'),
    'rnaseq_for_layer_nr_db_name'   => $self->o('dbowner') . '_' . $self->o('production_name') . '_star_rs_layer_nr_' . $self->o('release_number'),

    'scallop_initial_db_host'   => $self->o('databases_host'),
    'scallop_initial_db_port'   => $self->o('databases_port'),

    'scallop_blast_db_host'   => $self->o('databases_host'),
    'scallop_blast_db_port'   => $self->o('databases_port'),
    
    'pcp_db_host' => $self->o('databases_host'),
    'pcp_db_port' => $self->o('databases_port'),
    'pcp_db_name' => $self->o('dbowner').'_'.$self->o('production_name').'_pcp_'.$self->o('release_number'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'ensembl_release' => $ENV{ENSEMBL_RELEASE},    # this is the current release version on staging to be able to get the correct database

    databases_to_delete => ['rnaseq_for_layer_db','rnaseq_for_layer_nr_db','scallop_initial_db','scallop_blast_db','pcp_db','pcp_nr_db'],

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

    'min_toplevel_slice_length' => 250,
    rnaseq_merge_type => 'samtools',
    rnaseq_merge_threads => 12,

    # This is used for "messaging" other sub pipeline
    main_pipeline_url => undef,
    transcript_selection_url => undef,
    homology_rnaseq_url => undef,
    
    cpc2_output => catdir($self->o('output_path'),'cpc2_output'),
    rnasamba_output => catdir($self->o('output_path'),'rnasamba_output'),
    rna_samba_weights => '/nfs/production/flicek/ensembl/genebuild/rnasamba/full_length_weights.hdf5',

########################
# Executable paths
########################
    star_path       => catfile($self->o('binary_base'), 'STAR'),
    scallop_path    => catfile($self->o('binary_base'), 'scallop'),
    stringtie2_path => catfile($self->o('binary_base'), 'stringtie'),
    samtools_path   => catfile($self->o('binary_base'), 'samtools'), #You may need to specify the full path to the samtools binary
    picard_lib_jar  => catfile($self->o('linuxbrew_home_path'), 'Cellar', 'picard-tools', '2.6.0', 'libexec', 'picard.jar'), #You need to specify the full path to the picard library
    rnasamba => '/hps/software/users/ensembl/genebuild/singularity/rnasamba_latest.sif',
    cpc2 => '/hps/software/users/ensembl/genebuild/singularity/test_cpc2.sif',
    ensembl_analysis_scripts   => catdir($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts'),
    pcp_get_transcripts_script => catfile($self->o('ensembl_analysis_scripts'), 'pcp', 'get_transcripts.pl'),

    'blast_type'             => 'ncbi',                                                                         # It can be 'ncbi', 'wu', or 'legacy_ncbi'
    'uniprot_blast_exe_path' => catfile( $self->o('binary_base'), 'blastp' ),

    'summary_file_delimiter' => '\t',            # Use this option to change the delimiter for your summary data file
    'summary_csv_table'      => 'csv_data',
    'read_length_table'      => 'read_length',
    'rnaseq_data_provider'   => 'ENA',           #It will be set during the pipeline or it will use this value

    'rnaseq_dir' => catdir( $self->o('output_path'), 'rnaseq' ),
    'input_dir'  => catdir( $self->o('rnaseq_dir'),  'input' ),
    'output_dir' => catdir( $self->o('rnaseq_dir'),  'output' ),
    'merge_dir' => catdir( $self->o('rnaseq_dir'),  'merge' ),
    'pcp_dir' => catdir( $self->o('rnaseq_dir'),  'pcp' ),

    'rnaseq_ftp_base' => 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/',

    # Regular expression to allow FastQ files to be correctly paired,
    # for example: file_1.fastq and file_2.fastq could be paired using
    # the expression "\S+_(\d)\.\S+".  Need to identify the read number
    # in brackets; the name the read number (1, 2) and the
    # extension.
    pairing_regex => '\S+_(\d)\.\S+',

    # What Read group tag would you like to group your samples
    # by? Default = ID
    read_group_tag => 'SM',
    read_id_tag    => 'ID',

    use_threads => 3,

    star_threads    => 12,
    scallop_threads => 2,

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
    file_columns => [ 'SM', 'ID', 'is_paired', 'filename', 'is_mate_1', 'read_length', 'is_13plus', 'CN', 'PL', 'DS', 'url', 'md5sum' ],
    download_method => 'ftp',

    'filename_tag' => 'filename',    # For the analysis that creates star jobs, though I assume we should need to do it this way

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# No option below this mark should be modified
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#

########################
# db info
########################
    'rnaseq_for_layer_db' => {
      -dbname => $self->o('rnaseq_for_layer_db_name'),
      -host   => $self->o('rnaseq_for_layer_db_host'),
      -port   => $self->o('rnaseq_for_layer_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'rnaseq_for_layer_nr_db' => {
      -dbname => $self->o('rnaseq_for_layer_nr_db_name'),
      -host   => $self->o('rnaseq_for_layer_nr_db_host'),
      -port   => $self->o('rnaseq_for_layer_nr_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'scallop_initial_db' => {
      -dbname => $self->o('dbowner') . '_' . $self->o('production_name') . '_scallop_initial_' . $self->o('release_number'),
      -host   => $self->o('scallop_initial_db_host'),
      -port   => $self->o('scallop_initial_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'scallop_blast_db' => {
      -dbname => $self->o('dbowner') . '_' . $self->o('production_name') . '_scallop_blast_' . $self->o('release_number'),
      -host   => $self->o('scallop_blast_db_host'),
      -port   => $self->o('scallop_blast_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'pcp_db'=> {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_pcp_'.$self->o('release_number'),
      -host   => $self->o('pcp_db_host'),
      -port   => $self->o('pcp_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'pcp_nr_db'=> {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_pcp_nr_'.$self->o('release_number'),
      -host   => $self->o('pcp_db_host'),
      -port   => $self->o('pcp_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

  };
}

sub pipeline_create_commands {
  my ($self) = @_;

  my $tables;
  my %small_columns = (
    paired      => 1,
    read_length => 1,
    is_13plus   => 1,
    is_mate_1   => 1,
  );
  # We need to store the values of the csv file to easily process it. It will be used at different stages
  foreach my $key ( @{ $self->default_options->{'file_columns'} } ) {
    if ( exists $small_columns{$key} ) {
      $tables .= $key . ' SMALLINT UNSIGNED NOT NULL,'
    }
    elsif ( $key eq 'DS' or $key eq 'url') {
      $tables .= $key . ' VARCHAR(255) NOT NULL,'
    } elsif ($key eq 'SM' or $key eq 'CN') {
      $tables .= $key.' VARCHAR(138) NOT NULL,'
    } else {
      $tables .= $key . ' VARCHAR(50) NOT NULL,'
    }
  }
  $tables .= ' KEY(SM), KEY(ID)';

  return [
    # inheriting database and hive tables' creation
    @{ $self->SUPER::pipeline_create_commands },

    $self->db_cmd( 'CREATE TABLE ' . $self->o('summary_csv_table') . " ($tables)" ),

    $self->db_cmd( 'CREATE TABLE ' . $self->o('read_length_table') . ' (' .
        'fastq varchar(50) NOT NULL,' .
        'read_length int(50) NOT NULL,' .
        'PRIMARY KEY (fastq))' ),

    'mkdir -p ' . $self->o('rnaseq_dir'),
    'mkdir -p ' . $self->o('genome_dumps'),
    'mkdir -p ' . $self->o('rnasamba_output'),
    'mkdir -p ' . $self->o('cpc2_output'),
    'mkdir -p ' . $self->o('pcp_dir'),
  ];
}

sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{ $self->SUPER::pipeline_wide_parameters },
    genome_file          => $self->o('faidx_genome_file'),
    use_genome_flatfile => $self->o('use_genome_flatfile'),
    is_non_vert => $self->o('is_non_vert'),
  }
}

## See diagram for pipeline structure
sub pipeline_analyses {
  my ($self) = @_;

  my %commandline_params = (
    'ncbi'        => '-num_threads 3 -window_size 40',
    'wu'          => '-cpus 3 -hitdist 40',
    'legacy_ncbi' => '-a 3 -A 40',
  );

  return [

############################################################################
#
# RNA-seq analyses
#
############################################################################
    {
      -logic_name => 'download_rnaseq_csv',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadCsvENA',
      -rc_name => 'default',
      -parameters => {
        study_accession => $self->o('rnaseq_study_accession'),
        taxon_id => $self->o('taxon_id'),
        inputfile => $self->o('rnaseq_summary_file'),
        paired_end_only => $self->o('paired_end_only'),
        download_method => $self->o('download_method'),
      },
      -flow_into => {
        1 => ['download_genus_rnaseq_csv'],
      },
      -input_ids  => [{}],
    },

    {
      -logic_name => 'download_genus_rnaseq_csv',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadCsvENA',
      -rc_name => 'default',
      -parameters => {
        study_accession => $self->o('rnaseq_study_accession'),
        taxon_id => $self->o('genus_taxon_id'),
        inputfile => $self->o('rnaseq_summary_file_genus'),
        download_method => $self->o('download_method'),
      },
      -flow_into => {
        1 => ['checking_file_path'],
      },
    },

    {
      -logic_name => 'checking_file_path',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name    => 'default',
      -parameters => {
        cmd => 'EXIT_CODE=0; for F in ' . join( ' ',
          $self->o('uniprot_blast_exe_path')
          ) . '; do which "$F"; if [ "$?" == 1 ]; then EXIT_CODE=1;fi; done; '
          . 'if [ $EXIT_CODE -eq 1 ];then exit $EXIT_CODE;fi; '
          . 'for D in ' . join( ' ',
          $self->o('output_dir'),
          $self->o('input_dir'),
          $self->o('merge_dir'),
          ) . '; do mkdir -p "$D"; done; '
          . 'which lfs > /dev/null; if [ $? -eq 0 ]; then for D in ' . join( ' ',
          $self->o('output_dir'),
          $self->o('merge_dir'),
          $self->o('input_dir')
          ) . '; do lfs getdirstripe -q $D > /dev/null; if [ $? -eq 0 ]; then lfs setstripe -c -1 $D;fi;done;fi',
      },
      -flow_into => {
        '1->A' => [ 'create_fastq_download_jobs', 'index_rnaseq_genome_file' ],
        'A->1' => ['parse_summary_file'],
      },
    },

    {
      -logic_name => 'create_fastq_download_jobs',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        inputfile    => $self->o('rnaseq_summary_file'),
        column_names => $self->o('file_columns'),
        delimiter    => '\t',
      },
      -flow_into => {
        2 => ['download_RNASeq_fastqs'],
      },
    },

    {
      -logic_name => 'download_RNASeq_fastqs',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadData',
      -parameters => {
        output_dir => $self->o('input_dir'),
        download_method => $self->o('download_method'),
        uncompress => 0,
      },
      -flow_into => {
        2 => ['get_read_lengths'],
      },
      -analysis_capacity => 50,
      -max_retry_count => 12, #This is needed for big files > 10GB as there will be timeout and md5sum failures
    },

    {
      -logic_name => 'get_read_lengths',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCalculateReadLength',
      -parameters => {
        input_dir         => $self->o('input_dir'),
        read_length_table => $self->o('read_length_table'),
        _input_id_name => 'filename',
      },
    },

    {
      -logic_name => 'index_rnaseq_genome_file',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => '45GB_star',
      -parameters => {
        cmd => 'if [ ! -e "'.catfile($self->o('genome_dumps'), 'Genome').'" ]; then '.$self->o('star_path').' --runThreadN '.
          $self->o('star_threads').' --runMode genomeGenerate --outFileNamePrefix '.
          $self->o('genome_dumps').' --genomeDir '.
          $self->o('genome_dumps').' --genomeFastaFiles '.$self->o('faidx_genome_file').';fi',
      },
    },

    {
      -logic_name => 'parse_summary_file',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveParseCsvIntoTable',
      -rc_name    => 'default',
      -parameters => {
        column_names      => $self->o('file_columns'),
        sample_column     => $self->o('read_group_tag'),
        inputfile         => $self->o('rnaseq_summary_file'),
        delimiter         => $self->o('summary_file_delimiter'),
        csvfile_table     => $self->o('summary_csv_table'),
        pairing_regex     => $self->o('pairing_regex'),
        read_length_table => $self->o('read_length_table'),
      },
      -flow_into => {
        '2->A' => ['create_star_jobs'],
        'A->1' => ['scallopmerge'],
      },
    },

    {
      -logic_name => 'create_star_jobs',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateStarJobs',
      -parameters => {
        input_dir        => $self->o('input_dir'),
        sample_column    => $self->o('read_group_tag'),
        sample_id_column => $self->o('read_id_tag'),
        filename_column  => $self->o('filename_tag'),
        csvfile_table    => $self->o('summary_csv_table'),
        column_names     => $self->o('file_columns'),
      },
      -rc_name   => 'default',
      -flow_into => {
        2 => ['star'],
      },
    },

    {
      -logic_name => 'star',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveStar',
      -parameters => {
        disconnect_jobs    => 1,
        input_dir          => $self->o('input_dir'),
        output_dir         => $self->o('output_dir'),
        short_read_aligner => $self->o('star_path'),
        genome_dir         => catfile( $self->o('output_path'), 'genome_dumps' ),
        num_threads        => $self->o('star_threads'),
      },
      -flow_into => {
        2 => ['scallop'],
	ANYFAILURE => ['star_himem'],
      },
      -rc_name => '45GB_star',
    },

    {
      -logic_name => 'star_himem',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveStar',
      -parameters => {
        disconnect_jobs    => 1,
        input_dir          => $self->o('input_dir'),
        output_dir         => $self->o('output_dir'),
        short_read_aligner => $self->o('star_path'),
        genome_dir         => catfile( $self->o('output_path'), 'genome_dumps' ),
        num_threads        => $self->o('star_threads'),
      },
      -flow_into => {
        2 => ['scallop'],
      },
      -rc_name => '80GB_star',
    },

    {
      -logic_name => 'scallop',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::Scallop',
      -parameters => {
        output_dir => catdir( $self->o('output_dir'), 'scallop' ),
        scallop_path        => $self->o('scallop_path'),
        csv_summary_file       => $self->o('rnaseq_summary_file'),
        csv_summary_file_genus => $self->o('rnaseq_summary_file_genus'),
        num_threads            => $self->o('scallop_threads'),
      },
      -rc_name => '10GB_scallop',
      -flow_into => {
        MEMLIMIT => ['scallop_himem'],
      },
    },

    {
      -logic_name => 'scallop_himem',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::Scallop',
      -parameters => {
        output_dir => catdir( $self->o('output_dir'), 'scallop' ),
        scallop_path        => $self->o('scallop_path'),
        csv_summary_file       => $self->o('rnaseq_summary_file'),
        csv_summary_file_genus => $self->o('rnaseq_summary_file_genus'),
        num_threads            => $self->o('scallop_threads'),
      },
      -rc_name => '50GB_scallop',
      -flow_into => {
        MEMLIMIT => ['scallop_200GB'],
      },
    },

    {
      -logic_name => 'scallop_200GB',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::Scallop',
      -parameters => {
        output_dir => catdir( $self->o('output_dir'), 'scallop' ),
        scallop_path        => $self->o('scallop_path'),
        csv_summary_file       => $self->o('rnaseq_summary_file'),
        csv_summary_file_genus => $self->o('rnaseq_summary_file_genus'),
        num_threads            => $self->o('scallop_threads'),
      },
      -rc_name => '200GB_scallop',
    },

    {
      -logic_name => 'scallopmerge',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::Stringtie2Merge',
      -parameters => {
        input_gtf_dirs => [ catdir( $self->o('output_dir'), 'scallop' ) ],
        stringtie_merge_dir => catdir( $self->o('output_dir'), 'scallop', 'merge' ),
        stringtie2_path     => $self->o('stringtie2_path'),
        csv_summary_file    => $self->o('rnaseq_summary_file'),
        csv_summary_file_genus => $self->o('rnaseq_summary_file_genus'),
        num_threads            => $self->o('scallop_threads'),

      },
      -flow_into => {
        1 => ['create_scallop_initial_db'],
      },
      -rc_name => '10GB_scallop',
    },

    {
      -logic_name => 'create_scallop_initial_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('dna_db'),
        target_db   => $self->o('scallop_initial_db'),
        create_type => 'clone',
      },
      -max_retry_count => 0,
      -rc_name   => 'default',
      -flow_into => {
        1 => ['create_scallop_blast_db'],
      },
    },

    {
      -logic_name => 'create_scallop_blast_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('dna_db'),
        target_db   => $self->o('scallop_blast_db'),
        create_type => 'clone',
      },
      -max_retry_count => 0,
      -rc_name   => 'default',
      -flow_into => {
        '1->A' => ['generate_scallop_gtf_jobs'],
        'A->1' => ['star2introns'],
        1 => ['create_sample_jobs'],
      },
    },

    {
      -logic_name => 'create_sample_jobs',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        inputquery => 'SELECT DISTINCT('.$self->o('read_group_tag').'), "'.$self->o('rnaseq_data_provider').'" FROM '.$self->o('summary_csv_table'),
        column_names => ['sample_name', 'rnaseq_data_provider'],
      },
      -rc_name    => 'default',
      -priority => -2,
      -flow_into => {
        '2->A' => ['create_tissue_jobs'],
        'A->1' => ['merged_bam_file'],
      },
    },
    {
      -logic_name => 'create_tissue_jobs',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        inputquery => join(' ', 'SELECT', $self->o('read_group_tag'), ',', $self->o('read_id_tag'), 'FROM', $self->o('summary_csv_table'), 'WHERE', $self->o('read_group_tag'), '= "#sample_name#"'),
        column_names => [$self->o('read_group_tag'), $self->o('read_id_tag')],
      },
      -rc_name    => 'default',
      -priority => -2,
      -flow_into => {
        '2->A' => ['create_file_list'],
        'A->1' => ['merged_tissue_file'],
      },
    },
    {
      -logic_name => 'create_file_list',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        inputcmd => 'ls #input_dir#/#'.$self->o('read_id_tag').'#*.bam',
        column_names => ['filename'],
        input_dir => $self->o('output_dir'),
        fan_branch_code => 1,
      },
      -rc_name    => 'default',
      -priority => -2,
      -flow_into => {
          1 => ['?accu_name=filename&accu_address=[]&accu_input_variable=filename'],
      },
    },
    {
      -logic_name => 'merged_tissue_file',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMergeBamFiles',
      -parameters => {
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BamMergeStatic', $self->o('rnaseq_merge_type'))},
        # target_db is the database where we will write the files in the data_file table
        # You can use store_datafile => 0, if you don't want to store the output file
        target_db => $self->o('scallop_blast_db'),
        assembly_name => $self->o('assembly_name'),
        rnaseq_data_provider => $self->o('rnaseq_data_provider'),
        disconnect_jobs => 1,
        species => $self->o('species_name'),
        output_dir => $self->o('merge_dir'),
        input_dir => $self->o('output_dir'),
        samtools => $self->o('samtools_path'),
        picard_lib_jar => $self->o('picard_lib_jar'),
        use_threads => $self->o('rnaseq_merge_threads'),
        rename_file => 1,
      },
      -rc_name    => '3GB_rnaseq_multithread',
      -priority => -2,
      -flow_into => {
        1 => ['create_analyses_type_job', '?accu_name=filename&accu_address=[]&accu_input_variable=bam_file'],
      },
    },
    {
      -logic_name => 'create_analyses_type_job',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -rc_name    => 'default',
      -parameters => {
        inputlist => ['daf', 'ise'],
        column_names => ['type'],
        species => $self->o('species_name'),
      },
      -priority => -2,
      -flow_into => {
        2 => {'create_rnaseq_tissue_analyses' => {analyses => [{'-logic_name' => '#species#_#sample_name#_rnaseq_#type#'}]}},
      },
    },
    {
      -logic_name => 'create_rnaseq_tissue_analyses',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAddAnalyses',
      -rc_name    => 'default',
      -parameters => {
        source_type => 'list',
        target_db => $self->o('scallop_blast_db'),
      },
      -priority => -2,
    },
    {
      -logic_name => 'merged_bam_file',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMergeBamFiles',
      -parameters => {
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BamMergeStatic', $self->o('rnaseq_merge_type'))},
        # target_db is the database where we will write the files in the data_file table
        # You can use store_datafile => 0, if you don't want to store the output file
        target_db => $self->o('scallop_blast_db'),
        assembly_name => $self->o('assembly_name'),
        rnaseq_data_provider => $self->o('rnaseq_data_provider'),
        disconnect_jobs => 1,
        species => $self->o('species_name'),
        output_dir => $self->o('merge_dir'),
        input_dir => $self->o('merge_dir'),
        samtools => $self->o('samtools_path'),
        picard_lib_jar => $self->o('picard_lib_jar'),
        use_threads => $self->o('rnaseq_merge_threads'),
      },
      -rc_name    => '5GB_merge_multithread',
      -priority => -2,
      -flow_into => {
        1 => ['fan_merge_analyses'],
      },
    },

   {
      -logic_name => 'fan_merge_analyses',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'if [[ $(cut -f1 '.$self->o('rnaseq_summary_file')." | sort -u | wc -l) == 1 ]]; then exit 42; else exit 0;fi",
        return_codes_2_branches => {'42' => 2},
      },
      -rc_name    => 'default',
      -priority => -2,
      -flow_into  => {
  1 => ['create_merge_analyses_type_job'],
      },
    },

    {
      -logic_name => 'create_merge_analyses_type_job',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -rc_name    => 'default',
      -parameters => {
        inputlist => ['daf', 'ise'],
        column_names => ['type'],
        species => $self->o('species_name'),
      },
      -priority => -2,
      -flow_into => {
        2 => {'create_rnaseq_merge_analyses' => {analyses => [{'-logic_name' => '#species#_merged_rnaseq_#type#'}]}},
      },
    },
    {
      -logic_name => 'create_rnaseq_merge_analyses',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAddAnalyses',
      -rc_name    => 'default',
      -parameters => {
        source_type => 'list',
        target_db => $self->o('scallop_blast_db'),
      },
      -priority => -2,
    },
    {
      -logic_name => 'generate_scallop_gtf_jobs',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        inputcmd => 'ls #pattern#',
        pattern => catfile( $self->o('output_dir'), 'scallop', 'merge', '*.gtf'),
        column_names => ['filename'],
      },
      -flow_into => {
        '2->A' => ['load_scallop_transcripts'],
        'A->1' => ['create_logic_name_input_ids'],
      },
      -rc_name => '5GB',
    },

    {
      -logic_name => 'load_scallop_transcripts',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::LoadGTFBasic',
      -parameters => {
        target_db    => $self->o('scallop_initial_db'),
        loading_type => 'file',
        genome_file  => $self->o('faidx_genome_file'),
      },
      -rc_name   => '5GB',
    },

    {
      -logic_name => 'create_logic_name_input_ids',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -rc_name    => 'default',
      -parameters => {
        inputquery   => 'SELECT logic_name FROM analysis WHERE logic_name LIKE "#species_name#%rnaseq_gene"',
        column_names => ['logic_name'],
        db_conn      => $self->o('scallop_initial_db'),
        species_name => $self->o('species_name'),
      },
      -flow_into => {
        2 => ['create_slice_tissue_input_ids'],
      },
    },

    {
      -logic_name => 'create_slice_tissue_input_ids',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -rc_name    => 'default',
      -parameters => {
        iid_type              => 'slice',
        coord_system_name     => 'toplevel',
        slice                 => 1,
        include_non_reference => 0,
        top_level             => 1,
        feature_constraint    => 1,
        feature_type          => 'gene',
        target_db             => $self->o('scallop_initial_db'),
      },
      -flow_into => {
        2 => { 'create_gene_id_input_ids' => { iid => '#iid#', logic_name => '#logic_name#' } },
      },
    },

    {
      -logic_name => 'create_gene_id_input_ids',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -rc_name    => 'default',
      -parameters => {
        iid_type            => 'feature_id',
        coord_system_name   => 'toplevel',
        target_db           => $self->o('scallop_initial_db'),
        feature_type        => 'gene',
        batch_size          => 50,
        feature_logic_names => ['#logic_name#'],
      },
      -flow_into => {
        2 => { 'blast_scallop' => { iid => '#iid#', logic_name => '#logic_name#' } },
      },
    },

    {
      -logic_name => 'blast_scallop',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastRNASeqPep',
      -parameters => {
        source_db => $self->o('scallop_initial_db'),
        target_db => $self->o('scallop_blast_db'),
        dna_db    => $self->o('dna_db'),
        iid_type  => 'object_id',
        # path to index to fetch the sequence of the blast hit to calculate % coverage
        indicate_index => $self->o('protein_blast_index'),
        uniprot_index => [$self->o('protein_blast_db')],
        blast_program  => $self->o('uniprot_blast_exe_path'),
        %{ get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic', 'BlastGenscanPep', { BLAST_PARAMS => { -type => $self->o('blast_type') } } ) },
        commandline_params => $self->o('blast_type') eq 'wu' ? '-cpus=' . $self->o('use_threads') . ' -hitdist=40' : '-num_threads ' . $self->o('use_threads') . ' -window_size 40',
      },
      -flow_into => {
        '-1' => ['blast_scallop_longseq'],
        '2'  => ['blast_scallop_longseq'],
      },
      -rc_name => '3GB_multithread',
    },

    {
      -logic_name => 'blast_scallop_longseq',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastRNASeqPep',
      -parameters => {
        source_db => $self->o('scallop_initial_db'),
        target_db => $self->o('scallop_blast_db'),
        dna_db    => $self->o('dna_db'),
        iid_type  => 'object_id',
        # path to index to fetch the sequence of the blast hit to calculate % coverage
        indicate_index => $self->o('protein_blast_index'),
        uniprot_index  => [ $self->o('protein_blast_db') ],
        blast_program  => $self->o('uniprot_blast_exe_path'),
        %{ get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic', 'BlastGenscanPep', { BLAST_PARAMS => { -type => $self->o('blast_type') } } ) },
        commandline_params => $self->o('blast_type') eq 'wu' ? '-cpus=' . $self->o('use_threads') . ' -hitdist=40' : '-num_threads ' . $self->o('use_threads') . ' -window_size 40',
      },
      -rc_name => '10GB_multithread',
    },

    {
      -logic_name => 'star2introns',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveStar2Introns',
      -parameters => {
        star_junctions_dir => $self->o('output_dir'),
        intron_db          => $self->o('scallop_blast_db'),
        source_db          => $self->o('dna_db'),
        sample_column      => 'SM',
        sample_id_column   => 'ID',
        species            => $self->o('species_name'),
      },
      -rc_name   => '15GB',
      -flow_into => {
        '1' => ['copy_rnaseq_star_blast_db'],
      },
    },

    {
      -logic_name => 'copy_rnaseq_star_blast_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('scallop_blast_db'),
        target_db   => $self->o('rnaseq_for_layer_db'),
        create_type => 'copy',
      },
      -max_retry_count => 0,
      -rc_name   => 'default',
      -flow_into => {
        '1' => ['update_rnaseq_for_layer_biotypes_star'],
      },
    },

    {
      -logic_name => 'update_rnaseq_for_layer_biotypes_star',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('rnaseq_for_layer_db'),
        sql     => [
          'UPDATE gene SET biotype = "rnaseq_merged" WHERE source IN ("merged")',
          'UPDATE gene SET biotype = "rnaseq_tissue" WHERE biotype != "rnaseq_merged"',
          'UPDATE transcript JOIN gene USING(gene_id) SET transcript.biotype = gene.biotype',
        ],
      },
      -rc_name   => 'default',
      -flow_into => {
        '1' => ['classify_rnaseq_for_layer_models_star'],
      },
    },

    {
      -logic_name => 'classify_rnaseq_for_layer_models_star',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport',
      -parameters => {
        update_gene_biotype => 1,
        classification_type => 'standard',
        target_db           => $self->o('rnaseq_for_layer_db'),
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['rnaseq_for_layer_sanity_checks_star'],
      },

    },

    {
      -logic_name => 'rnaseq_for_layer_sanity_checks_star',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
      -parameters => {
        target_db                  => $self->o('rnaseq_for_layer_db'),
        sanity_check_type          => 'gene_db_checks',
        min_allowed_feature_counts => get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
          'gene_db_checks' )->{ $self->o('uniprot_set') }->{'rnaseq_blast'},
      },
      -flow_into => {
        1 => ['create_rnaseq_layer_nr_db_star'],
      },
    },

    {
      -logic_name => 'create_rnaseq_layer_nr_db_star',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('rnaseq_for_layer_db'),
        target_db   => $self->o('rnaseq_for_layer_nr_db'),
        create_type => 'copy',
      },
      -max_retry_count => 0,
      -rc_name   => 'default',
      -flow_into => {
        '1' => ['create_rnaseq_layer_nr_slices_star'],
      },
    },

    {
      -logic_name => 'create_rnaseq_layer_nr_slices_star',
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
        '2->A' => ['remove_redundant_rnaseq_layer_genes_star'],
        'A->1' => WHEN('#is_non_vert# eq "1"' => 'create_pcp_db', ELSE 'notification_pipeline_is_done',),
      },
    },

    {
      -logic_name => 'remove_redundant_rnaseq_layer_genes_star',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::RemoveRedundantGenes',
      -parameters => {
        target_db   => $self->o('rnaseq_for_layer_nr_db'),
        target_type => 'generic',
      },
      -rc_name => '5GB',
      -flow_into => {
      	 '2' => ['delete_short_reads'],
      },
    },

     { 
     -logic_name => 'delete_short_reads',
     -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
     -parameters => {
       cmd => 'rm '.$self->o('input_dir').'/*.gz',
     },
     -rc_name => 'default',
     },

    {
      -logic_name => 'create_pcp_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('scallop_initial_db'),
        target_db   => $self->o('pcp_db'),
        create_type => 'copy',
        force_drop  => 1,
      },
      -rc_name    => 'default',
      -flow_into  => {
        1 => ['dump_fasta'],
      },
    },

    {
      -logic_name => 'dump_fasta',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name    => '8GB',
      -parameters => {
        cmd => ' perl ' . $self->o('pcp_get_transcripts_script').
        ' -user ' . $self->o('user_r').
        ' -dna_user ' . $self->o('user_r').
        ' -dbname ' . $self->o('pcp_db_name').
        ' -dna_dbname ' . $self->o('dna_db_name').
        ' -host ' . $self->o('scallop_initial_db_host').
        ' -dna_host ' . $self->o('dna_db_host').
        ' -port ' . $self->o('scallop_initial_db_port').
        ' -dna_port ' . $self->o('dna_db_port').
        ' > ' . $self->o('cpc2_fasta_file')
      },
      -flow_into  => {
        '1->A' => ['run_cpc2','run_rnasamba'],
        'A->1' => ['impute_coding_genes']
      },
    },

    {
      -logic_name => 'run_cpc2',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name    => '8GB',
      -parameters => {
        cmd => 'singularity exec --bind '.$self->o('cpc2_output').' '.
        $self->o('cpc2').' python3 /CPC2_standalone-1.0.1/bin/CPC2.py'.
        ' -i '.$self->o('cpc2_fasta_file').' '.
        ' --ORF -o '.$self->o('cpc2_file').' '
      },
    },

    {
      -logic_name => 'run_rnasamba',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name    => '10GB',
      -parameters =>  {
        cmd => 'singularity exec --bind '.$self->o('rnasamba_output').' '.
        $self->o('rnasamba').' rnasamba classify '.
        $self->o('rnasamba_tsv_file').' '.
        $self->o('cpc2_fasta_file').' '.
        $self->o('rna_samba_weights')
      },
      -flow_into =>  {
        -1 => ['run_rnasamba_50GB'],
      },
    },

    {
      -logic_name => 'run_rnasamba_50GB',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name    => '50GB',
      -parameters =>  {
        cmd => 'singularity exec --bind '.$self->o('rnasamba_output').' '.
        $self->o('rnasamba').' rnasamba classify '.
        $self->o('rnasamba_tsv_file').' '.
        $self->o('cpc2_fasta_file').' '.
        $self->o('rna_samba_weights')
      },
      -flow_into =>  {
        -1 => ['run_rnasamba_100GB'],
      },
    },

    {
      -logic_name => 'run_rnasamba_100GB',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name    => '100GB',
      -parameters =>  {
        cmd => 'singularity exec --bind '.$self->o('rnasamba_output').' '.
        $self->o('rnasamba').' rnasamba classify '.
        $self->o('rnasamba_tsv_file').' '.
        $self->o('cpc2_fasta_file').' '.
        $self->o('rna_samba_weights')
      },
    },

    {
      -logic_name => 'impute_coding_genes',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . catfile($self->o('ensembl_analysis_scripts'), 'pcp', 'update_pcp_biotype.pl') .
        ' -user ' . $self->o('user').
        ' -pass ' . $self->o('password').
        ' -dbname  ' . $self->o('pcp_db_name').
        ' -port ' . $self->o('scallop_initial_db_port').
        ' -host ' . $self->o('scallop_initial_db_host').
        ' -cpc2 ' . $self->o('cpc2_txt_file').
        ' -rnas ' . $self->o('rnasamba_tsv_file'),
      },
      -rc_name    => '1GB',
      -flow_into => {
        '1'    => ['create_pcp_nr_db'],
      },
    },

    {
      -logic_name => 'create_pcp_nr_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('pcp_db'),
        target_db   => $self->o('pcp_nr_db'),
        create_type => 'copy',
        force_drop  => 1,
      },
      -rc_name    => 'default',
      -flow_into => {
        '1'    => ['remove_non_pcp_biotypes'],
      },
    },

    {
      -logic_name => 'remove_non_pcp_biotypes',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('pcp_nr_db'),
        sql => [
          'DELETE FROM gene WHERE biotype <> "pcp_protein_coding"'
        ],
      },
      -rc_name    => 'default',
      -flow_into => {
        '1' => ['create_pcp_nr_slices'],
      },
    },

    {
      -logic_name => 'create_pcp_nr_slices',
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
      -rc_name    => '2GB',
      -flow_into => {
        '2->A' => ['remove_redundant_pcp_genes'],
        'A->1' => ['notification_pipeline_is_done_pcp'],
      },
    },

    {
      -logic_name => 'remove_redundant_pcp_genes',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::RemoveRedundantGenes',
      -parameters => {
        dna_db => $self->o('dna_db'),
        target_db   => $self->o('pcp_nr_db'),
        target_type => 'generic',
      },
      -rc_name => '5GB',
    },

    {
      -logic_name => 'notification_pipeline_is_done',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::MessagePipeline',
      -parameters => {
        messages   => [
        {
          url => $self->o('homology_rnaseq_url'),
          logic_name => 'genblast_rnaseq_support',
          param => 'intron_db',
          data => $self->o('rnaseq_for_layer_nr_db'),
        },
        {
          url => $self->o('homology_rnaseq_url'),
          logic_name => 'genblast_rnaseq_support_himem',
          param => 'intron_db',
          data => $self->o('rnaseq_for_layer_nr_db'),
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'create_toplevel_slices',
          param => 'feature_dbs',
          data => $self->o('rnaseq_for_layer_nr_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'split_slices_on_intergenic',
          param => 'input_gene_dbs',
          data => $self->o('rnaseq_for_layer_nr_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'layer_annotation',
          param => 'SOURCEDB_REFS',
          data => $self->o('rnaseq_for_layer_nr_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'run_utr_addition',
          param => 'donor_dbs',
          data => $self->o('rnaseq_for_layer_nr_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'run_utr_addition_10GB',
          param => 'donor_dbs',
          data => $self->o('rnaseq_for_layer_nr_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'run_utr_addition_30GB',
          param => 'donor_dbs',
          data => $self->o('rnaseq_for_layer_nr_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'utr_memory_failover',
          param => 'donor_dbs',
          data => $self->o('rnaseq_for_layer_nr_db'),
          update => 1,
        },
        {
          url => $self->o('main_pipeline_url'),
          logic_name => 'initialise_rnaseq_db',
          param => 'rnaseq_blast_db',
          data => $self->o('scallop_blast_db'),
        },
        {
          url => $self->o('main_pipeline_url'),
          logic_name => 'initialise_rnaseq_db',
          param => 'rnaseq_refine_db',
          data => $self->o('scallop_blast_db'),
        }],
        tweak_script => catfile($self->o('enscode_root_dir'), 'ensembl-hive', 'scripts', 'tweak_pipeline.pl'),
      },
      -rc_name    => 'default',
    },

    {
      -logic_name => 'notification_pipeline_is_done_pcp',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::MessagePipeline',
      -parameters => {
        messages   => [
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'split_slices_on_intergenic',
          param => 'input_gene_dbs',
          data => $self->o('pcp_nr_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'layer_annotation',
          param => 'SOURCEDB_REFS',
          data => $self->o('pcp_nr_db'),
          update => 1,
        },
        ],
        tweak_script => catfile($self->o('enscode_root_dir'), 'ensembl-hive', 'scripts', 'tweak_pipeline.pl'),
      },
      -rc_name   => 'default',
      -flow_into => {
        '1' => ['notification_pipeline_is_done'],
      }
    }
  ];
}

sub resource_classes {
  my $self = shift;

  return {
    '1GB'     => { LSF => $self->lsf_resource_builder( 'production', 1000 ) },
    '2GB'     => { LSF => $self->lsf_resource_builder( 'production', 2000 ) },
    '5GB'     => { LSF => $self->lsf_resource_builder( 'production', 5000 ) },
    '8GB'     => { LSF => $self->lsf_resource_builder( 'production', 8000 ) },
    '10GB'     => { LSF => $self->lsf_resource_builder( 'production', 10000 ) },
    '50GB'     => { LSF => $self->lsf_resource_builder( 'production', 50000 ) },
    '100GB'     => { LSF => $self->lsf_resource_builder( 'production', 100000 ) },
    'default' => { LSF => $self->lsf_resource_builder( 'production', 900 ) },
    '3GB_multithread'     => { LSF => $self->lsf_resource_builder( 'production', 2900, undef, undef, $self->default_options->{use_threads} ) },
    '3GB_rnaseq_multithread'     => { LSF => $self->lsf_resource_builder( 'production', 2900, undef, undef, $self->default_options->{rnaseq_merge_threads} ) },
    '5GB_merge_multithread'     => { LSF => $self->lsf_resource_builder( 'production', 5000, undef, undef, $self->default_options->{rnaseq_merge_threads} ) },
    '10GB_multithread' => { LSF => $self->lsf_resource_builder( 'production', 10000, undef, undef, $self->default_options->{use_threads} ) },
    '45GB_star'    => { LSF => $self->lsf_resource_builder( 'production', 45000, undef, undef, ( $self->default_options->{'star_threads'} + 1 ) ) },
    '80GB_star'    => { LSF => $self->lsf_resource_builder( 'production', 80000, undef, undef, ( $self->default_options->{'star_threads'} + 1 ) ) },
    '10GB_scallop' => { LSF => $self->lsf_resource_builder( 'production', 10000, undef, undef, $self->default_options->{'scallop_threads'} ) },
    '50GB_scallop' => { LSF => $self->lsf_resource_builder( 'production', 50000, undef, undef, $self->default_options->{'scallop_threads'} ) },
    '200GB_scallop' => { LSF => $self->lsf_resource_builder( 'production', 200000, undef, undef, $self->default_options->{'scallop_threads'} ) },
    '15GB'     => { LSF => $self->lsf_resource_builder( 'production', 15000 ) },
    }
}

1;
