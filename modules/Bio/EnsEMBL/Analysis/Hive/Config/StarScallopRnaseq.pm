
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

package Bio::EnsEMBL::Analysis::Hive::Config::StarScallopRnaseq;

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
    'user_r'                    => '',                                                                                # read only db user
    'user'                      => '',                                                                                # write db user
    'password'                  => '',                                                                                # password for write db user
    'server_set'                => '',                                                                                # What server set to user, e.g. set1
    'pipe_db_server'            => '',                                                                                # host for pipe db
    'databases_server'          => '',                                                                                # host for general output dbs
    'dna_db_server'             => '',                                                                                # host for dna db
    'pipe_db_port'              => '',                                                                                # port for pipeline host
    'databases_port'            => '',                                                                                # port for general output db host
    'dna_db_port'               => '',                                                                                # port for dna db host
    'rnaseq_summary_file'       => '' || catfile( $self->o('rnaseq_dir'), $self->o('species_name') . '.csv' ),        # Set this if you have a pre-existing cvs file with the expected columns
    'rnaseq_summary_file_genus' => '' || catfile( $self->o('rnaseq_dir'), $self->o('species_name') . '_gen.csv' ),    # Set this if you have a pre-existing genus level cvs file with the expected columns
    'release_number' => '' || $self->o('ensembl_release'),
    'species_name'        => '',                                                                                      # e.g. mus_musculus
    'production_name'     => '',                                                                                      # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
    'taxon_id'            => '',                                                                                      # should be in the assembly report file
    'uniprot_set'         => '',                                                                                      # e.g. mammals_basic, check UniProtCladeDownloadStatic.pm module in hive config dir for suitable set,
    'output_path'         => '',                                                                                      # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses
    'assembly_name'       => '',                                                                                      # Name (as it appears in the assembly report file)
    'uniprot_version'     => 'uniprot_2019_04',                                                                       # What UniProt data dir to use for various analyses
    'paired_end_only'     => '1',                                                                                     # Will only use paired-end rnaseq data if 1

    # Keys for custom loading, only set/modify if that's what you're doing
    'protein_blast_db' => '' || catfile( $self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), 'PE12_vertebrata' ),    # Blast database for comparing the final models to.
    'protein_blast_index' => '' || catdir( $self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), 'PE12_vertebrata_index' ),    # Indicate Index for the blast database.
    'protein_entry_loc' => catfile( $self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), 'entry_loc' ),                       # Used by genscan blasts and optimise daf/paf. Don't change unless you know what you're doing

########################
# Pipe and ref db info
########################

    'pipe_db_name' => $self->o('dbowner') . '_' . $self->o('production_name') . '_pipe_' . $self->o('release_number'),
    'dna_db_name'  => $self->o('dbowner') . '_' . $self->o('production_name') . '_core_' . $self->o('release_number'),

    'reference_db_name'   => $self->o('dna_db_name'),
    'reference_db_server' => $self->o('dna_db_server'),
    'reference_db_port'   => $self->o('dna_db_port'),

    'star_rnaseq_for_layer_db_server' => $self->o('databases_server'),
    'star_rnaseq_for_layer_db_port'   => $self->o('databases_port'),

    'scallop_initial_db_server' => $self->o('databases_server'),
    'scallop_initial_db_port'   => $self->o('databases_port'),

    'scallop_blast_db_server' => $self->o('databases_server'),
    'scallop_blast_db_port'   => $self->o('databases_port'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'ensembl_release' => $ENV{ENSEMBL_RELEASE},    # this is the current release version on staging to be able to get the correct database

    databases_to_delete => ['star_rnaseq_for_layer_db','star_rnaseq_for_layer_nr_db','scallop_initial_db','scallop_blast_db'],

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

    'min_toplevel_slice_length' => 250,

########################
# Extra db settings
########################

    'num_tokens' => 10,

########################
# Executable paths
########################
    'star_path'       => '/homes/fergal/bin/STAR',
    'scallop_path'    => '/homes/fergal/bin/scallop',
    'stringtie2_path' => '/homes/fergal/bin/stringtie2',

    'blast_type'             => 'ncbi',                                                                         # It can be 'ncbi', 'wu', or 'legacy_ncbi'
    'uniprot_blast_exe_path' => catfile( $self->o('binary_base'), 'blastp' ),
    'bwa_path'               => catfile( $self->o('software_base_path'), 'opt', 'bwa-051mt', 'bin', 'bwa' ),    #You may need to specify the full path to the bwa binary

    'summary_file_delimiter' => '\t',            # Use this option to change the delimiter for your summary data file
    'summary_csv_table'      => 'csv_data',
    'read_length_table'      => 'read_length',
    'rnaseq_data_provider'   => 'ENA',           #It will be set during the pipeline or it will use this value

    'rnaseq_dir' => catdir( $self->o('output_path'), 'rnaseq' ),
    'input_dir'  => catdir( $self->o('rnaseq_dir'),  'input' ),
    'output_dir' => catdir( $self->o('rnaseq_dir'),  'output' ),

    'rnaseq_ftp_base' => 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/',

    'use_ucsc_naming' => 0,

    # If your reads are unpaired you may want to run on slices to avoid
    # making overlong rough models.  If you want to do this, specify a
    # slice length here otherwise it will default to whole chromosomes.
    slice_length => 10000000,

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
    file_columns => [ 'SM', 'ID', 'is_paired', 'filename', 'is_mate_1', 'read_length', 'is_13plus', 'CN', 'PL', 'DS' ],

    'filename_tag' => 'filename',    # For the analysis that creates star jobs, though I assume we should need to do it this way

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# No option below this mark should be modified
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#

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

    'star_rnaseq_for_layer_db' => {
      -dbname => $self->o('dbowner') . '_' . $self->o('production_name') . 'star_rs_layer_' . $self->o('release_number'),
      -host   => $self->o('star_rnaseq_for_layer_db_server'),
      -port   => $self->o('star_rnaseq_for_layer_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'star_rnaseq_for_layer_nr_db' => {
      -dbname => $self->o('dbowner') . '_' . $self->o('production_name') . '_star_rs_layer_nr_' . $self->o('release_number'),
      -host   => $self->o('star_rnaseq_for_layer_db_server'),
      -port   => $self->o('star_rnaseq_for_layer_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'scallop_initial_db' => {
      -dbname => $self->o('dbowner') . '_' . $self->o('production_name') . '_scallop_initial_' . $self->o('release_number'),
      -host   => $self->o('scallop_initial_db_server'),
      -port   => $self->o('scallop_initial_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'scallop_blast_db' => {
      -dbname => $self->o('dbowner') . '_' . $self->o('production_name') . '_scallop_blast_' . $self->o('release_number'),
      -host   => $self->o('scallop_blast_db_server'),
      -port   => $self->o('scallop_blast_db_port'),
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
    elsif ( $key eq 'DS' ) {
      $tables .= $key . ' VARCHAR(255) NOT NULL,'
    }
    else {
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
  ];
}

sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{ $self->SUPER::pipeline_wide_parameters },
    wide_ensembl_release => $self->o('ensembl_release'),
    genome_file          => $self->o('faidx_genome_file'),
    }
}

=head2 create_header_line

 Arg [1]    : Arrayref String, it will contains the values of 'file_columns'
 Example    : create_header_line($self->o('file_columns');
 Description: It will create a RG line using only the keys present in your csv file
 Returntype : String representing the RG line in a BAM file
 Exceptions : None


=cut

sub create_header_line {
  my ($items) = shift;

  my @read_tags = qw(ID SM DS CN DT FO KS LB PG PI PL PM PU);
  my $read_line = '@RG';
  foreach my $rt (@read_tags) {
    $read_line .= "\t$rt:#$rt#" if ( grep( $rt eq $_, @$items ) );
  }
  return $read_line . "\n";
}

## See diagram for pipeline structure
sub pipeline_analyses {
  my ($self) = @_;

  my %commandline_params = (
    'ncbi'        => '-num_threads 3 -window_size 40',
    'wu'          => '-cpus 3 -hitdist 40',
    'legacy_ncbi' => '-a 3 -A 40',
  );
  my $header_line = create_header_line( $self->default_options->{'file_columns'} );

  return [

############################################################################
#
# RNA-seq analyses
#
############################################################################
    {
      -logic_name => 'create_star_rnaseq_for_layer_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('dna_db'),
        target_db   => $self->o('star_rnaseq_for_layer_db'),
        create_type => 'clone',
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['fan_star_rnaseq_for_layer_db'],
      },
    },

    {
      -logic_name => 'fan_star_rnaseq_for_layer_db',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'if [ "#skip_rnaseq#" -ne "0" ]; then exit 42; else exit 0;fi',
        return_codes_2_branches => { '42' => 2 },
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['checking_file_path'],
      },
    },

    {
      -logic_name => 'checking_file_path',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name    => '1GB',
      -parameters => {
        cmd => 'EXIT_CODE=0; for F in ' . join( ' ',
          $self->o('bwa_path'),
          $self->o('uniprot_blast_exe_path')
          ) . '; do which "$F"; if [ "$?" == 1 ]; then EXIT_CODE=1;fi; done; '
          . 'if [ $EXIT_CODE -eq 1 ];then exit $EXIT_CODE;fi; '
          . 'for D in ' . join( ' ',
          $self->o('output_dir'),
          $self->o('input_dir'),
          ) . '; do mkdir -p "$D"; done; '
          . 'which lfs > /dev/null; if [ $? -eq 0 ]; then for D in ' . join( ' ',
          $self->o('output_dir'),
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
        2 => { 'download_RNASeq_fastqs' => { 'iid' => '#filename#' } },
      },
    },

    {
      -logic_name => 'download_RNASeq_fastqs',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadRNASeqFastqs',
      -parameters => {
        ftp_base_url => $self->o('rnaseq_ftp_base'),
        input_dir    => $self->o('input_dir'),
      },
      -flow_into => {
        1 => ['get_read_lengths'],
      },
      -analysis_capacity => 50,
    },

    {
      -logic_name => 'get_read_lengths',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCalculateReadLength',
      -parameters => {
        input_dir         => $self->o('input_dir'),
        read_length_table => $self->o('read_length_table'),
      },
      -flow_into => {
        1 => ['split_fastq_files'],
      },
    },

    {
      -logic_name => 'index_rnaseq_genome_file',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name    => '5GB',
      -parameters => {
        cmd => 'if [ ! -e "' . $self->o('faidx_genome_file') . '.ann" ]; then ' . $self->o('bwa_path') . ' index -a bwtsw ' . $self->o('faidx_genome_file') . ';fi',
      },
    },

    {
      -logic_name => 'parse_summary_file',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveParseCsvIntoTable',
      -rc_name    => '1GB',
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
        use_threading    => $self->o('use_threads'),
      },
      -rc_name   => '1GB',
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
      },
      -rc_name => '45GB_star',
    },

    {
      -logic_name => 'scallop',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::Scallop',
      -parameters => {
        output_dir => catdir( $self->o('output_dir'), 'scallop' ),
        stringtie2_path        => $self->o('scallop_path'),
        csv_summary_file       => $self->o('rnaseq_summary_file'),
        csv_summary_file_genus => $self->o('rnaseq_summary_file_genus'),
        num_threads            => $self->o('scallop_threads'),
      },
      -rc_name => '10GB_scallop',
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
      -rc_name   => 'default',
      -flow_into => {
        '1->A' => ['generate_scallop_gtf_jobs'],
        'A->1' => ['star2introns'],
      },
    },

    {
      -logic_name => 'generate_scallop_gtf_jobs',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::GenerateGTFLoadingJobs',
      -parameters => {
        gtf_dir => catdir( $self->o('output_dir'), 'scallop', 'merge' ),
      },
      -flow_into => {
        2 => ['load_scallop_transcripts'],
      },
      -rc_name => '5GB',
    },

    {
      -logic_name => 'load_scallop_transcripts',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::LoadGTFBasic',
      -parameters => {
        target_db    => $self->o('scallop_initial_db'),
        loading_type => 'range',
        genome_file  => $self->o('faidx_genome_file'),
        logic_name   => 'scallop',
        module       => 'Stringtie2',
      },
      -rc_name   => '5GB',
      -flow_into => {
        1 => 'create_slice_tissue_input_ids',
      },
    },

    {
      -logic_name => 'create_slice_tissue_input_ids',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -rc_name    => '1GB',
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
      -rc_name    => '1GB',
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
        uniprot_index  => [ $self->o('protein_blast_db') ],
        blast_program  => $self->o('uniprot_blast_exe_path'),
        %{ get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic', 'BlastGenscanPep', { BLAST_PARAMS => { -type => $self->o('blast_type') } } ) },
        commandline_params => $self->o('blast_type') eq 'wu' ? '-cpus=' . $self->o('use_threads') . ' -hitdist=40' : '-num_threads ' . $self->o('use_threads') . ' -window_size 40',
      },
      -flow_into => {
        '-1' => ['blast_scallop_longseq'],
        '2'  => ['blast_scallop_longseq'],
      },
      -rc_name => 'blast',
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
      -rc_name => 'blast10GB',
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
      -rc_name   => 'default',
      -flow_into => {
        '1' => ['copy_rnaseq_star_blast_db'],
      },
    },

    {
      -logic_name => 'copy_rnaseq_star_blast_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('scallop_blast_db'),
        target_db   => $self->o('star_rnaseq_for_layer_db'),
        create_type => 'copy',
        force_drop  => 1,
      },
      -rc_name   => 'default',
      -flow_into => {
        '1' => ['update_rnaseq_for_layer_biotypes_star'],
      },
    },

    {
      -logic_name => 'update_rnaseq_for_layer_biotypes_star',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('star_rnaseq_for_layer_db'),
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
        target_db           => $self->o('star_rnaseq_for_layer_db'),
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
        target_db                  => $self->o('star_rnaseq_for_layer_db'),
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
        source_db   => $self->o('star_rnaseq_for_layer_db'),
        target_db   => $self->o('star_rnaseq_for_layer_nr_db'),
        create_type => 'copy',
        force_drop  => 1,
      },
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
        '2' => ['remove_redundant_rnaseq_layer_genes_star'],
      },
    },

    {
      -logic_name => 'remove_redundant_rnaseq_layer_genes_star',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::RemoveRedundantGenes',
      -parameters => {
        target_db   => $self->o('star_rnaseq_for_layer_nr_db'),
        target_type => 'generic',
      },
      -rc_name => '5GB',
    },

  ];
}

sub resource_classes {
  my $self = shift;

  return {
    '1GB'     => { LSF => $self->lsf_resource_builder( 'production-rh74', 1000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    '2GB'     => { LSF => $self->lsf_resource_builder( 'production-rh74', 2000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    '5GB'     => { LSF => $self->lsf_resource_builder( 'production-rh74', 5000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    'default' => { LSF => $self->lsf_resource_builder( 'production-rh74', 900,  [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    'blast'     => { LSF => $self->lsf_resource_builder( 'production-rh74', 2900,  [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], undef, 3 ) },
    'blast10GB' => { LSF => $self->lsf_resource_builder( 'production-rh74', 10000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], undef, 3 ) },
    '45GB_star'    => { LSF => $self->lsf_resource_builder( 'production-rh74', 45000, undef, undef, ( $self->default_options->{'star_threads'} + 1 ) ) },
    '10GB_scallop' => { LSF => $self->lsf_resource_builder( 'production-rh74', 10000, undef, undef, ( $self->default_options->{'scallop_threads'} ) ) },
    }
}

sub hive_capacity_classes {
  my $self = shift;

  return {
    'hc_very_low' => 35,
    'hc_low'      => 200,
    'hc_medium'   => 500,
    'hc_high'     => 1000,
  };
}

sub check_file_in_ensembl {
  my ( $self, $file_path ) = @_;
  push @{ $self->{'_ensembl_file_paths'} }, $file_path;
  return $self->o('enscode_root_dir') . '/' . $file_path;
}

1;
