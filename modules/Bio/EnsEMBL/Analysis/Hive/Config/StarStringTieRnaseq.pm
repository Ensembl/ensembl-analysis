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

package Bio::EnsEMBL::Analysis::Hive::Config::StarStringTieRnaseq;

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
    'dbowner'                   => '' || $ENV{EHIVE_USER} || $ENV{USER},
    'pipeline_name'             => '' || $self->o('production_name').'_'.$self->o('ensembl_release'),
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
    'registry_host'             => '', # host for registry db
    'registry_port'             => '', # port for registry db
    'registry_db'               => '', # name for registry db
    'rnaseq_summary_file'       => '' || catfile($self->o('rnaseq_dir'), $self->o('species_name').'.csv'), # Set this if you have a pre-existing cvs file with the expected columns
    'star_rnaseq_summary_file'  => '' || catfile($self->o('rnaseq_dir'), 'star_'.$self->o('species_name').'.csv'),
    'rnaseq_summary_file_genus' => '' || catfile($self->o('rnaseq_dir'), $self->o('species_name').'_gen.csv'), # Set this if you have a pre-existing genus level cvs file with the expected columns
    'release_number'            => '' || $self->o('ensembl_release'),
    'species_name'              => '', # e.g. mus_musculus
    'production_name'           => '', # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
    'taxon_id'                  => '', # should be in the assembly report file
    'uniprot_set'               => '', # e.g. mammals_basic, check UniProtCladeDownloadStatic.pm module in hive config dir for suitable set,
    'output_path'               => '', # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses
    'assembly_name'             => '', # Name (as it appears in the assembly report file)
    'registry_file'             => '' || catfile($self->o('output_path'), "Databases.pm"), # Path to databse registry for LastaZ and Production sync
    'stable_id_prefix'          => '', # e.g. ENSPTR. When running a new annotation look up prefix in the assembly registry db
    'use_genome_flatfile'       => '1',# This will read sequence where possible from a dumped flatfile instead of the core db
    'skip_rnaseq'               => '0', # Will skip rnaseq analyses if 1
    'uniprot_version'           => 'uniprot_2019_04', # What UniProt data dir to use for various analyses
    'paired_end_only'           => '1', # Will only use paired-end rnaseq data if 1

    # Keys for custom loading, only set/modify if that's what you're doing
    'protein_blast_db'             => '' || catfile($self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), 'PE12_vertebrata'), # Blast database for comparing the final models to.
    'protein_blast_index'          => '' || catdir($self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), 'PE12_vertebrata_index'), # Indicate Index for the blast database.
    'protein_entry_loc'            => catfile($self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), 'entry_loc'), # Used by genscan blasts and optimise daf/paf. Don't change unless you know what you're doing

########################
# Pipe and ref db info
########################

    'pipe_db_name'                  => $self->o('dbowner').'_'.$self->o('production_name').'_pipe_'.$self->o('release_number'),
    'dna_db_name'                   => $self->o('dbowner').'_'.$self->o('production_name').'_core_'.$self->o('release_number'),

    'reference_db_name'            => $self->o('dna_db_name'),
    'reference_db_server'          => $self->o('dna_db_server'),
    'reference_db_port'            => $self->o('dna_db_port'),

    'rnaseq_for_layer_db_server'   => $self->o('databases_server'),
    'rnaseq_for_layer_db_port'     => $self->o('databases_port'),

    'star_rnaseq_for_layer_db_server'   => $self->o('databases_server'),
    'star_rnaseq_for_layer_db_port'     => $self->o('databases_port'),

    'rnaseq_db_server'             => $self->o('databases_server'),
    'rnaseq_db_port'               => $self->o('databases_port'),

    'rnaseq_rough_db_server'       => $self->o('databases_server'),
    'rnaseq_rough_db_port'         => $self->o('databases_port'),

    'rnaseq_refine_db_server'       => $self->o('databases_server'),
    'rnaseq_refine_db_port'         => $self->o('databases_port'),

    'rnaseq_blast_db_server'       => $self->o('databases_server'),
    'rnaseq_blast_db_port'         => $self->o('databases_port'),

    'stringtie_initial_db_server'  => $self->o('databases_server'),
    'stringtie_initial_db_port'    => $self->o('databases_port'),

    'stringtie_blast_db_server'    => $self->o('databases_server'),
    'stringtie_blast_db_port'      => $self->o('databases_port'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'ensembl_release'              => $ENV{ENSEMBL_RELEASE}, # this is the current release version on staging to be able to get the correct database

    databases_to_delete => ['reference_db', 'genewise_db', 'projection_db', 'selected_projection_db', 'layering_db', 'utr_db', 'genebuilder_db', 'pseudogene_db', 'ncrna_db', 'final_geneset_db', 'refseq_db', 'cdna2genome_db', 'rnaseq_blast_db', 'rnaseq_refine_db', 'rnaseq_rough_db', 'otherfeatures_db', 'rnaseq_db'],#, 'projection_realign_db'

########################
# BLAST db paths
########################
    'base_blast_db_path'        => $ENV{BLASTDB_DIR},

######################################################
#
# Mostly constant settings
#
######################################################

    genome_dumps                  => catdir($self->o('output_path'), 'genome_dumps'),
    # This one is used by most analyses that run against a genome flatfile like exonerate, genblast etc. Has slice name style headers. Is softmasked
    softmasked_genome_file        => catfile($self->o('genome_dumps'), $self->o('species_name').'_softmasked_toplevel.fa'),
    # This one is used in replacement of the dna table in the core db, so where analyses override slice->seq. Has simple headers with just the seq_region name. Also used by bwa in the RNA-seq analyses. Not masked
    faidx_genome_file             => catfile($self->o('genome_dumps'), $self->o('species_name').'_toplevel.fa'),
    # This one is a cross between the two above, it has the seq_region name header but is softmasked. It is used by things that would both want to skip using the dna table and also want to avoid the repeat_feature table, e.g. bam2introns
    faidx_softmasked_genome_file  => catfile($self->o('genome_dumps'), $self->o('species_name').'_softmasked_toplevel.fa.reheader'),

    'min_toplevel_slice_length'   => 250,

########################
# Extra db settings
########################

    'num_tokens' => 10,

########################
# Executable paths
########################
    'star_path'              => '/homes/fergal/bin/STAR',
    'stringtie2_path'        => '/homes/fergal/bin/stringtie',

    'blast_type' => 'ncbi', # It can be 'ncbi', 'wu', or 'legacy_ncbi'
    'uniprot_blast_exe_path' => catfile($self->o('binary_base'), 'blastp'),
    'exonerate_path'         => catfile($self->o('software_base_path'), 'opt', 'exonerate09', 'bin', 'exonerate'),
    samtools_path => catfile($self->o('binary_base'), 'samtools'), #You may need to specify the full path to the samtools binary
    picard_lib_jar => catfile($self->o('software_base_path'), 'Cellar', 'picard-tools', '2.6.0', 'libexec', 'picard.jar'), #You need to specify the full path to the picard library
    bwa_path => catfile($self->o('software_base_path'), 'opt', 'bwa-051mt', 'bin', 'bwa'), #You may need to specify the full path to the bwa binary
    refine_ccode_exe => catfile($self->o('binary_base'), 'RefineSolexaGenes'), #You may need to specify the full path to the RefineSolexaGenes binary

# RNA-seq pipeline stuff
    # You have the choice between:
    #  * using a csv file you already created
    #  * using a study_accession like PRJEB19386
    #  * using the taxon_id of your species
    # 'rnaseq_summary_file' should always be set. If 'taxon_id' or 'study_accession' are not undef
    # they will be used to retrieve the information from ENA and to create the csv file. In this case,
    # 'file_columns' and 'summary_file_delimiter' should not be changed unless you know what you are doing
    'study_accession'     => '',
    'max_reads_per_split' => 2500000, # This sets the number of reads to split the fastq files on
    'max_total_reads'     => 200000000, # This is the total number of reads to allow from a single, unsplit file

    'summary_file_delimiter' => '\t', # Use this option to change the delimiter for your summary data file
    'summary_csv_table' => 'csv_data',
    'read_length_table' => 'read_length',
    'rnaseq_data_provider' => 'ENA', #It will be set during the pipeline or it will use this value

    'rnaseq_dir'    => catdir($self->o('output_path'), 'rnaseq'),
    'input_dir'     => catdir($self->o('rnaseq_dir'),'input'),
    'output_dir'    => catdir($self->o('rnaseq_dir'),'output'),
    'merge_dir'     => catdir($self->o('rnaseq_dir'),'merge'),
    'sam_dir'       => catdir($self->o('rnaseq_dir'),'sams'),
    'header_file'   => catfile($self->o('output_dir'), '#'.$self->o('read_id_tag').'#_header.h'),

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
    
    # Regular expressions for splitting the fastq files
    split_paired_regex   => '(\S+)(\_\d\.\S+)',
    split_single_regex  => '([^.]+)(\.\S+)',

    # Do you want to make models for the each individual sample as well
    # as for the pooled samples (1/0)?
    single_tissue => 1,

    # What Read group tag would you like to group your samples
    # by? Default = ID
    read_group_tag => 'SM',
    read_id_tag => 'ID',

    use_threads => 3,
    rnaseq_merge_threads => 12,
    rnaseq_merge_type => 'samtools',
    read_min_paired => 50,
    read_min_mapped => 50,
    other_isoforms => 'other', # If you don't want isoforms, set this to undef
    maxintron => 200000,

    star_threads => 12,
    stringtie_threads => 2,

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

    'filename_tag'   => 'filename', # For the analysis that creates star jobs, though I assume we should need to do it this way

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# No option below this mark should be modified
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#

##################################
# Memory settings for the analyses
##################################
    'default_mem'          => '900',

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

    'rnaseq_for_layer_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_rs_layer_'.$self->o('release_number'),
      -host   => $self->o('rnaseq_for_layer_db_server'),
      -port   => $self->o('rnaseq_for_layer_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'star_rnaseq_for_layer_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'star_rs_layer_'.$self->o('release_number'),
      -host   => $self->o('rnaseq_for_layer_db_server'),
      -port   => $self->o('rnaseq_for_layer_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'rnaseq_for_layer_nr_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_rs_layer_nr_'.$self->o('release_number'),
      -host   => $self->o('rnaseq_for_layer_db_server'),
      -port   => $self->o('rnaseq_for_layer_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'star_rnaseq_for_layer_nr_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_star_rs_layer_nr_'.$self->o('release_number'),
      -host   => $self->o('rnaseq_for_layer_db_server'),
      -port   => $self->o('rnaseq_for_layer_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

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

    'rnaseq_rough_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_rough_'.$self->o('release_number'),
      -host   => $self->o('rnaseq_rough_db_server'),
      -port   => $self->o('rnaseq_rough_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'stringtie_initial_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_stringtie_initial_'.$self->o('release_number'),
      -host   => $self->o('stringtie_initial_db_server'),
      -port   => $self->o('stringtie_initial_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'stringtie_blast_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_stringtie_blast_'.$self->o('release_number'),
      -host   => $self->o('stringtie_blast_db_server'),
      -port   => $self->o('stringtie_blast_db_port'),
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
        paired => 1,
        read_length => 1,
        is_13plus => 1,
        is_mate_1 => 1,
        );
    # We need to store the values of the csv file to easily process it. It will be used at different stages
    foreach my $key (@{$self->default_options->{'file_columns'}}) {
        if (exists $small_columns{$key}) {
            $tables .= $key.' SMALLINT UNSIGNED NOT NULL,'
        }
        elsif ($key eq 'DS') {
            $tables .= $key.' VARCHAR(255) NOT NULL,'
        }
        else {
            $tables .= $key.' VARCHAR(50) NOT NULL,'
        }
    }
    $tables .= ' KEY(SM), KEY(ID)';

    return [
    # inheriting database and hive tables' creation
      @{$self->SUPER::pipeline_create_commands},

      $self->db_cmd('CREATE TABLE '.$self->o('summary_csv_table')." ($tables)"),

      $self->db_cmd('CREATE TABLE '.$self->o('read_length_table').' ('.
                    'fastq varchar(50) NOT NULL,'.
                    'read_length int(50) NOT NULL,'.
                    'PRIMARY KEY (fastq))'),

      'mkdir -p '.$self->o('rnaseq_dir'),
      'mkdir -p '.$self->o('genome_dumps'),
    ];
}


sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
    skip_rnaseq => $self->o('skip_rnaseq'),
    load_toplevel_only => $self->o('load_toplevel_only'),
    wide_ensembl_release => $self->o('ensembl_release'),
    use_genome_flatfile  => $self->o('use_genome_flatfile'),
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
        $read_line .= "\t$rt:#$rt#" if (grep($rt eq $_, @$items));
    }
    return $read_line."\n";
}

## See diagram for pipeline structure
sub pipeline_analyses {
    my ($self) = @_;

    my %commandline_params = (
      'ncbi' => '-num_threads 3 -window_size 40',
      'wu' => '-cpus 3 -hitdist 40',
      'legacy_ncbi' => '-a 3 -A 40',
      );
    my $header_line = create_header_line($self->default_options->{'file_columns'});

    return [

############################################################################
#
# RNA-seq analyses
#
############################################################################
      {
        -logic_name => 'create_rnaseq_for_layer_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
          source_db => $self->o('dna_db'),
          target_db => $self->o('rnaseq_for_layer_db'),
          create_type => 'clone',
        },
        -rc_name    => 'default',
        -flow_into => {
          '1->A' => ['fan_rnaseq_for_layer_db'],
          'A->1' => ['create_long_read_final_db'],
        },
      },


      {
        -logic_name => 'fan_rnaseq_for_layer_db',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'if [ "#skip_rnaseq#" -ne "0" ]; then exit 42; else exit 0;fi',
                         return_codes_2_branches => {'42' => 2},
                       },
        -rc_name => 'default',
        -flow_into  => {
          1 => ['checking_file_path'],
        },
      },


      {
        -logic_name => 'checking_file_path',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -rc_name => '1GB',
        -parameters => {
          cmd => 'EXIT_CODE=0; for F in '.join (' ',
              $self->o('bwa_path'),
              $self->o('samtools_path'),
              $self->o('exonerate_path'),
              $self->o('uniprot_blast_exe_path')
              ).'; do which "$F"; if [ "$?" == 1 ]; then EXIT_CODE=1;fi; done; '
            .'if [ $EXIT_CODE -eq 1 ];then exit $EXIT_CODE;fi; '
            .'for D in '.join(' ',
              $self->o('output_dir'),
              $self->o('input_dir'),
              $self->o('merge_dir'),
              $self->o('sam_dir')
              ).'; do mkdir -p "$D"; done; '
            .'which lfs > /dev/null; if [ $? -eq 0 ]; then for D in '.join(' ',
              $self->o('output_dir'),
              $self->o('input_dir'),
              $self->o('merge_dir')
              ).'; do lfs getdirstripe -q $D > /dev/null; if [ $? -eq 0 ]; then lfs setstripe -c -1 $D;fi;done;fi',
        },
        -flow_into => {
          '1->A' => ['create_fastq_download_jobs','index_rnaseq_genome_file'],
          'A->1' => ['create_rough_db'],
        },
      },


      {
        -logic_name => 'create_fastq_download_jobs',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
        -parameters => {
          inputfile => $self->o('rnaseq_summary_file'),
          column_names => $self->o('file_columns'),
          delimiter => '\t',
        },
        -flow_into => {
          2 => {'download_RNASeq_fastqs' => {'iid' => '#filename#'}},
        },
      },


      {
        -logic_name => 'download_RNASeq_fastqs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadRNASeqFastqs',
        -parameters =>{
          ftp_base_url => $self->o('rnaseq_ftp_base'),
          input_dir => $self->o('input_dir'),
        },
        -flow_into => {
          1 => ['get_read_lengths'],
        },
        -analysis_capacity => 50,
      },


      {
        -logic_name => 'get_read_lengths',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCalculateReadLength',
        -parameters =>{
          input_dir => $self->o('input_dir'),
          read_length_table => $self->o('read_length_table'),
        },
        -flow_into => {
          1 => ['split_fastq_files'],
        },
      },


      {
        -logic_name => 'split_fastq_files',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::SplitFastQFiles',
        -parameters => {
          'max_reads_per_split' => $self->o('max_reads_per_split'),
          'max_total_reads'     => $self->o('max_total_reads'),
          'rnaseq_summary_file' => $self->o('rnaseq_summary_file'),
          'fastq_dir'           => $self->o('input_dir'),
	  'pairing_regex'       => $self->o('split_paired_regex'),
	  'single_regex'        => $self->o('split_single_regex'),
        },
      },


     {
        -logic_name => 'index_rnaseq_genome_file',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -rc_name => '5GB',
        -parameters => {
          cmd => 'if [ ! -e "'.$self->o('faidx_genome_file').'.ann" ]; then '.$self->o('bwa_path').' index -a bwtsw '.$self->o('faidx_genome_file').';fi',
        },
      },


      {
        -logic_name => 'create_rough_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
          source_db => $self->o('dna_db'),
          target_db => $self->o('rnaseq_rough_db'),
          create_type => 'clone',
        },
        -rc_name => '1GB',
        -flow_into => {
          1 => ['backup_original_csv'],
        },
      },


      {
        -logic_name => 'backup_original_csv',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'cp '.$self->o('rnaseq_summary_file').' '.$self->o('rnaseq_summary_file').'_orig_bak',
        },
        -flow_into => {
          1 => ['create_updated_csv'],
        },
      },


      {
        -logic_name => 'create_updated_csv',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'cat '.catfile($self->o('input_dir'), '*_new.csv').' > '.$self->o('rnaseq_summary_file'),
        },
        -flow_into => {
          1 => ['parse_summary_file_bwa', 'parse_summary_file_star'],
        },
      },

      {
        -logic_name => 'parse_summary_file_bwa',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveParseCsvIntoTable',
        -rc_name => '1GB',
        -parameters => {
          column_names => $self->o('file_columns'),
          sample_column => $self->o('read_group_tag'),
          inputfile => $self->o('rnaseq_summary_file'),
          delimiter => $self->o('summary_file_delimiter'),
          csvfile_table => $self->o('summary_csv_table'),
          pairing_regex => $self->o('pairing_regex'),
          read_length_table => $self->o('read_length_table'),
        },
        -flow_into => {
          '2->A' => [ 'create_tissue_jobs'],
          'A->1' => [ 'merged_bam_file' ],
        },
      },
      {
        -logic_name => 'create_tissue_jobs',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
        -parameters => {
          inputquery => join(' ', 'SELECT', $self->o('read_group_tag'), ',', $self->o('read_id_tag'), ', is_paired, CN', 'FROM', $self->o('summary_csv_table'), 'WHERE', $self->o('read_group_tag'), '= "#sample_name#"'),
          column_names => [$self->o('read_group_tag'), $self->o('read_id_tag'), 'is_paired', 'rnaseq_data_provider'],
        },
        -rc_name    => '1GB',
        -flow_into => {
          '2->A' => ['create_bwa_jobs'],
          'A->1' => ['merged_tissue_file'],
        },
      },
      {
        -logic_name => 'create_bwa_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateBWAJobs',
        -parameters => {
          sample_column => $self->o('read_group_tag'),
          sample_id_column => $self->o('read_id_tag'),
          csvfile_table => $self->o('summary_csv_table'),
          column_names => $self->o('file_columns'),
          use_threading => $self->o('use_threads'),
        },
        -rc_name    => '1GB',
        -flow_into => {
          '2->A' => ['bwa', 'create_header_files'],
          'A->1' => ['bwa2bam'],
        },
      },
      {
        -logic_name => 'create_header_files',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -rc_name => '1GB',
        -parameters => {
          cmd => 'if [ ! -e "'.$self->o('header_file').'" ]; then printf "'.$header_line.'" > '.$self->o('header_file').'; fi',
        },
      },
      {
        -logic_name => 'bwa',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBWA',
        -parameters => {
          disconnect_jobs => 1,
          short_read_aligner => $self->o('bwa_path'),
          input_dir => $self->o('input_dir'),
          genome_file => $self->o('faidx_genome_file'),
          output_dir => $self->o('output_dir'),
        },
        -flow_into => {
          1 => [ ':////accu?fastq=[]' ],
          -1 => [ 'bwa_20GB' ],
          -2 => [ 'bwa_20GB' ],
        },
        -rc_name    => '10GB_multithread',
      },
      {
        -logic_name => 'bwa_20GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBWA',
        -can_be_empty => 1,
        -parameters => {
          disconnect_jobs => 1,
          short_read_aligner => $self->o('bwa_path'),
          input_dir => $self->o('input_dir'),
          genome_file => $self->o('faidx_genome_file'),
          output_dir => $self->o('output_dir'),
        },
        -flow_into => {
          1 => [ ':////accu?fastq=[]' ],
        },
        -rc_name    => '20GB_multithread',
      },
      {
        -logic_name => 'bwa2bam',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBWA2BAM',
        -parameters => {
          sampe_options => '-A -a '.$self->o('maxintron'),
          samse_options => '',
          min_paired => $self->o('read_min_paired'),
          min_mapped => $self->o('read_min_mapped'),
          header_file => $self->o('header_file'),
          bam_prefix => $self->o('read_id_tag'),
          email => $self->o('email_address'),
          disconnect_jobs => 1,
          short_read_aligner => $self->o('bwa_path'),
          input_dir => $self->o('input_dir'),
          genome_file => $self->o('faidx_genome_file'),
          output_dir => $self->o('output_dir'),
          samtools => $self->o('samtools_path'),
        },
        -flow_into => {
          1 => [ ':////accu?filename=[]' ],
          -1 => {'bwa2bam_20GB' => { fastq => '#fastq#', is_paired => '#is_paired#', $self->o('read_id_tag') => '#'.$self->o('read_id_tag').'#'}},
          -2 => {'bwa2bam_20GB' => { fastq => '#fastq#', is_paired => '#is_paired#', $self->o('read_id_tag') => '#'.$self->o('read_id_tag').'#'}},
        },
        -rc_name    => '10GB',
      },
      {
        -logic_name => 'bwa2bam_20GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBWA2BAM',
        -can_be_empty => 1,
        -parameters => {
          sampe_options => '-A -a '.$self->o('maxintron'),
          samse_options => '',
          min_paired => $self->o('read_min_paired'),
          min_mapped => $self->o('read_min_mapped'),
          header_file => $self->o('header_file'),
          bam_prefix => $self->o('read_id_tag'),
          email => $self->o('email_address'),
          disconnect_jobs => 1,
          short_read_aligner => $self->o('bwa_path'),
          input_dir => $self->o('input_dir'),
          genome_file => $self->o('faidx_genome_file'),
          output_dir => $self->o('output_dir'),
          samtools => $self->o('samtools_path'),
        },
        -flow_into => {
          1 => [ ':////accu?filename=[]' ],
        },
        -rc_name    => '20GB',
      },

      {
        -logic_name => 'merged_tissue_file',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMergeBamFiles',
        -parameters => {
          %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BamMergeStatic', $self->o('rnaseq_merge_type'))},
          # target_db is the database where we will write the files in the data_file table
          # You can use store_datafile => 0, if you don't want to store the output file
          target_db => $self->o('rnaseq_rough_db'),
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
        -rc_name    => '3GB_merged_multithread',
        -flow_into => {
          1 => ['create_analyses_type_job', '?accu_name=filename&accu_address=[]&accu_input_variable=alignment_bam_file'],
        },
      },
      {
        -logic_name => 'create_analyses_type_job',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
        -rc_name    => '1GB',
        -parameters => {
          inputlist => ['gene', 'daf', 'ise'],
          column_names => ['type'],
          species => $self->o('species_name'),
        },
        -flow_into => {
          2 => {'create_rnaseq_tissue_analyses' => {analyses => [{'-logic_name' => '#species#_#sample_name#_rnaseq_#type#'}]}},
        },
      },
      {
        -logic_name => 'create_rnaseq_tissue_analyses',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAddAnalyses',
        -rc_name    => '1GB',
        -parameters => {
          source_type => 'list',
          target_db => $self->o('rnaseq_rough_db'),
        },
      },
      {
        -logic_name => 'merged_bam_file',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMergeBamFiles',
        -parameters => {
          %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BamMergeStatic', $self->o('rnaseq_merge_type'))},
          # target_db is the database where we will write the files in the data_file table
          # You can use store_datafile => 0, if you don't want to store the output file
          target_db => $self->o('rnaseq_rough_db'),
          assembly_name => $self->o('assembly_name'),
          rnaseq_data_provider => $self->o('rnaseq_data_provider'),
          disconnect_jobs => 1,
          alignment_bam_file => catfile($self->o('merge_dir'), '#assembly_name#.#rnaseq_data_provider#.merged.1.bam'),
          species => $self->o('species_name'),
          output_dir => $self->o('merge_dir'),
          input_dir => $self->o('merge_dir'),
          samtools => $self->o('samtools_path'),
          picard_lib_jar => $self->o('picard_lib_jar'),
          use_threads => $self->o('rnaseq_merge_threads'),
        },
        -rc_name    => '5GB_merged_multithread',
        -flow_into => {
          1 => ['fan_merge_analyses'],
          2 => ['create_header_intron'],
        },
      },

     {
        -logic_name => 'fan_merge_analyses',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
	  cmd => "if [[ \$(cut -d\$'\\t' -f1 ".$self->o('rnaseq_summary_file')." | sort | uniq | wc -l) == 1 ]]; then exit 42; else exit 0;fi",
          return_codes_2_branches => {'42' => 2},
        },
        -rc_name    => 'default',
        -flow_into  => {
	  1 => ['create_merge_analyses_type_job'],
        },
      },

      {
        -logic_name => 'create_merge_analyses_type_job',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
        -rc_name    => '1GB',
        -parameters => {
          inputlist => ['gene', 'daf', 'ise'],
          column_names => ['type'],
          species => $self->o('species_name'),
        },
        -flow_into => {
          2 => {'create_rnaseq_merge_analyses' => {analyses => [{'-logic_name' => '#species#_merged_rnaseq_#type#'}]}},
        },
      },
      {
        -logic_name => 'create_rnaseq_merge_analyses',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAddAnalyses',
        -rc_name    => '1GB',
        -parameters => {
          source_type => 'list',
          target_db => $self->o('rnaseq_rough_db'),
        },
      },
      {
        -logic_name => 'create_header_intron',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -rc_name    => '1GB',
        -parameters => {
          cmd => $self->o('samtools_path').' view -H #filename# | grep -v @SQ | grep -v @HD > '.catfile($self->o('rnaseq_dir'),'merged_header.h'),
        },
        -flow_into => {
          '1->A' => [ 'create_toplevel_input_ids'],
          'A->1' => ['sam2bam'],
        },
      },


      {
        -logic_name => 'create_toplevel_input_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -rc_name    => '1GB',
        -parameters => {
          iid_type => 'slice',
          batch_slice_ids => 1,
          batch_target_size => 1000000,
          target_db => $self->o('rnaseq_rough_db'),
        },
        -flow_into => {
          2 => {'create_overlapping_slices' => {'iid' => '#iid#', alignment_bam_file => '#filename#'}},
        },
      },

      {
        -logic_name => 'create_overlapping_slices', # Hopefully this can be removed when the new Bam2Genes module is ready
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -rc_name    => '1GB',
        -parameters => {
          iid_type => 'split_slice',
          slice_size => 5000000,
          slice_overlaps => 2500000,
          target_db => $self->o('rnaseq_rough_db'),
        },
        -flow_into => {
          '2->A' => {'split_on_low_coverage' => {'iid' => '#iid#', alignment_bam_file => '#alignment_bam_file#'}},
          'A->1' => ['check_and_delete_broken_duplicated'],
        },
      },

      {
        -logic_name => 'split_on_low_coverage',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::CalculateLowCoverageSlices',
        -parameters => {
          'disconnect_jobs' => 1,
          'dna_db'   => $self->o('dna_db'),
        },

        -rc_name   => '10GB',
        -flow_into => {
          2 => ['rough_transcripts'],
          -1 => ['split_on_low_coverage_20GB'],
        },

      },

      {
        -logic_name => 'split_on_low_coverage_20GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::CalculateLowCoverageSlices',
        -parameters => {
          'disconnect_jobs' => 1,
          'dna_db'   => $self->o('dna_db'),
        },

        -rc_name          => '20GB',
        -flow_into => {
          2 => ['rough_transcripts'],
        },

      },


      {
        -logic_name => 'rough_transcripts',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBam2Genes',
        -parameters => {
          logic_name => 'rough_transcripts',
          target_db    => $self->o('rnaseq_rough_db'),
          dna_db    => $self->o('dna_db'),
          min_length => 300,
          min_exons  =>   1,
          max_intron_length => $self->o('maxintron'),
          min_single_exon_length => 1000,
          min_span   =>   1.5,
          paired => $self->o('paired_end_only'),
          use_ucsc_naming => $self->o('use_ucsc_naming'),
        },
        -rc_name    => '2GB',
        -flow_into => {
          -1 => {'rough_transcripts_10GB' => {'iid' => '#iid#', alignment_bam_file => '#alignment_bam_file#'}},
          -2 => {'failed_rough_transcript_jobs' => {'iid' => '#iid#', alignment_bam_file => '#alignment_bam_file#'}},
          -3 => {'failed_rough_transcript_jobs' => {'iid' => '#iid#', alignment_bam_file => '#alignment_bam_file#'}},
        },
      },


      {
        -logic_name => 'rough_transcripts_10GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBam2Genes',
        -can_be_empty => 1,
        -parameters => {
          logic_name => 'rough_transcripts',
          target_db    => $self->o('rnaseq_rough_db'),
          dna_db    => $self->o('dna_db'),
          min_length => 300,
          min_exons  =>   1,
          max_intron_length => $self->o('maxintron'),
          min_single_exon_length => 1000,
          min_span   =>   1.5,
          paired => $self->o('paired_end_only'),
          use_ucsc_naming => $self->o('use_ucsc_naming'),
        },
        -rc_name    => '10GB',
        -flow_into => {
          -1 => {'rough_transcripts_30GB' => {'iid' => '#iid#', alignment_bam_file => '#alignment_bam_file#'}},
          -2 => {'failed_rough_transcript_jobs' => {'iid' => '#iid#', alignment_bam_file => '#alignment_bam_file#'}},
          -3 => {'failed_rough_transcript_jobs' => {'iid' => '#iid#', alignment_bam_file => '#alignment_bam_file#'}},
        },
      },


      {
        -logic_name => 'rough_transcripts_30GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBam2Genes',
        -can_be_empty => 1,
        -parameters => {
          logic_name => 'rough_transcripts',
          target_db    => $self->o('rnaseq_rough_db'),
          dna_db    => $self->o('dna_db'),
          min_length => 300,
          min_exons  =>   1,
          max_intron_length => $self->o('maxintron'),
          min_single_exon_length => 1000,
          min_span   =>   1.5,
          paired => $self->o('paired_end_only'),
          use_ucsc_naming => $self->o('use_ucsc_naming'),
        },
        -rc_name    => '30GB',
        -flow_into => {
          -1 => {'failed_rough_transcript_jobs' => {'iid' => '#iid#', alignment_bam_file => '#alignment_bam_file#'}},
          -2 => {'failed_rough_transcript_jobs' => {'iid' => '#iid#', alignment_bam_file => '#alignment_bam_file#'}},
          -3 => {'failed_rough_transcript_jobs' => {'iid' => '#iid#', alignment_bam_file => '#alignment_bam_file#'}},
        },
      },


      {
        -logic_name => 'failed_rough_transcript_jobs',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name          => 'default',
        -can_be_empty  => 1,
      },


      {
        -logic_name => 'check_and_delete_broken_duplicated',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRemoveBrokenAndDuplicatedObjects',
        -parameters => {
          target_db => $self->o('rnaseq_rough_db'),
          check_support => 0,
        },
        -rc_name    => '5GB',
        -analysis_capacity => 1, # Because there is slice overlap, having parallel jobs can cause problems
        -flow_into => {
          1 => ['create_bam2introns_input_ids'],
        },
      },

      {
        -logic_name => 'create_bam2introns_input_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
          iid_type => 'slice_to_feature_ids',
          target_db => $self->o('rnaseq_rough_db'),
          feature_type => 'gene',
          logic_name => ['rough_transcripts'],
          use_stable_ids => 1,
          create_stable_ids => 1,
          stable_id_prefix => 'RNASEQ',
        },
        -rc_name    => '1GB',
        -batch_size => 100,
        -flow_into => {
          2 => {'bam2introns' => {iid => '#iid#', bam_file => '#alignment_bam_file#'}},
        },
      },

      {
        -logic_name => 'bam2introns',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBam2Introns',
        -parameters => {
          program_file => $self->o('exonerate_path'),
          source_db => $self->o('rnaseq_rough_db'),
          dna_db => $self->o('dna_db'),
          missmatch => 6,
          word_length => 10,
          saturate_threshold => 10000,
          mask => 1,
          percent_id => 97,
          coverage => 90,
          fullseq   => 1,
          max_transcript => 1000000,
          batch_size => 10000,
          maxintron => $self->o('maxintron'),
          use_ucsc_naming => $self->o('use_ucsc_naming'),
          output_dir => $self->o('sam_dir'),
          use_genome_flatfile => 1,
          flatfile_masked => 1,
          genome_file => $self->o('faidx_softmasked_genome_file'),
          timer => '15m',
        },
        -rc_name    => '5GB',
        -batch_size => 200,
        -analysis_capacity => 500,
        -flow_into => {
          1 => [':////accu?filename=[]'],
          2 => {'bam2introns' => {iid => '#iid#', bam_file => '#bam_file#'}},
          -1 => {'bam2introns_20GB' => {iid => '#iid#', bam_file => '#bam_file#'}},
          -2 => ['failed_bam2introns_jobs'],
          -3 => ['failed_bam2introns_jobs'],
        },
      },


      {
        -logic_name => 'bam2introns_20GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBam2Introns',
        -can_be_empty => 1,
        -parameters => {
          program_file => $self->o('exonerate_path'),
          source_db => $self->o('rnaseq_rough_db'),
          dna_db => $self->o('dna_db'),
          missmatch => 6,
          word_length => 10,
          saturate_threshold => 10000,
          mask => 1,
          percent_id => 97,
          coverage => 90,
          fullseq   => 1,
          max_transcript => 1000000,
          batch_size => 10000,
          maxintron => $self->o('maxintron'),
          use_ucsc_naming => $self->o('use_ucsc_naming'),
          output_dir => $self->o('sam_dir'),
          use_genome_flatfile => 1,
          flatfile_masked => 1,
          genome_file => $self->o('faidx_softmasked_genome_file'),
          timer => '30m',
        },
        -rc_name    => '20GB',
        -flow_into => {
          1 => [':////accu?filename=[]'],
          2 => {'bam2introns' => {iid => '#iid#', bam_file => '#bam_file#'}},
          -1 => {'bam2introns_30GB' => {iid => '#iid#', bam_file => '#bam_file#'}},
          -2 => ['failed_bam2introns_jobs'],
          -3 => ['failed_bam2introns_jobs'],
        },
      },

      {
        -logic_name => 'bam2introns_30GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBam2Introns',
        -can_be_empty => 1,
        -parameters => {
          program_file => $self->o('exonerate_path'),
          source_db => $self->o('rnaseq_rough_db'),
          dna_db => $self->o('dna_db'),
          missmatch => 6,
          word_length => 10,
          saturate_threshold => 10000,
          mask => 1,
          percent_id => 97,
          coverage => 90,
          fullseq   => 1,
          max_transcript => 1000000,
          batch_size => 10000,
          maxintron => $self->o('maxintron'),
          use_ucsc_naming => $self->o('use_ucsc_naming'),
          output_dir => $self->o('sam_dir'),
          use_genome_flatfile => 1,
          flatfile_masked => 1,
          genome_file => $self->o('faidx_softmasked_genome_file'),
          timer => '30m',
        },
        -rc_name    => '30GB',
        -flow_into => {
          1 => [':////accu?filename=[]'],
          2 => {'bam2introns' => {iid => '#iid#', bam_file => '#bam_file#'}},
          -1 => ['failed_bam2introns_jobs'],
          -2 => ['failed_bam2introns_jobs'],
          -3 => ['failed_bam2introns_jobs'],
        },
      },

      {
        -logic_name => 'failed_bam2introns_jobs',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name          => 'default',
        -can_be_empty  => 1,
      },

      {
        -logic_name => 'sam2bam',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSam2Bam',
        -parameters => {
          headerfile => catfile($self->o('rnaseq_dir'), 'merged_header.h'),
          disconnect_jobs => 1,
          samtools => $self->o('samtools_path'),
          intron_bam_file => catfile($self->o('output_dir'), 'introns'),
          genome_file => $self->o('faidx_genome_file'),
          use_threading => $self->o('use_threads'),
          sam_dir => $self->o('sam_dir'),
        },
        -rc_name    => '5GB_multithread',
        -flow_into => ['create_refine_db'],
      },


     {
        -logic_name => 'create_refine_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
          source_db => $self->o('rnaseq_rough_db'),
          target_db => $self->o('rnaseq_refine_db'),
          create_type => 'clone',
          extra_data_tables => ['data_file'],
        },
        -rc_name => '1GB',
        -flow_into => ['create_blast_db'],
      },

      {
        -logic_name => 'create_blast_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
          source_db => $self->o('rnaseq_refine_db'),
          target_db => $self->o('rnaseq_blast_db'),
          create_type => 'clone',
          extra_data_tables => ['data_file'],
        },
        -rc_name => '1GB',
        -flow_into => ['create_ccode_config'],
      },
      {
        -logic_name => 'create_ccode_config',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenerateRefineConfig',
        -parameters => {
          species => $self->o('species_name'),
          output_dir => $self->o('output_dir'),
          intron_bam_file => catfile($self->o('output_dir'), 'introns'),
          single_tissue => $self->o('single_tissue'),
          sample_column => $self->o('read_group_tag'),
          sample_id_column => $self->o('read_id_tag'),
          csvfile_table => $self->o('summary_csv_table'),
          source_db => $self->o('rnaseq_rough_db'),
          dna_db => $self->o('dna_db'),
          target_db => $self->o('rnaseq_refine_db'),
          # write the intron features into the OUTPUT_DB along with the models
          write_introns => 1,
          # maximum number of times to loop when building all possible paths through the transcript
          max_recursions => 10000000000000,
          # analysis logic_name for the dna_align_features to fetch from the INTRON_DB
          # If left blank all features will be fetched
          logicname => [],
          # logic name of the gene models to fetch
          model_ln  => '',
          # penalty for removing a retined intron
          retained_intron_penalty => 2,
          #Remove introns that overlap X introns
          filter_on_overlap => 0,
          # minimum size for an intron
          min_intron_size  => 30,
          max_intron_size  => $self->o('maxintron'),
          # biotype to give to single exon models if left blank single exons are ignored
          # minimum single exon size (bp)
          min_single_exon => 1000,
          # minimum percentage of single exon length that is coding
          single_exon_cds => 66,
          # Intron with most support determines the splice sites for an internal exon
          # lower scoring introns with different splice sites are rejected
          strict_internal_splice_sites => 1,
          # In some species alternate splice sites for end exons seem to be common
          strict_internal_end_exon_splice_sites => 1,
          # biotypes to give gene models if left blank these models will not get written to the output database
          # best score - model with most supporting intron features
          # all other possible models
          # max number of other models to make - blank = all
          other_num      => '10',
          # max number of other models to process - blank = all
          max_num      => '1000',
          other_isoforms => $self->o('other_isoforms'),
          # biotype to label bad models ( otherwise they are not written )
          # do you want to trim UTR
          trim_utr => 1,
          # config for trimming UTR
          max_3prime_exons => 2,
          max_3prime_length => 5000,
          max_5prime_exons => 3,
          max_5prime_length => 1000,
          # % of average intron score that a UTR intron must have
          reject_intron_cutoff => 5,
        },

        -rc_name          => '1GB',
        -flow_into => {
          '2->A' => ['create_ccode_input_ids'],
          'A->1' => ['copy_rnaseq_blast_db'],
        },
      },
      {
        -logic_name => 'create_ccode_input_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -rc_name    => '1GB',
        -parameters => {
          iid_type => 'slice',
          coord_system_name => 'toplevel',
          slice => 1,
          include_non_reference => 0,
          top_level => 1,
          feature_constraint => 1,
          feature_type => 'gene',
          target_db => $self->o('rnaseq_rough_db'),
        },
        -flow_into => {
          2 => {'refine_genes' => {iid => '#iid#', logic_name => '#logic_name#', config_file => '#config_file#'}},
        },
      },

      {
        -logic_name => 'refine_genes',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => $self->o('refine_ccode_exe').($self->o('use_ucsc_naming') ? ' -u ' : ' ').($self->o('use_threads') ? ' -t '.$self->o('use_threads').' ' : ' ').'-c #config_file# -i #iid# -l #logic_name# -v 0',
          return_codes_2_branches => {
            42 => 2,
          },
        },
        -rc_name => '2GB_multithread',
        -flow_into => {
          1 => ['create_gene_id_input_ids'],
          -1 => {'refine_genes_20GB' => {iid => '#iid#', config_file => '#config_file#', logic_name => '#logic_name#'}},
        },
      },
      {
        -logic_name => 'refine_genes_20GB',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => $self->o('refine_ccode_exe').($self->o('use_ucsc_naming') ? ' -u ' : ' ').($self->o('use_threads') ? ' -t '.$self->o('use_threads').' ' : ' ').'-c #config_file# -i #iid# -l #logic_name# -v 0',
          return_codes_2_branches => {
            42 => 2,
          },
        },
        -rc_name => '20GB_multithread',
        -flow_into => {
          1 => ['create_gene_id_input_ids'],
        },
      },


      {
        -logic_name => 'create_gene_id_input_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -rc_name    => '1GB',
        -parameters => {
          iid_type => 'feature_id',
          coord_system_name => 'toplevel',
          target_db => $self->o('rnaseq_refine_db'),
          feature_logic_names => ['#logic_name#'],
          feature_type => 'gene',
          batch_size => 50,
        },
        -flow_into => {
          2 => {'blast_rnaseq' => {iid => '#iid#', logic_name => '#logic_name#'}},
        },
      },


      {
        -logic_name => 'blast_rnaseq',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastRNASeqPep',
        -parameters => {
          source_db => $self->o('rnaseq_refine_db'),
          target_db => $self->o('rnaseq_blast_db'),
          dna_db => $self->o('dna_db'),
          iid_type => 'object_id',
          # path to index to fetch the sequence of the blast hit to calculate % coverage
          indicate_index => $self->o('protein_blast_index'),
          uniprot_index => [$self->o('protein_blast_db')],
          blast_program => $self->o('uniprot_blast_exe_path'),
          %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastGenscanPep', {BLAST_PARAMS => {-type => $self->o('blast_type')}})},
          commandline_params => $self->o('blast_type') eq 'wu' ? '-cpus='.$self->o('use_threads').' -hitdist=40' : '-num_threads '.$self->o('use_threads').' -window_size 40',
        },
        -rc_name => 'blast',
      },

      {
        -logic_name => 'copy_rnaseq_blast_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('rnaseq_blast_db'),
                         target_db => $self->o('rnaseq_for_layer_db'),
                         create_type => 'copy',
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '1' => ['update_rnaseq_for_layer_biotypes'],
                      },
      },


      {
        -logic_name => 'update_rnaseq_for_layer_biotypes',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('rnaseq_for_layer_db'),
          sql => [
            'UPDATE gene SET biotype = "rnaseq_merged" WHERE biotype IN ("best","single","other_merged")',
            'UPDATE gene SET biotype = "rnaseq_tissue" WHERE biotype != "rnaseq_merged"',
            'UPDATE transcript JOIN gene USING(gene_id) SET transcript.biotype = gene.biotype',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
                        '1' => ['remove_rnaseq_for_layer_daf_features'],
                      },
      },


      {
        -logic_name => 'remove_rnaseq_for_layer_daf_features',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('rnaseq_for_layer_db'),
          sql => [
            'TRUNCATE dna_align_feature',
            'DELETE transcript_supporting_feature FROM transcript_supporting_feature WHERE feature_type = "dna_align_feature"',
            'DELETE supporting_feature FROM supporting_feature WHERE feature_type = "dna_align_feature"',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
                        '1' => ['classify_rnaseq_for_layer_models'],
                      },
      },


      {
        -logic_name => 'classify_rnaseq_for_layer_models',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport',
        -parameters => {
                         update_gene_biotype => 1,
                         classification_type => 'standard',
                         target_db => $self->o('rnaseq_for_layer_db'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['rnaseq_for_layer_sanity_checks'],
                      },

      },


      {
        -logic_name => 'rnaseq_for_layer_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         target_db => $self->o('rnaseq_for_layer_db'),
                         sanity_check_type => 'gene_db_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'gene_db_checks')->{$self->o('uniprot_set')}->{'rnaseq_blast'},
                       },
        -flow_into => {
          1 => ['create_rnaseq_layer_nr_db'],
        },
        -rc_name    => '4GB',
      },


      {
        -logic_name => 'create_rnaseq_layer_nr_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('rnaseq_for_layer_db'),
                         target_db => $self->o('rnaseq_for_layer_nr_db'),
                         create_type => 'copy',
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '1' => ['create_rnaseq_layer_nr_slices'],
                      },
      },


      {
        -logic_name => 'create_rnaseq_layer_nr_slices',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db        => $self->o('dna_db'),
                         coord_system_name => 'toplevel',
                         iid_type => 'slice',
                         slice_size => 20000000,
                         include_non_reference => 0,
                         top_level => 1,
                         min_slice_length => $self->o('min_toplevel_slice_length'),
                         batch_slice_ids => 1,
                         batch_target_size => 20000000,
                       },
        -rc_name    => '2GB',
        -flow_into => {
                         '2'    => ['remove_redundant_rnaseq_layer_genes'],
                      },
      },


    {
      -logic_name => 'remove_redundant_rnaseq_layer_genes',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::RemoveRedundantGenes',
      -parameters => {
        target_db => $self->o('rnaseq_for_layer_nr_db'),
        target_type => 'generic',
      },
      -rc_name          => '5GB',
    },


# STAR rnaseq
	    {
        -logic_name => 'parse_summary_file_star',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveParseCsvIntoTable',
        -rc_name => '1GB',
        -parameters => {
          column_names => $self->o('file_columns'),
          sample_column => $self->o('read_group_tag'),
          inputfile => $self->o('star_rnaseq_summary_file'),
          delimiter => $self->o('summary_file_delimiter'),
          csvfile_table => $self->o('summary_csv_table'),
          pairing_regex => $self->o('pairing_regex'),
          read_length_table => $self->o('read_length_table'),
        },
        -flow_into => {
          '1->A' => ['create_star_jobs'],
          'A->1' => ['stringtie2merge'],
        },
      },


	    {
        -logic_name => 'create_star_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateStarJobs',
        -parameters => {
                          input_dir => $self->o('input_dir'),
          sample_column => $self->o('read_group_tag'),
          sample_id_column => $self->o('read_id_tag'),
          filename_column => $self->o('filename_tag'),
          csvfile_table => $self->o('summary_csv_table'),
          column_names => $self->o('file_columns'),
          use_threading => $self->o('use_threads'),
        },
        -rc_name    => '1GB',
        -flow_into => {
          2 => ['star'],
        },
      },

	    {
	     -logic_name => 'star',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveStar',
        -parameters => {
          disconnect_jobs => 1,
          input_dir => $self->o('input_dir'),
          output_dir => $self->o('output_dir'),
          short_read_aligner => $self->o('star_path'),
          genome_dir => catfile($self->o('output_path'),'genome_dumps'),
          num_threads => $self->o('star_threads'),
        },
        -flow_into => {
          2 => ['stringtie2'],
        },
        -rc_name    => '45GB_star',
      },


	    {
        -logic_name => 'stringtie2',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::Stringtie2',
         -parameters => {
           output_dir => catdir($self->o('output_dir'),'stringtie'),
           stringtie2_path        => $self->o('stringtie2_path'),
           csv_summary_file       => $self->o('star_rnaseq_summary_file'),
           csv_summary_file_genus => $self->o('rnaseq_summary_file_genus'),
           num_threads => $self->o('stringtie_threads'),
        },
        -rc_name    => '10GB_stringtie',
      },

	    {
        -logic_name => 'stringtie2merge',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::Stringtie2Merge',
         -parameters => {
           stringtie_gtf_dir => catdir($self->o('output_dir'),'stringtie'),
           stringtie_merge_dir => catdir($self->o('output_dir'),'stringtie','merge'),
           stringtie2_path        => $self->o('stringtie2_path'),
           csv_summary_file       => $self->o('star_rnaseq_summary_file'),
           csv_summary_file_genus => $self->o('rnaseq_summary_file_genus'),
           num_threads => $self->o('stringtie_threads'),

                         },
         -flow_into => {
          1 => ['create_stringtie_initial_db'],
        },
        -rc_name    => '10GB_stringtie',
      },

	    {
        -logic_name => 'create_stringtie_initial_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('stringtie_initial_db'),
                         create_type => 'clone',
                       },
        -rc_name    => 'default',
        -flow_into => {
          1 => ['create_stringtie_blast_db'],
        },
      },

	    {
        -logic_name => 'create_stringtie_blast_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('stringtie_blast_db'),
                         create_type => 'clone',
                       },
        -rc_name    => 'default',
        -flow_into => {
           '1->A' => ['generate_stringtie_gtf_jobs'],
           'A->1' => ['star2introns'],
        },
      },

	    {
        -logic_name => 'generate_stringtie_gtf_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::GenerateGTFLoadingJobs',
         -parameters => {
           gtf_dir => catdir($self->o('output_dir'),'stringtie','merge'),
        },
        -flow_into => {
          2 => ['load_stringtie_transcripts'],
        },
        -rc_name    => '5GB',
      },

	    {
        -logic_name => 'load_stringtie_transcripts',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::LoadGTFBasic',
        -parameters => {
          target_db => $self->o('stringtie_initial_db'),
          loading_type => 'range',
          genome_file => $self->o('faidx_genome_file'),
          logic_name  => 'stringtie2',
          module      => 'Stringtie2',
        },
        -rc_name    => '5GB',
        -flow_into => {
          1 => 'create_slice_tissue_input_ids',
        },
      },

	    {
        -logic_name => 'create_slice_tissue_input_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
	     -rc_name    => '1GB',
        -parameters => {
          iid_type => 'slice',
          coord_system_name => 'toplevel',
          slice => 1,
          include_non_reference => 0,
          top_level => 1,
          feature_constraint => 1,
          feature_type => 'gene',
          target_db => $self->o('stringtie_initial_db'),
        },
        -flow_into => {
          2 => {'create_gene_id_blast_stringtie_input_ids' => {iid => '#iid#', logic_name => '#logic_name#'}},
        },
      },

	    {
        -logic_name => 'create_gene_id_blast_stringtie_input_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -rc_name    => '1GB',
        -parameters => {
          iid_type => 'feature_id',
          coord_system_name => 'toplevel',
          target_db => $self->o('stringtie_initial_db'),
          feature_type => 'gene',
          batch_size => 50,
          feature_logic_names => ['#logic_name#'],
        },
        -flow_into => {
          2 => {'blast_stringtie' => {iid => '#iid#', logic_name => '#logic_name#'}},
        },
      },


	    {
        -logic_name => 'blast_stringtie',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastRNASeqPep',
        -parameters => {
          source_db => $self->o('stringtie_initial_db'),
          target_db => $self->o('stringtie_blast_db'),
          dna_db => $self->o('dna_db'),
          iid_type => 'object_id',
          # path to index to fetch the sequence of the blast hit to calculate % coverage
          indicate_index => $self->o('protein_blast_index'),
          uniprot_index => [$self->o('protein_blast_db')],
          blast_program => $self->o('uniprot_blast_exe_path'),
          %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastGenscanPep', {BLAST_PARAMS => {-type => $self->o('blast_type')}})},
          commandline_params => $self->o('blast_type') eq 'wu' ? '-cpus='.$self->o('use_threads').' -hitdist=40' : '-num_threads '.$self->o('use_threads').' -window_size 40',
        },
        -flow_into => {
         '-1' => ['blast_stringtie_longseq'],
          '2' => ['blast_stringtie_longseq'],
        },
        -rc_name => 'blast',
      },

	    {
        -logic_name => 'blast_stringtie_longseq',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastRNASeqPep',
        -parameters => {
          source_db => $self->o('stringtie_initial_db'),
          target_db => $self->o('stringtie_blast_db'),
          dna_db => $self->o('dna_db'),
                               iid_type => 'object_id',
                               # path to index to fetch the sequence of the blast hit to calculate % coverage
                               indicate_index => $self->o('protein_blast_index'),
          uniprot_index => [$self->o('protein_blast_db')],
          blast_program => $self->o('uniprot_blast_exe_path'),
          %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastGenscanPep', {BLAST_PARAMS => {-type => $self->o('blast_type')}})},
          commandline_params => $self->o('blast_type') eq 'wu' ? '-cpus='.$self->o('use_threads').' -hitdist=40' : '-num_threads '.$self->o('use_threads').' -window_size 40',
                             },
        -rc_name => 'blast10GB',
      },

	    {
        -logic_name => 'star2introns',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveStar2Introns',
        -parameters => {
                        star_junctions_dir => $self->o('output_dir'),
                        intron_db => $self->o('stringtie_blast_db'),
                        source_db => $self->o('dna_db'),
			sample_column => 'SM',
			sample_id_column => 'ID',
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '1' => ['copy_rnaseq_blast_db_star'],
                      },
      },

	    {
        -logic_name => 'copy_rnaseq_blast_db_star',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('stringtie_blast_db'),
                         target_db => $self->o('star_rnaseq_for_layer_db'),
                         create_type => 'copy',
                         force_drop => 1,
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '1' => ['update_rnaseq_for_layer_biotypes_star'],
                      },
      },

	    {
        -logic_name => 'update_rnaseq_for_layer_biotypes_star',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('star_rnaseq_for_layer_db'),
          sql => [
            'UPDATE gene SET biotype = "rnaseq_merged" WHERE source IN ("merged")',
            'UPDATE gene SET biotype = "rnaseq_tissue" WHERE biotype != "rnaseq_merged"',
            'UPDATE transcript JOIN gene USING(gene_id) SET transcript.biotype = gene.biotype',
          ],
        },
        -rc_name    => 'default',
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
                         target_db => $self->o('star_rnaseq_for_layer_db'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['rnaseq_for_layer_sanity_checks_star'],
                      },

      },

	    {
        -logic_name => 'rnaseq_for_layer_sanity_checks_star',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         target_db => $self->o('star_rnaseq_for_layer_db'),
                         sanity_check_type => 'gene_db_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'gene_db_checks')->{$self->o('uniprot_set')}->{'rnaseq_blast'},
                       },

	     -flow_into => {
          1 => ['create_rnaseq_layer_nr_db_star'],
        },
     },


	    {
	     -logic_name => 'create_rnaseq_layer_nr_db_star',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('star_rnaseq_for_layer_db'),
                         target_db => $self->o('star_rnaseq_for_layer_nr_db'),
                         create_type => 'copy',
                         force_drop => 1,
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '1' => ['create_rnaseq_layer_nr_slices_star'],
                      },
      },

     {
	     -logic_name => 'create_rnaseq_layer_nr_slices_star',
	     -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db        => $self->o('dna_db'),
                         coord_system_name => 'toplevel',
                         iid_type => 'slice',
                         slice_size => 20000000,
                         include_non_reference => 0,
                         top_level => 1,
                         min_slice_length => $self->o('min_toplevel_slice_length'),
                         batch_slice_ids => 1,
                         batch_target_size => 20000000,
                       },
        -rc_name    => '2GB',
        -flow_into => {
                         '2'    => ['remove_redundant_rnaseq_layer_genes_star'],
                      },
      },


     {
      -logic_name => 'remove_redundant_rnaseq_layer_genes_star',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::RemoveRedundantGenes',
      -parameters => {
        target_db => $self->o('star_rnaseq_for_layer_nr_db'),
        target_type => 'generic',
      },
      -rc_name          => '5GB',
    },

    ];
}


sub resource_classes {
  my $self = shift;

  return {
    '1GB' => { LSF => $self->lsf_resource_builder('production-rh74', 1000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '2GB' => { LSF => $self->lsf_resource_builder('production-rh74', 2000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '4GB' => { LSF => $self->lsf_resource_builder('production-rh74', 4000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '5GB' => { LSF => $self->lsf_resource_builder('production-rh74', 5000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '10GB' => { LSF => $self->lsf_resource_builder('production-rh74', 10000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '20GB' => { LSF => $self->lsf_resource_builder('production-rh74', 20000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '30GB' => { LSF => $self->lsf_resource_builder('production-rh74', 30000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'default' => { LSF => $self->lsf_resource_builder('production-rh74', 900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'blast' => { LSF => $self->lsf_resource_builder('production-rh74', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], undef, 3)},
    'blast10GB' => { LSF => $self->lsf_resource_builder('production-rh74', 10000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], undef, 3)},
    '2GB_multithread' => { LSF => $self->lsf_resource_builder('production-rh74', 2000, [$self->default_options->{'pipe_db_server'}], undef, $self->default_options->{'use_threads'})},
    '3GB_merged_multithread' => { LSF => $self->lsf_resource_builder('production-rh74', 3000, [$self->default_options->{'pipe_db_server'}], undef, $self->default_options->{'rnaseq_merge_threads'})},
    '5GB_merged_multithread' => { LSF => $self->lsf_resource_builder('production-rh74', 5000, [$self->default_options->{'pipe_db_server'}], undef, ($self->default_options->{'rnaseq_merge_threads'}))},
    '5GB_multithread' => { LSF => $self->lsf_resource_builder('production-rh74', 5000, [$self->default_options->{'pipe_db_server'}], undef, ($self->default_options->{'use_threads'}+1))},
    '10GB_multithread' => { LSF => $self->lsf_resource_builder('production-rh74', 10000, [$self->default_options->{'pipe_db_server'}], undef, ($self->default_options->{'use_threads'}+1))},
    '20GB_multithread' => { LSF => $self->lsf_resource_builder('production-rh74', 20000, [$self->default_options->{'pipe_db_server'}], undef, ($self->default_options->{'use_threads'}+1))},
    '45GB_star' => { LSF => $self->lsf_resource_builder('production-rh74', 45000, undef, undef, ($self->default_options->{'star_threads'}+1))},
    '10GB_stringtie' => { LSF => $self->lsf_resource_builder('production-rh74', 10000, undef, undef, ($self->default_options->{'stringtie_threads'}))},
  }
}

sub hive_capacity_classes {
  my $self = shift;

  return {
           'hc_very_low'    => 35,
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
