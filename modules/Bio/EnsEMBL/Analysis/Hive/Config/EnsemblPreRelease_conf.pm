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

package EnsemblPreRelease_conf;

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
    #BUSCO parameters
    'busco_singularity_image'   => '/hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/busco_v5.8.2_cv1.sif',
    'busco_download_path'       => '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/data/busco_data/data_odb12/',
    'helixer_singularity_image' => '/hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/helixer-docker_helixer_v0.3.5_cuda_12.2.2-cudnn8.sif',
    'gffread_path' => '/hps/software/users/ensembl/genebuild/genebuild_virtual_user/bin/gffread',
    'current_genebuild'            => 0,
    'cores'                        => 30,
    'num_threads'                  => 20,
    'gpu'                          => 'gpu:a100:2',
    'dbowner'                      => '' || $ENV{EHIVE_USER} || $ENV{USER},
    'base_output_dir'              => '',
    'init_config'               => '', #path for configuration file (custom loading)
    'override_clade'               => '', #optional, already defined in ProcessGCA
    'protein_file'                 => '', #optional, already defined in ProcessGCA
    'busco_protein_file'           => '', #optional, already defined in ProcessGCA
    'rfam_accessions_file'         => '', #optional, already defined in ProcessGCA
    'use_existing_short_read_dir'  => '', #path for esisting short read data
    'registry_file'                => catfile( $self->o('enscode_root_dir'),'ensembl-analysis/scripts/genebuild/gbiab/support_files/Databases.pm' ), # This should be the path to the pipeline's copy of the Databases.pm registry file, core adaptors will be written to it
    'generic_registry_file'        => '',                                                                                                                # Could use this to hold the path to ensembl-analysis/scripts/genebuild/gbiab/support_files/Databases.pm to copy as a generic registry
    'diamond_validation_db'        => '/hps/nobackup/flicek/ensembl/genebuild/blastdb/uniprot_euk_diamond/uniprot_euk.fa.dmnd',
    'validation_type'              => 'moderate',
    'release_number'               => '114' || $self->o('ensembl_release'),
    'production_name'              => '' || $self->o('species_name'),
    'pipeline_name'                => '' || $self->o('production_name') . $self->o('production_name_modifier'),
    'user_r'                       => 'ensro',                                                                                                                # read only db user
    'user'                         => 'ensadmin',                                                                                                                # write db user
    'password'                     => '',                                                                                                                # password for write db user
    'server_set'                   => '',                                                                                                                # What server set to user, e.g. set1
    'busco_input_file_stid'        => 'stable_id_to_dump.txt',
    'species_name'                 => '', #optional, already defined in ProcessGCA e.g. mus_musculus
    'taxon_id'                     => '', #optional, already defined in ProcessGCA, should be in the assembly report file
    'species_taxon_id'             => '' || $self->o('taxon_id'),                                                                                        # Species level id, could be different to taxon_id if we have a subspecies, used to get species level RNA-seq CSV data
    'genus_taxon_id'               => '' || $self->o('taxon_id'),                                                                                        # Genus level taxon id, used to get a genus level csv file in case there is not enough species level transcriptomic data
    'uniprot_set'                  => '', #optional, already defined in ProcessGCA e.g. mammals_basic, check UniProtCladeDownloadStatic.pm module in hive config dir for suitable set,
    'output_path'                  => '', #optional, already defined in ProcessGCA
    'assembly_name'                => '', #optional aleady defined in the registry
    'assembly_accession'           => '', #the pipeline is initialed via standalone job  # Versioned GCA assembly accession, e.g. GCA_001857705.1
    'stable_id_prefix'             => '', #optional, already defined in ProcessGCA
    'use_genome_flatfile'          => '1',# This will read sequence where possible from a dumped flatfile instead of the core db
    'species_url'                  => '' || $self->o('production_name') . $self->o('production_name_modifier'),                                          # sets species.url meta key
    'species_division'             => '', #optional, already defined in ProcessGCA # sets species.division meta key
    'stable_id_start'              => '', #optional, already defined in ProcessGCA When mapping is not required this is usually set to 0
    'mapping_required'             => '0',# If set to 1 this will run stable_id mapping sometime in the future. At the moment it does nothing
    'uniprot_version'              => 'uniprot_2021_04',                                                                                                 # What UniProt data dir to use for various analyses
    'production_name_modifier'     => '',                                                                                                                # Do not set unless working with non-reference strains, breeds etc. Must include _ in modifier, e.g. _hni for medaka strain HNI

    # Keys for custom loading, only set/modify if that's what you're doing
    'load_toplevel_only'        => '1',                                                                                                                  # This will not load the assembly info and will instead take any chromosomes, unplaced and unlocalised scaffolds directly in the DNA table
    'custom_toplevel_file_path' => '',                                                                                                                   # Only set this if you are loading a custom toplevel, requires load_toplevel_only to also be set to 2
    'repeatmodeler_library'     => '', #no needed, it can be an option for the anno command This should be the path to a custom repeat library, leave blank if none exists
    'base_blast_db_path'    => $ENV{BLASTDB_DIR},
    'protein_entry_loc'         => catfile( $self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), 'entry_loc' ),                       # Used by genscan blasts and optimise daf/paf. Don't change unless you know what you're doing

    'softmask_logic_names' => [],

    # busco threshold for the analysis that checks wether produce pre-release files or not!
    'busco_threshold' => 70, # If the busco score is above this threshold, the pre-release files will be produced
    'busco_lower_threshold' => 50, # If the busco score is above this threshod and the difference less than 'busco_difference_threshold', the pre-release files will be produced
    'busco_difference_threshold' => 10, # If the difference between the gene and protein busco score is less than this value, the pre-release files will be produced as long as the busco score is above 'busco_lower_threshold'
    
    
    #gff file dump options
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
# Pipe and ref db info
########################


    'provider_name' => 'Ensembl',
    'provider_url'  => 'www.ensembl.org',

    'pipe_db_name' => $self->o('dbowner') . '_' . $self->o('pipeline_name') . '_pipe_' . $self->o('release_number'),
    'dna_db_name' => $self->o('dbowner') . '_' . $self->o('production_name') . $self->o('production_name_modifier') . '_core_' . $self->o('release_number'),


    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'ensembl_release'      => $ENV{ENSEMBL_RELEASE},     # this is the current release version on staging to be able to get the correct database
    'production_db_server' => 'mysql-ens-meta-prod-1',
    'production_db_port'   => '4483',


    ensembl_analysis_script           => catdir( $self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts' ),
    load_optimise_script              => catfile( $self->o('ensembl_analysis_script'), 'genebuild', 'load_external_db_ids_and_optimize_af.pl' ),
    remove_duplicates_script          => catfile( $self->o('ensembl_analysis_script'), 'find_and_remove_duplicates.pl' ),
    ensembl_misc_script               => catdir( $self->o('enscode_root_dir'),        'ensembl',   'misc-scripts' ),
    meta_coord_script                 => catfile( $self->o('ensembl_misc_script'),     'meta_coord', 'update_meta_coord.pl' ),
    meta_levels_script                => catfile( $self->o('ensembl_misc_script'),     'meta_levels.pl' ),
    frameshift_attrib_script          => catfile( $self->o('ensembl_misc_script'),     'frameshift_transcript_attribs.pl' ),
    select_canonical_script           => catfile( $self->o('ensembl_misc_script'),     'canonical_transcripts', 'select_canonical_transcripts.pl' ),
    print_protein_script_path         => catfile( $self->o('ensembl_analysis_script'), 'genebuild', 'print_translations.pl' ),
    ensembl_gst_script                => catdir( $self->o('enscode_root_dir'), 'ensembl-genes', 'pipelines' , 'gene_symbol_classifier'  ),   
    gst_dump_proteins_script          => catfile( $self->o('ensembl_gst_script'), 'dump_protein_sequences.pl' ),
    gst_load_symbols_script          => catfile( $self->o('ensembl_gst_script'), 'load_gene_symbols.pl' ),	
    registry_status_update_script => catfile( $self->o('ensembl_analysis_script'), 'update_assembly_registry.pl' ),
    core_metadata_script     => catdir( $self->o('enscode_root_dir'), 'ensembl-genes', 'src', 'python', 'ensembl', 'genes', 'metadata', 'core_meta_data.py'),
    core_stats_script        => catdir( $self->o('enscode_root_dir'), 'ensembl-genes', 'src', 'perl', 'ensembl', 'genes', 'generate_species_homepage_stats.pl'),	

########################
# Extra db settings
########################

    'num_tokens'       => 10,

########################
# Executable paths
########################

    samtools_path            => catfile( $self->o('binary_base'),        'samtools' ),                                                    #You may need to specify the full path to the samtools binary

    'uniprot_table_name'          => 'uniprot_sequences',


# Best targetted stuff
    cdna_table_name                             => 'cdna_sequences',


# RNA-seq pipeline stuff
    # You have the choice between:
    #  * using a csv file you already created
    #  * using a study_accession like PRJEB19386
    #  * using the taxon_id of your species
    # 'rnaseq_summary_file' should always be set. If 'taxon_id' or 'study_accession' are not undef
    # they will be used to retrieve the information from ENA and to create the csv file. In this case,
    # 'file_columns' and 'summary_file_delimiter' should not be changed unless you know what you are doing
    'summary_csv_table'      => 'csv_data',
    'read_length_table'      => 'read_length',
    'rnaseq_data_provider'   => 'ENA',           #It will be set during the pipeline or it will use this value

    'rnaseq_dir'  => catdir( $self->o('output_path'), 'rnaseq' ),
    'input_dir'   => catdir( $self->o('rnaseq_dir'),  'input' ),

    'rnaseq_ftp_base' => 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/',

    'rnaseq_summary_file'          => '' || catfile( $self->o('rnaseq_dir'),    $self->o('species_name') . '.csv' ),                                     # Set this if you have a pre-existing cvs file with the expected columns
    'rnaseq_summary_file_genus'    => '' || catfile( $self->o('rnaseq_dir'),    $self->o('species_name') . '_gen.csv' ),                                 # Set this if you have a pre-existing genus level cvs file with the expected columns
    'long_read_dir'       => catdir( $self->o('output_path'),   'long_read' ),
    'long_read_summary_file'       => '' || catfile( $self->o('long_read_dir'), $self->o('species_name') . '_long_read.csv' ),                           # csv file for minimap2, should have 2 columns tab separated cols: sample_name\tfile_name
    'long_read_summary_file_genus' => '' || catfile( $self->o('long_read_dir'), $self->o('species_name') . '_long_read_gen.csv' ),                       # csv file for minimap2, should have 2 columns tab separated cols: sample_name\tfile_name
    'long_read_fastq_dir'          => '' || catdir( $self->o('long_read_dir'), 'input' ),


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
    file_columns      => [ 'SM',     'ID', 'is_paired', 'filename', 'is_mate_1', 'read_length', 'is_13plus', 'CN', 'PL', 'DS' ],
    long_read_columns => [ 'sample', 'filename' ],

########################
# Interproscan
########################
    'realign_table_name'               => 'projection_source_sequences',

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
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# No option below this mark should be modified
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
########################################################
# URLs for retrieving the INSDC contigs and RefSeq files
########################################################
    'ncbi_base_ftp'          => 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all',
    'insdc_base_ftp'         => $self->o('ncbi_base_ftp') . '/#expr(substr(#assembly_accession#, 0, 3))expr#/#expr(substr(#assembly_accession#, 4, 3))expr#/#expr(substr(#assembly_accession#, 7, 3))expr#/#expr(substr(#assembly_accession#, 10, 3))expr#/#assembly_accession#_#assembly_name#',
    'assembly_ftp_path'      => $self->o('insdc_base_ftp'),

########################
# db info
########################
    'pipe_db_server'               => $ENV{GBS4},                                                                                                        # host for pipe db
    'dna_db_server'                => $ENV{GBS2},                                                                                                        # host for dna db
    'pipe_db_port'                 => $ENV{GBP4},                                                                                                        # port for pipeline host
    'dna_db_port'                  => $ENV{GBP2},                                                                                                        # port for dna db host
    'registry_db_server'           => $ENV{GBS1},                                                                                                        # host for registry db
    'registry_db_port'             => $ENV{GBP1},                                                                                                        # port for registry db
    'registry_db_name'             => 'gb_assembly_registry', 

    'core_db' => {
      -dbname => $self->o('dna_db_name'),
      -host   => $self->o('dna_db_server'),
      -port   => $self->o('dna_db_port'),
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


    'pipe_db' => {
      -dbname => $self->o('pipe_db_name'),
      -host   => $self->o('pipe_db_server'),
      -port   => $self->o('pipe_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    #######################
    # Extra db settings
    ########################
    num_tokens => 10,

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
      $tables .= $key . ' SMALLINT UNSIGNED NOT NULL,';
    }
    elsif ( $key eq 'DS' ) {
      $tables .= $key . ' VARCHAR(255) NOT NULL,';
    }
    else {
      $tables .= $key . ' VARCHAR(50) NOT NULL,';
    }
  }
  $tables .= ' KEY(SM), KEY(ID)';


################
# LastZ
################

  my $second_pass = exists $self->{'_is_second_pass'};
  $self->{'_is_second_pass'} = $second_pass;
  return $self->SUPER::pipeline_create_commands if $self->can('no_compara_schema');
  my $pipeline_url = $self->pipeline_url();
  my $parsed_url   = $second_pass && Bio::EnsEMBL::Hive::Utils::URL::parse($pipeline_url);
  my $driver       = $second_pass ? $parsed_url->{'driver'} : '';

################
# /LastZ
################

  return [
    # inheriting database and hive tables' creation
    @{ $self->SUPER::pipeline_create_commands },

    $self->hive_data_table( 'protein', $self->o('uniprot_table_name') ),

    $self->hive_data_table( 'refseq', $self->o('cdna_table_name') ),

    $self->db_cmd( 'CREATE TABLE ' . $self->o('realign_table_name') . ' (' .
        'accession varchar(50) NOT NULL,' .
        'seq text NOT NULL,' .
        'PRIMARY KEY (accession))' ),

    $self->db_cmd( 'CREATE TABLE ' . $self->o('summary_csv_table') . " ($tables)" ),

    $self->db_cmd( 'CREATE TABLE ' . $self->o('read_length_table') . ' (' .
        'fastq varchar(50) NOT NULL,' .
        'read_length int(50) NOT NULL,' .
        'PRIMARY KEY (fastq))' ),

  ];
}


sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{ $self->SUPER::pipeline_wide_parameters },
    wide_ensembl_release => $self->o('ensembl_release'),
    load_toplevel_only => $self->o('load_toplevel_only'),
    skip_braker => 1, # default skip otherfeatures braker
  };
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

  my %genblast_params = (
    wu              => '-P wublast -gff -e #blast_eval# -c #blast_cov#',
    ncbi            => '-P blast -gff -e #blast_eval# -c #blast_cov# -W 3 -softmask -scodon 50 -i 30 -x 10 -n 30 -d 200000 -g T',
    wu_genome       => '-P wublast -gff -e #blast_eval# -c #blast_cov#',
    ncbi_genome     => '-P blast -gff -e #blast_eval# -c #blast_cov# -W 3 -softmask -scodon 50 -i 30 -x 10 -n 30 -d 200000 -g T',
    wu_projection   => '-P wublast -gff -e #blast_eval# -c #blast_cov# -n 100 -x 5 ',
    ncbi_projection => '-P blast -gff -e #blast_eval# -c #blast_cov# -W 3 -scodon 50 -i 30 -x 10 -n 30 -d 200000 -g T',
  );
  my %commandline_params = (
    'ncbi'        => '-num_threads 3 -window_size 40',
    'wu'          => '-cpus 3 -hitdist 40',
    'legacy_ncbi' => '-a 3 -A 40',
  );
  my $header_line = create_header_line( $self->default_options->{'file_columns'} );

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
        'num_threads'                 => $self->o('num_threads'),
        'dbowner'                     => $self->o('dbowner'),
        'core_db'                     => $self->o('core_db'),
        'ensembl_release'             => $self->o('ensembl_release'),
        'base_output_dir'             => $self->o('base_output_dir'),
        'registry_db'                 => $self->o('registry_db'),
        'enscode_root_dir'            => $self->o('enscode_root_dir'),
        'registry_file'               => $self->o('registry_file'),
        'diamond_validation_db'       => $self->o('diamond_validation_db'),
        'validation_type'             => $self->o('validation_type'),
        'use_existing_short_read_dir' => $self->o('use_existing_short_read_dir'),
        'override_clade'              => $self->o('override_clade'),
        'pipe_db'                     => $self->o('pipe_db'),
        'current_genebuild'           => $self->o('current_genebuild'),
	'init_config'     =>$self->o('init_config'),
        'assembly_accession'     =>$self->o('assembly_accession'),
   	'repeatmodeler_library' =>$self->o('repeatmodeler_library'),
   },
      -rc_name => 'default',

      -flow_into => {
        1 => ['update_annotation_tracking_started'],
      },
      -analysis_capacity => 1,
      -input_ids         => [
        #{'assembly_accession' => 'GCA_910591885.1'},
	  ],
    },
    
      {#we need to insert the script or command to update the annotation tracking in the new assembly registry
	  -logic_name => 'update_annotation_tracking_started',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
	  -parameters => {
		cmd => 'echo update command goes here',
	},
            -rc_name => 'default',
	    -flow_into       => { 1 => ['download_rnaseq_csv'], },
	},

    {
      -logic_name => 'download_rnaseq_csv',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name    => '1GB',
      -parameters => {
        cmd => 'python ' . catfile( $self->o('enscode_root_dir'), 'ensembl-genes', 'scripts','transcriptomic_data','get_transcriptomic_data.py' ) . ' -t #species_taxon_id# ' .'-f #rnaseq_summary_file# --read_type short -l 500' ,
        
      },
        -flow_into => {
        1 => ['download_genus_rnaseq_csv'],
      },
    },
  
    {
      -logic_name => 'download_genus_rnaseq_csv',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name    => '1GB',
      -parameters => {
        cmd => 'python ' . catfile( $self->o('enscode_root_dir'), 'ensembl-genes', 'scripts','transcriptomic_data','get_transcriptomic_data.py' ) . ' -t #genus_taxon_id# ' .'-f #rnaseq_summary_file_genus# --read_type short --tree -l 250' ,
       },
      -flow_into => {
        '1->A' => WHEN(
	    # if rnaseq_summary_file has at least 20 lines, i.e. 10 runs (2 read files per run) 
	    '[ $(wc -l < "#rnaseq_summary_file#") -ge 20 ]' => {'fan_short_read_download' => {'inputfile'  => '#rnaseq_summary_file#','input_dir'  => '#short_read_dir#',},},
	    # if rnaseq_summary_file has less than 20 lines BUT rnaseq_summary_file_genus has at least 10 lines, i.e. 5 runs, use that
	    '[ $(wc -l < "#rnaseq_summary_file#") -lt 20 ] && [ $(wc -l < "#rnaseq_summary_file_genus#") -ge 10 ]' => {'fan_short_read_download' => {'inputfile'  => '#rnaseq_summary_file_genus#','input_dir'  => '#short_read_dir#',},}
	    ),
        'A->1' => ['download_long_read_csv'],
      },
    },
      
    {
      -logic_name => 'fan_short_read_download',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd                     => 'if [ -s "#inputfile#" ]; then exit 0; else exit 42;fi',
        return_codes_2_branches => { '42' => 2 },
      },
      -rc_name   => 'default',
      -flow_into => {
	1 => { 'create_sr_fastq_download_jobs' => { 'inputfile' => '#inputfile#', 'input_dir' => '#input_dir#' } },
      },
    },
      
    {
      -logic_name => 'create_sr_fastq_download_jobs',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        column_names => $self->o('file_columns'),
        delimiter    => '\t',
      },
      -flow_into => {
              2 => { 'download_short_read_fastqs' => { 'iid' => '#filename#', 'input_dir' => '#input_dir#' }}, 
         },
    },


    {
      -logic_name => 'download_short_read_fastqs',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadRNASeqFastqs',
      -parameters => {
        ftp_base_url => $self->o('rnaseq_ftp_base'),
        input_dir    => $self->o('input_dir'),
      },
      -analysis_capacity => 50,
      -flow_into => {
        1 => ['download_long_read_csv'],
      },
    },


    {
      -logic_name => 'download_long_read_csv',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name    => '1GB',
      -parameters => {
        cmd => 'python ' . catfile( $self->o('enscode_root_dir'), 'ensembl-genes', 'scripts','transcriptomic_data','get_transcriptomic_data.py' ) . ' -t #species_taxon_id# ' .'-f #long_read_summary_file# --read_type long' ,
      },

      -flow_into => {
        '1->A' => { 'fan_long_read_download' => { 'inputfile' => '#long_read_summary_file#', 'input_dir' => '#long_read_dir#' } },
        'A->1' => ['create_core_db'],
      },
    },

    {
      -logic_name => 'fan_long_read_download',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd                     => 'if [ -s "#inputfile#" ]; then exit 0; else exit 42;fi',
        return_codes_2_branches => { '42' => 2 },
      },
      -flow_into => {
        1 => ['create_lr_fastq_download_jobs'],
      },
      -rc_name => 'default',
    },


    {
      -logic_name => 'create_lr_fastq_download_jobs',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
      -parameters => {
        column_names => $self->o('file_columns'),
        delimiter    => '\t',
      },
      -flow_into => {
        2 => { 'download_long_read_fastq' => { 'iid' => '#filename#', 'input_dir' => '#input_dir#' } },
      },
    },


    {
      -logic_name => 'download_long_read_fastq',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadRNASeqFastqs',
      -parameters => {
        ftp_base_url  => $self->o('rnaseq_ftp_base'),
        input_dir     => $self->o('long_read_fastq_dir'),
        samtools_path => $self->o('samtools_path'),
        decompress    => 1,
        create_faidx  => 1,
      },
      -rc_name           => '1GB',
      -analysis_capacity => 50,
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
      -rc_name => 'default',

      -flow_into => {
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
      -rc_name => 'default',

      -flow_into => {
	      # 1 => ['process_assembly_info'],
	      1 => WHEN ('#load_toplevel_only# == 1' => ['process_assembly_info'],
                        '#load_toplevel_only# == 2' => ['custom_load_toplevel']),
      },
    },
    ####
    # Loading custom assembly where the user provide a FASTA file, probably a repeat library
    ####
    {
      -logic_name => 'custom_load_toplevel',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.catfile($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'assembly_loading', 'load_seq_region.pl').
          ' -dbhost '.$self->o('dna_db_server').
          ' -dbuser '.$self->o('user').
          ' -dbpass '.$self->o('password').
          ' -dbport '.$self->o('dna_db_port').
          ' -dbname '.'#core_dbname#'.
          ' -coord_system_version '.$self->o('assembly_name').
          ' -default_version'.
          ' -coord_system_name primary_assembly'.
          ' -rank 1'.
          ' -fasta_file '. $self->o('custom_toplevel_file_path').
          ' -sequence_level'.
          ' -noverbose',
      },
      -rc_name => '4GB',
      -flow_into => {
        1 => ['custom_set_toplevel'],
      },
    },

    {
      -logic_name => 'custom_set_toplevel',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.catfile($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'assembly_loading', 'set_toplevel.pl').
          ' -dbhost '.$self->o('dna_db_server').
          ' -dbuser '.$self->o('user').
          ' -dbpass '.$self->o('password').
          ' -dbport '.$self->o('dna_db_port').
          ' -dbname '.'#core_dbname#',
      },
      -rc_name => 'default',
      -flow_into  => {
        1 => ['custom_add_meta_keys'],
      },
    },

    {
      -logic_name => 'custom_add_meta_keys',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => '#core_db#',
        sql => [
          'INSERT INTO meta (species_id,meta_key,meta_value) VALUES (1,"assembly.default","'.$self->o('assembly_name').'")',
          'INSERT INTO meta (species_id,meta_key,meta_value) VALUES (1,"assembly.name","'.$self->o('assembly_name').'")',
          'INSERT INTO meta (species_id,meta_key,meta_value) VALUES (1,"species.taxonomy_id","'.$self->o('taxon_id').'")',
        ],
      },
      -rc_name    => 'default',
      -flow_into => {
        1 => ['anno_load_meta_info'],
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
      -rc_name         => '8GB',
      -max_retry_count => 3,
      -flow_into       => {
        1 => ['check_load_meta_info'],
      },
    },

    {
      -logic_name => 'check_load_meta_info',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd                     => 'if [ -s "#rnaseq_summary_file#" ] || [ -s "#long_read_summary_file#" ]; then exit 0; else  exit 42;fi',
        return_codes_2_branches => { '42' => 2 },
      },
      -rc_name => 'default',
      -flow_into => {
        1 => ['anno_load_meta_info'],
        2 => ['helixer_load_meta_info'],
      },
    },
    
    {
      # Load some meta info and seq_region_synonyms
      -logic_name => 'anno_load_meta_info',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => '#core_db#',
        sql     => [
          'DELETE FROM meta WHERE meta_key="species.display_name"',
          'INSERT IGNORE INTO meta (species_id, meta_key, meta_value) VALUES ' .
            '(1, "annotation.provider_name", "Ensembl"),' .
            '(1, "annotation.provider_url", "www.ensembl.org"),' .
            '(1, "assembly.coverage_depth", "high"),' .
            '(1, "assembly.provider_name", NULL),' .
            '(1, "assembly.provider_url", NULL),' .
            '(1, "assembly.ucsc_alias", NULL),' .
            '(1, "species.stable_id_prefix", "#stable_id_prefix#"),' .
            '(1, "species.url", "#species_url#"),' .
            '(1, "species.display_name", "#species_display_name#"),' .
            '(1, "species.division", "#species_division#"),' .
            '(1, "species.strain", "#species_strain#"),' .
            '(1, "species.production_name", "#production_name#"),' .
            '(1, "strain.type", "#strain_type#"),' .
            '(1, "repeat.analysis", "repeatdetector"),' .
            '(1, "repeat.analysis", "dust"),' .
            '(1, "repeat.analysis", "trf"),' .
            '(1, "genebuild.initial_release_date", NULL),' .
            '(1, "genebuild.id", ' . $self->o('genebuilder_id') . '),' .
            '(1, "genebuild.method", "anno"),'.
	    '(1, "genebuild.method_display", "Ensembl Genebuild"),'.
        '(1, "species.annotation_source", "ensembl")'
        ],
      },
      -max_retry_count => 0,
      -rc_name         => 'default',
      -flow_into       => {
        1 => ['load_taxonomy_info'],
      },
    },
    {
      # Load some meta info and seq_region_synonyms
      -logic_name => 'helixer_load_meta_info',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => '#core_db#',
        sql     => [
          'DELETE FROM meta WHERE meta_key="species.display_name"',
          'INSERT IGNORE INTO meta (species_id, meta_key, meta_value) VALUES ' .
            '(1, "annotation.provider_name", "Ensembl"),' .
            '(1, "annotation.provider_url", "www.ensembl.org"),' .
            '(1, "assembly.coverage_depth", "high"),' .
            '(1, "assembly.provider_name", NULL),' .
            '(1, "assembly.provider_url", NULL),' .
            '(1, "assembly.ucsc_alias", NULL),' .
            '(1, "species.stable_id_prefix", "HELIXER#species_prefix#"),' .
            '(1, "species.url", "#species_url#"),' .
            '(1, "species.display_name", "#species_display_name#"),' .
            '(1, "species.division", "#species_division#"),' .
            '(1, "species.strain", "#species_strain#"),' .
            '(1, "species.production_name", "#production_name#"),' .
            '(1, "strain.type", "#strain_type#"),' .
            '(1, "repeat.analysis", "repeatdetector"),' .
            '(1, "repeat.analysis", "dust"),' .
            '(1, "repeat.analysis", "trf"),' .
            '(1, "genebuild.initial_release_date", NULL),' .
            '(1, "genebuild.id", ' . $self->o('genebuilder_id') . '),' .
            '(1, "genebuild.method", "helixer"),'.
	    '(1, "genebuild.method_display", "HELIXER"),'.
	    '(1, "species.annotation_source", "helixer")'
        ],
      },
      -max_retry_count => 0,
      -rc_name         => 'default',
      -flow_into       => {
        1 => ['load_taxonomy_info'],
      },
    },

    {
      -logic_name => 'load_taxonomy_info',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadTaxonomyInfo',
      -parameters => {
        'target_db'   => '#core_db#',
        'taxonomy_db' => $self->o('taxonomy_db'),
      },
      -rc_name => 'default',

      -flow_into => {
        1 => ['dump_toplevel_file'],    #['load_windowmasker_repeats'],# 'fan_refseq_import'],
      },
    },


    {
      -logic_name => 'dump_toplevel_file',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDumpGenome',
      -parameters => {
        'coord_system_name'  => 'toplevel',
        'target_db'          => '#core_db#',
        'output_path'        => '#output_path#',
        'enscode_root_dir'   => $self->o('enscode_root_dir'),
        'species_name'       => '#species_name#',
        'repeat_logic_names' => $self->o('softmask_logic_names'),    # This is emtpy as we just use masking present in downloaded file
      },
      -flow_into => {
        1 => ['reheader_toplevel_file'],
      },
      -rc_name => '3GB',
    },


    {
      -logic_name => 'reheader_toplevel_file',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . catfile( $self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'genebuild', 'convert_genome_dump.pl' ) .
          ' -conversion_type slice_name_to_seq_region_name' .
          ' -input_file #toplevel_genome_file#' .
          ' -output_file #reheadered_toplevel_genome_file#',
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['check_transcriptomic_data'],
      },
    },
    {
      -logic_name => 'check_transcriptomic_data',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd                     => 'if [ -s "#rnaseq_summary_file#" ] ||[ -s "#long_read_summary_file#" ]; then exit 0; else exit 42;fi',
        return_codes_2_branches => { '42' => 2 },
      },
      -flow_into => {
        1 => ['run_anno'],
        2 => ['run_anno_softmasking'],
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'run_anno',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
         cmd => 'python ' . catfile( $self->o('enscode_root_dir'), 'ensembl-anno', 'ensembl_anno.py' ) . ' #anno_commandline#',
      },
      -rc_name         => 'anno',
      -max_retry_count => 0,
      -flow_into       => {
        1 => ['update_biotypes_and_analyses']
      },
    },
    {
      -logic_name => 'run_anno_softmasking',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'python ' . catfile( $self->o('enscode_root_dir'), 'ensembl-anno', 'ensembl_anno.py' ) . ' #anno_red_commandline#;' .
          'cp #output_path#/red_output/mask_output/#species_name#_reheadered_toplevel.msk #output_path#/#species_name#_softmasked_toplevel.fa',
      },
      -rc_name         => 'anno',
      -max_retry_count => 0,
      -flow_into       => {
        1 => ['run_helixer'],
      },
    },
    {
      -logic_name => 'run_helixer',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
	cmd => 'mkdir -p #output_path#/helixer;' .
	  'singularity exec --nv '.  $self->o ('helixer_singularity_image') .  ' bash -c "'. 
	  'export PATH=\$PATH:/hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/HelixerPost/target/release/ && '. 
	  'Helixer.py --fasta-path #reheadered_toplevel_genome_file# --lineage fungi --gff-output-path #output_path#/helixer/#assembly_accession#_#species_name#.gff3";' .
	  $self->o('gffread_path').' #output_path#/helixer/#assembly_accession#_#species_name#.gff3 -T -o #output_path#/helixer/helixer.gtf;' .
          $self->o ('gffread_path').' #output_path#/helixer/#assembly_accession#_#species_name#.gff3 -g #reheadered_toplevel_genome_file# --adj-stop #output_path#/helixer/#assembly_accession#_#species_name#.gff3 -y #output_path#/helixer/helixer_proteins.fa;',  
      },
      -rc_name         => 'helixer',
      -max_retry_count => 0,
      -flow_into       => {
        1 => ['load_gtf_file'],
      },
    },
    {		   
      -logic_name => 'load_gtf_file',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . catfile( $self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'genebuild', 'helixer', 'parse_gtf.pl' ) .
          ' -dnahost ' . $self->o('dna_db_server') .
          ' -dnauser ' . $self->o('user_r') .
          ' -dnaport ' . $self->o('dna_db_port') .
          ' -dnadbname #core_dbname#' .
          ' -host ' . $self->o('dna_db_server') .
          ' -user ' . $self->o('user') .
          ' -pass ' . $self->o('password') .
          ' -port ' . $self->o('dna_db_port') .
          ' -dbname #core_dbname#' .
          ' -write' .
          ' -file #output_path#/helixer/helixer.gtf',
      },
      -rc_name         => 'default',
      -max_retry_count => 0,
      -flow_into       => {
        1 => ['update_biotypes_and_analyses'],
      },
    },
    {
      -logic_name => 'update_biotypes_and_analyses',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => '#core_db#',
        sql     => [
          'UPDATE analysis SET module=NULL',
          'UPDATE gene SET biotype = "protein_coding" WHERE biotype = "anno_protein_coding"',
          'UPDATE gene SET biotype = "lncRNA" WHERE biotype = "anno_lncRNA"',
          'UPDATE transcript JOIN gene USING(gene_id) SET transcript.biotype = gene.biotype',
          'UPDATE transcript JOIN gene USING(gene_id) SET transcript.analysis_id = gene.analysis_id',
          'UPDATE repeat_feature SET repeat_start = 1 WHERE repeat_start < 1',
          'UPDATE repeat_feature SET repeat_end = 1 WHERE repeat_end < 1',
          'UPDATE repeat_feature JOIN seq_region USING(seq_region_id) SET seq_region_end = length WHERE seq_region_end > length',
          'UPDATE gene SET display_xref_id=NULL',
          'UPDATE transcript SET display_xref_id=NULL',
        ],
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['delete_duplicate_genes'],
      },
    },

    {
      -logic_name => 'delete_duplicate_genes',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . $self->o('remove_duplicates_script') .
          ' -dbuser ' . $self->o('user') .
          ' -dbpass ' . $self->o('password') .
          ' -dbhost ' . $self->o( 'core_db', '-host' ) .
          ' -dbport ' . $self->o( 'core_db', '-port' ) .
          ' -dbname ' . '#core_dbname#'
      },
      -rc_name   => '5GB',
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
          ' -host ' . $self->o( 'core_db', '-host' ) .
          ' -port ' . $self->o( 'core_db', '-port' ) .
          ' -dbpattern ' . '#core_dbname#'
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
          ' -host ' . $self->o( 'core_db', '-host' ) .
          ' -port ' . $self->o( 'core_db', '-port' ) .
          ' -dbname ' . '#core_dbname#'
      },
      -rc_name   => 'default',
      -flow_into => { 1 => ['set_frameshift_introns'] },
    },

    {
      -logic_name => 'set_frameshift_introns',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . $self->o('frameshift_attrib_script') .
          ' -user ' . $self->o('user') .
          ' -pass ' . $self->o('password') .
          ' -host ' . $self->o( 'core_db', '-host' ) .
          ' -port ' . $self->o( 'core_db', '-port' ) .
          ' -dbpattern ' . '#core_dbname#'
      },
      -rc_name   => '10GB',
      -flow_into => { 1 => ['set_canonical_transcripts'] },
    },

    {
      -logic_name => 'set_canonical_transcripts',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . $self->o('select_canonical_script') .
          ' -dbuser ' . $self->o('user') .
          ' -dbpass ' . $self->o('password') .
          ' -dbhost ' . $self->o( 'core_db', '-host' ) .
          ' -dbport ' . $self->o( 'core_db', '-port' ) .
          ' -dbname ' . '#core_dbname#' .
          ' -coord toplevel -write'
      },
      -rc_name   => '10GB',
      -flow_into => { 1 => ['null_columns'] },
    },

    {
      -logic_name => 'null_columns',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => '#core_db#',
        sql     => [
          'UPDATE gene SET stable_id = NULL',
          'UPDATE transcript SET stable_id = NULL',
          'UPDATE translation SET stable_id = NULL',
          'UPDATE exon SET stable_id = NULL',
          'UPDATE protein_align_feature set external_db_id = NULL',
          'UPDATE dna_align_feature set external_db_id = NULL',
        ],
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['check_run_stable_ids'],
      },
    },

    {
      -logic_name => 'check_run_stable_ids',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd                     => 'if [ -s "#rnaseq_summary_file#" ] || [ -s "#long_read_summary_file#" ]; then exit 0; else exit 42;fi',
        return_codes_2_branches => { '42' => 2 },
      },
      -flow_into => {
        1 => ['anno_run_stable_ids'],
        2 => ['helixer_run_stable_ids'],
      },
      -rc_name => 'default',
    },
    {
      -logic_name => 'anno_run_stable_ids',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::SetStableIDs',
      -parameters => {
        enscode_root_dir => $self->o('enscode_root_dir'),
        mapping_required => 0,
        target_db        => '#core_db#',
        id_start         => '#stable_id_prefix#' . '#stable_id_start#',
        output_path      => '#output_path#',
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['final_meta_updates'],
      },
    },
    {
      -logic_name => 'helixer_run_stable_ids',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::SetStableIDs',
      -parameters => {
        enscode_root_dir => $self->o('enscode_root_dir'),
        mapping_required => 0,
        target_db        => '#core_db#',
        id_start         => 'HELIXER#species_prefix#' . '#stable_id_start#',
        output_path      => '#output_path#',
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['final_meta_updates'],
      },
    },

    {
      -logic_name => 'final_meta_updates',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => '#core_db#',
        sql     => [
          'INSERT IGNORE INTO meta (species_id, meta_key, meta_value) VALUES ' .
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
        db_conn => '#core_db#',
        sql     => [
          'TRUNCATE associated_xref',
          'TRUNCATE dependent_xref',
          'TRUNCATE identity_xref',
          'TRUNCATE object_xref',
          'TRUNCATE ontology_xref',
          'TRUNCATE xref',
	  'DELETE from meta where meta_key="species.strain_group"',
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
        input_db => '#core_db#',
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
        species => '#production_name#',
        group   => 'core',
      },
      -flow_into => {
        1 => ['run_busco_core_genome_mode'],
      },
      -rc_name => 'default_registry',
    },
    {
      -logic_name => 'run_busco_core_genome_mode',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'cd #output_path#; ' .
          'singularity exec ' . $self->o('busco_singularity_image') . ' busco -f -i #output_path#/#species_name#_reheadered_toplevel.fa  -m genome -l #busco_group# -c ' . $self->o('cores') . ' -o busco_core_genome_mode_output --offline --download_path ' . $self->o('busco_download_path') . ' ; ' .
          'rm -rf  #output_path#/busco_core_genome_mode_output/logs;' .
          'rm -rf  #output_path#/busco_core_genome_mode_output/busco_downloads;' .
          'rm -rf  #output_path#/busco_core_genome_mode_output/run*;' .
          'sed  -i "/genebuild/d"  #output_path#/busco_core_genome_mode_output/*.txt;' .
          'mv #output_path#/busco_core_genome_mode_output/*.txt #output_path#/busco_core_genome_mode_output/#species_strain_group#_genome_busco_short_summary.txt',
      },
      -rc_name   => '32GB',
      -flow_into => { 1 => ['fan_busco_output'] },
    },
    
     {
      -logic_name => 'fan_busco_output',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd                     => 'if [ -s "#rnaseq_summary_file#" ] || [ -s "#long_read_summary_file#" ]; then exit 0; else  exit 42;fi',
        return_codes_2_branches => { '42' => 2 },
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['create_busco_dirs'],
        2 => ['run_helixer_busco'],
      },
    },

   {
      -logic_name => 'run_helixer_busco',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',

      -parameters => {
        cmd => 'cd #output_path#/;' .
          'singularity exec ' . $self->o('busco_singularity_image') . ' busco -f -i #output_path#/helixer/helixer_proteins.fa  -m prot -l #busco_group# -c ' . $self->o('cores') . ' -o busco_core_protein_mode_output --offline --download_path ' . $self->o('busco_download_path') . ' ; ' .
	  'rm -rf  #output_path#/busco_core_protein_mode_output/logs;' .
	  'rm -rf  #output_path#/busco_core_protein_mode_output/busco_downloads;' .
	  'rm -rf  #output_path#/busco_core_protein_mode_output/run*;' .
	  'sed  -i "/genebuild/d"  #output_path#/busco_core_protein_mode_output/*.txt;' .
	  'mv #output_path#/busco_core_protein_mode_output/*.txt #output_path#/busco_core_protein_mode_output/#species_strain_group#_busco_short_summary.txt;',
      },
      -rc_name   => '32GB',
      -flow_into => {
        1 => ['gst_dump_protein_sequences'],
      },
    },


    {
      -logic_name => 'create_busco_dirs',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'mkdir -p #output_path#' . '/busco_score_data',
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['dump_canonical_stable_ids'],
      },
    },

    {
      -logic_name => 'dump_canonical_stable_ids',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
      -parameters => {
        db_conn     => '#core_db#',
        input_query => 'SELECT transcript.stable_id from gene, transcript ' .
          ' WHERE gene.gene_id = transcript.gene_id ' .
          ' AND gene.canonical_transcript_id = transcript.transcript_id ' .
          ' AND transcript.biotype = "protein_coding" ',
        command_out        => q( grep 'ENS' > #busco_input_file_m#),
        busco_input_file_m => catfile( '#output_path#', '/busco_score_data/', $self->o('busco_input_file_stid') ),
        prepend            => [ '-NB', '-q' ],
      },
      -rc_name   => '2GB',
      -flow_into => {
        1 => ['print_translations'],
      },
    },

    {
      -logic_name => 'print_translations',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . $self->o('print_protein_script_path') .
          ' --user=' . $self->o('user_r') .
          ' --host=' . $self->o( 'core_db', '-host' ) .
          ' --port=' . $self->o( 'core_db', '-port' ) .
          ' --dbname=' . '#core_dbname#' .
          ' --id_file=' . catfile( '#output_path#', '/busco_score_data/', $self->o('busco_input_file_stid') ) .
          ' --output_file=' . catfile( '#output_path#', '/busco_score_data/', 'canonical_proteins.fa' ),
      },
      -rc_name   => 'default',
      -flow_into => { 1 => ['run_busco_anno'] },
    },

    {
      -logic_name => 'run_busco_anno',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'cd #output_path#; ' .
	  'singularity exec ' . $self->o('busco_singularity_image') . ' busco -f -i #output_path#/busco_score_data/canonical_proteins.fa  -m prot -l #busco_group# -c ' . $self->o('cores') . ' -o busco_core_protein_mode_output --offline --download_path ' . $self->o('busco_download_path') . ' ; ' .
	  'rm -rf  #output_path#/busco_core_protein_mode_output/logs;' .
          'rm -rf  #output_path#/busco_core_protein_mode_output/busco_downloads;' .
          'rm -rf  #output_path#/busco_core_protein_mode_output/run*;' .
          'sed  -i "/genebuild/d"  #output_path#/busco_core_protein_mode_output/*.txt;' .
	  'mv #output_path#/busco_core_protein_mode_output/*.txt #output_path#/busco_core_protein_mode_output/#species_strain_group#_busco_short_summary.txt',
      },
      -rc_name   => '32GB',
      -flow_into => { 1 => ['gst_dump_protein_sequences'] },
    },
    
     {
      -logic_name => 'gst_dump_protein_sequences',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
	  cmd => 'perl ' . $self->o('gst_dump_proteins_script') . ' --group core --species #production_name# --registry #registry_file# --output_file #gst_dir#/#production_name#_protein_sequences.fa',
      },
	  -rc_name => 'default_registry',
	  -flow_into       => { 1 => ['gst_assign_gene_symbols'], },

    },

    {
      -logic_name => 'gst_assign_gene_symbols',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
	  cmd => 'singularity run --bind /hps/software/users/ensembl/genebuild/gene_symbol_classifier/data:/app/checkpoints --bind #gst_dir#:/app/data /hps/software/users/ensembl/genebuild/gene_symbol_classifier/singularity/gene_symbol_classifier_0.12.1.sif --checkpoint /app/checkpoints/mlp_10_min_frequency_2022-01-29_03.08.32.ckpt --sequences_fasta /app/data/#production_name#_protein_sequences.fa --scientific_name #species_name#',
      },
	  -rc_name => 'default_registry',
	  -flow_into       => { 1 => ['gst_filter_assignments'], },
    },

    {
     -logic_name => 'gst_filter_assignments',
     -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
     -parameters => {
	 cmd => 'singularity run --bind #gst_dir#:/app/data /hps/software/users/ensembl/genebuild/gene_symbol_classifier/singularity/gene_symbol_classifier_filter_0.3.0.sif --symbol_assignments /app/data/#production_name#_protein_sequences_symbols.csv --threshold 0.7',
     },
	 -rc_name => 'default_registry',
	 -flow_into       => { 1 => ['gst_load_gene_symbols'], },	 
    },

    {
     -logic_name => 'gst_load_gene_symbols',
     -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
     -parameters => {
	 cmd => 'perl ' . $self->o('gst_load_symbols_script') . ' --species #production_name# --group core --registry #registry_file# --symbol_assignments #gst_dir#/#production_name#_protein_sequences_symbols_filtered.csv --primary_ids_file /hps/software/users/ensembl/genebuild/gene_symbol_classifier/data/display_name_dbprimary_acc_105.dat --program_version 0.12.1',
     },
	 -rc_name => 'default_registry',
	 -flow_into       => { 1 => ['run_meta_updates'], },
    },
{
      -logic_name => 'run_meta_updates',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
          cmd => 'python ' . $self->o('core_metadata_script')  .  ' -o #output_path# -d #core_dbname# --team genebuild -s ' . $self->o('dna_db_server') . ' -p ' .$self->o('dna_db_port'),
      },
      -rc_name => '1GB',
      -flow_into       => { 1 => ['load_meta_updates'], },
    },

    {
      -logic_name => 'load_meta_updates',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
          cmd => '/hps/software/users/ensembl/ensw/mysql-cmds/ensembl/ensadmin/' . $self->o('dna_db_server') . ' #core_dbname# <#output_path#/#core_dbname#.sql',
      },
      -rc_name => 'default_registry',
      -flow_into       => { 1 => ['run_core_stats'], },
    },

    {
      -logic_name => 'run_core_stats',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
          cmd => 'perl ' . $self->o('core_stats_script')  .  ' -dbname #core_dbname# -host ' .  $self->o('dna_db_server') . ' -port ' .$self->o('dna_db_port') . ' -production_name #production_name# -output_dir #output_path#', 
      },
      -rc_name => '5GB',
      -flow_into       => { 1 => ['load_core_stats'], },
    },

    {
      -logic_name => 'load_core_stats',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
          cmd => '/hps/software/users/ensembl/ensw/mysql-cmds/ensembl/ensadmin/' . $self->o('dna_db_server') . ' #core_dbname# <#output_path#/stats_#core_dbname#.sql',
      },
      -rc_name => 'default_registry',
      -flow_into       => { 1 => ['pepstats'], },
    },

    {
	  -logic_name => 'pepstats',
	      -module     => 'Bio::EnsEMBL::Production::Pipeline::Production::PepStatsBatch',
	      -parameters => {
		  dbtype => 'core',
		  species => '#production_name#',
		  pepstats_binary => 'pepstats',
		  tmpdir => '#output_path#',
		  reg_conf => '#registry_file#',
	  },
	      -max_retry_count => 1,
	      -hive_capacity   => 50,
	      -rc_name => '50GB',
	      -flow_into       => { 1 => ['load_genome_busco_into_core'], }
      },

    {
        -logic_name => 'load_genome_busco_into_core',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
          cmd => 'python ' . catfile( $self->o('enscode_root_dir'), 'ensembl-genes','src','python','ensembl','genes','metrics','busco_metakeys_patch.py' ) .
           ' -db #core_dbname# -host ' .  $self->o('dna_db_server') . 
           ' -port ' . $self->o('dna_db_port') . 
           ' -user ' . $self->o('user') . 
           ' -password ' . $self->o('password') . 
           ' -assembly_id #assembly_id#' .
           ' -file #output_path#/busco_core_genome_mode_output/#species_strain_group#_genome_busco_short_summary.txt ' .
           ' -output_dir #output_path#/busco_core_genome_mode_output/ -run_query true',
    
    },
         -rc_name => 'default',
         -flow_into       => { 1 => ['load_protein_busco_into_core'], },
    }, 
    {
        -logic_name => 'load_protein_busco_into_core',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
          cmd => 'python ' .  catfile( $self->o('enscode_root_dir'), 'ensembl-genes','src','python','ensembl','genes','metrics','busco_metakeys_patch.py' ) .
           ' -db #core_dbname# -host ' .  $self->o('dna_db_server') .
           ' -port ' . $self->o('dna_db_port') .
           ' -user ' . $self->o('user') .
           ' -password ' . $self->o('password') .
           ' -assembly_id #assembly_id#' .
           ' -file #output_path#/busco_core_protein_mode_output/#species_strain_group#_busco_short_summary.txt ' .
           ' -output_dir #output_path#/busco_core_protein_mode_output/ -run_query true',

    },
         -rc_name => 'default',
         -flow_into       => { 1 => ['check_busco_score'], },
    },
    {
        -logic_name => 'check_busco_score',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
          cmd => 'if python ' .  catfile( $self->o('enscode_root_dir'), 'ensembl-genes','src','python','ensembl','genes','metrics', 'check_busco_score.py' ) .
          ' --genome #output_path#/busco_core_genome_mode_output/#core_dbname#_busco_genome_metakey.json ' .
          ' --protein #output_path#/busco_core_protein_mode_output/#core_dbname#_busco_protein_metakey.json' . '; then exit 0; else exit 42; fi',
        return_codes_2_branches => { '42' => 2 },
    },
         -rc_name => 'default',
         -flow_into  => {
              1 => 'backbone_job_pipeline',
              0 => 'update_registry_as_check',
    }
    },
  {
    -logic_name => 'update_assembly_registry_status',
    -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -parameters => {
	cmd => 'perl ' . $self->o('registry_status_update_script') .
	    ' -user ' . $self->o('user') .
	    ' -pass ' . $self->o('password') .
	    ' -assembly_accession ' . '#assembly_accession#' .
	    ' -registry_host ' . $self->o('registry_db_server') .
	    ' -registry_port ' . $self->o('registry_db_port') .
	    ' -registry_db ' . $self->o('registry_db_name'),
    },
	-rc_name => 'default',
	-flow_into       => { 1 => ['update_annotation_tracking_complete'], },
  },

  {
    #we need to insert the script or command to update the annotation tracking in the new assembly registry
     # Update status to COMPLETE to indicate that the pipeline made it this far. Should be quickly updated 
     # to something like 'to be checked' or 'pre-released'
  -logic_name => 'update_annotation_tracking_complete',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
    cmd => 'echo update command goes here',
        },
      -rc_name => 'default',
      -flow_into       => { 1 => ['delete_short_reads'], },
    },
    
    {
    -logic_name => 'delete_short_reads',
    -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -parameters => {
      cmd => 'if [ -f ' . '#short_read_dir#' . '/*.gz ]; then rm ' . '#short_read_dir#' . '/*.gz; fi',
    },
    -rc_name => 'default',
    -flow_into       => { 1 => ['delete_long_reads'], },
    },
    {
    -logic_name => 'delete_long_reads',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'if [ -f ' . '#long_read_dir#' . '/* ]; then rm ' . '#long_read_dir#' . '/*; fi',
      },
      -rc_name => 'default',
      -flow_into       => { 1 => ['busco_check'], },
    },
    {
    #we need to insert the script or command to update the annotation tracking in the new assembly registry
    # Update status to COMPLETE to indicate that the pipeline made it this far. Should be quickly updated 
    # to something like 'to be checked' or 'pre-released'
    -logic_name => 'busco_check',
    -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -parameters => {
      cmd => 'echo ',
      threshold => $self->o('busco_threshold'),
    },
    -rc_name    => 'default',
    -flow_into  => {
      # maybe we need to parse this with a different analysis and for this to be a dummy only working on this
      1 => WHEN(
        '#some_parameter# >= #threshold#' => [ 'backbone_job_pipeline' ],
        '#some_parameter# < #threshold#'  => [ 'update_registry_as_check' ],
      ), # furthermore, we can make it so that the script "fails" if the busco check fails, and then redirect channel "-2" to update_registry_as_check, and normal channel "1" to backbone_job_pipeline
    },
  },
  {
    -logic_name => 'update_registry_as_check',
    -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    -parameters => {
      cmd => 'echo update command goes here',
    },
    -rc_name => 'default',
  },
  {
    -logic_name     => 'backbone_job_pipeline',
    -module         => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
    -hive_capacity  => -1,
    -flow_into      => {
      '1->A'  => ['gff3','gtf','softmasked_genome_copy'],
      'A->1'  => ['checksum_generator'],
    }
  },
  {
    -logic_name      => 'checksum_generator',
    -module        => 'Bio::EnsEMBL::Production::Pipeline::Common::ChksumGenerator',
    -parameters    => {
      dumps              => ['gff3','gtf','softmasked_genome_copy'],
      # skip_convert_fasta => $self->o('skip_convert_fasta')
    },
    -hive_capacity => 10,
    -flow_into => {
      '1' => ['sync'],
    },
  },
  {
    -logic_name => 'softmasked_genome_copy',
      -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters    => { cmd => 'mv #out_file#.sorted.gz #out_file#', },
      -hive_capacity => 10,
  },
  ### GTF
  { -logic_name      => 'gtf',
      -module        => 'Bio::EnsEMBL::Production::Pipeline::GTF::DumpFile',
      -parameters    => {
          gtf_to_genepred => $self->o('gtftogenepred_exe'),
          gene_pred_check => $self->o('genepredcheck_exe'),
          abinitio        => $self->o('abinitio'),
          gene            => $self->o('gene')
      },
      -hive_capacity => 50,
      -rc_name       => '2GB',
      -flow_into     => { '-1' => 'gtf_32GB', '1' => 'move_gtf'},
  },

  { -logic_name      => 'gtf_32GB',
      -module        => 'Bio::EnsEMBL::Production::Pipeline::GTF::DumpFile',
      -parameters    => {
          gtf_to_genepred => $self->o('gtftogenepred_exe'),
          gene_pred_check => $self->o('genepredcheck_exe'),
          abinitio        => $self->o('abinitio'),
          gene            => $self->o('gene')
      },
      -hive_capacity => 50,
      -rc_name       => '32GB',
      -flow_into     => 'move_gtf',
  },

  {
      -logic_name    => 'move_gtf',
      -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters    => { cmd => 'mv #out_file#.sorted.gz #out_file#', },
      -hive_capacity => 10,
  },


  ### GFF3
  { -logic_name      => 'gff3',
      -module        => 'Bio::EnsEMBL::Production::Pipeline::GFF3::DumpFile',
      -parameters    => {
          feature_type     => $self->o('feature_type'),
          per_chromosome   => $self->o('per_chromosome'),
          include_scaffold => $self->o('include_scaffold'),
          logic_name       => $self->o('logic_name'),
          db_type          => $self->o('db_type'),
          abinitio         => $self->o('abinitio'),
          gene             => $self->o('gene'),
          out_file_stem    => $self->o('out_file_stem'),
          xrefs            => $self->o('xrefs'),
          base_path        => $self->o('output_path'),
          species          => $self->o('production_name'),
          release          => $self->o('release_number'),
      },
      -hive_capacity => 50,
      -rc_name       => '2GB',
      -flow_into     => { '-1' => 'gff3_32GB', '1' => 'tidy_gff3', },
  },

  { -logic_name      => 'gff3_32GB',
      -module        => 'Bio::EnsEMBL::Production::Pipeline::GFF3::DumpFile',
      -parameters    => {
          feature_type     => $self->o('feature_type'),
          per_chromosome   => $self->o('per_chromosome'),
          include_scaffold => $self->o('include_scaffold'),
          logic_name       => $self->o('logic_name'),
          db_type          => $self->o('db_type'),
          abinitio         => $self->o('abinitio'),
          gene             => $self->o('gene'),
          out_file_stem    => $self->o('out_file_stem'),
          xrefs            => $self->o('xrefs'),
      },
      -hive_capacity => 50,
      -rc_name       => '32GB',
      -flow_into     => { '1' => 'tidy_gff3', },
  },

  ### GFF3:post-processing
  { -logic_name      => 'tidy_gff3',
      -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters    => { cmd => $self->o('gff3_tidy') . ' -gzip -o #out_file#.sorted.gz #out_file#', },
      -hive_capacity => 10,
      -batch_size    => 10,
      -flow_into     => 'move_gff3',
  },

  {
      -logic_name    => 'move_gff3',
      -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters    => { cmd => 'mv #out_file#.sorted.gz #out_file#', },
      -hive_capacity => 10,
      -flow_into     => 'validate_gff3',
  },

  {
      -logic_name    => 'validate_gff3',
      -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters    => { cmd => $self->o('gff3_validate') . ' #out_file#', },
      -hive_capacity => 10,
      -batch_size    => 10,
  },

  {
      -logic_name    => 'sync',
      -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters    => {
          cmd => 'echo "Syncing files to the output directory with rsync"',
      },
      -flow_into => {
          1 => ['update_registry_pre_release'],
      },
  },
  {
      -logic_name    => 'update_registry_pre_release',
      -module        => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters    => {
        cmd => 'echo update registry command goes here',
      }
  }
  ];
}

sub resource_classes {
  my $self = shift;

  return {
    #inherit other stuff from the base class
     %{ $self->SUPER::resource_classes() },
     'anno'             => {
     SLURM =>  $self->slurm_resource_builder(50000, '7-00:00:00', $self->default_options->{'num_threads'} ),
     },
     
     'helixer'       => {
     SLURM =>  $self->slurm_resource_builder(100000, '7-00:00:00',undef, $self->default_options->{'gpu'} ),
     },
     
     '32GB'           => {
     SLURM =>  $self->slurm_resource_builder(32000, '2-00:00:00',  $self->default_options->{'cores'} ),
    },
    };
    }

sub hive_capacity_classes {
  my $self = shift;

  return {
    'hc_low'    => 200,
    'hc_medium' => 500,
    'hc_high'   => 1000,
  };
}

sub check_file_in_ensembl {
  my ( $self, $file_path ) = @_;
  push @{ $self->{'_ensembl_file_paths'} }, $file_path;
  return $self->o('enscode_root_dir') . '/' . $file_path;
}

1;
