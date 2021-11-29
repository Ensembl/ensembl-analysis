
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

package braker_subpipeline;

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

    #BRAKER
    'augustus_config_path'  => '/nfs/production/flicek/ensembl/genebuild/ftricomi/augustus_config/config',
    'augustus_species_path' => '/nfs/production/flicek/ensembl/genebuild/ftricomi/augustus_config/config/species/',
    'cores'                 => 30,
    'genome_file'           => '',
    'braker_singularity_image' => '/hps/software/users/ensembl/genebuild/ftricomi/singularity/test-braker2_es_ep_etp.simg',
    'python_singularity_image' => '/hps/software/users/ensembl/genebuild/ftricomi/singularity/test_clean_gtf.sif',
    'agat_singularity_image'   => '/hps/software/users/ensembl/genebuild/ftricomi/singularity/test-agat.simg',
    'busco_dataset'            => 'lepidoptera_odb10',
    'busco_singularity_image'  => '/hps/software/users/ensembl/genebuild/ftricomi/singularity/busco-v5.1.2_cv1.simg',

    #Gbiab
    'num_threads' => 20,
    'base_output_dir'      => '',
    'protein_file'         => '',
    'busco_protein_file'   => '',
    'rfam_accessions_file' => '',
    'registry_file'        => 'Databases.pm',
    'release_number'       => '001' || $self->o('ensembl_release'),

    'dbowner'                      => '' || $ENV{EHIVE_USER} || $ENV{USER},
    'pipeline_name'                => 'braker_prod_test' || $self->o('production_name') . $self->o('production_name_modifier') . '_' . $self->o('ensembl_release'),
    'user_r'                       => 'ensro',                                                                                        # read only db user
    'user'                         => 'ensadmin',                                                                                     # write db user
    'password'                     => '',                                                                                             # password for write db user
    'server_set'                   => '',                                                                                             # What server set to user, e.g. set1
    'pipe_db_server'               => $ENV{GBS7},                                                                                     # host for pipe db
    'databases_server'             => $ENV{GBS5},                                                                                     # host for general output dbs
    'dna_db_server'                => $ENV{GBS6},                                                                                     # host for dna db
    'pipe_db_port'                 => $ENV{GBP7},                                                                                     # port for pipeline host
    'databases_port'               => $ENV{GBP5},                                                                                     # port for general output db host
    'dna_db_port'                  => $ENV{GBP6},                                                                                     # port for dna db host
    'registry_db_server'           => $ENV{GBS1},                                                                                     # host for registry db
    'registry_db_port'             => $ENV{GBP1},                                                                                     # port for registry db
    'registry_db_name'             => 'gb_assembly_registry',                                                                         # name for registry db
    'repbase_logic_name'           => '',                                                                                             # repbase logic name i.e. repeatmask_repbase_XXXX, ONLY FILL THE XXXX BIT HERE!!! e.g primates
    'repbase_library'              => '',                                                                                             # repbase library name, this is the actual repeat repbase library to use, e.g. "Mus musculus"
    'rnaseq_summary_file'          => '' || catfile( $self->o('rnaseq_dir'),    $self->o('species_name') . '.csv' ),                  # Set this if you have a pre-existing cvs file with the expected columns
    'rnaseq_summary_file_genus'    => '' || catfile( $self->o('rnaseq_dir'),    $self->o('species_name') . '_gen.csv' ),              # Set this if you have a pre-existing genus level cvs file with the expected columns
    'long_read_summary_file'       => '' || catfile( $self->o('long_read_dir'), $self->o('species_name') . '_long_read.csv' ),        # csv file for minimap2, should have 2 columns tab separated cols: sample_name\tfile_name
    'long_read_summary_file_genus' => '' || catfile( $self->o('long_read_dir'), $self->o('species_name') . '_long_read_gen.csv' ),    # csv file for minimap2, should have 2 columns tab separated cols: sample_name\tfile_name
    'long_read_fastq_dir'          => '' || catdir( $self->o('long_read_dir'), 'input' ),
    'species_name'                 => '',                                                                                             # e.g. mus_musculus
    'production_name'              => 'test_scale_braker' || $self->o('species_name'),                                                # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
    'taxon_id'                     => '',                                                                                             # should be in the assembly report file
    'species_taxon_id'             => '' || $self->o('taxon_id'),                                                                     # Species level id, could be different to taxon_id if we have a subspecies, used to get species level RNA-seq CSV data
    'genus_taxon_id'               => '' || $self->o('taxon_id'),                                                                     # Genus level taxon id, used to get a genus level csv file in case there is not enough species level transcriptomic data
    'uniprot_set'                  => '',                                                                                             # e.g. mammals_basic, check UniProtCladeDownloadStatic.pm module in hive config dir for suitable set,
    'output_path'                  => '',                                                                                             # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses
    'assembly_name'                => '',                                                                                             # Name (as it appears in the assembly report file)
    'assembly_accession'           => '',                                                                                             # Versioned GCA assembly accession, e.g. GCA_001857705.1
    'assembly_refseq_accession'    => '',                                                                                             # Versioned GCF accession, e.g. GCF_001857705.1
    'stable_id_prefix'             => '',                                                                                             # e.g. ENSPTR. When running a new annotation look up prefix in the assembly registry db
    'use_genome_flatfile'          => '1',                                                                                            # This will read sequence where possible from a dumped flatfile instead of the core db
    'species_url'                  => '' || $self->o('production_name') . $self->o('production_name_modifier'),                       # sets species.url meta key
    'species_division'             => '',                                                                                             # sets species.division meta key
    'stable_id_start'              => '',                                                                                             # When mapping is not required this is usually set to 0
    'mapping_required'             => '0',                                                                                            # If set to 1 this will run stable_id mapping sometime in the future. At the moment it does nothing
    'mapping_db'                   => '',                                                                                             # Tied to mapping_required being set to 1, we should have a mapping db defined in this case, leave undef for now
    'uniprot_version'              => 'uniprot_2019_04',                                                                              # What UniProt data dir to use for various analyses
    'vertrna_version'              => '136',                                                                                          # The version of VertRNA to use, should correspond to a numbered dir in VertRNA dir
    'paired_end_only'              => '1',                                                                                            # Will only use paired-end rnaseq data if 1
    'ig_tr_fasta_file'             => 'human_ig_tr.fa',                                                                               # What IMGT fasta file to use. File should contain protein segments with appropriate headers
    'mt_accession'                 => undef,                                                                                          # This should be set to undef unless you know what you are doing. If you specify an accession, then you need to add the parameters to the load_mitochondrion analysis
    'production_name_modifier'     => '',                                                                                             # Do not set unless working with non-reference strains, breeds etc. Must include _ in modifier, e.g. _hni for medaka strain HNI
    'compara_registry_file'        => '',

    # Keys for custom loading, only set/modify if that's what you're doing
    'skip_genscan_blasts'       => '1',
    'load_toplevel_only'        => '1',                                                                                                                # This will not load the assembly info and will instead take any chromosomes, unplaced and unlocalised scaffolds directly in the DNA table
    'custom_toplevel_file_path' => '',                                                                                                                 # Only set this if you are loading a custom toplevel, requires load_toplevel_only to also be set to 2
    'repeatmodeler_library'     => '',                                                                                                                 # This should be the path to a custom repeat library, leave blank if none exists
    'use_repeatmodeler_to_mask' => '0',                                                                                                                # Setting this will include the repeatmodeler library in the masking process
    'protein_blast_db'          => '' || catfile( $self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), 'PE12_vertebrata' ),         # Blast database for comparing the final models to.
    'protein_blast_index'       => '' || catdir( $self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), 'PE12_vertebrata_index' ),    # Indicate Index for the blast database.
    'protein_entry_loc'         => catfile( $self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), 'entry_loc' ),                     # Used by genscan blasts and optimise daf/paf. Don't change unless you know what you're doing

    'softmask_logic_names' => [],

########################
# Pipe and ref db info
########################

    'projection_source_db_name'         => 'homo_sapiens_core_98_38',                                                                                  # This is generally a pre-existing db, like the current human/mouse core for example
    'projection_source_db_server'       => 'mysql-ens-mirror-1',
    'projection_source_db_port'         => '4240',
    'projection_source_production_name' => 'homo_sapiens',

    # The following might not be known in advance, since the come from other pipelines
    # These values can be replaced in the analysis_base table if they're not known yet
    # If they are not needed (i.e. no projection or rnaseq) then leave them as is

    'provider_name' => 'Ensembl',
    'provider_url'  => 'www.ensembl.org',

    'pipe_db_name' => $self->o('dbowner') . '_' . $self->o('production_name') . $self->o('production_name_modifier') . '_pipe_' . $self->o('release_number'),
    'dna_db_name'  => $self->o('dbowner') . '_' . $self->o('production_name') . $self->o('production_name_modifier') . '_core_' . $self->o('release_number'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'ensembl_release'      => $ENV{ENSEMBL_RELEASE},     # this is the current release version on staging to be able to get the correct database
    'production_db_server' => 'mysql-ens-meta-prod-1',
    'production_db_port'   => '4483',

    databases_to_delete => [ 'reference_db', 'cdna_db', 'genblast_db', 'genewise_db', 'projection_db', 'selected_projection_db', 'layering_db', 'utr_db', 'genebuilder_db', 'pseudogene_db', 'ncrna_db', 'final_geneset_db', 'refseq_db', 'cdna2genome_db', 'rnaseq_blast_db', 'rnaseq_refine_db', 'rnaseq_rough_db', 'lincrna_db', 'otherfeatures_db', 'rnaseq_db' ],    #, 'projection_realign_db'

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
    # This one is used by most analyses that run against a genome flatfile like exonerate, genblast etc. Has slice name style headers. Is softmasked
    softmasked_genome_file => catfile( $self->o('genome_dumps'), $self->o('species_name') . '_softmasked_toplevel.fa' ),
    # This one is used in replacement of the dna table in the core db, so where analyses override slice->seq. Has simple headers with just the seq_region name. Also used by bwa in the RNA-seq analyses. Not masked
    faidx_genome_file => catfile( $self->o('genome_dumps'), $self->o('species_name') . '_toplevel.fa' ),
    # This one is a cross between the two above, it has the seq_region name header but is softmasked. It is used by things that would both want to skip using the dna table and also want to avoid the repeat_feature table, e.g. bam2introns

    ensembl_analysis_script => catdir( $self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts' ),

    ensembl_misc_script      => catdir( $self->o('enscode_root_dir'), 'ensembl', 'misc-scripts' ),
    meta_coord_script        => catfile( $self->o('ensembl_misc_script'), 'meta_coord', 'update_meta_coord.pl' ),
    meta_levels_script       => catfile( $self->o('ensembl_misc_script'), 'meta_levels.pl' ),
    frameshift_attrib_script => catfile( $self->o('ensembl_misc_script'), 'frameshift_transcript_attribs.pl' ),
    select_canonical_script  => catfile( $self->o('ensembl_misc_script'), 'canonical_transcripts', 'select_canonical_transcripts.pl' ),
    stable_id_convert_script => catfile( $self->o('ensembl_analysis_script'), 'genebuild', 'braker', 'stable_id_convert.pl' ),

########################
# Extra db settings
########################

    'num_tokens'       => 10,
    mysql_dump_options => '--max_allowed_packet=1000MB',

########################
# Executable paths
########################
    'blast_type'  => 'ncbi',                                            # It can be 'ncbi', 'wu', or 'legacy_ncbi'
    samtools_path => catfile( $self->o('binary_base'), 'samtools' ),    #You may need to specify the full path to the samtools binary

    'uniprot_genblast_batch_size' => 15,
    'uniprot_table_name'          => 'uniprot_sequences',

    cdna_table_name => 'cdna_sequences',

# RNA-seq pipeline stuff
    # You have the choice between:
    #  * using a csv file you already created
    #  * using a study_accession like PRJEB19386
    #  * using the taxon_id of your species
    # 'rnaseq_summary_file' should always be set. If 'taxon_id' or 'study_accession' are not undef
    # they will be used to retrieve the information from ENA and to create the csv file. In this case,
    # 'file_columns' and 'summary_file_delimiter' should not be changed unless you know what you are doing
    'summary_csv_table' => 'csv_data',
    'read_length_table' => 'read_length',

    'rnaseq_dir'  => catdir( $self->o('output_path'), 'rnaseq' ),
    'input_dir'   => catdir( $self->o('rnaseq_dir'),  'input' ),
    'output_dir'  => catdir( $self->o('rnaseq_dir'),  'output' ),
    'merge_dir'   => catdir( $self->o('rnaseq_dir'),  'merge' ),
    'sam_dir'     => catdir( $self->o('rnaseq_dir'),  'sams' ),
    'header_file' => catfile( $self->o('output_dir'), '#' . $self->o('read_id_tag') . '#_header.h' ),

    'rnaseq_ftp_base' => 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/',

    'long_read_dir'       => catdir( $self->o('output_path'),   'long_read' ),
    'long_read_fastq_dir' => catdir( $self->o('long_read_dir'), 'input' ),
    'use_ucsc_naming'     => 0,

    # What Read group tag would you like to group your samples
    # by? Default = ID
    read_group_tag => 'SM',
    read_id_tag    => 'ID',

    use_threads          => 3,
    rnaseq_merge_threads => 12,
    rnaseq_merge_type    => 'samtools',
    read_min_paired      => 50,
    read_min_mapped      => 50,
    other_isoforms       => 'other',      # If you don't want isoforms, set this to undef
    maxintron            => 200000,

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

    'realign_table_name' => 'projection_source_sequences',

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
    'refseq_base_ftp'        => $self->o('ncbi_base_ftp') . '/#expr(substr(#assembly_refseq_accession#, 0, 3))expr#/#expr(substr(#assembly_refseq_accession#, 4, 3))expr#/#expr(substr(#assembly_refseq_accession#, 7, 3))expr#/#expr(substr(#assembly_refseq_accession#, 10, 3))expr#/#assembly_refseq_accession#_#assembly_name#',
    'refseq_import_ftp_path' => $self->o('refseq_base_ftp') . '/#assembly_refseq_accession#_#assembly_name#_genomic.gff.gz',
    'refseq_mrna_ftp_path'   => $self->o('refseq_base_ftp') . '/#assembly_refseq_accession#_#assembly_name#_rna.fna.gz',
    'refseq_report_ftp_path' => $self->o('refseq_base_ftp') . '/#assembly_refseq_accession#_#assembly_name#_assembly_report.txt',


########################
# db info
########################
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
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -dbname => $self->o('registry_db_name'),
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
    species_name         => $self->o('species_name'),
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
        'num_threads'      => $self->o('num_threads'),
        'dbowner'          => $self->o('dbowner'),
        'core_db'          => $self->o('core_db'),
        'ensembl_release'  => $self->o('ensembl_release'),
        'base_output_dir'  => $self->o('base_output_dir'),
        'registry_db'      => $self->o('registry_db'),
        'enscode_root_dir' => $self->o('enscode_root_dir'),
        'registry_file'    => $self->o('registry_file'),
      },
      -rc_name => 'default',

      -flow_into => {
        1 => ['download_rnaseq_csv'],
      },
      -analysis_capacity => 1,
      -input_ids         => [
        #{ 'assembly_accession' => 'GCA_' },
      ],

    },

    {
      -logic_name => 'download_rnaseq_csv',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadCsvENA',
      -rc_name    => '1GB',
      -parameters => {
        # This is specifically for gbiab, as taxon_id is populated in the input id with
        # the actual taxon, had to add some code override. Definitely better solutions available,
        # one might be to just branch this off and then only pass the genus taxon id
        override_taxon_id => 1,
        taxon_id          => '#genus_taxon_id#',
        inputfile         => '#rnaseq_summary_file#',
      },

      -flow_into => {
        '1->A' => { 'fan_short_read_download' => { 'inputfile' => '#rnaseq_summary_file#', 'input_dir' => '#short_read_dir#' } },
        'A->1' => ['download_long_read_csv'],
#       1 => ['fan_rnaseq_data_available'],
#        1 => { 'fan_rnaseq_data_available' => { 'long_read_dir' => '#long_read_dir#', 'short_read_dir' => '#short_read_dir#','output_path' => '#output_path#', 'reheadered_toplevel_genome_file' => '#reheadered_toplevel_genome_file#', 'species_name' => '#species_name#' } },
#               1 => {'run_braker_ab_initio' => {'output_path' => '#output_path#', 'reheadered_toplevel_genome_file' => '#reheadered_toplevel_genome_file#', 'species_name' => '#species_name#'}},
      },
    },

    {
      -logic_name => 'fan_short_read_download',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd                     => 'if [ -e "#inputfile#" ]; then exit 0; else exit 42;fi',
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
        2 => { 'download_short_read_fastqs' => { 'iid' => '#filename#', 'input_dir' => '#input_dir#' } },
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
    },

    {
      -logic_name => 'download_long_read_csv',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadCsvENA',
      -rc_name    => '1GB',
      -parameters => {
        # This is specifically for gbiab, as taxon_id is populated in the input id with
        # the actual taxon, had to add some code override. Definitely better solutions available,
        # one might be to just branch this off and then only pass the genus taxon id
        override_taxon_id => 1,
        taxon_id          => '#genus_taxon_id#',
        read_type         => 'isoseq',
        inputfile         => '#long_read_summary_file#',
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
        cmd                     => 'if [ -e "#inputfile#" ]; then exit 0; else exit 42;fi',
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
        column_names => $self->o('long_read_columns'),
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
        1 => ['process_assembly_info'],
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
      -max_retry_count => 0,
      -flow_into       => {
        1 => ['load_meta_info'],
      },
    },

    {
      # Load some meta info and seq_region_synonyms
      -logic_name => 'load_meta_info',
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
            '(1, "species.strain_group", "#species_strain_group#"),' .
            '(1, "species.production_name", "#production_name#"),' .
            '(1, "strain.type", "#strain_type#"),' .
            '(1, "repeat.analysis", "repeatdetector"),' .
            '(1, "repeat.analysis", "dust"),' .
            '(1, "repeat.analysis", "trf"),' .
            '(1, "genebuild.initial_release_date", NULL),' .
            '(1, "genebuild.projection_source_db", NULL),' .
            '(1, "genebuild.id", ' . $self->o('genebuilder_id') . '),' .
            '(1, "genebuild.method", "full_genebuild")'
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
        1 => ['run_gbiab_softmasking'],
      },
    },

    {
      -logic_name => 'run_gbiab_softmasking',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'python3.7 ' . catfile( $self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'genebuild', 'gbiab', 'gbiab.py' ) . ' #gbiab_red_commandline#;' .
          'cp #output_path#/red_output/mask_output/#species_name#_reheadered_toplevel.msk #output_path#/#species_name#_softmasked_toplevel.fa',
      },
      -rc_name         => 'gbiab',
      -max_retry_count => 0,
      -flow_into       => {
        1 => ['run_braker_ep_mode'],
      },
    },
######################BRAKER###########
    #{
    #  -logic_name => 'fan_rnaseq_data_available',
    #  -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    #  -parameters => {
    #    cmd                     => 'if [ "$(ls -A #short_read_dir#)" && "$(ls -A #long_read_dir#)" ]; then exit 0; else exit 42;fi',
    #    return_codes_2_branches => { '42' => 2 },
    #  },
    #  -rc_name   => 'default',
    #  -flow_into => {
	      ##        1 => ['run_gbiab'],
	      #    2 => ['run_braker_ep_mode'],
	      #  },
	      # },
    #    {
    #      -logic_name => 'run_gbiab',
    #      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    #      -parameters => {
    #        cmd => 'python3.7 ' . catfile( $self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'genebuild', 'gbiab', 'gbiab.py' ) . ' #gbiab_commandline#',
    #      },
    #      -rc_name         => 'gbiab',
    #      -max_retry_count => 0,
    #      -flow_into       => {
    #        1 => ['fan_braker_etp_setup'],
    #      },
    #    },
    #    {
    #      -logic_name => 'fan_braker_etp_setup',
    #      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    #      -parameters => {
    #        cmd                     => 'if [ -d "' . $self->o('augustus_species_path') . '#species_name#" ]; then exit 0; else exit 42;fi',
    #        return_codes_2_branches => { '42' => 2 },
    #      },
    #      -rc_name   => 'default',
    #      -flow_into => {
    #        1      => ['run_braker_etp_mode'],
    #        '2->A' => ['run_braker_ab_initio'],
    #        'A->2' => ['run_braker_etp_mode'],
    #      },
    #    },
    # {
    #  -logic_name => 'fan_braker_ep_setup',
    #  -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
    #  -parameters => {
    #    cmd                     => 'if [ -d "' . $self->o('augustus_species_path') . '#species_name#" ]; then exit 0; else exit 42;fi',
    #    return_codes_2_branches => { '42' => 2 },
    #  },
    #  -rc_name   => 'default',
    #  -flow_into => {
    #    1      => ['run_braker_ep_mode'],
    #    '2->A' => ['run_braker_ab_initio'],
    #    'A->2' => ['run_braker_ep_mode'],
    #  },
    #},
    #{
    #  -logic_name => 'run_braker_ab_initio',
    #  -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',

    #  -parameters => {
    #    cmd => 'cd #output_path#;' .
    #      'singularity exec --bind #output_path#/:/data:rw  ' . $self->o('braker_singularity_image') . ' braker.pl --genome=#species_name#_softmasked_toplevel.fa --softmasking --esmode --species=#species_name#  --AUGUSTUS_CONFIG_PATH=' . $self->o('augustus_config_path') . ' --cores ' . $self->o('cores') . ';' .
    #      'mv braker braker_ab_initio;',
    #  },
    #  -rc_name         => 'braker32',
    #  -max_retry_count => 0,
    #  -flow_into       => {
#          1 => ['load_gtf_file'],
    #  },
    #},
    #    {
    #      -logic_name => 'run_braker_etp_mode',
    #      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',

    #      -parameters => {
    #        cmd => 'cd #output_path#;' .

    #          'singularity exec --bind #output_path#/:/data:rw  ' . $self->o('braker_singularity_image') . ' braker.pl --genome=#species_name#_softmasked_toplevel.fa --prot_seq=/data/#protein_file# --bam=$bam_file --softmasking --etpmode --species=#species_name#  --useexisting --AUGUSTUS_CONFIG_PATH=' . $self->o('augustus_config_path') . ' --cores ' . $self->o('cores') . ';',
#          'mv braker braker_etp_mode;',
    #      },
    #      -rc_name         => 'braker32',
    #      -max_retry_count => 0,
    #      -flow_into       => {
    #        1 => ['clean_gtf_file'],
    #      },
    #    },
    {
      -logic_name => 'run_braker_ep_mode',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',

      -parameters => {
        cmd => 'mkdir #output_path#/prothint;' .
          'cd #output_path#/prothint;' .
          'singularity exec --bind #output_path#/prothint/:/data:rw  ' . $self->o('braker_singularity_image') . ' prothint.py #output_path#/#species_name#_softmasked_toplevel.fa #protein_file# ;' .
          'cd #output_path#/;' .
          'singularity exec --bind #output_path#/:/data:rw  ' . $self->o('braker_singularity_image') . ' braker.pl --genome=#species_name#_softmasked_toplevel.fa --softmasking  --hints=/data/prothint/prothint_augustus.gff --prothints=/data/prothint/prothint.gff --evidence=/data/prothint/evidence.gff --epmode --species=#assembly_accession#_#species_name# --AUGUSTUS_CONFIG_PATH=' . $self->o('augustus_config_path') . ' --cores ' . $self->o('cores') . ';',
      },
      -rc_name         => 'braker32',
      -max_retry_count => 0,
      -flow_into       => {
        1 => ['load_gtf_file'],
      },
    },
    # {
    #  -logic_name => 'clean_gtf_file',
    #  -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',

    #  -parameters => {
    # #        cmd => 'python ' . catfile( $self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'genebuild', 'braker', 'process_braker_gtf.py' ) . ' ./braker/braker.gtf braker_gtf.gtf;' .
        #          'sed "s/ \+ /\t/g" braker_gtf.gtf > braker_gtf_new.gtf',
	#    cmd => 'singularity exec --bind #output_path#/braker/:/data:rw  ' . $self->o('python_singularity_image') . 'python ' . catfile( $self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'genebuild', 'braker', 'process_braker_gtf.py' ) . ' #output_path#/braker/braker.gtf #output_path#/braker/braker_gtf.gtf;',
	# },
	#  -rc_name         => 'braker32',
	#  -max_retry_count => 0,
	#  -flow_into       => {
	#   1 => ['load_gtf_file'],
	# },
	# },
    {
      -logic_name => 'load_gtf_file',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',

      -parameters => {
        cmd => 'perl ' . catfile( $self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'genebuild', 'braker', 'parse_gtf.pl' ) .
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
          ' -file #output_path#/braker/braker.gtf',
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
        1 => ['run_stable_ids'],
      },
    },

    {
      -logic_name => 'run_stable_ids',
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

    # {
    #-logic_name => 'update_stable_ids',
    #-module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::SystemCmd',
    #-parameters => {
    #  cmd => 'perl ' . $self->o('stable_id_convert_script') .
    #    ' -user ' . $self->o('user') .
    #    ' -pass ' . $self->o('password') .
    #    ' -host ' . $self->o('core_db', '-host') .
    #    ' -port ' . $self->o('core_db', '-port') .
    #    ' -dbname ' . '#core_dbname#'
    #},
    #-rc_name   => 'default',
    #-flow_into => {
    #  1 => ['final_meta_updates'],
    #},
    #},

    {
      -logic_name => 'final_meta_updates',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => '#core_db#',
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
        db_conn => '#core_db#',
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
      -rc_name   => 'default_registry',
      -flow_into => {
        1 => ['run_agat_protein_file'],
      },
    },
    {
      -logic_name => 'run_agat_protein_file',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',

      -parameters => {
        cmd => 'singularity exec --bind #output_path#/:/data:rw  ' . $self->o('agat_singularity_image') . ' agat_sp_extract_sequences.pl --gff /data/braker/braker.gtf -f  #output_path#/#species_name#_softmasked_toplevel.fa -p  -o  #output_path#/braker/braker_proteins.fa;',
      },
      -rc_name         => 'braker32',
      -max_retry_count => 0,
      -flow_into       => {
        1 => ['run_busco'],
      },
    },
    {
      -logic_name => 'run_busco',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',

      -parameters => {
        cmd => 'singularity exec ' . $self->o('busco_singularity_image') . ' busco -i #output_path#/braker/braker_proteins.fa -m prot -l #busco_group# -o output_busco_#assembly_accession#  ;',
      },
      -rc_name => 'braker32',
    },
  ];
}

sub resource_classes {
  my $self = shift;

  return {
    'default_registry' => { LSF => [ $self->lsf_resource_builder( 'production', 900, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ), '-reg_conf ' . $self->default_options->{'registry_file'} ] },
    'gbiab'            => { LSF => $self->lsf_resource_builder( 'production', 50000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], undef, $self->default_options->{'num_threads'} ) },
    'braker32'         => { LSF => $self->lsf_resource_builder( 'production', 32000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], undef, $self->default_options->{'cores'} ) },
    '1GB'              => { LSF => $self->lsf_resource_builder( 'production', 1000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    '2GB_lastz'        => { LSF => [ $self->lsf_resource_builder( 'production', 2000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ), '-reg_conf ' . $self->default_options->{compara_registry_file} ] },
    '2GB'              => { LSF => $self->lsf_resource_builder( 'production', 2000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    '3GB'              => { LSF => $self->lsf_resource_builder( 'production', 3000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    '8GB_lastz'        => { LSF => [ $self->lsf_resource_builder( 'production', 8000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ), '-reg_conf ' . $self->default_options->{compara_registry_file} ] },
    '8GB'              => { LSF => $self->lsf_resource_builder( 'production', 8000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    '10GB'             => { LSF => $self->lsf_resource_builder( 'production', 10000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    '15GB_lastz'       => { LSF => [ $self->lsf_resource_builder( 'production', 15000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ), '-reg_conf ' . $self->default_options->{compara_registry_file} ] },
    'default'          => { LSF => $self->lsf_resource_builder( 'production', 900, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    'blast'            => { LSF => $self->lsf_resource_builder( 'production', 2900, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], undef, 3 ) },
    'genblast'         => { LSF => $self->lsf_resource_builder( 'production', 3900, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'exonerate_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
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
