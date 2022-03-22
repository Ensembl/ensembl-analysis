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

package EnsemblAnnoBraker_conf;

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
    #BRAKER parameters
    'augustus_config_path'     => '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/augustus_config/config',
    'augustus_species_path'    => '/nfs/production/flicek/ensembl/genebuild/genebuild_virtual_user/augustus_config/config/species/',
    'braker_singularity_image' => '/hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/test-braker2_es_ep_etp.simg',
    'agat_singularity_image'   => '/hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/test-agat.simg',
    'busco_singularity_image' => '/hps/software/users/ensembl/genebuild/genebuild_virtual_user/singularity/busco-v5.1.2_cv1.simg',

    'run_braker'               => 1,

    'current_genebuild'            => 1,
    'cores'                        => 30,
    'num_threads'                  => 20,
    'dbowner'                      => '',                                                                                                                # || $ENV{EHIVE_USER} || $ENV{USER},
    'base_output_dir'              => '',
    'init_config'               => '', #path for configuration file (custom loading)
    'override_clade'               => '',
    'protein_file'                 => '',
    'busco_protein_file'           => '',
    'rfam_accessions_file'         => '',
    'use_existing_short_read_dir'  => '',                                                                                                                #'/nfs/production/flicek/ensembl/genebuild/ftricomi/bees_data/bombus_short_read_fastq/',
    'registry_file'                => catfile( $self->o('enscode_root_dir'), 'ensembl-analysis/scripts/genebuild/gbiab/support_files/Databases.pm' ),    #'/hps/nobackup/flicek/ensembl/genebuild/ftricomi/braker_subpipe/Databases.pm', # This should be the path to the pipeline's copy of the Databases.pm registry file, core adaptors will be written to it
    'generic_registry_file'        => '',                                                                                                                # Could use this to hold the path to ensembl-analysis/scripts/genebuild/gbiab/support_files/Databases.pm to copy as a generic registry
    'diamond_validation_db'        => '/hps/nobackup/flicek/ensembl/genebuild/blastdb/uniprot_euk_diamond/uniprot_euk.fa.dmnd',
    'validation_type'              => 'moderate',
    'release_number'               => '' || $self->o('ensembl_release'),
    'production_name'              => '' || $self->o('species_name'),
    'pipeline_name'                => '' || $self->o('production_name') . $self->o('production_name_modifier') . '_' . $self->o('ensembl_release'),
    'user_r'                       => '',                                                                                                                # read only db user
    'user'                         => '',                                                                                                                # write db user
    'password'                     => '',                                                                                                                # password for write db user
    'server_set'                   => '',                                                                                                                # What server set to user, e.g. set1
    'busco_input_file_stid'        => 'stable_id_to_dump.txt',
    'pipe_db_server'               => $ENV{GBS7},                                                                                                        # host for pipe db
    'databases_server'             => $ENV{GBS5},                                                                                                        # host for general output dbs
    'dna_db_server'                => $ENV{GBS6},                                                                                                        # host for dna db
    'pipe_db_port'                 => $ENV{GBP7},                                                                                                        # port for pipeline host
    'databases_port'               => $ENV{GBP5},                                                                                                        # port for general output db host
    'dna_db_port'                  => $ENV{GBP6},                                                                                                        # port for dna db host
    'registry_db_server'           => $ENV{GBS1},                                                                                                        # host for registry db
    'registry_db_port'             => $ENV{GBP1},                                                                                                        # port for registry db
    'registry_db_name'             => 'gb_assembly_registry',                                                                                            # name for registry db
    'repbase_logic_name'           => '',                                                                                                                # repbase logic name i.e. repeatmask_repbase_XXXX, ONLY FILL THE XXXX BIT HERE!!! e.g primates
    'repbase_library'              => '',                                                                                                                # repbase library name, this is the actual repeat repbase library to use, e.g. "Mus musculus"
    'rnaseq_summary_file'          => '' || catfile( $self->o('rnaseq_dir'),    $self->o('species_name') . '.csv' ),                                     # Set this if you have a pre-existing cvs file with the expected columns
    'rnaseq_summary_file_genus'    => '' || catfile( $self->o('rnaseq_dir'),    $self->o('species_name') . '_gen.csv' ),                                 # Set this if you have a pre-existing genus level cvs file with the expected columns
    'long_read_summary_file'       => '' || catfile( $self->o('long_read_dir'), $self->o('species_name') . '_long_read.csv' ),                           # csv file for minimap2, should have 2 columns tab separated cols: sample_name\tfile_name
    'long_read_summary_file_genus' => '' || catfile( $self->o('long_read_dir'), $self->o('species_name') . '_long_read_gen.csv' ),                       # csv file for minimap2, should have 2 columns tab separated cols: sample_name\tfile_name
    'long_read_fastq_dir'          => '' || catdir( $self->o('long_read_dir'), 'input' ),
    'species_name'                 => '',                                                                                                                # e.g. mus_musculus
    'taxon_id'                     => '',                                                                                                                # should be in the assembly report file
    'species_taxon_id'             => '' || $self->o('taxon_id'),                                                                                        # Species level id, could be different to taxon_id if we have a subspecies, used to get species level RNA-seq CSV data
    'genus_taxon_id'               => '' || $self->o('taxon_id'),                                                                                        # Genus level taxon id, used to get a genus level csv file in case there is not enough species level transcriptomic data
    'uniprot_set'                  => '',                                                                                                                # e.g. mammals_basic, check UniProtCladeDownloadStatic.pm module in hive config dir for suitable set,
    'output_path'                  => '',                                                                                                                # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses
    'assembly_name'                => '',                                                                                                                # Name (as it appears in the assembly report file)
    'assembly_accession'           => '',                                                                                                                # Versioned GCA assembly accession, e.g. GCA_001857705.1
    'assembly_refseq_accession'    => '',                                                                                                                # Versioned GCF accession, e.g. GCF_001857705.1
    'stable_id_prefix'             => '',                                                                                                                # e.g. ENSPTR. When running a new annotation look up prefix in the assembly registry db
    'use_genome_flatfile'          => '1',                                                                                                               # This will read sequence where possible from a dumped flatfile instead of the core db
    'species_url'                  => '' || $self->o('production_name') . $self->o('production_name_modifier'),                                          # sets species.url meta key
    'species_division'             => '',                                                                                                                # sets species.division meta key
    'stable_id_start'              => '',                                                                                                                # When mapping is not required this is usually set to 0
    'mapping_required'             => '0',                                                                                                               # If set to 1 this will run stable_id mapping sometime in the future. At the moment it does nothing
    'mapping_db'                   => '',                                                                                                                # Tied to mapping_required being set to 1, we should have a mapping db defined in this case, leave undef for now
    'uniprot_version'              => 'uniprot_2019_04',                                                                                                 # What UniProt data dir to use for various analyses
    'vertrna_version'              => '136',                                                                                                             # The version of VertRNA to use, should correspond to a numbered dir in VertRNA dir
    'paired_end_only'              => '1',                                                                                                               # Will only use paired-end rnaseq data if 1
    'ig_tr_fasta_file'             => 'human_ig_tr.fa',                                                                                                  # What IMGT fasta file to use. File should contain protein segments with appropriate headers
    'mt_accession'                 => undef,                                                                                                             # This should be set to undef unless you know what you are doing. If you specify an accession, then you need to add the parameters to the load_mitochondrion analysis
    'production_name_modifier'     => '',                                                                                                                # Do not set unless working with non-reference strains, breeds etc. Must include _ in modifier, e.g. _hni for medaka strain HNI
    'compara_registry_file'        => '',

    # Keys for custom loading, only set/modify if that's what you're doing
    'skip_genscan_blasts'       => '1',
    'load_toplevel_only'        => '1',                                                                                                                  # This will not load the assembly info and will instead take any chromosomes, unplaced and unlocalised scaffolds directly in the DNA table
    'custom_toplevel_file_path' => '',                                                                                                                   # Only set this if you are loading a custom toplevel, requires load_toplevel_only to also be set to 2
    'repeatmodeler_library'     => '',                                                                                                                   # This should be the path to a custom repeat library, leave blank if none exists
    'use_repeatmodeler_to_mask' => '0',                                                                                                                  # Setting this will include the repeatmodeler library in the masking process
    'protein_blast_db'          => '' || catfile( $self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), 'PE12_vertebrata' ),           # Blast database for comparing the final models to.
    'protein_blast_index'       => '' || catdir( $self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), 'PE12_vertebrata_index' ),      # Indicate Index for the blast database.
    'protein_entry_loc'         => catfile( $self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_version'), 'entry_loc' ),                       # Used by genscan blasts and optimise daf/paf. Don't change unless you know what you're doing

    'softmask_logic_names' => [],
    'desired_slice_length' => 10000000,
    'store_rejected'       => 0,


########################
## Small ncRNAs params
#########################
    'mirBase_fasta' => 'all_mirnas.fa',                                                                                                                  # What mirBase file to use. It is currently best to use on with the most appropriate set for your species
    'rfc_scaler'    => 'filter_dafs_rfc_scaler_human.pkl',
    'rfc_model'     => 'filter_dafs_rfc_model_human.pkl',

    # Clade-based filtering on rfam accessions
    # Rfam db details should stay constant but check periodically
    'rfam_user'   => 'rfamro',
    'rfam_dbname' => 'Rfam',
    'rfam_host'   => 'mysql-rfam-public.ebi.ac.uk',
    'rfam_port'   => 4497,

    'rfam_path'        => catfile( $self->o('base_blast_db_path'), 'ncrna', 'Rfam_14.0' ),
    'rfam_seeds'       => $self->o('rfam_path') . "/Rfam.seed",
    'rfam_cm'          => $self->o('rfam_path') . "/Rfam.cm",
    'filtered_rfam_cm' => $self->o('output_path') . '/Rfam.cm',
    'clade'            => $self->o('repbase_logic_name'),


########################
# Pipe and ref db info
########################

    'projection_source_db_name'         => 'homo_sapiens_core_98_38',    # This is generally a pre-existing db, like the current human/mouse core for example
    'projection_source_db_server'       => 'mysql-ens-mirror-1',
    'projection_source_db_port'         => '4240',
    'projection_source_production_name' => 'homo_sapiens',

    'compara_db_name'   => 'leanne_ensembl_compara_95',
    'compara_db_server' => 'mysql-ens-genebuild-prod-5',
    'compara_db_port'   => 4531,

    # The following might not be known in advance, since the come from other pipelines
    # These values can be replaced in the analysis_base table if they're not known yet
    # If they are not needed (i.e. no projection or rnaseq) then leave them as is
    'projection_lastz_db_name'   => $self->o('pipe_db_name'),
    'projection_lastz_db_server' => $self->o('pipe_db_server'),
    'projection_lastz_db_port'   => $self->o('pipe_db_port'),

    'provider_name' => 'Ensembl',
    'provider_url'  => 'www.ensembl.org',

    'pipe_db_name' => $self->o('dbowner') . '_' . $self->o('production_name') . $self->o('production_name_modifier') . '_pipe_' . $self->o('release_number'),
    #'pipe_db_name'                  => $self->o('pipeline_name') . 'f_pipe_' . '105',
    'dna_db_name' => $self->o('dbowner') . '_' . $self->o('production_name') . $self->o('production_name_modifier') . '_core_' . $self->o('release_number'),

    'core_db_name'   => $self->o('dna_db_name'),
    'core_db_server' => $self->o('dna_db_server'),
    'core_db_port'   => $self->o('dna_db_port'),



    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'ensembl_release'      => $ENV{ENSEMBL_RELEASE},     # this is the current release version on staging to be able to get the correct database
    'production_db_server' => 'mysql-ens-meta-prod-1',
    'production_db_port'   => '4483',

########################
# BLAST db paths
########################
    'base_blast_db_path'    => $ENV{BLASTDB_DIR},
    'vertrna_blast_db_path' => catfile( $self->o('base_blast_db_path'), 'vertrna', $self->o('vertrna_version'), 'embl_vertrna-1' ),
    'unigene_blast_db_path' => catfile( $self->o('base_blast_db_path'), 'unigene', 'unigene' ),
    'ncrna_blast_path'      => catfile( $self->o('base_blast_db_path'), 'ncrna',   'ncrna_2016_05' ),
    'mirna_blast_path'      => catfile( $self->o('base_blast_db_path'), 'ncrna',   'mirbase_22' ),
    'ig_tr_blast_path'      => catfile( $self->o('base_blast_db_path'), 'ig_tr_genes' ),

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
    faidx_softmasked_genome_file             => catfile( $self->o('genome_dumps'), $self->o('species_name') . '_softmasked_toplevel.fa.reheader' ),
    'primary_assembly_dir_name'              => 'Primary_Assembly',
    'refseq_cdna_calculate_coverage_and_pid' => '0',
    'contigs_source'                         => 'ena',

    full_repbase_logic_name => "repeatmask_repbase_" . $self->o('repbase_logic_name'),

    'utr_biotype_priorities' => {
      'rnaseq' => 2,
      'cdna'   => 1,
    },

    'cleaning_blessed_biotypes' => {
      'pseudogene'           => 1,
      'processed_pseudogene' => 1,
      'IG_C_gene'            => 1,
      'IG_V_gene'            => 1,
      'TR_C_gene'            => 1,
      'TR_D_gene'            => 1,
      'TR_V_gene'            => 1,
      'lncRNA'               => 1,
    },

    'min_toplevel_slice_length' => 250,

    'repeatmodeler_logic_name' => 'repeatmask_repeatmodeler',
    'homology_models_path'     => catdir( $self->o('output_path'), 'homology_models' ),

    ncrna_dir       => catdir( $self->o('output_path'), 'ncrna' ),
    targetted_path  => catdir( $self->o('output_path'), 'targetted' ),
    cdna_file       => catfile( $self->o('targetted_path'), 'cdnas' ),
    annotation_file => $self->o('cdna_file') . '.annotation',

    ensembl_analysis_script           => catdir( $self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts' ),
    remove_duplicates_script_path     => catfile( $self->o('ensembl_analysis_script'), 'find_and_remove_duplicates.pl' ),
    flag_potential_pseudogenes_script => catfile( $self->o('ensembl_analysis_script'), 'genebuild', 'flag_potential_pseudogenes.pl' ),
    load_optimise_script              => catfile( $self->o('ensembl_analysis_script'), 'genebuild', 'load_external_db_ids_and_optimize_af.pl' ),
    prepare_cdnas_script              => catfile( $self->o('ensembl_analysis_script'), 'genebuild', 'prepare_cdnas.pl' ),
    load_fasta_script_path            => catfile( $self->o('ensembl_analysis_script'), 'genebuild', 'load_fasta_to_db_table.pl' ),
    loading_report_script             => catfile( $self->o('ensembl_analysis_script'), 'genebuild', 'report_genome_prep_stats.pl' ),
    refseq_synonyms_script_path       => catfile( $self->o('ensembl_analysis_script'), 'refseq',    'load_refseq_synonyms.pl' ),
    refseq_import_script_path         => catfile( $self->o('ensembl_analysis_script'), 'refseq',    'parse_ncbi_gff3.pl' ),
    sequence_dump_script              => catfile( $self->o('ensembl_analysis_script'), 'sequence_dump.pl' ),
    remove_duplicates_script          => catfile( $self->o('ensembl_analysis_script'), 'find_and_remove_duplicates.pl' ),
    sncrna_analysis_script            => catdir( $self->o('ensembl_analysis_script'), 'genebuild', 'sncrna' ),
    ensembl_misc_script               => catdir( $self->o('enscode_root_dir'),        'ensembl',   'misc-scripts' ),
    clean_lncrna_script               => catfile( $self->o('ensembl_analysis_script'), 'genebuild',  'gbiab', 'support_scripts_perl', 'clean_lncrna.pl' ),
    repeat_types_script               => catfile( $self->o('ensembl_misc_script'),     'repeats',    'repeat-types.pl' ),
    meta_coord_script                 => catfile( $self->o('ensembl_misc_script'),     'meta_coord', 'update_meta_coord.pl' ),
    meta_levels_script                => catfile( $self->o('ensembl_misc_script'),     'meta_levels.pl' ),
    frameshift_attrib_script          => catfile( $self->o('ensembl_misc_script'),     'frameshift_transcript_attribs.pl' ),
    select_canonical_script           => catfile( $self->o('ensembl_misc_script'),     'canonical_transcripts', 'select_canonical_transcripts.pl' ),
    assembly_name_script              => catfile( $self->o('ensembl_analysis_script'), 'update_assembly_name.pl' ),
    print_protein_script_path         => catfile( $self->o('ensembl_analysis_script'), 'genebuild', 'print_translations.pl' ),

    registry_status_update_script => catfile( $self->o('ensembl_analysis_script'), 'update_assembly_registry.pl' ),
    rnaseq_daf_introns_file       => catfile( $self->o('output_dir'),              'rnaseq_daf_introns.dat' ),

    # Genes biotypes to ignore from the final db when copying to core
    copy_biotypes_to_ignore => {
      'low_coverage' => 1,
      'CRISPR'       => 1,
    },
########################
# Extra db settings
########################

    'num_tokens'       => 10,
    mysql_dump_options => '--max_allowed_packet=1000MB',

########################
# Executable paths
########################
    'minimap2_genome_index' => $self->o('faidx_genome_file') . '.mmi',
    'minimap2_path'         => '/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/bin/minimap2',
    'paftools_path'         => '/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/bin/paftools.js',
    'minimap2_batch_size'   => '5000',

    'blast_type'             => 'ncbi',                                                                                                   # It can be 'ncbi', 'wu', or 'legacy_ncbi'
    'dust_path'              => catfile( $self->o('binary_base'),        'dustmasker' ),
    'trf_path'               => catfile( $self->o('binary_base'),        'trf' ),
    'eponine_java_path'      => catfile( $self->o('binary_base'),        'java' ),
    'eponine_jar_path'       => catfile( $self->o('software_base_path'), 'opt', 'eponine', 'libexec', 'eponine-scan.jar' ),
    'cpg_path'               => catfile( $self->o('binary_base'),        'cpg_lh' ),
    'trnascan_path'          => catfile( $self->o('binary_base'),        'tRNAscan-SE' ),
    'repeatmasker_path'      => catfile( $self->o('binary_base'),        'RepeatMasker' ),
    'genscan_path'           => catfile( $self->o('binary_base'),        'genscan' ),
    'genscan_matrix_path'    => catfile( $self->o('software_base_path'), 'share', 'HumanIso.smat' ),
    'uniprot_blast_exe_path' => catfile( $self->o('binary_base'),        'blastp' ),
    'blastn_exe_path'        => catfile( $self->o('binary_base'),        'blastn' ),
    'vertrna_blast_exe_path' => catfile( $self->o('binary_base'),        'tblastn' ),
    'unigene_blast_exe_path' => catfile( $self->o('binary_base'),        'tblastn' ),
    genewise_path            => catfile( $self->o('binary_base'),        'genewise' ),
    'exonerate_path'         => catfile( $self->o('software_base_path'), 'opt', 'exonerate09', 'bin', 'exonerate' ),
    'cmsearch_exe_path'      => catfile( $self->o('software_base_path'), 'bin', 'cmsearch' ),                                             # #'opt', 'infernal10', 'bin', 'cmsearch'),
    indicate_path            => catfile( $self->o('binary_base'),        'indicate' ),
    pmatch_path              => catfile( $self->o('binary_base'),        'pmatch' ),
    exonerate_annotation     => catfile( $self->o('binary_base'),        'exonerate' ),
    samtools_path            => catfile( $self->o('binary_base'),        'samtools' ),                                                    #You may need to specify the full path to the samtools binary
    picard_lib_jar           => catfile( $self->o('software_base_path'), 'Cellar', 'picard-tools', '2.6.0', 'libexec', 'picard.jar' ),    #You need to specify the full path to the picard library
    bwa_path                 => catfile( $self->o('software_base_path'), 'opt',    'bwa-051mt',    'bin',   'bwa' ),                      #You may need to specify the full path to the bwa binary
    refine_ccode_exe         => catfile( $self->o('binary_base'),        'RefineSolexaGenes' ),                                           #You may need to specify the full path to the RefineSolexaGenes binary
    interproscan_exe         => catfile( $self->o('binary_base'),        'interproscan.sh' ),
    bedtools                 => catfile( $self->o('binary_base'),        'bedtools' ),
    bedGraphToBigWig         => catfile( $self->o('binary_base'),        'bedGraphToBigWig' ),
    'cesar_path'             => catdir( $self->o('software_base_path'), 'opt', 'cesar', 'bin' ),

    'uniprot_genblast_batch_size' => 15,
    'uniprot_table_name'          => 'uniprot_sequences',

    'genblast_path'               => catfile( $self->o('binary_base'), 'genblast' ),
    'genblast_eval'               => $self->o('blast_type') eq 'wu' ? '1e-20' : '1e-1',
    'genblast_cov'                => '0.5',
    'genblast_pid'                => '30',
    'genblast_max_rank'           => '5',
    'genblast_flag_small_introns' => 1,
    'genblast_flag_subpar_models' => 0,

    'ig_tr_table_name'        => 'ig_tr_sequences',
    'ig_tr_genblast_cov'      => '0.8',
    'ig_tr_genblast_pid'      => '70',
    'ig_tr_genblast_eval'     => '1',
    'ig_tr_genblast_max_rank' => '5',
    'ig_tr_batch_size'        => 10,

    'exonerate_cdna_pid' => '95',    # Cut-off for percent id
    'exonerate_cdna_cov' => '50',    # Cut-off for coverage

    'cdna_selection_pid' => '97',    # Cut-off for percent id for selecting the cDNAs
    'cdna_selection_cov' => '90',    # Cut-off for coverage for selecting the cDNAs

# Best targetted stuff
    exonerate_logic_name => 'exonerate',
    ncbi_query           => '((txid' . $self->o('taxon_id') . '[Organism:noexp]+AND+biomol_mrna[PROP]))  NOT "tsa"[Properties]',

    cdna_table_name                             => 'cdna_sequences',
    target_exonerate_calculate_coverage_and_pid => 0,
    exonerate_protein_pid                       => 95,
    exonerate_protein_cov                       => 50,
    cdna2genome_region_padding                  => 2000,
    exonerate_max_intron                        => 200000,

    best_targetted_min_coverage => 50,    # This is to avoid having models based on fragment alignment and low identity
    best_targetted_min_identity => 50,    # This is to avoid having models based on fragment alignment and low identity


# RNA-seq pipeline stuff
    # You have the choice between:
    #  * using a csv file you already created
    #  * using a study_accession like PRJEB19386
    #  * using the taxon_id of your species
    # 'rnaseq_summary_file' should always be set. If 'taxon_id' or 'study_accession' are not undef
    # they will be used to retrieve the information from ENA and to create the csv file. In this case,
    # 'file_columns' and 'summary_file_delimiter' should not be changed unless you know what you are doing
    'study_accession'     => '',
    'max_reads_per_split' => 2500000,      # This sets the number of reads to split the fastq files on
    'max_total_reads'     => 200000000,    # This is the total number of reads to allow from a single, unsplit file

    'summary_file_delimiter' => '\t',            # Use this option to change the delimiter for your summary data file
    'summary_csv_table'      => 'csv_data',
    'read_length_table'      => 'read_length',
    'rnaseq_data_provider'   => 'ENA',           #It will be set during the pipeline or it will use this value

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
    split_paired_regex => '(\S+)(\_\d\.\S+)',
    split_single_regex => '([^.]+)(\.\S+)',

    # Do you want to make models for the each individual sample as well
    # as for the pooled samples (1/0)?
    single_tissue => 1,

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

# lincRNA pipeline stuff
    'lncrna_dir'               => catdir( $self->o('output_path'), 'lincrna' ),
    'file_translations'        => catfile( $self->o('lncrna_dir'), 'hive_dump_translations.fasta' ),
    'file_for_length'          => catfile( $self->o('lncrna_dir'), 'check_lincRNA_length.out' ),                              # list of genes that are smaller than 200bp, if any
    'file_for_biotypes'        => catfile( $self->o('lncrna_dir'), 'check_lincRNA_need_to_update_biotype_antisense.out' ),    # mysql queries that will apply or not in your dataset (check update_database) and will update biotypes
    'file_for_introns_support' => catfile( $self->o('lncrna_dir'), 'check_lincRNA_Introns_supporting_evidence.out' ),         # for debug
    biotype_output             => 'rnaseq',
    lincrna_protein_coding_set => [
      'rnaseq_merged_1',
      'rnaseq_merged_2',
      'rnaseq_merged_3',
      'rnaseq_merged_4',
      'rnaseq_merged_5',
      'rnaseq_tissue_1',
      'rnaseq_tissue_2',
      'rnaseq_tissue_3',
      'rnaseq_tissue_4',
      'rnaseq_tissue_5',
    ],

########################
# SPLIT PROTEOME File
########################
    'max_seqs_per_file'       => 20,
    'max_seq_length_per_file' => 20000,                                 # Maximum sequence length in a file
    'max_files_per_directory' => 1000,                                  # Maximum number of files in a directory
    'max_dirs_per_directory'  => $self->o('max_files_per_directory'),

########################
# FINAL Checks parameters - Update biotypes to lincRNA, antisense, sense, problem ...
########################

    update_database => 'yes',                                           # Do you want to apply the suggested biotypes? yes or no.

########################
# Interproscan
########################
    required_externalDb              => '',
    interproscan_lookup_applications => [
      'PfamA',
    ],
    required_externalDb => [],
    pathway_sources     => [],
    required_analysis   => [
      {
        'logic_name'    => 'pfam',
        'db'            => 'Pfam',
        'db_version'    => '31.0',
        'ipscan_name'   => 'Pfam',
        'ipscan_xml'    => 'PFAM',
        'ipscan_lookup' => 1,
      },
    ],




# Max internal stops for projected transcripts
    'projection_pid'                        => '50',
    'projection_cov'                        => '50',
    'projection_max_internal_stops'         => '1',
    'projection_calculate_coverage_and_pid' => '1',

    'projection_lincrna_percent_id'    => 90,
    'projection_lincrna_coverage'      => 90,
    'projection_pseudogene_percent_id' => 60,
    'projection_pseudogene_coverage'   => 75,
    'projection_ig_tr_percent_id'      => 70,
    'projection_ig_tr_coverage'        => 90,
    'projection_exonerate_padding'     => 5000,

    'realign_table_name'               => 'projection_source_sequences',
    'max_projection_structural_issues' => 1,

## Add in genewise path and put in matching code
    'genewise_pid'                        => '50',
    'genewise_cov'                        => '50',
    'genewise_region_padding'             => '50000',
    'genewise_calculate_coverage_and_pid' => '1',

########################
# Misc setup info
########################
    'repeatmasker_engine' => 'crossmatch',
    'masking_timer_long'  => '5h',
    'masking_timer_short' => '2h',

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

##################################
# Memory settings for the analyses
##################################
    'default_mem'          => '900',
    'genblast_mem'         => '1900',
    'genblast_retry_mem'   => '4900',
    'genewise_mem'         => '3900',
    'genewise_retry_mem'   => '5900',
    'refseq_mem'           => '9900',
    'projection_mem'       => '1900',
    'layer_annotation_mem' => '3900',
    'genebuilder_mem'      => '1900',

########################
# LastZ
########################

    'compara_master'                    => 'compara_master',
    'compara_conf_file'                 => '',
    'compara_innodb_schema'             => 1,
    'compara_genome_db_update_path'     => catfile( $self->o('enscode_root_dir'), '/ensembl-compara/scripts/pipeline/update_genome.pl' ),
    'compara_mlss_script_path'          => catfile( $self->o('enscode_root_dir'), '/ensembl-compara/scripts/pipeline/create_mlss.pl' ),
    'compara_mlss_reg_conf_path'        => catfile( $self->o('enscode_root_dir'), '/ensembl-compara/scripts/pipeline/production_reg_ensembl_conf.pl' ),
    'compara_populate_new_database_exe' => catfile( $self->o('enscode_root_dir'), 'ensembl-compara/scripts/pipeline/populate_new_database.pl' ),
    'compara_only_cellular_component'   => undef,
    'compara_dump_dir'                  => catdir( $self->o('output_path'), 'lastz' ),

    'mlss_id_list'       => undef,
    'compara_collection' => '',

    'compara_ref_species'     => $self->o('projection_source_production_name'),
    'compara_non_ref_species' => $self->o('production_name'),
    'only_cellular_component' => undef,                                           # Do we load *all* the dnafrags or only the ones from a specific cellular-component ?
    'mix_cellular_components' => 0,                                               # Do we try to allow the nuclear genome vs MT, etc ?
    'dump_min_nib_size'       => 11500000,
    'dump_min_chunk_size'     => 1000000,
    'dump_min_chunkset_size'  => 1000000,
    'quick'                   => 1,
    'default_chunks'          => {
      'reference' => {
        'homo_sapiens' => {
          'chunk_size'            => 30000000,
          'overlap'               => 0,
          'include_non_reference' => -1,                              #1  => include non_reference regions (eg human assembly patches)
                                                                      #0  => do not include non_reference regions
                                                                      #-1 => auto-detect (only include non_reference regions if the non-reference species is high-coverage
                                                                      #ie has chromosomes since these analyses are the only ones we keep up-to-date with the patches-pipeline)
          'masking_options'       => '{default_soft_masking => 1}',
          # if you have a specific selection of repeat elements for the masking
          #'masking_options_file' => $self->check_file_in_ensembl('ensembl-compara/scripts/pipeline/human36.spec'),
        },
        #non human example
        'default' => {
          'chunk_size'      => 10000000,
          'overlap'         => 0,
          'masking_options' => '{default_soft_masking => 1}'
        },
      },
      'non_reference' => {
        'chunk_size'      => 10100000,
        'group_set_size'  => 10100000,
        'overlap'         => 100000,
        'masking_options' => '{default_soft_masking => 1}'
      },
    },

    'compara_window_size'             => 10000,
    'filter_duplicates_rc_name'       => '2GB_lastz',
    'filter_duplicates_himem_rc_name' => '8GB_lastz',

    #
    #Default pair_aligner
    #
    'pair_aligner_method_link' => [ 1001, 'LASTZ_RAW' ],
    'pair_aligner_logic_name'  => 'LastZ',
    'pair_aligner_module'      => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::LastZ',

    'pair_aligner_options' => {
      default => 'T=1 L=3000 H=2200 O=400 E=30 --ambiguous=iupac',                                                                                                          # ensembl genomes settings
      7742    => 'T=1 K=3000 L=3000 H=2200 O=400 E=30 --ambiguous=iupac',                                                                                                   # vertebrates - i.e. ensembl-specific
      9526    => 'T=1 K=5000 L=5000 H=3000 M=10 O=400 E=30 Q=' . $self->check_file_in_ensembl('ensembl-compara/scripts/pipeline/primate.matrix') . ' --ambiguous=iupac',    # primates
      33554   => 'T=1 K=5000 L=5000 H=3000 M=10 O=400 E=30 --ambiguous=iupac',                                                                                              # carnivora
      3913    => 'T=1 L=3000 H=2200 O=400 E=30 --ambiguous=iupac --matchcount=1000',
      4070    => 'T=1 L=3000 H=2200 O=400 E=30 --ambiguous=iupac --matchcount=1000',
    },

    #
    #Default chain
    #
    'chain_input_method_link'  => [ 1001, 'LASTZ_RAW' ],
    'chain_output_method_link' => [ 1002, 'LASTZ_CHAIN' ],

    #linear_gap=>medium for more closely related species, 'loose' for more distant
    'linear_gap' => 'medium',

    'chain_parameters' => { 'max_gap' => '50', 'linear_gap' => $self->o('linear_gap'), 'faToNib' => $self->o('faToNib_exe'), 'lavToAxt' => $self->o('lavToAxt_exe'), 'axtChain' => $self->o('axtChain_exe'), 'max_blocks_for_chaining' => 100000 },

    #
    #Default patch_alignments
    #
    'patch_alignments' => 0,    #set to 1 to align the patches of a species to many other species

    #
    #Default net
    #
    'net_input_method_link'  => [ 1002, 'LASTZ_CHAIN' ],
    'net_output_method_link' => [ 16,   'LASTZ_NET' ],
    'net_ref_species'        => $self->o('compara_ref_species'),                                 #default to ref_species
    'net_parameters'         => { 'max_gap' => '50', 'chainNet' => $self->o('chainNet_exe') },
    'bidirectional'          => 0,

    #
    #Default healthcheck
    #
    'previous_db'               => 'compara_prev',
    'prev_release'              => 0,                # 0 is the default and it means "take current release number and subtract 1"
    'max_percent_diff'          => 20,
    'max_percent_diff_patches'  => 99.99,
    'do_pairwise_gabs'          => 1,
    'do_compare_to_previous_db' => 0,

    'compara_bed_dir'     => $self->o('compara_dump_dir') . '/bed_dir',
    'compara_feature_dir' => $self->o('compara_dump_dir') . '/feature_dumps',

    #
    #Default pairaligner config
    #
    'skip_pairaligner_stats' => 1,    #skip this module if set to 1

    'pair_aligner_method_link' => [ 1001, 'LASTZ_RAW' ],
    'pair_aligner_logic_name'  => 'LastZ',
    'pair_aligner_module'      => 'Bio::EnsEMBL::Compara::RunnableDB::PairAligner::LastZ',
    'chain_input_method_link'  => [ 1001, 'LASTZ_RAW' ],
    'chain_output_method_link' => [ 1002, 'LASTZ_CHAIN' ],
    'linear_gap'               => 'medium',
    'net_input_method_link'    => [ 1002, 'LASTZ_CHAIN' ],
    'net_output_method_link'   => [ 16,   'LASTZ_NET' ],

    # Capacities
    'pair_aligner_analysis_capacity'  => 700,
    'pair_aligner_batch_size'         => 40,
    'chain_hive_capacity'             => 200,
    'chain_batch_size'                => 10,
    'net_hive_capacity'               => 300,
    'net_batch_size'                  => 10,
    'filter_duplicates_hive_capacity' => 200,
    'filter_duplicates_batch_size'    => 10,

    # LastZ is used to align the genomes
    'pair_aligner_exe'             => $self->o('lastz_exe'),
    'cellar_dir'                   => '/nfs/software/ensembl/RHEL7-JUL2017-core2/linuxbrew/Cellar/',
    'lastz_exe'                    => catfile( $self->o('cellar_dir'),       'lastz/1.04.00/bin/lastz' ),
    'axtChain_exe'                 => catfile( $self->o('cellar_dir'),       'kent/v335_1/bin/axtChain' ),
    'chainNet_exe'                 => catfile( $self->o('cellar_dir'),       'kent/v335_1/bin/chainNet' ),
    'faToNib_exe'                  => catfile( $self->o('cellar_dir'),       'kent/v335_1/bin/faToNib' ),
    'lavToAxt_exe'                 => catfile( $self->o('cellar_dir'),       'kent/v335_1/bin/lavToAxt' ),
    'compare_beds_exe'             => catfile( $self->o('enscode_root_dir'), 'ensembl-compara/scripts/pipeline/compare_beds.pl' ),
    'create_pair_aligner_page_exe' => catfile( $self->o('enscode_root_dir'), 'ensembl-compara/scripts/report/create_pair_aligner_page.pl' ),
    'dump_features_exe'            => catfile( $self->o('enscode_root_dir'), 'ensembl-compara/scripts/dumps/DumpMultiAlign.pl' ),

    #gene_db_host => $self->o('data_db_server'),
    #gene_db_port => $self->o('data_db_port'),
    #gene_db_name => $self->o('dbowner').'_'.$self->o('production_name').'_gene_'.$self->o('release_number'),

    clean_utr_db_host => $self->o('dna_db_server'),
    clean_utr_db_port => $self->o('dna_db_port'),
    clean_utr_db_name => $self->o('dbowner') . '_' . $self->o('production_name') . '_clean_utr_' . $self->o('release_number'),

    otherfeatures_db_host => $self->o('dna_db_server'),
    otherfeatures_db_port => $self->o('dna_db_port'),
    otherfeatures_db_name => $self->o('dbowner') . '_' . $self->o('production_name') . '_otherfeatures_' . $self->o('release_number'),

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
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -dbname => $self->o('registry_db_name'),
      -driver => $self->o('hive_driver'),
    },

    'clean_utr_db' => {
      -dbname => $self->o('clean_utr_db_name'),
      -host   => $self->o('clean_utr_db_host'),
      -port   => $self->o('clean_utr_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
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
    'otherfeatures_db' => {
      -dbname => $self->o('otherfeatures_db_name'),
      -host   => $self->o('otherfeatures_db_host'),
      -port   => $self->o('otherfeatures_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },
    databases_to_delete => [ 'reference_db', 'cdna_db', 'genblast_db', 'genewise_db', 'projection_db', 'selected_projection_db', 'layering_db', 'utr_db', 'genebuilder_db', 'pseudogene_db', 'ncrna_db', 'final_geneset_db', 'refseq_db', 'cdna2genome_db', 'rnaseq_blast_db', 'rnaseq_refine_db', 'rnaseq_rough_db', 'lincrna_db', 'rnaseq_db', 'clean_utr_db' ],

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
        'clean_utr_db'                => $self->o('clean_utr_db'),
        'otherfeatures_db'            => $self->o('otherfeatures_db'),
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
        1 => ['download_rnaseq_csv'],
      },
      -analysis_capacity => 1,
      -input_ids         => [
        #{'assembly_accession' => 'GCA_910591885.1'},
        #	{'assembly_accession' => 'GCA_905333015.1'},
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
        input_dir         => $self->o('use_existing_short_read_dir'),
      },

      -flow_into => {
        '1->A' => { 'fan_short_read_download' => { 'inputfile' => '#rnaseq_summary_file#', 'input_dir' => '#short_read_dir#' } },
        'A->1' => ['download_long_read_csv'],
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
          ' -dbname '.'#core_db#'.
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
          ' -dbname '.'#core_db#',
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
      -max_retry_count => 0,
      -flow_into       => {
        1 => ['check_load_meta_info'],
      },
    },
    {
      -logic_name => 'check_load_meta_info',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd                     => 'if [ -e "#rnaseq_summary_file#" ] || [ -e "#long_read_summary_file#" ]; then exit 0; else  exit 42;fi',
        return_codes_2_branches => { '42' => 2 },
      },
      -flow_into => {
        1 => ['anno_load_meta_info'],
        2 => ['braker_load_meta_info'],
      },
      -rc_name => 'default',
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
            '(1, "species.strain_group", "#species_strain_group#"),' .
            '(1, "species.production_name", "#production_name#"),' .
            '(1, "strain.type", "#strain_type#"),' .
            '(1, "repeat.analysis", "repeatdetector"),' .
            '(1, "repeat.analysis", "dust"),' .
            '(1, "repeat.analysis", "trf"),' .
            '(1, "genebuild.initial_release_date", NULL),' .
            '(1, "genebuild.projection_source_db", NULL),' .
            '(1, "genebuild.id", ' . $self->o('genebuilder_id') . '),' .
            '(1, "genebuild.method", "anno")'
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
      -logic_name => 'braker_load_meta_info',
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
            '(1, "species.stable_id_prefix", "BRAKER#species_prefix#"),' .
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
            '(1, "genebuild.method", "braker")'
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
        cmd                     => 'if [ -e "#rnaseq_summary_file#" ] || [ -e "#long_read_summary_file#" ]; then exit 0; else exit 42;fi',
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
        cmd => 'python3.7 ' . catfile( $self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'genebuild', 'gbiab', 'gbiab.py' ) . ' #anno_commandline#',
      },
      -rc_name         => 'anno',
      -max_retry_count => 0,
      -flow_into       => {
        1 => ['update_biotypes_and_analyses']
      },
    },

    {
      -logic_name => 'update_biotypes_and_analyses',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => '#core_db#',
        sql     => [
          'UPDATE analysis SET module=NULL',
          'UPDATE gene SET biotype = "protein_coding" WHERE biotype = "gbiab_protein_coding"',
          'UPDATE gene SET biotype = "lncRNA" WHERE biotype = "gbiab_lncRNA"',
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
        cmd                     => 'if [ -e "#rnaseq_summary_file#" ] || [ -e "#long_read_summary_file#" ]; then exit 0; else  exit 42;fi',
        return_codes_2_branches => { '42' => 2 },
      },
      -flow_into => {
        1 => ['anno_run_stable_ids'],
        2 => ['braker_run_stable_ids'],
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
      -logic_name => 'braker_run_stable_ids',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::SetStableIDs',
      -parameters => {
        enscode_root_dir => $self->o('enscode_root_dir'),
        mapping_required => 0,
        target_db        => '#core_db#',
        id_start         => 'BRAKER#species_prefix#' . '#stable_id_start#',
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
        1 => ['fan_busco_output'],
      },
      -rc_name => 'default_registry',
    },
    {
      -logic_name => 'fan_busco_output',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        #cmd                     => 'if '.$self->o('run_braker').'==1] ; then exit 0; else exit 42;fi',
        cmd                     => 'if [ -e "#rnaseq_summary_file#" ] || [ -e "#long_read_summary_file#" ]; then exit 0; else  exit 42;fi',
        return_codes_2_branches => { '42' => 2 },
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['create_busco_dirs'],
        2 => ['run_agat_protein_file'],
      },
    },
    {
      -logic_name => 'run_agat_protein_file',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',

      -parameters => {
        cmd => 'singularity exec --bind #output_path#/:/data:rw  ' . $self->o('agat_singularity_image') . ' agat_sp_extract_sequences.pl --gff /data/braker/braker.gtf -f  #output_path#/#species_name#_softmasked_toplevel.fa -p  -o  #output_path#/braker/braker_proteins.fa;',
      },
      -rc_name         => 'braker',
      -max_retry_count => 0,
      -flow_into       => {
        1 => ['run_busco_braker'],
      },
    },
    {
      -logic_name => 'run_busco_braker',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',

      -parameters => {
        cmd => 'cd #output_path#/;' .
          'singularity exec ' . $self->o('busco_singularity_image') . ' busco -i #output_path#/braker/braker_proteins.fa -m prot -l #busco_group# -o output_busco_#assembly_accession#  ;',
      },
      -rc_name   => 'braker',
      -flow_into => {
        1 => ['update_assembly_registry_status'],
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
        cmd => 'cd #output_path#; singularity exec ' .
          $self->o('busco_singularity_image') .
          ' busco -f -i #output_path#/busco_score_data/canonical_proteins.fa -m prot -l ' . '#busco_group#' .
          ' -o busco_score_output;',
      },
      -rc_name   => 'busco',
      -flow_into => { 1 => ['fan_otherfeatures_db'] },
    },
    {
      -logic_name => 'fan_otherfeatures_db',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd                     => 'if [ -e "#rnaseq_summary_file#" ] || [ -e "#long_read_summary_file#" ]; then exit 0; else  exit 42;fi',
        return_codes_2_branches => { '42' => 2 },
      },
      -flow_into => {
        1 => ['create_otherfeatures_db'],

        2 => ['update_assembly_registry_status'],
      },
      -rc_name => 'default',
    },
    {
      # Creates a reference db for each species
      -logic_name => 'create_otherfeatures_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        'source_db'        => '#core_db#',
        'target_db'        => '#otherfeatures_db#',
        'enscode_root_dir' => $self->o('enscode_root_dir'),
        'create_type'      => 'clone',
      },
      -rc_name => 'default',

      -flow_into => {
        1 => ['populate_production_tables_otherfeatures'],
      },
    },
    {
      # Load production tables into each reference
      -logic_name => 'populate_production_tables_otherfeatures',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HivePopulateProductionTables',
      -parameters => {
        'target_db'        => '#otherfeatures_db#',
        'output_path'      => '#output_path#',
        'enscode_root_dir' => $self->o('enscode_root_dir'),
        'production_db'    => $self->o('production_db'),
      },
      -rc_name => 'default',

      -flow_into => {
        1 => ['run_braker_ep_mode_otherfeatures'],
      },
    },
    {
      -logic_name => 'run_braker_ep_mode_otherfeatures',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',

      -parameters => {
        cmd => 'cp #output_path#/red_output/mask_output/#species_name#_reheadered_toplevel.msk #output_path#/#species_name#_softmasked_toplevel.fa;' .
          'mkdir #output_path#/prothint;' .
          'cd #output_path#/prothint;' .
          'singularity exec --bind #output_path#/prothint/:/data:rw  ' . $self->o('braker_singularity_image') . ' prothint.py #output_path#/#species_name#_softmasked_toplevel.fa #protein_file# ;' .
          'cd #output_path#/;' .
          'singularity exec --bind #output_path#/:/data:rw  ' . $self->o('braker_singularity_image') . ' braker.pl --genome=#species_name#_softmasked_toplevel.fa --softmasking  --hints=/data/prothint/prothint_augustus.gff --prothints=/data/prothint/prothint.gff --evidence=/data/prothint/evidence.gff --epmode --species=#assembly_accession#_#species_name# --AUGUSTUS_CONFIG_PATH=' . $self->o('augustus_config_path') . ' --cores ' . $self->o('cores') . ';',
      },
      -rc_name         => 'braker',
      -max_retry_count => 0,
      -flow_into       => {
        1 => ['load_gtf_file_in_otherfeatures_db'],
      },
    },
    {
      -logic_name => 'load_gtf_file_in_otherfeatures_db',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',

      -parameters => {
        cmd => 'perl ' . catfile( $self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'genebuild', 'braker', 'parse_gtf.pl' ) .
          ' -dnahost ' . $self->o('otherfeatures_db_host') .
          ' -dnauser ' . $self->o('user_r') .
          ' -dnaport ' . $self->o('otherfeatures_db_port') .
          ' -dnadbname #otherfeatures_dbname#' .
          ' -host ' . $self->o('otherfeatures_db_host') .
          ' -user ' . $self->o('user') .
          ' -pass ' . $self->o('password') .
          ' -port ' . $self->o('otherfeatures_db_port') .
          ' -dbname #otherfeatures_dbname#' .
          ' -write' .
          ' -file #output_path#/braker/braker.gtf',
      },
      -rc_name         => 'default',
      -max_retry_count => 0,
      -flow_into       => {
        1 => ['update_otherfeatures_db'],
      },
    },
    {
      -logic_name => 'update_otherfeatures_db',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => '#otherfeatures_db#',
        sql     => [
          'DELETE FROM analysis_description WHERE analysis_id IN (SELECT analysis_id FROM analysis WHERE logic_name IN' .
            ' ("dust","repeatdetector","trf","cpg","eponine"))',
          'DELETE FROM analysis WHERE logic_name IN' .
            ' ("dust","repeatdetector","trf","cpg","eponine")',
          'DELETE FROM meta WHERE meta_key LIKE "%.level"',
	  'DELETE FROM meta WHERE meta_key LIKE "sample.%"',
          'DELETE FROM meta WHERE meta_key LIKE "assembly.web_accession%"',
          'DELETE FROM meta WHERE meta_key LIKE "removed_evidence_flag.%"',
          'DELETE FROM meta WHERE meta_key LIKE "marker.%"',
          'DELETE FROM meta WHERE meta_key IN' .
            ' ("repeat.analysis","genebuild.method","genebuild.last_geneset_update","genebuild.projection_source_db","genebuild.start_date")',
          'INSERT INTO meta (species_id,meta_key,meta_value) VALUES (1,"genebuild.last_otherfeatures_update",NOW())',
	  'UPDATE meta set meta_value="BRAKER#species_prefix#" where meta_key="species.stable_id_prefix"',
          'UPDATE transcript JOIN transcript_supporting_feature USING(transcript_id)'.
              ' JOIN dna_align_feature ON feature_id = dna_align_feature_id SET stable_id = hit_name',  
      ],
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['otherfeatures_set_meta_coords'],
      },
    },
    {
      -logic_name => 'otherfeatures_set_meta_coords',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . $self->o('meta_coord_script') .
          ' -user ' . $self->o('user') .
          ' -pass ' . $self->o('password') .
          ' -host ' . $self->o( 'otherfeatures_db', '-host' ) .
          ' -port ' . $self->o( 'otherfeatures_db', '-port' ) .
          ' -dbpattern ' . '#otherfeatures_dbname#'
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['otherfeatures_set_meta_levels'],
      },
    },

    {
      -logic_name => 'otherfeatures_set_meta_levels',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . $self->o('meta_levels_script') .
          ' -user ' . $self->o('user') .
          ' -pass ' . $self->o('password') .
          ' -host ' . $self->o( 'otherfeatures_db', '-host' ) .
          ' -port ' . $self->o( 'otherfeatures_db', '-port' ) .
          ' -dbname ' . '#otherfeatures_dbname#'
      },
      -rc_name   => 'default',
      -flow_into => { 1 => ['otherfeatures_set_frameshift_introns'] },
    },

    {
      -logic_name => 'otherfeatures_set_frameshift_introns',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . $self->o('frameshift_attrib_script') .
          ' -user ' . $self->o('user') .
          ' -pass ' . $self->o('password') .
          ' -host ' . $self->o( 'otherfeatures_db', '-host' ) .
          ' -port ' . $self->o( 'otherfeatures_db', '-port' ) .
          ' -dbpattern ' . '#otherfeatures_dbname#'
      },
      -rc_name   => '10GB',
      -flow_into => { 1 => ['otherfeatures_set_canonical_transcripts'] },
    },
    {
      -logic_name => 'otherfeatures_set_canonical_transcripts',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . $self->o('select_canonical_script') .
          ' -dbuser ' . $self->o('user') .
          ' -dbpass ' . $self->o('password') .
          ' -dbhost ' . $self->o( 'otherfeatures_db', '-host' ) .
          ' -dbport ' . $self->o( 'otherfeatures_db', '-port' ) .
          ' -dbname ' . '#otherfeatures_dbname#' .
          ' -dnadbuser ' . $self->o('user_r') .
          ' -dnadbhost ' . $self->o( 'core_db', '-host' ) .
          ' -dnadbport ' . $self->o( 'core_db', '-port' ) .
          ' -dnadbname ' . '#core_dbname#' .
          ' -coord toplevel -write'
      },
      -rc_name   => '10GB',
      -flow_into => { 1 => ['otherfeatures_null_columns'] },
    },

    {
      -logic_name => 'otherfeatures_null_columns',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => '#otherfeatures_db#',
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
        1 => ['otherfeatures_braker_run_stable_ids'],
      },
    },

    {
      -logic_name => 'otherfeatures_braker_run_stable_ids',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::SetStableIDs',
      -parameters => {
        enscode_root_dir => $self->o('enscode_root_dir'),
        mapping_required => 0,
        target_db        => '#otherfeatures_db#',
        id_start         => 'BRAKER#species_prefix#' . '#stable_id_start#',
        output_path      => '#output_path#',
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['load_external_db_ids_and_optimise_otherfeatures'],
      },
    },
    {
      -logic_name => 'load_external_db_ids_and_optimise_otherfeatures',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . $self->o('load_optimise_script') .
          ' -output_path ' . catdir( '#output_path#', 'optimise_otherfeatures' ) .
          ' -uniprot_filename ' . $self->o('protein_entry_loc') .
          ' -dbuser ' . $self->o('user') .
          ' -dbpass ' . $self->o('password') .
          ' -dbport ' . $self->o( 'otherfeatures_db', '-port' ) .
          ' -dbhost ' . $self->o( 'otherfeatures_db', '-host' ) .
          ' -dbname ' . '#otherfeatures_dbname#' .
          ' -prod_dbuser ' . $self->o('user_r') .
          ' -prod_dbhost ' . $self->o( 'production_db', '-host' ) .
          ' -prod_dbname ' . $self->o( 'production_db', '-dbname' ) .
          ' -prod_dbport ' . $self->o( 'production_db', '-port' ) .
          ' -verbose'
      },
      -max_retry_count => 0,
      -rc_name         => '4GB',
      -flow_into       => {
        1 => ['otherfeatures_final_cleaning'],
      },
    },

    {
      -logic_name => 'otherfeatures_final_cleaning',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => '#otherfeatures_db#',
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
        1 => ['otherfeatures_add_placeholder_sample_location'],
      },

    },
    {
      -logic_name => 'otherfeatures_add_placeholder_sample_location',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAddPlaceholderLocation',
      -parameters => {
        input_db => '#otherfeatures_db#',
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['otherfeatures_populate_analysis_descriptions'],
      },
    },

    {
      -logic_name => 'otherfeatures_populate_analysis_descriptions',
      -module     => 'Bio::EnsEMBL::Production::Pipeline::ProductionDBSync::PopulateAnalysisDescription',
      -parameters => {
        species => '#production_name#',
        group   => 'otherfeatures',
      },
      -flow_into => {
        1 => ['otherfeatures_sanity_checks'],
      },
      -rc_name => 'default_registry',
    },
    {
      -logic_name => 'otherfeatures_sanity_checks',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
      -parameters => {
        target_db                  => '#otherfeatures_db#',
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
        input_db => '#otherfeatures_db#',
        species  => '#species_name#',
        group    => 'otherfeatures_handover',
      },
      -max_retry_count => 0,
      -rc_name         => '4GB',
      -flow_into       => { 1 => ['update_assembly_registry_status'], },
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
     },    

    {
      -logic_name => 'run_anno_softmasking',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'python3.7 ' . catfile( $self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'genebuild', 'gbiab', 'gbiab.py' ) . ' #anno_red_commandline#;' .
          'cp #output_path#/red_output/mask_output/#species_name#_reheadered_toplevel.msk #output_path#/#species_name#_softmasked_toplevel.fa',
      },
      -rc_name         => 'anno',
      -max_retry_count => 0,
      -flow_into       => {
        1 => ['run_braker_ep_mode'],
      },
    },
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
      -rc_name         => 'braker',
      -max_retry_count => 0,
      -flow_into       => {
        1 => ['load_gtf_file'],
      },
    },
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
  ];
}

sub resource_classes {
  my $self = shift;

  return {
    'default_registry' => { LSF => [ $self->lsf_resource_builder( 'production', 900, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ), '-reg_conf ' . $self->default_options->{'registry_file'} ] },
    'anno'             => { LSF => $self->lsf_resource_builder( 'production', 50000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], undef, $self->default_options->{'num_threads'} ) },
    'braker'           => { LSF => $self->lsf_resource_builder( 'production', 32000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], undef, $self->default_options->{'cores'} ) },
    'busco'            => { LSF => $self->lsf_resource_builder( 'production', 32000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], undef, $self->default_options->{'cores'} ) },
    '1GB'              => { LSF => $self->lsf_resource_builder( 'production', 1000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    '2GB'              => { LSF => $self->lsf_resource_builder( 'production', 2000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    '3GB'              => { LSF => $self->lsf_resource_builder( 'production', 3000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    '4GB'              => { LSF => $self->lsf_resource_builder( 'production', 4000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'refseq_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    '5GB'              => { LSF => $self->lsf_resource_builder( 'production', 5000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    '8GB'              => { LSF => $self->lsf_resource_builder( 'production', 8000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    '10GB'             => { LSF => $self->lsf_resource_builder( 'production', 10000, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
    'default'          => { LSF => $self->lsf_resource_builder( 'production', 900, [ $self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'} ], [ $self->default_options->{'num_tokens'} ] ) },
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
