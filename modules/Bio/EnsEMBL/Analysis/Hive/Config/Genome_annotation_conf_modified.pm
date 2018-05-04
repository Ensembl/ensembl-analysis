=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

package Genome_annotation_conf;

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
    'pipeline_name'             => '', # What you want hive to call the pipeline, not the db name itself
    'user_r'                    => '', # read only db user
    'user'                      => '', # write db user
    'password'                  => '', # password for write db user
    'pipe_db_server'            => '', # host for pipe db
    'databases_server'          => '', # host for general output dbs
    'dna_db_server'             => '', # host for dna db
    'pipe_db_port'              => '', # port for pipeline host
    'databases_port'            => '', # port for general output db host
    'dna_db_port'               => '', # prot for dna db host
    'repbase_logic_name'        => '', # repbase logic name i.e. repeatmask_repbase_XXXX, ONLY FILL THE XXXX BIT HERE!!! e.g primates
    'repbase_library'           => '', # repbase library name, this is the actual repeat repbase library to use, e.g. "Mus musculus"
    'repeatmodeler_library'     => '', # This should be the path to a custom repeat library, leave blank if none exists
    'release_number'            => '' || $self->o('ensembl_release'),
    'species_name'              => '', # e.g. mus_musculus
    'production_name'           => '', # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
    'taxon_id'                  => '', # should be in the assembly report file
    'uniprot_set'               => '', # e.g. mammals_basic, check UniProtCladeDownloadStatic.pm module in hive config dir for suitable set,
    'output_path'               => '', # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses
    'wgs_id'                    => '', # Can be found in assembly report file on ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/
    'assembly_name'             => '', # Name (as it appears in the assembly report file)
    'assembly_accession'        => '', # Versioned GCA assembly accession, e.g. GCA_001857705.1
    'assembly_refseq_accession' => '', # Versioned GCF accession, e.g. GCF_001857705.1
    'stable_id_prefix'          => '', # e.g. ENSPTR. When running a new annotation look up prefix in the assembly registry db
    'stable_id_start'           => '0', # When mapping is not required this is usually set to 0
    'load_toplevel_only'        => 1, # This will not load the assembly info and will instead take any chromosomes, unplaced and unlocalised scaffolds directly in the DNA table
    'skip_projection'           => 0, # Will skip projection process if 1
    'skip_rnaseq'               => 0, # Will skip rnaseq analyses if 1
    'skip_ncrna'                => 0, # Will skip ncrna process if 1
    'skip_cleaning'             => 0, # Will skip the cleaning phase, will keep more genes/transcripts but some lower quality models may be kept
    'mapping_required'          => 0, # If set to 1 this will run stable_id mapping sometime in the future. At the moment it does nothing
    'mapping_db'                => undef, # Tied to mapping_required being set to 1, we should have a mapping db defined in this case, leave undef for now
    'uniprot_db_dir'            => 'uniprot_2018_01', # What UniProt data dir to use for various analyses
    'vertrna_version'           => 134, # The version of VertRNA to use, should correspond to a numbered dir in VertRNA dir
    'mirBase_fasta'             => 'all_mirnas.fa', # What mirBase file to use. It is currently best to use on with the most appropriate set for your species
    'rfc_scaler'                => 'filter_dafs_rfc_scaler_human.pkl',
    'rfc_model'                 => 'filter_dafs_rfc_model_human.pkl',
    'ig_tr_fasta_file'          => 'human_ig_tr.fa', # What IMGT fasta file to use. File should contain protein segments with appropriate headers
    'mt_accession'              => undef, # This should be set to undef unless you know what you are doing

########################
# Pipe and ref db info
########################

    'projection_source_db_name'    => 'homo_sapiens_core_91_38', # This is generally a pre-existing db, like the current human/mouse core for example
    'projection_source_db_server'  => 'mysql-ensembl-mirror',
    'projection_source_db_port'    => 4240,

    # The following might not be known in advance, since the come from other pipelines
    # These values can be replaced in the analysis_base table if they're not known yet
    # If they are not needed (i.e. no projection or rnaseq) then leave them as is
    'projection_lastz_db_name'     => 'PROJECTION_LASTZ_DBNAME',
    'projection_lastz_db_server'   => 'PROJECTION_LASTZ_SERVER',
    'projection_lastz_db_port'     => 'PROJECTION_LASTZ_PORT',
    'rnaseq_refine_db_name'        => 'RNASEQ_REFINE_DBNAME',
    'rnaseq_refine_db_server'      => 'RNASEQ_REFINE_SERVER',
    'rnaseq_refine_db_port'        => 'RNASEQ_REFINE_PORT',
    'rnaseq_blast_db_name'         => 'RNASEQ_BLAST_DBNAME',
    'rnaseq_blast_db_server'       => 'RNASEQ_BLAST_SERVER',
    'rnaseq_blast_db_port'         => 'RNASEQ_BLAST_PORT',

    'provider_name'                => 'Ensembl',
    'provider_url'                 => 'www.ensembl.org',

    'pipe_db_name'                  => $self->o('dbowner').'_'.$self->o('production_name').'_pipe_'.$self->o('release_number'),
    'dna_db_name'                   => $self->o('dbowner').'_'.$self->o('production_name').'_core_'.$self->o('release_number'),

    'reference_db_name'            => $self->o('dna_db_name'),
    'reference_db_server'          => $self->o('dna_db_server'),
    'reference_db_port'            => $self->o('dna_db_port'),

    'cdna_db_server'               => $self->o('databases_server'),
    'cdna_db_port'                 => $self->o('databases_port'),

    cdna2genome_db_server          => $self->o('databases_server'),
    cdna2genome_db_port            => $self->o('databases_port'),

    'genblast_db_server'           => $self->o('databases_server'),
    'genblast_db_port'             => $self->o('databases_port'),

    'genblast_select_db_server'    => $self->o('databases_server'),
    'genblast_select_db_port'      => $self->o('databases_port'),

    'ig_tr_db_server'              => $self->o('databases_server'),
    'ig_tr_db_port'                => $self->o('databases_port'),

    'genewise_db_server'           => $self->o('databases_server'),
    'genewise_db_port'             => $self->o('databases_port'),

    'projection_coding_db_server'  => $self->o('databases_server'),
    'projection_coding_db_port'    => $self->o('databases_port'),

    'projection_realign_db_server' => $self->o('databases_server'),
    'projection_realign_db_port'   => $self->o('databases_port'),

    'projection_lincrna_db_server' => $self->o('databases_server'),
    'projection_lincrna_db_port'   => $self->o('databases_port'),

    'projection_pseudogene_db_server' => $self->o('databases_server'),
    'projection_pseudogene_db_port'   => $self->o('databases_port'),

    'rnaseq_for_layer_db_server'   => $self->o('databases_server'),
    'rnaseq_for_layer_db_port'     => $self->o('databases_port'),

    'rnaseq_final_db_server'       => $self->o('databases_server'),
    'rnaseq_final_db_port'         => $self->o('databases_port'),

    'layering_db_server'           => $self->o('databases_server'),
    'layering_db_port'             => $self->o('databases_port'),

    'utr_db_server'                => $self->o('databases_server'),
    'utr_db_port'                  => $self->o('databases_port'),

    'genebuilder_db_server'        => $self->o('databases_server'),
    'genebuilder_db_port'          => $self->o('databases_port'),

    'pseudogene_db_server'         => $self->o('databases_server'),
    'pseudogene_db_port'           => $self->o('databases_port'),

    'ncrna_db_server'              => $self->o('databases_server'),
    'ncrna_db_port'                => $self->o('databases_port'),

    'final_geneset_db_server'      => $self->o('databases_server'),
    'final_geneset_db_port'        => $self->o('databases_port'),

    'refseq_db_server'             => $self->o('databases_server'),
    'refseq_db_port'               => $self->o('databases_port'),

    'killlist_db_server'           => $self->o('databases_server'),
    'killlist_db_port'             => $self->o('databases_port'),

    'otherfeatures_db_server'      => $self->o('databases_server'),
    'otherfeatures_db_port'        => $self->o('databases_port'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'ensembl_release'              => $ENV{ENSEMBL_RELEASE}, # this is the current release version on staging to be able to get the correct database
    'staging_1_db_server'          => 'mysql-ens-sta-1',
    'staging_1_port'               => 4519,


    databases_to_delete => ['reference_db', 'cdna_db', 'genblast_db', 'genewise_db', 'projection_coding_db', 'projection_realign_db', 'layering_db', 'utr_db', 'genebuilder_db', 'pseudogene_db', 'ncrna_db', 'final_geneset_db', 'refseq_db', 'cdna2genome_db'],

########################
# BLAST db paths
########################
    'base_blast_db_path'        => $ENV{BLASTDB_DIR},
    'uniprot_entry_loc'         => catfile($self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_db_dir'), 'entry_loc'),
    'uniprot_blast_db_path'     => catfile($self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_db_dir'), 'uniprot_vertebrate'),
    'vertrna_blast_db_path'     => catfile($self->o('base_blast_db_path'), 'vertrna', $self->o('vertrna_version'), 'embl_vertrna-1'),
    'unigene_blast_db_path'     => catfile($self->o('base_blast_db_path'), 'unigene', 'unigene'),
    'ncrna_blast_path'          => catfile($self->o('base_blast_db_path'), 'ncrna', 'ncrna_2016_05'),
    'mirna_blast_path'          => catfile($self->o('base_blast_db_path'), 'ncrna', 'mirbase_22'),
    'ig_tr_blast_path'          => catfile($self->o('base_blast_db_path'), 'ig_tr_genes'),

######################################################
#
# Mostly constant settings
#
######################################################

    'genome_file'               => catfile($self->o('output_path'), 'genome_dumps', $self->o('species_name').'_softmasked_toplevel.fa'),
    'primary_assembly_dir_name' => 'Primary_Assembly',
    'refseq_cdna_calculate_coverage_and_pid' => '0',
    'contigs_source'            => 'ena',



    'layering_input_gene_dbs' => [
                                   $self->o('genblast_db'),
                                   $self->o('rnaseq_for_layer_db'),
                                   $self->o('projection_coding_db'),
                                   $self->o('ig_tr_db'),
                                 ],

    utr_donor_dbs => [
      $self->o('cdna_db'),
      $self->o('rnaseq_for_layer_db'),
    ],

    utr_acceptor_dbs => [
      $self->o('layering_db'),
    ],

    'utr_biotype_priorities'  => {
                                   'rnaseq' => 2,
                                   'cdna' => 1,
                                 },

    'cleaning_blessed_biotypes' => {
                                     'pseudogene' => 1,
                                     'processed_pseudogene' => 1,
                                     'IG_C_gene' => 1,
                                     'IG_V_gene' => 1,
                                     'TR_C_gene' => 1,
                                     'TR_D_gene' => 1,
                                     'TR_V_gene' => 1,
                                   },

    'min_toplevel_slice_length'   => 250,

    'repeat_logic_names'          => ['repeatmask_repbase_'.$self->o('repbase_logic_name'),'dust'],
    'homology_models_path'        => catdir($self->o('output_path'),'homology_models'),

    targetted_path => catdir($self->o('output_path'), 'targetted'),
    cdna_file      => catfile($self->o('targetted_path'), 'cdnas'),
    annotation_file => $self->o('cdna_file').'.annotation',

    ensembl_analysis_script       => catdir($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts'),
    remove_duplicates_script_path => catfile($self->o('ensembl_analysis_script'), 'find_and_remove_duplicates.pl'),
    load_optimise_script          => catfile($self->o('ensembl_analysis_script'), 'load_external_db_ids_and_optimize_af.pl'),
    prepare_cdnas_script          => catfile($self->o('ensembl_analysis_script'), 'genebuild', 'prepare_cdnas.pl'),
    load_fasta_script_path        => catfile($self->o('ensembl_analysis_script'), 'genebuild', 'load_fasta_to_db_table.pl'),
    loading_report_script         => catfile($self->o('ensembl_analysis_script'), 'genebuild', 'report_genome_prep_stats.pl'),
    refseq_synonyms_script_path   => catfile($self->o('ensembl_analysis_script'), 'refseq', 'load_refseq_synonyms.pl'),
    refseq_import_script_path     => catfile($self->o('ensembl_analysis_script'), 'refseq', 'parse_ncbi_gff3.pl'),
    'mirna_analysis_script'       => catdir($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'genebuild', 'sncrna'),

    ensembl_misc_script        => catdir($self->o('enscode_root_dir'), 'ensembl', 'misc-scripts'),
    repeat_types_script        => catfile($self->o('ensembl_misc_script'), 'repeats', 'repeat-types.pl'),
    meta_coord_script          => catfile($self->o('ensembl_misc_script'), 'meta_coord', 'update_meta_coord.pl'),
    meta_levels_script         => catfile($self->o('ensembl_misc_script'), 'meta_levels.pl'),
    frameshift_attrib_script   => catfile($self->o('ensembl_misc_script'), 'frameshift_transcript_attribs.pl'),
    select_canonical_script    => catfile($self->o('ensembl_misc_script'), 'select_canonical_transcripts.pl'),
    generate_stable_ids_script => catfile($self->o('ensembl_misc_script'), 'generate_stable_ids.pl'),

    rnaseq_stable_id_file => catfile($self->o('output_path'), 'rnaseq_stable_ids.sql'),
    rnaseq_daf_introns_file => catfile($self->o('output_path'), 'rnaseq_daf_introns.dat'),

########################
# Extra db settings
########################

    'num_tokens' => 10,

########################
# Executable paths
########################
    'blast_type' => 'ncbi', # It can be 'ncbi', 'wu', or 'legacy_ncbi'
    'dust_path' => catfile($self->o('binary_base'), 'dustmasker'),
    'trf_path' => catfile($self->o('binary_base'), 'trf'),
    'eponine_java_path' => catfile($self->o('binary_base'), 'java'),
    'eponine_jar_path' => catfile($self->o('software_base_path'), 'opt', 'eponine', 'libexec', 'eponine-scan.jar'),
    'cpg_path' => catfile($self->o('binary_base'), 'cpg_lh'),
    'trnascan_path' => catfile($self->o('binary_base'), 'tRNAscan-SE'),
    'repeatmasker_path' => catfile($self->o('binary_base'), 'RepeatMasker'),
    'genscan_path' => catfile($self->o('binary_base'), 'genscan'),
    'genscan_matrix_path' => catfile($self->o('software_base_path'), 'share', 'HumanIso.smat'),
    'uniprot_blast_exe_path' => catfile($self->o('binary_base'), 'blastp'),
    'blastn_exe_path' => catfile($self->o('binary_base'), 'blastn'),
    'vertrna_blast_exe_path' => catfile($self->o('binary_base'), 'tblastn'),
    'unigene_blast_exe_path' => catfile($self->o('binary_base'), 'tblastn'),
    genewise_path => catfile($self->o('binary_base'), 'genewise'),
    'exonerate_path'         => catfile($self->o('software_base_path'), 'opt', 'exonerate09', 'bin', 'exonerate'),
    'cmsearch_exe_path'    => catfile($self->o('software_base_path'), 'opt', 'infernal10', 'bin', 'cmsearch'),
    indicate_path  => catfile($self->o('binary_base'), 'indicate'),
    pmatch_path  => catfile($self->o('binary_base'), 'pmatch'),
    exonerate_annotation => catfile($self->o('binary_base'), 'exonerate'),

    'uniprot_genblast_batch_size' => 15,
    'uniprot_table_name'          => 'uniprot_sequences',

    'genblast_path'     => catfile($self->o('binary_base'), 'genblast'),
    'genblast_eval'     => $self->o('blast_type') eq 'wu' ? '1e-20' : '1e-1',
    'genblast_cov'      => '0.5',
    'genblast_pid'      => '50',
    'genblast_max_rank' => '5',
    'genblast_flag_small_introns' => 1,
    'genblast_flag_subpar_models' => 1,

    'ig_tr_table_name'    => 'ig_tr_sequences',
    'ig_tr_genblast_cov'  => '0.8',
    'ig_tr_genblast_pid'  => '70',
    'ig_tr_genblast_eval' => '1',
    'ig_tr_genblast_max_rank' => '10',
    'ig_tr_batch_size'    => 10,

    'exonerate_cdna_pid' => '95', # Cut-off for percent id
    'exonerate_cdna_cov' => '50', # Cut-off for coverage

    'cdna_selection_pid' => '97', # Cut-off for percent id for selecting the cDNAs
    'cdna_selection_cov' => '90', # Cut-off for coverage for selecting the cDNAs

# Best targetted stuff
    exonerate_logic_name => 'exonerate',
    ncbi_query => 'txid'.$self->o('taxon_id').'[Organism:noexp]+AND+biomol_mrna[PROP]',

    cdna_table_name    => 'cdna_sequences',
    target_exonerate_calculate_coverage_and_pid => 0,
    exonerate_protein_pid => 95,
    exonerate_protein_cov => 50,
    cdna2genome_region_padding => 2000,
    exonerate_max_intron => 200000,

    best_targetted_min_coverage => 50, # This is to avoid having models based on fragment alignment and low identity
    best_targetted_min_identity => 50, # This is to avoid having models based on fragment alignment and low identity


# Max internal stops for projected transcripts
    'projection_pid'                        => '50',
    'projection_cov'                        => '50',
    'projection_max_internal_stops'         => '1',
    'projection_calculate_coverage_and_pid' => '1',

    'projection_lincrna_percent_id'         => 90,
    'projection_lincrna_coverage'           => 90,
    'projection_pseudogene_percent_id'      => 60,
    'projection_pseudogene_coverage'        => 75,
    'projection_ig_tr_percent_id'           => 70,
    'projection_ig_tr_coverage'             => 90,
    'projection_exonerate_padding'          => 5000,

    'realign_table_name'                    => 'projection_source_sequences',
    'max_projection_structural_issues'      => 1,

## Add in genewise path and put in matching code
    'genewise_pid'                        => '50',
    'genewise_cov'                        => '50',
    'genewise_region_padding'             => '50000',
    'genewise_calculate_coverage_and_pid' => '1',

########################
# Misc setup info
########################
    'repeatmasker_engine' => 'crossmatch',

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# No option below this mark should be modified
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
########################################################
# URLs for retrieving the INSDC contigs and RefSeq files
########################################################
    'ncbi_base_ftp'           => 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all',
    'insdc_base_ftp'          => $self->o('ncbi_base_ftp').'/#expr(substr(#assembly_accession#, 0, 3))expr#/#expr(substr(#assembly_accession#, 4, 3))expr#/#expr(substr(#assembly_accession#, 7, 3))expr#/#expr(substr(#assembly_accession#, 10, 3))expr#/#assembly_accession#_#assembly_name#',
    'assembly_ftp_path'       => $self->o('insdc_base_ftp'),
    'refseq_base_ftp'         => $self->o('ncbi_base_ftp').'/#expr(substr(#assembly_refseq_accession#, 0, 3))expr#/#expr(substr(#assembly_refseq_accession#, 4, 3))expr#/#expr(substr(#assembly_refseq_accession#, 7, 3))expr#/#expr(substr(#assembly_refseq_accession#, 10, 3))expr#/#assembly_refseq_accession#_#assembly_name#',
    'refseq_import_ftp_path'  => $self->o('refseq_base_ftp').'/#assembly_refseq_accession#_#assembly_name#_genomic.gff.gz',
    'refseq_mrna_ftp_path'    => $self->o('refseq_base_ftp').'/#assembly_refseq_accession#_#assembly_name#_rna.fna.gz',

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

    'cdna_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_cdna_'.$self->o('release_number'),
      -host   => $self->o('cdna_db_server'),
      -port   => $self->o('cdna_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },


    'genblast_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_genblast_'.$self->o('release_number'),
      -host   => $self->o('genblast_db_server'),
      -port   => $self->o('genblast_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },


    'genblast_select_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_gensel_'.$self->o('release_number'),
      -host   => $self->o('genblast_select_db_server'),
      -port   => $self->o('genblast_select_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },


    'ig_tr_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_igtr_'.$self->o('release_number'),
      -host   => $self->o('ig_tr_db_server'),
      -port   => $self->o('ig_tr_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    cdna2genome_db => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_cdna2genome_'.$self->o('release_number'),
      -host   => $self->o('cdna2genome_db_server'),
      -port   => $self->o('cdna2genome_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'genewise_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_genewise_'.$self->o('release_number'),
      -host   => $self->o('genewise_db_server'),
      -port   => $self->o('genewise_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'best_targeted_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_bt_'.$self->o('release_number'),
      -host   => $self->o('genewise_db_server'),
      -port   => $self->o('genewise_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'projection_coding_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_proj_coding_'.$self->o('release_number'),
      -host   => $self->o('projection_coding_db_server'),
      -port   => $self->o('projection_coding_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'projection_realign_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_realign_'.$self->o('release_number'),
      -host   => $self->o('projection_realign_db_server'),
      -port   => $self->o('projection_realign_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'projection_lincrna_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_proj_linc_'.$self->o('release_number'),
      -host   => $self->o('projection_lincrna_db_server'),
      -port   => $self->o('projection_lincrna_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'projection_pseudogene_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_proj_pseudo_'.$self->o('release_number'),
      -host   => $self->o('projection_pseudogene_db_server'),
      -port   => $self->o('projection_pseudogene_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'projection_source_db' => {
      -dbname => $self->o('projection_source_db_name'),
      -host   => $self->o('projection_source_db_server'),
      -port   => $self->o('projection_source_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -driver => $self->o('hive_driver'),
    },

    'projection_lastz_db' => {
      -dbname => $self->o('projection_lastz_db_name'),
      -host   => $self->o('projection_lastz_db_server'),
      -port   => $self->o('projection_lastz_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -driver => $self->o('hive_driver'),
    },

    'rnaseq_for_layer_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_rnaseq_layer_'.$self->o('release_number'),
      -host   => $self->o('rnaseq_for_layer_db_server'),
      -port   => $self->o('rnaseq_for_layer_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'rnaseq_final_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_rnaseq_final_'.$self->o('release_number'),
      -host   => $self->o('rnaseq_final_db_server'),
      -port   => $self->o('rnaseq_final_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'rnaseq_refine_db' => {
      -dbname => $self->o('rnaseq_refine_db_name'),
      -host   => $self->o('rnaseq_refine_db_server'),
      -port   => $self->o('rnaseq_refine_db_port'),
      -user   => $self->o('user_r'),
      -driver => $self->o('hive_driver'),
    },

    'rnaseq_blast_db' => {
      -dbname => $self->o('rnaseq_blast_db_name'),
      -host   => $self->o('rnaseq_blast_db_server'),
      -port   => $self->o('rnaseq_blast_db_port'),
      -user   => $self->o('user_r'),
      -driver => $self->o('hive_driver'),
    },

    'layering_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_layer_'.$self->o('release_number'),
      -host   => $self->o('layering_db_server'),
      -port   => $self->o('layering_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'utr_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_utr_'.$self->o('release_number'),
      -host   => $self->o('utr_db_server'),
      -port   => $self->o('utr_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'genebuilder_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_gbuild_'.$self->o('release_number'),
      -host   => $self->o('genebuilder_db_server'),
      -port   => $self->o('genebuilder_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'pseudogene_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_pseudo_'.$self->o('release_number'),
      -host   => $self->o('pseudogene_db_server'),
      -port   => $self->o('pseudogene_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'ncrna_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_ncrna_'.$self->o('release_number'),
      -host   => $self->o('ncrna_db_server'),
      -port   => $self->o('ncrna_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'final_geneset_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_final_'.$self->o('release_number'),
      -host   => $self->o('final_geneset_db_server'),
      -port   => $self->o('final_geneset_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'refseq_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_refseq_'.$self->o('release_number'),
      -host   => $self->o('refseq_db_server'),
      -port   => $self->o('refseq_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'killlist_db' => {
      -dbname => $self->o('killlist_db_name'),
      -host   => $self->o('killlist_db_server'),
      -port   => $self->o('killlist_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -driver => $self->o('hive_driver'),
    },

    'production_db' => {
      -host   => $self->o('staging_1_db_server'),
      -port   => $self->o('staging_1_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -dbname => 'ensembl_production_'.$self->o('ensembl_release'),
      -driver => $self->o('hive_driver'),
    },

    'taxonomy_db' => {
      -host   => $self->o('staging_1_db_server'),
      -port   => $self->o('staging_1_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -dbname => 'ncbi_taxonomy',
      -driver => $self->o('hive_driver'),
    },

    'otherfeatures_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('production_name').'_otherfeatures_'.$self->o('release_number'),
      -host   => $self->o('otherfeatures_db_server'),
      -port   => $self->o('otherfeatures_db_port'),
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
      @{$self->SUPER::pipeline_create_commands},

      $self->hive_data_table('protein', $self->o('uniprot_table_name')),

      $self->hive_data_table('refseq', $self->o('cdna_table_name')),

      $self->db_cmd('CREATE TABLE '.$self->o('realign_table_name').' ('.
                    'accession varchar(50) NOT NULL,'.
                    'seq text NOT NULL,'.
                    'PRIMARY KEY (accession))'),
    ];
}


sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
    skip_projection => $self->o('skip_projection'),
    skip_rnaseq => $self->o('skip_rnaseq'),
    skip_ncrna => $self->o('skip_ncrna'),
    load_toplevel_only => $self->o('load_toplevel_only'),
  }
}

## See diagram for pipeline structure
sub pipeline_analyses {
    my ($self) = @_;

    my %genblast_params = (
      wu    => '-P wublast -gff -e #blast_eval# -c #blast_cov#',
      ncbi  => '-P blast -gff -e #blast_eval# -c #blast_cov# -W 3 -rl 5000',
      wu_genome    => '-P wublast -gff -e #blast_eval# -c #blast_cov#',
      ncbi_genome  => '-P blast -gff -e #blast_eval# -c #blast_cov# -W 3 -rl 5000 -softmask',
      wu_projection    => '-P wublast -gff -e #blast_eval# -c #blast_cov# -n 100 -rl 5000 -x 5 ',
      ncbi_projection  => '-P blast -gff -e #blast_eval# -c #blast_cov# -W 3 -rl 5000',
      );
    my %commandline_params = (
      'ncbi' => '-num_threads 3 -window_size 40',
      'wu' => '-cpus 3 -hitdist 40',
      'legacy_ncbi' => '-a 3 -A 40',
      );

    return [


###############################################################################
#
# ASSEMBLY LOADING ANALYSES
#
###############################################################################

      {
        # Creates a reference db for each species
        -logic_name => 'create_core_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         'target_db'        => $self->o('reference_db'),
                         'enscode_root_dir' => $self->o('enscode_root_dir'),
                         'create_type'      => 'core_only',
                       },
        -rc_name    => 'default',
        -flow_into  => {
                         1 => ['populate_production_tables'],
                       },
        -input_ids  => [
          {
            assembly_name => $self->o('assembly_name'),
            assembly_accession => $self->o('assembly_accession'),
            assembly_refseq_accession => $self->o('assembly_refseq_accession'),
          },
        ],

      },

      {
        # Load production tables into each reference
        -logic_name => 'populate_production_tables',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HivePopulateProductionTables',
        -parameters => {
                         'target_db'        => $self->o('reference_db'),
                         'output_path'      => $self->o('output_path'),
                         'enscode_root_dir' => $self->o('enscode_root_dir'),
                         'production_db'    => $self->o('production_db'),
                       },
        -rc_name    => 'default',
        -flow_into  => {
                         1 => WHEN ('#load_toplevel_only#' => ['process_assembly_info'],
                              ELSE ['download_assembly_info']),
                       },
      },


      {
        # Download the files and dir structure from the NCBI ftp site. Uses the link to a species in the ftp_link_file
        -logic_name => 'process_assembly_info',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProcessAssemblyReport',
        -parameters => {
          full_ftp_path => $self->o('assembly_ftp_path'),
          output_path   => $self->o('output_path'),
          target_db     => $self->o('reference_db'),
        },
        -rc_name    => 'default',
        -flow_into  => {
          '2->A' => ['download_assembly_fasta'],
          'A->1' => ['load_toplevel_sequences'],
        },
      },


      {
        # Download the files and dir structure from the NCBI ftp site. Uses the link to a species in the ftp_link_file
        -logic_name => 'download_assembly_fasta',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadData',
        -parameters => {
          output_dir => $self->o('output_path'),
          download_method => 'ftp',
        },
        -rc_name    => 'default',
      },


      {
        # Load the contigs into each reference db
        -logic_name => 'load_toplevel_sequences',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadGenomeSequences',
        -parameters => {
          target_db   => $self->o('reference_db'),
          output_path => $self->o('output_path'),
        },
        -rc_name    => '4GB',
        -flow_into  => {
                         1 => ['load_meta_info'],
                       },
      },


      {
        # Load some meta info and seq_region_synonyms
        -logic_name => 'load_meta_info',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('reference_db'),
          sql => [
            'INSERT INTO meta (species_id, meta_key, meta_value) VALUES '.
              '(1, "assembly.coverage_depth", "high"),'.
              '(1, "genebuild.id", '.$self->o('genebuilder_id').'),'.
              '(1, "genebuild.method", "full_genebuild"),'.
              '(1, "genebuild.projection_source_db", "'.$self->o('projection_source_db_name').'"),'.
              '(1, "provider.name", "'.$self->o('provider_name').'"),'.
              '(1, "provider.url", "'.$self->o('provider_url').'"),'.
              '(1, "species.production_name", "'.$self->o('production_name').'"),'.
              '(1, "repeat.analysis", "repeatmask_repbase_'.$self->o('repbase_logic_name').'"),'.
              '(1, "repeat.analysis", "dust"),'.
              '(1, "repeat.analysis", "trf")',
          ],
        },
        -rc_name    => 'default',
        -flow_into  => {
          1 => ['load_taxonomy_info'],
                       },
      },

      {
#       Download the files and dir structure from the NCBI ftp site. Uses the link to a species in the ftp_link_file
        -logic_name => 'download_assembly_info',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveDownloadNCBIFtpFiles',
        -parameters => {
          'full_ftp_path'             => $self->o('assembly_ftp_path'),
          'output_path'               => $self->o('output_path'),
          'primary_assembly_dir_name' => $self->o('primary_assembly_dir_name'),
        },
        -rc_name    => 'default',
        -flow_into  => {
          1 => ['find_contig_accessions'],
        },
      },

      {
# Get the prefixes for all contigs from the AGP files
        -logic_name => 'find_contig_accessions',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveFindContigAccessions',
        -parameters => {
          'contigs_source'            => $self->o('contigs_source'),
          'wgs_id'                    => $self->o('wgs_id'),
          'output_path'               => $self->o('output_path'),
          'primary_assembly_dir_name' => $self->o('primary_assembly_dir_name'),
        },
        -rc_name => 'default',
        -flow_into => {
          1 => ['download_contigs'],
        },
      },

      {
# Download contig from NCBI
        -logic_name => 'download_contigs',
        -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveDownloadContigs',
        -parameters => {
          'contigs_source' => $self->o('contigs_source'),
          'wgs_id' => $self->o('wgs_id'),
          'output_path' => $self->o('output_path'),
          'primary_assembly_dir_name' => $self->o('primary_assembly_dir_name'),
        },
        -rc_name => 'default',
        -flow_into => {
          1 => ['load_contigs'],
        },
      },


      {
# Load the contigs into each reference db
        -logic_name => 'load_contigs',
        -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadSeqRegions',
        -parameters => {
          'coord_system_version' => $self->o('assembly_name'),
          'target_db' => $self->o('reference_db'),
          'output_path' => $self->o('output_path'),
          'enscode_root_dir' => $self->o('enscode_root_dir'),
          'primary_assembly_dir_name' => $self->o('primary_assembly_dir_name'),
        },
        -rc_name => '4GB',
        -flow_into => {
          1 => ['load_assembly_info'],
        },
      },


      {
# Load the AGP files
        -logic_name => 'load_assembly_info',
        -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadAssembly',
        -parameters => {
          'target_db' => $self->o('reference_db'),
          'output_path' => $self->o('output_path'),
          'primary_assembly_dir_name' => $self->o('primary_assembly_dir_name'),
        },
        -rc_name => 'default',
        -flow_into => {
          1 => ['set_toplevel'],
        },
      },


      {
# Set the toplevel
        -logic_name => 'set_toplevel',
        -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveSetAndCheckToplevel',
        -parameters => {
          'target_db' => $self->o('reference_db'),
          'output_path' => $self->o('output_path'),
          'enscode_root_dir' => $self->o('enscode_root_dir'),
          'primary_assembly_dir_name' => $self->o('primary_assembly_dir_name'),
        },
        -rc_name => '1.5GB',
        -flow_into => {
          1 => ['load_meta_info_full'],
        },
      },


      {
# Load some meta info and seq_region_synonyms
        -logic_name => 'load_meta_info_full',
        -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveSetMetaAndSeqRegionSynonym',
        -parameters => {
          'target_db' => $self->o('reference_db'),
          'output_path' => $self->o('output_path'),
          'enscode_root_dir' => $self->o('enscode_root_dir'),
          'primary_assembly_dir_name' => $self->o('primary_assembly_dir_name'),
          'meta_key_list' => {
            'assembly.accession' => $self->o('assembly_accession'),
            'assembly.coverage_depth' => 'high',
            'assembly.default' => $self->default_options->{'assembly_name'},
            'assembly.name' => $self->default_options->{'assembly_name'},
            'assembly.web_accession_source' => 'NCBI',
            'assembly.web_accession_type' => 'GenBank Assembly ID',
            'genebuild.id' => $self->o('genebuilder_id'),
            'genebuild.method' => 'full_genebuild',
            'genebuild.projection_source_db' => $self->default_options->{'projection_source_db_name'},
            'provider.name' => $self->o('provider_name'),
            'provider.url' => $self->o('provider_url'),
            'repeat.analysis' => ['repeatmask_repbase_'.$self->o('repbase_logic_name'), 'dust', 'trf'],
            'species.production_name' => $self->o('production_name'),
            'species.taxonomy_id' => $self->default_options->{'taxon_id'},
          }
        },
        -rc_name => 'default',
        -flow_into => {
          1 => ['load_taxonomy_info'],
        },
      },

      {
        -logic_name => 'load_taxonomy_info',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadTaxonomyInfo',
        -parameters => {
                         'target_db'        => $self->o('reference_db'),
                         'enscode_root_dir' => $self->o('enscode_root_dir'),
                         'production_db'    => $self->o('production_db'),
                         'taxonomy_db'      => $self->o('taxonomy_db'),
                       },
        -rc_name    => 'default',

        -flow_into  => {
                          1 => ['load_mitochondrion', 'fan_refseq_import'],
                       },
      },

###############################################################################
#
# REFSEQ ANNOTATION
#
###############################################################################
      {
        -logic_name => 'fan_refseq_import',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'if [ -n "#assembly_refseq_accession#" ]; then exit 0; else exit 42;fi',
                         return_codes_2_branches => {'42' => 2},
                       },
        -rc_name => 'default',
        -flow_into  => {
                          1 => ['load_refseq_synonyms'],
                       },
      },

      {
        -logic_name => 'load_refseq_synonyms',
        -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadRefSeqSynonyms',
        -parameters => {
                         'target_db'        => $self->o('reference_db'),
                       },
        -rc_name => 'default',
        -flow_into  => {
                          1 => ['download_refseq_gff'],
                       },
      },


      {
        -logic_name => 'download_refseq_gff',
        -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadData',
        -parameters => {
          output_dir => catdir($self->o('output_path'), 'refseq_import'),
          url => $self->o('refseq_import_ftp_path'),
          download_method => 'ftp',
        },
        -rc_name => 'default',
        -flow_into  => {
                         1 => ['create_refseq_db'],
                       },
      },


      {
        -logic_name => 'create_refseq_db',
        -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('refseq_db'),
                         create_type => 'clone',
                       },
       -rc_name => 'default',
       -flow_into  => {
                          1 => ['load_refseq_gff'],
                      },
      },

      {
        -logic_name => 'load_refseq_gff',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('refseq_import_script_path').
                                      ' -dnahost '.$self->o('dna_db','-host').
                                      ' -dnadbname '.$self->o('dna_db','-dbname').
                                      ' -dnaport '.$self->o('dna_db','-port').
                                      ' -user '.$self->o('user').
                                      ' -pass '.$self->o('password').
                                      ' -host '.$self->o('refseq_db','-host').
                                      ' -port '.$self->o('refseq_db','-port').
                                      ' -dbname '.$self->o('refseq_db','-dbname').
                                      ' -infile '.$self->o('output_path').'/refseq_import/#assembly_refseq_accession#_#assembly_name#_genomic.gff',
                       },
        -rc_name => 'refseq_import',
      },


      {
        -logic_name => 'load_mitochondrion',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadMitochondrion',
        -parameters => {
          'target_db'        => $self->o('reference_db'),
          'output_path'      => $self->o('output_path'),
          'enscode_root_dir' => $self->o('enscode_root_dir'),
          'species_name'     => $self->o('species_name'),
          'mt_accession'     => $self->o('mt_accession'),
        },
        -rc_name    => 'default',
        -flow_into  => {
          '1->A' => ['create_10mb_slice_ids'],
          'A->1' => ['genome_prep_sanity_checks'],
        },
      },

###############################################################################
#
# REPEATMASKER ANALYSES
#
###############################################################################
      {
        # Create 10mb toplevel slices, these will be split further for repeatmasker
        -logic_name => 'create_10mb_slice_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db        => $self->o('dna_db'),
                         coord_system_name => 'toplevel',
                         iid_type => 'slice',
                         slice_size => 10000000,
                         include_non_reference => 0,
                         top_level => 1,
                         min_slice_length => $self->o('min_toplevel_slice_length'),
                         batch_slice_ids => 1,
                         batch_target_size => 10000000,
                       },
        -flow_into => {
                         '2'    => ['semaphore_10mb_slices'],
                      },
      },


      {
         # Wait for repeatmasker to complete all the sub slices for a 10mb slice before passing to dust
        -logic_name => 'semaphore_10mb_slices',
         -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
         -parameters => {},
         -flow_into => {
                        '1->A' => ['create_repeatmasker_slices'],
                        'A->1' => ['run_dust'],
                       },
      },


      {
        -logic_name => 'create_repeatmasker_slices',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('dna_db'),
                         iid_type => 'rebatch_and_resize_slices',
                         slice_size => 1000000,
                         batch_target_size => 500000,
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '2' => ['run_repeatmasker','fan_repeatmodeler'],
                      },
      },


      {
        -logic_name => 'run_repeatmasker',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveRepeatMasker',
        -parameters => {
                         timer_batch => '5h',
                         target_db => $self->o('reference_db'),
                         logic_name => 'repeatmask_repbase_'.$self->o('repbase_logic_name'),
                         module => 'HiveRepeatMasker',
                         repeatmasker_path => $self->o('repeatmasker_path'),
                         commandline_params => '-nolow -species "'.$self->o('repbase_library').'" -engine "'.$self->o('repeatmasker_engine').'"',
                       },
        -rc_name    => 'repeatmasker',
        -flow_into => {
                        '-1' => ['rebatch_repeatmasker'],
                        '-2' => ['rebatch_repeatmasker'],
                      },
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
      },


      {
        -logic_name => 'rebatch_repeatmasker',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('dna_db'),
                         iid_type => 'rebatch_and_resize_slices',
                         slice_size => 100000,
                         batch_target_size => 10000,
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '2' => ['run_repeatmasker_small_batch'],
                      },
        -can_be_empty  => 1,
      },


      {
        -logic_name => 'run_repeatmasker_small_batch',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveRepeatMasker',
        -parameters => {
                         timer_batch => '2h',
                         target_db => $self->o('reference_db'),
                         logic_name => 'repeatmask_repbase_'.$self->o('repbase_logic_name'),
                         module => 'HiveRepeatMasker',
                         repeatmasker_path => $self->o('repeatmasker_path'),
                         commandline_params => '-nolow -species "'.$self->o('repbase_library').'" -engine "'.$self->o('repeatmasker_engine').'"',
                       },
        -rc_name    => 'repeatmasker_rebatch',
        -flow_into => {
                         -1 => ['failed_repeatmasker_batches'],
                         -2 => ['failed_repeatmasker_batches'],
                      },
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -can_be_empty  => 1,
      },


      {
        -logic_name => 'failed_repeatmasker_batches',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name          => 'default',
        -can_be_empty  => 1,
      },


      {
        -logic_name => 'fan_repeatmodeler',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'if [ -n "'.$self->o('repeatmodeler_library').'" ]; then exit 0; else exit 42;fi',
          return_codes_2_branches => {'42' => 2},
        },
        -rc_name    => 'default',
        -flow_into  => { '1' => ['run_repeatmsker_repeatmodeler'] },
      },


      {
        -logic_name => 'run_repeatmsker_repeatmodeler',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveRepeatMasker',
        -parameters => {
                         timer_batch => '3h',
                         target_db => $self->o('reference_db'),
                         logic_name => 'repeatmask_repeatmodeler',
                         module => 'HiveRepeatMasker',
                         repeatmasker_path => $self->o('repeatmasker_path'),
                         commandline_params => '-nolow -lib "'.$self->o('repeatmodeler_library').'" -engine "'.$self->o('repeatmasker_engine').'"',
                       },
        -rc_name    => 'repeatmasker',
        -flow_into => {
                        '-1' => ['rebatch_repeatmasker_repeatmodeler'],
                        '-2' => ['rebatch_repeatmasker_repeatmodeler'],
                      },
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
      },


      {
        -logic_name => 'rebatch_repeatmasker_repeatmodeler',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('dna_db'),
                         iid_type => 'rebatch_and_resize_slices',
                         slice_size => 100000,
                         batch_target_size => 10000,
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '2' => ['run_repeatmasker_repeatmodeler_small_batch'],
                      },
        -can_be_empty  => 1,
      },


      {
        -logic_name => 'run_repeatmasker_repeatmodeler_small_batch',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveRepeatMasker',
        -parameters => {
                         timer_batch => '1h',
                         target_db => $self->o('reference_db'),
                         logic_name => 'repeatmask_repeatmodeler',
                         module => 'HiveRepeatMasker',
                         repeatmasker_path => $self->o('repeatmasker_path'),
                         commandline_params => '-nolow -lib "'.$self->o('repeatmodeler_library').'" -engine "'.$self->o('repeatmasker_engine').'"',
                       },
        -rc_name    => 'repeatmasker_rebatch',
        -flow_into => {
                         -1 => ['failed_repeatmasker_repeatmodeler_batches'],
                         -2 => ['failed_repeatmasker_repeatmodeler_batches'],
                      },
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -can_be_empty  => 1,
      },


      {
        -logic_name => 'failed_repeatmasker_repeatmodeler_batches',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name          => 'default',
        -can_be_empty  => 1,
      },



      {
        # Set the toplevel
        -logic_name => 'dump_softmasked_toplevel',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDumpGenome',
        -parameters => {
                         'coord_system_name'    => 'toplevel',
                         'target_db'            => $self->o('reference_db'),
                         'output_path'          => $self->o('output_path')."/genome_dumps/",
                         'enscode_root_dir'     => $self->o('enscode_root_dir'),
                         'species_name'         => $self->o('species_name'),
                         'repeat_logic_names'   => $self->o('repeat_logic_names'),
                       },
        -input_ids => [{}],
        -wait_for => ['run_dust'],
        -flow_into => {
          1 => ['format_softmasked_toplevel'],
        },
        -rc_name    => 'default_himem',
      },


      {
        # This should probably be a proper module
        -logic_name => 'format_softmasked_toplevel',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         'cmd'    => 'if [ "'.$self->o('blast_type').'" = "ncbi" ]; then convert2blastmask -in '.catfile($self->o('output_path'), 'genome_dumps', $self->o('species_name')).'_softmasked_toplevel.fa -parse_seqids -masking_algorithm repeatmasker -masking_options "repeatmasker, default" -outfmt maskinfo_asn1_bin -out '.catfile($self->o('output_path'), 'genome_dumps', $self->o('species_name')).'_softmasked_toplevel.fa.asnb;makeblastdb -in '.catfile($self->o('output_path'), 'genome_dumps', $self->o('species_name')).'_softmasked_toplevel.fa -dbtype nucl -parse_seqids -mask_data '.catfile($self->o('output_path'), 'genome_dumps', $self->o('species_name')).'_softmasked_toplevel.fa.asnb -title "'.$self->o('species_name').'"; else xdformat -n '.catfile($self->o('output_path'), 'genome_dumps', $self->o('species_name')).'_softmasked_toplevel.fa;fi',
                       },
        -rc_name    => 'default_himem',
      },

###############################################################################
#
# SIMPLE FEATURE AND OTHER REPEAT ANALYSES
#
###############################################################################

      {
        # Run dust
        -logic_name => 'run_dust',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveDust',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         logic_name => 'dust',
                         module => 'HiveDust',
                         dust_path => $self->o('dust_path'),
                       },
        -rc_name    => 'simple_features',
        -flow_into => {
                         1 => ['run_trf'],
                         -1 => ['run_trf'],
		         -2 => ['run_trf'],
                      },
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -batch_size => 20,
      },


      {
        # Run TRF
        -logic_name => 'run_trf',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveTRF',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         logic_name => 'trf',
                         module => 'HiveTRF',
                         trf_path => $self->o('trf_path'),
                       },
        -rc_name    => 'simple_features',
        -flow_into => {
                         1 => ['run_eponine'],
                        -1 => ['run_eponine'],
                        -2 => ['run_eponine'],
                      },
       -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
       -batch_size => 20,
      },


      {
        # Run eponine
        -logic_name => 'run_eponine',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveEponine',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         logic_name => 'eponine',
                         module => 'HiveEponine',
                         eponine_path => $self->o('eponine_java_path'),
                         commandline_params => '-epojar => '.$self->o('eponine_jar_path').', -threshold => 0.999',
                       },
        -rc_name    => 'simple_features',
        -flow_into => {
                         1 => ['run_cpg'],
                         -1 => ['run_cpg'],
                         -2 => ['run_cpg'],
                      },
       -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
       -batch_size => 20,
      },


      {
        # Run CPG
        -logic_name => 'run_cpg',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveCPG',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         logic_name => 'cpg',
                         module => 'HiveCPG',
                         cpg_path => $self->o('cpg_path'),
                       },
        -rc_name    => 'simple_features',
        -flow_into => {
                        1 => ['run_trnascan'],
                        -1 => ['run_trnascan'],
                        -2 => ['run_trnascan'],
                      },
       -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
       -batch_size => 20,
      },


      {
        # Run tRNAscan
        -logic_name => 'run_trnascan',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveTRNAScan',
        -parameters => {
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('reference_db'),
                         logic_name => 'trnascan',
                         module => 'HiveTRNAScan',
                         trnascan_path => $self->o('trnascan_path'),
                       },
        -rc_name    => 'simple_features',
        -flow_into => {
                         1 => ['create_genscan_slices'],
                         -1 => ['create_genscan_slices'],
                         -2 => ['create_genscan_slices'],

                      },
         -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
         -batch_size => 20,
      },

###############################################################################
#
# GENSCAN ANALYSIS
#
##############################################################################
      {
         # Genscan has issues with large slices, so instead rebatch the 10mb slices
         # into 1mb slices with a target batch size of 10mb
        -logic_name => 'create_genscan_slices',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('dna_db'),
                         iid_type => 'rebatch_and_resize_slices',
                         slice_size => 1000000,
                         batch_target_size => 10000000,
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '2' => ['run_genscan'],
                      },
      },


      {
        # Run genscan, uses 1mb slices from repeatmasker. Flows into create_prediction_transcript_ids which
        # then takes these 1mb slices and converts them into individual prediction transcript input ids based
        # on the dbID of each feature generate by this analysis
        -logic_name => 'run_genscan',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveGenscan',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         logic_name => 'genscan',
                         module => 'HiveGenscan',
                         genscan_path => $self->o('genscan_path'),
                         genscan_matrix_path => $self->o('genscan_matrix_path'),
                         repeat_masking_logic_names => ['repeatmask_repbase_'.$self->o('repbase_logic_name')],
                       },
        -rc_name    => 'genscan',
        -flow_into => {
                        # No need to semaphore the jobs with issues as the blast analyses work off prediction transcript
                        # ids from slices that genscan succeeds on. So passing small slices in and ignore failed slices is fine
                        1 => ['create_prediction_transcript_ids'],
                        -1 => ['decrease_genscan_slice_size'],
                        -2 => ['decrease_genscan_slice_size'],
                        -3 => ['decrease_genscan_slice_size'],
                      },
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -batch_size => 20,
      },


      {
        -logic_name => 'decrease_genscan_slice_size',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db        => $self->o('dna_db'),
                         iid_type => 'rebatch_and_resize_slices',
                         slice_size => 100000,
                         batch_target_size => 100000,
                       },
        -flow_into => {
                        '2' => ['run_genscan_short_slice'],
                      },
        -rc_name    => 'default',
        -can_be_empty  => 1,
      },


      {
        -logic_name => 'run_genscan_short_slice',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveGenscan',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         logic_name => 'genscan',
                         module => 'HiveGenscan',
                         genscan_path => $self->o('genscan_path'),
                         genscan_matrix_path => $self->o('genscan_matrix_path'),
                         repeat_masking_logic_names => ['repeatmask_repbase_'.$self->o('repbase_logic_name')],
                       },
        -rc_name    => 'genscan',
        -flow_into => {
                        1 => ['create_prediction_transcript_ids'],
                        -1 => ['failed_genscan_slices'],
                        -2 => ['failed_genscan_slices'],
                        -3 => ['failed_genscan_slices'],
                      },
        -rc_name    => 'genscan_short',
        -can_be_empty  => 1,
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
      },


      {
        -logic_name => 'failed_genscan_slices',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name          => 'default',
        -can_be_empty  => 1,
      },


      {
        # Create input ids for prediction transcripts. Takes a slice as an input id and converts it
        # to a set of input ids that are individual dbIDs for the prediction transcripts. This avoids empty slices
        # being submitted as jobs and also means one feature corresponds to one job. Each species flows into this
        # independantly with 1mb slices. Also has the advantage that downstream analyses can start working as soon
        # as a single slice is finished
        -logic_name => 'create_prediction_transcript_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('dna_db'),
                         feature_type => 'prediction_transcript',
                         iid_type => 'slice_to_feature_ids',
                         prediction_transcript_logic_names => ['genscan'],
                         batch_size => 50,
                       },
        -flow_into => {
                        2 => ['run_uniprot_blast'],
                      },
        -rc_name    => 'default',
      },

##############################################################################
#
# BLAST ANALYSES
#
##############################################################################

      {
        # BLAST individual prediction transcripts against uniprot. The config settings are held lower in this
        # file in the master_config_settings sub
        -logic_name => 'run_uniprot_blast',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlastGenscan',
        -parameters => {
                         timer => '30m',
                         sequence_type => 'peptide',
                         prediction_transcript_db => $self->o('dna_db'),
                         target_db => $self->o('reference_db'),
                         repeat_masking_logic_names => ['repeatmask_repbase_'.$self->o('repbase_logic_name')], # not sure if this is used
                         blast_db_path => $self->o('uniprot_blast_db_path'),
                         blast_exe_path => $self->o('uniprot_blast_exe_path'),
                         commandline_params => $commandline_params{$self->o('blast_type')},
                         iid_type => 'feature_id',
                         logic_name => 'uniprot',
                         module => 'HiveBlastGenscan',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic',
                                                 'BlastGenscanPep',
                                                 {BLAST_PARAMS => {type => $self->o('blast_type')}})},
                      },
        -flow_into => {
                        1 => ['run_vertrna_blast'],
                        -1 => ['split_blast_jobs'],
                        -2 => ['split_blast_jobs'],
                      },
        -rc_name    => 'blast',
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -batch_size => 20,
      },


      {
        # Only do the split on the uniprot ones as they're the only ones that seem to have any issues
        -logic_name => 'split_blast_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         iid_type => 'rechunk',
                         batch_size => 1,
                       },
        -rc_name      => 'default',
        -can_be_empty  => 1,
        -flow_into => {
                        2 => ['run_uniprot_blast_retry'],
                      },
      },


      {
        # BLAST individual prediction transcripts against uniprot. The config settings are held lower in this
        # file in the master_config_settings sub
        -logic_name => 'run_uniprot_blast_retry',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlastGenscan',
        -parameters => {
                         timer => '45m',
                         sequence_type => 'peptide',
                         prediction_transcript_db => $self->o('dna_db'),
                         target_db => $self->o('reference_db'),
                         repeat_masking_logic_names => ['repeatmask_repbase_'.$self->o('repbase_logic_name')], # not sure if this is used
                         blast_db_path => $self->o('uniprot_blast_db_path'),
                         blast_exe_path => $self->o('uniprot_blast_exe_path'),
                         commandline_params => $commandline_params{$self->o('blast_type')},
                         iid_type => 'feature_id',
                         logic_name => 'uniprot',
                         module => 'HiveBlastGenscan',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic',
                                                 'BlastGenscanPep',
                                                 {BLAST_PARAMS => {type => $self->o('blast_type')}})},
                      },
        -can_be_empty  => 1,
        -flow_into => {
                        1 => ['run_vertrna_blast'],
                        -1 => ['run_vertrna_blast'],
                        -2 => ['run_vertrna_blast'],
                      },
        -rc_name    => 'blast_retry',
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -batch_size => 20,
      },


      {
        # BLAST individual prediction transcripts against vertRNA. The config settings are held lower in this
        # file in the master_config_settings sub
        -logic_name => 'run_vertrna_blast',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlastGenscan',
        -parameters => {
                         timer => '10m',
                         sequence_type => 'dna',
                         prediction_transcript_db => $self->o('dna_db'),
                         target_db => $self->o('reference_db'),
                         repeat_masking_logic_names => ['repeatmask_repbase_'.$self->o('repbase_logic_name')], # not sure if this is used
                         blast_db_path => $self->o('vertrna_blast_db_path'),
                         blast_exe_path => $self->o('vertrna_blast_exe_path'),
                         commandline_params => $commandline_params{$self->o('blast_type')},
                         iid_type => 'feature_id',
                         logic_name => 'vertrna',
                         module => 'HiveBlastGenscan',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic',
                                                 'BlastGenscanVertRNA',
                                                 {BLAST_PARAMS => {type => $self->o('blast_type')}})},
                      },
        -flow_into => {
                        1 => ['run_unigene_blast'],
                        -1 => ['run_unigene_blast'],
                        -2 => ['run_unigene_blast'],
                      },
        -rc_name    => 'blast',
       -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
       -batch_size => 20,
      },


      {
        # BLAST individual prediction transcripts against unigene. The config settings are held lower in this
        # file in the master_config_settings sub
        -logic_name => 'run_unigene_blast',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlastGenscan',
        -parameters => {
                         timer => '10m',
                         sequence_type => 'dna',
                         prediction_transcript_db => $self->o('dna_db'),
                         target_db => $self->o('reference_db'),
                         repeat_masking_logic_names => ['repeatmask_repbase_'.$self->o('repbase_logic_name')], # not sure if this is used
                         blast_db_path => $self->o('unigene_blast_db_path'),
                         blast_exe_path => $self->o('unigene_blast_exe_path'),
                         commandline_params => $commandline_params{$self->o('blast_type')},
                         iid_type => 'feature_id',
                         logic_name => 'unigene',
                         module => 'HiveBlastGenscan',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic',
                                                 'BlastGenscanUnigene',
                                                 {BLAST_PARAMS => {type => $self->o('blast_type')}})},
                       },
        -flow_into => {
                        -1 => ['failed_blast_job'],
                        -2 => ['failed_blast_job'],
                      },
        -rc_name    => 'blast',
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -batch_size => 20,
      },


      {
        -logic_name => 'failed_blast_job',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name          => 'default',
        -can_be_empty  => 1,
        -failed_job_tolerance => 100,
      },

      {
        -logic_name => 'genome_prep_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         target_db => $self->o('dna_db'),
                         sanity_check_type => 'genome_preparation_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'genome_preparation_checks')->{$self->default_options->{'uniprot_set'}},
                       },

        -flow_into =>  {
                         1 => ['backup_core_db'],
                       },
        -rc_name    => '16GB',
     },

     {
        # Creates a reference db for each species
        -logic_name => 'backup_core_db',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::DatabaseDumper',
        -parameters => {
                         src_db_conn => $self->o('dna_db'),
                         output_file => catfile($self->o('output_path'), 'core_bak.sql.gz'),
                       },
        -rc_name    => 'default',
        -flow_into => { 1 => ['assembly_loading_report'] },
      },

      {
        -logic_name => 'assembly_loading_report',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('loading_report_script').
                                ' -user '.$self->o('user_r').
                                ' -host '.$self->o('dna_db','-host').
                                ' -port '.$self->o('dna_db','-port').
                                ' -dbname '.$self->o('dna_db','-dbname').
                                ' -report_type assembly_loading'.
                                ' > '.$self->o('output_path').'/loading_report.txt'
                       },
         -rc_name => 'default',
         -flow_into => { 1 => ['set_repeat_types','email_loading_report'] },
      },

      {
        -logic_name => 'email_loading_report',
        -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::TextfileByEmail',
        -parameters => {
                         email => $self->o('email_address'),
                         subject => 'AUTOMATED REPORT: assembly loading and feature annotation for '.$self->o('dna_db','-dbname').' completed',
                         text => 'Assembly loading and feature annotation have completed for '.$self->o('dna_db','-dbname').". Basic stats can be found below",
                         file => $self->o('output_path').'/loading_report.txt',
                       },
        -rc_name => 'default',
        -flow_into => {
          1 =>['repeat_masker_coverage'],
        },
      },
      {
        -logic_name => 'repeat_masker_coverage',
        -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRepeatCoverage',
        -parameters => {
          source_db => $self->o('reference_db'),
          repeat_logic_names => $self->o('repeat_logic_names'),
          coord_system_version => '#assembly_name#',
          email => $self->o('email_address'),
        },
        -rc_name => 'default',
      },

     {
        -logic_name => 'set_repeat_types',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('repeat_types_script').
                                ' -user '.$self->o('reference_db', '-user').
                                ' -pass '.$self->o('reference_db', '-pass').
                                ' -host '.$self->o('reference_db','-host').
                                ' -port '.$self->o('reference_db','-port').
                                ' -dbpattern '.$self->o('reference_db','-dbname')
                       },
         -rc_name => 'default',
         -flow_into => { 1 => ['create_genblast_output_db'] },
      },


######################################################################################
#
# Protein models (genblast and genewise)
#
######################################################################################

      {
        -logic_name => 'create_genblast_output_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('genblast_db'),
                         create_type => 'clone',
                       },
        -rc_name    => 'default',
        -wait_for => ['format_softmasked_toplevel'],
        -flow_into => { 1 => ['create_genblast_select_output_db'] },
      },


      {
        -logic_name => 'create_genblast_select_output_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('genblast_select_db'),
                         create_type => 'clone',
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '1->A' => ['download_uniprot_files'],
                        'A->1' => ['classify_genblast_models'],
                      },
      },


      {

        -logic_name => 'download_uniprot_files',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadUniProtFiles',
        -parameters => {
                         multi_query_download => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::UniProtCladeDownloadStatic', $self->default_options->{'uniprot_set'}),
                         taxon_id => $self->o('taxon_id'),
                         output_path => $self->o('homology_models_path'),
                       },
        -rc_name          => 'default',
        -flow_into => {
                        '2->A' => ['process_uniprot_files'],
                        'A->1' => ['generate_genblast_jobs'],
                      },
      },

      {
        -logic_name => 'process_uniprot_files',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProcessUniProtFiles',
        -parameters => {
                         killlist_type => 'protein',
                         killlist_db => $self->o('killlist_db'),
                         sequence_table_name => $self->o('uniprot_table_name'),
                      },
        -rc_name => 'default',
      },



      {
        -logic_name => 'generate_genblast_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         iid_type => 'sequence_accession',
                         batch_size => $self->o('uniprot_genblast_batch_size'),
                         sequence_table_name => $self->o('uniprot_table_name'),
                       },
        -rc_name      => 'default_himem',
        -flow_into => {
                        2 => ['genblast'],
                        1 => ['create_seleno_homology_jobs'],
                      },
      },

      {
        -logic_name => 'create_seleno_homology_jobs',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
        -parameters => {
          inputquery => 'SELECT accession FROM '.$self->o('uniprot_table_name').' WHERE source_db = "seleno"',
          column_names => ['iid'],
        },
        -rc_name          => 'default',
        -flow_into => {
          2 => ['process_homology_selenocysteine'],
        },
      },

      {
        -logic_name => 'process_homology_selenocysteine',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSelenocysteineFinder',
        -parameters => {
          target_db => $self->o('genblast_db'),
          dna_db => $self->o('dna_db'),
          genome => $self->o('genome_file'),
          exonerate => $self->o('exonerate_path'),
          genewise => $self->o('genewise_path'),
          iid_type => 'db_seq',
          sequence_table_name => $self->o('uniprot_table_name'),
          biotype => 'seleno_other',
          missmatch_allowed => 10,
        },
        -rc_name          => 'default',
      },



      {
        -logic_name => 'genblast',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
        -parameters => {
                         iid_type => 'db_seq',
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('genblast_db'),
                         logic_name => 'genblast',
                         module => 'HiveGenblast',
                         genblast_path => $self->o('genblast_path'),
                         genblast_db_path => $self->o('genome_file'),
                         commandline_params => $genblast_params{$self->o('blast_type').'_genome'},
                         sequence_table_name => $self->o('uniprot_table_name'),
                         max_rank => $self->o('genblast_max_rank'),
                         genblast_pid => $self->o('genblast_pid'),
                         timer => '2h',
                         blast_eval => $self->o('genblast_eval'),
                         blast_cov  => $self->o('genblast_cov'),
                         flag_small_introns => $self->o('genblast_flag_small_introns'),
                         flag_subpar_models => $self->o('genblast_flag_subpar_models'),
                       },
        -rc_name    => 'genblast',
        -flow_into => {
                        -1 => ['split_genblast_jobs'],
                        -2 => ['split_genblast_jobs'],
                        -3 => ['split_genblast_jobs'],
                      },
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
      },

      {
        -logic_name => 'split_genblast_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         iid_type => 'rechunk',
                         batch_size => 1,
                       },
        -rc_name      => 'default',
        -can_be_empty  => 1,
        -flow_into => {
                        2 => ['genblast_retry'],
                      },
      },

      {
        -logic_name => 'genblast_retry',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
        -parameters => {
                         iid_type => 'db_seq',
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('genblast_db'),
                         logic_name => 'genblast',
                         module => 'HiveGenblast',
                         genblast_path => $self->o('genblast_path'),
                         genblast_db_path => $self->o('genome_file'),
                         commandline_params => $genblast_params{$self->o('blast_type').'_genome'},
                         sequence_table_name => $self->o('uniprot_table_name'),
                         max_rank => $self->o('genblast_max_rank'),
                         genblast_pid => $self->o('genblast_pid'),
                         timer => '1h',
                         blast_eval => $self->o('genblast_eval'),
                         blast_cov  => $self->o('genblast_cov'),
                         flag_small_introns => $self->o('genblast_flag_small_introns'),
                         flag_subpar_models => $self->o('genblast_flag_subpar_models'),
                       },
        -rc_name          => 'genblast_retry',
        -can_be_empty  => 1,
        -failed_job_tolerance => 100,
        -flow_into => {
                        -1 => ['failed_genblast_proteins'],
                        -2 => ['failed_genblast_proteins'],
                        -3 => ['failed_genblast_proteins'],
                      },
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
      },

      {
        -logic_name => 'failed_genblast_proteins',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name          => 'default',
        -can_be_empty  => 1,
        -failed_job_tolerance => 100,
      },


      {
        -logic_name => 'classify_genblast_models',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport',
        -parameters => {
                         classification_type => 'standard',
                         update_gene_biotype => 1,
                         target_db => $self->o('genblast_db'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['genblast_sanity_checks'],
                      },
      },


      {
        -logic_name => 'genblast_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         target_db => $self->o('genblast_db'),
                         sanity_check_type => 'gene_db_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'gene_db_checks')->{$self->default_options->{'uniprot_set'}}->{'genblast'},
                       },

        -rc_name    => '4GB',
        -flow_into => { 1 => ['create_genblast_slices'] },
      },


      {
        -logic_name => 'create_genblast_slices',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('dna_db'),
                         iid_type => 'slice',
                         slice_size => 10000000,
                         coord_system_name => 'toplevel',
                         include_non_reference => 0,
                         top_level => 1,
                         # These options will create only slices that have a gene on the slice in one of the feature dbs
                         feature_constraint => 1,
                         feature_type => 'gene',
                         feature_dbs => [$self->o('genblast_db')],
                      },
        -flow_into => {
                       '2->A' => ['genblast_select'],
                       'A->1' => ['update_genblast_select_biotypes'],
                      },

        -rc_name    => 'default',
      },


      {
        -logic_name => 'genblast_select',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.catfile($self->o('ensembl_analysis_script'), 'genebuild', 'pick_best_alt_transcripts.pl').
                                     ' -slice_name #iid#'.
                                     ' -source_host '.$self->o('genblast_db','-host').
                                     ' -source_user '.$self->o('user_r').
                                     ' -source_port '.$self->o('genblast_db','-port').
                                     ' -source_dbname '.$self->o('genblast_db','-dbname').
                                     ' -dna_user '.$self->o('user_r').
                                     ' -dna_host '.$self->o('dna_db','-host').
                                     ' -dna_port '.$self->o('dna_db','-port').
                                     ' -dna_dbname '.$self->o('dna_db','-dbname').
                                     ' -out_user '.$self->o('user').
                                     ' -out_pass '.$self->o('password').
                                     ' -out_host '.$self->o('genblast_select_db','-host').
                                     ' -out_port '.$self->o('genblast_select_db','-port').
                                     ' -out_dbname '.$self->o('genblast_select_db','-dbname')
                      },
        -rc_name => '4GB',
     },


      {
        -logic_name => 'update_genblast_select_biotypes',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('genblast_select_db'),
          sql => [
            'UPDATE gene set biotype="genblast_select"',' UPDATE transcript set biotype="genblast_select"'
          ],
        },
        -rc_name    => 'default',
        -flow_into  => {
                         1 => ['classify_genblast_select_models'],
                       },
      },


      {
        -logic_name => 'classify_genblast_select_models',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport',
        -parameters => {
                         classification_type => 'standard',
                         update_gene_biotype => 1,
                         target_db => $self->o('genblast_select_db'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['genblast_select_sanity_checks'],
                      },
      },


      {
        -logic_name => 'genblast_select_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         target_db => $self->o('genblast_db'),
                         sanity_check_type => 'gene_db_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'gene_db_checks')->{$self->default_options->{'uniprot_set'}}->{'genblast_select'},
                       },

        -rc_name    => '4GB',
        -flow_into => {
                        '1->A' => ['create_cdna_db'],
                        'A->1' => ['create_ig_tr_db'],
                      },

      },


######################################################################################
#
# cDNA alignment
#
######################################################################################

    {
      -logic_name => 'create_cdna_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db => $self->o('dna_db'),
        target_db => $self->o('cdna_db'),
        create_type => 'clone',
      },
      -rc_name    => 'default',
      -flow_into => {
        '1->A' => ['create_genewise_db', 'download_mRNA', 'download_selenocysteines'],
        'A->1' => ['create_besttargetted_db'],
      },
    },

    {
      -logic_name => 'create_genewise_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db => $self->o('dna_db'),
        target_db => $self->o('genewise_db'),
        create_type => 'clone',
      },
      -rc_name    => 'default',
      -flow_into => {
        '1->A' => ['download_uniprot_self', 'download_refseq_self'],
        'A->1' => ['load_self'],
      },
    },
    {

      -logic_name => 'download_refseq_self',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadNCBImRNA',
      -parameters => {
        query => '#species_name#[Organism] AND RefSeq[Filter]',
        species_name => $self->o('species_name'),
        output_file => catfile($self->o('targetted_path'), 'ncbi_self.fa'),
        ncbidb => 'protein',
        _branch_to_flow_to => 1,
      },
      -rc_name          => 'default',
      -flow_into => {
        1 => [':////accu?iid=[]'],
      },
    },
    {

      -logic_name => 'download_uniprot_self',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadUniProtFiles',
      -parameters => {
        taxon_id => $self->o('taxon_id'),
        multi_query_download => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::UniProtCladeDownloadStatic', 'self_isoforms_12'),
        output_path => $self->o('targetted_path'),
        _branch_to_flow_to => 1,
      },
      -rc_name          => 'default',
      -flow_into => {
        1 => [':////accu?iid=[]'],
      },
    },

    {

      -logic_name => 'load_self',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProcessUniProtFiles',
      -parameters => {
        output_file => catfile($self->o('targetted_path'), 'proteome.fa'),
        skip_Xs => 5,
        killlist_db => $self->o('killlist_db'),
        killlist_type => 'protein',
        sequence_table_name => $self->o('uniprot_table_name'),
      },
      -rc_name          => 'default',
      -flow_into => {
        '1' => ['indicate_proteome', 'generate_exonerate_jobs'],
      },
    },
    {
      -logic_name => 'generate_exonerate_jobs',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        iid_type => 'sequence_accession',
        sequence_table_name => $self->o('uniprot_table_name'),
        batch_size => 1,
        constraint => 'biotype IN ("self_pe12_sp", "self_pe12_tr", "self_isoforms_12_sp", "self_isoforms_12_tr", "ncbi_self_refseq")',
      },
      -rc_name      => 'default',
      -flow_into => {
        2 => ['targetted_exonerate'],
      },
    },
    {

      -logic_name => 'targetted_exonerate',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
      -parameters => {
        iid_type => 'db_seq',
        sequence_table_name => $self->o('uniprot_table_name'),
        dna_db => $self->o('dna_db'),
        target_db => $self->o('genewise_db'),
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_protein')},
        genome_file      => $self->o('genome_file'),
        exonerate_path   => $self->o('exonerate_path'),
        repeat_libraries => $self->o('repeat_logic_names'),
        calculate_coverage_and_pid => $self->o('target_exonerate_calculate_coverage_and_pid'),
      },
      -rc_name          => 'exonerate',
    },


    {
      -logic_name => 'indicate_proteome',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => '#indicate_path# -d #indicate_dir# -f #proteome# -i #proteome_index# -p singleWordParser',
        indicate_path => $self->o('indicate_path'),
        proteome => 'proteome.fa',
        indicate_dir => $self->o('targetted_path'),
        proteome_index => catdir($self->o('targetted_path'), 'proteome_index'),
      },
      -rc_name => 'default',
      -flow_into => {
        '1' => ['generate_pmatch_jobs'],
      },
    },


    {
      -logic_name => 'generate_pmatch_jobs',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db        => $self->o('dna_db'),
        coord_system_name => 'toplevel',
        iid_type => 'slice',
        top_level => 1,
      },
      -rc_name      => 'default',
      -flow_into => {
        '2->A' => ['pmatch'],
        'A->1' => ['bestpmatch'],
      },
    },
    {

      -logic_name => 'pmatch',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivePmatch',
      -parameters => {
        target_db => $self->o('genewise_db'),
        dna_db => $self->o('dna_db'),
        PROTEIN_FILE => catfile($self->o('targetted_path'), 'proteome.fa'),
        MIN_COVERAGE => 25,
        BINARY_LOCATION => $self->o('pmatch_path'),
        REPEAT_MASKING => [],
        MAX_INTRON_LENGTH => 50000,
        OPTIONS => '-T 20', # set threshold to 14 for more sensitive search
      },
      -rc_name          => 'default',
    },
    {

      -logic_name => 'bestpmatch',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBestPmatch',
      -parameters => {
        source_db => $self->o('genewise_db'),
        target_db => $self->default_options->{'genewise_db'},
        PMATCH_LOGIC_NAME => ['pmatch'],
        MIN_COVERAGE => 50,
      },
      -rc_name          => 'default',
      -flow_into => {
        1 => ['generate_targetted_jobs'],
      },
    },
    {
      -logic_name => 'generate_targetted_jobs',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        iid_type => 'feature_region',
        feature_type => 'protein_align_feature',
        target_db        => $self->o('genewise_db'),
        logic_name => ['bestpmatch'],
        coord_system_name => 'toplevel',
      },
      -rc_name      => 'default',
      -flow_into => {
        2 => ['targetted_genewise_gtag', 'targetted_genewise_nogtag', 'targetted_exo'],
      },
    },
    {

      -logic_name => 'targetted_genewise_gtag',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastMiniGenewise',
      -parameters => {
        source_db => $self->o('genewise_db'),
        target_db => $self->o('genewise_db'),
        killlist_db => $self->o('killlist_db'),
        dna_db => $self->o('dna_db'),
        gtag => 0, # 0 is for gtag, 1 is for non canonical
        biotype => 'gw_gtag',
        max_intron_length => 200000,
        seqfetcher_index => [catfile($self->o('targetted_path'), 'proteome_index')],
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::GeneWiseStatic', 'targetted_genewise')},
      },
      -rc_name          => 'exonerate',
    },
    {

      -logic_name => 'targetted_genewise_nogtag',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastMiniGenewise',
      -parameters => {
        source_db => $self->o('genewise_db'),
        target_db => $self->o('genewise_db'),
        killlist_db => $self->o('killlist_db'),
        dna_db => $self->o('dna_db'),
        gtag => 1, # 0 is for gtag, 1 is for non canonical
        biotype => 'gw_nogtag',
        max_intron_length => 200000,
        seqfetcher_index => [catfile($self->o('targetted_path'), 'proteome_index')],
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::GeneWiseStatic', 'targetted_genewise')},
      },
      -rc_name          => 'exonerate',
    },
    {

      -logic_name => 'targetted_exo',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerateForGenewise',
      -parameters => {
        source_db => $self->o('genewise_db'),
        target_db => $self->o('genewise_db'),
        killlist_db => $self->o('killlist_db'),
        dna_db => $self->o('dna_db'),
        biotype => 'gw_exo',
        seqfetcher_index => [catfile($self->o('targetted_path'), 'proteome_index')],
        max_intron_length => 700000,
        program_file => $self->o('exonerate_path'),
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::GeneWiseStatic', 'targetted_exonerate')},
      },
      -rc_name          => 'exonerate',
    },
    {
      -logic_name => 'download_mRNA',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadNCBImRNA',
      -parameters => {
        output_file => $self->o('cdna_file'),
        filetype => 'gb',
        query => $self->o('ncbi_query'),
      },
      -rc_name    => 'default',
      -flow_into => {
        2 => ['prepare_cdna'],
      },
    },
    {
      -logic_name => 'prepare_cdna',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('prepare_cdnas_script').
          ' -killdbnames '.$self->o('killlist_db','-dbname').
          ' -killdbhost '.$self->o('killlist_db','-host').
          ' -killdbuser '.$self->o('killlist_db','-user').
          ' -killdbport '.$self->o('killlist_db','-port').
          ($self->o('killlist_db','-pass') ? ' -killdbpass '.$self->o('killlist_db','-pass') : '').
          ' -infile #sequence_file#'.
          ' -outfile #sequence_file#.clipped'.
          ($self->o('taxon_id') ? ' -tax_id '.$self->o('taxon_id') : '').
          ' -nomole',
        sequence_file => $self->o('cdna_file'),
      },
      -rc_name    => 'default',
      -flow_into => {
        1 => ['load_cdna_file'],
      },
    },
    {
      -logic_name => 'load_cdna_file',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadmRNAs',
      -parameters => {
        sequence_table_name => $self->o('cdna_table_name'),
        filetype => 'fasta',
        sequence_file => $self->o('cdna_file').'.clipped',
      },
      -flow_into => {
        '2->A' => ['exonerate'],
        'A->1' => ['prepare_cdna2genome'],
      },
      -rc_name    => 'default',
    },

    {
      -logic_name => 'exonerate',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
      -rc_name    => 'exonerate',
      -parameters => {
        iid_type => 'db_seq',
        sequence_table_name => $self->o('cdna_table_name'),
        dna_db => $self->o('dna_db'),
        target_db => $self->o('cdna_db'),
        logic_name => $self->o('exonerate_logic_name'),
        genome_file      => $self->o('genome_file'),
        exonerate_path   => $self->o('exonerate_path'),
        repeat_libraries => $self->o('repeat_logic_names'),
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','cdna_est2genome')},
        exonerate_cdna_pid => 50,
        exonerate_cdna_cov => 50,
        calculate_coverage_and_pid => 0,
      },
      -batch_size => 100,
      -flow_into => {
        -1 => ['exonerate_retry'],
      },
      -batch_size => 100,
      -failed_job_tolerance => 5,
    },

    {
      -logic_name => 'exonerate_retry',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
      -rc_name    => 'exonerate_6G',
      -parameters => {
        iid_type => 'db_seq',
        sequence_table_name => $self->o('cdna_table_name'),
        dna_db => $self->o('dna_db'),
        target_db => $self->o('cdna_db'),
        logic_name => $self->o('exonerate_logic_name'),
        genome_file      => $self->o('genome_file'),
        exonerate_path   => $self->o('exonerate_path'),
        repeat_libraries => $self->o('repeat_logic_names'),
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','cdna_est2genome')},
        exonerate_cdna_pid => 50,
        exonerate_cdna_cov => 50,
        calculate_coverage_and_pid => 0,
      },
      -batch_size => 100,
      -failed_job_tolerance => 100,
      -can_be_empty => 1,
    },
    {
      -logic_name => 'prepare_cdna2genome',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('prepare_cdnas_script').
          ' -killdbnames '.$self->o('killlist_db','-dbname').
          ' -killdbhost '.$self->o('killlist_db','-host').
          ' -killdbuser '.$self->o('killlist_db','-user').
          ' -killdbport '.$self->o('killlist_db','-port').
          ($self->o('killlist_db','-pass') ? ' -killdbpass '.$self->o('killlist_db','-pass') : '').
          ' -infile #sequence_file#'.
          ' -outfile #sequence_file#.cds'.
          ' -annotation '.$self->o('annotation_file').
          ($self->o('taxon_id') ? ' -tax_id '.$self->o('taxon_id') : '').
          ' -nomole',
        sequence_file => $self->o('cdna_file'),
      },
      -rc_name    => 'default',
      -flow_into => {
        1 => ['create_cdna2genome_db'],
      },
    },
    {
      -logic_name => 'create_cdna2genome_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db => $self->o('dna_db'),
        target_db => $self->o('cdna2genome_db'),
        create_type => 'clone',
      },
      -rc_name    => 'default',
      -flow_into => {
        '1' => ['create_cdna_toplevel_slices'],
      },
    },
    {
      -logic_name => 'create_cdna_toplevel_slices',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db => $self->o('cdna_db'),
        iid_type => 'slice',
        coord_system_name => 'toplevel',
        include_non_reference => 0,
        top_level => 1,
        feature_constraint => 1,
        feature_type => 'gene',
      },
      -flow_into => {
        '2' => ['apply_threshold'],
      },
      -rc_name    => 'default',
    },
    {
      -logic_name => 'apply_threshold',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSelectGeneOnFilter',
      -parameters => {
        dna_db => $self->o('dna_db'),
        source_db => $self->o('cdna_db'),
        logic_name => $self->o('species_name').'_cdna',
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','cdna_est2genome')},
        exonerate_cdna_pid => $self->o('cdna_selection_pid'),
        exonerate_cdna_cov => $self->o('cdna_selection_cov'),
      },
      -rc_name => 'default',
      -analysis_capacity => 5,
      -batch_size => 10,
      -flow_into => {
        '1' => ['create_cdna2genome_slices'],
      },
    },
    {
      -logic_name => 'create_cdna2genome_slices',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db => $self->o('cdna_db'),
        iid_type => 'feature_region',
        feature_type => 'gene',
        logic_name => [$self->o('exonerate_logic_name')],
        coord_system_name => 'toplevel',
        include_non_reference => 0,
        top_level => 1,
        use_annotation => 1,
# These options will create only slices that have a gene on the slice in one of the feature dbs
        annotation_file => $self->o('annotation_file'),
        region_padding => $self->o('cdna2genome_region_padding'),
      },
      -flow_into => {
        '2->A' => ['cdna2genome'],
        'A->1' => ['internal_stop'],
      },

      -rc_name    => 'default',
    },
    {
      -logic_name => 'cdna2genome',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2GenesRegion',
      -rc_name    => 'exonerate',
      -parameters => {
        iid_type => 'db_seq',
        dna_db => $self->o('dna_db'),
        sequence_table_name => $self->o('cdna_table_name'),
        source_db => $self->o('cdna_db'),
        target_db => $self->o('cdna2genome_db'),
        logic_name => 'cdna2genome',
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','cdna2genome')},
        calculate_coverage_and_pid => 1,
        exonerate_path => $self->o('exonerate_annotation'),
        annotation_file => $self->o('annotation_file'),
        SOFT_MASKED_REPEATS => $self->o('repeat_logic_names'),
      },
      -batch_size => 10,
      -flow_into => {
        '-1' => ['cdna2genome_himem'],
      },
    },

    {
      -logic_name => 'cdna2genome_himem',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2GenesRegion',
      -rc_name    => 'exonerate_6G',
      -parameters => {
        iid_type => 'db_seq',
        dna_db => $self->o('dna_db'),
        sequence_table_name => $self->o('cdna_table_name'),
        source_db => $self->o('cdna_db'),
        target_db => $self->o('cdna2genome_db'),
        logic_name => 'cdna2genome',
        module     => 'HiveExonerate2Genes',
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','cdna2genome')},
        calculate_coverage_and_pid => 1,
        program => $self->o('exonerate_annotation'),
        annotation_file => $self->o('annotation_file'),
        repeat_libraries => $self->o('repeat_logic_names'),
      },
      -batch_size => 10,
    },
    {
      -logic_name => 'internal_stop',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveInternalStopFix',
      -parameters => {
        dna_db => $self->o('dna_db'),
        source_db => $self->o('cdna2genome_db'),
        edited_biotype => 'edited',
        stop_codon_biotype => 'stop_codon',
        logic_name => 'cdna2genome',
        biotype => undef,
        source => undef,
      },
      -rc_name    => 'default',
    },
    {

      -logic_name => 'download_selenocysteines',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadUniProtFiles',
      -parameters => {
        taxon_id => $self->o('taxon_id'),
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::UniProtCladeDownloadStatic', 'selenocysteine')},
        output_path => $self->o('targetted_path'),
      },
      -rc_name          => 'default',
      -flow_into => {
        2 => ['load_selenocysteine'],
      },
    },
    {

      -logic_name => 'load_selenocysteine',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProcessUniProtFiles',
      -parameters => {
        sequence_table_name => $self->o('uniprot_table_name'),
      },
      -rc_name          => 'default',
      -flow_into => {
        2 => ['process_selenocysteine'],
      },
    },
    {

      -logic_name => 'process_selenocysteine',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSelenocysteineFinder',
      -parameters => {
        target_db => $self->o('genewise_db'),
        dna_db => $self->o('dna_db'),
        genome => $self->o('genome_file'),
        biotype => 'seleno_self',
        exonerate => $self->o('exonerate_path'),
        genewise => $self->o('genewise_path'),
        iid_type => 'db_seq',
        sequence_table_name => $self->o('uniprot_table_name'),
      },
      -rc_name          => 'default',
    },
    {
      -logic_name => 'create_besttargetted_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db => $self->o('dna_db'),
        target_db => $self->o('best_targeted_db'),
        create_type => 'clone',
      },
      -rc_name    => 'default',
      -flow_into => {
        1 => ['generate_besttargetted_index'],
      },
    },
    {
      -logic_name => 'generate_besttargetted_index',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenerateBestTargettedIndex',
      -parameters => {
        source_db        => $self->o('cdna2genome_db'),
        seqfetcher_index => [catfile($self->o('targetted_path'), 'proteome_index')],
        fasta_filename => catfile($self->o('targetted_path'), 'best_targetted.fa'),
        email => $self->o('email_address'),
        genbank_file => $self->o('cdna_file'),
        files => [catfile($self->o('targetted_path'), 'proteome.fa')],
      },
      -rc_name      => 'default',
      -flow_into => {
        1 => ['indicate_BT'],
      },
    },
    {
      -logic_name => 'indicate_BT',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => '#indicate_path# -d #indicate_dir# -f #proteome# -i #indicate_dir#/best_targetted_index -p singleWordParser -M BTMultiParser',
        indicate_path => $self->o('indicate_path'),
        proteome => 'best_targetted.fa',
        indicate_dir => $self->o('targetted_path'),
      },
      -rc_name => 'default',
      -flow_into => {
        1 => ['generate_besttargetted_jobs'],
      },
    },
    {
      -logic_name => 'generate_besttargetted_jobs',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db        => $self->o('genewise_db'),
        coord_system_name => 'toplevel',
        iid_type => 'slice',
        feature_constraint => 1,
        feature_type => 'gene',
        top_level => 1,
        feature_dbs => [$self->default_options->{'genewise_db'}, $self->o('cdna2genome_db')],
      },
      -rc_name      => 'default',
      -flow_into => {
        2 => ['best_targetted'],
      },
    },
    {

      -logic_name => 'best_targetted',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBestTargetted',
      -parameters => {
        target_db => $self->o('best_targeted_db'),
        dna_db => $self->o('dna_db'),
        source_db => { protein_db => $self->o('genewise_db'), cdna2genome_db => $self->o('cdna2genome_db')},
        SEQFETCHER_DIR => [catfile($self->o('targetted_path'), 'proteome_index'),
        catfile($self->o('targetted_path'), 'best_targetted_index')],
        SEQFETCHER_OBJECT => 'Bio::EnsEMBL::Analysis::Tools::SeqFetcher::OBDAIndexSeqFetcher',
        INPUT_DATA_FROM_DBS  => {
          protein_db => ['seleno_self', 'gw_gtag', 'gw_nogtag', 'gw_exo', 'targetted_exonerate'],
          cdna2genome_db => ['cdna2genome', 'edited'],
        },
        BIOTYPES => ['seleno_self', 'cdna2genome', 'edited', 'gw_gtag', 'gw_nogtag', 'gw_exo', 'targetted_exonerate'], # sorted list, preferred is first
        protein_min_coverage => $self->o('best_targetted_min_coverage'),
        protein_min_identity => $self->o('best_targetted_min_identity'),
      },
      -rc_name          => 'exonerate',
    },


######################################################################################
#
# IG and TR genes
#
######################################################################################

      {
        -logic_name => 'create_ig_tr_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                     source_db => $self->o('dna_db'),
                     target_db => $self->o('ig_tr_db'),
                     create_type => 'clone',
                   },
         -flow_into => {
                         1 => ['load_ig_tr_seqs'],
                       },
      },


      {
        -logic_name => 'load_ig_tr_seqs',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('load_fasta_script_path')
                                ." -dbhost ".$self->o('pipeline_db','-host')
                                ." -dbuser ".$self->o('pipeline_db','-user')
                                ." -dbpass ".$self->o('pipeline_db','-pass')
                                ." -dbname ".$self->o('pipeline_db','-dbname')
                                ." -dbport ".$self->o('pipeline_db','-port')
                                ." -fasta_file ".$self->o('ig_tr_blast_path')."/".$self->o('ig_tr_fasta_file')
                                ." -sequence_table_name ".$self->o('ig_tr_table_name')
                                ." -create_table 1",
                       },

         -rc_name => 'default',
         -flow_into => {
                         1 => ['generate_ig_tr_jobs'],
                       },
      },


      {
        -logic_name => 'generate_ig_tr_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         iid_type => 'sequence_accession',
                         batch_size => $self->o('ig_tr_batch_size'),
                         sequence_table_name => $self->o('ig_tr_table_name'),
                       },
        -rc_name      => 'default',
        -flow_into => {
                        '2->A' => ['ig_tr_genblast'],
                        'A->1' => ['update_ig_tr_hitnames'],
                      },
      },


      {
        -logic_name => 'ig_tr_genblast',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
        -parameters => {
                         iid_type => 'db_seq',
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('ig_tr_db'),
                         logic_name => 'ig_tr_gene',
                         module => 'HiveGenblast',
                         genblast_path => $self->o('genblast_path'),
                         genblast_db_path => $self->o('genome_file'),
                         commandline_params => $genblast_params{$self->o('blast_type').'_genome'},
                         sequence_table_name => $self->o('ig_tr_table_name'),
                         max_rank => $self->o('ig_tr_genblast_max_rank'),
                         genblast_pid => $self->o('ig_tr_genblast_pid'),
                         timer => '2h',
                         blast_eval => $self->o('ig_tr_genblast_eval'),
                         blast_cov  => $self->o('ig_tr_genblast_cov'),
                       },
        -rc_name    => 'genblast',
        -flow_into => {
                        -1 => ['split_ig_tr_genblast_jobs'],
                        -2 => ['split_ig_tr_genblast_jobs'],
                        -3 => ['split_ig_tr_genblast_jobs'],
                      },
      },


      {
        -logic_name => 'split_ig_tr_genblast_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         iid_type => 'rechunk',
                         batch_size => 1,
                       },
        -rc_name      => 'default',
        -can_be_empty  => 1,
        -flow_into => {
                        2 => ['ig_tr_genblast_retry'],
                      },
      },


      {
        -logic_name => 'ig_tr_genblast_retry',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
        -parameters => {
                         iid_type => 'db_seq',
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('ig_tr_db'),
                         logic_name => 'genblast',
                         module => 'HiveGenblast',
                         genblast_path => $self->o('genblast_path'),
                         genblast_db_path => $self->o('genome_file'),
                         commandline_params => $genblast_params{$self->o('blast_type').'_genome'},
                         sequence_table_name => $self->o('ig_tr_table_name'),
                         max_rank => $self->o('genblast_max_rank'),
                         genblast_pid => $self->o('genblast_pid'),
                         timer => '1h',
                         blast_eval => $self->o('ig_tr_genblast_eval'),
                         blast_cov  => $self->o('ig_tr_genblast_cov'),
                       },
        -rc_name          => 'genblast_retry',
        -can_be_empty  => 1,
        -failed_job_tolerance => 100,
        -flow_into => {
                        -1 => ['failed_ig_tr_genblast_proteins'],
                        -2 => ['failed_ig_tr_genblast_proteins'],
                        -3 => ['failed_ig_tr_genblast_proteins'],
                      },
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
      },


      {
        -logic_name => 'failed_ig_tr_genblast_proteins',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name          => 'default',
        -can_be_empty  => 1,
        -failed_job_tolerance => 100,
      },


      {
        -logic_name => 'update_ig_tr_hitnames',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          sql => 'UPDATE protein_align_feature set hit_name = "IMGT"',
          db_conn => $self->o('ig_tr_db'),
        },
        -rc_name    => 'default',
        -flow_into => {
                        '1' => ['ig_tr_sanity_checks'],
                      },
      },


      {
        -logic_name => 'ig_tr_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         target_db => $self->o('ig_tr_db'),
                         sanity_check_type => 'gene_db_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'gene_db_checks')->{$self->default_options->{'uniprot_set'}}->{'ig_tr'},
                       },

        -rc_name    => '1.5GB',
        -flow_into => {
          1 => ['create_projection_coding_db'],
        },
      },

########################################################################
#
# Projection analyses
#
########################################################################
      {
        -logic_name => 'create_projection_coding_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
          source_db   => $self->o('dna_db'),
          target_db   => $self->o('projection_coding_db'),
          create_type => 'clone',
        },
        -rc_name    => 'default',
        -flow_into  => {
          '1->A' => ['fan_projection'],
          'A->1' => ['cluster_ig_tr_genes'],
        },
      },


      {
        -logic_name => 'fan_projection',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'if [ "#skip_projection#" -ne "0" ]; then exit 42; else exit 0;fi',
          return_codes_2_branches => {'42' => 2},
        },
        -rc_name    => 'default',
        -flow_into  => {
          '1->A' => ['create_projection_coding_input_ids'],
          'A->1' => ['classify_projection_coding_models'],
        },
      },


      {
        -logic_name => 'create_projection_coding_input_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
          target_db           => $self->o('projection_source_db'),
          iid_type            => 'feature_id',
          feature_type        => 'transcript',
          feature_restriction => 'projection',
          biotypes            => {
            'protein_coding' => 1,
          },
          batch_size          => 500,
        },
        -rc_name    => '4GB',
        -flow_into => {
          2 => ['project_coding_transcripts'],
        },
      },


      {
        -logic_name => 'project_coding_transcripts',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveWGA2GenesDirect',
        -parameters => {
          logic_name => 'project_transcripts',
          module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveWGA2GenesDirect',
          source_dna_db         => $self->default_options()->{'projection_source_db'},
          target_dna_db         => $self->o('dna_db'),
          source_transcript_db  => $self->default_options()->{'projection_source_db'},
          target_transcript_db   => $self->o('projection_coding_db'),
          compara_db            => $self->o('projection_lastz_db'),
          method_link_type => 'LASTZ_NET',
          max_exon_readthrough_dist => 15,
          TRANSCRIPT_FILTER => {
            OBJECT     => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
            PARAMETERS => {
              -coverage => $self->o('projection_cov'),
              -percent_id => $self->o('projection_pid'),
            },
          },
          iid_type => 'feature_id',
          feature_type => 'transcript',
          calculate_coverage_and_pid => $self->o('projection_calculate_coverage_and_pid'),
          max_internal_stops => $self->o('projection_max_internal_stops'),
          timer => '30m',
        },
        -flow_into => {
          -3 => ['failed_projection_coding_jobs'],
        },
        -rc_name    => 'project_transcripts',
        -batch_size => 100,
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
      },


      {
        -logic_name => 'failed_projection_coding_jobs',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {},
        -rc_name    => 'default',
        -can_be_empty  => 1,
        -failed_job_tolerance => 100,
      },


      {
        -logic_name => 'classify_projection_coding_models',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport',
        -parameters => {
          skip_analysis => $self->o('skip_projection'),
          classification_type => 'standard',
          update_gene_biotype => 1,
          target_db => $self->o('projection_coding_db'),
        },
        -rc_name    => 'default',
        -flow_into => {
          1 => ['fix_unaligned_protein_hit_names'],
#          When the realign part is fix, the line above needs to be deleted and the lines below uncommented
#          All analysis below create_projection_realign_db need to be uncommented too
#          '1->A' => ['fix_unaligned_protein_hit_names'],
#          'A->1' => ['create_projection_realign_db'],
        },
      },


      {
        -logic_name => 'fix_unaligned_protein_hit_names',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn    => $self->o('projection_coding_db'),
          sql => [
            'UPDATE protein_align_feature JOIN transcript_supporting_feature ON feature_id = protein_align_feature_id'.
              ' JOIN transcript USING(transcript_id) SET hit_name = stable_id',
            'UPDATE protein_align_feature JOIN supporting_feature ON feature_id = protein_align_feature_id'.
              ' JOIN exon_transcript USING(exon_id) JOIN transcript USING(transcript_id) SET hit_name = stable_id',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
          1 => ['projection_coding_sanity_checks'],
        },
      },


      {
        -logic_name => 'projection_coding_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
          target_db => $self->o('projection_coding_db'),
          sanity_check_type => 'gene_db_checks',
          min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
            'gene_db_checks')->{$self->default_options->{'uniprot_set'}}->{'projection_coding'},
        },
        -rc_name    => '4GB',
      },


#      When the realign part is fixed, uncomment the whole block
#      {
#        -logic_name => 'create_projection_realign_db',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
#        -parameters => {
#          source_db => $self->o('dna_db'),
#          target_db => $self->o('projection_realign_db'),
#          create_type => 'clone',
#        },
#        -rc_name    => 'default',
#        -flow_into => {
#          '1->A' => ['create_ids_for_evaluate_projection'],
#          'A->1' => ['update_realign_biotypes'],
#        },
#      },
#
#
#      {
#        -logic_name => 'create_ids_for_evaluate_projection',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
#        -parameters => {
#          target_db => $self->o('projection_coding_db'),
#          iid_type => 'feature_id',
#          feature_type => 'transcript',
#          batch_size => 500,
#        },
#        -rc_name    => 'default',
#        -flow_into => {
#          2 => ['evaluate_coding_transcripts'],
#        },
#      },
#
#
#      {
#        -logic_name => 'evaluate_coding_transcripts',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProjectionRealignment',
#        -parameters => {
#          dna_db => $self->o('dna_db'),
#          projection_db => $self->o('projection_coding_db'),
#          projection_source_db => $self->o('projection_source_db'),
#          projection_realign_db => $self->o('projection_realign_db'),
#          protein_table_name => $self->o('realign_table_name'),
#          max_feature_issues => $self->o('max_projection_structural_issues'),
#        },
#        -rc_name    => 'default',
#        -flow_into  => {
#          2 => ['realign_projection'],
#        },
#      },
#
#
#      {
#        -logic_name => 'realign_projection',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
#        -parameters => {
#          iid_type => 'projection_transcript_id',
#          sequence_table_name => 'projection_source_sequences',
#          projection_padding => 0,
#          dna_db => $self->o('dna_db'),
#          target_db => $self->o('projection_realign_db'),
#          projection_db => $self->o('projection_coding_db'),
#          logic_name => 'genblast',
#          module => 'HiveGenblast',
#          genblast_path => $self->o('genblast_path'),
#          commandline_params => $genblast_params{$self->o('blast_type').'_projection'},
#          query_seq_dir => $self->o('homology_models_path'),
#          max_rank => 1,
#          genblast_pid => $self->o('genblast_pid'),
#          blast_eval => $self->o('genblast_eval'),
#          blast_cov => $self->o('genblast_cov'),
#          timer => '30m',
#        },
#        -rc_name    => 'default',
#        -flow_into  => {
#          -1 => ['failed_realignments'],
#          -2 => ['failed_realignments'],
#          -3 => ['failed_realignments'],
#        },
#      },
#
#
#      {
#        -logic_name => 'failed_realignments',
#        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
#        -parameters => {},
#        -rc_name    => 'default',
#        -can_be_empty  => 1,
#      },
#
#
#      {
#        -logic_name => 'update_realign_biotypes',
#        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
#        -parameters => {
#          db_conn    => $self->o('projection_realign_db'),
#          sql => [
#            'UPDATE gene SET biotype = "realign"',
#            'UPDATE transcript SET biotype = "realign"',
#          ],
#        },
#        -rc_name    => 'default',
#        -flow_into => {
#          1 => ['classify_realigned_models'],
#        },
#      },
#
#
#      {
#        -logic_name => 'classify_realigned_models',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport',
#        -parameters => {
#          classification_type => 'standard',
#          update_gene_biotype => 1,
#          target_db => $self->o('projection_realign_db'),
#        },
#        -rc_name    => 'default',
#        -flow_into => {
#          1 => ['realign_sanity_checks'],
#        },
#      },
#
#
#      {
#        -logic_name => 'realign_sanity_checks',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
#        -parameters => {
#                         target_db => $self->o('projection_realign_db'),
#                         sanity_check_type => 'gene_db_checks',
#                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
#                                                                             'gene_db_checks')->{$self->default_options->{'uniprot_set'}}->{'realign'},
#                       },
#
#        -rc_name    => '4GB',
#        -flow_into => {
#                       1 => ['create_projection_lincrna_db'],
#                      },
#      },
#
#
#      {
#        -logic_name => 'create_projection_lincrna_db',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
#        -parameters => {
#                         source_db => $self->o('dna_db'),
#                         target_db => $self->o('projection_lincrna_db'),
#                         create_type => 'clone',
#                       },
#        -flow_into => {
#                       1 => ['create_projection_pseudogene_db'],
#                      },
#      },
#
#
#      {
#        -logic_name => 'create_projection_pseudogene_db',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
#        -parameters => {
#                         source_db => $self->o('dna_db'),
#                         target_db => $self->o('projection_pseudogene_db'),
#                         create_type => 'clone',
#                       },
#
#        -flow_into => {
#                        1 => ['generate_projection_lincrna_ids','generate_projection_pseudogene_ids','generate_projection_ig_tr_ids'],
#                      },
#      },
#
#
#      {
#        -logic_name => 'generate_projection_lincrna_ids',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
#        -parameters => {
#                         target_db => $self->o('projection_source_db'),
#                         iid_type => 'feature_id',
#                         feature_type => 'transcript',
#                         feature_restriction => 'biotype',
#                         biotypes => {
#                                       'lincRNA' => 1,
#                                     },
#                         batch_size => 100,
#        },
#        -flow_into => {
#                       '2->A' => ['run_project_lincrnas'],
#                       'A->1' => ['project_lincrna_sanity_checks'],
#                      },
#        -rc_name    => 'default',
#      },
#
#
#      {
#        -logic_name => 'run_project_lincrnas',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProjectionExonerate',
#        -parameters => {
#                         'logic_name'                 => 'project_lincrna',
#                         'source_dna_db'              => $self->o('projection_source_db'),
#                         'target_dna_db'              => $self->o('dna_db'),
#                         'source_db'                  => $self->o('projection_source_db'),
#                         'target_db'                  => $self->o('projection_lincrna_db'),
#                         'compara_db'                 => $self->o('projection_lastz_db'),
#                         'method_link_type'           => 'LASTZ_NET',
#                         'exon_region_padding'        => $self->o('projection_exonerate_padding'),
#                         'exonerate_path'             => $self->o('exonerate_path'),
#                         'exonerate_coverage'         => $self->o('projection_lincrna_coverage'),
#                         'exonerate_percent_id'       => $self->o('projection_lincrna_percent_id'),
#                         'calculate_coverage_and_pid' => 0,
#                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_projection_dna')},
#			},
#        -rc_name    => 'default',
#        -hive_capacity => 900,
#      },
#
#
#      {
#        -logic_name => 'project_lincrna_sanity_checks',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
#        -parameters => {
#                         target_db => $self->o('projection_lincrna_db'),
#                         sanity_check_type => 'gene_db_checks',
#                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
#                                                                             'gene_db_checks')->{$self->default_options->{'uniprot_set'}}->{'projection_lincrna'},
#                       },
#        -rc_name    => '4GB',
#      },
#
#
#      {
#        -logic_name => 'generate_projection_pseudogene_ids',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
#        -parameters => {
#                         target_db => $self->o('projection_source_db'),
#                         iid_type => 'feature_id',
#                         feature_type => 'transcript',
#                         feature_restriction => 'biotype',
#                         biotypes => {
#                                       'pseudogene' => 1,
#                                       'processed_pseudogene' => 1,
#                                       'unprocessed_pseudogene' => 1,
#                                     },
#                         batch_size => 100,
#        },
#        -flow_into => {
#                       '2->A' => ['run_project_pseudogenes'],
#                       'A->1' => ['project_pseudogene_sanity_checks'],
#        },
#        -rc_name    => 'default',
#      },
#
#
#      {
#        -logic_name => 'run_project_pseudogenes',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProjectionExonerate',
#        -parameters => {
#                         'logic_name'                 => 'project_pseudogene',
#                         'source_dna_db'              => $self->o('projection_source_db'),
#                         'target_dna_db'              => $self->o('dna_db'),
#                         'source_db'                  => $self->o('projection_source_db'),
#                         'target_db'                  => $self->o('projection_pseudogene_db'),
#                         'compara_db'                 => $self->o('projection_lastz_db'),
#                         'method_link_type'           => 'LASTZ_NET',
#                         'exon_region_padding'        => $self->o('projection_exonerate_padding'),
#                         'exonerate_path'             => $self->o('exonerate_path'),
#                         'exonerate_coverage'         => $self->o('projection_pseudogene_coverage'),
#                         'exonerate_percent_id'       => $self->o('projection_pseudogene_percent_id'),
#                         'calculate_coverage_and_pid' => 0,
#                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_projection_pseudogene')},
#			},
#        -rc_name    => 'default',
#        -hive_capacity => 900,
#      },
#
#
#      {
#        -logic_name => 'project_pseudogene_sanity_checks',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
#        -parameters => {
#                         target_db => $self->o('projection_pseudogene_db'),
#                         sanity_check_type => 'gene_db_checks',
#                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
#                                                                             'gene_db_checks')->{$self->default_options->{'uniprot_set'}}->{'projection_pseudogene'},
#                       },
#        -rc_name    => '4GB',
#      },
#
#
#      {
#        -logic_name => 'generate_projection_ig_tr_ids',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
#        -parameters => {
#                         target_db => $self->o('projection_source_db'),
#                         iid_type => 'feature_id',
#                         feature_type => 'transcript',
#                         feature_restriction => 'biotype',
#                         biotypes => {
#                                       'IG_V_gene' => 1,
#                                       'IG_D_gene' => 1,
#                                       'IG_J_gene' => 1,
#                                       'IG_C_gene' => 1,
#                                       'TR_V_gene' => 1,
#                                       'TR_J_gene' => 1,
#                                       'TR_D_gene' => 1,
#                                       'TR_C_gene' => 1,
#                                     },
#                         batch_size => 100,
#        },
#        -flow_into => {
#                       '2->A' => ['run_project_ig_tr'],
#                       'A->1' => ['project_ig_tr_sanity_checks'],
#        },
#        -rc_name    => 'default',
#      },
#
#
#      {
#        -logic_name => 'run_project_ig_tr',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProjectionExonerate',
#        -parameters => {
#                         'logic_name'                 => 'project_ig_tr',
#                         'source_dna_db'              => $self->o('projection_source_db'),
#                         'target_dna_db'              => $self->o('dna_db'),
#                         'source_db'                  => $self->o('projection_source_db'),
#                         'target_db'                  => $self->o('ig_tr_db'),
#                         'compara_db'                 => $self->o('projection_lastz_db'),
#                         'method_link_type'           => 'LASTZ_NET',
#                         'exon_region_padding'        => $self->o('projection_exonerate_padding'),
#                         'exonerate_path'             => $self->o('exonerate_path'),
#                         'exonerate_coverage'         => $self->o('projection_ig_tr_coverage'),
#                         'exonerate_percent_id'       => $self->o('projection_ig_tr_percent_id'),
#                         'calculate_coverage_and_pid' => 1,
#                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_projection_ig_tr_protein')},
#			},
#        -rc_name    => 'default',
#        -hive_capacity => 900,
#      },
#
#
#      {
#        -logic_name => 'project_ig_tr_sanity_checks',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
#        -parameters => {
#                         target_db => $self->o('ig_tr_db'),
#                         sanity_check_type => 'gene_db_checks',
#                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
#                                                                             'gene_db_checks')->{$self->default_options->{'uniprot_set'}}->{'projection_ig_tr'},
#                       },
#        -rc_name    => '1.5GB',
#      },


      {
        -logic_name => 'cluster_ig_tr_genes',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCollapseIGTR',
        -parameters => {
                         target_db => $self->o('ig_tr_db'),
                         dna_db => $self->o('dna_db'),
                         logic_name => 'ig_tr_gene',
                         logic_names_to_cluster => ['ig_tr_gene','ig_tr_gene_not_best'],
                       },
        -rc_name    => 'genblast',
        -flow_into => {
                        1 => ['update_ig_tr_biotypes'],
                      },
      },


      {
        -logic_name => 'update_ig_tr_biotypes',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('ig_tr_db'),
          sql => [
            'UPDATE transcript JOIN analysis USING(analysis_id) SET biotype = CONCAT(biotype, "_pre_collapse")'.
              ' WHERE logic_name != "ig_tr_collapse"',
            'UPDATE gene JOIN transcript USING(gene_id) SET gene.biotype = transcript.biotype',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
          1 => ['create_rnaseq_for_layer_db'],
        },
      },

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
          'A->1' => ['create_layering_output_db'],
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
          1 => ['copy_rnaseq_blast_db'],
        },
      },


      {
        -logic_name => 'copy_rnaseq_blast_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('rnaseq_blast_db'),
                         target_db => $self->o('rnaseq_for_layer_db'),
                         create_type => 'copy',
                         force_drop => 1,
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
                         skip_analysis => $self->o('skip_rnaseq'),
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
                         skip_check => $self->o('skip_rnaseq'),
                         target_db => $self->o('rnaseq_for_layer_db'),
                         sanity_check_type => 'gene_db_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'gene_db_checks')->{$self->default_options->{'uniprot_set'}}->{'rnaseq_blast'},
                       },

        -rc_name    => '4GB',
      },

############################################################################
#
# Finalisation analyses
#
############################################################################
      {
        -logic_name => 'create_layering_output_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('layering_db'),
                         create_type => 'clone',
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['create_utr_db'],
                      },
      },


      {
        -logic_name => 'create_utr_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('utr_db'),
                         create_type => 'clone',
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['create_genebuilder_db'],
                      },
      },


      {
        -logic_name => 'create_genebuilder_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('genebuilder_db'),
                         create_type => 'clone',
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['create_toplevel_slices'],
                      },
      },


      {
        -logic_name => 'create_toplevel_slices',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('dna_db'),
                         iid_type => 'slice',
                         coord_system_name => 'toplevel',
                         include_non_reference => 0,
                         top_level => 1,
                         # These options will create only slices that have a gene on the slice in one of the feature dbs
                         feature_constraint => 1,
                         feature_type => 'gene',
                         feature_dbs => [$self->o('genblast_db'),$self->o('projection_realign_db'),$self->o('rnaseq_for_layer_db')],
                      },
        -flow_into => {
                       '2->A' => ['layer_annotation'],
                       'A->1' => ['layer_annotation_sanity_checks'],
                      },

        -rc_name    => 'default',
      },


      {
        -logic_name => 'layer_annotation',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLayerAnnotation',
        -parameters => {
                         dna_db     => $self->o('dna_db'),
                         logic_name => 'layer_annotation',
                         module     => 'HiveLayerAnnotation',
                         TARGETDB_REF => $self->o('layering_db'),
                         SOURCEDB_REFS => $self->o('layering_input_gene_dbs'),

                         # Filtering is using done at the exon-overlap level
                         # When no FILTER exists in this file, this is the default behaviour
                         # If you would like to filter in a different way, please specify filter
                         #FILTER => 'Bio::EnsEMBL::Analysis::Tools::GenomeOverlapFilter',
                         #FILTER => 'Bio::EnsEMBL::Analysis::Tools::AllExonOverlapFilter',
                         FILTER => 'Bio::EnsEMBL::Analysis::Tools::CodingExonOverlapFilter',
                         # ordered list of annotation layers. Genes from lower layers
                         # are only retained if they do not "interfere" with genes from
                         # higher layers. Genes in "Discard" layers are when assessing
                         # interference, but are not written to the final database
                         LAYERS => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::LayerAnnotationStatic', $self->default_options->{'uniprot_set'}, undef, 'ARRAY'),
                       },
        -rc_name    => 'layer_annotation',
        -flow_into  => {
                         '1->A' => ['split_slices_on_intergenic'],
                         'A->1' => ['genebuilder'],
                       },
      },

      {
        -logic_name => 'split_slices_on_intergenic',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveFindIntergenicRegions',
        -parameters => {
                         dna_db => $self->o('dna_db'),
                         input_gene_dbs => $self->o('utr_acceptor_dbs'),
                         iid_type => 'slice',
                       },
        -batch_size => 100,
        -hive_capacity => $self->hive_capacity_classes->{'hc_medium'},
        -rc_name    => 'transcript_finalisation',
        -flow_into => {
                        2 => ['run_utr_addition'],
                      },
      },


      {
        -logic_name => 'run_utr_addition',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveUTRAddition',
        -parameters => {
                         logic_name => 'utr_addition',
                         dna_db => $self->o('dna_db'),
                         donor_dbs => $self->o('utr_donor_dbs'),
                         acceptor_dbs => $self->o('utr_acceptor_dbs'),
                         utr_biotype_priorities => $self->o('utr_biotype_priorities'),
                         target_db => $self->o('utr_db'),
                         iid_type => 'slice',
                       },
        -batch_size => 20,
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -rc_name    => 'transcript_finalisation',
     },


     {
        -logic_name => 'genebuilder',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGeneBuilder',
        -parameters => {
                         layering_output_db => $self->o('utr_db'),
                         genebuilder_output_db => $self->o('genebuilder_db'),
                         dna_db     => $self->o('dna_db'),
                         logic_name => 'ensembl',
                         module     => 'HiveGeneBuilder',
                         INPUT_GENES => {
                           'input_db' => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::GenebuilderStatic',
                                                               $self->default_options->{'uniprot_set'}, undef, 'ARRAY'),
                         },
                         OUTPUT_BIOTYPE => 'protein_coding',
                         MAX_TRANSCRIPTS_PER_CLUSTER => 10,
                         MIN_SHORT_INTRON_LEN => 7, #introns shorter than this seem
                         #to be real frame shifts and shoudn't be ignored
                         MAX_SHORT_INTRON_LEN => 15,
                         BLESSED_BIOTYPES => {
                                              'ccds_gene' => 1,
                                              'IG_C_gene' => 1,
                                              'IG_J_gene' => 1,
                                              'IG_V_gene' => 1,
                                              'IG_D_gene' => 1,
                                              'TR_C_gene' => 1,
                                              'TR_J_gene' => 1,
                                              'TR_V_gene' => 1,
                                              'TR_D_gene' => 1,
                                             },
                         #the biotypes of the best genes always to be kept
                         MAX_EXON_LENGTH => 20000,
                         #if the coding_only flag is set to 1, the transcript clustering into genes is done over coding exons only
                         # the current standard way is to cluster only on coding exons
                         CODING_ONLY => 1,
                       },
        -rc_name    => 'transcript_finalisation',
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
      },


      {
        -logic_name => 'layer_annotation_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         skip_rnaseq => $self->o('skip_rnaseq'),
                         skip_projection => $self->o('skip_projection'),
                         target_db => $self->o('layering_db'),
                         sanity_check_type => 'gene_db_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'gene_db_checks')->{$self->default_options->{'uniprot_set'}}->{'layer'},
                       },

        -rc_name    => '4GB',
        -flow_into => {
                        1 => ['genebuilder_sanity_checks'],
                      },
      },


      {
        -logic_name => 'genebuilder_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         skip_rnaseq => $self->o('skip_rnaseq'),
                         skip_projection => $self->o('skip_projection'),
                         target_db => $self->o('genebuilder_db'),
                         sanity_check_type => 'gene_db_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'gene_db_checks')->{$self->default_options->{'uniprot_set'}}->{'genebuilder'},
                       },

        -rc_name    => '4GB',
        -flow_into => {
                        1 => ['restore_ig_tr_biotypes'],
                      },
      },


      {
        -logic_name => 'restore_ig_tr_biotypes',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('genebuilder_db'),
          sql => [
            'UPDATE gene JOIN transcript USING(gene_id) SET gene.biotype = transcript.biotype'.
              ' WHERE transcript.biotype LIKE "IG\_%" OR transcript.biotype LIKE "TR\_%"',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['create_pseudogene_db'],
                      },
      },


      {
        -logic_name => 'create_pseudogene_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('pseudogene_db'),
                         create_type => 'clone',
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['pseudogenes'],
                      },
      },


      {
        -logic_name => 'pseudogenes',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivePseudogenes',
        -parameters => {
                         single_multi_file => 1,
                         output_path => $self->o('output_path').'/pseudogenes/',
                         input_gene_db => $self->o('genebuilder_db'),
                         repeat_db => $self->default_options->{'dna_db'},
                         output_db => $self->o('pseudogene_db'),
                         dna_db => $self->default_options->{'dna_db'},
                         logic_name => 'pseudogenes',
                         module     => 'HivePseudogenes',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::PseudoGeneStatic','pseudogenes')},
                       },

	     -rc_name    => 'default_himem',
	     -flow_into => {
			    1 => ['format_blast_db'],
                      },
      },


      {
        -logic_name => 'format_blast_db',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'if [ "'.$self->o('blast_type').
                                '" = "ncbi" ];then makeblastdb -dbtype nucl -in '.
			        $self->o('output_path').'/pseudogenes/all_multi_exon_genes.fasta;'.
	                        ' else xdformat -n '.$self->o('output_path').'/pseudogenes/all_multi_exon_genes.fasta;fi'
		       },
         -rc_name => 'default',
         -flow_into => {
                         1 => ['spliced_elsewhere'],
                       },
      },


      {
        -logic_name => 'spliced_elsewhere',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSplicedElsewhere',
        -parameters => {
                         multi_exon_db_path => $self->o('output_path').'/pseudogenes/',
                         input_gene_db => $self->o('genebuilder_db'),
                         repeat_db => $self->default_options->{'dna_db'},
                         output_db => $self->o('pseudogene_db'),
                         dna_db => $self->default_options->{'dna_db'},
                         logic_name => 'spliced_elsewhere',
                         module     => 'HiveSplicedElsewhere',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::PseudoGeneStatic','pseudogenes')},
                       },
        -rc_name          => 'default_himem',
        -flow_into => {
                        1 => ['create_final_geneset_db'],
                      },
      },


      {
        -logic_name => 'create_final_geneset_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('pseudogene_db'),
                         target_db => $self->o('final_geneset_db'),
                         create_type => 'copy',
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '1' => ['update_rnaseq_ise_logic_names'],
                      },
      },


      {
        -logic_name => 'update_rnaseq_ise_logic_names',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('final_geneset_db'),
          sql => [
            'UPDATE analysis SET logic_name = CONCAT(logic_name, "_ise") WHERE logic_name LIKE "%\_rnaseq"',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['run_cleaner'],
                      },
      },


      {
        -logic_name => 'run_cleaner',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCleanGeneset',
        -parameters => {
                         skip_analysis => $self->o('skip_cleaning'),
                         input_db => $self->o('final_geneset_db'),
                         dna_db => $self->o('dna_db'),
                         output_path => $self->o('output_path').'/clean_genes/',
                         blessed_biotypes => $self->o('cleaning_blessed_biotypes'),
                         flagged_redundancy_coverage_threshold => 95,
                         general_redundancy_coverage_threshold => 95,
                       },
        -rc_name    => '4GB',
        -flow_into => {
                        '1' => ['delete_flagged_transcripts'],
                      },
      },


      {
       -logic_name => 'delete_flagged_transcripts',
       -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDeleteTranscripts',
       -parameters => {
                        skip_analysis => $self->o('skip_cleaning'),
                        dbhost => $self->o('final_geneset_db','-host'),
                        dbname => $self->o('final_geneset_db','-dbname'),
                        dbuser => $self->o('user'),
                        dbpass => $self->o('password'),
                        dbport => $self->o('final_geneset_db','-port'),
                        transcript_ids_file => catfile($self->o('output_path'), 'clean_genes', 'transcript_ids_to_remove.txt'),
                        delete_transcripts_path => catdir($self->o('ensembl_analysis_script'), 'genebuild'),
                        delete_genes_path => catdir($self->o('enscode_root_dir'), 'genebuild'),
                        delete_transcripts_script_name => 'delete_transcripts.pl',
                        delete_genes_script_name => 'delete_genes.pl',
                        output_path => catdir($self->o('output_path'), 'clean_genes'),
                        output_file_name => 'delete_transcripts.out',
                      },
        -max_retry_count => 0,
        -flow_into => {
          '1->A' => ['fan_ncrna'],
          'A->1' => ['dummy_projected_lincrna_pseudos'],
        },
     },

############################################################################
#
# ncRNA pipeline
#
############################################################################
      {
        -logic_name => 'fan_ncrna',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'if [ "#skip_ncrna#" -ne "0" ]; then exit 42; else exit 0;fi',
                         return_codes_2_branches => {'42' => 2},
                       },
        -rc_name => 'default',
        -flow_into  => {
          1 => ['create_ncrna_db'],
        },
      },


      {
        -logic_name => 'create_ncrna_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('ncrna_db'),
                         create_type => 'clone',
                       },
        -rc_name    => 'default',
        -flow_into => {
                          '1->A' => ['create_small_rna_slice_ids'],
                          'A->1' => ['ncrna_sanity_checks'],
                      },

      },


      {
        -logic_name => 'create_small_rna_slice_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         coord_system_name => 'toplevel',
                         iid_type => 'slice',
                         slice_size => 1000000,
                         top_level => 1,
                         target_db => $self->o('dna_db'),
                         batch_slice_ids => 1,
                         batch_target_size => 2000000,
                      },
        -flow_into => {
                       '2->A' => ['mirna_blast','rfam_blast'],
                       'A->1' => ['filter_ncrnas'],
                      },
        -rc_name    => 'default',
      },


      {
        -logic_name => 'rfam_blast',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastRfam',
        -parameters => {
                         repeat_logic_names => ['dust'],
                         repeat_db => $self->o('dna_db'),
                         output_db => $self->o('ncrna_db'),
                         dna_db => $self->o('dna_db'),
                         logic_name => 'rfamblast',
                         module     => 'HiveBlastRfam',
                         blast_db_path => $self->o('ncrna_blast_path')."/filtered.fasta",
                         blast_exe_path => $self->o('blastn_exe_path'),
                         commandline_params => ' -num_threads 3 -word_size 12 -num_alignments 5000  -num_descriptions 5000 -max_hsps 1 ',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastRFam', {BLAST_PARAMS => {type => $self->o('blast_type')}})},
                         timer => '3h',
                       },
       -flow_into => {
                       '-1' => ['rebatch_rfam'],
                       '-2' => ['rebatch_rfam'],
                     },
       -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
       -rc_name    => 'rfam_blast',
      },


      {
        -logic_name => 'rebatch_rfam',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('dna_db'),
                         iid_type => 'rebatch_and_resize_slices',
                         slice_size => 100000,
                         batch_target_size => 100000,
                       },
       -flow_into => {
                       '2' => ['rfam_blast_retry'],
                     },
       -rc_name    => 'default',
       -can_be_empty  => 1,
      },


      {
        -logic_name => 'rfam_blast_retry',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastRfam',
        -parameters => {
                         repeat_logic_names => ['dust'],
                         repeat_db => $self->o('dna_db'),
                         output_db => $self->o('ncrna_db'),
                         dna_db => $self->o('dna_db'),
                         logic_name => 'rfamblast',
                         module     => 'HiveBlastRfam',
                         blast_db_path => $self->o('ncrna_blast_path')."/filtered.fasta",
                         blast_exe_path => $self->o('blastn_exe_path'),
                         commandline_params => ' -num_threads 3 -word_size 12 -num_alignments 5000  -num_descriptions 5000 -max_hsps 1 ',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastRFam', {BLAST_PARAMS => {type => $self->o('blast_type')}})},
                         timer => '1h',
                       },
       -flow_into => {
                       '-1' => ['failed_rfam_blast_job'],
                       '-2' => ['failed_rfam_blast_job'],
                     },
       -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
       -rc_name    => 'rfam_blast_retry',
       -can_be_empty  => 1,
      },


      {
        -logic_name => 'failed_rfam_blast_job',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name          => 'default',
        -can_be_empty  => 1,
      },


      {
        -logic_name => 'mirna_blast',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastmiRNA',
        -parameters => {
                         repeat_logic_names => ['dust'],
                         repeat_db => $self->o('dna_db'),
                         output_db => $self->o('ncrna_db'),
                         dna_db => $self->o('dna_db'),
                         logic_name => 'blastmirna',
                         module     => 'HiveBlastmiRNA',
                         blast_db_path => $self->o('mirna_blast_path') . '/' . $self->o('mirBase_fasta') ,
                         blast_exe_path => $self->o('blastn_exe_path'),
                         commandline_params => ' -num_threads 3 ',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastmiRBase', {BLAST_PARAMS => {type => $self->o('blast_type')}})},
                         timer => '2h',
                       },

        -flow_into => {
                        '-1' => ['rebatch_mirna'],
                        '-2' => ['rebatch_mirna'],
                      },
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -rc_name    => 'blast',
      },


      {
        -logic_name => 'rebatch_mirna',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('dna_db'),
                         iid_type => 'rebatch_and_resize_slices',
                         slice_size => 100000,
                         batch_target_size => 100000,
                       },
       -flow_into => {
                       '2' => ['mirna_blast_retry'],
                     },
       -rc_name    => 'default',
       -can_be_empty  => 1,
      },


      {
        -logic_name => 'mirna_blast_retry',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastmiRNA',
        -parameters => {
                         repeat_logic_names => ['dust'],
                         repeat_db => $self->o('dna_db'),
                         output_db => $self->o('ncrna_db'),
                         dna_db => $self->o('dna_db'),
                         logic_name => 'blastmirna',
                         module     => 'HiveBlastmiRNA',
                         blast_db_path => $self->o('mirna_blast_pathh') . '/' . $self->o('mirBase_fasta') ,
                         blast_exe_path => $self->o('blastn_exe_path'),
                         commandline_params => ' -num_threads 3 ',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastmiRBase', {BLAST_PARAMS => {type => $self->o('blast_type')}})},
                         timer => '1h',
                       },

        -flow_into => {
                        '-1' => ['failed_mirna_blast_job'],
                        '-2' => ['failed_mirna_blast_job'],
                      },
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -rc_name    => 'blast_retry',
      },


      {
        -logic_name => 'failed_mirna_blast_job',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name          => 'default',
        -can_be_empty  => 1,
      },


      {
        -logic_name => 'filter_ncrnas',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveFilterncRNAs',
        -parameters => {
                         output_db => $self->o('ncrna_db'),
                         dna_db => $self->o('dna_db'),
                         logic_name => 'filter_ncrnas',
                         module     => 'HiveFilterncRNAs',
                       },
        -rc_name    => 'filter',
        -flow_into => {
                        '2' => ['run_mirna'],
                        '3' => ['run_infernal'],
                      },
      },


      {
        -logic_name => 'run_mirna',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivemiRNA',
        -parameters => {
                         output_db => $self->o('ncrna_db'),
                         dna_db => $self->o('dna_db'),
                         logic_name => 'ncrna',
                         module     => 'HivemiRNA',
                         blast_db_dir_path => $self->o('mirna_blast_path').'/all_mirnas.embl',
                         output_dir => $self->o('output_path'),
                       },
        -batch_size => 20,
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -rc_name    => 'filter',
        -flow_into  => { '1->A' => ['dump_features', 'dump_genome'], 'A->1' => ['filter_mirnas']},
      },

      {
        -logic_name => 'dump_features',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => "perl " . $self->o('mirna_analysis_script') . "/dump_prefilter_features.pl " .
                                $self->o('dbowner').'_'.$self->o('production_name').'_core_'.$self->o('release_number') . " " .
                                $self->o('dna_db_server') . " " .
                                $self->o('dna_db_port') . " " .
                                $self->o('user_r') . " " .
                                $self->o('output_path') . " blastmirna",
                        },
         -rc_name   => 'filter',
      },
      {
        -logic_name => 'dump_genome',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                          cmd => "perl " . $self->o('mirna_analysis_script') . "/sequence_dump.pl -dbhost " .
                                  $self->o('dna_db_server') . " -dbname  " .
                                  $self->o('dbowner').'_'.$self->o('production_name').'_core_'.$self->o('release_number') . "-dbport " .
                                  $self->o('dna_db_port') . " -dbuser " .
                                  $self->o('user_r') . " -coord_system_name toplevel -mask -mask_repeat RepeatMask -output_dir " .
                                  $self->o('output_path') . " -filename genome.fasta -softmask -onefile -header rnaseq",
                      },
         -rc_name   => 'filter',           
                                  
      },

      {
        -logic_name => 'filter_mirnas',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                          cmd => "sh " . $self->o('mirna_analysis_script') . "/FilterMiRNAs.sh -d " .
                                  $self->o('output_path') . "/blastmirna_dafs.bed -r " .
                                  $self->o('output_path') . "/repeats.bed -g " .
                                  $self->o('output_path') . "/genome.fasta -w " .
                                  $self->o('output_path') . " -m " .
                                  $self->o('mirna_blast_path') . "/rfc_filters/" . $self->o('rfc_model') . " -s " .
                                  $self->o('mirna_blast_path') . "/rfc_filters/" . $self->o('rfc_scaler') . " -c " .
                                  $self->o('ncrna_db_server') . ":" . $self->o('ncrna_db_port') . ":" . 
                                  $self->o('dbowner').'_'.$self->o('production_name').'_ncrna_'.$self->o('release_number') . ":" .
                                  $self->o('user') . ":" . $self->o('password'),
                        },
        -rc_name   => 'filter',
        
      },


      {
        -logic_name => 'run_infernal',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveInfernal',
        -parameters => {
                         output_db => $self->o('ncrna_db'),
                         dna_db => $self->o('dna_db'),
                         logic_name => 'ncrna',
                         module     => 'HiveInfernal',
                         cmsearch_exe_path => $self->o('cmsearch_exe_path'),
                         blast_db_dir_path => $self->o('ncrna_blast_path').'/',
                       },
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -rc_name    => 'transcript_finalisation',
      },


      {
        -logic_name => 'ncrna_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         target_db => $self->o('ncrna_db'),
                         sanity_check_type => 'gene_db_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'gene_db_checks')->{$self->default_options->{'uniprot_set'}}->{'ncrna'},
                       },

        -rc_name    => '4GB',
        -flow_into => {
          1 => ['transfer_ncrnas'],
        },
      },


     {
       -logic_name => 'transfer_ncrnas',
       -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
       -parameters => {
                        cmd => 'perl '.catfile($self->o('ensembl_analysis_script'), 'genebuild', 'copy_genes.pl').
                                     ' -sourcehost '.$self->o('ncrna_db','-host').
                                     ' -sourceuser '.$self->o('user_r').
                                     ' -sourceport '.$self->o('ncrna_db','-port').
                                     ' -sourcedbname '.$self->o('ncrna_db','-dbname').
                                     ' -dnauser '.$self->o('user_r').
                                     ' -dnahost '.$self->o('dna_db','-host').
                                     ' -dnaport '.$self->o('dna_db','-port').
                                     ' -dnadbname '.$self->o('dna_db','-dbname').
                                     ' -targetuser '.$self->o('user').
                                     ' -targetpass '.$self->o('password').
                                     ' -targethost '.$self->o('final_geneset_db','-host').
                                     ' -targetport '.$self->o('final_geneset_db','-port').
                                     ' -targetdbname '.$self->o('final_geneset_db','-dbname').
                                     ' -all'
                      },
        -rc_name => 'default',
     },

############################################################################
#
# Transfer projected lincRNA and pseudogenes
#
############################################################################
     {
       -logic_name => 'dummy_projected_lincrna_pseudos',
       -module => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
       -parameters => {},
       -rc_name => 'default',
       -flow_into  => {
         '1->A' => ['fan_transfer_projected_genes'],
         'A->1' => ['delete_duplicate_genes'],
       },
     },


     {
       -logic_name => 'fan_transfer_projected_genes',
       -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
       -parameters => {
         cmd => 'if [ "#skip_projection#" -ne -1 ]; then exit 42; else exit 0;fi', # To enable this part, replace -1 with 0
         return_codes_2_branches => {'42' => 2},
       },
       -rc_name => 'default',
       -flow_into  => {
         1 => ['create_projected_lincrna_ids_to_copy'],
       },
     },


     {
       -logic_name => 'create_projected_lincrna_ids_to_copy',
       -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
       -parameters => {
                        target_db    => $self->o('projection_lincrna_db'),
                        iid_type     => 'feature_id',
                        feature_type => 'gene',
                        batch_size   => 500,
                     },
       -flow_into => {
         '2->A' => ['transfer_projected_lincrnas'],
         'A->1' => ['create_projected_pseudogene_ids_to_copy'],
       },
       -rc_name    => 'default',
     },


     {
       -logic_name => 'transfer_projected_lincrnas',
       -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCopyGenes',
       -parameters => {
         copy_genes_directly => 1,
         source_db => $self->o('projection_lincrna_db'),
         dna_db => $self->o('dna_db'),
         target_db => $self->o('final_geneset_db'),
         filter_on_overlap => 1,
         filter_on_strand  => 0,
         overlap_filter_type => 'genomic_overlap',
         filter_against_biotypes => {
           protein_coding => 1,
           IG_C_gene => 1,
           IG_V_gene => 1,
           TR_C_gene => 1,
           TR_J_gene => 1,
           TR_V_gene => 1,
           processed_pseudogene => 1,
           pseudogene => 1,
         },
       },
       -rc_name    => 'default',
     },


     {
       -logic_name => 'create_projected_pseudogene_ids_to_copy',
       -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
       -parameters => {
                        target_db    => $self->o('projection_pseudogene_db'),
                        iid_type     => 'feature_id',
                        feature_type => 'gene',
                        batch_size   => 500,
                      },
       -flow_into => {
                       '2' => ['transfer_projected_pseudogenes'],
                     },
       -rc_name    => 'default',
     },


      {
        -logic_name => 'transfer_projected_pseudogenes',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCopyGenes',
        -parameters => {
                         copy_genes_directly => 1,
                         source_db => $self->o('projection_pseudogene_db'),
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('final_geneset_db'),
                         filter_on_overlap => 1,
                         filter_on_strand  => 1,
                         overlap_filter_type => 'genomic_overlap',
                         filter_against_biotypes => {
                                                      protein_coding => 1,
                                                      IG_C_gene => 1,
                                                      IG_V_gene => 1,
                                                      TR_C_gene => 1,
                                                      TR_J_gene => 1,
                                                      TR_V_gene => 1,
                                                      processed_pseudogene => 1,
                                                      pseudogene => 1,
                                                      lincRNA => 1,
                                                    },
                       },
        -rc_name    => 'default',
      },


     {
       -logic_name => 'delete_duplicate_genes',
       -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
       -parameters => {
                        cmd => "perl ".$self->o('remove_duplicates_script_path')
                               ." -dbhost ".$self->o('final_geneset_db','-host')
                               ." -dbuser ".$self->o('user')
                               ." -dbpass ".$self->o('password')
                               ." -dbname ".$self->o('final_geneset_db','-dbname')
                               ." -dbport ".$self->o('final_geneset_db','-port')
                               ." -dnadbhost ".$self->o('dna_db','-host')
                               ." -dnadbuser ".$self->o('user_r')
                               ." -dnadbname ".$self->o('dna_db','-dbname')
                               ." -dnadbport ".$self->o('dna_db','-port'),
                     },
        -max_retry_count => 0,
        -rc_name => 'default',
        -flow_into => {
                        '1' => ['final_db_sanity_checks'],
                      },
      },


      {
        -logic_name => 'final_db_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         skip_rnaseq => $self->o('skip_rnaseq'),
                         skip_projection => $self->o('skip_projection'),
                         target_db => $self->o('final_geneset_db'),
                         sanity_check_type => 'gene_db_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'gene_db_checks')->{$self->default_options->{'uniprot_set'}}->{'final'},
                       },

        -rc_name    => '4GB',
        -flow_into => {
                        '1->A' => ['create_gene_ids_to_copy'],
                        'A->1' => ['update_biotypes_and_analyses'],
                      },
      },


      {
        -logic_name => 'create_gene_ids_to_copy',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db    => $self->o('final_geneset_db'),
                         iid_type     => 'feature_id',
                         feature_type => 'gene',
                         batch_size   => 500,
                      },
        -flow_into => {
                       '2' => ['copy_genes_to_core'],
                      },

        -rc_name    => 'default',
      },


      {
        -logic_name => 'copy_genes_to_core',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCopyGenes',
        -parameters => {
                         copy_genes_directly => 1,
                         source_db => $self->o('final_geneset_db'),
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('reference_db'),
                       },
        -rc_name    => 'default',
      },


      {
        -logic_name => 'update_biotypes_and_analyses',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('reference_db'),
          sql => [
            'UPDATE gene SET biotype = "protein_coding" WHERE biotype = "ensembl"',
            'UPDATE gene SET analysis_id = (SELECT analysis_id FROM analysis WHERE logic_name = "ensembl")'.
              ' WHERE analysis_id IN'.
              ' (SELECT analysis_id FROM analysis'.
              ' WHERE logic_name IN ("spliced_elsewhere","pseudogenes","genblast","project_transcripts","project_pseudogene","project_lincrna"))',
            'UPDATE transcript JOIN gene USING(gene_id) SET transcript.biotype = gene.biotype',
            'UPDATE transcript JOIN gene USING(gene_id) SET transcript.analysis_id = gene.analysis_id',
            'INSERT IGNORE into analysis (created,logic_name,db) VALUES (NOW(),"other_protein","uniprot")',
            'INSERT IGNORE into analysis (created,logic_name,db) VALUES (NOW(),"projected_transcript","'.$self->o('projection_source_db_name').'")',
            'UPDATE protein_align_feature SET analysis_id ='.
              '(SELECT analysis_id FROM analysis WHERE logic_name = "projected_transcript") WHERE analysis_id ='.
              '(SELECT analysis_id FROM analysis WHERE logic_name = "project_transcripts")',
            'UPDATE protein_align_feature SET analysis_id ='.
              '(SELECT analysis_id FROM analysis WHERE logic_name = "other_protein") WHERE analysis_id NOT IN'.
              '(SELECT analysis_id FROM analysis WHERE logic_name IN ("uniprot","projected_transcript"))',
            'UPDATE dna_align_feature SET analysis_id ='.
              '(SELECT analysis_id FROM analysis WHERE logic_name = "projected_transcript") WHERE analysis_id IN'.
              '(SELECT analysis_id FROM analysis WHERE logic_name IN ("project_lincrna","project_pseudogene"))',
            'UPDATE repeat_feature SET seq_region_start = 1 WHERE seq_region_start < 1',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['set_meta_coords'],
                      },
      },


      {
        -logic_name => 'set_meta_coords',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('meta_coord_script').
                                ' -user '.$self->o('user').
                                ' -pass '.$self->o('password').
                                ' -host '.$self->o('reference_db','-host').
                                ' -port '.$self->o('reference_db','-port').
                                ' -dbpattern '.$self->o('reference_db','-dbname')
                       },
        -rc_name => 'default',
        -flow_into => {
                        1 => ['set_meta_levels'],
                      },
      },


      {
        -logic_name => 'set_meta_levels',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('meta_levels_script').
                                ' -user '.$self->o('user').
                                ' -pass '.$self->o('password').
                                ' -host '.$self->o('reference_db','-host').
                                ' -port '.$self->o('reference_db','-port').
                                ' -dbname '.$self->o('reference_db','-dbname')
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['set_frameshift_introns'] },
      },

      {
        -logic_name => 'set_frameshift_introns',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('frameshift_attrib_script').
                                ' -user '.$self->o('user').
                                ' -pass '.$self->o('password').
                                ' -host '.$self->o('reference_db','-host').
                                ' -port '.$self->o('reference_db','-port').
                                ' -dbpattern '.$self->o('reference_db','-dbname')
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['set_canonical_transcripts'] },
      },

      {
        -logic_name => 'set_canonical_transcripts',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('select_canonical_script').
                                ' -dbuser '.$self->o('user').
                                ' -dbpass '.$self->o('password').
                                ' -dbhost '.$self->o('reference_db','-host').
                                ' -dbport '.$self->o('reference_db','-port').
                                ' -dbname '.$self->o('reference_db','-dbname').
                                ' -coord toplevel -write'
                       },
        -rc_name => '1.5GB',
        -flow_into => { 1 => ['null_columns'] },
      },


      {
        -logic_name => 'null_columns',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('reference_db'),
          sql => [
            'UPDATE transcript set stable_id = NULL',
            'UPDATE translation set stable_id = NULL',
            'UPDATE exon set stable_id = NULL',
            'UPDATE protein_align_feature set external_db_id = NULL',
            'UPDATE dna_align_feature set external_db_id = NULL',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['run_stable_ids'],
                      },
      },


      {
        -logic_name => 'run_stable_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::SetStableIDs',
        -parameters => {
                         enscode_root_dir => $self->o('enscode_root_dir'),
                         mapping_required => $self->o('mapping_required'),
                         target_db => $self->o('reference_db'),
                         mapping_db => $self->o('mapping_db'),
                         id_start => $self->o('stable_id_prefix').$self->o('stable_id_start'),
                         output_path => $self->o('output_path'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['backup_core_db_pre_optimise'],
                      },
      },


      {
        # Creates a reference db for each species
        -logic_name => 'backup_core_db_pre_optimise',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::DatabaseDumper',
        -parameters => {
                         src_db_conn => $self->o('dna_db'),
                         output_file => catfile($self->o('output_path'), 'core_post_stable_idsbak.sql.gz'),
                       },
        -rc_name    => 'default',
        -flow_into => { 1 => ['load_external_db_ids_and_optimise_af_tables'] },
      },


      {
        -logic_name => 'load_external_db_ids_and_optimise_af_tables',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('load_optimise_script').
                                ' -output_path '.catdir($self->o('output_path'), 'optimise').
                                ' -uniprot_filename '.$self->o('uniprot_entry_loc').
                                ' -dbuser '.$self->o('user').
                                ' -dbpass '.$self->o('password').
                                ' -dbport '.$self->o('reference_db','-port').
                                ' -dbhost '.$self->o('reference_db','-host').
                                ' -dbname '.$self->o('reference_db','-dbname').
                                ' -prod_dbuser '.$self->o('user_r').
                                ' -prod_dbhost '.$self->o('production_db','-host').
                                ' -prod_dbname '.$self->o('production_db','-dbname').
                                ' -prod_dbport '.$self->o('production_db','-port').
                                ' -ise -core'
                       },
        -max_retry_count => 1,
        -rc_name => '8GB',
        -flow_into => {
                        1 => ['clean_unused_analyses'],
                      },
      },


      {
        -logic_name => 'clean_unused_analyses',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('reference_db'),
          sql => [
            'DELETE FROM analysis WHERE logic_name IN'.
              ' ("spliced_elsewhere","pseudogenes","genblast","project_pseudogene",'.
              ' "project_lincrna","project_transcripts","ig_tr_collapse")',
            'DELETE FROM ad USING analysis_description ad LEFT JOIN analysis a ON ad.analysis_id = a.analysis_id'.
              ' WHERE a.analysis_id IS NULL',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['drop_backup_tables_job'],
                      },
      },


      {
        -logic_name => 'drop_backup_tables_job',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
        -parameters => {
          inputquery => 'SHOW TABLES LIKE "%bak%"',
          column_names => ['table'],
          db_conn =>$self->o('reference_db'),
        },
        -rc_name    => 'default',
        -flow_into => {
          '2->A' => ['drop_backup_tables'],
          'A->1' => ['core_gene_set_sanity_checks'],
        },
      },


      {
        -logic_name => 'drop_backup_tables',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          sql => 'DROP TABLE #table#',
          db_conn =>$self->o('reference_db'),
        },
        -rc_name    => 'default',
      },


      {
        -logic_name => 'core_gene_set_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         skip_rnaseq => $self->o('skip_rnaseq'),
                         skip_projection => $self->o('skip_projection'),
                         target_db => $self->o('reference_db'),
                         sanity_check_type => 'gene_db_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'gene_db_checks')->{$self->default_options->{'uniprot_set'}}->{'core'},
                       },

        -rc_name    => '4GB',
        -flow_into => {
                        1 => ['core_healthchecks'],
                      },
      },


      {
        -logic_name => 'core_healthchecks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveHealthcheck',
        -parameters => {
                         input_db         => $self->o('reference_db'),
                         species          => $self->o('species_name'),
                         group            => 'core_handover',
                         enscode_root_dir => $self->o('enscode_root_dir'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['fan_otherfeatures_db','fan_rnaseq_final_db'],
                      },
      },


      {
        -logic_name => 'fan_otherfeatures_db',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'if [ -n "#assembly_refseq_accession#" ]; then exit 0; else exit 42;fi',
                         return_codes_2_branches => {'42' => 2},
                       },
        -rc_name => 'default',
        -flow_into  => {
                          1 => ['create_otherfeatures_db'],
                       },
      },


      {
        -logic_name => 'create_otherfeatures_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('cdna_db'),
                         target_db => $self->o('otherfeatures_db'),
                         create_type => 'copy',
                       },
        -rc_name    => 'default',

        -flow_into => {
                        '1->A' => ['create_refseq_import_ids_to_copy'],
                        'A->1' => ['update_otherfeatures_db'],
                      },
      },


      {
        -logic_name => 'create_refseq_import_ids_to_copy',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db    => $self->o('refseq_db'),
                         iid_type     => 'feature_id',
                         feature_type => 'gene',
                         batch_size   => 500,
                      },
        -flow_into => {
                       '2' => ['copy_refseq_genes_to_otherfeatures'],
                      },

        -rc_name    => 'default',
      },


      {
        -logic_name => 'copy_refseq_genes_to_otherfeatures',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCopyGenes',
        -parameters => {
                         copy_genes_directly => 1,
                         source_db => $self->o('refseq_db'),
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('otherfeatures_db'),
                       },
        -rc_name    => 'default',
      },


      {
        -logic_name => 'update_otherfeatures_db',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('otherfeatures_db'),
          sql => [
            'DELETE analysis_description FROM analysis_description join analysis using(analysis_id)'.
              ' WHERE logic_name NOT IN ("refseq_import","cdna_alignment")',
            'DELETE FROM analysis WHERE logic_name NOT IN ("refseq_import","cdna_alignment")',
            'DELETE FROM meta WHERE meta_key LIKE "%.level"',
            'DELETE FROM meta WHERE meta_key LIKE "provider.%"',
            'DELETE FROM meta WHERE meta_key LIKE "assembly.web_accession%"',
            'DELETE FROM meta WHERE meta_key LIKE "removed_evidence_flag.%"',
            'DELETE FROM meta WHERE meta_key LIKE "marker.%"',
            'DELETE FROM meta WHERE meta_key = "repeat.analysis"',
            'DELETE FROM meta WHERE meta_key IN'.
              ' ("genebuild.last_geneset_update","genebuild.method","genebuild.projection_source_db","genebuild.start_date")',
            'INSERT INTO meta (species_id,meta_key,meta_value) VALUES (1,"genebuild.last_otherfeatures_update",NOW())',
            'UPDATE transcript JOIN transcript_supporting_feature USING(transcript_id)'.
              ' JOIN dna_align_feature ON feature_id = dna_align_feature_id SET stable_id = hit_name',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['set_otherfeatures_meta_coords'],
                      },
      },


      {
        -logic_name => 'set_otherfeatures_meta_coords',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('meta_coord_script').
                                ' -user '.$self->o('user').
                                ' -pass '.$self->o('password').
                                ' -host '.$self->o('otherfeatures_db','-host').
                                ' -port '.$self->o('otherfeatures_db','-port').
                                ' -dbpattern '.$self->o('otherfeatures_db','-dbname')
                       },
        -rc_name => 'default',
        -flow_into => {
                        1 => ['set_otherfeatures_meta_levels'],
                      },
      },


      {
        -logic_name => 'set_otherfeatures_meta_levels',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('meta_levels_script').
                                ' -user '.$self->o('user').
                                ' -pass '.$self->o('password').
                                ' -host '.$self->o('otherfeatures_db','-host').
                                ' -port '.$self->o('otherfeatures_db','-port').
                                ' -dbname '.$self->o('otherfeatures_db','-dbname')
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['set_otherfeatures_frameshift_introns'] },
      },


      {
        -logic_name => 'set_otherfeatures_frameshift_introns',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('frameshift_attrib_script').
                                ' -user '.$self->o('user').
                                ' -pass '.$self->o('password').
                                ' -host '.$self->o('otherfeatures_db','-host').
                                ' -port '.$self->o('otherfeatures_db','-port').
                                ' -dbpattern '.$self->o('otherfeatures_db','-dbname')
                       },
        -rc_name => '4GB',
        -flow_into => { 1 => ['set_otherfeatures_canonical_transcripts'] },
      },


      {
        -logic_name => 'set_otherfeatures_canonical_transcripts',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('select_canonical_script').
                                ' -dbuser '.$self->o('user').
                                ' -dbpass '.$self->o('password').
                                ' -dbhost '.$self->o('otherfeatures_db','-host').
                                ' -dbport '.$self->o('otherfeatures_db','-port').
                                ' -dbname '.$self->o('otherfeatures_db','-dbname').
                                ' -dnadbuser '.$self->o('user_r').
                                ' -dnadbhost '.$self->o('dna_db','-host').
                                ' -dnadbport '.$self->o('dna_db','-port').
                                ' -dnadbname '.$self->o('dna_db','-dbname').
                                ' -coord toplevel -write'
                       },
        -rc_name => '1.5GB',
        -flow_into => { 1 => ['populate_production_tables_otherfeatures'] },
      },


      {
        -logic_name => 'populate_production_tables_otherfeatures',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HivePopulateProductionTables',
        -parameters => {
                         'target_db'        => $self->o('otherfeatures_db'),
                         'output_path'      => $self->o('output_path'),
                         'enscode_root_dir' => $self->o('enscode_root_dir'),
                         'production_db'    => $self->o('production_db'),
                       },
        -rc_name    => 'default',
        -flow_into  => {
                         1 => ['null_otherfeatures_columns'],
                       },
      },


      {
        -logic_name => 'null_otherfeatures_columns',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
        -parameters => {
          db_conn => $self->o('otherfeatures_db'),
          sql => [
            'UPDATE dna_align_feature SET external_db_id = NULL',
          ],
        },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['load_external_db_ids_and_optimise_otherfeatures'],
                      },
      },


      {
        -logic_name => 'load_external_db_ids_and_optimise_otherfeatures',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('load_optimise_script').
                                ' -output_path '.$self->o('output_path').'/optimise_otherfeatures/'.
                                ' -uniprot_filename '.$self->o('uniprot_entry_loc').
                                ' -dbuser '.$self->o('user').
                                ' -dbpass '.$self->o('password').
                                ' -dbport '.$self->o('otherfeatures_db','-port').
                                ' -dbhost '.$self->o('otherfeatures_db','-host').
                                ' -dbname '.$self->o('otherfeatures_db','-dbname').
                                ' -prod_dbuser '.$self->o('user_r').
                                ' -prod_dbhost '.$self->o('production_db','-host').
                                ' -prod_dbname '.$self->o('production_db','-dbname').
                                ' -prod_dbport '.$self->o('production_db','-port').
                                ' -verbose'
                       },
        -max_retry_count => 1,
        -rc_name => '4GB',
        -flow_into => {
                        1 => ['otherfeatures_sanity_checks'],
                      },
      },


      {
        -logic_name => 'otherfeatures_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         target_db => $self->o('otherfeatures_db'),
                         sanity_check_type => 'gene_db_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'gene_db_checks')->{'otherfeatures'},
                       },

        -rc_name    => '4GB',
        -flow_into => {
                        1 => ['otherfeatures_healthchecks'],
                      },
      },


      {
        -logic_name => 'otherfeatures_healthchecks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveHealthcheck',
        -parameters => {
			 input_db         => $self->o('reference_db'),
			 species          => $self->o('species_name'),
                         group            => 'otherfeatures_handover',
                         enscode_root_dir => $self->o('enscode_root_dir'),
                       },

        -rc_name    => '4GB',
      },


      {
        -logic_name => 'fan_rnaseq_final_db',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'if [ "#skip_rnaseq#" -ne "0" ]; then exit 42; else exit 0;fi',
                         return_codes_2_branches => {'42' => 2},
                       },
        -rc_name => 'default',
        -flow_into  => {
                          1 => ['create_rnaseq_final_db'],
                       },

      },


      {
        -logic_name => 'create_rnaseq_final_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('rnaseq_blast_db'),
                         target_db => $self->o('rnaseq_final_db'),
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
                         db_conn => $self->o('rnaseq_final_db'),
                         sql => [
                          'TRUNCATE dna_align_feature',
                          'DELETE FROM transcript_supporting_feature WHERE feature_type = "dna_align_feature"',
                          'DELETE FROM supporting_feature WHERE feature_type = "dna_align_feature"',
                          'INSERT IGNORE INTO analysis (logic_name, module, created, db_version) VALUES ("other_protein", "HiveBlastRNAseq", NOW(), "#uniprot_version#")',
                          'UPDATE protein_align_feature paf, analysis a SET paf.analysis_id = a.analysis_id WHERE a.logic_name = "other_protein"',
                          'DELETE FROM meta WHERE meta_key LIKE "provider\.%"',
                          'DELETE FROM meta WHERE meta_key LIKE "assembly.web_accession%"',
                          'DELETE FROM meta WHERE meta_key LIKE "removed_evidence_flag\.%"',
                          'DELETE FROM meta WHERE meta_key LIKE "marker\.%"',
                          'DELETE FROM meta WHERE meta_key IN ("genebuild.method","genebuild.projection_source_db","genebuild.start_date","repeat.analysis")',
                          'UPDATE gene JOIN transcript USING(gene_id) SET canonical_transcript_id = transcript_id',
                          'UPDATE transcript JOIN translation USING(transcript_id) SET canonical_translation_id = translation_id',
                         ],
                         uniprot_version => $self->o('uniprot_db_dir'),
                       },
        -rc_name    => 'default',

        -flow_into => {
                        '1' => ['generate_rnaseq_stable_ids'],
                      },
      },


      {
        -logic_name => 'generate_rnaseq_stable_ids',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('generate_stable_ids_script').
                                ' -user '.$self->o('rnaseq_final_db','-user').
                                ' -pass '.$self->o('rnaseq_final_db','-pass').
                                ' -host '.$self->o('rnaseq_final_db','-host').
                                ' -port '.$self->o('rnaseq_final_db','-port').
                                ' -dbname '.$self->o('rnaseq_final_db','-dbname').
                                ' -start RNASEQ > #query_file#',
                          query_file => $self->o('rnaseq_stable_id_file'),
                       },
        -rc_name => 'default',
        -flow_into => {
                        1 => ['load_rnaseq_stable_ids'],
                      },
      },


      {
        -logic_name => 'load_rnaseq_stable_ids',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
        -parameters => {
                         db_conn => $self->o('rnaseq_final_db'),
                         input_file => $self->o('rnaseq_stable_id_file'),
                       },
        -rc_name => 'default',
        -flow_into => {
                        1 => ['populate_production_tables_rnaseq_final'],
                      },
      },


      {
        -logic_name => 'populate_production_tables_rnaseq_final',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HivePopulateProductionTables',
        -parameters => {
                         'target_db'        => $self->o('rnaseq_final_db'),
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
                         input_query => 'SELECT * FROM dna_align_feature',
                         command_out => q(sort -nk2 -nk3 -nk4 | sed 's/NULL/\\N/g'),
                         output_file => $self->o('rnaseq_daf_introns_file'),
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
                         db_conn => $self->o('rnaseq_final_db'),
                         input_file => $self->o('rnaseq_daf_introns_file'),
                       },
        -rc_name => 'default',
        -flow_into => {
                        1 => ['set_rnaseq_final_meta_coords'],
                      },
      },


      {
        -logic_name => 'set_rnaseq_final_meta_coords',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('meta_coord_script').
                                ' -user '.$self->o('user').
                                ' -pass '.$self->o('password').
                                ' -host '.$self->o('rnaseq_final_db','-host').
                                ' -port '.$self->o('rnaseq_final_db','-port').
                                ' -dbpattern '.$self->o('rnaseq_final_db','-dbname')
                       },
        -rc_name => 'default',
        -flow_into => {
                        1 => ['set_rnaseq_final_meta_levels'],
                      },
      },


      {
        -logic_name => 'set_rnaseq_final_meta_levels',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('meta_levels_script').
                                ' -user '.$self->o('user').
                                ' -pass '.$self->o('password').
                                ' -host '.$self->o('rnaseq_final_db','-host').
                                ' -port '.$self->o('rnaseq_final_db','-port').
                                ' -dbname '.$self->o('rnaseq_final_db','-dbname')
                       },
        -rc_name => 'default',
        -flow_into => { 1 => ['optimise_rnaseq_final'] },
      },


      {
        -logic_name => 'optimise_rnaseq_final',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('load_optimise_script').
                                ' -output_path '.$self->o('output_path').'/optimise_rnaseq_final/'.
                                ' -uniprot_filename '.$self->o('uniprot_entry_loc').
                                ' -dbuser '.$self->o('user').
                                ' -dbpass '.$self->o('password').
                                ' -dbport '.$self->o('rnaseq_final_db','-port').
                                ' -dbhost '.$self->o('rnaseq_final_db','-host').
                                ' -dbname '.$self->o('rnaseq_final_db','-dbname').
                                ' -prod_dbuser '.$self->o('user_r').
                                ' -prod_dbhost '.$self->o('production_db','-host').
                                ' -prod_dbname '.$self->o('production_db','-dbname').
                                ' -prod_dbport '.$self->o('production_db','-port').
                                ' -nodaf -ise'
                       },
        -max_retry_count => 1,
        -rc_name => '4GB',
        -flow_into => {
                        1 => ['rnaseq_final_gene_sanity_checks'],
                      },
      },


      {
        -logic_name => 'rnaseq_final_gene_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         target_db => $self->o('rnaseq_final_db'),
                         sanity_check_type => 'gene_db_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'gene_db_checks')->{'rnaseq_final'},
                       },

        -rc_name    => '4GB',
        -flow_into => {
                        1 => ['rnaseq_healthchecks'],
                      },
      },


      {
        -logic_name => 'rnaseq_healthchecks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveHealthcheck',
        -parameters => {
			 input_db         => $self->o('rnaseq_final_db'),
                         dna_db           => $self->o('dna_db'),
                         refine_db        => $self->o('rnaseq_refine_db'),
                         blast_db         => $self->o('rnaseq_blast_db'),
			 species          => $self->o('species_name'),
                         group            => 'rnaseq_handover',
                         enscode_root_dir => $self->o('enscode_root_dir'),
                       },

        -rc_name    => '4GB',
      },

    ];
}


sub resource_classes {
  my $self = shift;

  return {
    '1.5GB' => { LSF => $self->lsf_resource_builder('production-rh7', 1500, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '4GB' => { LSF => $self->lsf_resource_builder('production-rh7', 4000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '8GB' => { LSF => $self->lsf_resource_builder('production-rh7', 8000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    '16GB' => { LSF => $self->lsf_resource_builder('production-rh7', 16000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'default' => { LSF => $self->lsf_resource_builder('production-rh7', 900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'default_himem' => { LSF => $self->lsf_resource_builder('production-rh7', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'repeatmasker' => { LSF => $self->lsf_resource_builder('production-rh7', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'repeatmasker_rebatch' => { LSF => $self->lsf_resource_builder('production-rh7', 5900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'simple_features' => { LSF => $self->lsf_resource_builder('production-rh7', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'genscan' => { LSF => $self->lsf_resource_builder('production-rh7', 3900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'refseq_cdna' => { LSF => $self->lsf_resource_builder('production-rh7', 4000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'refseq_cdna_retry' => { LSF => $self->lsf_resource_builder('production-rh7', 6000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'genscan_short' => { LSF => $self->lsf_resource_builder('production-rh7', 5900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'blast' => { LSF => $self->lsf_resource_builder('production-rh7', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], undef, 3)},
    'blast_retry' => { LSF => $self->lsf_resource_builder('production-rh7', 5900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], undef, 3)},
    'rfam_blast' => { LSF => $self->lsf_resource_builder('production-rh7', 4000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], undef, 3)},
    'rfam_blast_retry' => { LSF => $self->lsf_resource_builder('production-rh7', 6000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], undef, 3)},
    'genblast' => { LSF => $self->lsf_resource_builder('production-rh7', 1900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genblast_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'genblast_retry' => { LSF => $self->lsf_resource_builder('production-rh7', 4900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genblast_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'project_transcripts' => { LSF => $self->lsf_resource_builder('production-rh7', 1900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'projection_coding_db_server'}, $self->default_options->{'projection_lastz_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'refseq_import' => { LSF => $self->lsf_resource_builder('production-rh7', 9900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'refseq_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'layer_annotation' => { LSF => $self->lsf_resource_builder('production-rh7', 3900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genblast_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'genebuilder' => { LSF => $self->lsf_resource_builder('production-rh7', 1900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genblast_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'transcript_finalisation' => { LSF => $self->lsf_resource_builder('production-rh7', 1900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genblast_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'filter' => { LSF => $self->lsf_resource_builder('production-rh7', 4900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genblast_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'exonerate' => { LSF => $self->lsf_resource_builder('production-rh7', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'exonerate_6G' => { LSF => $self->lsf_resource_builder('production-rh7', 5900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
  }
}

sub hive_capacity_classes {
  my $self = shift;

  return {
           'hc_low'    => 200,
           'hc_medium' => 500,
           'hc_high'   => 1000,
         };
}

1;
