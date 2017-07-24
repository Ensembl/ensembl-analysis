=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

package Genome_annotation_static_conf;

use strict;
use warnings;
use File::Spec::Functions;

use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_analysis_settings);
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
    use_tokens => 0,
########################
# Misc setup info
########################
    'email_address'             => '' || $ENV{HIVE_EMAIL}, #!!!!!!!!!!!
    'genebuilder_id'            => '' || $ENV{GENEBUILDER_ID}, #!!!!!!!!!!!
    'release_number'            => '', #!!!!!!!!!!! (you can put whatever here if you're not sure)
    'enscode_root_dir'          => '', #!!!!!!!!!!! git repo checkouts
    'repbase_logic_name'        => '', #!!!!!!!!!!!!!!!!! repbase logic name i.e. repeatmask_repbase_XXXX, ONLY FILL THE XXXX BIT HERE!!! e.g primates
    'repbase_library'           => '', #!!!!!!!!!!!!!!!!! repbase library name, this is the actual repeat repbase library to use, e.g. "Mus musculus"
    'species_name'              => '', #!!!!!!!!!!!!!!!!! e.g. mus_musculus
    'production_name'           => '',
    'taxon_id'                  => '', #!!!!!!!!!!!!!!!!! should be in the assembly report file
    'uniprot_set'               => '', #!!!!!!!!!!!!!!!!! Check sub uniprot_clade_download below for suitable set
    'output_path'               => '', # Lustre output dir
    'wgs_id'                    => '', #!!!!!!!!!!!!!! Can be found in assembly report file on ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/
    'assembly_name'             => '', #!!!!!!!!!!!!!! Name (as it appears in the assembly report file)
    'assembly_accession'        => '', #!!!!!!!!!!!!!! GCA
    'assembly_refseq_accession' => '', #!!!!!!!!!!!!!! GCF
    'mt_accession'              => undef, # This should be set to undef unless you know what you are doing
    'skip_projection'           => 0, # Will ignore projection biotypes when making gene set
    'skip_rnaseq'               => 0, # Will ignore rnaseq biotypes when making gene set
    'skip_mito'                 => 0,
    'mapping_required'          => 0,
    'mapping_db'                => '',
    'stable_id_prefix'          => '',
    'stable_id_start'           => '0',
    'uniprot_db_dir'            => 'uniprot_2017_04', # e.g. 'uniprot_2017_04'
    'mirBase_fasta'             => 'human_mirnas.fa',
    'ig_tr_fasta_file'          => 'human_ig_tr.fa',

########################
# Pipe and ref db info
########################
    'pipeline_name'                => '', #!!!!!!!!!!! What you want hive to call the pipeline, not the db name itself
    'user_r'                       => '', #!!!!!!!!!!!
    'user'                         => '', #!!!!!!!!!!!
    'password'                     => '', #!!!!!!!!!!!

    'pipe_db_server'               => '', #!!!!!!!!!!!
    'databases_server'             => '', #!!!!!!!!!!!
    'dna_db_server'                => '', #!!!!!!!!!!!

    'pipe_db_port'                 => '', #!!!!!!!!!!!
    'databases_port'               => '', #!!!!!!!!!!!
    'dna_db_port'                  => '', #!!!!!!!!!!!

    'projection_source_db_name'    => '', # This is generally a pre-existing db, like the current human/mouse core for example
    'projection_source_db_server'  => 'mysql-ensembl-mirror',
    'projection_source_db_port'    => 4240,

    'pipe_dbname' => $self->o('dbowner').'_'.$self->o('pipeline_name').'_pipe_'.$self->o('release_number'),
    'dna_dbname'                   => $self->o('dbowner').'_'.$self->o('production_name').'_core_'.$self->o('release_number'),
    'port'                         => $self->o('pipe_db_port'),

    'reference_dbname'             => $self->o('dna_dbname'),
    'reference_db_server'          => $self->o('dna_db_server'),
    'reference_db_port'            => $self->o('dna_db_port'),

    'cdna_db_server'               => $self->o('databases_server'),
    'cdna_db_port'                 => $self->o('databases_port'),

    'genblast_db_server'           => $self->o('databases_server'),
    'genblast_db_port'             => $self->o('databases_port'),

    'ig_tr_db_server'              => $self->o('databases_server'),
    'ig_tr_db_port'                => $self->o('databases_port'),

    'genewise_db_server'           => $self->o('databases_server'),
    'genewise_db_port'             => $self->o('databases_port'),

    'projection_db_server'         => $self->o('databases_server'),
    'projection_db_port'           => $self->o('databases_port'),

    'projection_realign_db_server' => $self->o('databases_server'),
    'projection_realign_db_port'   => $self->o('databases_port'),

    'projection_lastz_db_server'   => $self->o('databases_server'),
    'projection_lastz_db_port'     => $self->o('databases_port'),

    'rnaseq_dbname'                => $self->o('dbowner').'_'.$self->o('species_name').'_rnaseq_'.$self->o('release_number'),
    'rnaseq_db_server'             => $self->o('databases_server'),
    'rnaseq_db_port'               => $self->o('databases_port'),

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

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'staging_1_db_server'          => 'mysql-ens-sta-1',
    'staging_1_port'               => 4519,


    databases_to_delete => ['reference_db', 'cdna_db', 'genblast_db', 'genewise_db', 'projection_db', 'projection_realign_db', 'layering_db', 'utr_output_db', 'genebuilder_db', 'pseudogene_db', 'ncrna_db', 'final_geneset_db', 'refseq_db'],


########################
# BLAST db paths
########################
    'base_blast_db_path'        => '' || $ENV{BLASTDB_DIR},
    'uniprot_entry_loc'         => catfile($self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_db_dir'), '/entry_loc'),
    'uniprot_blast_db_path'     => catfile($self->o('base_blast_db_path'), 'uniprot', $self->o('uniprot_db_dir'), '/uniprot_vertebrate'),
    'vertrna_blast_db_path'     => catfile($self->o('base_blast_db_path'), 'vertrna', '131/embl_vertrna-1'),
    'unigene_blast_db_path'     => catfile($self->o('base_blast_db_path'), 'unigene', 'unigene'),
    'mito_index_path'           => undef, # Set this path only if you don't want to use the GCF report and if you haven't set 'mt_accession' '/nfs/production/panda/ensembl/genebuild/blastdb/refseq_mitochondria_set/mito_index.txt',
    'ncrna_blast_path'          => catfile($self->o('base_blast_db_path'), 'ncrna', 'ncrna_2016_05'),
    'ig_tr_blast_path'          => catfile($self->o('base_blast_db_path'), 'ig_tr_genes'),

######################################################
#
# Mostly constant settings
#
######################################################

    'genome_file'               => catfile($self->o('output_path'), 'genome_dumps', $self->o('species_name').'_softmasked_toplevel.fa'),
    'primary_assembly_dir_name' => 'Primary_Assembly',
    'refseq_cdna_calculate_coverage_and_pid' => '0',
    'refseq_cdna_table_name'    => 'refseq_sequences',
    'refseq_cdna_batch_size'    => 10,
    'contigs_source'            => 'ena',



    'layering_input_gene_dbs' => [
                                   $self->o('genblast_db'),
                                   $self->o('rnaseq_db'),
                                   $self->o('projection_db'),
                                   $self->o('ig_tr_db'),
                                 ],

    'utr_gene_dbs' => {
                        'cdna_db'       => $self->o('cdna_db'),
                        'rnaseq_utr_db' => $self->o('rnaseq_db'),
                        'no_utr_db'     => $self->o('layering_db'),
                      },

    'utr_biotype_priorities'  => {
                                   'rnaseq' => 1,
                                   'cdna' => 1,
                                   'cdna_predicted' => 2,
                                 },

    'cleaning_blessed_biotypes' => {'pseudogene' => 1, 'processed_pseudogene' => 1},

    'min_toplevel_slice_length'   => 250,

    'repeat_logic_names'          => ['repeatmask_repbase_'.$self->o('repbase_logic_name'),'dust'],
    'homology_models_path'        => catdir($self->o('output_path'),'homology_models'),
    'load_fasta_script_path'        => catfile($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'genebuild', 'load_fasta_to_db_table.pl'),
    'clone_db_script_path'      => catfile($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'clone_database.ksh'),
    'refseq_synonyms_script_path' => catfile($self->o('enscode_root_dir'), 'ensembl-pipeline', 'scripts', 'refseq_import', 'load_refseq_synonyms.pl'),
    'refseq_import_script_path'   => catfile($self->o('enscode_root_dir'), 'ensembl-pipeline', 'scripts', 'refseq_import', 'parse_ncbi_gff3.pl'),
    'loading_report_script'     => catfile($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'genebuild', 'report_genome_prep_stats.pl'),
    'remove_duplicates_script_path' => catfile($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts','find_and_remove_duplicates.pl'),
    repeat_types_script => catfile($self->o('enscode_root_dir'), 'ensembl', 'misc-scripts', 'repeats', 'repeat-types.pl'),
    meta_coord_script => catfile($self->o('enscode_root_dir'), 'ensembl', 'misc-scripts', 'meta_coord', 'update_meta_coord.pl'),
    meta_levels_script => catfile($self->o('enscode_root_dir'), 'ensembl', 'misc-scripts', 'meta_levels.pl'),

########################
# Extra db settings
########################

    'num_tokens' => 10,

########################
# Executable paths
########################
    'blast_type' => 'ncbi', # It can be 'ncbi', 'wu', or 'legacy_ncbi'
    'binary_base' => '/nfs/software/ensembl/RHEL7/linuxbrew/bin',
    'dust_path' => catfile($self->o('binary_base'), 'dustmasker'),
    'trf_path' => catfile($self->o('binary_base'), 'trf'),
    'eponine_java_path' => '/nfs/software/ensembl/RHEL7/jenv/shims/java',
    'eponine_jar_path' => '/nfs/software/ensembl/RHEL7/linuxbrew/Cellar/eponine/1.0/libexec/eponine-scan.jar',
    'cpg_path' => catfile($self->o('binary_base'), 'cpg_lh'),
    'trnascan_path' => catfile($self->o('binary_base'), 'tRNAscan-SE'),
    'repeatmasker_path' => catfile($self->o('binary_base'), 'RepeatMasker'),
    'genscan_path' => catfile($self->o('binary_base'), 'genscan'),
    'genscan_matrix_path' => '/nfs/software/ensembl/RHEL7/linuxbrew/share/HumanIso.smat',
    'uniprot_blast_exe_path' => catfile($self->o('binary_base'), 'blastp'),
    'blastn_exe_path' => catfile($self->o('binary_base'), 'blastn'),
    'vertrna_blast_exe_path' => catfile($self->o('binary_base'), 'tblastn'),
    'unigene_blast_exe_path' => catfile($self->o('binary_base'), 'tblastn'),
    'exonerate_path'         => '/nfs/software/ensembl/RHEL7/linuxbrew/opt/exonerate09/bin/exonerate',
    'cmsearch_exe_path'    => '/nfs/production/panda/ensembl/leanne/programs/infernal/infernal-1.0/src/cmsearch',

    'uniprot_genblast_batch_size' => 15,
    'uniprot_table_name'          => 'uniprot_sequences',

    'genblast_path'     => catfile($self->o('binary_base'), 'genblast'),
    'genblast_eval'     => $self->o('blast_type') eq 'wu' ? '1e-20' : '1e-1',
    'genblast_cov'      => '0.5',
    'genblast_pid'      => '50',
    'genblast_max_rank' => '5',

    'ig_tr_table_name'    => 'ig_tr_sequences',
    'ig_tr_genblast_cov'  => '0.8',
    'ig_tr_genblast_pid'  => '70',
    'ig_tr_genblast_eval' => '1',
    'ig_tr_genblast_max_rank' => '10',
    'ig_tr_batch_size'    => 10,

    'exonerate_cdna_pid' => '95', # Cut-off for percent id
    'exonerate_cdna_cov' => '50', # Cut-off for coverage

# Max internal stops for projected transcripts
    'projection_pid'                        => '50',
    'projection_cov'                        => '50',
    'projection_max_internal_stops'         => '1',
    'projection_calculate_coverage_and_pid' => '1',

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
    'killlist_dbname' => 'gb_kill_list',

    'reference_db' => {
      -dbname => $self->o('reference_dbname'),
      -host   => $self->o('reference_db_server'),
      -port   => $self->o('reference_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'cdna_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('species_name').'_cdna_'.$self->o('release_number'),
      -host   => $self->o('cdna_db_server'),
      -port   => $self->o('cdna_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },


    'genblast_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('species_name').'_genblast_'.$self->o('release_number'),
      -host   => $self->o('genblast_db_server'),
      -port   => $self->o('genblast_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'ig_tr_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('species_name').'_igtr_'.$self->o('release_number'),
      -host   => $self->o('ig_tr_db_server'),
      -port   => $self->o('ig_tr_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'genewise_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('species_name').'_genewise_'.$self->o('release_number'),
      -host   => $self->o('genewise_db_server'),
      -port   => $self->o('genewise_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'projection_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('species_name').'_proj_'.$self->o('release_number'),
      -host   => $self->o('projection_db_server'),
      -port   => $self->o('projection_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'projection_realign_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('species_name').'_realign_'.$self->o('release_number'),
      -host   => $self->o('projection_realign_db_server'),
      -port   => $self->o('projection_realign_db_port'),
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
      -dbname => $self->o('dbowner').'_'.$self->o('species_name').'_lastz_'.$self->o('release_number'),
      -host   => $self->o('projection_lastz_db_server'),
      -port   => $self->o('projection_lastz_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -driver => $self->o('hive_driver'),
    },

    'rnaseq_db' => {
      -dbname => $self->o('rnaseq_dbname'),
      -host   => $self->o('rnaseq_db_server'),
      -port   => $self->o('rnaseq_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'layering_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('species_name').'_layer_'.$self->o('release_number'),
      -host   => $self->o('layering_db_server'),
      -port   => $self->o('layering_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'utr_output_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('species_name').'_utr_'.$self->o('release_number'),
      -host   => $self->o('utr_db_server'),
      -port   => $self->o('utr_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'genebuilder_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('species_name').'_gbuild_'.$self->o('release_number'),
      -host   => $self->o('genebuilder_db_server'),
      -port   => $self->o('genebuilder_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'pseudogene_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('species_name').'_pseudo_'.$self->o('release_number'),
      -host   => $self->o('pseudogene_db_server'),
      -port   => $self->o('pseudogene_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'ncrna_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('species_name').'_ncrna_'.$self->o('release_number'),
      -host   => $self->o('ncrna_db_server'),
      -port   => $self->o('ncrna_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'final_geneset_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('species_name').'_final_'.$self->o('release_number'),
      -host   => $self->o('final_geneset_db_server'),
      -port   => $self->o('final_geneset_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'refseq_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('species_name').'_refseq_'.$self->o('release_number'),
      -host   => $self->o('refseq_db_server'),
      -port   => $self->o('refseq_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'killlist_db' => {
      -dbname => $self->o('killlist_dbname'),
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
      -dbname => 'ensembl_production',
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
    };
}

sub pipeline_create_commands {
    my ($self) = @_;
    return [
    # inheriting database and hive tables' creation
      @{$self->SUPER::pipeline_create_commands},

      $self->db_cmd('CREATE TABLE '.$self->o('uniprot_table_name').' ('.
                    'accession varchar(50) NOT NULL,'.
                    'source_db varchar(50) NOT NULL,'.
                    'pe_level varchar(50) NOT NULL,'.
                    'biotype varchar(255) NOT NULL,'.
                    'group_name varchar(255) NOT NULL,'.
                    'seq text NOT NULL,'.
                    'PRIMARY KEY (accession))'),

      $self->db_cmd('CREATE TABLE '.$self->o('refseq_cdna_table_name').' ('.
                    'accession varchar(50) NOT NULL,'.
                    'source_db varchar(50) NOT NULL,'.
                    'biotype varchar(25) NOT NULL,'.
                    'date varchar(50) NOT NULL,'.
                    'seq text NOT NULL,'.
                    'PRIMARY KEY (accession))'),
    ];
}


## See diagram for pipeline structure
sub pipeline_analyses {
    my ($self) = @_;

    my %genblast_params = (
      wu    => '-P wublast -gff -e '.$self->o('genblast_eval').' -c '.$self->o('genblast_cov'),
      ncbi  => '-P blast -gff -e '.$self->o('genblast_eval').' -c '.$self->o('genblast_cov').' -W 3 -rl 5000',
      ig_tr => '-P blast -gff -e '.$self->o('ig_tr_genblast_eval').' -c '.$self->o('ig_tr_genblast_cov').' -W 3 -rl 5000',
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
        # Download the files and dir structure from the NCBI ftp site. Uses the link to a species in the ftp_link_file
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
        -input_ids  => [{assembly_name => $self->o('assembly_name'), assembly_accession => $self->o('assembly_accession'), assembly_refseq_accession => $self->o('assembly_refseq_accession')}],
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
        -rc_name    => 'default',
        -flow_into  => {
                         1 => ['download_contigs'],
                       },
      },

            {
        # Download contig from NCBI
        -logic_name => 'download_contigs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveDownloadContigs',
        -parameters => {
                         'contigs_source'            => $self->o('contigs_source'),
                         'wgs_id'                    => $self->o('wgs_id'),
                         'output_path'               => $self->o('output_path'),
                         'primary_assembly_dir_name' => $self->o('primary_assembly_dir_name'),
                       },
        -rc_name    => 'default',
        -flow_into  => {
                         1 => ['create_core_db'],
                       },
      },


      {
        # Creates a reference db for each species
        -logic_name => 'create_core_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         'target_db'        => $self->o('reference_db'),
                         'user_w'           => $self->o('user'),
                         'pass_w'           => $self->o('password'),
                         'enscode_root_dir' => $self->o('enscode_root_dir'),
                         'create_type'      => 'core_only',
                       },
        -rc_name    => 'default',
        -flow_into  => {
                         1 => ['populate_production_tables'],
                       },

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
                         1 => ['load_contigs'],
                       },
      },


      {
        # Load the contigs into each reference db
        -logic_name => 'load_contigs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadSeqRegions',
        -parameters => {
                         'coord_system_version'      => $self->o('assembly_name'),
                         'target_db'                 => $self->o('reference_db'),
                         'output_path'               => $self->o('output_path'),
                         'enscode_root_dir'          => $self->o('enscode_root_dir'),
                         'primary_assembly_dir_name' => $self->o('primary_assembly_dir_name'),
                       },
        -rc_name    => 'default',
        -flow_into  => {
                         1 => ['load_assembly_info'],
                       },
      },


      {
        # Load the AGP files
        -logic_name => 'load_assembly_info',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadAssembly',
        -parameters => {
                         'target_db'                 => $self->o('reference_db'),
                         'output_path'               => $self->o('output_path'),
                         'enscode_root_dir'          => $self->o('enscode_root_dir'),
                         'primary_assembly_dir_name' => $self->o('primary_assembly_dir_name'),
                       },

        -rc_name    => 'default',
        -flow_into  => {
                         1 => ['set_toplevel'],
                       },
      },


      {
        # Set the toplevel
        -logic_name => 'set_toplevel',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveSetAndCheckToplevel',
        -parameters => {
                         'target_db'            => $self->o('reference_db'),
                         'output_path'          => $self->o('output_path'),
                         'enscode_root_dir'     => $self->o('enscode_root_dir'),
                         'primary_assembly_dir_name' => $self->o('primary_assembly_dir_name'),
                       },
        -rc_name    => '1.5GB',
        -flow_into  => {
                         1 => ['load_meta_info'],
                       },
      },


      {
        # Load some meta info and seq_region_synonyms
        -logic_name => 'load_meta_info',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveSetMetaAndSeqRegionSynonym',
        -parameters => {
                         'target_db'                 => $self->o('reference_db'),
                         'output_path'               => $self->o('output_path'),
                         'enscode_root_dir'          => $self->o('enscode_root_dir'),
                         'primary_assembly_dir_name' => $self->o('primary_assembly_dir_name'),
                         'meta_key_list' => {
                                              'assembly.accession'                      => $self->o('assembly_accession'),
                                              'assembly.coverage_depth'                 => 'high',
                                              'assembly.default'                        => $self->default_options->{'assembly_name'},
                                              'assembly.name'                           => $self->default_options->{'assembly_name'},
                                              'assembly.web_accession_source'           => 'NCBI',
                                              'assembly.web_accession_type'             => 'GenBank Assembly ID',
                                              'genebuild.id'                          => $self->o('genebuilder_id'),
                                              'genebuild.method'                        => 'full_genebuild',
                                              'genebuild.projection_source_db'          => $self->default_options->{'projection_source_db_name'},
                                              'provider.name'                           => 'Ensembl',
                                              'provider.url'                            => 'www.ensembl.org',
                                              'removed_evidence_flag.ensembl_dbversion' => $self->default_options->{'projection_source_db_name'},
                                              'removed_evidence_flag.uniprot_dbversion' => $self->o('uniprot_db_dir'),
                                              'repeat.analysis'                         => 'repeatmask_repbase_'.$self->o('repbase_logic_name'),
                                              'repeat.analysis'                         => 'dust',
                                              'repeat.analysis'                         => 'trf',
                                              'species.production_name'                 => $self->o('production_name'),
                                              'species.taxonomy_id'                     => $self->default_options->{'taxon_id'},
                                            }
                       },
        -rc_name    => 'default',
        -flow_into  => {
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
                          1 => ['load_mitochondrion', 'fan_refseq'],
                       },
      },

      {
        -logic_name => 'fan_refseq',
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
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('refseq_synonyms_script_path').
                                      ' -species '.$self->o('species_name').
                                      ' -dbuser '.$self->o('user').
                                      ' -dbpass '.$self->o('password').
                                      ' -dbhost '.$self->o('reference_db','-host').
                                      ' -dbport '.$self->o('reference_db','-port').
                                      ' -dbname '.$self->o('reference_db','-dbname').
                                      ' -workdir '.$self->o('output_path').'/refseq_import/',
                       },
        -rc_name => 'default',
        -flow_into  => {
                          1 => ['download_refseq_gff'],
                       },
      },

      {
        -logic_name => 'load_mitochondrion',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadMitochondrion',
        -parameters => {
                         'skip_analysis'             => $self->o('skip_mito'),
                         'target_db'                 => $self->o('reference_db'),
                         'output_path'               => $self->o('output_path'),
                         'enscode_root_dir'          => $self->o('enscode_root_dir'),
                         'mito_index_path'           => $self->o('mito_index_path'),
                         'species_name'              => $self->o('species_name'),
                         'mt_accession'              => $self->o('mt_accession'),
                      },
        -rc_name    => 'default',

        -flow_into => {
                        '1->A' => ['create_10mb_slice_ids'],
                        'A->1' => ['genome_prep_sanity_checks'],
                      },

      },


###############################################################################
#
# REFSEQ ANNOTATION
#
###############################################################################

      {
        -logic_name => 'download_refseq_gff',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'mkdir -p '.catdir($self->o('output_path'), 'refseq_import').';'.
                                'wget -O '.catfile($self->o('output_path'), 'refseq_import', 'refseq.gff.gz').' '.$self->o('refseq_import_ftp_path').';'.
                                'gunzip '.catfile($self->o('output_path'), 'refseq_import', 'refseq.gff.gz'),
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
                         script_path => $self->o('clone_db_script_path'),
                         user_r => $self->o('user_r'),
                         user => $self->o('user'),
                         pass_w => $self->o('password'),
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
                                      ' -infile '.$self->o('output_path').'/refseq_import/refseq.gff',
                       },
        -rc_name => 'refseq_import',
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
                        '2' => ['semaphore_10mb_slices'],
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
                        '2' => ['run_repeatmasker'],
                      },
      },


      {
        -logic_name => 'run_repeatmasker',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveRepeatMasker',
        -parameters => {
                         timer_batch => '3h',
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
                         timer_batch => '1h',
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
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         'source_db'        => $self->o('dna_db'),
                         'user_w'           => $self->o('user'),
                         'pass_w'           => $self->o('password'),
                         'create_type'      => 'backup',
                         'output_path'      => $self->o('output_path')."/core_db_bak/",
                         'backup_name'      => 'core_bak.sql',
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
         -flow_into => { 1 => ['create_cdna_db'] },
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
                         script_path => $self->o('clone_db_script_path'),
                         user_r => $self->o('user_r'),
                         user => $self->o('user'),
                         pass_w => $self->o('password'),
                       },
        -rc_name    => 'default',
        -wait_for => ['format_softmasked_toplevel'],
        -flow_into => {
                        '1->A' => ['fan_refseq_2'],
                        'A->1' => ['create_genblast_output_db'],
                      },
      },


      {
        -logic_name => 'fan_refseq_2',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'if [ -n "#assembly_refseq_accession#" ]; then exit 0; else exit 42;fi',
                         return_codes_2_branches => {'42' => 2},
                       },
        -rc_name => 'default',
        -flow_into  => {
                          1 => ['get_refseq_cdnas'],
                       },
      },


      {
        -logic_name => 'get_refseq_cdnas',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadRefSeqmRNA',
        -parameters => {
                         output_path => $self->o('output_path')."/cdnas/",
                         ftp_path    => $self->o('refseq_mrna_ftp_path'),
                         min_seq_length => 60,
                       },
        -flow_into => {
                        1 => ['generate_refseq_cdna_jobs'],
                      },
        -rc_name    => 'default',
      },


      {
        -logic_name => 'generate_refseq_cdna_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         iid_type => 'sequence_accession',
                         batch_size => $self->o('refseq_cdna_batch_size'),
                         sequence_table_name => $self->o('refseq_cdna_table_name'),
                       },
        -rc_name      => 'default_himem',
        -flow_into => {
                        2 => ['refseq_cdna'],
                      },
      },


      {
        -logic_name => 'refseq_cdna',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
        -parameters => {
                         iid_type => 'db_seq',
                         sequence_table_name => $self->o('refseq_cdna_table_name'), # there is a problem here. I add that, but it should be query_table_name or sequence_table_name

                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('cdna_db'),
                         logic_name => 'refseq_cdna',
                         module     => 'HiveExonerate2Genes',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_cov_per_sub')},
                         GENOMICSEQS         => $self->o('genome_file'),
                         PROGRAM             => $self->o('exonerate_path'),
                         SOFT_MASKED_REPEATS => ['repeatmask_repbase_'.$self->o('repbase_logic_name'),'dust'],
                         query_seq_dir => $self->o('output_path'),
                         calculate_coverage_and_pid => $self->o('refseq_cdna_calculate_coverage_and_pid'),
                         exonerate_cdna_pid => $self->o('exonerate_cdna_pid'),
                         exonerate_cdna_cov => $self->o('exonerate_cdna_cov'),
                      },
        -flow_into => {
                        -1 => ['split_refseq_cdna_jobs'],
                        -2 => ['split_refseq_cdna_jobs'],
                      },
        -batch_size => 20,
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -rc_name    => 'refseq_cdna',
        -failed_job_tolerance => 0.5,
     },


     {
       -logic_name => 'split_refseq_cdna_jobs',
       -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
       -parameters => {
                        iid_type => 'rechunk',
                        batch_size => 1,
                      },
       -rc_name      => 'default',
       -can_be_empty  => 1,
       -flow_into => {
                       2 => ['refseq_cdna_retry'],
                     },
     },


     {
        -logic_name => 'refseq_cdna_retry',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
        -rc_name    => 'genblast_retry',
        -parameters => {
                         iid_type => 'db_seq',
                         sequence_table_name => $self->o('refseq_cdna_table_name'), # there is a problem here. I add that, but it should be query_table_name or sequence_table_name
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('cdna_db'),
                         logic_name => 'refseq_cdna',
                         module     => 'HiveExonerate2Genes',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_cov_per_sub')},
                         GENOMICSEQS         => $self->o('genome_file'),
                         PROGRAM             => $self->o('exonerate_path'),
                         SOFT_MASKED_REPEATS => ['repeatmask_repbase_'.$self->o('repbase_logic_name'),'dust'],
                         query_seq_dir => $self->o('output_path'),
                         calculate_coverage_and_pid => $self->o('refseq_cdna_calculate_coverage_and_pid'),
                         exonerate_cdna_pid => $self->o('exonerate_cdna_pid'),
                         exonerate_cdna_cov => $self->o('exonerate_cdna_cov'),
                      },
        -flow_into => {
                        -1 => ['failed_refseq_cdnas'],
                        -2 => ['failed_refseq_cdnas'],
                      },
        -batch_size => 20,
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -failed_job_tolerance => 100,
        -can_be_empty => 1,
        -rc_name    => 'refseq_cdna_retry',
      },


      {
        -logic_name => 'failed_refseq_cdnas',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {},
        -rc_name     => 'default',
        -can_be_empty  => 1,
        -failed_job_tolerance => 100,
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
                         script_path => $self->o('clone_db_script_path'),
                         user_r => $self->o('user_r'),
                         user => $self->o('user'),
                         pass_w => $self->o('password'),
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
                      },
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
                         commandline_params => $genblast_params{$self->o('blast_type')},
                         sequence_table_name => $self->o('uniprot_table_name'),
                         max_rank => $self->o('genblast_max_rank'),
                         genblast_pid => $self->o('genblast_pid'),
                         timer => '2h',
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
                         commandline_params => $genblast_params{$self->o('blast_type')},
                         sequence_table_name => $self->o('uniprot_table_name'),
                         max_rank => $self->o('genblast_max_rank'),
                         genblast_pid => $self->o('genblast_pid'),
                         timer => '1h',
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
        -flow_into => {
                        1 => ['create_ig_tr_db'],
                      },

      },


      {
        -logic_name => 'create_ig_tr_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                     source_db => $self->o('dna_db'),
                     target_db => $self->o('ig_tr_db'),
                     create_type => 'clone',
                     script_path => $self->o('clone_db_script_path'),
                     user_r => $self->o('user_r'),
                     user => $self->o('user'),
                     pass_w => $self->o('password'),
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
                                ." -create_table 1"
                       },

         -rc_name => 'default',
         -flow_into => {
                         '1->A' => ['generate_ig_tr_jobs'],
                         'A->1' => ['cluster_ig_tr_genes'],
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
                        2 => ['ig_tr_genblast'],
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
                         commandline_params => $genblast_params{$self->o('blast_type')},
                         sequence_table_name => $self->o('ig_tr_table_name'),
                         max_rank => $self->o('ig_tr_genblast_max_rank'),
                         genblast_pid => $self->o('ig_tr_genblast_pid'),
                         timer => '2h',
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
                         commandline_params => $genblast_params{$self->o('blast_type')},
                         sequence_table_name => $self->o('uniprot_table_name'),
                         max_rank => $self->o('genblast_max_rank'),
                         genblast_pid => $self->o('genblast_pid'),
                         timer => '1h',
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
        -logic_name => 'cluster_ig_tr_genes',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCollapseIGTR',
        -parameters => {
                         target_db => $self->o('ig_tr_db'),
                         dna_db => $self->o('dna_db'),
                         logic_name => 'ig_tr_gene',
                         logic_names_to_cluster => ['genblast','genblast_not_best'],
                       },
        -rc_name    => 'genblast',
        -flow_into => {
                        1 => ['update_ig_tr_biotypes'],
                      },
      },

      {
        -logic_name => 'update_ig_tr_biotypes',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'mysql '.
                                ' -u'.$self->o('user').
                                ' -p'.$self->o('password').
                                ' -h'.$self->o('ig_tr_db','-host').
                                ' -P'.$self->o('ig_tr_db','-port').
                                ' -D'.$self->o('ig_tr_db','-dbname').
                                ' -e \'UPDATE gene join analysis using(analysis_id) set biotype=concat(biotype,"_not_best")'.
                                ' where logic_name != "ig_tr_gene"; UPDATE transcript join gene using(gene_id) set transcript.biotype=gene.biotype;\'',
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['ig_tr_sanity_checks'],
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
                         script_path => $self->o('clone_db_script_path'),
                         user_r => $self->o('user_r'),
                         user => $self->o('user'),
                         pass_w => $self->o('password'),
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
                         dna_db => $self->default_options->{'dna_db'},
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
                         dna_db => $self->default_options->{'dna_db'},
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
                         dna_db => $self->default_options->{'dna_db'},
                         logic_name => 'blastmirna',
                         module     => 'HiveBlastmiRNA',
                         blast_db_path => $self->o('ncrna_blast_path') . '/' . $self->o('mirBase_fasta') ,
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
                         dna_db => $self->default_options->{'dna_db'},
                         logic_name => 'blastmirna',
                         module     => 'HiveBlastmiRNA',
                         blast_db_path => $self->o('ncrna_blast_path') . '/' . $self->o('mirBase_fasta') ,
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
                         blast_db_dir_path => $self->o('ncrna_blast_path').'/all_mirnas.embl',
                       },
        -batch_size => 20,
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -rc_name    => 'filter',
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
                        1 => ['create_projection_output_db'],
                      },
      },

########################################################################
#
# Projection analyses
#
########################################################################

      {
        -logic_name => 'create_projection_output_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('projection_db'),
                         create_type => 'clone',
                         script_path => $self->o('clone_db_script_path'),
                         user_r => $self->o('user_r'),
                         user => $self->o('user'),
                         pass_w => $self->o('password'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '1->A' => ['create_projection_input_ids'],
                        'A->1' => ['classify_projection_models'],
                      },
      },

      {
        -logic_name => 'create_projection_input_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         skip_analysis => $self->o('skip_projection'),
                         target_db => $self->o('projection_source_db'),
                         iid_type => 'feature_id',
                         feature_type => 'transcript',
                         feature_restriction => 'projection',
                       },

        -flow_into => {
                        2 => ['project_transcripts'],
                      },

         -rc_name    => 'default',
      },

      {
        -logic_name => 'project_transcripts',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveWGA2GenesDirect',
        -parameters => {
                         logic_name => 'project_transcripts',
                         module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveWGA2GenesDirect',
                         QUERY_CORE_DB => $self->default_options()->{'projection_source_db'},
                         QUERY_CORE_DNA_DB => $self->default_options()->{'projection_source_db'},
                         TARGET_CORE_DB => $self->o('projection_db'),
                         TARGET_CORE_DNA_DB => $self->o('dna_db'),
                         COMPARA_DB => $self->o('projection_lastz_db'),
                         INPUT_METHOD_LINK_TYPE => 'LASTZ_NET',
                         MAX_EXON_READTHROUGH_DIST => 15,
                         TRANSCRIPT_FILTER => {
                           OBJECT     => 'Bio::EnsEMBL::Analysis::Tools::ClusterFilter',
                           PARAMETERS => {
                             -coverage => $self->o('projection_cov'),
                             -percent_id => $self->o('projection_pid'),
                             -max_editable_stops => $self->o('projection_max_internal_stops'),
                             -best_in_genome => 1,
                           },
                         },
                         iid_type => 'feature_id',
                         feature_type => 'transcript',
                         calculate_coverage_and_pid => $self->o('projection_calculate_coverage_and_pid'),
                         max_internal_stops => $self->o('projection_max_internal_stops'),
                       },

        -flow_into => {
                        -3 => ['failed_projection_jobs'],
                      },
        -rc_name          => 'project_transcripts',
        -batch_size => 100,
        -failed_job_tolerance => 0.5,
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
      },

      {
        -logic_name => 'failed_projection_jobs',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name          => 'default',
        -can_be_empty  => 1,
        -failed_job_tolerance => 100,
      },

      {
        -logic_name => 'classify_projection_models',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport',
        -parameters => {
                         skip_analysis => $self->o('skip_projection'),
                         classification_type => 'standard',
                         update_gene_biotype => 1,
                         target_db => $self->o('projection_db'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['fix_unalgined_protein_hit_names'],
                      },
      },

      {
        -logic_name => 'fix_unalgined_protein_hit_names',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
        -parameters => {
                               'db_conn'    => $self->o('projection_db'),
                                input_query => 'UPDATE protein_align_feature join transcript_supporting_feature on feature_id = protein_align_feature_id '.
                                               'join transcript using(transcript_id) SET hit_name = stable_id; '.
                                               'UPDATE protein_align_feature join supporting_feature on feature_id = protein_align_feature_id '.
                                               'join exon_transcript using(exon_id) join transcript using(transcript_id) SET hit_name = stable_id;',
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['projection_sanity_checks'],
                      },
      },

      {
        -logic_name => 'projection_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         skip_check => $self->o('skip_projection'),
                         target_db => $self->o('projection_db'),
                         sanity_check_type => 'gene_db_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'gene_db_checks')->{$self->default_options->{'uniprot_set'}}->{'projection'},
                       },

        -rc_name    => '4GB',
        -flow_into => {
                        1 => ['create_realign_output_db'],
                      },
      },

      {
        -logic_name => 'create_realign_output_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('projection_realign_db'),
                         create_type => 'clone',
                         script_path => $self->o('clone_db_script_path'),
                         user_r => $self->o('user_r'),
                         user => $self->o('user'),
                         pass_w => $self->o('password'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                       '1->A' => ['create_ids_for_evaluate_projection'],
                       'A->1' => ['classify_realigned_models'],
                      },
      },

      {
        -logic_name => 'create_ids_for_evaluate_projection',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         skip_analysis => $self->o('skip_projection'),
                         target_db => $self->o('projection_db'),
                         iid_type => 'feature_id',
                         feature_type => 'transcript',
                       },
        -flow_into => {
                        2 => ['evaluate_projection_transcripts'],
                      },
        -rc_name    => 'default',
      },


      {
        -logic_name => 'evaluate_projection_transcripts',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProjectionRealignment',
        -parameters => {
                         dna_db => $self->o('dna_db'),
                         projection_db => $self->o('projection_db'),
                         projection_source_db => $self->o('projection_source_db'),
                         projection_realign_db => $self->o('projection_realign_db'),
                         protein_table_name => 'projection_protein_sequences',
                       },
        -rc_name    => 'default',
        -flow_into  => {
                         2 => ['realign_projection'],
                       },
      },

      {
        -logic_name => 'realign_projection',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
        -parameters => {
                         iid_type => 'projection_transcript_id',
                         sequence_table_name => 'projection_protein_sequences',
                         projection_padding => 0,
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('projection_realign_db'),
                         projection_db => $self->o('projection_db'),
                         logic_name => 'genblast',
                         module => 'HiveGenblast',
                         genblast_path => $self->o('genblast_path'),
                         commandline_params => ' -P wublast -gff -e '.$self->o('genblast_eval').' -c '.$self->o('genblast_cov').' -n 100 -rl 5000 -x 5 ',
                         query_seq_dir => $self->o('homology_models_path'),
                         max_rank => 1,
                         genblast_pid => $self->o('genblast_pid'),
                       },
        -rc_name    => 'default',
        -flow_into  => {
                         -1 => ['failed_realignments'],
                         -2 => ['failed_realignments'],
                         -3 => ['failed_realignments'],
                       },
      },


      {
        -logic_name => 'failed_realignments',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name          => 'default',
        -can_be_empty  => 1,
        -failed_job_tolerance => 100,
      },

     {
        -logic_name => 'classify_realigned_models',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport',
        -parameters => {
                         classification_type => 'standard',
                         update_gene_biotype => 1,
                         target_db => $self->o('projection_realign_db'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['realign_sanity_checks'],
                      },
      },

      {
        -logic_name => 'realign_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         skip_check => $self->o('skip_projection'),
                         target_db => $self->o('projection_realign_db'),
                         sanity_check_type => 'gene_db_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'gene_db_checks')->{$self->default_options->{'uniprot_set'}}->{'realign'},
                       },

        -rc_name    => '4GB',
        -flow_into => {
                        1 => ['create_rnaseq_db'],
                      },
      },

############################################################################
#
# RNA-seq analyses
#
############################################################################

      {
        -logic_name => 'create_rnaseq_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('dna_db'),
                         target_db => $self->o('rnaseq_db'),
                         create_type => 'clone',
                         script_path => $self->o('clone_db_script_path'),
                         user_r => $self->o('user_r'),
                         user => $self->o('user'),
                         pass_w => $self->o('password'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['classify_rnaseq_models'],
                      },
      },


      {
        -logic_name => 'classify_rnaseq_models',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport',
        -parameters => {
                         skip_analysis => $self->o('skip_rnaseq'),
                         update_gene_biotype => 1,
                         classification_type => 'standard',
                         target_db => $self->o('rnaseq_db'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['rnaseq_sanity_checks'],
                      },

      },

      {
        -logic_name => 'rnaseq_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         skip_check => $self->o('skip_rnaseq'),
                         target_db => $self->o('rnaseq_db'),
                         sanity_check_type => 'gene_db_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'gene_db_checks')->{$self->default_options->{'uniprot_set'}}->{'rnaseq_blast'},
                       },

        -rc_name    => '4GB',
        -flow_into => {
                        1 => ['create_layering_output_db'],
                      },
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
                         script_path => $self->o('clone_db_script_path'),
                         user_r => $self->o('user_r'),
                         user => $self->o('user'),
                         pass_w => $self->o('password'),
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
                         target_db => $self->o('utr_output_db'),
                         create_type => 'clone',
                         script_path => $self->o('clone_db_script_path'),
                         user_r => $self->o('user_r'),
                         user => $self->o('user'),
                         pass_w => $self->o('password'),
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
                         script_path => $self->o('clone_db_script_path'),
                         user_r => $self->o('user_r'),
                         user => $self->o('user'),
                         pass_w => $self->o('password'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '1->A' => ['create_toplevel_slices'],
                        'A->1' => ['layer_annotation_sanity_checks'],
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
                         feature_dbs => [$self->o('genblast_db'),$self->o('projection_realign_db'),$self->o('rnaseq_db')],
                      },
        -flow_into => {
                       '2' => ['layer_annotation'],
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
                         input_gene_dbs => $self->o('utr_gene_dbs'),
                         iid_type => 'slice',
                       },
        -batch_size => 100,
        -hive_capacity => $self->hive_capacity_classes->{'hc_medium'},
        -rc_name    => 'transcript_finalisation',
        -flow_into => {
                        2 => ['cluster_input_genes'],
                      },
      },


      {
        -logic_name => 'cluster_input_genes',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClusterSourceGenes',
        -parameters => {
                         logic_name => 'cluster_input_genes',
                         dna_db => $self->o('dna_db'),
                         input_gene_dbs => $self->o('utr_gene_dbs'),
                         allowed_input_sets => undef,
                         iid_type => 'slice',
                       },
        -batch_size => 10,
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
                         input_gene_dbs => $self->default_options()->{'utr_gene_dbs'},
                         utr_biotype_priorities => $self->o('utr_biotype_priorities'),
                         utr_output_db => $self->default_options()->{'utr_output_db'},
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
                         layering_output_db => $self->o('utr_output_db'),
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
                         script_path => $self->o('clone_db_script_path'),
                         user_r => $self->o('user_r'),
                         user => $self->o('user'),
                         pass_w => $self->o('password'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '1->A' => ['create_pseudogene_ids'],
                        'A->1' => ['create_final_geneset_db'],
                      },
      },

            {
        -logic_name => 'create_pseudogene_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('genebuilder_db'),
                         iid_type => 'feature_id',
                         feature_type => 'gene',
                      },
        -flow_into => {
                       2 => ['pseudogenes'],
                      },
        -rc_name    => 'default',
      },

            {
        -logic_name => 'pseudogenes',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivePseudogenes',
        -parameters => {
                         input_gene_db => $self->o('genebuilder_db'),
                         repeat_db => $self->default_options->{'dna_db'},
                         output_db => $self->o('pseudogene_db'),
                         dna_db => $self->default_options->{'dna_db'},
                         logic_name => 'pseudogenes',
                         module     => 'HivePseudogenes',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::PseudoGeneStatic','pseudogenes')},
                         PS_MULTI_EXON_DIR       => $self->o('output_path').'/pseudogenes/multi_exon_dir/',
                       },
        -batch_size => 20,
        -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
        -rc_name    => 'transcript_finalisation',
        -flow_into => {
                       2 => ['spliced_elsewhere'],
                      },
      },


      {
        -logic_name => 'concat_blast_db',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'for i in '.$self->o('output_path').'/pseudogenes/multi_exon_dir/multi_exon_seq*.fasta;'.
                                'do cat $i >> '.$self->o('output_path').'/pseudogenes/multi_exon_dir/all_multi_exon_genes.fasta;'.
                                'done'
                       },
         -rc_name => 'default',
         -wait_for => ['pseudogenes'],
         -input_ids => [{}],
         -flow_into => { 1 => ['format_blast_db'] },
      },


      {
        -logic_name => 'format_blast_db',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'if [ "'.$self->o('blast_type').'" = "ncbi" ];then makeblastdb -dbtype nucl -in '.$self->o('output_path').'/pseudogenes/multi_exon_dir/all_multi_exon_genes.fasta; else xdformat -n '.$self->o('output_path').'/pseudogenes/multi_exon_dir/all_multi_exon_genes.fasta;fi'
                       },
         -rc_name => 'default',
      },


      {
        -logic_name => 'spliced_elsewhere',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSplicedElsewhere',
        -parameters => {
                         input_gene_db => $self->o('genebuilder_db'),
                         repeat_db => $self->default_options->{'dna_db'},
                         output_db => $self->o('pseudogene_db'),
                         dna_db => $self->default_options->{'dna_db'},
                         logic_name => 'spliced_elsewhere',
                         module     => 'HiveSplicedElsewhere',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::PseudoGeneStatic','pseudogenes')},
                         PS_MULTI_EXON_DIR       => $self->o('output_path').'/pseudogenes/multi_exon_dir/',
                       },
        -rc_name          => 'transcript_finalisation',
        -can_be_empty  => 1,
        -wait_for => ['format_blast_db'],
      },


      {
        -logic_name => 'create_final_geneset_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('pseudogene_db'),
                         target_db => $self->o('final_geneset_db'),
                         create_type => 'copy',
                         script_path => $self->o('clone_db_script_path'),
                         user_r => $self->o('user_r'),
                         user_w => $self->o('user'),
                         pass_w => $self->o('password'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '1' => ['run_cleaner'],
                      },
      },


      {
        -logic_name => 'run_cleaner',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCleanGeneset',
        -parameters => {
                         input_db => $self->o('final_geneset_db'),
                         dna_db => $self->o('dna_db'),
                         output_path => $self->o('output_path').'/clean_genes/',
                         blessed_biotypes => $self->o('cleaning_blessed_biotypes'),
                         flagged_redundancy_coverage_threshold => 95,
                         general_redundancy_coverage_threshold => 95,
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '1' => ['delete_flagged_transcripts'],
                      },
      },


      {
       -logic_name => 'delete_flagged_transcripts',
       -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDeleteTranscripts',
       -parameters => {
                        dbhost => $self->o('final_geneset_db','-host'),
                        dbname => $self->o('final_geneset_db','-dbname'),
                        dbuser => $self->o('user'),
                        dbpass => $self->o('password'),
                        dbport => $self->o('final_geneset_db','-port'),
                        transcript_ids_file => $self->o('output_path').'/clean_genes/transcript_ids_to_remove.txt',
                        delete_transcripts_path => $self->o('enscode_root_dir').'/ensembl-analysis/scripts/genebuild/',
                        delete_genes_path => $self->o('enscode_root_dir').'/ensembl-analysis/scripts/genebuild/',
                        delete_transcripts_script_name => 'delete_transcripts.pl',
                        delete_genes_script_name => 'delete_genes.pl',
                        output_path => $self->o('output_path').'/clean_genes/',
                        output_file_name => 'delete_transcripts.out',
                      },
        -max_retry_count => 0,
        -flow_into => {
                        '1' => ['transfer_ncrnas'],
                      },
     },


     {
       -logic_name => 'transfer_ncrnas',
       -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
       -parameters => {
                        cmd => 'perl '.$self->o('enscode_root_dir').'/ensembl-analysis/scripts/genebuild/copy_genes.pl'.
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
        -flow_into => {
                        '1' => ['delete_duplicate_genes'],
                      },
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
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'mysql '.
                                ' -u'.$self->o('user').
                                ' -p'.$self->o('password').
                                ' -h'.$self->o('reference_db','-host').
                                ' -P'.$self->o('reference_db','-port').
                                ' -D'.$self->o('reference_db','-dbname').
                                ' -e \'UPDATE gene set biotype="protein_coding" where biotype="ensembl";'.
                                ' UPDATE gene set analysis_id=(SELECT analysis_id from analysis where logic_name="ensembl")'.
                                ' WHERE analysis_id in (SELECT analysis_id from analysis where logic_name in ("spliced_elsewhere","pseudogenes","genblast"));'.
                                ' UPDATE transcript join gene using(gene_id) set transcript.biotype=gene.biotype;'.
                                ' UPDATE transcript join gene using(gene_id) set transcript.analysis_id=gene.analysis_id;'.
                                ' INSERT into analysis (created,logic_name,db) values (now(),"other_protein","uniprot");'.
                                ' UPDATE protein_align_feature set analysis_id='.
                                '(SELECT analysis_id from analysis where logic_name="other_protein") where analysis_id='.
                                '(SELECT analysis_id from analysis where logic_name="genblast");\'',
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
                         cmd => 'perl '.$self->o('enscode_root_dir').'/ensembl/misc-scripts/meta_coord/update_meta_coord.pl'.
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
                         cmd => 'perl '.$self->o('enscode_root_dir').'/ensembl/misc-scripts/meta_levels.pl'.
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
                         cmd => 'perl '.$self->o('enscode_root_dir').'/ensembl/misc-scripts/frameshift_transcript_attribs.pl'.
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
                         cmd => 'perl '.$self->o('enscode_root_dir').'/ensembl/misc-scripts/canonical_transcripts/select_canonical_transcripts.pl'.
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
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'mysql '.
                                ' -u'.$self->o('user').
                                ' -p'.$self->o('password').
                                ' -h'.$self->o('reference_db','-host').
                                ' -P'.$self->o('reference_db','-port').
                                ' -D'.$self->o('reference_db','-dbname').
                                ' -e \'UPDATE gene set stable_id=NULL;'.
                                ' UPDATE transcript set stable_id=NULL;'.
                                ' UPDATE translation set stable_id=NULL;'.
                                ' UPDATE exon set stable_id=NULL;'.
                                ' UPDATE protein_align_feature set external_db_id=NULL;'.
                                ' UPDATE dna_align_feature set external_db_id=NULL;\'',
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
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         'source_db'        => $self->o('dna_db'),
                         'user_w'           => $self->o('user'),
                         'pass_w'           => $self->o('password'),
                         'create_type'      => 'backup',
                         'output_path'      => $self->o('output_path')."/core_db_bak/",
                         'backup_name'      => 'core_bak_post_stable_id.sql',
                       },
        -rc_name    => 'default',
        -flow_into => { 1 => ['load_external_db_ids_and_optimise_af_tables'] },
      },


      {
        -logic_name => 'load_external_db_ids_and_optimise_af_tables',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('enscode_root_dir').'/ensembl-analysis/scripts/genebuild/load_external_db_ids_and_optimize_af.pl'.
                                ' -output_path '.$self->o('output_path').'/optimise/'.
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
                                ' -verbose'
                       },
        -max_retry_count => 1,
        -rc_name => '8GB',
        -flow_into => {
                        1 => ['clean_unused_analyses'],
                      },
      },


      {
        -logic_name => 'clean_unused_analyses',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'mysql '.
                                ' -u'.$self->o('user').
                                ' -p'.$self->o('password').
                                ' -h'.$self->o('reference_db','-host').
                                ' -P'.$self->o('reference_db','-port').
                                ' -D'.$self->o('reference_db','-dbname').
                                ' -e \'DELETE from analysis_description where analysis_id in '.
                                '(SELECT analysis_id from analysis where logic_name in ("spliced_elsewhere","pseudogenes","genblast"));'.
                                ' DELETE from analysis where logic_name in ("spliced_elsewhere","pseudogenes","genblast");\'',
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['update_rnaseq_analyses'],
                      },
      },


      {
        -logic_name => 'update_rnaseq_analyses',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'mysql '.
                                ' -u'.$self->o('user').
                                ' -p'.$self->o('password').
                                ' -h'.$self->o('reference_db','-host').
                                ' -P'.$self->o('reference_db','-port').
                                ' -D'.$self->o('reference_db','-dbname').
                                ' -e \'UPDATE analysis set logic_name=replace(logic_name,"_rnaseq","_introns")'.
                                ' WHERE logic_name like "'.$self->o('species_name').'%_rnaseq";\'',
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['drop_backup_tables'],
                      },
      },


      {
        -logic_name => 'drop_backup_tables',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'mysql '.
                                ' -u'.$self->o('user').
                                ' -p'.$self->o('password').
                                ' -h'.$self->o('reference_db','-host').
                                ' -P'.$self->o('reference_db','-port').
                                ' -D'.$self->o('reference_db','-dbname').
                                ' -e \'show tables like "%bak%";\' '.
                                '| xargs -I "@@" echo '.
                                'mysql '.
                                ' -u'.$self->o('user').
                                ' -p'.$self->o('password').
                                ' -h'.$self->o('reference_db','-host').
                                ' -P'.$self->o('reference_db','-port').
                                ' -D'.$self->o('reference_db','-dbname').
                                ' -e \'DROP TABLE @@";\'',
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['core_gene_set_sanity_checks'],
                      },
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
                        1 => ['core_general_sanity_checks'],
                      },
      },


      {
        -logic_name => 'core_general_sanity_checks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
        -parameters => {
                         skip_rnaseq => $self->o('skip_rnaseq'),
                         skip_projection => $self->o('skip_projection'),
                         target_db => $self->o('reference_db'),
                         sanity_check_type => 'final_core_checks',
                         min_allowed_feature_counts => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
                                                                             'final_core_checks')->{$self->default_options->{'uniprot_set'}},
                       },

        -rc_name    => '4GB',
        -flow_into => {
                        1 => ['run_healthchecks'],
                      },
      },


      {
        -logic_name => 'run_healthchecks',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveHealthcheck',
        -parameters => {
                         input_db         => $self->o('reference_db'),
                         species          => $self->o('species_name'),
                         group            => 'GenebuildHandover',
                         enscode_root_dir => $self->o('enscode_root_dir'),
                       },
        -rc_name    => 'default',
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
    'simple_features' => { LSF => $self->lsf_resource_builder('production-rh7', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]).' -W 20'},
    'genscan' => { LSF => $self->lsf_resource_builder('production-rh7', 3900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]).' -W 120'},
    'refseq_cdna' => { LSF => $self->lsf_resource_builder('production-rh7', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]).' -W 120'},
    'refseq_cdna_retry' => { LSF => $self->lsf_resource_builder('production-rh7', 5900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]).' -W 120'},
    'genscan_short' => { LSF => $self->lsf_resource_builder('production-rh7', 5900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]).' -W 60'},
    'blast' => { LSF => $self->lsf_resource_builder('production-rh7', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], undef, 3)},
    'blast_retry' => { LSF => $self->lsf_resource_builder('production-rh7', 5900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], undef, 3)},
    'rfam_blast' => { LSF => $self->lsf_resource_builder('production-rh7', 4000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], undef, 3)},
    'rfam_blast_retry' => { LSF => $self->lsf_resource_builder('production-rh7', 6000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], undef, 3)},
    'genblast' => { LSF => $self->lsf_resource_builder('production-rh7', 1900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genblast_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'genblast_retry' => { LSF => $self->lsf_resource_builder('production-rh7', 4900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genblast_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'project_transcripts' => { LSF => $self->lsf_resource_builder('production-rh7', 1900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'projection_db_server'}, $self->default_options->{'projection_lastz_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]).' -W 120'},
    'refseq_import' => { LSF => $self->lsf_resource_builder('production-rh7', 9900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'refseq_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'layer_annotation' => { LSF => $self->lsf_resource_builder('production-rh7', 3900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genblast_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]).' -W 200'},
    'genebuilder' => { LSF => $self->lsf_resource_builder('production-rh7', 1900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genblast_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]).' -W 200'},
    'transcript_finalisation' => { LSF => $self->lsf_resource_builder('production-rh7', 1900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genblast_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'filter' => { LSF => $self->lsf_resource_builder('production-rh7', 4900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genblast_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
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
