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
    'email_address'             => '', #!!!!!!!!!!!
    'genebuilder_id'            => '', #!!!!!!!!!!!
    'release_number'            => '', #!!!!!!!!!!! (you can put whatever here if you're not sure)
    'enscode_root_dir'          => '', #!!!!!!!!!!! git repo checkouts
    'repeatmasker_library'      => '', #!!!!!!!!!!!!!!!!! repbase library name e.g. rodents
    'repeatmasker_species'      => '', #!!!!!!!!!!!!!!!!! repbase library name e.g. rodents, without '_'
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

########################
# Pipe and ref db info
########################
    'pipeline_name'                => '', #!!!!!!!!!!! What you want hive to call the pipeline, not the db name itself
    'user_r'                       => '', #!!!!!!!!!!!
    'user'                         => '', #!!!!!!!!!!!
    'password'                     => '', #!!!!!!!!!!!
    'port'                         => '', #!!!!!!!!!!!

    'pipe_db_server'               => '',
    'pipe_db_port'                 => '',

    'dna_dbname'                   => $self->o('dbowner').'_'.$self->o('production_name').'_core_'.$self->o('release_number'),
    'dna_db_server'                => '',
    'dna_db_port'                  => $self->o('port'),

    'reference_dbname'             => $self->o('dna_dbname'),
    'reference_db_server'          => '',
    'reference_db_port'            => $self->o('dna_db_port'),

    'cdna_db_server'               => '',
    'cdna_db_port'                 => $self->o('port'),

    'genblast_db_server'           => '',
    'genblast_db_port'             => $self->o('port'),

    'genewise_db_server'           => '',
    'genewise_db_port'             => $self->o('port'),

    'projection_db_server'         => '',
    'projection_db_port'           => $self->o('port'),

    'projection_realign_db_server' => '',
    'projection_realign_db_port'   => $self->o('port'),

    'projection_source_db_name'    => '', # This is generally a pre-existing db, like the current human/mouse core for example
    'projection_source_db_server'  => 'mysql-ensembl-mirror',
    'projection_source_db_port'    => 4240,

    'projection_lastz_db_server'   => '',
    'projection_lastz_db_port'     => $self->o('port'),

    'rnaseq_dbname'                => $self->o('dbowner').'_'.$self->o('species_name').'_rnaseq_'.$self->o('release_number'),
    'rnaseq_db_server'             => '',
    'rnaseq_db_port'               => $self->o('port'),

    'layering_db_server'           => '',
    'layering_db_port'             => $self->o('port'),

    'utr_db_server'                => '',
    'utr_db_port'                  => $self->o('port'),

    'genebuilder_db_server'        => '',
    'genebuilder_db_port'          => $self->o('port'),

    'pseudogene_db_server'         => '',
    'pseudogene_db_port'           => $self->o('port'),

    'ncrna_db_server'              => '',
    'ncrna_db_port'                => $self->o('port'),

    'final_geneset_db_server'      => '',
    'final_geneset_db_port'        => $self->o('port'),

    'refseq_db_server'             => '',
    'refseq_db_port'               => $self->o('port'),

    'killlist_db_server'           => '',
    'killlist_db_port'             => $self->o('port'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'staging_1_db_server'          => 'mysql-ens-sta-1',
    'staging_1_port'               => 4519,


    databases_to_delete => ['reference_db', 'cdna_db', 'genblast_db', 'genewise_db', 'projection_db', 'projection_realign_db', 'layering_db', 'utr_output_db', 'genebuilder_db', 'pseudogene_db', 'ncrna_db', 'final_geneset_db', 'refseq_db'],


########################
# BLAST db paths
########################
    'uniprot_blast_db_path'     => '/nfs/production/panda/ensembl/genebuild/blastdb/uniprot/uniprot_2016_10/uniprot_vertebrate',
    'vertrna_blast_db_path'     => '/nfs/production/panda/ensembl/genebuild/blastdb/vertrna/129/embl_vertrna-1',
    'unigene_blast_db_path'     => '/nfs/production/panda/ensembl/genebuild/blastdb/unigene/unigene',
    'mito_index_path'           => '/nfs/production/panda/ensembl/genebuild/blastdb/refseq_mitochondria_set/mito_index.txt',
    'ncrna_blast_path'          => '/nfs/production/panda/ensembl/genebuild/blastdb/ncrna_2016_05/',
    'mirBase_fasta'             => 'mouse_mirnas.fa',

######################################################
#
# Mostly constant settings
#
######################################################

    'genome_file'               => catfile($self->o('output_path'), 'genome_dumps', $self->o('species_name').'_softmasked_toplevel.fa'),
    'primary_assembly_dir_name' => 'Primary_Assembly',
    'refseq_cdna_calculate_coverage_and_pid' => '0',
    'refseq_cdna_table_name'    => 'refseq_sequences',
    'contigs_source'            => 'NCBI',



    'layering_input_gene_dbs' => [
                                   $self->o('genblast_db'),
                                   $self->o('rnaseq_db'),
                                   $self->o('projection_realign_db'),
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


    'min_toplevel_slice_length'   => 0,

    'repeat_logic_names'          => ['repeatmasker_repbase_'.$self->o('repeatmasker_library'),'dust'],
    'homology_models_path'        => catdir($self->o('output_path'),'homology_models'),
    'clone_db_script_path'        => catfile($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'clone_database.ksh'),
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
    'dust_path' => catfile($self->o('binary_base'), 'tcdust'),
    'trf_path' => catfile($self->o('binary_base'), 'trf'),
    'eponine_java_path' => '/nfs/software/ensembl/RHEL7/jenv/shims/java',
    'eponine_jar_path' => '/nfs/software/ensembl/RHEL7/linuxbrew/Cellar/eponine/1.0/libexec/eponine-scan.jar',
    'cpg_path' => catfile($self->o('binary_base'), 'cpg_lh'),
    'trnascan_path' => catfile($self->o('binary_base'), 'tRNAscan-SE'),
    'repeatmasker_path' => catfile($self->o('binary_base'), 'RepeatMasker'),
    'genscan_path' => catfile($self->o('binary_base'), 'genscan'),
    'genscan_matrix_path' => '/nfs/software/ensembl/RHEL7/linuxbrew/share/HumanIso.smat',
    'uniprot_blast_exe_path' => catfile($self->o('binary_base'), 'blastp'),
    'blastn_exe_path' => catfile($self->o('binary_base'), 'blastp'),
    'vertrna_blast_exe_path' => catfile($self->o('binary_base'), 'tblastn'),
    'unigene_blast_exe_path' => catfile($self->o('binary_base'), 'tblastn'),
    'exonerate_path'         => '/nfs/software/ensembl/RHEL7/linuxbrew/opt/exonerate09/bin/exonerate',
    'cmsearch_exe_path'    => catfile($self->o('binary_base'), 'cmsearch'),

    'uniprot_index_name'          => 'uniprot_index',
    'uniprot_db_name'             => 'uniprot_db',
    'uniprot_query_dir_name'      => 'uniprot_temp',
    'uniprot_genblast_batch_size' => 5,
    'uniprot_table_name'          => 'uniprot_sequences',

    'genblast_path'     => catfile($self->o('binary_base'), 'genblast'),
    'genblast_eval'     => '1e-20',
    'genblast_cov'      => '0.5',
    'genblast_pid'      => '50',
    'genblast_max_rank' => '5',

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
    'assembly_ftp_path'       => $self->o('insdc_base_ftp').'/#assembly_accession#_#assembly_name#_assembly_structure/',
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

#sub pipeline_create_commands {
#    my ($self) = @_;
#    return [
    # inheriting database and hive tables' creation
#      @{$self->SUPER::pipeline_create_commands},

#      $self->db_cmd('CREATE TABLE '.$self->o('uniprot_table_name').' ('.
#                    'accession varchar(50) NOT NULL,'.
#                    'source_db varchar(50) NOT NULL,'.
#                    'pe_level varchar(50) NOT NULL,'.
#                    'biotype varchar(255) NOT NULL,'.
#                    'group_name varchar(255) NOT NULL,'.
#                    'seq text NOT NULL,'.
#                    'PRIMARY KEY (accession))'),

#      $self->db_cmd('CREATE TABLE '.$self->o('refseq_cdna_table_name').' ('.
#                    'accession varchar(50) NOT NULL,'.
#                    'source_db varchar(50) NOT NULL,'.
#                    'biotype varchar(25) NOT NULL,'.
#                    'date varchar(50) NOT NULL,'.
#                    'seq text NOT NULL,'.
#                    'PRIMARY KEY (accession))'),
#    ];
#}
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
     $self->db_cmd('CREATE TABLE projection_protein_sequences ('.
                   'accession varchar(50) NOT NULL,'.
                   'seq text NOT NULL,'.
                   'PRIMARY KEY (accession))'),

       ];
}

## See diagram for pipeline structure
sub pipeline_analyses {
    my ($self) = @_;

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
#                         'coord_system_version'      => $self->o('assembly_name'),
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
        -rc_name    => 'default',
        -flow_into  => {
                         1 => ['load_meta_info'],
                       },
      },


            {
        # Load some meta info and seq_region_synonyms
        -logic_name => 'load_meta_info',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveSetMetaAndSeqRegionSynonym',
        -parameters => {
                         'taxon_id'                  => $self->o('taxon_id'),
                         'production_name'			=> $self->o('production_name'),
                         'genebuilder_id'            => $self->o('genebuilder_id'),
                         'target_db'                 => $self->o('reference_db'),
                         'taxonomy_db'                 => $self->o('taxonomy_db'),
                         'production_db'                 => $self->o('production_db'),
                         'output_path'               => $self->o('output_path'),
                         'enscode_root_dir'          => $self->o('enscode_root_dir'),
                         'primary_assembly_dir_name' => $self->o('primary_assembly_dir_name'),
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
                          1 => ['load_mitochondrion', 'download_refseq_gff'],
                       },
      },


      {
        # Load the AGP files
        -logic_name => 'load_mitochondrion',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadMitochondrion',
        -parameters => {
                         'target_db'                 => $self->o('reference_db'),
                         'output_path'               => $self->o('output_path'),
                         'enscode_root_dir'          => $self->o('enscode_root_dir'),
                         'mito_index_path'           => $self->o('mito_index_path'),
                         'species_name'              => $self->o('species_name'),
                         'chromosomes_present'       => $self->o('chromosomes_present'), # may need to remove later
                         'mt_accession'              => $self->o('mt_accession'),
                      },
        -rc_name    => 'default',

        -flow_into => {
                        # FirstEF is left out here as it runs on repeats. Repeatmasker analyses flow into it
                        '1->A' => ['create_1mb_slice_ids'],
                        'A->1' => ['backup_core_db'],
                      },

      },



      {
        # Create 1mb toplevel slices, each species flow into this independantly
        -logic_name => 'create_1mb_slice_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db        => $self->o('dna_db'),
                         coord_system_name => 'toplevel',
                         iid_type => 'slice',
                         slice_size => 1000000,
                         include_non_reference => 0,
                         top_level => 1,
                         min_slice_length => $self->o('min_toplevel_slice_length'),
#                         batch_slice_ids => 1,
#                         batch_target_size => 1000000,
                       },
        -flow_into => {
                        '2' => ['run_repeatmasker'],
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
                                      ' -dnadbport '.$self->o('dna_db','-port').
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
        -logic_name => 'run_repeatmasker',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveRepeatMasker',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         logic_name => 'repeatmasker_repbase_'.$self->o('repeatmasker_library'),
                         module => 'HiveRepeatMasker',
                         repeatmasker_path => $self->o('repeatmasker_path'),
                         commandline_params => '-nolow -species "'.$self->o('repeatmasker_species').'" -engine "'.$self->o('repeatmasker_engine').'"',
                       },
        -rc_name    => 'repeatmasker',
        -flow_into => {
                        -1 => ['run_repeatmasker_himem'],
                        1 => ['run_dust'],
                      },
        -hive_capacity => 900,
      },

            {
        -logic_name => 'run_repeatmasker_himem',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveRepeatMasker',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         logic_name => 'repeatmasker_repbase_'.$self->o('repeatmasker_library'),
                         module => 'HiveRepeatMasker',
                         repeatmasker_path => $self->o('repeatmasker_path'),
                         commandline_params => '-nolow -species "'.$self->o('repeatmasker_species').'" -engine "'.$self->o('repeatmasker_engine').'"',
                       },
        -rc_name    => 'repeatmasker_himem',
        -flow_into => {
                         1 => ['run_dust'],
                      },
        -hive_capacity => 900,
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
        -wait_for => ['run_dust','run_repeatmasker','run_repeatmasker_himem'],
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
                      },
        -hive_capacity => 900,
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
                      },
       -hive_capacity => 900,
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
#                         1 => ['run_firstef'],
                         1 => ['run_cpg'],
                      },
       -hive_capacity => 900,
       -batch_size => 20,
      },

#            {
#        # Run FirstEF. This differs slightly from the other in that it uses repeatmasking so in pipeline terms
#        # it is set to run after repeatmasker instead of in parallel, unlike the others
#        -logic_name => 'run_firstef',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveFirstEF',
#        -parameters => {
#                         target_db => $self->o('dna_db'),
#                         logic_name => 'firstef',
#                         module => 'HiveFirstEF',
#                         firstef_path => $self->o('firstef_path'),
#                         repeat_masking_logic_names => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
#                         commandline_params => '-repeatmasked',
#                       },
#        -rc_name    => 'simple_features',
#        -flow_into => {
#                         1 => ['run_cpg'],
#                      },
#       -hive_capacity => 900,
#       -batch_size => 20,
#      },

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
                      },
       -hive_capacity => 900,
       -batch_size => 20,
      },

            {
        # Run tRNAscan
        -logic_name => 'run_trnascan',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveTRNAScan',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         logic_name => 'trnascan',
                         module => 'HiveTRNAScan',
                         trnascan_path => $self->o('trnascan_path'),
                       },
        -rc_name    => 'simple_features',
        -flow_into => {
                         1 => ['run_genscan'],
                      },
       -hive_capacity => 900,
       -batch_size => 20,
      },


###############################################################################
#
# GENSCAN ANALYSIS
#
##############################################################################

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
                         repeat_masking_logic_names => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
                       },
        -rc_name    => 'genscan',
        -flow_into => {
                        1 => ['create_prediction_transcript_ids'],
                        -1 => ['decrease_genscan_slice_size'],
                        -2 => ['decrease_genscan_slice_size'],
                        -3 => ['decrease_genscan_slice_size'],
                      },
       -hive_capacity => 900,
       -batch_size => 20,
      },



            {
        # Create 1mb toplevel slices, each species flow into this independantly
        -logic_name => 'decrease_genscan_slice_size',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         iid_type => 'split_slice',
                         slice_size => 100000,
                       },
        -flow_into => {
                        2 => ['run_genscan_short_slice'],
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
                         repeat_masking_logic_names => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
                       },
        -rc_name    => 'genscan',
        -flow_into => {
                        1 => ['create_prediction_transcript_ids'],
                        -1 => ['failed_genscan_slices'],
                        -2 => ['failed_genscan_slices'],
                        -3 => ['decrease_genscan_slice_size'],
                      },
        -rc_name    => 'genscan_short',
        -can_be_empty  => 1,
        -hive_capacity => 900,
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
        # Create input ids for individual prediction transcripts. Takes a slice as an input id and converts it
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
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlastGenscanPep',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         blast_db_path => $self->o('uniprot_blast_db_path'),
                         blast_exe_path => $self->o('uniprot_blast_exe_path'),
                         commandline_params => $commandline_params{$self->o('blast_type')},
                         repeat_masking_logic_names => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
                         prediction_transcript_logic_names => ['genscan'],
                         iid_type => 'feature_id',
                         logic_name => 'uniprot',
                         module => 'HiveBlastGenscanPep',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastGenscanPep', {BLAST_PARAMS => {type => $self->o('blast_type')}})},
                      },
        -flow_into => {
                        1 => ['run_vertrna_blast'],
                        -1 => ['failed_blast_job'],
                        -2 => ['failed_blast_job'],
                      },
        -rc_name    => 'blast',
        -failed_job_tolerance => 0.5,
        -hive_capacity => 900,
        -batch_size => 20,
      },

      {
        # BLAST individual prediction transcripts against vertRNA. The config settings are held lower in this
        # file in the master_config_settings sub
        -logic_name => 'run_vertrna_blast',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlastGenscanDNA',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         blast_db_path => $self->o('vertrna_blast_db_path'),
                         blast_exe_path => $self->o('vertrna_blast_exe_path'),
                         commandline_params => $commandline_params{$self->o('blast_type')},
                         repeat_masking_logic_names => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
                         prediction_transcript_logic_names => ['genscan'],
                         iid_type => 'feature_id',
                         logic_name => 'vertrna',
                         module => 'HiveBlastGenscanDNA',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastGenscanVertRNA', {BLAST_PARAMS => {type => $self->o('blast_type')}})},
                      },
        -flow_into => {
                        1 => ['run_unigene_blast'],
                        -1 => ['failed_blast_job'],
                        -2 => ['failed_blast_job'],
                      },
        -rc_name    => 'blast',
       -failed_job_tolerance => 0.5,
       -hive_capacity => 900,
       -batch_size => 20,
      },


      {
        # BLAST individual prediction transcripts against unigene. The config settings are held lower in this
        # file in the master_config_settings sub
        -logic_name => 'run_unigene_blast',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlastGenscanDNA',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         blast_db_path => $self->o('unigene_blast_db_path'),
                         blast_exe_path => $self->o('unigene_blast_exe_path'),
                         commandline_params => $commandline_params{$self->o('blast_type')},
                         prediction_transcript_logic_names => ['genscan'],
                         iid_type => 'feature_id',
                         repeat_masking_logic_names => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
                         logic_name => 'unigene',
                         module => 'HiveBlastGenscanDNA',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastGenscanUnigene', {BLAST_PARAMS => {type => $self->o('blast_type')}})},
                      },
        -flow_into => {
                        -1 => ['failed_blast_job'],
                        -2 => ['failed_blast_job'],
                      },
        -rc_name    => 'blast',
        -failed_job_tolerance => 0.5,
        -hive_capacity => 900,
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
        # Creates a reference db for each species
        -logic_name => 'backup_core_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         'source_db'        => $self->o('dna_db'),
                         'user'           	=> $self->o('user_r'),
                         'user_w' 			=> $self->o('user'),
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
         -flow_into => { 1 => ['set_meta_coords'] },
      },

     {
        -logic_name => 'set_meta_coords',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('meta_coord_script').
                                ' -user '.$self->o('reference_db', '-user').
                                ' -pass '.$self->o('reference_db', '-pass').
                                ' -host '.$self->o('reference_db','-host').
                                ' -port '.$self->o('reference_db','-port').
                                ' -dbpattern '.$self->o('reference_db','-dbname')
                       },
         -rc_name => 'default',
         -flow_into => { 1 => ['set_meta_levels'] },
      },

      {
        -logic_name => 'set_meta_levels',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('meta_levels_script').
                                ' -user '.$self->o('reference_db', '-user').
                                ' -pass '.$self->o('reference_db', '-pass').
                                ' -host '.$self->o('reference_db','-host').
                                ' -port '.$self->o('reference_db','-port').
                                ' -dbname '.$self->o('reference_db','-dbname')
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
                        '1->A' => ['get_refseq_cdnas'],
                        'A->1' => ['create_genblast_output_db'],
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
                        2 => ['refseq_cdna'],
                      },
        -rc_name    => 'default',
      },

      {
        -logic_name => 'refseq_cdna',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
        -rc_name    => 'genblast',
        -parameters => {
                         iid_type => 'db_seq',
                         query_table_name => $self->o('refseq_cdna_table_name'),
                         sequence_table_name => $self->o('refseq_cdna_table_name'), # there is a problem here. I add that, but it should be query_table_name or sequence_table_name

                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('cdna_db'),
                         logic_name => 'refseq_cdna',
                         module     => 'HiveExonerate2Genes',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_cov_per_sub')},
                         GENOMICSEQS         => $self->o('genome_file'),
                         PROGRAM             => $self->o('exonerate_path'),
                         SOFT_MASKED_REPEATS => ['repeatmasker_repbase_'.$self->o('repeatmasker_library'),'dust'],
                         query_seq_dir => $self->o('output_path'),
                         calculate_coverage_and_pid => $self->o('refseq_cdna_calculate_coverage_and_pid'),
                         exonerate_cdna_pid => $self->o('exonerate_cdna_pid'),
                         exonerate_cdna_cov => $self->o('exonerate_cdna_cov'),
                      },
        -flow_into => {
                        -1 => ['refseq_cdna_retry'],
                        -2 => ['failed_refseq_cdnas'],
                      },
        -batch_size => 20,
        -hive_capacity => 900,
        -failed_job_tolerance => 0.5,
     },

     {
        -logic_name => 'refseq_cdna_retry',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
        -rc_name    => 'genblast_retry',
        -parameters => {
                         iid_type => 'db_seq',
                         query_table_name => $self->o('refseq_cdna_table_name'),
                         sequence_table_name => $self->o('refseq_cdna_table_name'), # there is a problem here. I add that, but it should be query_table_name or sequence_table_name

                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('cdna_db'),
                         logic_name => 'refseq_cdna',
                         module     => 'HiveExonerate2Genes',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_cov_per_sub')},
                         GENOMICSEQS         => $self->o('genome_file'),
                         PROGRAM             => $self->o('exonerate_path'),
                         SOFT_MASKED_REPEATS => ['repeatmasker_repbase_'.$self->o('repeatmasker_library'),'dust'],
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
        -hive_capacity => 900,
        -failed_job_tolerance => 100,
        -can_be_empty => 1,
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
                        1 => ['process_uniprot_files'],
                      },
      },

            {
        -logic_name => 'process_uniprot_files',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProcessUniProtFiles',
        -parameters => {
                         uniprot_db_name => $self->o('uniprot_db_name'),
                         uniprot_index_name => $self->o('uniprot_index_name'),
                         dest_dir   => $self->o('homology_models_path'),
                         killlist_type => 'protein',
                         killlist_db => $self->o('killlist_db'),
                      },
        -rc_name => 'default',
        -flow_into => {
                       '2->A' => ['load_uniprot_seqs'],
                       'A->1' => ['generate_genblast_jobs'],
                      },
      },

      {
        -logic_name => 'load_uniprot_seqs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadSeqsPipeDB',
        -parameters => {
                         uniprot_db_path => catdir($self->o('homology_models_path'), $self->o('uniprot_db_name')),
                         uniprot_index_path => catdir($self->o('homology_models_path'), $self->o('uniprot_index_name')),
                       },
        -rc_name          => 'default',
        -hive_capacity => 100,
      },



      {
        -logic_name => 'generate_genblast_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         iid_type => 'uniprot_accession',
                         uniprot_batch_size => $self->o('uniprot_genblast_batch_size'),
                         uniprot_table_name => $self->o('uniprot_table_name'),
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
                         commandline_params => ' -P blast -gff -e '.$self->o('genblast_eval').' -c '.$self->o('genblast_cov').' ',
                         query_seq_dir => catdir($self->o('homology_models_path'), $self->o('uniprot_query_dir_name')),
                         sequence_table_name => $self->o('uniprot_table_name'),
                         max_rank => $self->o('genblast_max_rank'),
                         genblast_pid => $self->o('genblast_pid'),
                       },
        -rc_name    => 'genblast',
        -flow_into => {
                        -1 => ['split_genblast_jobs'],
                        -2 => ['split_genblast_jobs'],
                        -3 => ['split_genblast_jobs'],
                      },
        -failed_job_tolerance => 0.5,
        -hive_capacity => 900,
      },

      {
        -logic_name => 'split_genblast_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         iid_type => 'rechunk_uniprot_accession',
                         uniprot_batch_size => 1,
                         uniprot_table_name => $self->o('uniprot_table_name'),
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
                         commandline_params => ' -P wublast -gff -e '.$self->o('genblast_eval').' -c '.$self->o('genblast_cov').' ',
                         query_seq_dir => catdir($self->o('homology_models_path'), $self->o('uniprot_query_dir_name')),
                         sequence_table_name => $self->o('uniprot_table_name'),
                         max_rank => $self->o('genblast_max_rank'),
                         genblast_pid => $self->o('genblast_pid'),
                       },
        -rc_name          => 'genblast_retry',
        -can_be_empty  => 1,
        -failed_job_tolerance => 100,
        -flow_into => {
                        -1 => ['failed_genblast_proteins'],
                        -2 => ['failed_genblast_proteins'],
                        -3 => ['failed_genblast_proteins'],
                      },
        -hive_capacity => 900,
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
                        '1->A' => ['create_200kb_slice_ids'],
                        'A->1' => ['create_projection_output_db'],
                      },

      },

#      {
#        -logic_name => 'ncrna_semaphore',
#        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
#        -parameters => {
#                       },
#        -rc_name          => 'default',
#        -flow_into => {
#                        '1->A' => ['create_200kb_slice_ids'],
#                        'A->1' => ['filter_ncrnas'],
#                      },
#      },


            {
        -logic_name => 'create_200kb_slice_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         coord_system_name => 'toplevel',
                         iid_type => 'slice',
                         slice_size => 200000,
                         top_level => 1,
                         target_db => $self->o('dna_db'),
                         batch_slice_ids => 1,
                         batch_target_size => 600000,
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
#                         commandline_params => ' W=12 B=10000 V=10000 -hspmax 0 -gspmax 0 -kap ',
                         commandline_params => ' -word_size 12 -num_alignments 10000  -num_descriptions 10000 -max_hsps 0 ',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastRFam', {BLAST_PARAMS => {type => $self->o('blast_type')}})},
                       },
       -hive_capacity => 900,
       -rc_name    => 'blast',
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
                         blast_db_path => $self->o('ncrna_blast_path') . '/' . $self->o('mirBase_fasta') ,
                         blast_exe_path => $self->o('blastn_exe_path'),
                         commandline_params => ' ',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastmiRBase', {BLAST_PARAMS => {type => $self->o('blast_type')}})},
                       },
        -hive_capacity => 900,
        -rc_name    => 'blast',
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
        -hive_capacity => 900,
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
        -hive_capacity => 900,
        -rc_name    => 'transcript_finalisation',
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
                         target_db => $self->o('projection_source_db'),
                         iid_type => 'feature_id',
                         feature_type => 'transcript',
                         feature_restriction => 'protein_coding',
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
        -hive_capacity => 900,
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
                         classification_type => 'standard',
                         update_gene_biotype => 1,
                         target_db => $self->o('projection_db'),
                       },
        -rc_name    => 'default',
        -flow_into => {
        	           1=>['create_realign_db'],
           
                      },
      },

      {
        -logic_name => 'create_realign_db',
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
                       'A->1' => ['update_realigned_biotypes'],
                      },      
      },


      {
        -logic_name => 'create_ids_for_evaluate_projection',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
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
        -logic_name => 'update_realigned_biotypes',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                                'cmd'   => 'mysql -NB -u'.$self->o('user').
                                               ' -p'.$self->o('password').
                                               ' -h'.$self->o('projection_realign_db','-host').
                                               ' -P'.$self->o('projection_realign_db','-port').
                                               ' -D'.$self->o('projection_realign_db','-dbname').
                                               ' -e "UPDATE gene SET biotype =\'realign\' ;  UPDATE transcript SET biotype =\'realign\'; "  '                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['classify_realigned_models'],
                      },
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
                        1 => ['classify_rnaseq_models'],
                      },
      },

############################################################################
#
# RNA-seq analyses
#
############################################################################


      {
        -logic_name => 'classify_rnaseq_models',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport',
        -parameters => {
                         update_gene_biotype => 1,
                         classification_type => 'standard',
                         target_db => $self->o('rnaseq_db'),
                       },
        -rc_name    => 'default',
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
                        'A->1' => ['create_pseudogene_db'],
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
                         feature_dbs => [$self->o('genblast_db'),$self->o('projection_db'),$self->o('rnaseq_db')],
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
                         LAYERS => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::LayerAnnotationStatic', $self->default_options->{'uniprot_set'}, undef, 'ARRAY') ,
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
        -hive_capacity => 200,
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
        -batch_size => 20,
        -hive_capacity => 200,
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
        -hive_capacity => 900,
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
                           'input_db' => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::GenebuilderStatic',$self->default_options->{'uniprot_set'}, undef, 'ARRAY'),
                         },
                         OUTPUT_BIOTYPE => 'ensembl',
                         MAX_TRANSCRIPTS_PER_CLUSTER => 10,
                         MIN_SHORT_INTRON_LEN => 7, #introns shorter than this seem
                         #to be real frame shifts and shoudn't be ignored
                         MAX_SHORT_INTRON_LEN => 15,
                         BLESSED_BIOTYPES => {
                                              'ccds_gene' => 1,
                                              'Blessed_UTR_Genes' => 1,
                                             },
                         #the biotypes of the best genes always to be kept
                         MAX_EXON_LENGTH => 20000,
                         #if the coding_only flag is set to 1, the transcript clustering into genes is done over coding exons only
                         # the current standard way is to cluster only on coding exons
                         CODING_ONLY => 1,
                       },
        -rc_name    => 'transcript_finalisation',
        -hive_capacity => 900,
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
                         repeat_db => $self->o('dna_db'),
                         output_db => $self->o('pseudogene_db'),
                         dna_db => $self->o('dna_db'),
                         logic_name => 'pseudogenes',
                         module     => 'HivePseudogenes',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::PseudoGeneStatic','pseudogenes')},
                         PS_MULTI_EXON_DIR       => $self->o('output_path').'/pseudogenes/multi_exon_dir/',
                       },
        -batch_size => 20,
        -hive_capacity => 900,
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
                         cmd => 'xdformat -n '.$self->o('output_path').'/pseudogenes/multi_exon_dir/all_multi_exon_genes.fasta'
                       },
         -rc_name => 'default',
      },

            {
        -logic_name => 'spliced_elsewhere',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSplicedElsewhere',
        -parameters => {
                         input_gene_db => $self->o('genebuilder_db'),
                         repeat_db => $self->o('dna_db'),
                         output_db => $self->o('pseudogene_db'),
                         dna_db => $self->o('dna_db'),
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
                         user => $self->o('user'),
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
      },



    ];
}


sub resource_classes {
  my $self = shift;

  return {
    'default' => { LSF => $self->lsf_resource_builder('production-rh7', 900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'default_himem' => { LSF => $self->lsf_resource_builder('production-rh7', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'repeatmasker' => { LSF => $self->lsf_resource_builder('production-rh7', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'repeatmasker_himem' => { LSF => $self->lsf_resource_builder('production-rh7', 5900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'simple_features' => { LSF => $self->lsf_resource_builder('production-rh7', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'genscan' => { LSF => $self->lsf_resource_builder('production-rh7', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]).' -W 180'},
    'genscan_short' => { LSF => $self->lsf_resource_builder('production-rh7', 5900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]).' -W 60'},
    'blast' => { LSF => $self->lsf_resource_builder('production-rh7', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], undef, 3)},
    'blast_himem' => { LSF => $self->lsf_resource_builder('production-rh7', 5900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], undef, 3)},
    'blast' => { LSF => $self->lsf_resource_builder('production-rh7', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], undef, 3)},
    'genblast' => { LSF => $self->lsf_resource_builder('production-rh7', 1900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genblast_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]).' -W 120'},
    'genblast_retry' => { LSF => $self->lsf_resource_builder('production-rh7', 4900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genblast_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]).' -W 120'},
    'genewise' => { LSF => $self->lsf_resource_builder('production-rh7', 3900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genewise_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]).' -W 120'},
    'genewise_retry' => { LSF => $self->lsf_resource_builder('production-rh7', 5900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genewise_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]).' -W 120'},
    'project_transcripts' => { LSF => $self->lsf_resource_builder('production-rh7', 1900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'projection_db_server'}, $self->default_options->{'projection_lastz_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]).' -W 120'},
    'refseq_import' => { LSF => $self->lsf_resource_builder('production-rh7', 9900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'refseq_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'layer_annotation' => { LSF => $self->lsf_resource_builder('production-rh7', 3900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genblast_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]).' -W 200'},
    'genebuilder' => { LSF => $self->lsf_resource_builder('production-rh7', 1900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genblast_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]).' -W 200'},
    'transcript_finalisation' => { LSF => $self->lsf_resource_builder('production-rh7', 1900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genblast_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'filter' => { LSF => $self->lsf_resource_builder('production-rh7', 4900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'genblast_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
  }
}

1;
