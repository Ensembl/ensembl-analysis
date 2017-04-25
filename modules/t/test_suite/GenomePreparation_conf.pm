=head1 LICENSE

Copyright [2017] The EMBL-European Bioinformatics Institute

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

package GenomePreparation_conf;

use strict;
use warnings;
use feature 'say';

use File::Spec::Functions;

use Bio::EnsEMBL::ApiVersion qw/software_version/;
use base ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

sub default_options {
  my ($self) = @_;
  return {
    # inherit other stuff from the base class
    %{ $self->SUPER::default_options() },

###########################################
# User parameters: you fill in stuff here # 
###########################################

# Assembly directory
#'assembly_data_directory' => $ENV{PWD}.'/data/assembly', # Fill in!!!!

# Meta directory
#'meta_data_directory'     => $ENV{PWD}.'/data/meta_data', # Fill in!!!!

# Repeatmasker parameters
'repeatmasker_library_path'       => $ENV{PWD}.'/data/repeats/rodent_repeats.fa', # Fill in!!!!

#    'genebuilder_id'            => '' || $ENV{GENEBUILDER_ID}, #!!!!!!!!!!!
    'release_number'            => '', #!!!!!!!!!!! (you can put whatever here if you're not sure)
    'enscode_root_dir'          => $ENV{PWD}.'/../../../..', #!!!!!!!!!!! git repo checkouts
    'repeatmasker_library'      => 'rat', #!!!!!!!!!!!!!!!!! repbase library name e.g. rodents
    'repeatmasker_species'      => 'rat', #!!!!!!!!!!!!!!!!! repbase library name e.g. rodents, without '_'
    'species_name'              => 'rattus_norvegicus', #!!!!!!!!!!!!!!!!! e.g. mus_musculus
    'production_name'           => 'rattus_norvegicus',
#    'taxon_id'                  => '', #!!!!!!!!!!!!!!!!! should be in the assembly report file
    'output_path'               => $ENV{PWD}.'/data/assembly', # Lustre output dir
    'assembly_name'             => 'Rnor_workshop', #!!!!!!!!!!!!!! Name (as it appears in the assembly report file)
    'assembly_accession'        => '', #!!!!!!!!!!!!!! GCA
    'assembly_refseq_accession' => '', #!!!!!!!!!!!!!! GCF
    'mt_accession'              => undef, # This should be set to undef unless you know what you are doing

########################
# Pipe and ref db info
########################
    'pipeline_name'                => '', #!!!!!!!!!!! What you want hive to call the pipeline, not the db name itself
    'user_r'                       => $ENV{EHIVE_ROUSER}, #!!!!!!!!!!!
    'databases_port'               => '', #!!!!!!!!!!!
    'databases_server'             => '', #!!!!!!!!!!!

    'pipe_db_server'               => $self->o('host'),
    'pipe_db_port'                  => $self->o('port'), #!!!!!!!!!!!
    'dna_dbname'                   => $self->o('dbowner').'_'.$self->o('production_name').'_core_'.$self->o('release_number'),
    'dna_db_server'                => $self->o('host'), #!!!!!!!!!!!
    'dna_db_port'                  => $self->o('port'), #!!!!!!!!!!!

    'reference_dbname'             => $self->o('dna_dbname'),
    'reference_db_server'          => $self->o('dna_db_server'),
    'reference_db_port'            => $self->o('dna_db_port'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'staging_1_db_server'          => 'mysql-ens-sta-1',
    'staging_1_port'               => 4519,


    databases_to_delete => ['reference_db'],


#########################
## BLAST db paths
#########################
#    'base_blast_db_path'        => '' || $ENV{BLASTDB_DIR},
#    'uniprot_blast_db_path'     => catfile($self->o('base_blast_db_path'), 'uniprot', 'uniprot_2016_10/uniprot_vertebrate'),
#    'vertrna_blast_db_path'     => catfile($self->o('base_blast_db_path'), 'vertrna', '130/embl_vertrna-1'),
#    'unigene_blast_db_path'     => catfile($self->o('base_blast_db_path'), 'unigene', 'unigene'),
    'mito_index_path'           => undef, # Set this path only if you don't want to use the GCF report and if you haven't set 'mt_accession' '/nfs/production/panda/ensembl/genebuild/blastdb/refseq_mitochondria_set/mito_index.txt',
#    'ncrna_blast_path'          => catfile($self->o('base_blast_db_path'), 'ncrna_2016_05'),
#    'mirBase_fasta'             => 'mouse_mirnas.fa',

######################################################
#
# Mostly constant settings
#
######################################################

    'genome_file'               => catfile($self->o('output_path'), 'genome_dumps', $self->o('species_name').'_softmasked_toplevel.fa'),
    'primary_assembly_dir_name' => undef,
    'refseq_cdna_calculate_coverage_and_pid' => '0',
    'refseq_cdna_table_name'    => 'refseq_sequences',
    'contigs_source'            => 'NCBI',

    'min_toplevel_slice_length'   => 0,

    'repeat_logic_names'          => ['repeatmasker_repbase_'.$self->o('repeatmasker_library'),'dust'],

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
    'software_base_path' => '' || $ENV{LINUXBREW_HOME},
    'binary_base' => catdir($self->o('software_base_path'), 'bin'),
    'dust_path' => catfile($self->o('binary_base'), 'dustmasker'),
    'trf_path' => catfile($self->o('binary_base'), 'trf'),
    'eponine_java_path' => '/nfs/software/ensembl/RHEL7/jenv/shims/java',
    'eponine_jar_path' => catfile($self->o('software_base_path'), 'Cellar', 'eponine', '1.0', 'libexec', 'eponine-scan.jar'),
    'cpg_path' => catfile($self->o('binary_base'), 'cpg_lh'),
    'trnascan_path' => catfile($self->o('binary_base'), 'tRNAscan-SE'),
#    'repeatmasker_path' => catfile($self->o('binary_base'), 'RepeatMasker'),
    'repeatmasker_path' => catfile('/nfs/production/panda/ensembl/thibaut/linuxbrew/bin', 'RepeatMasker'),
    'genscan_path' => catfile($self->o('binary_base'), 'genscan'),
    'genscan_matrix_path' => catfile($self->o('software_base_path'), 'share', 'HumanIso.smat'),
    'uniprot_blast_exe_path' => catfile($self->o('binary_base'), 'blastp'),
    'blastn_exe_path' => catfile($self->o('binary_base'), 'blastp'),
    'vertrna_blast_exe_path' => catfile($self->o('binary_base'), 'tblastn'),
    'unigene_blast_exe_path' => catfile($self->o('binary_base'), 'tblastn'),
    'exonerate_path'         => catfile($self->o('software_base_path'), 'opt', 'exonerate09', 'bin', 'exonerate'),
    'cmsearch_exe_path'    => catfile($self->o('binary_base'), 'cmsearch'),

    'uniprot_genblast_batch_size' => 5,
    'uniprot_table_name'          => 'uniprot_sequences',

    'genblast_path'     => catfile($self->o('binary_base'), 'genblast'),
    'genblast_eval'     => $self->o('blast_type') eq 'wu' ? '1e-20' : '1e-1',
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
    'repeatmasker_engine' => 'rmblast',

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# No option below this mark should be modified
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
########################################################
# URLs for retrieving the INSDC contigs and RefSeq files
########################################################
    'ncbi_base_ftp'           => 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all',
    'insdc_base_ftp'          => $self->o('ncbi_base_ftp').'/#expr(substr(#assembly_accession#, 0, 3))expr#/#expr(substr(#assembly_accession#, 4, 3))expr#/#expr(substr(#assembly_accession#, 7, 3))expr#/#expr(substr(#assembly_accession#, 10, 3))expr#/#assembly_accession#_#assembly_name#',
    'assembly_ftp_path'       => $self->o('insdc_base_ftp').'/#assembly_accession#_#assembly_name#_assembly_structure',
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
    'production_dbname' => 'ensembl_production',

    'reference_db' => {
      -dbname => $self->o('reference_dbname'),
      -host   => $self->o('reference_db_server'),
      -port   => $self->o('reference_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'production_db' => {
      -host   => $self->o('staging_1_db_server'),
      -port   => $self->o('staging_1_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -dbname => $self->o('production_dbname'),
      -driver => $self->o('hive_driver'),
    },

#    'taxonomy_db' => {
#      -host   => $self->o('staging_1_db_server'),
#      -port   => $self->o('staging_1_port'),
#      -user   => $self->o('user_r'),
#      -pass   => $self->o('password_r'),
#      -dbname => 'ncbi_taxonomy',
#      -driver => $self->o('hive_driver'),
#    },
    };
}

# This will run commands when beginning the pipeline. In this case we are not running any
# extra commands, so we have not added anything in here
sub pipeline_create_commands {
    my ($self) = @_;
    return [
    # inheriting database and hive tables' creation
      @{$self->SUPER::pipeline_create_commands},
    ];
}


# These are the analyses that will be run by the pipeline. Each analysis is an anonymous hash reference
# Each analysis can be uniquely identified by its logic_name. The general structure is as follows:
#
#
#
#

sub pipeline_analyses {
    my ($self) = @_;

    return [

## my config
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
        -input_ids => [{}],
        -max_retry_count => 0,
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
        -max_retry_count => 0,
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
                         chromosome_scaffold => 'chromosome_scaffold.agp',
                         chromosome_contig => 'chromosome_contig.agp',
                         scaffold_contig => 'scaffold_contig.agp',
                       },
        -rc_name    => 'default',
        -max_retry_count => 0,
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
                         chromosome_scaffold => 'chromosome_scaffold.agp',
                         chromosome_contig => 'chromosome_contig.agp',
                         scaffold_contig => 'scaffold_contig.agp',
                       },

        -rc_name    => 'default',
        -max_retry_count => 0,
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
        -max_retry_count => 0,
        -flow_into  => {
#                         1 => ['load_meta_info'],
                         1 => ['load_mitochondrion'],
                       },
      },


#            {
#              -logic_name => 'load_meta_info',
#              -module     => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
#              -parameters => {
#                               db_conn => $self->o('core_db'),
#                               input_file => $self->o('meta_data_directory')."/meta_data.sql",
#              },
#              -flow_into => {
#                              '1' => ['create_200kb_slice_ids'],
#              },
#
#              -rc_name    => 'local',
#            },

#            {
#        # Load some meta info and seq_region_synonyms
#        -logic_name => 'load_meta_info',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveSetMetaAndSeqRegionSynonym',
#        -parameters => {
#                         'taxon_id'                  => $self->o('taxon_id'),
#                         'genebuilder_id'            => $self->o('genebuilder_id'),
#                         'target_db'                 => $self->o('reference_db'),
#                         'output_path'               => $self->o('output_path'),
#                         'enscode_root_dir'          => $self->o('enscode_root_dir'),
#                         'primary_assembly_dir_name' => $self->o('primary_assembly_dir_name'),
#                         'production_name'           => $self->o('production_name'),
#                       },
#        -rc_name    => 'default',
#        -flow_into  => {
#                          1 => ['load_taxonomy_info'],
#                       },
#      },
#      {
#        -logic_name => 'load_taxonomy_info',
#        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadTaxonomyInfo',
#        -parameters => {
#                         'target_db'        => $self->o('reference_db'),
#                         'enscode_root_dir' => $self->o('enscode_root_dir'),
#                         'production_db'    => $self->o('production_db'),
#                         'taxonomy_db'      => $self->o('taxonomy_db'),
#                       },
#        -rc_name    => 'default',
#
#        -flow_into  => {
#                          1 => ['load_mitochondrion'],
#                       },
#      },

      {
        -logic_name => 'load_mitochondrion',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadMitochondrion',
        -parameters => {
                         'target_db'                 => $self->o('reference_db'),
                         'output_path'               => $self->o('output_path'),
                         'enscode_root_dir'          => $self->o('enscode_root_dir'),
                         'mito_index_path'           => $self->o('mito_index_path'),
                         'species_name'              => $self->o('species_name'),
                         'mt_accession'              => $self->o('mt_accession'),
                      },
        -rc_name    => 'default',
        -max_retry_count => 0,

        -flow_into => {
                        1 => ['create_toplevel_slice_ids'],
                      },

      },

      {
        -logic_name => 'create_toplevel_slice_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db        => $self->o('dna_db'),
                         coord_system_name => 'toplevel',
                         iid_type => 'slice',
                         include_non_reference => 0,
                         top_level => 1,
                       },
        -flow_into => {
                        2 => ['create_5mb_slice_ids'],
                      },
      },

      {
        # Create 10mb toplevel slices, each species flow into this independantly
        -logic_name => 'create_5mb_slice_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db        => $self->o('dna_db'),
                         coord_system_name => 'toplevel',
                         iid_type => 'split_slice',
                         slice_size => 5000000,
                         include_non_reference => 0,
                         top_level => 1,
                         min_slice_length => $self->o('min_toplevel_slice_length'),
#                         batch_slice_ids => 1,
#                         batch_target_size => 1000000,
                       },
        -flow_into => {
                        '2->A' => ['run_repeatmasker'],
                        'A->1' => ['run_dust'],
                      },

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
                         commandline_params => '-nolow '.(defined($self->o('repeatmasker_library_path')) ? '-lib '.$self->o('repeatmasker_library_path') : '-species "'.$self->o('repeatmasker_species').'"').' -engine "'.$self->o('repeatmasker_engine').'"',
                       },
        -rc_name    => 'repeatmasker',
        -max_retry_count => 0,
        -flow_into => {
                        -1 => ['run_repeatmasker_himem'],
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
                         commandline_params => '-nolow '.(defined($self->o('repeatmasker_library_path')) ? '-lib '.$self->o('repeatmasker_library_path') : '-species "'.$self->o('repeatmasker_species').'"').' -engine "'.$self->o('repeatmasker_engine').'"',
                       },
        -rc_name    => 'repeatmasker_himem',
        -max_retry_count => 0,
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
        -max_retry_count => 0,
      },
      {
        # This should probably be a proper module
        -logic_name => 'format_softmasked_toplevel',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         'cmd'    => 'if [ "'.$self->o('blast_type').'" = "ncbi" ]; then convert2blastmask -in '.catfile($self->o('output_path'), 'genome_dumps', $self->o('species_name')).'_softmasked_toplevel.fa -parse_seqids -masking_algorithm repeatmasker -masking_options "repeatmasker, default" -outfmt maskinfo_asn1_bin -out '.catfile($self->o('output_path'), 'genome_dumps', $self->o('species_name')).'_softmasked_toplevel.fa.asnb;makeblastdb -in '.catfile($self->o('output_path'), 'genome_dumps', $self->o('species_name')).'_softmasked_toplevel.fa -dbtype nucl -parse_seqids -mask_data '.catfile($self->o('output_path'), 'genome_dumps', $self->o('species_name')).'_softmasked_toplevel.fa.asnb -title "'.$self->o('species_name').'"; else xdformat -n '.catfile($self->o('output_path'), 'genome_dumps', $self->o('species_name')).'_softmasked_toplevel.fa;fi',
                       },
        -rc_name    => 'default_himem',
        -max_retry_count => 0,
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
        -max_retry_count => 0,
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
        -max_retry_count => 0,
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
        -max_retry_count => 0,
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
        -max_retry_count => 0,
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
        -max_retry_count => 0,
        -flow_into => {
                         1 => ['create_1mb_slice_ids'],
                      },
       -hive_capacity => 900,
       -batch_size => 20,
      },

      {
        # Create 1mb toplevel slices, each species flow into this independantly
        -logic_name => 'create_1mb_slice_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db        => $self->o('dna_db'),
                         coord_system_name => 'toplevel',
                         iid_type => 'split_slice',
                         slice_size => 1000000,
                         include_non_reference => 0,
                         top_level => 1,
                         min_slice_length => $self->o('min_toplevel_slice_length'),
#                         batch_slice_ids => 1,
#                         batch_target_size => 1000000,
                       },
        -flow_into => {
                        '2' => ['run_genscan'],
                      },

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
        -max_retry_count => 0,
       -hive_capacity => 900,
       -batch_size => 20,
      },

#            {
#              # Create 200kb slices, each species flow into this independantly
#              -logic_name => 'create_200kb_slice_ids',
#              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
#              -parameters => {
#                               target_db => $self->o('core_db'),
#                               coord_system_name => 'toplevel',
#                               iid_type => 'slice',
#                               slice_size => 200000,
#                               include_non_reference => 0,
#                               top_level => 1,
#                             },
#              -flow_into => {
#                              '2' => ['run_repeatmasker'],
#                            },
#
#            },

    ];
}

sub pipeline_wide_parameters {
    my ($self) = @_;

    return {
        %{ $self->SUPER::pipeline_wide_parameters() },  # inherit other stuff from the base class
    };
}


sub resource_classes {
    my $self = shift;
    return {
      'local' => {'LOCAL' => ''},
    'default' => { LSF => $self->lsf_resource_builder('production-rh7', 900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'default_himem' => { LSF => $self->lsf_resource_builder('production-rh7', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'repeatmasker' => { LSF => $self->lsf_resource_builder('production-rh7', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'repeatmasker_himem' => { LSF => $self->lsf_resource_builder('production-rh7', 5900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'simple_features' => { LSF => $self->lsf_resource_builder('production-rh7', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}])},
    'genscan' => { LSF => $self->lsf_resource_builder('production-rh7', 2900, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}], [$self->default_options->{'num_tokens'}]).' -W 180'},
    }
}

1;
