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

package Genome_annotation_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');
use Bio::EnsEMBL::ApiVersion qw/software_version/;

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
'farm_user_name'       => '', # for ref db prefix
'genebuilder_id'       => '', # for meta table
'enscode_root_dir'     => '', # git repo checkouts
'repeatmasker_library' => '', # repbase library to use
'species_name'         => '',
'taxon_id'             => '',
'uniprot_set'          => '',

'input_gene_dbs' => [
                      $self->o('genewise_db'),
                      $self->o('rnaseq_db'),
                      $self->o('projection_db'),
                    ],

########################
# Pipe and ref db info
########################
'pipeline_name'                  => $self->o('farm_user_name').'_'.$self->o('species_name').'_pipe',

'pipe_db_server'                 => '', # NOTE! used to generate tokens in the resource_classes sub below
'reference_db_server'            => '', # NOTE! used to generate tokens in the resource_classes sub below
'dna_db_server'                  => '',
'genblast_db_server'             => '', # NOTE! used to generate tokens in the resource_classes sub below
'genewise_db_server'             => '', # NOTE! used to generate tokens in the resource_classes sub below
'projection_db_server'           => '',
'projection_source_db_name'      => '', # This is generally a pre-existing db, like the current human core for example
'projection_source_db_server'    => '',
'projection_lastz_db_server'     => '',
'rnaseq_db_server'               => '',
'layering_db_server'             => '',
'genebuilder_db_server'          => '',
'refseq_db_server'               => '',
'killlist_db_server'             => '',
'user_r'                         => '',
'user_w'                         => '',
'password'                       => '',
'port'                           => '',


'output_path'               => '',
'genome_file'               => '',

'primary_assembly_dir_name' => 'Primary_Assembly',

'contigs_source'            => '',
'wgs_id'                    => '',
'assembly_name'             => '',
'assembly_accession'        => '',
'ftp_path'                  => '',
'full_ftp_path'             => '',
'refseq_ftp_path'           => '',
'chromosomes_present'       => '',


########################
# BLAST db paths
########################
'uniprot_blast_db_path'     => '',
'vertrna_blast_db_path'     => '',
'unigene_blast_db_path'     => '',
'mito_index_path'           => '/data/blastdb/Ensembl/refseq_mitochondria_set/mito_index.txt',

######################################################
#
# Mostly constant settings
#
######################################################

'min_toplevel_slice_length'   => 0,

'repeat_logic_names'          => ['repeatmasker_repbase_'.$self->o('repeatmasker_library'),'dust'],
'homology_models_path'        => $self->o('output_path').'/homology_models',
'clone_db_script_path'        => $self->o('enscode_root_dir').'/ensembl-analysis/scripts/clone_database.ksh',
'refseq_synonyms_script_path' => $self->o('enscode_root_dir').'/ensembl-pipeline/scripts/refseq_import/load_refseq_synonyms.pl',
'refseq_import_script_path'   => $self->o('enscode_root_dir').'/ensembl-pipeline/scripts/refseq_import/parse_ncbi_gff3.pl',

########################
# Extra db settings
########################

'driver' => 'mysql',
'num_tokens' => 10,

########################
# Executable paths
########################
'dust_path' => '/software/ensembl/genebuild/usrlocalensemblbin/tcdust',
'trf_path' => '/software/ensembl/genebuild/usrlocalensemblbin/trf',
'eponine_path' => '/software/jdk1.6.0_14/bin/java',
'firstef_path' => '/software/ensembl/genebuild/usrlocalensemblbin/firstef',
'cpg_path' => '/software/ensembl/genebuild/usrlocalensemblbin/cpg',
'trnascan_path' => '/software/ensembl/genebuild/usrlocalensemblbin/tRNAscan-SE',
'repeatmasker_path' => '/software/ensembl/bin/RepeatMasker_4_0_5/RepeatMasker',
'genscan_path' => '/software/ensembl/genebuild/usrlocalensemblbin/genscan',
'uniprot_blast_exe_path' => '/software/ensembl/genebuild/usrlocalensemblbin/wublastp',
'vertrna_blast_exe_path' => '/software/ensembl/genebuild/usrlocalensemblbin/wutblastn',
'unigene_blast_exe_path' => '/software/ensembl/genebuild/usrlocalensemblbin/wutblastn',


'uniprot_index_name'         => 'uniprot_index',
'uniprot_db_name'            => 'uniprot_db',
'uniprot_query_dir_name'     => 'uniprot_temp',
'uniprot_genblast_batch_size' => 5,
'uniprot_table_name'         => 'uniprot_sequences',

'genblast_path'              => 'genblast',
'genblast_eval'              => '1e-20',
'genblast_cov'               => '0.5',
'genblast_pid'               => '50',
'genblast_max_rank'          => '5',

# Max internal stops for projected transcripts
'projection_pid'                        => '50',
'projection_cov'                        => '50',
'projection_max_internal_stops'         => '1',
'projection_calculate_coverage_and_pid' => '1',

## Add in genewise path and put in matching code
'genewise_pid'              => '50',
'genewise_cov'              => '50',
'genewise_region_padding'   => '50000',
'genewise_calculate_coverage_and_pid' => '1',

'default_mem'                => '900',
'genblast_mem'               => '1900',
'genblast_retry_mem'         => '4900',
'genewise_mem'               => '3900',
'genewise_retry_mem'         => '5900',
'refseq_mem'                 => '9900',
'projection_mem'             => '1900',
'layer_annotation_mem'       => '3900',
'genebuilder_mem'            => '1900',


########################
# Misc setup info
########################
'repeatmasker_engine' => 'wublast',
'contigs_source' => 'ncbi',
'primary_assembly_dir_name' => 'Primary_Assembly',

########################
# db info
########################
'user' => 'ensro',
'killlist_dbname' => 'gb_kill_list',

'ensembl_release' => $ENV{ENSEMBL_RELEASE},

'pipeline_db' => {
  -dbname => $self->o('farm_user_name').'_'.$self->o('species_name').'_pipe',
  -host   => $self->o('pipe_db_server'),
  -port   => $self->o('port'),
  -user   => $self->o('user_w'),
  -pass   => $self->o('password'),
  -driver => $self->o('driver'),
},

# NOTE! the dbname for each species is generated in the pipeline itself by setup_assembly_loading_pipeline
'reference_db' => {
  -dbname => $self->o('farm_user_name').'_'.$self->o('species_name').'_core',
  -host   => $self->o('reference_db_server'),
  -port   => $self->o('port'),
  -user   => $self->o('user_w'),
  -pass   => $self->o('password'),
},

'dna_db' => {
  -dbname => $self->o('farm_user_name').'_'.$self->o('species_name').'_core',
  -host   => $self->o('reference_db_server'),
  -port   => $self->o('port'),
  -user   => $self->o('user_r'),
},

'genblast_db' => {
  -dbname => $self->o('farm_user_name').'_'.$self->o('species_name').'_genblast',
  -host   => $self->o('genblast_db_server'),
  -port   => $self->o('port'),
  -user   => $self->o('user_w'),
  -pass   => $self->o('password'),
},

'genewise_db' => {
  -dbname => $self->o('farm_user_name').'_'.$self->o('species_name').'_genewise',
  -host   => $self->o('genewise_db_server'),
  -port   => $self->o('port'),
  -user   => $self->o('user_w'),
  -pass   => $self->o('password'),
},

'projection_db' => {
  -dbname => $self->o('farm_user_name').'_'.$self->o('species_name').'_proj',
  -host   => $self->o('projection_db_server'),
  -port   => $self->o('port'),
  -user   => $self->o('user_w'),
  -pass   => $self->o('password'),
},

'projection_source_db' => {
  -dbname => $self->o('projection_source_db_name'),
  -host   => $self->o('projection_source_db_server'),
  -port   => $self->o('port'),
  -user   => $self->o('user_r'),
},

'projection_lastz_db' => {
  -dbname => $self->o('farm_user_name').'_'.$self->o('species_name').'_lastz',
  -host   => $self->o('projection_lastz_db_server'),
  -port   => $self->o('port'),
  -user   => $self->o('user_r'),
},

'rnaseq_db' => {
  -dbname => $self->o('farm_user_name').'_'.$self->o('species_name').'_rnaseq',
  -host   => $self->o('rnaseq_db_server'),
  -port   => $self->o('port'),
  -user   => $self->o('user_w'),
  -pass   => $self->o('password'),
},

'layering_db' => {
  -dbname => $self->o('farm_user_name').'_'.$self->o('species_name').'_layer',
  -host   => $self->o('layering_db_server'),
  -port   => $self->o('port'),
  -user   => $self->o('user_w'),
  -pass   => $self->o('password'),
},

'genebuilder_db' => {
  -dbname => $self->o('farm_user_name').'_'.$self->o('species_name').'_gbuild',
  -host   => $self->o('genebuilder_db_server'),
  -port   => $self->o('port'),
  -user   => $self->o('user_w'),
  -pass   => $self->o('password'),
},

'refseq_db' => {
  -dbname => $self->o('farm_user_name').'_'.$self->o('species_name').'_refseq',
  -host   => $self->o('refseq_db_server'),
  -port   => $self->o('port'),
  -user   => $self->o('user_w'),
  -pass   => $self->o('password'),
},

'killlist_db' => {
                   -dbname    => $self->o('killlist_dbname'),
                   -host      => $self->o('killlist_db_server'),
                   -port      => $self->o('port'),
                   -user      => $self->o('user_r'),
},

'production_db' => {
  -host   => 'ens-staging1',
  -port   => 3306,
  -user   => 'ensro',
  -dbname => 'ensembl_production_'.$self->o('ensembl_release'),
},

'taxonomy_db' => {
  -host   => 'ens-livemirror',
  -port   => 3306,
  -user   => 'ensro',
  -dbname => 'ncbi_taxonomy',
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
    ];
}


## See diagram for pipeline structure
sub pipeline_analyses {
    my ($self) = @_;

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
                         'full_ftp_path'             => $self->o('full_ftp_path'),
                         'output_path'               => $self->o('output_path'),
                         'primary_assembly_dir_name' => $self->o('primary_assembly_dir_name'),
                       },
        -rc_name    => 'default',
        -flow_into  => {
                         1 => ['find_contig_accessions'],
                       },
        -input_ids  => [{}],
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
                         1 => ['populate_production_tables'],
                       },
      },

            {
        # Creates a reference db for each species
        -logic_name => 'create_core_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         'target_db'        => $self->o('reference_db'),
                         'user_w'           => $self->o('user_w'),
                         'pass_w'           => $self->o('password'),
                         'enscode_root_dir' => $self->o('enscode_root_dir'),
                         'create_type'      => 'core_only',
                       },
        -rc_name    => 'default',
        -input_ids => [{}],

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
         -wait_for => ['create_core_db'],
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
                         'chromosomes_present'       => $self->o('chromosomes_present'),
                         'genebuilder_id'            => $self->o('genebuilder_id'),
                         'target_db'                 => $self->o('reference_db'),
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
                       },
        -rc_name    => 'default',

        -flow_into  => {
                          1 => ['load_mitochondrion'],
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
                         'chromosomes_present'       => $self->o('chromosomes_present'),
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
                         target_db        => $self->o('reference_db'),
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
                         commandline_params => '-nolow -species "'.$self->o('repeatmasker_library').'" -engine "'.$self->o('repeatmasker_engine').'"',
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
                         commandline_params => '-nolow -species "'.$self->o('repeatmasker_library').'" -engine "'.$self->o('repeatmasker_engine').'"',
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
                         'species_name'         => $self->o('species_name')."_".$self->o('strain'),
                         'repeat_logic_names'   => $self->o('repeat_logic_names'),
                       },
        -input_ids => [{}],
        -wait_for => ['run_dust','run_repeatmasker','run_repeatmasker_himem'],
        -rc_name    => 'default',
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
                         eponine_path => $self->o('eponine_path'),
                         commandline_params => '-epojar=> /software/ensembl/genebuild/usrlocalensembllib/eponine-scan.jar, -threshold=> 0.999',
                       },
        -rc_name    => 'simple_features',
        -flow_into => {
                         1 => ['run_firstef'],
                      },
       -hive_capacity => 900,
       -batch_size => 20,
      },

            {
        # Run FirstEF. This differs slightly from the other in that it uses repeatmasking so in pipeline terms
        # it is set to run after repeatmasker instead of in parallel, unlike the others
        -logic_name => 'run_firstef',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveFirstEF',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         logic_name => 'firstef',
                         module => 'HiveFirstEF',
                         firstef_path => $self->o('firstef_path'),
                         repeat_masking_logic_names => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
                         commandline_params => '-repeatmasked',
                       },
        -rc_name    => 'simple_features',
        -flow_into => {
                         1 => ['run_cpg'],
                      },
       -hive_capacity => 900,
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
                         repeat_masking_logic_names => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
                       },
        -rc_name    => 'genscan',
        -flow_into => {
                        1 => ['create_prediction_transcript_ids'],
                        -1 => ['decrease_genscan_slice_size'],
                        -2 => ['decrease_genscan_slice_size'],
                      },
       -hive_capacity => 900,
       -batch_size => 20,
      },



            {
        # Create 1mb toplevel slices, each species flow into this independantly
        -logic_name => 'decrease_genscan_slice_size',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         split_slice => 1,
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
                         repeat_masking_logic_names => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
                       },
        -rc_name    => 'genscan',
        -flow_into => {
                        1 => ['create_prediction_transcript_ids'],
                        -1 => ['failed_genscan_slices'],
                        -2 => ['failed_genscan_slices'],
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
                         target_db => $self->o('reference_db'),
                         feature_type => 'prediction_transcript',
                         slice_to_feature_ids => 1,
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
                         commandline_params => '-cpus 3 -hitdist 40',
                         repeat_masking_logic_names => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
                         prediction_transcript_logic_names => ['genscan'],
                         iid_type => 'feature_id',
                         logic_name => 'uniprot',
                         module => 'HiveBlastGenscanPep',
                         config_settings => $self->get_config_settings('HiveBlast','HiveBlastGenscanPep'),
                      },
        -flow_into => {
                        -1 => ['run_uniprot_blast_himem'],
                        1 => ['run_vertrna_blast'],
                      },
        -rc_name    => 'blast',
        -failed_job_tolerance => 0.5,
        -hive_capacity => 900,
        -batch_size => 20,
      },

            {
        -logic_name => 'run_uniprot_blast_himem',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlastGenscanPep',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         blast_db_path => $self->o('uniprot_blast_db_path'),
                         blast_exe_path => $self->o('uniprot_blast_exe_path'),
                         commandline_params => '-cpus 3 -hitdist 40',
                         repeat_masking_logic_names => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
                         prediction_transcript_logic_names => ['genscan'],
                         iid_type => 'feature_id',
                         logic_name => 'uniprot',
                         module => 'HiveBlastGenscanPep',
                         config_settings => $self->get_config_settings('HiveBlast','HiveBlastGenscanPep'),
                      },
        -rc_name    => 'blast_himem',
        -can_be_empty  => 1,
        -flow_into => {
                        1 => ['run_vertrna_blast'],
                      },
       -failed_job_tolerance => 100,
       -hive_capacity => 900,
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
                         commandline_params => '-cpus 3 -hitdist 40',
                         repeat_masking_logic_names => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
                         prediction_transcript_logic_names => ['genscan'],
                         iid_type => 'feature_id',
                         logic_name => 'vertrna',
                         module => 'HiveBlastGenscanDNA',
                         config_settings => $self->get_config_settings('HiveBlast','HiveBlastGenscanVertRNA'),
                      },
        -flow_into => {
                        -1 => ['run_vertrna_blast_himem'],
                        1 => ['run_unigene_blast'],
                      },
        -rc_name    => 'blast',
       -failed_job_tolerance => 0.5,
       -hive_capacity => 900,
       -batch_size => 20,
      },

            {
        -logic_name => 'run_vertrna_blast_himem',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlastGenscanDNA',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         blast_db_path => $self->o('vertrna_blast_db_path'),
                         blast_exe_path => $self->o('vertrna_blast_exe_path'),
                         commandline_params => '-cpus 3 -hitdist 40',
                         repeat_masking_logic_names => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
                         prediction_transcript_logic_names => ['genscan'],
                         iid_type => 'feature_id',
                         logic_name => 'vertrna',
                         module => 'HiveBlastGenscanDNA',
                         config_settings => $self->get_config_settings('HiveBlast','HiveBlastGenscanVertRNA'),
                      },
        -rc_name    => 'blast_himem',
        -can_be_empty  => 1,
        -flow_into => {
                        1 => ['run_unigene_blast'],
                      },
       -failed_job_tolerance => 100,
       -hive_capacity => 900,
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
                         commandline_params => '-cpus 3 -hitdist 40',
                         prediction_transcript_logic_names => ['genscan'],
                         iid_type => 'feature_id',
                         repeat_masking_logic_names => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
                         logic_name => 'unigene',
                         module => 'HiveBlastGenscanDNA',
                         config_settings => $self->get_config_settings('HiveBlast','HiveBlastGenscanUnigene'),
                      },
        -flow_into => {
                        -1 => ['run_unigene_blast_himem'],
                      },
        -rc_name    => 'blast',
        -failed_job_tolerance => 0.5,
        -hive_capacity => 900,
        -batch_size => 20,
      },

      {
        -logic_name => 'run_unigene_blast_himem',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlastGenscanDNA',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         blast_db_path => $self->o('unigene_blast_db_path'),
                         blast_exe_path => $self->o('unigene_blast_exe_path'),
                         commandline_params => '-cpus 3 -hitdist 40',
                         prediction_transcript_logic_names => ['genscan'],
                         iid_type => 'feature_id',
                         repeat_masking_logic_names => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
                         logic_name => 'unigene',
                         module => 'HiveBlastGenscanDNA',
                         config_settings => $self->get_config_settings('HiveBlast','HiveBlastGenscanUnigene'),
                      },
        -rc_name    => 'blast_himem',
        -can_be_empty  => 1,
        -failed_job_tolerance => 100,
        -hive_capacity => 900,
      },

     {
        # Creates a reference db for each species
        -logic_name => 'backup_core_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         'source_db'        => $self->o('reference_db'),
                         'user_w'           => $self->o('user_w'),
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
                                ' -host '.$self->o('reference_db','-host').
                                ' -port '.$self->o('reference_db','-port').
                                ' -dbname '.$self->o('reference_db','-dbname').
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
                         subject => 'AUTOMATED REPORT: assembly loading and feature annotation for '.$self->o('reference_db','-dbname').' completed',
                         text => 'Assembly loading and feature annotation have completed for '.$self->o('reference_db','-dbname').". Basic stats can be found below",
                         file => $self->o('output_path').'/loading_report.txt',
                       },
        -rc_name => 'default',
      },

     {
        -logic_name => 'set_repeat_types',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('enscode_root_dir').'/ensembl/misc-scripts/repeats/repeat-types.pl'.
                                ' -user '.$self->o('user_w').
                                ' -pass '.$self->o('password').
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
                         cmd => 'perl '.$self->o('enscode_root_dir').'/ensembl/misc-scripts/meta_coord/update_meta_coord.pl'.
                                ' -user '.$self->o('user_w').
                                ' -pass '.$self->o('password').
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
                         cmd => 'perl '.$self->o('enscode_root_dir').'/ensembl/misc-scripts/meta_levels.pl'.
                                ' -user '.$self->o('user_w').
                                ' -pass '.$self->o('password').
                                ' -host '.$self->o('reference_db','-host').
                                ' -port '.$self->o('reference_db','-port').
                                ' -dbname '.$self->o('reference_db','-dbname')
                       },
         -rc_name => 'default',
         -flow_into => { 1 => ['create_genewise_output_db'] },
      },


######################################################################################
#
# Protein models (genblast and genewise)
#
######################################################################################

      {
        -logic_name => 'create_genewise_output_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('reference_db'),
                         target_db => $self->o('genewise_db'),
                         create_type => 'clone',
                         script_path => $self->o('clone_db_script_path'),
                         user_r => $self->o('user_r'),
                         user_w => $self->o('user_w'),
                         pass_w => $self->o('password'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['create_genblast_output_db'],
                      },

      },


      {
        -logic_name => 'create_genblast_output_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('reference_db'),
                         target_db => $self->o('genblast_db'),
                         create_type => 'clone',
                         script_path => $self->o('clone_db_script_path'),
                         user_r => $self->o('user_r'),
                         user_w => $self->o('user_w'),
                         pass_w => $self->o('password'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['download_uniprot_files'],
                      },
      },

      {

        -logic_name => 'download_uniprot_files',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadUniProtFiles',
        -parameters => {
                         multi_query_download => $self->uniprot_clade_download(),
                       },
        -rc_name          => 'default',
        -input_ids  => [ {} ],
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
                       1 => ['load_uniprot_seqs'],
                      },
      },

      {
        -logic_name => 'load_uniprot_seqs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadSeqsPipeDB',
        -parameters => {
                         uniprot_db_path => $self->o('homology_models_path').'/'.$self->o('uniprot_db_name'),
                         uniprot_index_path => $self->o('homology_models_path').'/'.$self->o('uniprot_index_name'),
                       },
        -rc_name          => 'default',
        -flow_into => {
                       1 => ['generate_genblast_jobs'],
                      },
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
                        1 => ['genblast'],
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
                         commandline_params => ' -P wublast -gff -e '.$self->o('genblast_eval').' -c '.$self->o('genblast_cov').' ',
                         query_seq_dir => $self->o('homology_models_path').'/'.$self->o('uniprot_query_dir_name'),
                         max_rank => $self->o('genblast_max_rank'),
                         genblast_pid => $self->o('genblast_pid'),
                       },
        -rc_name    => 'genblast',
        -wait_for => ['create_genblast_output_db',],
        -flow_into => {
                        -1 => ['split_genblast_jobs'],
                        -2 => ['split_genblast_jobs'],
                      },
        -failed_job_tolerance => 50,
      },

      {
        -logic_name => 'split_genblast_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         rechunk_uniprot_accession => 1,
                         uniprot_batch_size => 1,
                         uniprot_table_name => $self->o('uniprot_table_name'),
                       },
        -rc_name      => 'default',
        -can_be_empty  => 1,
        -flow_into => {
                        1 => ['genblast_retry'],
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
                         query_seq_dir => $self->o('homology_models_path').'/'.$self->o('uniprot_query_dir_name'),
                       },
        -rc_name          => 'genblast_retry',
        -can_be_empty  => 1,
        -failed_job_tolerance => 100,
        -flow_into => {
                        -2 => ['failed_genblast_proteins'],
                      },

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
                       },
        -rc_name    => 'default',
        -input_ids  => [{"iid" => $self->o('genblast_db')}],
        -wait_for => ['genblast','genblast_retry'],
        -flow_into => {
                        1 => ['create_genewise_ids'],
                      },
      },


      {
        -logic_name => 'create_genewise_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('genblast_db'),
                         feature_id => 1,
                         feature_type => 'transcript',
                         feature_logic_names => ['genblast','genblast_not_best'],
                       },
        -flow_into => {
                        1 => ['genewise'],
                      },
        -rc_name    => 'default_himem',
        -flow_into => {
                        1 => ['genewise'],
                      },
      },



      {
        -logic_name => 'genewise',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenewise',
        -parameters => {
                         iid_type => 'feature_id',
                         feature_type => 'transcript',
                         transcript_db => $self->o('genblast_db'),
                         use_genblast_best_in_genome => 0,
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('genewise_db'),
                         logic_name => 'genewise',
                         module     => 'HiveGenewise',
                         config_settings => $self->get_config_settings('genewise_protein','genewise'),
                         genewise_cov => $self->o('genewise_cov'),
                         genewise_pid => $self->o('genewise_pid'),
                         calculate_coverage_and_pid => $self->o('genewise_calculate_coverage_and_pid'),
                      },
        -flow_into => {
                        -1 => ['genewise_retry'],
                        -2 => ['failed_genewise_proteins'],
                      },
        -rc_name => 'genewise',
        -batch_size => 10,
        -failed_job_tolerance => 50,
      },

      {
        -logic_name => 'genewise_retry',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenewise',
        -parameters => {
                         iid_type => 'feature_id',
                         feature_type => 'transcript',
                         transcript_db => $self->o('genblast_db'),
                         use_genblast_best_in_genome => 0,
                         dna_db => $self->o('dna_db'),
                         target_db => $self->o('genewise_db'),
                         logic_name => 'mini_genewise',
                         module     => 'HiveGenewise',
                         config_settings => $self->get_config_settings('genewise_protein','genewise'),
                         genewise_cov => $self->o('genewise_cov'),
                         genewise_pid => $self->o('genewise_pid'),
                         calculate_coverage_and_pid => $self->o('genewise_calculate_coverage_and_pid'),
                      },
        -flow_into => {
                        -2 => ['failed_genewise_proteins'],
                      },
        -rc_name => 'genewise_retry',
        -batch_size => 10,
        -failed_job_tolerance => 50,
        -can_be_empty  => 1,
      },

      {
        -logic_name => 'classify_genewise_models',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport',
        -parameters => {
                         classification_type => 'standard',
                         update_gene_biotype => 1,
                       },
        -rc_name    => 'default',
        -input_ids  => [{"iid" => $self->o('genewise_db')}],
        -wait_for => ['genewise','genewise_retry'],
      },

      {
        -logic_name => 'failed_genewise_proteins',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name          => 'default',
        -can_be_empty  => 1,
        -failed_job_tolerance => 100,
      },

########################################################################
#
# Projection analyses
#
########################################################################

      {
        -logic_name => 'wait_for_lastz_db',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name    => 'default',
        -flow_into  => {
                         1 => ['create_projection_output_db'],
                       },
      },

      {
        -logic_name => 'create_projection_output_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('reference_db'),
                         target_db => $self->o('projection_db'),
                         create_type => 'clone',
                         script_path => $self->o('clone_db_script_path'),
                         user_r => $self->o('user_r'),
                         user_w => $self->o('user_w'),
                         pass_w => $self->o('password'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['create_projection_input_ids'],
                      },
      },

      {
        -logic_name => 'create_projection_input_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('projection_source_db'),
                         feature_id => 1,
                         feature_type => 'transcript',
                         feature_restriction => 'protein_coding',
                       },

        -flow_into => {
                        1 => ['project_transcripts'],
                      },

         -rc_name    => 'default',
      },

      {
        -logic_name => 'project_transcripts',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveWGA2GenesDirect',
        -parameters => {
                         logic_name => 'project_transcripts',
                         module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveWGA2GenesDirect',
                         config_settings => $self->get_config_settings('projection','project_transcripts'),
                         iid_type => 'feature_id',
                         feature_type => 'transcript',
                         calculate_coverage_and_pid => $self->o('projection_calculate_coverage_and_pid'),
                         max_internal_stops => $self->o('projection_max_internal_stops'),
                       },

        -rc_name          => 'project_transcripts',
        -batch_size => 100,
        -failed_job_tolerance => 0.5,
      },

      {
        -logic_name => 'classify_projection_models',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport',
        -parameters => {
                         classification_type => 'standard',
                         update_gene_biotype => 1,
                       },
        -rc_name    => 'default',
        -input_ids  => [{"iid" => $self->o('projection_db')}],
        -wait_for => ['project_transcripts'],
      },

############################################################################
#
# RNA-seq analyses
#
############################################################################

      {
        -logic_name => 'wait_for_rnaseq_db',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name    => 'default',
        # An input id will be inserted for this by running a script when both the RNA-seq and
        # lastz db are ready. This is cos these things can't be directly linked into the pipeline yet
      },

      {
        -logic_name => 'classify_rnaseq_models',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport',
        -parameters => {
                         update_gene_biotype => 1,
                         classification_type => 'standard',
                       },
        -rc_name    => 'default',
        -input_ids  => [{"iid" => $self->o('rnaseq_db')}],
        -wait_for => ['wait_for_rnaseq_db'],
      },














































      {
        -logic_name => 'load_refseq_synonyms',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'perl '.$self->o('refseq_synonyms_script_path').
                                      ' -species '.$self->o('species_name').
                                      ' -dbuser '.$self->o('user_w').
                                      ' -dbpass '.$self->o('password').
                                      ' -dbhost '.$self->o('reference_db','-host').
                                      ' -dbport '.$self->o('reference_db','-port').
                                      ' -dbname '.$self->o('reference_db','-dbname').
                                      ' -workdir '.$self->o('output_path').'/refseq_import/',
                       },
        -rc_name => 'default',

      },

      {
        -logic_name => 'download_refseq_gff',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'mkdir -p '.$self->o('output_path').'/refseq_import/;'.
                                'wget -O '.$self->o('output_path').'/refseq_import/refseq.gff.gz'." ".$self->o('refseq_ftp_path').';'.
                                'gunzip '.$self->o('output_path').'/refseq_import/refseq.gff.gz',
                       },
        -rc_name => 'default',
        -flow_into  => {
                         1 => ['create_refseq_db'],
                       },
##download_refseq_gff##        -input_ids => [{}],
        -wait_for => ['load_refseq_synonyms'],
      },


      {
        -logic_name => 'create_refseq_db',
        -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('reference_db'),
                         target_db => $self->o('refseq_db'),
                         create_type => 'clone',
                         script_path => $self->o('clone_db_script_path'),
                         user_r => $self->o('user_r'),
                         user_w => $self->o('user_w'),
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
                                      ' -user '.$self->o('user_w').
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
# ASSEMBLY LOADING ANALYSES
#
###############################################################################

      {
        # Download the files and dir structure from the NCBI ftp site. Uses the link to a species in the ftp_link_file
        -logic_name => 'download_assembly_info',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveDownloadNCBIFtpFiles',
        -parameters => {
                         'full_ftp_path'             => $self->o('full_ftp_path'),
                         'output_path'               => $self->o('output_path'),
                         'primary_assembly_dir_name' => $self->o('primary_assembly_dir_name'),
                       },
        -rc_name    => 'default',
        -flow_into  => {
                         1 => ['find_contig_accessions'],
                       },
        -input_ids  => [{}],
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
                         1 => ['populate_production_tables'],
                       },
      },

      {
        # Creates a reference db for each species
        -logic_name => 'create_core_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         'target_db'        => $self->o('reference_db'),
                         'user_w'           => $self->o('user_w'),
                         'pass_w'           => $self->o('password'),
                         'enscode_root_dir' => $self->o('enscode_root_dir'),
                         'create_type'      => 'core_only',
                       },
        -rc_name    => 'default',
        -input_ids => [{}],

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
         -wait_for => ['create_core_db'],
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
                         'chromosomes_present'       => $self->o('chromosomes_present'),
                         'genebuilder_id'            => $self->o('genebuilder_id'),
                         'target_db'                 => $self->o('reference_db'),
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
        # I've commented this out at the moment, because for some species it fails because of the species vs subspecies name
        # problem. Really it should just be added with failure tolerance until the loading script is fixed
        -logic_name => 'load_taxonomy_info',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadTaxonomyInfo',
        -parameters => {
                         'target_db'        => $self->o('reference_db'),
                         'enscode_root_dir' => $self->o('enscode_root_dir'),
                       },
        -rc_name    => 'default',
        -flow_into  => {
                         1 => ['load_mitochondrion'],
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
                         'chromosomes_present'       => $self->o('chromosomes_present'),
                      },
        -rc_name    => 'default',
        -flow_into  => {
                         1 => ['create_1mb_slice_ids','load_refseq_synonyms'],
                       },
      },


      {
        # Create 1mb toplevel slices, each species flow into this independantly
        -logic_name => 'create_1mb_slice_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db        => $self->o('reference_db'),
                         coord_system_name => 'toplevel',
                         slice => 1,
                         slice_size => 1000000,
                         include_non_reference => 0,
                         top_level => 1,
                         min_slice_length => $self->o('min_toplevel_slice_length'),
                       },
        -flow_into => {
                        # FirstEF is left out here as it runs on repeats. Repeatmasker analyses flow into it
                        4 => ['run_repeatmasker'],
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
                         commandline_params => '-nolow -species "'.$self->o('repeatmasker_library').'" -engine "'.$self->o('repeatmasker_engine').'"',
                       },
        -rc_name    => 'repeatmasker',
        -flow_into => {
                        -1 => ['run_repeatmasker_himem'],
                        4 => ['run_dust'],
                      },
        -hive_capacity => 500,
      },

      {
        -logic_name => 'run_repeatmasker_himem',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveRepeatMasker',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         logic_name => 'repeatmasker_repbase_'.$self->o('repeatmasker_library'),
                         module => 'HiveRepeatMasker',
                         repeatmasker_path => $self->o('repeatmasker_path'),
                         commandline_params => '-nolow -species "'.$self->o('repeatmasker_library').'" -engine "'.$self->o('repeatmasker_engine').'"',
                       },
        -rc_name    => 'repeatmasker_himem',
        -flow_into => {
                         4 => ['run_dust'],
                      },

        -can_be_empty  => 1,
        -hive_capacity => 500,
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
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['create_genewise_output_db'],
                      },
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
                         4 => ['run_trf'],
                      },
        -batch_size => 100,
        -hive_capacity => 500,
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
                         4 => ['run_eponine'],
                      },
        -batch_size => 100,
        -hive_capacity => 500,
      },

      {
        # Run eponine
        -logic_name => 'run_eponine',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveEponine',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         logic_name => 'eponine',
                         module => 'HiveEponine',
                         eponine_path => $self->o('eponine_path'),
                         commandline_params => '-epojar=> /software/ensembl/genebuild/usrlocalensembllib/eponine-scan.jar, -threshold=> 0.999',
                       },
        -rc_name    => 'simple_features',
        -flow_into => {
                         4 => ['run_firstef'],
                      },
        -batch_size => 100,
        -hive_capacity => 500,
      },

      {
        # Run FirstEF. This differs slightly from the other in that it uses repeatmasking so in pipeline terms
        # it is set to run after repeatmasker instead of in parallel, unlike the others
        -logic_name => 'run_firstef',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveFirstEF',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         logic_name => 'firstef',
                         module => 'HiveFirstEF',
                         firstef_path => $self->o('firstef_path'),
                         repeat_masking_logic_names => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
                         commandline_params => '-repeatmasked',
                       },
        -rc_name    => 'simple_features',
        -flow_into => {
                         4 => ['run_cpg'],
                      },
        -batch_size => 100,
        -hive_capacity => 500,
      },

      {
        # Run CPG
        -logic_name => 'run_cpg',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveCPG',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         logic_name => 'run_cpg',
                         module => 'HiveCPG',
                         cpg_path => $self->o('cpg_path'),
                       },
        -rc_name    => 'simple_features',
        -flow_into => {
                        4 => ['run_trnascan'],
                      },
        -batch_size => 100,
        -hive_capacity => 500,
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
                         4 => ['run_genscan'],
                      },
        -batch_size => 100,
        -hive_capacity => 500,
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
                         repeat_masking_logic_names => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
                       },
        -rc_name    => 'genscan',
        -flow_into => {
                        4 => ['create_prediction_transcript_ids'],
                        -1 => ['run_genscan_himem'],
                      },
        -batch_size => 10,
        -hive_capacity => 500,
      },

      {
        -logic_name => 'run_genscan_himem',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveGenscan',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         logic_name => 'genscan',
                         module => 'HiveGenscan',
                         genscan_path => $self->o('genscan_path'),
                         repeat_masking_logic_names => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
                       },
        -rc_name    => 'genscan_himem',
        -flow_into => {
                        4 => ['create_prediction_transcript_ids'],
                      },
        -can_be_empty  => 1,
        -hive_capacity => 500,
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
                         target_db => $self->o('reference_db'),
                         feature_type => 'prediction_transcript',
                         slice_to_feature_ids => 1,
                         prediction_transcript_logic_names => ['genscan'],
                       },
        -flow_into => {
                        4 => ['run_uniprot_blast'],
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
                         commandline_params => '-cpus 3 -hitdist 40',
                         repeat_masking_logic_names => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
                         prediction_transcript_logic_names => ['genscan'],
                         iid_type => 'feature_id',
                         logic_name => 'uniprot',
                         module => 'HiveBlastGenscanPep',
                         config_settings => $self->get_config_settings('HiveBlast','HiveBlastGenscanPep'),
                      },
        -flow_into => {
                        -1 => ['run_uniprot_blast_himem'],
                        4 => ['run_vertrna_blast'],
                      },
        -rc_name    => 'blast',
        -batch_size => 10,
        -hive_capacity => 500,
      },

      {
        -logic_name => 'run_uniprot_blast_himem',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlastGenscanPep',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         blast_db_path => $self->o('uniprot_blast_db_path'),
                         blast_exe_path => $self->o('uniprot_blast_exe_path'),
                         commandline_params => '-cpus 3 -hitdist 40',
                         repeat_masking_logic_names => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
                         prediction_transcript_logic_names => ['genscan'],
                         iid_type => 'feature_id',
                         logic_name => 'uniprot',
                         module => 'HiveBlastGenscanPep',
                         config_settings => $self->get_config_settings('HiveBlast','HiveBlastGenscanPep'),
                      },
        -rc_name    => 'blast_himem',
        -can_be_empty  => 1,
        -flow_into => {
                        4 => ['run_vertrna_blast'],
                      },
        -hive_capacity => 500,
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
                         commandline_params => '-cpus 3 -hitdist 40',
                         repeat_masking_logic_names => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
                         prediction_transcript_logic_names => ['genscan'],
                         iid_type => 'feature_id',
                         logic_name => 'vertrna',
                         module => 'HiveBlastGenscanDNA',
                         config_settings => $self->get_config_settings('HiveBlast','HiveBlastGenscanVertRNA'),
                      },
        -flow_into => {
                        -1 => ['run_vertrna_blast_himem'],
                        4 => ['run_unigene_blast'],
                      },
        -rc_name    => 'blast',
        -batch_size => 10,
        -hive_capacity => 500,
      },

      {
        -logic_name => 'run_vertrna_blast_himem',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlastGenscanDNA',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         blast_db_path => $self->o('vertrna_blast_db_path'),
                         blast_exe_path => $self->o('vertrna_blast_exe_path'),
                         commandline_params => '-cpus 3 -hitdist 40',
                         repeat_masking_logic_names => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
                         prediction_transcript_logic_names => ['genscan'],
                         iid_type => 'feature_id',
                         logic_name => 'vertrna',
                         module => 'HiveBlastGenscanDNA',
                         config_settings => $self->get_config_settings('HiveBlast','HiveBlastGenscanVertRNA'),
                      },
        -rc_name    => 'blast_himem',
        -can_be_empty  => 1,
        -flow_into => {
                        4 => ['run_unigene_blast'],
                      },
        -hive_capacity => 500,
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
                         commandline_params => '-cpus 3 -hitdist 40',
                         prediction_transcript_logic_names => ['genscan'],
                         iid_type => 'feature_id',
                         repeat_masking_logic_names => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
                         logic_name => 'unigene',
                         module => 'HiveBlastGenscanDNA',
                         config_settings => $self->get_config_settings('HiveBlast','HiveBlastGenscanUnigene'),
                      },
        -flow_into => {
                        -1 => ['run_unigene_blast_himem'],
                      },
        -rc_name    => 'blast',
        -batch_size => 10,
        -hive_capacity => 500,
      },

      {
        -logic_name => 'run_unigene_blast_himem',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlastGenscanDNA',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         blast_db_path => $self->o('unigene_blast_db_path'),
                         blast_exe_path => $self->o('unigene_blast_exe_path'),
                         commandline_params => '-cpus 3 -hitdist 40',
                         prediction_transcript_logic_names => ['genscan'],
                         iid_type => 'feature_id',
                         repeat_masking_logic_names => ['repeatmasker_repbase_'.$self->o('repeatmasker_library')],
                         logic_name => 'unigene',
                         module => 'HiveBlastGenscanDNA',
                         config_settings => $self->get_config_settings('HiveBlast','HiveBlastGenscanUnigene'),
                      },
        -rc_name    => 'blast_himem',
        -can_be_empty  => 1,
        -hive_capacity => 500,
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
                         source_db => $self->o('reference_db'),
                         target_db => $self->o('layering_db'),
                         create_type => 'clone',
                         script_path => $self->o('clone_db_script_path'),
                         user_r => $self->o('user_r'),
                         user_w => $self->o('user_w'),
                         pass_w => $self->o('password'),
                       },
        -rc_name    => 'default',
        -wait_for => ['classify_genewise_models','classify_projection_models','classify_rnaseq_models'],
        -input_ids => [{}],
        -flow_into => {
                        1 => ['create_genebuilder_output_db'],
                      },
      },

      {
        -logic_name => 'create_genebuilder_output_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('reference_db'),
                         target_db => $self->o('genebuilder_db'),
                         create_type => 'clone',
                         script_path => $self->o('clone_db_script_path'),
                         user_r => $self->o('user_r'),
                         user_w => $self->o('user_w'),
                         pass_w => $self->o('password'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        1 => ['create_toplevel_slices'],
                      },
      },

      {
        # Create toplevel slices, each species flow into this independantly
        -logic_name => 'create_toplevel_slices',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('reference_db'),
                         coord_system_name => 'toplevel',
                         slice => 1,
                         include_non_reference => 0,
                         top_level => 1,
                         min_slice_length => $self->o('min_toplevel_slice_length'),
                         # These options will create only slices that have a gene on the slice in one of the feature dbs
                         feature_constraint => 1,
                         feature_type => 'gene',
                         feature_dbs => [$self->o('genewise_db'),$self->o('projection_db'),$self->o('rnaseq_db')],
                       },
        -flow_into => {
                        1 => ['layer_annotation'],
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
                         config_settings => $self->get_config_settings('layer_annotation','layers'),
                       },
        -rc_name    => 'layer_annotation',
        -flow_into  => {  1 => ['genebuilder'],
                       },
      },

      {
        -logic_name => 'genebuilder',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGeneBuilder',
        -parameters => {
                         layering_output_db => $self->o('layering_db'),
                         genebuilder_output_db => $self->o('genebuilder_db'),
                         dna_db     => $self->o('dna_db'),
                         logic_name => 'ensembl',
                         module     => 'HiveGeneBuilder',
                         config_settings => $self->get_config_settings('genebuilder','genebuilder_set'),
                       },
        -rc_name    => 'genebuilder',
      },

    ];
}



sub pipeline_wide_parameters {
    my ($self) = @_;

    return {
        %{ $self->SUPER::pipeline_wide_parameters() },  # inherit other stuff from the base class
    };
}

# override the default method, to force an automatic loading of the registry in all workers
#sub beekeeper_extra_cmdline_options {
#    my $self = shift;
#    return "-reg_conf ".$self->o("registry");
#}

sub resource_classes {
  my $self = shift;

  # Note that this builds resource requests off some of the variables at the top. Works off the idea
  # that the references all get put on one server and the pipe db is on another
  my $reference_db_server = $self->default_options()->{'reference_db_server'};
  my $pipe_db_server = $self->default_options()->{'pipe_db_server'};

  my $dna_db_server = $self->default_options()->{'dna_db_server'};
  my $genblast_db_server = $self->default_options()->{'genblast_db_server'};
  my $genewise_db_server = $self->default_options()->{'genewise_db_server'};
  my $refseq_db_server = $self->default_options()->{'refseq_db_server'};
  my $projection_db_server = $self->default_options()->{'projection_db_server'};
  my $projection_lastz_db_server = $self->default_options()->{'projection_lastz_db_server'};
  my $killlist_db_server = $self->default_options()->{'killlist_db_server'};

  my $default_mem = $self->default_options()->{'default_mem'};
  my $genblast_mem = $self->default_options()->{'genblast_mem'};
  my $genblast_retry_mem = $self->default_options()->{'genblast_retry_mem'};
  my $genewise_mem = $self->default_options()->{'genewise_mem'};
  my $genewise_retry_mem = $self->default_options()->{'genewise_retry_mem'};
  my $projection_mem = $self->default_options()->{'projection_mem'};
  my $refseq_mem = $self->default_options()->{'refseq_mem'};
  my $layer_annotation_mem = $self->default_options()->{'layer_annotation_mem'};
  my $genebuilder_mem = $self->default_options()->{'genebuilder_mem'};

  my $reference_db_server_number;
  my $pipe_db_server_number;
  my $dna_db_server_number;
  my $genblast_db_server_number;
  my $genewise_db_server_number;
  my $projection_db_server_number;
  my $projection_lastz_db_server_number;
  my $refseq_db_server_number;
  my $killlist_db_server_number;


  my $num_tokens = $self->default_options()->{'num_tokens'};

  unless($reference_db_server =~ /(\d+)$/) {
    die "Failed to parse the server number out of the reference db server name. This is needed for setting tokens\n".
        "reference_db_server: ".$reference_db_server;
  }

  $reference_db_server_number = $1;

  unless($pipe_db_server =~ /(\d+)$/) {
    die "Failed to parse the server number out of the pipeline db server name. This is needed for setting tokens\n".
        "pipe_db_server: ".$pipe_db_server;
  }

  $pipe_db_server_number = $1;


  unless($dna_db_server =~ /(\d+)$/) {
    die "Failed to parse the server number out of the pipeline db server name. This is needed for setting tokens\n".
        "dna_db_server: ".$dna_db_server;
  }

  $dna_db_server_number = $1;

  unless($genblast_db_server =~ /(\d+)$/) {
    die "Failed to parse the server number out of the genblast db server name. This is needed for setting tokens\n".
        "genblast_output_db_server: ".$genblast_db_server;
  }

  $genblast_db_server_number = $1;

  unless($genewise_db_server =~ /(\d+)$/) {
    die "Failed to parse the server number out of the genewise db server name. This is needed for setting tokens\n".
        "exonerate_output_db_server: ".$genewise_db_server;
  }

  $genewise_db_server_number = $1;


  unless($projection_db_server =~ /(\d+)$/) {
    die "Failed to parse the server number out of the genewise db server name. This is needed for setting tokens\n".
        "exonerate_output_db_server: ".$projection_db_server;
  }

  $projection_db_server_number = $1;

  unless($projection_lastz_db_server =~ /(\d+)$/) {
    die "Failed to parse the server number out of the genewise db server name. This is needed for setting tokens\n".
        "exonerate_output_db_server: ".$projection_lastz_db_server;
  }

  $projection_lastz_db_server_number = $1;


  unless($refseq_db_server =~ /(\d+)$/) {
    die "Failed to parse the server number out of the genewise db server name. This is needed for setting tokens\n".
        "refseq_db_server: ".$refseq_db_server;
  }

  $refseq_db_server_number = $1;

  unless($killlist_db_server=~ /(\d+)$/) {
    die "Failed to parse the server number out of the pipeline db server name. This is needed for setting tokens\n".
        "killlist_db_server: ".$killlist_db_server;
  }

  $killlist_db_server_number = $1;


  unless($num_tokens) {
    die "num_tokens is uninitialised or zero. num_tokens needs to be present in default_options and not zero\n".
        "num_tokens: ".$pipe_db_server;
  }



  return {
    'default' => { LSF => '-q normal -M900 -R"select[mem>900] rusage[mem=900,myens_build'.
                          $reference_db_server_number.'tok='.$num_tokens.',myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },
    'default_himem' => { LSF => '-q normal -M2900 -R"select[mem>2900] rusage[mem=2900,myens_build'.
                          $reference_db_server_number.'tok='.$num_tokens.',myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },
    'repeatmasker' => { LSF => '-q normal -M2900 -R"select[mem>2900] rusage[mem=2900,myens_build'.
                                $reference_db_server_number.'tok='.$num_tokens.',myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },
    'repeatmasker_himem' => { LSF => '-q normal -M5900 -R"select[mem>5900] rusage[mem=5900,myens_build'.
                                     $reference_db_server_number.'tok='.$num_tokens.',myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },
    'simple_features' => { LSF => '-q normal -M2900 -R"select[mem>2900] rusage[mem=2900,myens_build'.
                                  $reference_db_server_number.'tok='.$num_tokens.',myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },
    'genscan' => { LSF => '-q normal -M2900 -R"select[mem>2900] rusage[mem=2900,myens_build'.
                          $reference_db_server_number.'tok='.$num_tokens.',myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },
    'genscan_himem' => { LSF => '-q normal -M5900 -R"select[mem>5900] rusage[mem=5900,myens_build'.
                          $reference_db_server_number.'tok='.$num_tokens.',myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },
    'blast' => { LSF => '-q normal -M2900 -n 3 -R "select[mem>2900] rusage[mem=2900,myens_build'.
                        $reference_db_server_number.'tok='.$num_tokens.',myens_build'.$pipe_db_server_number.'tok='.$num_tokens.'] span[hosts=1]"' },
    'blast_himem' => { LSF => '-q normal -M5900 -n 3 -R "select[mem>5900] rusage[mem=5900,myens_build'.
                        $reference_db_server_number.'tok='.$num_tokens.',myens_build'.$pipe_db_server_number.'tok='.$num_tokens.'] span[hosts=1]"' },
    'genblast' => { LSF => '-q normal -W 120 -M'.$genblast_mem.' -R"select[mem>'.$genblast_mem.'] '.
                             'rusage[mem='.$genblast_mem.','.
                             'myens_build'.$genblast_db_server_number.'tok='.$num_tokens.','.
                             'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
                             'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },


    'genblast_retry' => { LSF => '-q normal -W 120 -M'.$genblast_retry_mem.' -R"select[mem>'.$genblast_retry_mem.'] '.
                                   'rusage[mem='.$genblast_retry_mem.','.
                                   'myens_build'.$genblast_db_server_number.'tok='.$num_tokens.','.
                                   'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
                                   'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

    'genewise' => { LSF => '-q normal -W 120 -M'.$genewise_mem.' -R"select[mem>'.$genewise_mem.'] '.
                              'rusage[mem='.$genewise_mem.','.
                              'myens_build'.$genewise_db_server_number.'tok='.$num_tokens.','.
                              'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
                              'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

    'genewise_retry' => { LSF => '-q normal -W 120 -M'.$genewise_retry_mem.' -R"select[mem>'.$genewise_retry_mem.'] '.
                                    'rusage[mem='.$genewise_retry_mem.','.
                                    'myens_build'.$genewise_db_server_number.'tok='.$num_tokens.','.
                                    'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
                                    'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

    'project_transcripts' => { LSF => '-q normal -W 120 -M'.$projection_mem.' -R"select[mem>'.$projection_mem.'] '.
                              'rusage[mem='.$projection_mem.','.
                              'myens_build'.$projection_db_server_number.'tok='.$num_tokens.','.
                              'myens_build'.$projection_lastz_db_server_number.'tok='.$num_tokens.','.
                              'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
                              'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

    'refseq_import' => { LSF => '-q normal -M'.$refseq_mem.' -R"select[mem>'.$refseq_mem.'] '.
                                    'rusage[mem='.$refseq_mem.','.
                                    'myens_build'.$refseq_db_server_number.'tok='.$num_tokens.','.
                                    'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
                                    'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

    'layer_annotation' => { LSF => '-q normal -W 200 -M'.$layer_annotation_mem.' -R"select[mem>'.$layer_annotation_mem.'] '.
                             'rusage[mem='.$layer_annotation_mem.','.
                             'myens_build'.$genblast_db_server_number.'tok='.$num_tokens.','.
                             'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
                             'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

    'genebuilder' => { LSF => '-q normal -W 200 -M'.$genebuilder_mem.' -R"select[mem>'.$genebuilder_mem.'] '.
                             'rusage[mem='.$genebuilder_mem.','.
                             'myens_build'.$genblast_db_server_number.'tok='.$num_tokens.','.
                             'myens_build'.$dna_db_server_number.'tok='.$num_tokens.','.
                             'myens_build'.$pipe_db_server_number.'tok='.$num_tokens.']"' },

  }
}

sub get_config_settings {

   # This is a helper sub created to access parameters that historically were held in separate configs in the
   # old pipeline. These are now stored in the master_config_settings sub below this one. In the analyses hashes
   # earlier in the config sets of these param can be requested and stored in the config_settings hash which
   # is them passed in as a parameter to the analysis. The converted analysis modules have code to take the key
   # value pairs from the config_settings hash and assign the values to the getter/setter sub associated with the
   # key.

   # Shift in the group name (a hash that has a collection of logic name hashes and a default hash)
   # Shift in the logic name of the specific analysis
   my $self = shift;
   my $config_group = shift;
   my $config_logic_name = shift;

   # And additional hash keys will be stored in here
   my @additional_configs = @_;

   # Return a ref to the master hash for the group using the group name
   my $config_group_hash = $self->master_config_settings($config_group);
   unless(defined($config_group_hash)) {
     die "You have asked for a group name in master_config_settings that doesn't exist. Group name:\n".$config_group;
   }
   # Final hash is the hash reference that gets returned. It is important to note that the keys added have
   # priority based on the call to this subroutine, with priority from left to right. Keys assigned to
   # $config_logic_name will have most priority, then keys in any additional hashes, then keys from the
   # default hash. A default hash key will never override a $config_logic_name key
   my $final_hash;

   # Add keys from the logic name hash
   my $config_logic_name_hash = $config_group_hash->{$config_logic_name};
   unless(defined($config_logic_name_hash)) {
     die "You have asked for a logic name hash that doesn't exist in the group you specified.\n".
         "Group name:\n".$config_group."\nLogic name:\n".$config_logic_name;
   }

   $final_hash = $self->add_keys($config_logic_name_hash,$final_hash);

   # Add keys from any additional hashes passed in, keys that are already present will not be overriden
   foreach my $additional_hash (@additional_configs) {
     my $config_additional_hash = $config_group_hash->{$additional_hash};
     $final_hash = $self->add_keys($config_additional_hash,$final_hash);
   }

   # Default is always loaded and has the lowest key value priority
   my $config_default_hash = $config_group_hash->{'Default'};
   $final_hash = $self->add_keys($config_default_hash,$final_hash);

   return($final_hash);
}

sub add_keys {
  my ($self,$hash_to_add,$final_hash) = @_;

  foreach my $key (keys(%$hash_to_add)) {
    unless(exists($final_hash->{$key})) {
      $final_hash->{$key} = $hash_to_add->{$key};
    }
  }

  return($final_hash);
}

sub master_config_settings {

  my ($self,$config_group) = @_;
  my $master_config_settings = {

  HiveBlast => {
    Default => {
      BLAST_PARSER => 'Bio::EnsEMBL::Analysis::Tools::BPliteWrapper',
      PARSER_PARAMS => {
        -regex => '^(\w+)',
        -query_type => undef,
        -database_type => undef,
      },
      BLAST_FILTER => undef,
      FILTER_PARAMS => {},
      BLAST_PARAMS => {
        -unknown_error_string => 'FAILED',
        -type => 'wu',
      }
    },

    HiveBlastGenscanPep => {
      BLAST_PARSER => 'Bio::EnsEMBL::Analysis::Tools::FilterBPlite',
      PARSER_PARAMS => {
                         -regex => '^(\w+\W\d+)',
                         -query_type => 'pep',
                         -database_type => 'pep',
                         -threshold_type => 'PVALUE',
                         -threshold => 0.01,
                       },
      BLAST_FILTER => 'Bio::EnsEMBL::Analysis::Tools::FeatureFilter',
      FILTER_PARAMS => {
                         -min_score => 200,
                         -prune => 1,
                       },
    },

    HiveBlastGenscanVertRNA => {
      BLAST_PARSER => 'Bio::EnsEMBL::Analysis::Tools::FilterBPlite',
      PARSER_PARAMS => {
                         -regex => '^(\w+\W\d+)',
                         -query_type => 'pep',
                         -database_type => 'dna',
                         -threshold_type => 'PVALUE',
                         -threshold => 0.001,
                       },
      BLAST_FILTER => 'Bio::EnsEMBL::Analysis::Tools::FeatureFilter',
      FILTER_PARAMS => {
                         -prune => 1,
                       },
    },

    HiveBlastGenscanUnigene => {
      BLAST_PARSER => 'Bio::EnsEMBL::Analysis::Tools::FilterBPlite',
      PARSER_PARAMS => {
                         -regex => '\/ug\=([\w\.]+)',
                         -query_type => 'pep',
                         -database_type => 'dna',
                         -threshold_type => 'PVALUE',
                         -threshold => 0.001,
                       },
      BLAST_FILTER => 'Bio::EnsEMBL::Analysis::Tools::FeatureFilter',
      FILTER_PARAMS => {
                         -prune => 1,
                       },
      },
    },

    genewise_protein => {

      Default => {},

      genewise => {
        GENEWISE_PARAMETERS => {
                                  # pass parameters go genewise here, i.e. -program =>"/usr/local/ensembl/bin/genewiseXXX"
                                  # for more options which can be passed see Runnable/Genewise.pm
                                  -endbias => 1,
                                  -matrix => 'BLOSUM80.bla',
                                  -gap => 20,
                                  -extension => 8,
                                  -splice_model => 0
                               },
        MINIGENEWISE_PARAMETERS => {
                                     -terminal_padding => 20000,
                                     -exon_padding => 200,
                                     -minimum_intron => 1000,
                                   },

                   FILTER_PARAMS => {
                           -max_exon_length => '20000',
                           -multi_exon_min_coverage => '25',
                           -single_exon_min_coverage => '80',
                           -max_intron_length => '100000',
                           -min_split_coverage => 95,
                           -max_low_complexity => 101,
                         },

                   LIMIT_TO_FEATURE_RANGE => 1,
        FEATURE_RANGE_PADDING => 20000,
      },
    },

  projection => {
      Default => {},
      project_transcripts => {
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
      },
    },

    layer_annotation => {
      Default => {
        TARGETDB_REF => $self->o('layering_db'),
        SOURCEDB_REFS => $self->o('input_gene_dbs'),
        # ordered list of annotation layers. Genes from lower layers
        # are only retained if they do not "interfere" with genes from
        # higher layers. Genes in "Discard" layers are when assessing
        # interference, but are not written to the final database

        # Filtering is using done at the exon-overlap level
        # When no FILTER exists in this file, this is the default behaviour

        # If you would like to filter in a different way, please specify filter
        #FILTER => 'Bio::EnsEMBL::Analysis::Tools::GenomeOverlapFilter',
        #FILTER => 'Bio::EnsEMBL::Analysis::Tools::AllExonOverlapFilter',
        FILTER => 'Bio::EnsEMBL::Analysis::Tools::CodingExonOverlapFilter',
      },

      layers => {
        LAYERS => $self->layering_set(),
      },
    },

    genebuilder => {
      Default => {
      },
      genebuilder_set => {
        INPUT_GENES => {
          'input_db' => $self->genebuilder_set(),
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
    },

  };

  return($master_config_settings->{$config_group});

}

sub uniprot_clade_download {
  my ($self) = @_;

  my $clade = $self->default_options()->{'uniprot_set'};
  my $output_path = $self->o('homology_models_path');

  my $taxon_ids = {
                   'human_taxon_id'    => '9606',
                   'mammals_taxon_id'  => '40674',
                   'mouse_taxon_id'    => '10090',
                   'primates_taxon_id'  => '9443',
                   'rodents_taxon_id'  => '9989',
                   'vert_taxon_id'     => '7742',
                 };

  if($clade eq 'primates_basic') {
    return ({
              self_pe12 =>{
                            file_name => 'self_pe12.fasta',
                            taxon_id  => $self->o('taxon_id'),
                            dest_dir  => $output_path,
                            compress  => 0,
                            pe_level  => [1,2],
                          },

              human_pe12 => {
                              file_name => 'human_pe12.fasta',
                              taxon_id  => $taxon_ids->{'human_taxon_id'},
                              dest_dir  => $output_path,
                              compress  => 0,
                              pe_level  => [1,2],
                            },

              primates_pe12 => {
                             file_name => 'primate_pe12.fasta',
                             taxon_id  => $taxon_ids->{'primates_taxon_id'},
                             exclude_id => [$self->o('taxon_id'),$taxon_ids->{'human_taxon_id'}],
                             dest_dir  => $output_path,
                             compress  => 0,
                             pe_level  => [1,2],
                           },

                primates_pe3 => {
                             file_name => 'primate_pe3.fasta',
                             taxon_id  => $taxon_ids->{'primates_taxon_id'},
                             dest_dir  => $output_path,
                             compress  => 0,
                             pe_level  => [3],
                           },

                primates_pe45 => {
                             file_name => 'primate_pe45.fasta',
                             taxon_id  => $taxon_ids->{'primates_taxon_id'},
                             dest_dir  => $output_path,
                             compress  => 0,
                             pe_level  => [4,5],
                           },

               mammals_pe12 => {
                                 file_name  => 'mammal_pe12.fasta',
                                 taxon_id   => $taxon_ids->{'mammals_taxon_id'},
                                 exclude_id => [$taxon_ids->{'primates_taxon_id'}],
                                 dest_dir   => $output_path,
                                 compress   => 0,
                                 pe_level   => [1,2],
                               },

               vert_pe12 => {
                              file_name  => 'vert_pe12.fasta',
                              taxon_id   => $taxon_ids->{'vert_taxon_id'},
                              exclude_id => [$taxon_ids->{'mammals_taxon_id'}],
                              dest_dir   => $output_path,
                              compress   => 0,
                              pe_level   => [1,2],
                            },

             });
  } elsif($clade eq 'rodents_basic') {
    return ({
              self_pe12 =>{
                            file_name => 'self_pe12.fasta',
                            taxon_id  => $self->o('taxon_id'),
                            dest_dir  => $output_path,
                            compress  => 0,
                            pe_level  => [1,2],
                          },

              human_pe12 => {
                              file_name => 'human_pe12.fasta',
                              taxon_id  => $taxon_ids->{'human_taxon_id'},
                              dest_dir  => $output_path,
                              compress  => 0,
                              pe_level  => [1,2],
                            },

              mouse_pe12 => {
                              file_name  => 'mouse_pe12.fasta',
                              taxon_id   => $taxon_ids->{'mouse_taxon_id'},
                              dest_dir   => $output_path,
                              compress   => 0,
                              pe_level   => [1,2],
                            },

              rodents_pe12 => {
                                file_name => 'rodent_pe12.fasta',
                                taxon_id  => $taxon_ids->{'rodents_taxon_id'},
                                exclude_id => [$taxon_ids->{'mouse_taxon_id'},$self->o('taxon_id')],
                                dest_dir  => $output_path,
                                compress  => 0,
                                pe_level  => [1,2],
                              },

              rodents_pe3 => {
                               file_name => 'rodent_pe3.fasta',
                               taxon_id  => $taxon_ids->{'rodents_taxon_id'},
                               dest_dir  => $output_path,
                               compress  => 0,
                               pe_level  => [3],
                             },

              rodents_pe45 => {
                                file_name => 'rodent_pe45.fasta',
                                taxon_id  => $taxon_ids->{'rodents_taxon_id'},
                                dest_dir  => $output_path,
                                compress  => 0,
                                pe_level  => [4,5],
                              },


               mammals_pe12 => {
                                 file_name  => 'mammal_pe12.fasta',
                                 taxon_id   => $taxon_ids->{'mammals_taxon_id'},
                                 exclude_id => [$taxon_ids->{'rodents_taxon_id'},$taxon_ids->{'human_taxon_id'}],
                                 dest_dir   => $output_path,
                                 compress   => 0,
                                 pe_level   => [1,2],
                               },

               vert_pe12 => {
                              file_name  => 'vert_pe12.fasta',
                              taxon_id   => $taxon_ids->{'vert_taxon_id'},
                              exclude_id => [$taxon_ids->{'mammals_taxon_id'}],
                              dest_dir   => $output_path,
                              compress   => 0,
                              pe_level   => [1,2],
                            },

             });
  } else {
    die "Unknown clade selected for UniProt protein download: ".$clade;
  }

}

sub layering_set {
  my ($self) = @_;

  my $clade = $self->default_options()->{'uniprot_set'};

  if($clade eq 'primates_basic') {
    return ([
            {
              ID         => 'LAYER1',
              BIOTYPES   => [
                             'projection_95',
                             'projection_80',
                             'rnaseq_95',
                             'rnaseq_80',
                             'self_pe12_sp_95',
                             'self_pe12_tr_95',
                             'self_pe12_sp_80',
                             'self_pe12_tr_80',
                             'human_pe12_sp_95',
                             'human_pe12_tr_95',
                             'primates_pe12_sp_95',
                             'primates_pe12_tr_95',
                             'mammals_pe12_sp_95',
                             'mammals_pe12_tr_95',
                            ],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER2',
              BIOTYPES   => [
                             'human_pe12_sp_80',
                             'human_pe12_tr_80',
                             'primates_pe12_sp_80',
                             'primates_pe12_tr_80',
                             'mammals_pe12_sp_80',
                             'mammals_pe12_tr_80',
                            ],
              FILTER_AGAINST => ['LAYER1'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER3',
              BIOTYPES   => [
                             'primates_pe3_sp_95',
                             'vert_pe12_sp_95',
                             'vert_pe12_tr_95',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER4',
              BIOTYPES   => [
                             'primates_pe3_sp_80',
                             'vert_pe12_sp_80',
                             'vert_pe12_tr_80',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER5',
              BIOTYPES   => [
                              'projection_50',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4'],
              DISCARD    => 0,
            },


    ]);
  } elsif($clade eq 'rodents_basic') {
    return ([
            {
              ID         => 'LAYER1',
              BIOTYPES   => ['projection_95','projection_80',
                            'rnaseq_95','rnaseq_80',
                            ],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER2',
              BIOTYPES   => [ 'self_pe12_sp_95','self_pe12_sp_80',
                              'mouse_pe12_sp_95','mouse_pe12_sp_80',
                            ],
              FILTER_AGAINST => ['LAYER1'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER3',
              BIOTYPES   => [ 'self_pe12_tr_95','self_pe12_tr_80',
                              'mouse_pe12_tr_95','mouse_pe12_tr_80',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER4',
              BIOTYPES   => ['rodents_pe12_sp_95','rodents_pe12_sp_80'
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER5',
              BIOTYPES   => ['rodents_pe12_tr_95','rodents_pe12_tr_80',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4'],
              DISCARD    => 0,
            },


            {
              ID         => 'LAYER6',
              BIOTYPES   => ['human_pe12_sp_95','human_pe12_sp_80',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER7',
              BIOTYPES   => ['human_pe12_tr_95','human_pe12_tr_80',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER8',
              BIOTYPES   => [
                              'mammals_pe12_sp_95','mammals_pe12_sp_80'
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER7'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER9',
              BIOTYPES   => [
                              'rodents_pe3_sp_95','rodents_pe3_sp_80'
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER8'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER10',
              BIOTYPES   => [
                              'mammals_pe12_tr_95','mammals_pe12_tr_80'
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER8','LAYER9'],
              DISCARD    => 0,
            },


            {
              ID         => 'LAYER11',
              BIOTYPES   => [
                              'rodents_pe3_tr_95','rodents_pe3_tr_80'
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER8','LAYER9','LAYER10'],
              DISCARD    => 0,
            },


            {
              ID         => 'LAYER12',
              BIOTYPES   => [
                              'vert_pe12_sp_95','vert_pe12_sp_80'
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER8','LAYER9','LAYER10','LAYER11'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER13',
              BIOTYPES   => [
                              'vert_pe12_tr_95','vert_pe12_tr_80'
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER8','LAYER9','LAYER10','LAYER11','LAYER12'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER14',
              BIOTYPES   => [
                              'projection_50'
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER8','LAYER9','LAYER10','LAYER11','LAYER12','LAYER13'],
              DISCARD    => 0,
            },
            ]);
  } else {
    die "Unknown clade selected for layering: ".$clade;
  }

}

sub genebuilder_set {
  my ($self) = @_;

  my $clade = $self->default_options()->{'uniprot_set'};

  if($clade eq 'primates_baisc') {
    return ([
              'projection_80',
              'projection_95',
              'projection_50',
              'self_pe12_sp_95',
              'self_pe12_sp_80',
              'self_pe12_tr_95',
              'self_pe12_tr_80',
              'rnaseq_95',
              'rnaseq_80',
              'vert_pe12_sp_95',
              'vert_pe12_sp_80',
              'vert_pe12_tr_80',
              'vert_pe12_tr_95',
              'primates_pe12_sp_80',
              'primates_pe12_tr_80',
              'primates_pe12_tr_95',
              'primates_pe12_sp_95',
              'human_pe12_sp_80',
              'human_pe12_sp_95',
              'human_pe12_tr_80',
              'human_pe12_tr_95',
              'primates_pe3_sp_95',
              'primates_pe3_tr_80',
              'primates_pe3_tr_95',
              'primates_pe3_sp_80',
              'mammals_pe12_sp_80',
              'mammals_pe12_sp_95',
              'mammals_pe12_tr_95',
              'mammals_pe12_tr_80',
    ]);
  } elsif($clade eq 'rodents_basic') {
    return ([
              'projection_80',
              'projection_95',
              'projection_50',
              'self_pe12_sp_95',
              'self_pe12_sp_80',
              'self_pe12_tr_95',
              'self_pe12_tr_80',
              'rnaseq_95',
              'rnaseq_80',
              'mouse_pe12_sp_95',
              'mouse_pe12_sp_80',
              'mouse_pe12_tr_95',
              'mouse_pe12_tr_80',
              'vert_pe12_sp_95',
              'vert_pe12_sp_80',
              'vert_pe12_tr_80',
              'vert_pe12_tr_95',
              'rodents_pe12_sp_80',
              'rodents_pe12_tr_80',
              'rodents_pe12_tr_95',
              'rodents_pe12_sp_95',
              'human_pe12_sp_80',
              'human_pe12_sp_95',
              'human_pe12_tr_80',
              'human_pe12_tr_95',
              'rodents_pe3_sp_95',
              'rodents_pe3_tr_80',
              'rodents_pe3_tr_95',
              'rodents_pe3_sp_80',
              'rodents_pe45_sp_95',
              'rodents_pe45_tr_80',
              'rodents_pe45_tr_95',
              'rodents_pe45_sp_80',
              'mammals_pe12_sp_80',
              'mammals_pe12_sp_95',
              'mammals_pe12_tr_95',
              'mammals_pe12_tr_80',
            ]);
  } else {
    die "Unknown clade selected for genebuilding: ".$clade;
  }

}

1;
