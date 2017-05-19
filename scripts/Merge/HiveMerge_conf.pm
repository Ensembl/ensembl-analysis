# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2017] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

package HiveMerge_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::Version 2.4;
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;           # Allow this particular config to use conditional dataflow

use parent ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

use File::Spec::Functions;
use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_analysis_settings);
use Env qw(ENSCODE);

sub default_options {
  my ($self) = @_;

  return {
    # Inherit other stuff from the parent class
    %{$self->SUPER::default_options()},

    # Name of the pipeline
    'pipeline_name' => 'merge',

    # Set to 1 if you want to delete the vega and core databases.
    # Only useful if you want to start the whole pipeline again.
    # You also need to specify the driver if you want the DROP to work,
    # the easiest is to set the driver in your database hash like this:
    # -driver => $self->o('pipeline_db', '-driver'),
    'drop_databases' => 0,

    # If you are working on human and mouse you need the CCDS. If you are doing the merge
    # for any other species which does NOT have CCDS, set process_ccds to 0
    'process_ccds' => 1,

    # If you are working on human and mouse and there is a patch update to be applied, set this to 1
    # to load the new patches and run the patches annotation pipeline on patches and set the patches_ftp_dir
    # to the corresponding ftp path on the NCBI FTP server where the PATCHES are located
    'patch_update' => 1,
    'patches_ftp_dir' => 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_00000XXXX.X_GRCXXX.pX/GCA_00000XXXX.X_GRCXXX.pX_assembly_structure/PATCHES/alt_scaffolds',

    # users and passwords for read-only and write access to the genebuild MySQL servers
    'pass_r' => '',
    'user_r' => '',
    'pass_w' => '',
    'user_w' => '',
    'user' => $self->o('user_w'),
    'password' => $self->o('pass_w'),

    # output directory. The merge will write log files here. The directory must not already exist.
    'output_dir' => '/hps/nobackup/production/ensembl/...',

    # name of the directory where the vega checks reports will be written. This directory must not be neither output_dir nor a subdirectory of output_dir.
    'reports_dir' => '/hps/nobackup/production/ensembl/...',

    # email to send the vega checks reports to
    'vega_checks_reports_email' => '@ebi.ac.uk',

    # email to send the missing CCDS report to
    'ccds_report_email' => '@ebi.ac.uk',

    # name of the files containing the temporary original vega and ensembl database dumps
    'vega_tmp_filename' => 'vegadump.tmp',
    'ensembl_tmp_filename' => 'ensembldump.tmp',

    # name of the file containing the temporary previous core database dump
    'prevcore_tmp_filename' => 'prevcoredump.tmp',

    # name of the file containing the list of processed gene ids after running the merge
    'processed_genes_filename' => 'havana_merge_list_processed_genes.ids',

    # name of the file containing the list of unprocessed gene ids to be copied from the Ensembl to the merge db
    # after running the merge
    'unprocessed_genes_filename' => 'havana_merge_list_unprocessed_genes.ids',

    # name of the file containing the list of genes from the Havana db to be merged
    'vega_genes_for_merge_filename' => 'vega_genes_for_merge.ids',

    # name of the file containing the list of genes from the merge db to be copied into the core db
    # (it will be all of them but this allows us to chunk the file and speed up the process by running multiple copy_genes processes in parallel)
    'merge_genes_for_copy_filename' => 'merge_genes_for_copy.ids',

    # name of the file containing the list of core db tables which will be dumped from the previous core db
    # and loaded into the new core db
    'core_db_tables_filename' => 'core_db_tables_for_dump.list',

    # name of the file containing the list of genes from the core db to be deleted just after its creation
    'core_genes_for_deletion_filename' => 'core_genes_for_deletion.ids',

    # name of the file containing the dump of the core database after the merge
    'backup_dump_core_db_after_merge_filename' => 'backup_dump_core_db_after_merge.sql',

    # name of the files containing the list of missing CCDS before copying the missing CCDS into the Ensembl db
    'ensembl_missing_ccds_filename1' => 'ensembl_missing_ccds_1.txt',
    'vega_missing_ccds_filename1' => 'vega_missing_ccds_1.txt',
    'missing_ccds_filename1' => 'missing_ccds_1.txt',
    
    # name of the files containing the list of missing CCDS after copying the missing CCDS into the Ensembl db
    # this file contents should always be empty after running the corresponding analysis
    'ensembl_missing_ccds_filename2' => 'ensembl_missing_ccds_2.txt',
    'vega_missing_ccds_filename2' => 'vega_missing_ccds_2.txt',
    'missing_ccds_filename2' => 'missing_ccds_2.txt',

    # The number of jobs in the job array for the merge script. The workload will be evenly
    # distributed over these jobs no matter what number of jobs you put here.
    'njobs' => '45',#'75',

    # The maximum number of consecutive jobs for the merge code to run at any point in time.
    # A number between 10 and 20 seems to be optimal.
    'concurrent' => '15',

    # assembly path without patch update extension required by some scripts
    'assembly_path' => 'GRCh38',

    # full path to uniprot file to be used for the script which assigns the external DB IDs and optimises the align feature tables
    'uniprot_file' => '/hps/nobackup/production/ensembl/genebuild/blastdb/uniprot/uniprot_20XX_XX/entry_loc',

    # full path to your local copy of the ensembl-analysis git repository
    'ensembl_analysis_base' => '$ENSCODE/ensembl-analysis',

    # database parameters
    'default_port' => 4527, # default port to be used for all databases except the original vega db provided by the Vega team
    'default_host' => 'genebuildXX', # default host to be used for the pipeline, vega and core databases. Used to build the resource requirements.

    'original_ensembl_host' => '', # vega database provided by the Vega team
    'original_ensembl_port' => '',
    'original_ensembl_name' => '',

    'original_vega_host' => '', # vega database provided by the Vega team
    'original_vega_port' => '',
    'original_vega_name' => '',

    'ensembl_host' => '', # ensembl database to be used for the merge
    'ensembl_name' => '',
    'ensembl_port' => '',

    'ccds_host' => '', # ccds database containing the CCDS gene set
    'ccds_name' => '',
    'ccds_port' => '',

    'prevcore_host' => '', # previous core database (available on staging or live)
    'prevcore_name' => '',
    'prevcore_port' => '',

    'production_host' => '', # production database (available on staging)
    'production_name' => '',
    'production_port' => '',

    'pipe_dbname' => '', # pipeline db, automatically created

    'vega_name' => '', # vega database to be used for the merge

    'core_name' => '', # core database which will contain the merged gene set at the end of the process

    'vega_host' => $self->o('default_host'),
    'vega_port' => $self->o('default_port'),

    'core_host' => $self->o('default_host'),
    'core_port' => $self->o('default_port'),
    
    'pipe_db_server' => $self->o('default_host'),
    'pipe_port' => $self->o('default_port'),

    ##################
    # patch annotation
    ##################

    'repeatmasker_library' => '',     # 'human' or 'mouse'
    'species_name' => '',             # 'homo_sapiens' or 'mus_musculus'
    
    # BLAST db paths
    'uniprot_blast_db_path'     => '/hps/nobackup/production/ensembl/genebuild/blastdb/uniprot/uniprot_20XX_XX/uniprot_vertebrate',
    'uniprot_index'             => '/hps/nobackup/production/ensembl/genebuild/blastdb/uniprot/uniprot_20XX_XX/entry_loc',
    'vertrna_blast_db_path'     => '/hps/nobackup/production/ensembl/genebuild/blastdb/vertrna/XXX/embl_vertrna-1',
    'unigene_blast_db_path'     => '/hps/nobackup/production/ensembl/genebuild/blastdb/unigene/unigene',
    
    # alignment annotation on patches
    'taxon_id'             => XXXX, # 9606 for human
    'clone_db_script_path' => "$ENSCODE/ensembl-analysis/scripts/clone_database.ksh",
    'genblast_name'        => 'XXX_genblast_XXX', # genblast database
    'genblast_host'        => $self->o('default_host'),
    'exonerate_settings'   => 'exonerate_protein_XXX_patch', # exonerate settings to use from ExonerateStatic
    'exonerate_name'        => 'XXX_exonerate_XXX', # genblast database
    'exonerate_host'        => $self->o('default_host'),
    'cdna_name'             => '', # latest cdna db on live-mirror
    'layering_name'        => 'XXX_layering_XXX',
    'layering_host'        => $self->o('default_host'),
    'utr_name'        => 'XXX_utr_XXX',
    'utr_host'        => $self->o('default_host'),
    'genebuilder_name'        => 'XXX_genebuilder_XXX',
    'genebuilder_host'        => $self->o('default_host'),
    'pseudogene_name'        => 'XXX_pseudogene_XXX',
    'pseudogene_host'        => $self->o('default_host'),
    'patch_geneset_name'        => 'XXX_patch_XXX',
    'patch_geneset_host'        => $self->o('default_host'),
    'homology_models_path' => $self->o('output_dir').'/homology_models/',

    ##
    # You shouldn't need to change anything below this line
    ##

    'repeatmasker_engine' => 'crossmatch',    
    'repeat_logic_names' => ['repeatmask_repbase_'.$self->o('repeatmasker_library'),
                             'dust'],
    'repeat_masking_logic_names' => ['repeatmask_repbase_'.$self->o('repeatmasker_library')],
    'uniprot_set'                => 'self_patch',

    # executable paths
    'repeatmasker_path' => '/nfs/software/ensembl/RHEL7/linuxbrew/Cellar/repeatmasker/4.0.5/bin/RepeatMasker',
    'dust_path' => '/nfs/software/ensembl/RHEL7/linuxbrew/bin/tcdust',
    'trf_path' => '/nfs/software/ensembl/RHEL7/linuxbrew/bin/trf',
    'eponine_java_path' => '/nfs/software/ensembl/RHEL7/jenv/shims/java',
    'eponine_jar_path' => '/nfs/software/ensembl/RHEL7/linuxbrew/Cellar/eponine/1.0/libexec/eponine-scan.jar',
    'cpg_path' => '/nfs/software/ensembl/RHEL7/linuxbrew/bin/cpg_lh',
    'trnascan_path' => '/nfs/software/ensembl/RHEL7/linuxbrew/bin/tRNAscan-SE',
    'genscan_path' => '/nfs/software/ensembl/RHEL7/linuxbrew/bin/genscan',
    'genscan_matrix_path' => '/nfs/software/ensembl/RHEL7/linuxbrew/share/HumanIso.smat',
    'uniprot_blast_exe_path' => '/nfs/software/ensembl/RHEL7/linuxbrew/bin/blastp',
    'vertrna_blast_exe_path' => '/nfs/software/ensembl/RHEL7/linuxbrew/bin/tblastn',
    'unigene_blast_exe_path' => '/nfs/software/ensembl/RHEL7/linuxbrew/bin/tblastn',
    'exonerate_path'         => '/nfs/software/ensembl/RHEL7/linuxbrew/opt/exonerate09/bin/exonerate',
    'cmsearch_exe_path'    => '/nfs/software/ensembl/RHEL7/linuxbrew/Cellar/infernal10/1.0/bin/cmsearch',
    'genblast_path'        => '/nfs/software/ensembl/RHEL7/linuxbrew/bin/genblast',

    'layering_input_gene_dbs' => [$self->o('exonerate_db')],
    'utr_gene_dbs' => {
                        'cdna_db'       => $self->o('cdna_db'),
                        'no_utr_db'     => $self->o('layering_db'),
                      },

    'utr_biotype_priorities'  => {'cdna_update' => 1},

    'cleaning_blessed_biotypes' => {'pseudogene' => 1,'processed_pseudogene' => 1},

    'blast_type'                 => 'ncbi',
    'uniprot_index_name'         => 'uniprot_index',
    'uniprot_db_name'            => 'uniprot_db',
    'uniprot_query_dir_name'     => 'uniprot_temp',
    'uniprot_genblast_batch_size' => 50,
    'uniprot_table_name'         => 'uniprot_sequences',
    
    'genblast_eval'              => '1e-10',
    'genblast_cov'               => '0.5',
    'genblast_pid'               => '50',
    'genblast_max_rank'          => '5',

    'exonerate_pid'              => '80', # Cut-off for percent id
    'exonerate_cov'              => '80', # Cut-off for coverage
    'exonerate_calculate_coverage_and_pid' => '1',
    'exonerate_region_padding' => 15000,#1000

    'killlist_name'        => 'gb_kill_list', # genblast database
    'killlist_host'        => 'mysql-ens-genebuild-prod-6.ebi.ac.uk',
    'killlist_port'        => '4532',

    'genome_file' => $self->o('output_dir')."/genome_dumps/".$self->o('species_name')."_softmasked_toplevel.fa",

    ## Required by HiveBaseConfig_conf
    #core_db (not required in the merge process, leave blank)
    'core_dbname' => '',
    'core_db_server' => '',
    'port' => '',
    #'user_r' => '',
    #dna_db (not required in the merge process, leave blank)
    'dna_dbname' => '',
    'dna_db_server' => '',
    #'port' => '',
    #'user_r' => '',

    # original Ensembl database (used in previous release)
    'original_ensembl_db' => {
                    -host      => $self->o('original_ensembl_host'),
                    -port      => $self->o('original_ensembl_port'),
                    -user      => $self->o('user_r'),
                    -pass      => $self->o('pass_r'),
                    -dbname    => $self->o('original_ensembl_name'),
    },

    # vega database provided by the Vega team
    'original_vega_db' => {
                    -host      => $self->o('original_vega_host'),
                    -port      => $self->o('original_vega_port'),
                    -user      => $self->o('user_r'),
                    -pass      => $self->o('pass_r'),
                    -dbname    => $self->o('original_vega_name'),
    },

    # ensembl database to be used for the merge
    'ensembl_db' => {
                    -host      => $self->o('ensembl_host'),
                    -port      => $self->o('default_port'),
                    -user      => $self->o('user_w'),
                    -pass      => $self->o('pass_w'),
                    -dbname    => $self->o('ensembl_name'),
    },

    # ccds database containing the CCDS gene set
    'ccds_db' => {
                    -host      => $self->o('ccds_host'),
                    -port      => $self->o('default_port'),
                    -user      => $self->o('user_r'),
                    -pass      => $self->o('pass_r'),
                    -dbname    => $self->o('ccds_name'),
    },

    # previous core database (available on ens-staging or ens-livemirror)
    'prevcore_db' => {
                    -host      => $self->o('prevcore_host'),
                    -port      => $self->o('default_port'),
                    -user      => $self->o('user_r'),
                    -pass      => $self->o('pass_r'),
                    -dbname    => $self->o('prevcore_name'),
    },
    'db_conn' => 'mysql://'.$self->o('user_r').'@'.$self->o('prevcore_host').'/'.$self->o('prevcore_name'),

    # pipeline db, pipeline will create this automatically
    'pipeline_db' => {
                    -host      => $self->o('default_host'),
                    -port      => $self->o('default_port'),
                    -user      => $self->o('user_w'),
                    -pass      => $self->o('pass_w'),
                    -dbname    => $self->o('pipe_dbname'),
                    -driver    => "mysql",
    },

    # vega database to be used for the merge
    'vega_db' => {
                    -host      => $self->o('vega_host'),
                    -port      => $self->o('default_port'),
                    -user      => $self->o('user_w'),
                    -pass      => $self->o('pass_w'),
                    -dbname    => $self->o('vega_name'),
                    -driver    => $self->o('pipeline_db', '-driver'),
    },

    # core database
    'core_db' => {
                    -host      => $self->o('core_host'),
                    -port      => $self->o('default_port'),
                    -user      => $self->o('user_w'),
                    -pass      => $self->o('pass_w'),
                    -dbname    => $self->o('core_name'),
                    -driver    => $self->o('pipeline_db', '-driver'),
    },
    
    # core database read-only access
    'ro_core_db' => {
                    -host      => $self->o('core_host'),
                    -port      => $self->o('default_port'),
                    -user      => $self->o('user_r'),
                    -pass      => $self->o('pass_r'),
                    -dbname    => $self->o('core_name'),
                    -driver    => $self->o('pipeline_db', '-driver'),
    },

    # genblast database
    'genblast_db' => {
                    -host      => $self->o('genblast_host'),
                    -port      => $self->o('default_port'),
                    -user      => $self->o('user_w'),
                    -pass      => $self->o('pass_w'),
                    -dbname    => $self->o('genblast_name'),
                    -driver    => $self->o('pipeline_db', '-driver'),
    },
    
    # exonerate database
    'exonerate_db' => {
                    -host      => $self->o('exonerate_host'),
                    -port      => $self->o('default_port'),
                    -user      => $self->o('user_w'),
                    -pass      => $self->o('pass_w'),
                    -dbname    => $self->o('exonerate_name'),
                    -driver    => $self->o('pipeline_db', '-driver'),
    },
    
    # patch layering database
    'layering_db' => {
                    -host      => $self->o('layering_host'),
                    -port      => $self->o('default_port'),
                    -user      => $self->o('user_w'),
                    -pass      => $self->o('pass_w'),
                    -dbname    => $self->o('layering_name'),
                    -driver    => $self->o('pipeline_db', '-driver'),
    },
    
    # patch utr database
    'utr_db' => {
                    -host      => $self->o('utr_host'),
                    -port      => $self->o('default_port'),
                    -user      => $self->o('user_w'),
                    -pass      => $self->o('pass_w'),
                    -dbname    => $self->o('utr_name'),
                    -driver    => $self->o('pipeline_db', '-driver'),
    },
    
    # patch pseudogene database
    'pseudogene_db' => {
                    -host      => $self->o('pseudogene_host'),
                    -port      => $self->o('default_port'),
                    -user      => $self->o('user_w'),
                    -pass      => $self->o('pass_w'),
                    -dbname    => $self->o('pseudogene_name'),
                    -driver    => $self->o('pipeline_db', '-driver'),
    },
    
    # patch genebuilder database
    'genebuilder_db' => {
                    -host      => $self->o('genebuilder_host'),
                    -port      => $self->o('default_port'),
                    -user      => $self->o('user_w'),
                    -pass      => $self->o('pass_w'),
                    -dbname    => $self->o('genebuilder_name'),
                    -driver    => $self->o('pipeline_db', '-driver'),
    },
    
    # patch geneset database
    'patch_geneset_db' => {
                    -host      => $self->o('patch_geneset_host'),
                    -port      => $self->o('default_port'),
                    -user      => $self->o('user_w'),
                    -pass      => $self->o('pass_w'),
                    -dbname    => $self->o('patch_geneset_name'),
                    -driver    => $self->o('pipeline_db', '-driver'),
    },
    
    # kill list database
    'killlist_db' => {
                    -host      => $self->o('killlist_host'),
                    -port      => $self->o('default_port'),
                    -user      => $self->o('user_r'),
                    -pass      => $self->o('pass_r'),
                    -dbname    => $self->o('killlist_name'),
                    -driver    => $self->o('pipeline_db', '-driver'),
    },

    # cdna database
    'cdna_db' => {
                   -host   => 'ens-livemirror',
                   -port   => $self->o('default_port'),
                   -user   => $self->o('user_r'),
                   -pass   => $self->o('pass_r'),
                   -dbname => $self->o('cdna_name'),
                   -driver => $self->o('pipeline_db', '-driver'),
    },

    # production database
    'production_db' => {
                    -host      => $self->o('production_host'),
                    -port      => $self->o('default_port'),
                    -user      => $self->o('user_r'),
                    -pass      => $self->o('pass_r'),
                    -dbname    => $self->o('production_name'),
    },

  };
}

sub resource_classes {
      my $self = shift;
          return {
            'default' => { 'LSF' => $self->lsf_resource_builder('normal',900,[$self->default_options->{'default_host'}]) },
            'normal_1500' => { 'LSF' => $self->lsf_resource_builder('normal',1500,[$self->default_options->{'default_host'}]) },
            'normal_1900' => { 'LSF' => $self->lsf_resource_builder('normal',1900,[$self->default_options->{'default_host'}]) },
            'normal_1900_120' => { 'LSF' => "-W 120 ".$self->lsf_resource_builder('normal',1900,[$self->default_options->{'default_host'}]) },
            'normal_2900' => { 'LSF' => $self->lsf_resource_builder('normal',2900,[$self->default_options->{'default_host'}]) },
            'normal_2900_180' => { 'LSF' => "-W 180 ".$self->lsf_resource_builder('normal',2900,[$self->default_options->{'default_host'}]) },
            'normal_2900_n3' => { 'LSF' => $self->lsf_resource_builder('normal',2900,[$self->default_options->{'default_host'}],undef,3) },
            'normal_3900_200' => { 'LSF' => "-W 200 ".$self->lsf_resource_builder('normal',3900,[$self->default_options->{'default_host'}]) },
            'normal_5900_60' => { 'LSF' => "-W 60 ".$self->lsf_resource_builder('normal',5900,[$self->default_options->{'default_host'}]) },
            'normal_4600' => { 'LSF' => $self->lsf_resource_builder('normal',4600,[$self->default_options->{'default_host'}]) },
            'normal_4900_120' => { 'LSF' => "-W 120 ".$self->lsf_resource_builder('normal',4900,[$self->default_options->{'default_host'}]) },
            'normal_5900' => { 'LSF' => $self->lsf_resource_builder('normal',5900,[$self->default_options->{'default_host'}]) },
            'normal_7900' => { 'LSF' => $self->lsf_resource_builder('normal',7900,[$self->default_options->{'default_host'}]) },
            'normal_12000' => { 'LSF' => $self->lsf_resource_builder('normal',12000,[$self->default_options->{'default_host'}]) },
            'basement_1500' => { 'LSF' => $self->lsf_resource_builder('basement',1500,[$self->default_options->{'default_host'}]) },
            'local' => {'LOCAL' => ''},   
                }
}

sub pipeline_create_commands {
  my ($self) = @_;

  my $create_commands = $self->SUPER::pipeline_create_commands;
  if ($self->o('drop_databases')) {
    foreach my $dbname ('vega_db','core_db') {
      push(@$create_commands,$self->db_cmd('DROP DATABASE IF EXISTS', $self->dbconn_2_url($dbname)));
    }
  }
  
  push (@$create_commands,$self->db_cmd('CREATE TABLE '.$self->o('uniprot_table_name').' ('.
                                        'accession varchar(50) NOT NULL,'.
                                        'source_db varchar(50) NOT NULL,'.
                                        'pe_level varchar(50) NOT NULL,'.
                                        'biotype varchar(255) NOT NULL,'.
                                        'group_name varchar(255) NOT NULL,'.
                                        'seq text NOT NULL,'.
                                        'PRIMARY KEY (accession))'));
  
  return $create_commands;
}

sub pipeline_analyses {
  my ($self) = @_;

  return [
            {
              -logic_name => 'create_reports_dir',
              -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                               'cmd'   => 'mkdir -p '.$self->o('reports_dir').";".
                                          'mkdir -p '.$self->o('output_dir')
                             },
              -flow_into => { 1 => ['create_vega_db'] },
              -rc_name => 'default',
              -input_ids => [ {} ],
            },

            {
              -logic_name => 'create_ensembl_db',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                                'cmd'   => 'mysql -NB -u'.$self->o('user_w').
                                               ' -p'.$self->o('pass_w').
                                               ' -h'.$self->o('ensembl_db','-host').
                                               ' -P'.$self->o('ensembl_db','-port').
                                               ' -e "CREATE DATABASE '.$self->o('ensembl_db','-dbname').'"'
                             },
              -input_ids => [ {} ],
              -rc_name => 'default',
              -flow_into => { 1 => ['list_ensembl_db_tables'] },
            },

            {
              -logic_name => 'list_ensembl_db_tables',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
              -parameters => {
                               'db_conn'    => $self->o('original_ensembl_db'),
                               'inputquery' => 'SHOW TABLE STATUS',
                             },
              -rc_name => 'default',
              -flow_into => { '2->A' => { 'parallel_dump_ensembl_db' => { 'table_name' => '#Name#' }, },
                              'A->1' => ['ensembl_db_creation_completed'] },
            },
            
            {
              -logic_name => 'ensembl_db_creation_completed',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                               'cmd'   => 'echo "Ensembl db creation completed."'
                             },
              -rc_name => 'default',
            },

            {
              -logic_name => 'parallel_dump_ensembl_db',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
              -parameters => {
                               'db_conn'       => $self->o('original_ensembl_db'),
                               'output_file'   => $self->o('output_dir').'/'.$self->o('original_ensembl_db','-dbname').'_'.'#table_name#.sql',
                               'executable'    => 'mysqldump',
                               'append'        => ['#table_name#'],
                             },
               -analysis_capacity => 200,
               -hive_capacity => 200,
               -max_retry_count => 2,
               -rc_name => 'normal_1500',
               -flow_into => { 1 => ['parallel_load_ensembl_db'] },
            },

            {
              -logic_name => 'parallel_load_ensembl_db',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
              -parameters => {
                                 db_conn => $self->o('ensembl_db'),
                                 input_file => $self->o('output_dir').'/'.$self->o('original_ensembl_db','-dbname').'_'.'#table_name#.sql',
                             },
              -rc_name => 'normal_1500',
            },         

            {
              -logic_name => 'create_vega_db',
              -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
              -parameters => {
              	                create_type => 'copy',
                                source_db => $self->o('original_vega_db'),
                                target_db => $self->o('vega_db'),
                                pass_w => $self->o('pass_w'),
                                user_w => $self->o('user_w'),
                                db_dump_file => $self->o('output_dir').'/'.$self->o('vega_tmp_filename'),
                                #ignore_dna => 1, # the vega db does not have any dna
                             },
             -flow_into => { 1 => ['list_toplevel_for_vega_checks_before'] },
            },

            {
              -logic_name => 'create_core_db',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                                'cmd'   => 'mysql -NB -u'.$self->o('user_w').
                                               ' -p'.$self->o('pass_w').
                                               ' -h'.$self->o('core_db','-host').
                                               ' -P'.$self->o('core_db','-port').
                                               ' -e "CREATE DATABASE '.$self->o('core_db','-dbname').'"'
                             },
              -input_ids => [ {} ],
              -rc_name => 'default',
              -flow_into => { 1 => ['list_core_db_tables'] },
            },

            {
              -logic_name => 'list_core_db_tables',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
              -parameters => {
                               'db_conn'    => $self->o('prevcore_db'),
                               'inputquery' => 'SHOW TABLE STATUS',
                             },
              -rc_name => 'default',
              -flow_into => { '2->A' => { 'parallel_dump_core_db' => { 'table_name' => '#Name#' }, },
                              'A->1' => ['list_core_genes'] },
            },

            {
              -logic_name => 'parallel_dump_core_db',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
              -parameters => {
                               'db_conn'       => $self->o('prevcore_db'),
                               'output_file'   => $self->o('output_dir').'/#table_name#.sql',
                               'executable'    => 'mysqldump',
                               'append'        => ['#table_name#'],
                             },
               -analysis_capacity => 200,
               -hive_capacity => 200,
               -max_retry_count => 2,
               -rc_name => 'normal_1500',
               -flow_into => { 1 => ['parallel_load_core_db'] },
            },

            {
              -logic_name => 'parallel_load_core_db',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
              -parameters => {
                                 db_conn => $self->o('core_db'),
                                 input_file => $self->o('output_dir').'/#table_name#.sql',
                             },
              -rc_name => 'normal_1500',
              #-flow_into => { 1 => ['list_core_genes'] },
            },

            {
              -logic_name => 'list_core_genes',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::DbCmd',
              -parameters => {
                               db_conn     => $self->o('core_db'),
                               append      => ['-NB'],
                               input_query => 'SELECT gene_id FROM gene g,seq_region sr WHERE g.seq_region_id=sr.seq_region_id AND name <> "MT"',
                               output_file => $self->o('output_dir').'/'.$self->o('core_genes_for_deletion_filename'),
                             },
              -rc_name => 'default',
              -flow_into => { 1 => ['chunk_core_genes'] },
            },

            {

              -logic_name => 'chunk_core_genes',
              -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::FileFactory',
              -parameters => {
                               inputfile => $self->o('output_dir').'/'.$self->o('core_genes_for_deletion_filename'),
                               output_dir => $self->o('output_dir'),
                               output_prefix => $self->o('core_genes_for_deletion_filename')."_chunk_",
                             },

              -flow_into => { '2->A' => [ 'delete_core_genes' ],
                              'A->1' => [ 'core_sql_truncates' ],
                            },
              -rc_name => 'default',
            },

            {
              -logic_name => 'delete_core_genes',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                               cmd => 'perl '.$self->o('ensembl_analysis_base').'/scripts/genebuild/delete_genes.pl'
                                     ." -dbhost ".$self->o('core_db','-host')
                                     ." -dbuser ".$self->o('core_db','-user')
                                     ." -dbpass ".$self->o('core_db','-pass')
                                     ." -dbname ".$self->o('core_db','-dbname')
                                     ." -dbport ".$self->o('core_db','-port')
                                     ." -idfile #file#"

                             },
               -analysis_capacity => 25,
               -hive_capacity => 25,
               -max_retry_count => 2,
               -rc_name => 'normal_1500',
            },
            {
              -logic_name => 'core_sql_truncates',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
              -parameters => {
                               db_conn => $self->o('core_db'),
                               sql => [ 'TRUNCATE xref',
                                        'TRUNCATE object_xref',
                                        'TRUNCATE external_synonym',
                                        'TRUNCATE dependent_xref',
                                        'TRUNCATE interpro',
                                        'TRUNCATE identity_xref',
                                        'TRUNCATE ontology_xref',
                                        'DELETE FROM unmapped_object WHERE type LIKE "xref"'
                                       ],
                             },
              -max_retry_count => 0,
              -rc_name => 'default',
            },

            {
              -logic_name => 'list_toplevel_for_vega_checks_before',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
              -parameters => {
                               db_conn => $self->o('vega_db'),
                               inputquery => 'SELECT sr.name FROM seq_region sr, seq_region_attrib sra, attrib_type at WHERE sr.seq_region_id = sra.seq_region_id AND sr.name NOT LIKE "LRG\_%" AND sra.attrib_type_id = at.attrib_type_id AND code = "toplevel"',
                               column_names => ['chromosome'],
                             },
              -flow_into => { '2->A' => [ 'vega_checks_before' ],
                              'A->1' => [ 'vega_checks_before_concat' ],
                            },
              -rc_name => 'default',
            },

            {
              -logic_name => 'vega_checks_before',
              -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveVegaChecks',
              -parameters => {
              	               dbname => $self->o('vega_db','-dbname'),
                               dbhost => $self->o('vega_db','-host'),
                               dnadbname => $self->o('ensembl_db','-dbname'),
                               dnadbhost => $self->o('ensembl_db','-host'),
                               coord_system => 'toplevel',
                               path => $self->o('assembly_path'),
                               sql_output => $self->o('output_dir').'/vega_checks_before_#chromosome#.sql',
                               dbtype => '', # can be 'vega' or '' (empty string)
                               port => '3306',
                               user => $self->o('user_w'),
                               pass => $self->o('pass_w'),
                               #chromosome => '',
                               write => 1,
                               affix => 0, # perform the checks by using the biotypes with or without the prefixes and suffixes like weird_, _Ens, _hav, ... ; without affixes by default
                               biotypes_extension => 0,
                               stdout_file => $self->o('reports_dir').'/vega_checks_before_#chromosome#.out',
                               stderr_file => $self->o('reports_dir').'/vega_checks_before_#chromosome#.err',
                             },
              -hive_capacity    => 30,
              -analysis_capacity => 30,
              -wait_for => ['ensembl_db_creation_completed'],
            },
            {
              -logic_name => 'vega_checks_before_concat',
              -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                               'cmd'   => 'cat '.$self->o('reports_dir').'/vega_checks_before_*.out > '.
                                                 $self->o('reports_dir').'/vega_checks_before.out'
                             },
              -flow_into => { 1 => ['vega_checks_before_report'] },
              -rc_name => 'default',
            },
            {
              -logic_name => 'vega_checks_before_report',
              -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::TextfileByEmail',
              -parameters => {
                               email => $self->o('vega_checks_reports_email'),
                               subject => 'AUTOMATED REPORT: vega biotype combinations',
                               text => 'Please find below the list of not allowed gene and transcript biotype combinations BEFORE the merge found in the vega database '.$self->o('vega_db','-dbname').'. Please note that any change applied to the list of allowed biotype combinations will take effect in a future release (not the next release).',
                               file => $self->o('reports_dir').'/vega_checks_before.out',
                               command => q{egrep '(Unknown)|(not allowed\.)' | awk '{print $9,$18}' | sort | uniq -c | sort -nr | sed s'/\. run/ (UNKNOWN gene biotype)/g'},
                             },
              -flow_into => { 1 => ['prepare_vega_db'] },
              -rc_name => 'default',
            },

            {
              -logic_name => 'prepare_vega_db',
              -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveVegaPreparation',
              -parameters => {
                               output_path => $self->o('output_dir').'/vega_preparation/',
                               dbhost => $self->o('vega_db','-host'),
                               dbname => $self->o('vega_db','-dbname'),
                               dbuser => $self->o('user_w'),
                               dbpass => $self->o('pass_w'),
                               dbport => $self->o('vega_db','-port'),
                               dnadbhost => $self->o('ensembl_db','-host'),
                               dnadbport => $self->o('ensembl_db','-port'),
                               dnadbname => $self->o('ensembl_db','-dbname'),
                               check_vega_met_stop_dir => $self->o('ensembl_analysis_base').'/scripts/Merge',
                               skip => 0,
                               only => 0,
                             },
              -flow_into => { 1 => WHEN ('#process_ccds#' => ['set_ccds_biotype'],
                                   ELSE                      ['list_vega_genes_for_merge']) },
            },

            {
              -logic_name => 'set_ccds_biotype',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
              -parameters => {
                               db_conn => $self->o('ensembl_db'),
                               sql => ['UPDATE transcript SET biotype="ccds" WHERE stable_id LIKE "CCDS%"'],
                             },
              -flow_into => { 1 => ['delete_ccds'] },
            },

            {
              -logic_name => 'delete_ccds',
              -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDeleteTranscripts',
              -parameters => {
                               biotype => 'ccds',
                               dbhost => $self->o('ensembl_db','-host'),
                               dbname => $self->o('ensembl_db','-dbname'),
                               dbuser => $self->o('user_w'),
                               dbpass => $self->o('pass_w'),
                               dbport => $self->o('ensembl_db','-port'),
                               delete_transcripts_path => $self->o('ensembl_analysis_base').'/scripts/genebuild/',
                               delete_genes_path => $self->o('ensembl_analysis_base').'/scripts/genebuild/',
                               delete_transcripts_script_name => 'delete_transcripts.pl',
                               delete_genes_script_name => 'delete_genes.pl',
                               output_path => $self->o('output_dir').'/delete_ccds/',
                               output_file_name => 'delete_ccds.out',
                               email => $self->o('vega_checks_reports_email'),
                               from => 'ensembl-genebuild@ebi.ac.uk'
                             },
              -max_retry_count => 0,
              -flow_into => { 1 => ['list_ensembl_toplevel_for_update_ccds_labels'] },
            },

            {
              -logic_name => 'list_ensembl_toplevel_for_update_ccds_labels',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
              -parameters => {
                               db_conn => $self->o('ensembl_db'),
                               inputquery => 'SELECT sr.name FROM seq_region sr, seq_region_attrib sra, attrib_type at WHERE sr.seq_region_id = sra.seq_region_id AND sr.name NOT LIKE "LRG\_%" AND sra.attrib_type_id = at.attrib_type_id AND code = "toplevel"',
                               column_names => ['chromosome'],
                             },
              -flow_into => { '2->A' => [ 'update_ensembl_ccds_labels' ],
                              'A->1' => [ 'update_ensembl_ccds_labels_concat' ],
                            },
              -rc_name => 'default',
            },

            {
              -logic_name => 'update_ensembl_ccds_labels',
              -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveUpdateCCDSLabels',
              -parameters => {
                               ccds_dbname => $self->o('ccds_db','-dbname'),
                               ccds_host => $self->o('ccds_db','-host'),
                               ccds_user => $self->o('ccds_db','-user'),
                               output_dbname => $self->o('ensembl_db','-dbname'),
                               output_host => $self->o('ensembl_db','-host'),
                               output_user => $self->o('user_w'),
                               output_pass => $self->o('pass_w'),
                               dna_dbname => $self->o('ensembl_db','-dbname'),
                               dna_host => $self->o('ensembl_db','-host'),
                               dna_user => $self->o('user_r'),
                               dna_pass => $self->o('pass_r'),
                               assembly_path => $self->o('assembly_path'),
                               reports_dir => $self->o('reports_dir'),
                               output_filename => $self->o('ensembl_missing_ccds_filename1'), # this will get an extra extension corresponding to the chromosome name
                               #chromosome => '',
                             },
              -hive_capacity    => 30,
              -analysis_capacity => 30,
            },
            {
              -logic_name => 'update_ensembl_ccds_labels_concat',
              -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                               'cmd'   => 'cat '.$self->o('reports_dir').'/'.$self->o('ensembl_missing_ccds_filename1').'.* > '.
                                                 $self->o('reports_dir').'/'.$self->o('ensembl_missing_ccds_filename1')
                             },
              -flow_into => { 1 => ['list_vega_toplevel_for_update_ccds_labels'] },
              -rc_name => 'default',
            },

            {
              -logic_name => 'list_vega_toplevel_for_update_ccds_labels',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
              -parameters => {
                               db_conn => $self->o('vega_db'),
                               inputquery => 'SELECT sr.name FROM seq_region sr, seq_region_attrib sra, attrib_type at WHERE sr.seq_region_id = sra.seq_region_id AND sr.name NOT LIKE "LRG\_%" AND sra.attrib_type_id = at.attrib_type_id AND code = "toplevel"',
                               column_names => ['chromosome'],
                             },
              -flow_into => { '2->A' => [ 'update_vega_ccds_labels' ],
                              'A->1' => [ 'update_vega_ccds_labels_concat' ],
                            },
              -rc_name => 'default',
            },

            {
              -logic_name => 'update_vega_ccds_labels',
              -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveUpdateCCDSLabels',
              -parameters => {
                               ccds_dbname => $self->o('ccds_db','-dbname'),
                               ccds_host => $self->o('ccds_db','-host'),
                               ccds_user => $self->o('ccds_db','-user'),
                               output_dbname => $self->o('vega_db','-dbname'),
                               output_host => $self->o('vega_db','-host'),
                               output_user => $self->o('user_w'),
                               output_pass => $self->o('pass_w'),
                               dna_dbname => $self->o('ensembl_db','-dbname'),
                               dna_host => $self->o('ensembl_db','-host'),
                               dna_user => $self->o('user_r'),
                               dna_pass => $self->o('pass_r'),
                               assembly_path => $self->o('assembly_path'),
                               reports_dir => $self->o('reports_dir'),
                               output_filename => $self->o('vega_missing_ccds_filename1'), # this will get an extra extension corresponding to the chromosome name
                               #chromosome => '',
                             },
              -hive_capacity    => 30,
              -analysis_capacity => 30,
            },

            {
              -logic_name => 'update_vega_ccds_labels_concat',
              -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                               'cmd'   => 'cat '.$self->o('reports_dir').'/'.$self->o('vega_missing_ccds_filename1').'.* > '.
                                                 $self->o('reports_dir').'/'.$self->o('vega_missing_ccds_filename1')
                             },
              -flow_into => { 1 => ['create_missing_ccds_report1'] },
              -rc_name => 'default',
            },

            {
              -logic_name => 'create_missing_ccds_report1',
              -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                               'cmd'   => 'grep -Fx -f '.$self->o('reports_dir').'/'.$self->o('ensembl_missing_ccds_filename1').' '.
                                                 $self->o('reports_dir').'/'.$self->o('vega_missing_ccds_filename1').' > '.
                                                 $self->o('reports_dir').'/'.$self->o('missing_ccds_filename1')
                             },
              -flow_into => { 1 => ['email_missing_ccds_report1'] },
              -rc_name => 'default',
            },

            {
              -logic_name => 'email_missing_ccds_report1',
              -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::TextfileByEmail',
              -parameters => {
                               email => $self->o('ccds_report_email'),
                               subject => 'AUTOMATED REPORT: missing CCDS before copying missing CCDS',
                               text => 'Please find below the list of missing CCDS models before copying the missing CCDS into the Ensembl database.',
                               file => $self->o('reports_dir').'/'.$self->o('missing_ccds_filename1'),
                               command => q{cat},
                             },
              -flow_into => { 1 => ['copy_missing_ccds'] },
              -rc_name => 'default',
            },

            {
              -logic_name => 'copy_missing_ccds',
              -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCopyGenes',
              -parameters => {
                               copy_genes_path => $self->o('ensembl_analysis_base').'/scripts/genebuild/',
                               copy_genes_script_name => 'copy_genes.pl',

                               # copy_genes.pl script parameters
                               sourcehost => $self->o('ccds_db','-host'),
                               sourceuser => $self->o('ccds_db','-user'),
                               sourceport => $self->o('ccds_db','-port'),
                               sourcepass => $self->o('ccds_db','-pass'),
                               sourcedbname => $self->o('ccds_db','-dbname'),
                               outhost => $self->o('ensembl_db','-host'),
                               outuser => $self->o('ensembl_db','-user'),
                               outpass => $self->o('ensembl_db','-pass'),
                               outdbname => $self->o('ensembl_db','-dbname'),
                               outport => $self->o('ensembl_db','-port'),
                               dnahost => $self->o('ensembl_db','-host'),
                               dnadbname => $self->o('ensembl_db','-dbname'),
                               dnauser => $self->o('user_r'),
                               dnaport => $self->o('ensembl_db','-port'),
                               logic => 'ensembl',
                               source => 'ensembl',
                               biotype => 'protein_coding',
                               stable_id => 1,
                               merge => 1,
                               file => $self->o('reports_dir').'/'.$self->o('missing_ccds_filename1'),
                             },
               -flow_into => { 1 => ['list_ensembl_toplevel_for_update_ccds_labels_after'] },
               -analysis_capacity => 20,
               -hive_capacity => 20,
               -max_retry_count => 0,
            },

            {
              -logic_name => 'list_ensembl_toplevel_for_update_ccds_labels_after',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
              -parameters => {
                               db_conn => $self->o('ensembl_db'),
                               inputquery => 'SELECT sr.name FROM seq_region sr, seq_region_attrib sra, attrib_type at WHERE sr.seq_region_id = sra.seq_region_id AND sr.name NOT LIKE "LRG\_%" AND sra.attrib_type_id = at.attrib_type_id AND code = "toplevel"',
                               column_names => ['chromosome'],
                             },
              -flow_into => { '2->A' => [ 'update_ensembl_ccds_labels_after' ],
                              'A->1' => [ 'update_ensembl_ccds_labels_concat_after' ],
                            },
              -rc_name => 'default',
            },

            {
              -logic_name => 'update_ensembl_ccds_labels_after',
              -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveUpdateCCDSLabels',
              -parameters => {
                               ccds_dbname => $self->o('ccds_db','-dbname'),
                               ccds_host => $self->o('ccds_db','-host'),
                               ccds_user => $self->o('ccds_db','-user'),
                               output_dbname => $self->o('ensembl_db','-dbname'),
                               output_host => $self->o('ensembl_db','-host'),
                               output_user => $self->o('user_w'),
                               output_pass => $self->o('pass_w'),
                               dna_dbname => $self->o('ensembl_db','-dbname'),
                               dna_host => $self->o('ensembl_db','-host'),
                               dna_user => $self->o('user_r'),
                               dna_pass => $self->o('pass_r'),
                               assembly_path => $self->o('assembly_path'),
                               reports_dir => $self->o('reports_dir'),
                               output_filename => $self->o('ensembl_missing_ccds_filename2'), # this will get an extra extension corresponding to the chromosome name
                               #chromosome => '',
                             },
              -hive_capacity    => 30,
              -analysis_capacity => 30,
            },
            {
              -logic_name => 'update_ensembl_ccds_labels_concat_after',
              -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                               'cmd'   => 'cat '.$self->o('reports_dir').'/'.$self->o('ensembl_missing_ccds_filename2').'.* > '.
                                                 $self->o('reports_dir').'/'.$self->o('ensembl_missing_ccds_filename2')
                             },
              -flow_into => { 1 => ['list_vega_toplevel_for_update_ccds_labels_after'] },
              -rc_name => 'default',
            },

            {
              -logic_name => 'list_vega_toplevel_for_update_ccds_labels_after',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
              -parameters => {
                               db_conn => $self->o('vega_db'),
                               inputquery => 'SELECT sr.name FROM seq_region sr, seq_region_attrib sra, attrib_type at WHERE sr.seq_region_id = sra.seq_region_id AND sr.name NOT LIKE "LRG\_%" AND sra.attrib_type_id = at.attrib_type_id AND code = "toplevel"',
                               column_names => ['chromosome'],
                             },
              -flow_into => { '2->A' => [ 'update_vega_ccds_labels_after' ],
                              'A->1' => [ 'update_vega_ccds_labels_concat_after' ],
                            },
              -rc_name => 'default',
            },

            {
              -logic_name => 'update_vega_ccds_labels_after',
              -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveUpdateCCDSLabels',
              -parameters => {
                               ccds_dbname => $self->o('ccds_db','-dbname'),
                               ccds_host => $self->o('ccds_db','-host'),
                               ccds_user => $self->o('ccds_db','-user'),
                               output_dbname => $self->o('vega_db','-dbname'),
                               output_host => $self->o('vega_db','-host'),
                               output_user => $self->o('user_w'),
                               output_pass => $self->o('pass_w'),
                               dna_dbname => $self->o('ensembl_db','-dbname'),
                               dna_host => $self->o('ensembl_db','-host'),
                               dna_user => $self->o('user_r'),
                               dna_pass => $self->o('pass_r'),
                               assembly_path => $self->o('assembly_path'),
                               reports_dir => $self->o('reports_dir'),
                               output_filename => $self->o('vega_missing_ccds_filename2'), # this will get an extra extension corresponding to the chromosome name
                               #chromosome => '',
                             },
              -hive_capacity    => 30,
              -analysis_capacity => 30,
            },

            {
              -logic_name => 'update_vega_ccds_labels_concat_after',
              -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                               'cmd'   => 'cat '.$self->o('reports_dir').'/'.$self->o('vega_missing_ccds_filename2').'.* > '.
                                                 $self->o('reports_dir').'/'.$self->o('vega_missing_ccds_filename2')
                             },
              -flow_into => { 1 => ['create_missing_ccds_report2'] },
              -rc_name => 'default',
            },

            {
              -logic_name => 'create_missing_ccds_report2',
              -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                               'cmd'   => 'grep -Fx -f '.$self->o('reports_dir').'/'.$self->o('ensembl_missing_ccds_filename2').' '.
                                                 $self->o('reports_dir').'/'.$self->o('vega_missing_ccds_filename2').' | cat > '.
                                                 $self->o('reports_dir').'/'.$self->o('missing_ccds_filename2')
                             },
              -flow_into => { 1 => ['email_missing_ccds_report2'] },
              -rc_name => 'default',
            },

            {
              -logic_name => 'email_missing_ccds_report2',
              -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::TextfileByEmail',
              -parameters => {
                               email => $self->o('ccds_report_email'),
                               subject => 'AUTOMATED REPORT: missing CCDS after copying missing CCDS',
                               text => 'Please find below the list of missing CCDS models after copying the missing CCDS into the Ensembl database. NOTE THIS LIST SHOULD BE EMPTY. OTHERWISE THERE IS A PROBLEM.',
                               file => $self->o('reports_dir').'/'.$self->o('missing_ccds_filename2'),
                               command => q{cat},
                             },
              -flow_into => { 1 => ['list_vega_genes_for_merge'] },
              -rc_name => 'default',
            },

            {
              -logic_name => 'list_vega_genes_for_merge',
              -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                               'cmd'   => 'mysql -NB -u'.$self->o('user_r').
                                               ' -h'.$self->o('vega_db','-host').
                                               ' -D'.$self->o('vega_db','-dbname').
                                               ' -P'.$self->o('vega_db','-port').
                                               ' -e"SELECT gene_id from gene;" > '.
                                               $self->o('output_dir').'/'.$self->o('vega_genes_for_merge_filename')
                             },
              -flow_into => { 1 => ['chunk_vega_genes_for_merge'] },
              -rc_name => 'default',
            },

            {
              -logic_name => 'chunk_vega_genes_for_merge',
              -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::FileFactory',
              -parameters => {
                               inputfile => $self->o('output_dir').'/'.$self->o('vega_genes_for_merge_filename'),
                               output_dir => $self->o('output_dir'),
                               output_prefix => $self->o('vega_genes_for_merge_filename')."_chunk_",
                             },
              -flow_into => { '2->A' => ['havana_merge'],
                              'A->1' => ['havana_merge_list_processed_genes'],
                            },
              -rc_name => 'default',
            },

            {
              -logic_name => 'havana_merge',
              -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMerge',
              -parameters => {
                               ensembl_analysis_base => $self->o('ensembl_analysis_base'),
                               host_secondary => $self->o('ensembl_db','-host'),
                               user_secondary => $self->o('user_r'),
                               password_secondary => $self->o('pass_r'),
                               database_secondary => $self->o('ensembl_db','-dbname'),
                               host_primary => $self->o('vega_db','-host'),
                               user_primary => $self->o('user_r'),
                               password_primary =>$self->o('pass_r'),
                               database_primary => $self->o('vega_db','-dbname'),
                               host_output => $self->o('core_db','-host'),
                               user_output =>$self->o('user_w'),
                               password_output => $self->o('pass_w'),
                               database_output => $self->o('core_db','-dbname'),
                               secondary_include => '',
                               secondary_exclude => '',
                               primary_include => '',
                               primary_exclude => '',

                               # Tagging:  Will be used as suffix for logic names ("_tag") and for
                               # source.  With the default settings, merged genes and transcripts will
                               # get the source "secondary_primary".

                               secondary_tag => 'ensembl',
                               primary_tag => 'havana',

                               # Xrefs:  The format is a comma-separated list of
                               # "db_name,db_display_name,type"

                               primary_gene_xref => 'OTTG,Havana gene,ALT_GENE',
                               primary_transcript_xref => 'OTTT,Havana transcript,ALT_TRANS',
                               primary_translation_xref => 'OTTP,Havana translation,MISC',

                               # as the chunks (and a job per chunk) are created in the step before,
                               # these parameters would define how many jobs per chunk we want, just 1 as we don't want chunks of chunks
                               # and we cannot use the LSF job index on the ehive to create chunks of chunks here anyway
                               njobs => 1, #$self->o('njobs'),
                               job => 1,   #$LSB_JOBINDEX

                               #file => '', this parameter will come from 'chunk_genes_for_merge' output, see FileFactory.pm
                             },
              -rc_name => 'normal_1500',
              -analysis_capacity => $self->o('njobs'),
              -hive_capacity => $self->o('njobs'),
              -max_retry_count => 0,
              -wait_for => ['core_sql_truncates'],
            },

            {
              -logic_name => 'havana_merge_list_processed_genes',
              -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                               'cmd'   => "awk '\$1 == ".'"'."PROCESSED".'"'." {print \$2}' ".$self->o('output_dir')."/*merge-run*.out ".
                               " | sort -u -n > ".$self->o('output_dir').'/'.$self->o('processed_genes_filename')
                             },
              -rc_name => 'default',
              -flow_into => { 1 => ['havana_merge_list_unprocessed_genes'] },
            },

            {
              -logic_name => 'havana_merge_list_unprocessed_genes',
              -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveListUnprocessed',
              -parameters => {
              	               processed_genes_filename => $self->o('processed_genes_filename'),
              	               output_dir => $self->o('output_dir'),
              	               output_file => $self->o('unprocessed_genes_filename'),
                               host_secondary => $self->o('ensembl_db','-host'),
                               user_secondary => $self->o('user_r'),
                               password_secondary => $self->o('pass_r'),
                               database_secondary => $self->o('ensembl_db','-dbname'),
                               secondary_include => '',
                               secondary_exclude => '',
                             },
              -rc_name => 'default',
              #-hive_capacity    => 100,
              -flow_into => { 1 => ['chunk_unprocessed_genes'] },
            },

            {
              -logic_name => 'chunk_unprocessed_genes',
              -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::FileFactory',
              -parameters => {
                               inputfile => $self->o('output_dir').'/'.$self->o('unprocessed_genes_filename'),
                               output_dir => $self->o('output_dir'),
                               output_prefix => $self->o('unprocessed_genes_filename')."_chunk_",
                             },
              -flow_into => { '2->A' => [ 'copy_unprocessed_genes' ],
                              'A->1' => [ 'set_ncrna' ],
                            },
              -rc_name => 'default',
            },

            {
              -logic_name => 'copy_unprocessed_genes',
              -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCopyGenes',
              -parameters => {
                               copy_genes_path => $self->o('ensembl_analysis_base').'/scripts/genebuild/',
                               copy_genes_script_name => 'copy_genes.pl',

                               # copy_genes.pl script parameters
                               logic => 'ensembl',
                               sourcehost => $self->o('ensembl_db','-host'),
                               sourceuser => $self->o('ensembl_db','-user'),
                               sourceport => $self->o('ensembl_db','-port'),
                               sourcepass => $self->o('ensembl_db','-pass'),
                               sourcedbname => $self->o('ensembl_db','-dbname'),
                               outhost => $self->o('core_db','-host'),
                               outuser => $self->o('core_db','-user'),
                               outpass => $self->o('core_db','-pass'),
                               outdbname => $self->o('core_db','-dbname'),
                               outport => $self->o('core_db','-port'),
                               dnahost => $self->o('ensembl_db','-host'),
                               dnadbname => $self->o('ensembl_db','-dbname'),
                               dnauser => $self->o('user_r'),
                               dnaport => $self->o('ensembl_db','-port'),
                               #file => $self->o('output_dir').$self->o('unprocessed_genes_filename'),
                             },
               -analysis_capacity => 20,
               -hive_capacity => 20,
               -max_retry_count => 0,
            },

            {
              -logic_name => 'set_ncrna',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
              -parameters => {
                               db_conn => $self->o('core_db'),
                               sql => [ 'INSERT IGNORE analysis(logic_name) VALUES("ncrna")',
                                        'UPDATE gene g,analysis a SET g.analysis_id=(SELECT analysis_id FROM analysis where logic_name="ncrna")
                                                                  WHERE g.analysis_id=a.analysis_id
                                                                         AND a.logic_name="ensembl"
                                                                         AND g.biotype in (
                                                                                           "miRNA",
                                                                                           "misc_RNA",
                                                                                           "ribozyme",
                                                                                           "rRNA",
                                                                                           "scaRNA",
                                                                                           "snoRNA",
                                                                                           "snRNA",
                                                                                           "sRNA"
                                                                                          )',
                                        'UPDATE transcript t,analysis a SET t.analysis_id=(SELECT analysis_id FROM analysis where logic_name="ncrna")
                                                                   WHERE t.analysis_id=a.analysis_id
                                                                         AND a.logic_name="ensembl"
                                                                         AND t.biotype in (
                                                                                           "miRNA",
                                                                                           "misc_RNA",
                                                                                           "ribozyme",
                                                                                           "rRNA",
                                                                                           "scaRNA",
                                                                                           "snoRNA",
                                                                                           "snRNA",
                                                                                           "sRNA"
                                                                                          )'
                                       ],
                             },
              -flow_into => { 1 => ['set_igtr_analysis_biotypes'] },
            },

            {
              -logic_name => 'set_igtr_analysis_biotypes',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
              -parameters => {
                               db_conn => $self->o('core_db'),
                               sql => [ 'UPDATE gene g SET g.biotype=REPLACE(g.biotype,"_","_V_")
                                                                     WHERE (g.description LIKE "%variable%" OR
                                                                            g.description LIKE "%variant%" OR
                                                                            g.description LIKE "%ighv%" OR
                                                                            g.description LIKE "%immunog.%V.%pseudogene%")
                                                                           AND
                                                                           (g.biotype LIKE "IG\_%" OR
                                                                            g.biotype LIKE "TR\_%")
                                         ',
                                         'UPDATE gene g SET g.biotype=REPLACE(g.biotype,"_","_C_")
                                                                     WHERE (g.description LIKE "%constant%")
                                                                           AND
                                                                           (g.biotype LIKE "IG\_%" OR
                                                                            g.biotype LIKE "TR\_%")
                                         ',
                                         'UPDATE gene g SET g.biotype=REPLACE(g.biotype,"_","_J_")
                                                                     WHERE (g.description LIKE "%joining%")
                                                                           AND
                                                                           (g.description NOT LIKE "%constant%")
                                                                           AND
                                                                           (g.biotype LIKE "IG\_%" OR
                                                                            g.biotype LIKE "TR\_%")
                                         ',
                                         'UPDATE gene g SET g.biotype=REPLACE(g.biotype,"_","_D_")
                                                                     WHERE (g.description LIKE "%diversity%")
                                                                           AND
                                                                           (g.biotype LIKE "IG\_%" OR
                                                                            g.biotype LIKE "TR\_%")
                                         ',
                                         'UPDATE transcript t,gene g SET t.biotype=g.biotype
                                                                     WHERE t.gene_id=g.gene_id
                                                                           AND
                                                                           (g.biotype LIKE "IG\_%" OR
                                                                            g.biotype LIKE "TR\_%")
                                                                           AND
                                                                           (t.biotype LIKE "IG\_%" OR
                                                                            t.biotype LIKE "TR\_%")
                                         ',
                                         'UPDATE gene g SET analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name="ensembl_ig_gene")
                                                      WHERE (g.biotype LIKE "IG\_%" OR
                                                             g.biotype LIKE "TR\_%")
                                                            AND
                                                            analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name="ensembl")
                                         ',
                                         'UPDATE transcript SET analysis_id=(select analysis_id FROM analysis WHERE logic_name="ensembl_ig_gene")
                                                            WHERE (biotype LIKE "IG\_%" OR
                                                                   biotype LIKE "TR\_%")
                                                                  AND
                                                                  analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name="ensembl")
                                         ',
                                         'UPDATE gene g SET analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name="havana_ig_gene")
                                                      WHERE (g.biotype LIKE "IG\_%" OR
                                                             g.biotype LIKE "TR\_%")
                                                            AND
                                                            analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name="havana")
                                         ',
                                         'UPDATE transcript SET analysis_id=(select analysis_id FROM analysis WHERE logic_name="havana_ig_gene")
                                                            WHERE (biotype LIKE "IG\_%" OR
                                                                   biotype LIKE "TR\_%")
                                                                  AND
                                                                  analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name="havana")
                                         ',
                                         'UPDATE gene g SET analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name="ensembl_havana_ig_gene")
                                                      WHERE (g.biotype LIKE "IG\_%" OR
                                                             g.biotype LIKE "TR\_%")
                                                            AND
                                                            analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name="ensembl_havana")
                                         ',
                                         'UPDATE transcript SET analysis_id=(select analysis_id FROM analysis WHERE logic_name="ensembl_havana_ig_gene")
                                                            WHERE (biotype LIKE "IG\_%" OR
                                                                   biotype LIKE "TR\_%")
                                                                  AND
                                                                  analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name="ensembl_havana")
                                         ',
                                         'UPDATE gene g SET analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name="proj_ensembl_ig_gene")
                                                      WHERE (g.biotype LIKE "IG\_%" OR
                                                             g.biotype LIKE "TR\_%")
                                                            AND
                                                            analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name="proj_ensembl")
                                         ',
                                         'UPDATE transcript SET analysis_id=(select analysis_id FROM analysis WHERE logic_name="proj_ensembl_ig_gene")
                                                            WHERE (biotype LIKE "IG\_%" OR
                                                                   biotype LIKE "TR\_%")
                                                                  AND
                                                                  analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name="proj_ensembl")
                                         ',
                                         'UPDATE gene g SET analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name="proj_havana_ig_gene")
                                                      WHERE (g.biotype LIKE "IG\_%" OR
                                                             g.biotype LIKE "TR\_%")
                                                            AND
                                                            analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name="proj_havana")
                                         ',
                                         'UPDATE transcript SET analysis_id=(select analysis_id FROM analysis WHERE logic_name="proj_havana_ig_gene")
                                                            WHERE (biotype LIKE "IG\_%" OR
                                                                   biotype LIKE "TR\_%")
                                                                  AND
                                                                  analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name="proj_havana")
                                         ',
                                         'UPDATE gene g SET analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name="proj_ensembl_havana_ig_gene")
                                                      WHERE (g.biotype LIKE "IG\_%" OR
                                                             g.biotype LIKE "TR\_%")
                                                            AND
                                                            analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name="proj_ensembl_havana")
                                         ',
                                         'UPDATE transcript SET analysis_id=(select analysis_id FROM analysis WHERE logic_name="proj_ensembl_havana_ig_gene")
                                                            WHERE (biotype LIKE "IG\_%" OR
                                                                   biotype LIKE "TR\_%")
                                                                  AND
                                                                  analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name="proj_ensembl_havana")
                                         '
                                       ],
                             },
              -max_retry_count => 0,
              -flow_into => { 1 => ['list_toplevel_for_vega_checks_after'] },
            },

            {

              -logic_name => 'list_toplevel_for_vega_checks_after',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
              -parameters => {
                               db_conn => $self->o('core_db'),
                               inputquery => 'SELECT sr.name FROM seq_region sr, seq_region_attrib sra, attrib_type at WHERE sr.seq_region_id = sra.seq_region_id AND sr.name NOT LIKE "LRG\_%" AND sra.attrib_type_id = at.attrib_type_id AND code = "toplevel"',
                               column_names => ['chromosome'],
                             },
              -flow_into => { '2->A' => [ 'vega_checks_after' ],
                              'A->1' => [ 'vega_checks_after_concat' ],
                            },
              -rc_name => 'default',
            },

            {
              -logic_name => 'vega_checks_after_concat',
              -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                               'cmd'   => 'cat '.$self->o('reports_dir').'/vega_checks_after_*.out > '.
                                                 $self->o('reports_dir').'/vega_checks_after.out'
                             },
              -flow_into => { 1 => ['vega_checks_after_report'] },
              -rc_name => 'default',
            },

            {
              -logic_name => 'vega_checks_after',
              -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveVegaChecks',
              -parameters => {
                               dbname => $self->o('core_db','-dbname'),
                               dbhost => $self->o('core_db','-host'),
                               dnadbname => $self->o('ensembl_db','-dbname'),
                               dnadbhost => $self->o('ensembl_db','-host'),
                               coord_system => 'toplevel',
                               path => $self->o('assembly_path'),
                               sql_output => $self->o('output_dir').'/vega_checks_after_#chromosome#.sql',
                               dbtype => '', # can be 'vega' or '' (empty string)
                               port => '3306',
                               user => $self->o('user_w'),
                               pass => $self->o('pass_w'),
                               #chromosome => '',
                               write => 1,
                               affix => 1, # perform the checks by using the biotypes with or without the prefixes and suffixes like weird_, _Ens, _hav, ... ; with affixes by default
                               biotypes_extension => 1,
                               stdout_file => $self->o('reports_dir').'/vega_checks_after_#chromosome#.out',
                               stderr_file => $self->o('reports_dir').'/vega_checks_after_#chromosome#.err',
                             },
              -hive_capacity    => 30,
              -analysis_capacity => 30,
            },

            {
              -logic_name => 'vega_checks_after_report',
              -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::TextfileByEmail',
              -parameters => {
                               email => $self->o('vega_checks_reports_email'),
                               subject => 'AUTOMATED REPORT: merge biotype combinations',
                               text => 'Please find below the list of not allowed gene and transcript biotype combinations AFTER the merge found in the core database '.$self->o('core_db','-dbname').'. Any artifact transcript listed below will be deleted. Please note that any change applied to the list of allowed biotype combinations will take effect in a future release (not the next release).',
                               file => $self->o('reports_dir').'/vega_checks_after.out',
                               command => q{egrep '(Unknown)|(not allowed\.)' | awk '{print $9,$18}' | sort | uniq -c | sort -nr | sed s'/\. run/ (UNKNOWN gene biotype)/g'},
                             },
              -flow_into => { 1 => ['delete_artifacts'] },
              -rc_name => 'default',
            },
            {
              -logic_name => 'delete_artifacts',
              -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDeleteTranscripts',
              -parameters => {
              	               biotype => 'artifact',
                               dbhost => $self->o('core_db','-host'),
                               dbname => $self->o('core_db','-dbname'),
                               dbuser => $self->o('user_w'),
                               dbpass => $self->o('pass_w'),
                               dbport => $self->o('core_db','-port'),
                               delete_transcripts_path => $self->o('ensembl_analysis_base').'/scripts/genebuild/',
                               delete_genes_path => $self->o('ensembl_analysis_base').'/scripts/genebuild/',
                               delete_transcripts_script_name => 'delete_transcripts.pl',
                               delete_genes_script_name => 'delete_genes.pl',
                               output_path => $self->o('output_dir').'/delete_artifacts/',
                               output_file_name => 'delete_artifacts.out',
                               email => $self->o('vega_checks_reports_email'),
                               from => 'ensembl-genebuild@ebi.ac.uk'
                             },
              -max_retry_count => 0,
              -flow_into => { 1 => ['set_temp_stable_ids'] },
            },

            {
              # make sure that all my transcript have stable IDs so that the CCDS can be added as supporting features later on
              -logic_name => 'set_temp_stable_ids',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
              -parameters => {
                               db_conn => $self->o('core_db'),
                               sql => [ "UPDATE transcript SET stable_id=CONCAT('TEMPSID',transcript_id) WHERE stable_id IS NULL" ],
                             },
               -max_retry_count => 3,
               -rc_name => 'default',
               -flow_into => { 1 => ['list_toplevel'] },
            },


            {
              -logic_name => 'list_toplevel',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
              -parameters => {
                               db_conn => $self->o('core_db'),
                               inputquery => 'SELECT sr.name FROM seq_region sr, seq_region_attrib sra, attrib_type at WHERE sr.seq_region_id = sra.seq_region_id AND sr.name NOT LIKE "LRG\_%" AND sra.attrib_type_id = at.attrib_type_id AND code = "toplevel"',
                               column_names => ['chr'],
                             },
               -flow_into => { '2->A' => ['alternative_atg_attributes'],

                              'A->1' => WHEN ('#process_ccds#' => ['ccds_sql_updates'],
                                        ELSE                      ['set_frameshift_transcript_attributes']),
                            },
              -rc_name => 'default',
            },

            {
              -logic_name => 'alternative_atg_attributes',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                               cmd => 'perl '.$self->o('ensembl_analysis_base').'/scripts/Merge/look_for_upstream_atg.pl'.
                                      ' -dbuser '.$self->o('user_w').
                                      ' -dbpass '.$self->o('pass_w').
                                      ' -dbhost '.$self->o('core_db','-host').
                                      ' -dbport '.$self->o('core_db','-port').
                                      ' -dbname '.$self->o('core_db','-dbname').
                                      ' -dna_user '.$self->o('user_r').
                                      ' -dna_host '.$self->o('core_db','-host').
                                      ' -dna_port '.$self->o('core_db','-port').
                                      ' -dna_dbname '.$self->o('core_db','-dbname').
                                      ' -upstream_dist 200'.
                                      ' -chromosomes #chr#'.
                                      ' -genetypes protein_coding'.
                                      ' -coord_system toplevel'.
                                      ' -path '.$self->o('assembly_path')
                             },
               -analysis_capacity => 25,
               -hive_capacity => 25,
               -max_retry_count => 2,
               -rc_name => 'normal_1500',
               -flow_into => { 1 => WHEN ('!#process_ccds#' => ['set_canonical_transcripts_without_ccds'])},
            },

            {
              -logic_name => 'set_canonical_transcripts_without_ccds',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                               cmd => 'perl '.$self->o('ensembl_analysis_base').'/../ensembl/misc-scripts/canonical_transcripts/select_canonical_transcripts.pl'.
                                      ' -dbuser '.$self->o('user_w').
                                      ' -dbpass '.$self->o('pass_w').
                                      ' -dbhost '.$self->o('core_db','-host').
                                      ' -dbport '.$self->o('core_db','-port').
                                      ' -dbname '.$self->o('core_db','-dbname').
                                      ' -seq_region_name #chr#'.
                                      ' -coord toplevel'.
                                      ' -write'
                             },
               -analysis_capacity => 25,
               -hive_capacity => 25,
               -max_retry_count => 2,
               -rc_name => 'normal_1500',
            },

            {
              -logic_name => 'ccds_sql_updates',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
              -parameters => {
                               db_conn => $self->o('core_db'),
                               sql => [ 'UPDATE gene SET biotype="protein_coding" WHERE biotype="ccds_gene"',
                                        'UPDATE transcript SET biotype="protein_coding" WHERE biotype="ccds_gene"',
                                        'INSERT IGNORE analysis(created,logic_name) VALUES(now(),"ccds")',
                                        'UPDATE dna_align_feature SET analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name="ccds")
                                                                WHERE analysis_id=(SELECT analysis_id FROM analysis WHERE logic_name="ccds_gene")
                                                                      OR dna_align_feature.hit_name LIKE "ccds%"',
                                        'UPDATE gene SET source="ensembl" WHERE source="ccds"',
                                        'DELETE FROM analysis WHERE logic_name="ccds_gene"',
                                       ],
                             },
              -max_retry_count => 0,
              -flow_into => { 1 => ['prepare_lincrnas'] },
            },

            {
              -logic_name => 'prepare_lincrnas',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
              -parameters => {
                               db_conn => $self->o('core_db'),
                               sql => [ "UPDATE gene SET biotype = 'new_lincRNA' WHERE biotype = 'lincRNA'" ],
                             },
               -max_retry_count => 3,
               -rc_name => 'default',
               -flow_into => { 1 => ['transfer_lincrnas'] },
            },

            {
              -logic_name => 'transfer_lincrnas',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                               cmd => 'perl '.$self->o('ensembl_analysis_base').'/scripts/Merge/transfer_lincRNAs_to_merged_gene_set.pl'.
                                      ' -dbname '.$self->o('prevcore_db','-dbname').
                                      ' -dbhost '.$self->o('prevcore_db','-host').
                                      ' -dbport '.$self->o('prevcore_db','-port').
                                      ' -dbuser '.$self->o('user_r').
                                      ' -newdbname '.$self->o('core_db','-dbname').
                                      ' -newdbhost '.$self->o('core_db','-host').
                                      ' -newdbuser '.$self->o('user_w').
                                      ' -newdbpass '.$self->o('pass_w').
                                      ' -vegadbname '.$self->o('vega_db','-dbname').
                                      ' -vegadbhost '.$self->o('vega_db','-host').
                                      ' -vegadbport '.$self->o('vega_db','-port').
                                      ' -vegadbuser '.$self->o('user_r').
                                      ' -coordsystem chromosome'.
                                      ' -path '.$self->o('assembly_path').
                                      ' -write '.
                                      ' -verbose'
                             },
               -max_retry_count => 0,
               -rc_name => 'normal_1500',
               -flow_into => { 1 => ['set_lincrna_biotypes'] },
            },

            {
              -logic_name => 'set_lincrna_biotypes',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
              -parameters => {
                               db_conn => $self->o('core_db'),
                               sql => [ "UPDATE gene SET biotype = 'lincRNA' WHERE biotype = 'new_lincRNA'" ],
                             },
               -max_retry_count => 3,
               -rc_name => 'default',
               -flow_into => { 1 => ['backup_dump_core_db_after_merge'] },
            },

            {
              -logic_name => 'set_frameshift_transcript_attributes',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                               cmd => 'perl '.$self->o('ensembl_analysis_base').'/../ensembl/misc-scripts/frameshift_transcript_attribs.pl'.
                                      ' -dbuser '.$self->o('user_w').
                                      ' -dbpass '.$self->o('pass_w').
                                      ' -dbhost '.$self->o('core_db','-host').
                                      ' -dbport '.$self->o('core_db','-port').
                                      ' -dbpattern '.$self->o('core_db','-dbname')
                             },
               -analysis_capacity => 25,
               -hive_capacity => 25,
               -max_retry_count => 2,
               -rc_name => 'normal_4600',
               -flow_into => { 1 => ['set_repeat_types_after_merge'] },
            },

            {
              -logic_name => 'set_repeat_types_after_merge',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                               cmd => 'perl '.$self->o('ensembl_analysis_base').'/../ensembl/misc-scripts/repeats/repeat-types.pl'.
                                      ' -user '.$self->o('user_w').
                                      ' -pass '.$self->o('pass_w').
                                      ' -host '.$self->o('core_db','-host').
                                      ' -port '.$self->o('core_db','-port').
                                      ' -dbpattern '.$self->o('core_db','-dbname')
                             },
               -analysis_capacity => 25,
               -hive_capacity => 25,
               -max_retry_count => 2,
               -rc_name => 'default',
               -flow_into => { 1 => ['load_external_db_ids_and_optimise_af_tables'] },
            },

            {
              -logic_name => 'load_external_db_ids_and_optimise_af_tables',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                               cmd => 'perl '.$self->o('ensembl_analysis_base').'/scripts/genebuild/load_external_db_ids_and_optimize_af.pl'.
                                      ' -output_path '.$self->o('output_dir').'/optimise'.
                                      ' -uniprot_filename '.$self->o('uniprot_file').
                                      ' -dbuser '.$self->o('user_w').
                                      ' -dbpass '.$self->o('pass_w').
                                      ' -dbhost '.$self->o('core_db','-host').
                                      ' -dbname '.$self->o('core_db','-dbname').
                                      ' -prod_dbuser '.$self->o('user_r').
                                      ' -prod_dbhost '.$self->o('production_db','-host').
                                      ' -prod_dbname '.$self->o('production_db','-dbname').
                                      ' -prod_dbport '.$self->o('production_db','-port').
                                      ' -verbose'
                             },
               -analysis_capacity => 25,
               -hive_capacity => 25,
               -max_retry_count => 2,
               -rc_name => 'normal_12000',
            },
            
            {
              -logic_name => 'backup_dump_core_db_after_merge',
              -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
              -parameters => {
                                create_type => 'backup',
                                source_db => $self->o('core_db'),
                                backup_name => $self->o('backup_dump_core_db_after_merge_filename'),
                                pass_w => $self->o('pass_w'),
                                user_w => $self->o('user_w'),
                                output_path => $self->o('output_dir'),
                             },
              -flow_into => { 1 => WHEN ('#patch_update#' => ['stop_for_sid_hcs_and_gencode_qc'])},
            },
            
            {
              -logic_name => 'stop_for_sid_hcs_and_gencode_qc',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                               'cmd'   => 'echo "Stop"'
                             },
              -rc_name => 'default',
              -wait_for => ['wait_forever'],
              -flow_into => { 1 => ['load_patches'] },
            },
            
            {
              -logic_name => 'wait_forever',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                               'cmd'   => 'echo "Wait forever"'
                             },
              -rc_name => 'default',
            },
            
            {
              # after gencode qc, comment all the above analyses and uncomment the following line, init the pipeline again and resume beekeeper
              #-input_ids => [ {} ],
              -logic_name => 'load_patches',
              -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadPatches',
              -parameters => {
                                dbhost => $self->o('core_db','-host'),
                                dbname => $self->o('core_db','-dbname'),
                                dbuser => $self->o('user_w'),
                                dbpass => $self->o('pass_w'),
                                ftp_path => $self->o('patches_ftp_dir'),
                                output_path => $self->o('output_dir')."/load_patches/",
                                cs_version => $self->o('assembly_path'),
                             },
              -max_retry_count => 0,
              -rc_name => 'normal_1500',
              -flow_into => { 1 => ['update_meta_coords'] },
            },
            
            {
               -logic_name => 'update_meta_coords',
               -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
               -parameters => {
                                 cmd => 'perl $ENSCODE/ensembl/misc-scripts/meta_coord/update_meta_coord.pl'.
                                      ' -user '.$self->o('user_w').
                                      ' -pass '.$self->o('pass_w').
                                      ' -host '.$self->o('core_db','-host').
                                      ' -port '.$self->o('core_db','-port').
                                      ' -dbpattern '.$self->o('core_db','-dbname')
                              },
               -rc_name => 'default',
               -flow_into => { 1 => ['list_chr_for_nt_transform_check'] },
            },
            
            {
              -logic_name => 'list_chr_for_nt_transform_check',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
              -parameters => {
                               db_conn => $self->o('core_db'),
                               inputquery => 'SELECT sr.name FROM seq_region sr, seq_region_attrib sra, attrib_type at, coord_system cs WHERE cs.coord_system_id=sr.coord_system_id AND sr.seq_region_id = sra.seq_region_id AND sr.name NOT LIKE "LRG\_%" AND sra.attrib_type_id = at.attrib_type_id AND code = "toplevel" AND sr.name LIKE "CHR%" and cs.version="'.$self->o('assembly_path').'"',
                               column_names => ['chromosome_name'],
                             },
              -flow_into => { '2->A' => [ 'nt_transform_check' ],
                              'A->1' => [ 'start_patch_annotation' ],
                            },
              -rc_name => 'default',
            },

            {
              -logic_name => 'nt_transform_check',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                               cmd => 'perl '.$self->o('ensembl_analysis_base').'/scripts/assembly_patches/nt_transform_check.pl'.
                                      ' -dbname '.$self->o('core_db','-dbname').
                                      ' -h '.$self->o('core_db','-host').
                                      ' -port '.$self->o('core_db','-port').
                                      ' -u '.$self->o('user_r').
                                      ' -central_coordsystem scaffold'.
                                      ' -path '.$self->o('assembly_path').
                                      ' -chromosome_name #chromosome_name#'.
                                      ' > '.$self->o('output_dir').'/nt_transform_check_multi_#chromosome_name#.out'
                             },
               -max_retry_count => 0,
               -rc_name => 'basement_1500',
            },

            {
               -logic_name => 'start_patch_annotation',
               -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
               -parameters => {},
               -rc_name       => 'default',
               -can_be_empty  => 1,
               -flow_into => { '1->A' => ['create_patch_slice_ids_1'],
                               'A->1' => ['set_repeat_types'],
                             },
            },

            # Create 1mb top-level slices
            {
              -logic_name => 'create_patch_slice_ids_1',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
              -parameters => {
                               target_db        => $self->o('core_db'),
                               coord_system_name => 'toplevel',
                               iid_type => 'patch_slice',
                               top_level => 1,
                             },
              -flow_into => {'2' => ['repeatmasker']},
            },

###############################################################################
#
# REPEATMASKER ANALYSES
#
###############################################################################

            {
              -logic_name => 'repeatmasker',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveRepeatMasker',
              -parameters => {
                               target_db => $self->o('core_db'),
                               logic_name => 'repeatmask_repbase_'.$self->o('repeatmasker_library'),
                               module => 'HiveRepeatMasker',
                               repeatmasker_path => $self->o('repeatmasker_path'),
                               commandline_params => '-nolow -species "'.$self->o('repeatmasker_library').'" -engine "'.$self->o('repeatmasker_engine').'"',
                             },
              -rc_name    => 'normal_2900',
              -flow_into => {
                            -1 => ['repeatmasker_himem'],
                              1 => ['dust'],
                            },
              -hive_capacity => 900,
            },

            {
              -logic_name => 'repeatmasker_himem',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveRepeatMasker',
              -parameters => {
                               target_db => $self->o('core_db'),
                               logic_name => 'repeatmask_repbase_'.$self->o('repeatmasker_library'),
                               module => 'HiveRepeatMasker',
                               repeatmasker_path => $self->o('repeatmasker_path'),
                               commandline_params => '-nolow -species "'.$self->o('repeatmasker_library').'" -engine "'.$self->o('repeatmasker_engine').'"',
                             },
              -rc_name    => 'normal_5900',
              -flow_into => {1 => ['dust']},
              -hive_capacity => 900,
              -can_be_empty  => 1,
            },

            {
              -logic_name => 'dump_softmasked_toplevel',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDumpGenome',
              -parameters => {
                               'coord_system_name'    => 'toplevel',
                               'target_db'            => $self->o('core_db'),
                               'output_path'          => $self->o('output_dir')."/genome_dumps/",
                               'enscode_root_dir'     => "$ENSCODE",
                               'species_name'         => $self->o('species_name'),
                               'repeat_logic_names'   => $self->o('repeat_logic_names'),
                               'patch_only'           => 1,
                             },
              -input_ids => [{}],
              -wait_for => ['dust','repeatmasker','repeatmasker_himem'],
              -rc_name    => 'default',
            },

###############################################################################
#
# SIMPLE FEATURE AND OTHER REPEAT ANALYSES
#
###############################################################################

            {
              -logic_name => 'dust',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveDust',
              -parameters => {
                               target_db => $self->o('core_db'),
                               logic_name => 'dust',
                               module => 'HiveDust',
                               dust_path => $self->o('dust_path'),
                             },
              -rc_name    => 'normal_2900',
              -flow_into => {1 => ['trf']},
              -hive_capacity => 900,
              -batch_size => 20,
            },

            {
              -logic_name => 'trf',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveTRF',
              -parameters => {
                               target_db => $self->o('core_db'),
                               logic_name => 'trf',
                               module => 'HiveTRF',
                               trf_path => $self->o('trf_path'),
                             },
              -rc_name    => 'normal_2900',
              -flow_into => {
                               1 => ['eponine'],
                            },
              -hive_capacity => 900,
              -batch_size => 20,
            },

            {
              -logic_name => 'eponine',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveEponine',
              -parameters => {
                               target_db => $self->o('core_db'),
                               logic_name => 'eponine',
                               module => 'HiveEponine',
                               eponine_path => $self->o('eponine_java_path'),
                               commandline_params => '-epojar => '.$self->o('eponine_jar_path').', -threshold => 0.999',
                             },
              -rc_name    => 'normal_2900',
              -flow_into => {1 => ['firstef']},
              -hive_capacity => 900,
              -batch_size => 20,
            },

            {
              -logic_name => 'firstef',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveFirstEF',
              -parameters => {
                               target_db => $self->o('core_db'),
                               logic_name => 'firstef',
                               module => 'HiveFirstEF',
                               firstef_path => $self->o('firstef_path'),
                               repeat_masking_logic_names => ['repeatmask_repbase_'.$self->o('repeatmasker_library')],
                               commandline_params => '-repeatmasked',
                             },
              -rc_name    => 'normal_2900',
              -flow_into => {1 => ['cpg']},
              -hive_capacity => 900,
              -batch_size => 20,
            },

            {
              -logic_name => 'cpg',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveCPG',
              -parameters => {
                               target_db => $self->o('core_db'),
                               logic_name => 'cpg',
                               module => 'HiveCPG',
                               cpg_path => $self->o('cpg_path'),
                             },
              -rc_name    => 'normal_2900',
              -flow_into => {1 => ['trnascan']},
              -hive_capacity => 900,
              -batch_size => 20,
            },

            {
              -logic_name => 'trnascan',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveTRNAScan',
              -parameters => {
                               target_db => $self->o('core_db'),
                               logic_name => 'trnascan',
                               module => 'HiveTRNAScan',
                               trnascan_path => $self->o('trnascan_path'),
                             },
              -rc_name    => 'normal_2900',
              -flow_into => {1 => ['genscan']},
              -hive_capacity => 900,
              -batch_size => 20,
            },


###############################################################################
#
# GENSCAN ANALYSIS
#
##############################################################################

# Run genscan, uses 1mb slices from repeatmasker. Flows into create_prediction_transcript_ids which
# then takes these 1mb slices and converts them into individual prediction transcript input ids based
# on the dbID of each feature generate by this analysis

            {
              -logic_name => 'genscan',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveGenscan',
              -parameters => {
                               target_db => $self->o('core_db'),
                               logic_name => 'genscan',
                               module => 'HiveGenscan',
                               genscan_path => $self->o('genscan_path'),
                               genscan_matrix_path => $self->o('genscan_matrix_path'),
                               repeat_masking_logic_names => ['repeatmask_repbase_'.$self->o('repeatmasker_library')],
                             },
              -rc_name    => 'normal_2900_180',
              -flow_into => {
                              1 => ['create_prediction_transcript_ids'],
                              -1 => ['decrease_genscan_slice_size'],
                              -2 => ['decrease_genscan_slice_size'],
                            },
              -hive_capacity => 900,
              -batch_size => 20,
            },

            # Create 1mb top-level slices
            {
              -logic_name => 'decrease_genscan_slice_size',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
              -parameters => {
                               iid_type => 'split_slice',
                               slice_size => 100000,
                             },
              -flow_into => {2 => ['genscan_short_slice']},
              -rc_name    => 'default',
              -can_be_empty  => 1,
            },

            {
              -logic_name => 'genscan_short_slice',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveGenscan',
              -parameters => {
                               target_db => $self->o('core_db'),
                               logic_name => 'genscan',
                               module => 'HiveGenscan',
                               genscan_path => $self->o('genscan_path'),
                               repeat_masking_logic_names => ['repeatmask_repbase_'.$self->o('repeatmasker_library')],
                             },
              -rc_name    => 'normal_2900_180',
              -flow_into => {
                              1 => ['create_prediction_transcript_ids'],
                              -1 => ['failed_genscan_slices'],
                              -2 => ['failed_genscan_slices'],
                            },
              -rc_name    => 'normal_5900_60',
              -can_be_empty  => 1,
              -hive_capacity => 900,
            },

            {
              -logic_name => 'failed_genscan_slices',
              -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
              -parameters => {},
              -rc_name          => 'default',
              -can_be_empty  => 1,
            },


            # Create input ids for individual prediction transcripts. Takes a slice as an input id and converts it
            # to a set of input ids that are individual dbIDs for the prediction transcripts. This avoids empty slices
            # being submitted as jobs and also means one feature corresponds to one job. Each species flows into this
            # independently with 1mb slices. Also has the advantage that downstream analyses can start working as soon
            # as a single slice is finished
            {
              -logic_name => 'create_prediction_transcript_ids',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
              -parameters => {
                               target_db => $self->o('core_db'),
                               feature_type => 'prediction_transcript',
                               iid_type => 'slice_to_feature_ids',
                               prediction_transcript_logic_names => ['genscan'],
                             },
              -flow_into => {2 => ['uniprot_blast']},
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
              -logic_name => 'uniprot_blast',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlastGenscanPep',
              -parameters => {
                               target_db => $self->o('core_db'),
                               blast_db_path => $self->o('uniprot_blast_db_path'),
                               blast_exe_path => $self->o('uniprot_blast_exe_path'),
                               commandline_params => '-cpus 3 -hitdist 40',
                               repeat_masking_logic_names => ['repeatmask_repbase_'.$self->o('repeatmasker_library')],
                               prediction_transcript_logic_names => ['genscan'],
                               iid_type => 'feature_id',
                               logic_name => 'uniprot',
                               module => 'HiveBlastGenscanPep',
                               %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastGenscanPep')},
                               #config_settings => $self->get_config_settings('HiveBlast','HiveBlastGenscanPep'),
                            },
              -flow_into => {
                              1 => ['vertrna_blast'],
                              -1 => ['failed_blast_job'],
                              -2 => ['failed_blast_job'],
                            },
              -rc_name    => 'normal_2900_n3',
              -failed_job_tolerance => 0.5,
              -hive_capacity => 900,
              -batch_size => 20,
            },

            {
              # BLAST individual prediction transcripts against vertRNA. The config settings are held lower in this
              # file in the master_config_settings sub
              -logic_name => 'vertrna_blast',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlastGenscanDNA',
              -parameters => {
                               target_db => $self->o('core_db'),
                               blast_db_path => $self->o('vertrna_blast_db_path'),
                               blast_exe_path => $self->o('vertrna_blast_exe_path'),
                               commandline_params => '-cpus 3 -hitdist 40',
                               repeat_masking_logic_names => ['repeatmask_repbase_'.$self->o('repeatmasker_library')],
                               prediction_transcript_logic_names => ['genscan'],
                               iid_type => 'feature_id',
                               logic_name => 'vertrna',
                               module => 'HiveBlastGenscanDNA',
                               %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastGenscanVertRNA')},
                               #config_settings => $self->get_config_settings('HiveBlast','HiveBlastGenscanVertRNA'),
                            },
              -flow_into => {
                              1 => ['unigene_blast'],
                              -1 => ['failed_blast_job'],
                              -2 => ['failed_blast_job'],
                            },
              -rc_name    => 'normal_2900_n3',
             -failed_job_tolerance => 0.5,
             -hive_capacity => 900,
             -batch_size => 20,
            },

            {
              # BLAST individual prediction transcripts against unigene. The config settings are held lower in this
              # file in the master_config_settings sub
              -logic_name => 'unigene_blast',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlastGenscanDNA',
              -parameters => {
                               target_db => $self->o('core_db'),
                               blast_db_path => $self->o('unigene_blast_db_path'),
                               blast_exe_path => $self->o('unigene_blast_exe_path'),
                               commandline_params => '-cpus 3 -hitdist 40',
                               prediction_transcript_logic_names => ['genscan'],
                               iid_type => 'feature_id',
                               repeat_masking_logic_names => ['repeatmask_repbase_'.$self->o('repeatmasker_library')],
                               logic_name => 'unigene',
                               module => 'HiveBlastGenscanDNA',
                               %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastGenscanUnigene')},
                               #config_settings => $self->get_config_settings('HiveBlast','HiveBlastGenscanUnigene'),
                            },
              -flow_into => {
                              -1 => ['failed_blast_job'],
                              -2 => ['failed_blast_job'],
                            },
              -rc_name    => 'normal_2900_n3',
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
              -logic_name => 'set_repeat_types',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                               cmd => 'perl $ENSCODE/ensembl/misc-scripts/repeats/repeat-types.pl'.
                                      ' -user '.$self->o('user_w').
                                      ' -pass '.$self->o('pass_w').
                                      ' -host '.$self->o('core_db','-host').
                                      ' -port '.$self->o('core_db','-port').
                                      ' -dbpattern '.$self->o('core_db','-dbname')
                             },
               -rc_name => 'default',
               -flow_into => { 1 => ['set_meta_coords'] },
            },

           {
              -logic_name => 'set_meta_coords',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                               cmd => 'perl $ENSCODE/ensembl/misc-scripts/meta_coord/update_meta_coord.pl'.
                                      ' -user '.$self->o('user_w').
                                      ' -pass '.$self->o('pass_w').
                                      ' -host '.$self->o('core_db','-host').
                                      ' -port '.$self->o('core_db','-port').
                                      ' -dbpattern '.$self->o('core_db','-dbname')
                             },
               -rc_name => 'default',
               -flow_into => { 1 => ['set_meta_levels'] },
            },

            {
              -logic_name => 'set_meta_levels',
              -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
              -parameters => {
                               cmd => 'perl $ENSCODE/ensembl/misc-scripts/meta_levels.pl'.
                                      ' -user '.$self->o('user_w').
                                      ' -pass '.$self->o('pass_w').
                                      ' -host '.$self->o('core_db','-host').
                                      ' -port '.$self->o('core_db','-port').
                                      ' -dbname '.$self->o('core_db','-dbname')
                             },
               -rc_name => 'default',
               -flow_into => { 1 => ['create_genblast_output_db'] },
            },

            {
              -logic_name => 'create_genblast_output_db',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
              -parameters => {
                               source_db => $self->o('ro_core_db'),
                               target_db => $self->o('genblast_db'),
                               create_type => 'clone',
                               script_path => $self->o('clone_db_script_path'),
                             },
              -rc_name    => 'default',
              -max_retry_count => 0,
              -flow_into => {
                              '1->A' => ['download_uniprot_files'],
                              'A->1' => ['classify_genblast_models'],
                            },
            },

            {
              -logic_name => 'download_uniprot_files',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadUniProtFiles',

              -parameters => {
                               multi_query_download => $self->uniprot_clade_download(),

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
                               uniprot_db_path => $self->o('homology_models_path').'/'.$self->o('uniprot_db_name'),
                               uniprot_index_path => $self->o('homology_models_path').'/'.$self->o('uniprot_index_name'),
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
              -rc_name      => 'normal_2900',
              -flow_into => {
                              2 => ['genblast'],
                            },

            },

            {
              -logic_name => 'genblast',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast',
              -parameters => {
                               iid_type => 'db_seq',
                               dna_db => $self->o('core_db'),
                               target_db => $self->o('genblast_db'),
                               logic_name => 'genblast',
                               module => 'HiveGenblast',
                               genblast_path => $self->o('genblast_path'),
                               genblast_db_path => $self->o('genome_file'),
                               commandline_params => ' -P wublast -gff -e '.$self->o('genblast_eval').' -c '.$self->o('genblast_cov').' ',
                               query_seq_dir => $self->o('homology_models_path').'/'.$self->o('uniprot_query_dir_name'),
                               sequence_table_name => $self->o('uniprot_table_name'),
                               max_rank => $self->o('genblast_max_rank'),
                               genblast_pid => $self->o('genblast_pid'),
                             },
              -rc_name    => 'normal_1900_120',
              -flow_into => {
                              -1 => ['split_genblast_jobs'],
                              -2 => ['split_genblast_jobs'],
                              -3 => ['split_genblast_jobs'],
                            },
              -failed_job_tolerance => 0.5,
-wait_for => ['dump_softmasked_toplevel'],
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
                               dna_db => $self->o('core_db'),
                               target_db => $self->o('genblast_db'),
                               logic_name => 'genblast',
                               module => 'HiveGenblast',
                               genblast_path => $self->o('genblast_path'),
                               genblast_db_path => $self->o('genome_file'),
                               commandline_params => ' -P wublast -gff -e '.$self->o('genblast_eval').' -c '.$self->o('genblast_cov').' ',
                               query_seq_dir => $self->o('homology_models_path').'/'.$self->o('uniprot_query_dir_name'),
                               sequence_table_name => $self->o('uniprot_table_name'),
                               max_rank => $self->o('genblast_max_rank'),
                               genblast_pid => $self->o('genblast_pid'),
                             },
              -rc_name          => 'normal_4900_120',
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
              -parameters => {},
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
                              1 => ['create_exonerate_output_db'],
                            },

            },

            {
              -logic_name => 'create_exonerate_output_db',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
              -parameters => {
                               source_db => $self->o('ro_core_db'),
                               target_db => $self->o('exonerate_db'),
                               create_type => 'clone',
                               script_path => $self->o('clone_db_script_path'),
                               user_r => $self->o('user_r'),
                               user_w => $self->o('user_w'),
                               pass_w => $self->o('pass_w'),
                             },
              -rc_name    => 'default',
              -flow_into => {
                              '1->A' => ['create_exonerate_ids'],
                              'A->1' => ['classify_exonerate_models'],
                            },
            },

            {
              -logic_name => 'create_exonerate_ids',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
              -parameters => {
                               target_db => $self->o('genblast_db'),
                               iid_type => 'feature_id',
                               feature_type => 'transcript',
                               feature_logic_names => ['genblast','genblast_not_best'],
                             },
              -flow_into => {
                              2 => ['exonerate'],
                            },
              -rc_name    => 'default',
            },

            {
              -logic_name => 'exonerate',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes',
              -rc_name    => 'normal_1900_120',
              -parameters => {
                               sequence_table_name => $self->o('uniprot_table_name'),
                               iid_type => 'feature_id',
                               feature_type => 'transcript',
                               transcript_db => $self->o('genblast_db'),
                               region_padding => $self->o('exonerate_region_padding'),
                               use_genblast_best_in_genome => 0,
                               dna_db => $self->o('core_db'),
                               target_db => $self->o('exonerate_db'),
                               logic_name => 'exonerate',
                               module     => 'HiveExonerate2Genes',
                               config_settings => $self->get_config_settings('exonerate_protein','exonerate'),
                               query_seq_dir => $self->o('homology_models_path').'/'.$self->o('uniprot_query_dir_name'),
                               calculate_coverage_and_pid => $self->o('exonerate_calculate_coverage_and_pid'),
                            },
            },

            {
              -logic_name => 'classify_exonerate_models',
              -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClassifyTranscriptSupport',
              -parameters => {
                               classification_type => 'standard',
                               update_gene_biotype => 1,
                               target_db => $self->o('exonerate_db'),
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
                         source_db => $self->o('ro_core_db'),
                         target_db => $self->o('layering_db'),
                         create_type => 'clone',
                         script_path => $self->o('clone_db_script_path'),
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
                         source_db => $self->o('ro_core_db'),
                         target_db => $self->o('utr_db'),
                         create_type => 'clone',
                         script_path => $self->o('clone_db_script_path'),
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
                         source_db => $self->o('ro_core_db'),
                         target_db => $self->o('genebuilder_db'),
                         create_type => 'clone',
                         script_path => $self->o('clone_db_script_path'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '1->A' => ['create_patch_slice_ids_3'],
                        'A->1' => ['create_pseudogene_db'],
                      },
      },

      {
        -logic_name => 'create_patch_slice_ids_3',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('core_db'),
                         iid_type => 'patch_slice',
                         coord_system_name => 'toplevel',
                         top_level => 1,
                         # These options will create only slices that have a gene on the slice in one of the feature dbs
                         feature_constraint => 1,
                         feature_type => 'gene',
                         feature_dbs => [$self->o('genblast_db')],
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
                         dna_db     => $self->o('core_db'),
                         logic_name => 'layer_annotation',
                         module     => 'HiveLayerAnnotation',
                         config_settings => $self->get_config_settings('layer_annotation','layers'),
                       },
        -rc_name    => 'normal_3900_200',
        -flow_into  => {
                         '1->A' => ['split_slices_on_intergenic'],
                         'A->1' => ['genebuilder'],
                       },
      },

      {
        -logic_name => 'split_slices_on_intergenic',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveFindIntergenicRegions',
        -parameters => {
                         dna_db => $self->o('core_db'),
                         input_gene_dbs => $self->o('utr_gene_dbs'),
                         iid_type => 'slice',
                       },
        -batch_size => 100,
        -hive_capacity => 200,
        -rc_name    => 'normal_1900',
        -flow_into => {
                        2 => ['cluster_input_genes'],
                      },
      },


      {
        -logic_name => 'cluster_input_genes',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveClusterSourceGenes',
        -parameters => {
                         logic_name => 'cluster_input_genes',
                         dna_db => $self->o('ro_core_db'),
                         input_gene_dbs => $self->o('utr_gene_dbs'),
                         allowed_input_sets => undef,
                         iid_type => 'slice',
                       },
        -batch_size => 20,
        -hive_capacity => 200,
        -rc_name    => 'normal_1900',
        -flow_into => {
                        2 => ['run_utr_addition'],
                      },

      },

      {
        -logic_name => 'run_utr_addition',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveUTRAddition',
        -parameters => {
                         logic_name => 'utr_addition',
                         dna_db => $self->o('ro_core_db'),
                         input_gene_dbs => $self->default_options()->{'utr_gene_dbs'},
                         utr_biotype_priorities => $self->o('utr_biotype_priorities'),
                         utr_output_db => $self->default_options()->{'utr_db'},
                         iid_type => 'slice',
                       },
        -batch_size => 20,
        -hive_capacity => 900,
        -rc_name    => 'normal_1900',
     },

     {
        -logic_name => 'genebuilder',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGeneBuilder',
        -parameters => {
                         layering_output_db => $self->o('utr_db'),
                         genebuilder_output_db => $self->o('genebuilder_db'),
                         dna_db     => $self->o('ro_core_db'),
                         logic_name => 'ensembl',
                         module     => 'HiveGeneBuilder',
                         config_settings => $self->get_config_settings('genebuilder','genebuilder_set'),
                       },
        -rc_name    => 'normal_1900',
        -hive_capacity => 900,
     },

     {
        -logic_name => 'create_pseudogene_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('ro_core_db'),
                         target_db => $self->o('pseudogene_db'),
                         create_type => 'clone',
                         script_path => $self->o('clone_db_script_path'),
                       },
        -rc_name    => 'default',
        -flow_into => {
                        '1->A' => ['create_pseudogene_ids'],
                        'A->1' => ['create_patch_geneset_db'],
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
                         repeat_db => $self->default_options()->{'ro_core_db'},
                         output_db => $self->o('pseudogene_db'),
                         dna_db => $self->o('ro_core_db'),
                         logic_name => 'pseudogenes',
                         module     => 'HivePseudogenes',
                         config_settings => $self->get_config_settings('pseudogenes','pseudogenes_set'),
                       },
        -batch_size => 20,
        -hive_capacity => 900,
        -rc_name    => 'normal_1900',
        -flow_into => {
                       2 => ['spliced_elsewhere'],
                      },
     },

     {
        -logic_name => 'concat_blast_db',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
                         cmd => 'for i in '.$self->o('output_dir').'/pseudogenes/multi_exon_dir/multi_exon_seq*.fasta;'.
                                'do cat $i >> '.$self->o('output_dir').'/pseudogenes/multi_exon_dir/all_multi_exon_genes.fasta;'.
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
                         cmd => 'xdformat -n '.$self->o('output_dir').'/pseudogenes/multi_exon_dir/all_multi_exon_genes.fasta'
                       },
         -rc_name => 'default',
     },

     {
        -logic_name => 'spliced_elsewhere',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSplicedElsewhere',
        -parameters => {
                         input_gene_db => $self->o('genebuilder_db'),
                         repeat_db => $self->o('core_db'),
                         output_db => $self->o('pseudogene_db'),
                         dna_db => $self->o('core_db'),
                         logic_name => 'spliced_elsewhere',
                         module     => 'HiveSplicedElsewhere',
                         config_settings => $self->get_config_settings('pseudogenes','pseudogenes_set'),
                       },
        -rc_name          => 'normal_1900',
        -can_be_empty  => 1,
        -wait_for => ['format_blast_db'],
     },

     {
        -logic_name => 'create_patch_geneset_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('pseudogene_db'),
                         target_db => $self->o('patch_geneset_db'),
                         create_type => 'copy',
                         script_path => $self->o('clone_db_script_path'),
                         user_r => $self->o('user_r'),
                         user_w => $self->o('user_w'),
                         pass_w => $self->o('pass_w'),
                       },
        -rc_name    => 'default',
     },
  ];
}


sub pipeline_wide_parameters {
    my ($self) = @_;

      return {
            # Inherit other stuff from the parent class
                %{$self->SUPER::pipeline_wide_parameters()},
                'process_ccds' => $self->o('process_ccds'),
      };
}

sub uniprot_clade_download {
  my ($self) = @_;

  my $clade = $self->default_options()->{'uniprot_set'};
  my $output_path = $self->o('homology_models_path');

  my $taxon_ids = {
                   'human_taxon_id'    => '9606',
                   'mouse_taxon_id'    => '10090',
                   'rodents_taxon_id'  => '9989',
                 };

  if ($clade eq 'human_patch') {
    return ({
              human_pe12 => {
                              file_name => 'human_pe12.fasta',
                              taxon_id  => $taxon_ids->{'human_taxon_id'},
                              dest_dir  => $output_path,
                              compress  => 0,
                              pe_level  => [1,2],
                            },
                            
              human_frag_pe12 => {
                                  file_name  => 'human_frag_pe12.fasta',
                                  taxon_id   => $taxon_ids->{'human_taxon_id'},
                                  dest_dir   => $output_path,
                                  compress   => 0,
                                  fragment  => 1,
                                  pe_level   => [1,2],
                                 },
             });
  } elsif ($clade eq 'mouse_patch') {
    return ({
              mouse_pe12 => {
                              file_name => 'mouse_pe12.fasta',
                              taxon_id  => $taxon_ids->{'mouse_taxon_id'},
                              dest_dir  => $output_path,
                              compress  => 0,
                              pe_level  => [1,2],
                            },
                            
              mouse_frag_pe12 => {
                                  file_name  => 'mouse_frag_pe12.fasta',
                                  taxon_id   => $taxon_ids->{'mouse_taxon_id'},
                                  dest_dir   => $output_path,
                                  compress   => 0,
                                  fragment  => 1,
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

  if ($clade eq 'human_patch') {
return ([
            {
              ID         => 'LAYER1',
              BIOTYPES   => [
                             'human_pe12_sp_95',
                             'human_pe12_tr_95',
                            ],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER2',
              BIOTYPES   => [
                              'human_frag_pe12_sp_95',
                              'human_frag_pe12_tr_95',
                            ],
              FILTER_AGAINST => ['LAYER1'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER3',
              BIOTYPES   => [
                             'human_pe12_sp_80',
                             'human_pe12_tr_80',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2'],
              DISCARD    => 0,
            },
            
            {
              ID         => 'LAYER4',
              BIOTYPES   => [
                             'human_frag_pe12_sp_80',
                             'human_frag_pe12_tr_80',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'],
              DISCARD    => 0,
            },
    ]);

  } elsif ($clade eq 'mouse_patch') {
return ([
            {
              ID         => 'LAYER1',
              BIOTYPES   => [
                             'mouse_pe12_sp_95',
                             'mouse_pe12_tr_95',
                            ],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER2',
              BIOTYPES   => [
                              'mouse_frag_pe12_sp_95',
                              'mouse_frag_pe12_tr_95',
                            ],
              FILTER_AGAINST => ['LAYER1'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER3',
              BIOTYPES   => [
                             'mouse_pe12_sp_80',
                             'mouse_pe12_tr_80',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2'],
              DISCARD    => 0,
            },
            
            {
              ID         => 'LAYER4',
              BIOTYPES   => [
                             'mouse_frag_pe12_sp_80',
                             'mouse_frag_pe12_tr_80',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'],
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

  if ($clade eq 'human_patch') {
    return ([
              'human_pe12_sp_95',
              'human_pe12_sp_80',
              'human_pe12_tr_95',
              'human_pe12_tr_80',
              'human_frag_pe12_sp_95',
              'human_frag_pe12_tr_95',
              'human_frag_pe12_sp_80',
              'human_frag_pe12_tr_80',
    ]);
  } elsif ($clade eq 'mouse_patch') {
    return ([
              'mouse_pe12_sp_95',
              'mouse_pe12_sp_80',
              'mouse_pe12_tr_95',
              'mouse_pe12_tr_80',
              'mouse_frag_pe12_sp_95',
              'mouse_frag_pe12_tr_95',
              'mouse_frag_pe12_sp_80',
              'mouse_frag_pe12_tr_80',
    ]);
  } else {
    die "Unknown clade selected for genebuilding: ".$clade;
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

  exonerate_protein => {
    Default => {
                 IIDREGEXP           => '(\d+):(\d+)',
                 OPTIONS             => '--model protein2genome --forwardcoordinates FALSE --softmasktarget TRUE --exhaustive FALSE --bestn 1 --maxintron 50000 --minintron 20',
                 COVERAGE_BY_ALIGNED => 0,
                 QUERYTYPE           => 'protein',
                 GENOMICSEQS         => $self->o('genome_file'),
                 PROGRAM             => $self->o('exonerate_path'),
                 SOFT_MASKED_REPEATS => $self->o('repeat_masking_logic_names'),
               },

    exonerate => {
                   FILTER     
                                      => {
                     OBJECT                      => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
                     PARAMETERS                  => {
                       -coverage                 => $self->o('exonerate_cov'),
                       -percent_id               => $self->o('exonerate_pid'),
                       -best_in_genome           => 1,
                       -reject_processed_pseudos => 1,
                     },
                   },
                 },

    killlist_protein => {
                         KILLLISTDB          => $self->o('killlist_db'),
                         USE_KILL_LIST       => 1,
                         KILL_TYPE           => 'protein',
                         KILL_LIST_FILTER    => {
                                                  -only_mol_type        => 'protein',
                                                  -user_id              => undef,
                                                  -from_source_species  => undef,
                                                  -before_date          => undef,
                                                  -having_status        => undef,
                                                  -reasons              => [],
                                                  -for_analyses         => [],
                                                  -for_species          => [],
                                                  -for_external_db_ids  => [],
                                                },
                       },
    },

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

    layer_annotation => {
      Default => {
        TARGETDB_REF => $self->o('layering_db'),
        SOURCEDB_REFS => $self->o('layering_input_gene_dbs'),
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

    pseudogenes => {
      Default => {
               # you can set the input- and output database - the names should point to
               # keys in Database.pm
               PS_INPUT_DATABASE  => 'GENEBUILD_DB',
               PS_OUTPUT_DATABASE => 'PSEUDO_DB',

               # configs for the introns in repeats test

               # introns longer than the following are considered "real"
               PS_FRAMESHIFT_INTRON_LENGTH => 9,
               # total length of introns
               PS_MAX_INTRON_LENGTH   => '5000',
               # Types of repeats to run the anaysis with
               PS_REPEAT_TYPES =>  ['LINE','LTR','SINE'],
               # max percent coverage of the introns with the above repeats
               PS_MAX_INTRON_COVERAGE => '80',
               # max allowed exon coverage with the above repeats
               PS_MAX_EXON_COVERAGE   => '99',
               PS_NUM_FRAMESHIFT_INTRONS  => 1,
               PS_NUM_REAL_INTRONS  => 1,
               # biotype of genes to check
               PS_BIOTYPE  => 'protein_coding',

               # Blessed genes dont get called pseudogenes
               # Biotype is a transcript biotype
               BLESSED_BIOTYPES => { 'ccds_gene' => 1 },

               # configs for the spliced elsewhere tests
               # %ID of a tbalstx of the (presumed) retrotransposed query sequence to its
               # homolog that is spliced elsewhere in the genome. hits falling below
               # this cutoff are ignored (80%) is suggested
               PS_PERCENT_ID_CUTOFF   => 40,
               PS_P_VALUE_CUTOFF   => '1.0e-50',
               PS_RETOTRANSPOSED_COVERAGE   => 80,
               PS_ALIGNED_GENOMIC  => 100,
               # logic name to give to pseudogenes
               PS_PSEUDO_TYPE      => 'pseudogene',
               # if a gene is found to be a pseudogene, its gene biotype will be changed to
               # PS_PSEUDO_TYPE. By default, the biotype of its transcript will also be changed
               # to PS_PSEUDO_TYPE.  If you want to keep the original transcript biotype
               # instead (so you can keep track of what type of models actually got turned into a
               # pseudogene), set KEEP_TRANS_BIOTYPE to 1.

               KEEP_TRANS_BIOTYPE  => 0,

               # logic name to give genes with exons covered by repeats
               # if left blank they will just get deleted (recommended)
               PS_REPEAT_TYPE      => '',

               # analysis logic names to run over genes falling into these categories
               SINGLE_EXON      => 'spliced_elsewhere',
               INDETERMINATE    => '',
               RETROTRANSPOSED  => '',
               # if you dont wish to run further tests on retro transposed genes
               # What type would you like to give them?
               # previously set to 'retrotransposed', we change to 'processed_pseudogene' from e70 onwards.
               RETRO_TYPE       => 'processed_pseudogene',

               SPLICED_ELSEWHERE_LOGIC_NAME => 'spliced_elsewhere',
               PSILC_LOGIC_NAME => 'Psilc',
               # SPLICED ELSEWHERE SPECIFIC CONFIG
               # ratio of the spans of the retrotransposed gene vs its spliced homologue
               # spliced / retrotransposed
               # ie: 1 is the same length genes
               # many retrotransposed genes have a ratio > 10
               # used to make retrotransposition decision
               PS_SPAN_RATIO          => 3,
               # mimimum number of exons for the spliced gene to have
               PS_MIN_EXONS           => 4,
               # path of blast db of multi exon genes
               PS_MULTI_EXON_DIR       => "/path/to/my/blast/directory/" ,
               # Chunk size
               PS_CHUNK => '50',
               DEBUG => '1',
              },

              pseudogenes_set => {
                PS_INPUT_DATABASE  => 'GBUILD_DB',
                # biotype of genes to check
                PS_BIOTYPE  => 'protein_coding',
                # path of blast db of multi exon genes
                PS_MULTI_EXON_DIR       => $self->o('output_dir').'/pseudogenes/multi_exon_dir/',
              },

              spliced_elsewhere_set => {
                PS_INPUT_DATABASE  => 'GBUILD_DB',
                # biotype of genes to check
                PS_BIOTYPE  => 'ensembl_utr',
                # path of blast db of multi exon genes
                PS_MULTI_EXON_DIR       =>$self->o('output_dir').'/pseudogenes/multi_exon_dir/',
             },

      }, # end pseudogenes
  };

  return($master_config_settings->{$config_group});

}

1;
