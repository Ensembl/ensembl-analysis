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


package Hive_cDNA_update_conf;

use strict;
use warnings;

use File::Spec::Functions qw(catfile);

use Bio::EnsEMBL::ApiVersion qw/software_version/;
use parent ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;           # Allow this particular config to use conditional dataflow and INPUT_PLUS

sub default_options {
  my ($self) = @_;

  return {
    # inherit other stuff from the base class
    %{ $self->SUPER::default_options() },

##########################################################################
#                                                                        #
# CHANGE STUFF HERE                                                      #
#                                                                        #
##########################################################################
    
    'recipient_email'            => $ENV{HIVE_EMAIL}, # email address where reports will be sent

    'species'                    => '', # either mouse or human

    'ensembl_release' => $ENV{ENSEMBL_RELEASE}, # What release are you doing this for?

    'complete_update' => 0, # set this as 1 or 0. Usually we'll probably only be adding alignments from new cDNAs and removing retired ones so you want to set this to 0

    'cdna_batch_size' => 1,

    'pipeline_name' => $self->o('species').'_cdna_update_'.$self->o('ensembl_release'),

    # Database connection info:
    'pipe_dbname' => '',
    'pipe_db_server' => '',
    'pipe_db_port' => '',

    'dna_db_server' => '',
    'dna_dbname' => '',
    'dna_db_port' => '',

    'output_db_server' => '',
    'output_db_name' => '',
    'output_db_port' => '',

    # details of the last cdna db (eg. on livemirror)
    'old_cdna_db_name' => '',
    'output_path' => '', # output directory you want to place downloaded files and log files

    'refseq_path' => $self->o('output_path'), # Where you'd like to store the RefSeq data to be downloaded, usually just the same as output_path
    'refseq_file' => catfile($self->o('refseq_path'), 'refseq_'.$self->o('species').'.fa'), # The name of the refseq file, eg. refseq_mouse.fa

    'genome_file' => '', #The softmasked genome file for either human or mouse

    'repeat_masking_logic_names' => [''], # the repeatmask logic name(s) of the analyses used in the genebuild

    'refseq_version' => '', # look at the refseq home page to get the latest version, this will go into the cDNA database

    'exonerate_path'             => catfile($self->o('software_base'), 'opt', 'exonerate09', 'bin', 'exonerate'),

##########################################################################
#                                                                        #
# MOSTLY STAYS CONSTANT, MOSTLY                                          #
#                                                                        #
##########################################################################

    'default_queue'              => 'production-rh7',

    'exonerate_batch_size'       => '50',

    'exonerate_time_limit'       => '2h',

    'fastasplit_random_path'     => catfile($self->o('binary_base'), 'fastasplit_random'),

    'cdna_file_name'             => 'cdna_update',

    'user'                     => 'ensadmin',
    'password'                   => '',
    'user_r'                     => 'ensro',
    'password_r'                 => undef,

    'cdna_query_dir_name'        => 'cdna_temp',

    'refseq_ftp'                 => 'ftp://ftp.ncbi.nlm.nih.gov/refseq/release/vertebrate_mammalian/',

    'retired_cdnas_file'         => catfile('#wide_output_dir#', 'cdna_update.retired'),
    'genes_to_delete'            => catfile('#wide_output_dir#', 'genes_to_delete.ids'),

    'many_hits_dir'              => 'many_hits',
    'cdna_table_name'            => 'cdna_sequences',
    'cdna_small_batch_size'      => '1',
    'cdna_big_batch_size'        => '5',

    'default_mem'                => '900',
    'exonerate_mem'              => '3900',
    'exonerate_retry_mem'        => '5900',
    'exonerate_high_mem'         => '9900',
    'optimise_mem'               => '11900',
    'download_mem'               => '8000',

    'exonerate_pid'              => '97',
    'exonerate_cov'              => '90',
    'exonerate_version' => '0.9.0', # This is unlikely to change any time soon

    'polyA_script'               => catfile($self->o('enscode_root_dir').'', 'ensembl-analysis', 'scripts', 'cdna_update', 'polyA_clipping.pl'),
    'delete_genes_script'        => catfile($self->o('enscode_root_dir').'', 'ensembl-analysis', 'scripts', 'genebuild', 'delete_genes.pl'),
    'findN_script'               => catfile($self->o('enscode_root_dir').'', 'ensembl-analysis', 'scripts', 'cdna_update', 'find_N.pl'),
    'optimize_script'            => catfile($self->o('enscode_root_dir').'', 'ensembl-analysis', 'scripts', 'genebuild', 'load_external_db_ids_and_optimize_af.pl'),

    'analysis_desc_script'       => catfile($self->o('enscode_root_dir').'', 'ensembl-production', 'scripts', 'production_database', 'populate_analysis_description.pl'),
    'populate_production_script' => catfile($self->o('enscode_root_dir').'', 'ensembl-production', 'scripts', 'production_database', 'populate_production_db_tables.pl'),

    'meta_level_script'          => catfile($self->o('enscode_root_dir').'', 'ensembl', 'misc-scripts', 'meta_levels.pl'),
    'meta_coord_script'          => catfile($self->o('enscode_root_dir').'', 'ensembl', 'misc-scripts', 'meta_coord', 'update_meta_coord.pl'),

    'old_cdna_db_server' => 'mysql-ensembl-mirror',
    'old_cdna_db_port' => '4240',

    'killlist_db_name' => 'gb_kill_list',
    'killlist_db_server' => 'mysql-ens-genebuild-prod-6.ebi.ac.uk',
    'killlist_db_port' => '4532',

    'production_db_name' => 'ensembl_production_'.$self->o('ensembl_release'),
    'production_db_server' => 'mysql-ens-sta-1',
    'production_db_port' => '4519',

    'production_db' => {
      -dbname => $self->o('production_db_name'),
      -host => $self->o('production_db_server'),
      -port => $self->o('production_db_port'),
      -user => $self->o('user_r'),
      -pass => $self->o('password_r'),
    },
 
    'output_db' => {
      -dbname => $self->o('output_db_name'),
      -host => $self->o('output_db_server'),
      -port => $self->o('output_db_port'),
      -user => $self->o('user_w'),
      -pass => $self->o('password'),
      -pass => $self->o('password_r'),
    },

    'killlist_db' => {
      -dbname => $self->o('killlist_db_name'),
      -host => $self->o('killlist_db_server'),
      -port => $self->o('killlist_db_port'),
      -user => $self->o('user_r'),
      -pass => $self->o('password_r'),
    },

    'old_cdna_db' => {
      -dbname => $self->o('old_cdna_db_name'),
      -host => $self->o('old_cdna_db_server'),
      -port => $self->o('old_cdna_db_port'),
      -user => $self->o('user_r'),
      -pass => $self->o('password_r'),
    },
  };
}

sub pipeline_create_commands {
  my ($self) = @_;
  return [
  # inheriting database and hive tables' creation
    @{$self->SUPER::pipeline_create_commands},

    $self->db_cmd('CREATE TABLE '.$self->o('cdna_table_name').' ('.
      'accession varchar(50) NOT NULL,'.
      'seq text NOT NULL,'.
      'biotype varchar(50) NOT NULL,'.
      'PRIMARY KEY (accession))'),
  ];
}

## See diagram for pipeline structure
sub pipeline_analyses {
  my ($self) = @_;

  my %taxon_id = (
    human => 9606,
    mouse => 10090,
  );

  return [
    {
      # need to make sure the database actually copies as if it doesn't the job does appear to complete according to eHive
      -logic_name => 'create_output_db',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db => $self->o('old_cdna_db'),
        target_db => $self->o('output_db'),
        #create_type => WHEN ('#complete_update#' => 'clone', ELSE 'copy'),
        create_type => $self->o('complete_update') == 1 ? 'clone' : 'copy',
        script_path => '#wide_clone_script#',
        pass_w => $self->o('password'),
        user_w => $self->o('user'),
        user_r => $self->o('user_r'),
      },
      -rc_name => 'default',
      -input_ids => [{
        cdna_file => '#wide_output_dir#/'.$self->o('cdna_file_name'),
      }],
      -max_retry_count => 0,
      -flow_into => {
        1 => ['create_output_dir'],
#            1 => [ 'compare_cdna_files'],
#            1 => [ 'download_cdnas'],
      }
    },
    {
      -logic_name => 'copy_tables',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'mysqldump -h '.$self->o('old_cdna_db','-host').' -P '.$self->o('old_cdna_db','-port').' -u '.$self->o('old_cdna_db','-user').' '.$self->o('old_cdna_db','-dbname').' mapping_set karyotype seq_region_synonym > #wide_output_dir#/copied_tables.sql;'.
               'mysql -h '.$self->o('output_db','-host').' -P '.$self->o('output_db','-port').' -u '.$self->o('output_db','-user').' -p'.$self->o('output_db','-pass').' '.$self->o('output_db','-dbname').' < #wide_output_dir#/copied_tables.sql'
      },
      -max_retry_count => 0,
      -flow_into => {
        '1' => [ 'populate_production'],
      }
    },
    {
      -logic_name => 'create_output_dir',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'mkdir -p '.$self->o('output_path')
      },
      -max_retry_count => 0,
      -flow_into => {
        '1' => [ 'copy_tables'],
      }
    },
    {
      -logic_name => 'populate_production',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl #wide_production_script#'.
               ' -dp #wide_output_dir#'.
               ' -d '.$self->o('output_db','-dbname').
               ' -h '.$self->o('output_db','-host').
               ' -u '.$self->o('output_db','-user').
               ' -p '.$self->o('output_db','-pass').
               ' -P '.$self->o('output_db','-port').
               ' -md '.$self->o('production_db','-dbname').
               ' -mh '.$self->o('production_db','-host').
               ' -mu '.$self->o('production_db','-user').
               ' -mP '.$self->o('production_db','-port').
               ' -t external_db -t attrib_type -t misc_set -t unmapped_reason'
      },
      -max_retry_count => 0,
      -flow_into => {
 #     #      '2->A' => [ 'download_cdnas'],
 #     #      'A->1' => [ 'load_cdnas' ],
        1 => [ 'download_cdnas'],
#        1 => ['prepare_cdnas'],
      }
    },
    {
      -logic_name => 'download_cdnas',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadcDNAFiles',
      -parameters => {
        embl_sequences => {
          output_path => '#wide_output_dir#',
          output_file => $self->o('cdna_file_name'),
          species => '#wide_species#',
        },
        refseq_sequences => {
          refseq_ftp => $self->o('refseq_ftp'),
        },
      },
#      -wait_for => ['populate_production'],  
      -max_retry_count => 0,
      -rc_name => 'download',
#      -input_ids => [{}],
#      -flow_into => {
#        1 => ['prepare_cdnas'],
#      },
    },
    {
      -logic_name => 'prepare_cdnas',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivePreparecDNAs',
      -parameters => {
        prepare_seqs => {
          dest_dir => '#wide_output_dir#',
          embl_file => 'embl_' . $taxon_id{$self->o('species')} . '.fa',
          refseq_file => $self->o('refseq_file'),
          #gss_file => '#wide_gss_file#',
          killlist_type => 'cdna_update',
          killlist_db => $self->o('killlist_db'),
          polyA_script => '#wide_polyA#',
          cdna_file => $self->o('cdna_file_name'),
          species => '#wide_species#',
        },
      },
      -input_ids => [{}],
      -wait_for => ['download_cdnas'],
      -max_retry_count => 0,
      -rc_name => 'default',
      -flow_into => {
        1 => WHEN ('#complete_update#'  => 'load_cdnas', ELSE 'load_new_cdnas'),
      },
    },
    {
      # there should probably be a check here to make sure that we get roughly the number of retired sequences we expect
      -logic_name => 'load_new_cdnas',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadUpdatecDNAs',
      -parameters => {
        cdna_file => $self->o('output_path').'/'.$self->o('cdna_file_name').'.clipped',
        species => '#wide_species#',
        old_cdna_db => $self->o('old_cdna_db'),
        new_cdna_db => $self->o('output_db'),
        retire_gene_file => $self->o('genes_to_delete'),
        strategy => 'update', 
      },
      -max_retry_count => 0,
      -rc_name => 'default',
      -flow_into => {
        1 => ['delete_retired_genes'],
        '1->A' => [ 'generate_jobs' ],
        'A->1' => [ 'find_many_hits' ],
      },
    },
    {
      -logic_name => 'delete_retired_genes',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl #wide_delete_genes# -h '.$self->o('output_db','-host').' -u '.$self->o('output_db','-user').' -P '.$self->o('output_db','-port')   
               .' -p '.$self->o('output_db','-pass').' -D '.$self->o('output_db','-dbname').' -idfile '.$self->o('genes_to_delete')
      },
      -max_retry_count => 0,
    },
    {
      -logic_name => 'load_cdnas',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadUpdatecDNAs',
      -parameters => {
        cdna_file => $self->o('output_path').'/'.$self->o('cdna_file_name').'.clipped',
        species => '#wide_species#',
        old_cdna_db => $self->o('old_cdna_db'),
        new_cdna_db => $self->o('output_db'),
        species => '#wide_species#',
        strategy => 'complete',
      },
      -rc_name => 'default',
      #-input_ids => [{}],
      #-wait_for => ['prepare_cdnas'],
      -max_retry_count => 0,
      -flow_into => {
        '1->A' => [ 'generate_jobs' ],
        'A->1' => [ 'find_many_hits' ],
      },
    },
    {
      -logic_name => 'generate_jobs',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        cdna_accession => 1,
        #cdna_batch_size => WHEN ('#complete_update#'  => $self->o('cdna_small_batch_size'), ELSE $self->o('cdna_big_batch_size')),
        cdna_batch_size => $self->o('cdna_batch_size'),
        cdna_table_name => $self->o('cdna_table_name'),
      },
      -rc_name => 'default',
      -flow_into => {
        1 => ['exonerate'],
      },
    },
    {
      -logic_name => 'exonerate',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes_cdna',
      -parameters => {
        iid_type => 'db_seq',
        dna_db => $self->o('dna_db'),
        target_db => $self->o('output_db'),
        logic_name => 'cdna_update',
        timer => $self->o('exonerate_time_limit'),
        module => 'HiveExonerate2Genes_cdna',
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_cdnaupdate')},
        exonerate_pid => $self->o('exonerate_pid'),
        exonerate_cov => $self->o('exonerate_cov'),
        query_seq_dir => $self->o('output_path').'/'.$self->o('cdna_query_dir_name'),
        GENOMICSEQS => '#wide_genome_file#',
        PROGRAM => '#wide_exonerate#',
        SOFT_MASKED_REPEATS => $self->o('repeat_masking_logic_names'),
      },
      -flow_into => {
     #   -1 => ['split_exonerate_jobs'],
     #   -2 => ['split_exonerate_jobs'],
        -1 => ['exonerate_himem'],
        -2 => ['exonerate_second_run'],
        -3 => ['exonerate_second_run'],
        2 => ['exonerate_second_run'],
      },
      -rc_name => 'exonerate',
      -failed_job_tolerance => 50,
      #-batch_size => $self->o('exonerate_batch_size'),
    },
    #{
    #  -logic_name => 'split_exonerate_jobs',
    #  -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
    #  -parameters => {
    #    iid_type => 'rechunk',
    #    cdna_batch_size => 1,
    #    cdna_table_name => $self->o('cdna_table_name'),
    #  },
    #  -rc_name => 'default',
    #  -can_be_empty => 1,
    #  -flow_into => {
    #    1 => ['exonerate_himem'],
    #  },
    #},
    {
      -logic_name => 'exonerate_himem',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes_cdna',
      -parameters => {
        iid_type => 'db_seq',
        dna_db => $self->o('dna_db'),
        target_db => $self->o('output_db'),
        logic_name => 'cdna_update',
        timer => $self->o('exonerate_time_limit'),
        module => 'HiveExonerate2Genes_cdna',
        %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_cdnaupdate_loose')},
        exonerate_pid => $self->o('exonerate_pid'),
        exonerate_cov => $self->o('exonerate_cov'),
        query_seq_dir => '#wide_output_dir#/'.$self->o('cdna_query_dir_name'),
        GENOMICSEQS => '#wide_genome_file#',
        PROGRAM => '#wide_exonerate#',
        SOFT_MASKED_REPEATS => $self->o('repeat_masking_logic_names'),
      },
      -flow_into => {
        2 => ['exonerate_second_run'],
        -1 => ['exonerate_second_run'],
        -2 => ['exonerate_second_run'],
        -3 => ['exonerate_second_run'],
      },
      -rc_name => 'exonerate_himem',
      -can_be_empty => 1,
      -failed_job_tolerance => 100,
    },
    {
      # it may be better to carry out this analysis as a series of sql commands like handover_preparation
      # instead of creating a module. I'll leave it like this for now but I'll come back to this
      -logic_name => 'find_missing_cdnas',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMissingcDNAs',
      -parameters => {
        dest_dir => '#wide_output_dir#',
        query_db => $self->o('output_db'),
        cdna_file => '#wide_output_dir#/'.$self->o('cdna_file_name').'.clipped',
      },
      -rc_name => 'default',
      -max_retry_count => 0,
      -wait_for => ['find_many_hits'],
      -input_ids => [{}],
    },
    {
      -logic_name => 'exonerate_second_run',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes_cdna',
      -parameters => {
        iid_type => 'db_seq',
        dna_db => $self->o('dna_db'),
        target_db => $self->o('output_db'),
        logic_name => 'cdna_update',
        timer => $self->o('exonerate_time_limit'),
        module => 'HiveExonerate2Genes_cdna',
        config_settings => $self->get_config_settings('exonerate_cdna','exonerate_2'),
        #%{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::ExonerateStatic','exonerate_2')},
        query_seq_dir => '#wide_output_dir#/'.$self->o('cdna_query_dir_name'),
        stdout_file => '#wide_output_dir#/exonerate_2.out',
        GENOMICSEQS => '#wide_genome_file#',
        PROGRAM => '#wide_exonerate#',
        SOFT_MASKED_REPEATS => $self->o('repeat_masking_logic_names'),
      },
      -rc_name => 'exonerate_2',
      -can_be_empty => 1,
      -failed_job_tolerance => 100,
      -flow_into => {
        -1 => ['failed_exonerate_jobs'],
        -2 => ['failed_exonerate_jobs'],
        -3 => ['failed_exonerate_jobs'],
      },
    },
    {
      -logic_name => 'failed_exonerate_jobs',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -parameters => {
      },
      -rc_name => 'default',
      -can_be_empty => 1,
      -failed_job_tolerance => 100,
    },
    {
      # it may be better to carry out this analysis as a series of sql commands like handover_preparation
      # instead of creating a module. I'll leave it like this for now but I'll come back to this
      -logic_name => 'find_many_hits',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivecDNAManyHits',
      -parameters => {
        dest_dir => '#wide_output_dir#',
        query_db => $self->o('output_db'),
        file_dir => $self->o('many_hits_dir'),
      },
      -rc_name => 'default',
      -max_retry_count => 0,
      #-input_ids => [{}],
      -flow_into => {
        '1->A' => ['filter_output'],
        'A->1' => ['database_compare'],
      },
    },
    {
      -logic_name => 'filter_output',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'grep -r "rpp" #wide_output_dir#/exonerate_2*out | awk \'{split($0,a,":"); print a[2]}\' >> ' 
               . '#wide_output_dir#/failed_hits.out;'
               . ' grep -r "only" #wide_output_dir#/exonerate_2*out | awk \'{split($0,a,":"); print a[2]}\' >> '
               . '#wide_output_dir#/failed_hits.out;'
               . ' grep -r "reject" #wide_output_dir#/exonerate_2*out | awk \'{split($0,a,":"); print a[2] ": " a[3]}\' >> '
               . '#wide_output_dir#/failed_hits.out;'
               . ' grep -r "max_coverage" #wide_output_dir#/exonerate_2*out | awk \'{split($0,a,":"); print a[1]}\' >> '
               . '#wide_output_dir#/failed_hits.out'
      },
      -max_retry_count => 0,
      -flow_into => {
        1 => [ 'store_unmapped' ],
      },
    },
    {
      -logic_name => 'store_unmapped',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveStoreUnmappedcDNAs',
      -parameters => {
        #gss => '#wide_gss_file#',
        seq_file => '#wide_output_dir#/missing_cdnas.fasta',
        user => $self->o('user'),
        pass => $self->o('password'),
        host => $self->o('output_db_server'),
        port => $self->o('output_db_port'),
        dbname => $self->o('output_db_name'),
        species => '#wide_species#',
        vertrna => '#wide_output_dir#/embl_'.$taxon_id{$self->o('species')}.'.fa',
        refseq => '#wide_refseq_file#',
        infile => '#wide_output_dir#/failed_hits.out',
        findN_prog => '#wide_findN#',
        reasons => '#wide_output_dir#/unmapped_reasons.txt',
        pid => $self->o('exonerate_pid'),
    	cov => $self->o('exonerate_cov'), 
        outdir => '#wide_output_dir#',
     },
      -max_retry_count => 0,
    },
    {
      -logic_name => 'database_compare',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveComparecDNAdbs',
      -parameters => {
        old_cdna_db => $self->o('old_cdna_db'),
        new_cdna_db => $self->o('output_db'),
        output_file => '#wide_output_dir#/comparison.out',
      },
      #-input_ids => [{}],
      -failed_job_tolerance => 0,
      -flow_into => {
        1 => [ 'comparison_report','null_align_features' ],
      },
    },
    {
      -logic_name => 'comparison_report',
      -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::TextfileByEmail',
      -parameters => {
        email => $self->o('recipient_email'),
        subject => 'AUTOMATED REPORT: cDNA update database comparison',
        text => 'Please find below the counts for each toplevel seq_region for the current and the previous cDNA updates:',
        file => '#wide_output_dir#/comparison.out',,
      },
      -failed_job_tolerance => 0,
    },
    {
      -logic_name => 'null_align_features',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => 'mysql://'.$self->o('output_db','-user').':'.
                   $self->o('output_db','-pass').'@'.$self->o('output_db','-host').
                   ':'.$self->o('output_db','-port').'/'.$self->o('output_db','-dbname'),
        sql => [
          "UPDATE dna_align_feature SET external_db_id=NULL",
          "UPDATE protein_align_feature SET external_db_id=NULL",
        ],
      },
      -flow_into => {
        1 => [ 'load_xdbids' ],
      },
      -max_retry_count => 0,
    },
    {
      -logic_name => 'load_xdbids',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl #wide_optimize_script#'.
               ' -output_path #wide_output_dir#/optimise_daf_paf.1'.
               ' -dbname '.$self->o('output_db','-dbname').
               ' -dbhost '.$self->o('output_db','-host').
               ' -dbport '.$self->o('output_db','-port').
               ' -dbuser '.$self->o('output_db','-user').
               ' -dbpass '.$self->o('output_db','-pass').
               ' -prod_dbuser '.$self->o('production_db','-dbname').
               ' -prod_dbpass '.$self->o('production_db','-host').
               ' -prod_dbhost '.$self->o('production_db','-port').
               ' -prod_dbname '.$self->o('production_db','-user').
               ' -verbose -clean -no_external_db'
      },
      -rc_name => 'optimise',
      -max_retry_count => 0,
      -flow_into => {
        1 => [ 'populate_analysis_description' ],
      },
    },
    {
      -logic_name => 'populate_analysis_description',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl #wide_analysis_desc#'.
               ' -dp #wide_output_dir#'.
               ' -d '.$self->o('output_db','-dbname').
               ' -h '.$self->o('output_db','-host').
               ' -u '.$self->o('output_db','-user').
               ' -p '.$self->o('output_db','-pass').
               ' -P '.$self->o('output_db','-port').
               ' -s #wide_species# -t cdna'
      },
      -flow_into => {
        1 => [ 'handover_preparation' ],
      },
      -max_retry_count => 0,
    },
    {
      -logic_name => 'handover_preparation',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => 'mysql://'.$self->o('output_db','-user').':'.
                   $self->o('output_db','-pass').'@'.$self->o('output_db','-host').
                   ':'.$self->o('output_db','-port').'/'.$self->o('output_db','-dbname'),
        sql => [
          "DELETE FROM meta WHERE meta_key = 'progress_status'",
          "UPDATE analysis SET db_version = '".$self->o('refseq_version')."', db = 'RefSeq' WHERE logic_name = 'cdna_update'",
          "UPDATE analysis SET program_file= '".$self->o('exonerate_version')."' WHERE logic_name = 'cdna_update'", 
          "DELETE FROM analysis WHERE logic_name != 'cdna_update'",
          "DROP TABLE IF EXISTS analysis_description_bak",
          "DROP TABLE IF EXISTS attrib_type_bak",
          "DROP TABLE IF EXISTS external_db_bak",
          "DROP TABLE IF EXISTS misc_set_bak",
          "DROP TABLE IF EXISTS unmapped_reason_bak",
          "UPDATE gene SET biotype = 'cdna_update'",
          "UPDATE transcript SET biotype = 'cdna_update'",
          "UPDATE gene g, transcript t SET g.canonical_transcript_id = t.transcript_id WHERE g.gene_id = t.gene_id",
        ],
      },
      -flow_into => {
        1 => [ 'meta_levels' ],
      },
      -max_retry_count => 0,
    },
    {
      -logic_name => 'meta_levels',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('meta_level_script').
               ' --dbpattern '.$self->o('output_db','-dbname').
               ' --host '.$self->o('output_db','-host').
               ' --user '.$self->o('output_db','-user').
               ' --pass '.$self->o('output_db','-pass').
               ' --port '.$self->o('output_db','-port')
      },
      -flow_into => {
        1 => [ 'update_meta_coord' ],
      },
      -max_retry_count => 0,
    },
    {
      -logic_name => 'update_meta_coord',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl '.$self->o('meta_coord_script').
               ' --dbpattern '.$self->o('output_db','-dbname').
               ' --host '.$self->o('output_db','-host').
               ' --user '.$self->o('output_db','-user').
               ' --pass '.$self->o('output_db','-pass').
               ' --port '.$self->o('output_db','-port')
      },
      -max_retry_count => 0,
    },
  ];
}

sub pipeline_wide_parameters {
  my ($self) = @_;
  return {
    %{ $self->SUPER::pipeline_wide_parameters() },  # inherit other stuff from the base class
    complete_update => $self->o('complete_update'),
    wide_output_dir => $self->o('output_path'),
    wide_genome_file => $self->o('genome_file'),
    #wide_gss_file => $self->o('gss_file'),
    wide_refseq_file => $self->o('refseq_path') . $self->o('refseq_file'),
    wide_clone_script => $self->o('clone_db_script_path'),
    wide_optimize_script => $self->o('optimize_script'),
    wide_fastasplit => $self->o('fastasplit_random_path'),
    wide_exonerate => $self->o('exonerate_path'),
    wide_polyA => $self->o('polyA_script'),
    wide_delete_genes => $self->o('delete_genes_script'),
    wide_production_script => $self->o('populate_production_script'),
    wide_findN => $self->o('findN_script'),
    wide_analysis_desc => $self->o('analysis_desc_script'),
    wide_species => $self->o('species'),
  };
}

sub resource_classes {
  my $self = shift;

  my $default_queue = $self->default_options()->{'default_queue'}; 

  my $default_mem = $self->default_options()->{'default_mem'};
  my $exonerate_mem = $self->default_options()->{'exonerate_mem'};
  my $exonerate_retry_mem = $self->default_options()->{'exonerate_retry_mem'};
  my $exonerate_high_mem = $self->default_options()->{'exonerate_high_mem'};
  my $optimise_mem = $self->default_options()->{'optimise_mem'};
  my $download_mem = $self->default_options()->{'download_mem'};

  return {
    'default' => { LSF => '-q '.$default_queue.' -M'.$default_mem.' -R"select[mem>'.$default_mem.'] '.
                          'rusage[mem='.$default_mem.']"'},

    'exonerate' => { LSF => '-q '.$default_queue.' -M'.$exonerate_mem.' -R"select[mem>'.$exonerate_mem.'] '.
                            'rusage[mem='.$exonerate_mem.']"'},

    'exonerate_himem' => { LSF => '-q '.$default_queue.' -M'.$exonerate_high_mem.' -R"select[mem>'.$exonerate_high_mem.'] '.
                           'rusage[mem='.$exonerate_high_mem.']"'},

    'exonerate_2' => { LSF => '-q '.$default_queue.' -M'.$exonerate_retry_mem.' -R"select[mem>'.$exonerate_retry_mem.'] '.
                              'rusage[mem='.$exonerate_retry_mem.']"'},

    'optimise' => { LSF => '-q '.$default_queue.' -M'.$optimise_mem.' -R"select[mem>'.$optimise_mem.'] '.
                           'rusage[mem='.$optimise_mem.']"'},

    'download' => { LSF => '-q '.$default_queue.' -M'.$download_mem.' -R"select[mem>'.$download_mem.'] '.
                           'rusage[mem='.$download_mem.']"' },
  }
}


1;
