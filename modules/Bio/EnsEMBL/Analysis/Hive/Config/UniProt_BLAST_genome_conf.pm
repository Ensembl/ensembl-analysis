package Bio::EnsEMBL::Analysis::Hive::Config::UniProt_BLAST_genome_conf;

use warnings;
use strict;
use File::Spec::Functions;

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_analysis_settings);
use base ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

use Bio::EnsEMBL::ApiVersion qw/software_version/;

sub default_options {
  my ($self) = @_;
  return {
    # inherit other stuff from the base class
    %{ $self->SUPER::default_options() },

    ##########################
    # Required settings
    ##########################
    'blast_db_path'   => '', # this should be the location of a formatted/indexed/rm dump of the genome
    'output_path'     => '', # dir to download uniprot files to

    # read/write user and password for servers
    'user'            => '',
    'password'        => '',
    'user_r'          => '',

    # server connection info. The pipe and blast output dbs will be automatically created. The dna db should be an existing db (e.g. the human core)
    'pipe_db_name'    => '',
    'pipe_db_host'  => '',
    'pipe_db_port'    => '',
    'dna_db_name'     => '',
    'dna_db_host'   => '',
    'dna_db_port'     => '',
    'blast_db_name'   => '',
    'blast_db_host' => '',
    'blast_db_port'   => '',


    ##########################
    # Preset variables
    ##########################

    # UniProt settings
    'uniprot_set'          => 'havana_teleost_blast', # the UniProt set in the UniProtCladeDownloadStatic config to download
    'uniprot_blast_batch_size' => 10, # number of protein sequences per job
    'uniprot_table_name'   => 'uniprot_sequences', # table name to load the sequences into in the pipe db

    # Blast settings. The commandline is where to tweak blast parameters
    'blast_commandline' => ' -num_threads 3 -seg yes -soft_masking true -word_size 4 -threshold 20 -evalue 1e-2 -num_alignments 10000 ',
    'blast_type' => 'ncbi',
    'blast_exe_path' => catfile($self->o('binary_base'), 'tblastn'),

    load_optimise_script => catfile($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'genebuild', 'load_external_db_ids_and_optimize_af.pl'),

    'blast_db_user'   => $self->o('user'),
    'blast_db_pass'   => $self->o('password'),
    'blast_db_driver' => $self->o('hive_driver'),

    'killlist_db_name'   => 'gb_kill_list',
    'killlist_db_user'   => $self->o('user_r'),
    'killlist_db_pass'   => $self->o('password_r'),
    'killlist_db_driver' => $self->o('hive_driver'),

    'production_db_server' => 'mysql-ens-meta-prod-1',
    'production_db_port'   => '4483',

    'blast_db' => {
      -dbname => $self->o('blast_db_name'),
      -host   => $self->o('blast_db_host'),
      -port   => $self->o('blast_db_port'),
      -user   => $self->o('blast_db_user'),
      -pass   => $self->o('blast_db_pass'),
      -driver => $self->o('blast_db_driver'),
    },

    'killlist_db' => {
      -dbname => $self->o('killlist_db_name'),
      -host   => $self->o('killlist_db_host'),
      -port   => $self->o('killlist_db_port'),
      -user   => $self->o('killlist_db_user'),
      -pass   => $self->o('killlist_db_pass'),
      -driver => $self->o('killlist_db_driver'),
    },

    'production_db' => {
      -host   => $self->o('production_db_server'),
      -port   => $self->o('production_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -dbname => 'ensembl_production',
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
    ];
  }


## See diagram for pipeline structure
sub pipeline_analyses {
    my ($self) = @_;

    return [

     {
       -logic_name => 'create_uniprot_output_db',
       -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
       -parameters => {
                        source_db => $self->o('dna_db'),
                        target_db => $self->o('blast_db'),
                        create_type => 'clone',
                      },
       -rc_name    => 'default',
       -input_ids => [{}],
       -flow_into => {
                       1 => ['download_uniprot_files'],
                     },
     },


     {
       -logic_name => 'download_uniprot_files',
       -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadUniProtFiles',
       -parameters => {
                        multi_query_download => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::UniProtCladeDownloadStatic', $self->o('uniprot_set')),
                        output_path => $self->o('output_path'),
                      },
        -rc_name          => 'default',
        -flow_into => {
                        '2->A' => ['process_uniprot_files'],
                        'A->1' => ['generate_blast_jobs'],
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
        -logic_name => 'generate_blast_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         iid_type => 'sequence_accession',
                         batch_size => $self->o('uniprot_blast_batch_size'),
                         sequence_table_name => $self->o('uniprot_table_name'),
                       },
        -rc_name      => 'default',
        -flow_into => {
                        '2->A' => ['run_uniprot_blast'],
                        'A->1' => ['load_external_db_ids_and_optimise_af_tables'],
                      },
     },


      {
        -logic_name => 'run_uniprot_blast',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastPepToGenome',
        -parameters => {
                         sequence_table_name => $self->o('uniprot_table_name'),
                         sequence_type => 'peptide',
                         output_db => $self->o('blast_db'),
                         dna_db => $self->o('dna_db'),
                         logic_name => 'uniprot',
                         module     => 'HiveBlastPepToGenome',
                         blast_db_path => $self->o('blast_db_path'),
                         blast_exe_path => $self->o('blast_exe_path'),
                         commandline_params => $self->o('blast_commandline'),
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastUniProtToGenome', {BLAST_PARAMS => {type => $self->o('blast_type')}})},
                         timer => '2h',
                       },

        -flow_into => {
                        'MEMLIMIT' => ['resize_blast_jobs'],
                        'RUNLIMIT' => ['resize_blast_jobs'],
                      },
        -rc_name    => '4GB',
      },

     {
        -logic_name => 'resize_blast_jobs',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         iid_type => 'rechunk',
                         batch_size => 1,
                       },
        -rc_name      => 'default',
        -flow_into => {
                        2 => ['run_uniprot_blast_8GB'],
                      },
     },

      {
        -logic_name => 'run_uniprot_blast_8GB',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastPepToGenome',
        -parameters => {
                         sequence_table_name => $self->o('uniprot_table_name'),
                         sequence_type => 'peptide',
                         output_db => $self->o('blast_db'),
                         dna_db => $self->o('dna_db'),
                         logic_name => 'uniprot',
                         module     => 'HiveBlastPepToGenome',
                         blast_db_path => $self->o('blast_db_path'),
                         blast_exe_path => $self->o('blast_exe_path'),
                         commandline_params => $self->o('blast_commandline'),
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastUniProtToGenome', {BLAST_PARAMS => {type => $self->o('blast_type')}})},
                         timer => '2h',
                       },

        -flow_into => {
                        'MEMLIMIT' => ['failed_blast_jobs'],
                        'RUNLIMIT' => ['run_uniprot_blast_restrained'],
                      },
        -rc_name    => '8GB',
      },
      {
        -logic_name => 'run_uniprot_blast_restrained',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastPepToGenome',
        -parameters => {
                         sequence_table_name => $self->o('uniprot_table_name'),
                         sequence_type => 'peptide',
                         output_db => $self->o('blast_db'),
                         dna_db => $self->o('dna_db'),
                         logic_name => 'uniprot_low',
                         module     => 'HiveBlastPepToGenome',
                         blast_db_path => $self->o('blast_db_path'),
                         blast_exe_path => $self->o('blast_exe_path'),
                         commandline_params => $self->o('blast_commandline'),
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastUniProtToGenome', {BLAST_PARAMS => {type => $self->o('blast_type')}})},
                         timer => '4h',
                       },

        -flow_into => {
                        'MEMLIMIT' => ['failed_blast_jobs'],
                        'RUNLIMIT' => ['failed_blast_jobs'],
                      },
        -rc_name    => '8GB',
      },

      {
        -logic_name => 'failed_blast_jobs',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -parameters => {
                       },
        -rc_name          => 'default',
        -can_be_empty  => 1,
        -failed_job_tolerance => 100,
      },
      {
        -logic_name => 'load_external_db_ids_and_optimise_af_tables',
        -module => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters => {
          cmd => 'perl '.$self->o('load_optimise_script').
            ' -output_path '.catdir($self->o('output_path'), 'optimise').
            ' -uniprot_filename '.$self->o('protein_entry_loc').
            ' -dbuser '.$self->o('user').
            ' -dbpass '.$self->o('password').
            ' -dbport '.$self->o('blast_db','-port').
            ' -dbhost '.$self->o('blast_db','-host').
            ' -dbname '.$self->o('blast_db','-dbname').
            ' -prod_dbuser '.$self->o('user_r').
            ' -prod_dbhost '.$self->o('production_db','-host').
            ' -prod_dbname '.$self->o('production_db','-dbname').
            ' -prod_dbport '.$self->o('production_db','-port').
            ' -core'
        },
        -max_retry_count => 0,
        -rc_name => '8GB',
      },

    ];
  }


sub resource_classes {
    my $self = shift;
    return {
      'default' => { LSF => $self->lsf_resource_builder('production-rh74', 900, [$self->default_options->{'pipe_db_host'}, $self->default_options->{'dna_db_host'}], [$self->default_options->{'num_tokens'}])},
      '4GB' => { LSF => $self->lsf_resource_builder('production-rh74', 4000, [$self->default_options->{'pipe_db_host'}, $self->default_options->{'dna_db_host'},$self->default_options->{'blast_db_host'}], undef, 3)},
      '8GB' => { LSF => $self->lsf_resource_builder('production-rh74', 8000, [$self->default_options->{'pipe_db_host'}, $self->default_options->{'dna_db_host'},$self->default_options->{'blast_db_host'}], undef, 3)},
    }
  }

1;
