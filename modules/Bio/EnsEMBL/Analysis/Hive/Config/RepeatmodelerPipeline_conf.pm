package RepeatmodelerPipeline_conf;

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
     user => $ENV{'USER_R'},
     password => $ENV{'GBPASS'},
     user_w => $ENV{'GBUSER'},
     'pipeline_db' => {
	-dbname => 'repeatmodeler_test',
        -host   => $ENV{GBS1},
        -port   => $ENV{GBP1},
        -user   => $self->o('user_w'),
         -pass   => $self->o('password'),
	-driver => 'mysql',
     },

# Handy SQL for looking at failed or passed on repeatmodeler jobs:
# select count(*),substring(data, 12,15) as failed from job join analysis_data on analysis_data_id=replace(input_id,'_extended_data_id ','') 
# where analysis_id=(select analysis_id from analysis_base where logic_name='run_repeatmodeler') and status in ('FAILED','PASSED_ON') group by failed;


# Things to set
'farm_user_name'              => '' || $ENV{EHIVE_USER} || $ENV{USER},
'user_r'                      => 'ensro',
'user_w'                      => 'ensadmin',
'password'                    => 'ensembl',

'base_repeat_dir'             => '/hps/nobackup/flicek/ensembl/genebuild/custom_repeat_libraries/repeatmodeler/',
'repeatmodeler_path'          => '/hps/software/users/ensembl/genebuild/RepeatModeler-2.0.3/',
'public_ftp'                  => '/nfs/ftp/public/databases/ensembl/repeats/unfiltered_repeatmodeler/species/',
#'pipeline_name'               => 'repeatmodeler_test',
#'pipe_db_server'              => $ENV{GBS1},
#'port'                        => $ENV{GBP1},

'assembly_registry_dbname'    => 'gb_assembly_registry',
'assembly_registry_db_server' => $ENV{GBS1},
'assembly_registry_db_port'   => $ENV{GBP1},

# Variables (with suggested values)
'min_contig_n50'              => 100000, # min contig n50 to process
'assembly_group'              => 'dtol', # assembly group to process
'min_scaffold_n50'            => 0, # min scaffold n50 to process
'min_total_length'            => 50000000, # min amount of total sequence
'sleep_length_hours'          => 1, # How long to sleep between checking for new assemblies
'min_consensi_files'          => 1, # Min number of repeatmodeler runs for an assembly before proceeding
'repeatmodeler_run_count'     => 1, # Number of runs of repeatmodeler per assembly
'num_cores'                   => 10, # Number of cores for BLAST phase of repeatmodeler

# Mostly constant stuff
#'pipeline_name'               => 'vertebrate_repeatmodeler',
'ncbi_base_ftp'               => 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all',
'base_output_path'            => catfile($self->o('base_repeat_dir'),'pipeline_output_dir'),
'default_mem'                 => '2900',
'default_himem'               => '20000',
'default_himem_50'            => '50000',
'default_himem_75'            => '75000',
'datamover_5'                 => '5000',
'user'                        => $self->o('user_w'),
'dna_db_name'                 => '', # Leave blank, just needs to be present
'assembly_registry_db'        => {
                                   -dbname => $self->o('assembly_registry_dbname'),
                                   -host   => $self->o('assembly_registry_db_server'),
                                   -port   => $self->o('assembly_registry_db_port'),
                                   -user   => $self->o('user_w'),
				   -pass   => $self->o('password'),
                                   -driver => $self->o('hive_driver'),
                                 },




 } # end return
} # end default_options


sub pipeline_create_commands {
  my ($self) = @_;
  return [
    @{$self->SUPER::pipeline_create_commands},
    # inheriting database and hive tables' creation
    $self->db_cmd('CREATE TABLE run_records ('.
                  'accession varchar(50) NOT NULL,'.
                  'species_name varchar(100) NOT NULL,'.
                  'run_count int NOT NULL,'.
                  'status varchar(50) NOT NULL,'.
                  'PRIMARY KEY (accession))'),
  ];
} # end pipeline_create_commands



sub pipeline_analyses {
  my ($self) = @_;
  return [

    {
      -logic_name => 'seed_assembly_jobs',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSeedRepeatPipeline',
      -parameters => {
                       assembly_registry_db => $self->o('assembly_registry_db'),
                       base_repeat_dir      => $self->o('base_repeat_dir'),
                       min_contig_n50       => $self->o('min_contig_n50'),
                       assembly_group       => $self->o('assembly_group'),
                       min_scaffold_n50     => $self->o('min_scaffold_n50'),
                       min_total_length     => $self->o('min_total_length'),
                       sleep_length_hours   => $self->o('sleep_length_hours'),
                       max_version_only     => 1,
                     },
      -flow_into  => {
                       '2' => ['download_genome_fasta'],
                       '3' => ['seed_assembly_jobs'],
                     },
      -input_ids  => [{'iid' => 1}],
    },


    {
      -logic_name => 'download_genome_fasta',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDownloadGenomeFlatfiles',
      -parameters => {
                       assembly_registry_db => $self->o('assembly_registry_db'),
                       ncbi_base_ftp        => $self->o('ncbi_base_ftp'),
                       base_output_path     => $self->o('base_output_path'),
                     },
      -flow_into  => {
                        '1->A' => ['build_repeatmodeler_db'],
                        'A->1' => ['store_generated_libraries'],
                     },
      -analysis_capacity => 20,
    },


    {
      -logic_name => 'build_repeatmodeler_db',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
                      cmd => 'sleep 120; cd #path_to_genomic_fasta# && '.$self->o('repeatmodeler_path').'/BuildDatabase -name repeatmodeler_db -engine ncbi #path_to_genomic_fasta#/genomic.fna',
                     },
      -rc_name    => 'default_himem',
      -max_retry_count => 1,
      -flow_into  => {
                       '1' => ['generate_repeatmodeler_jobs'],
                     },
    },


    {
      -logic_name => 'generate_repeatmodeler_jobs',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenerateRepeatmodelerJobs',
      -parameters => {
                       repeatmodeler_run_count => $self->o('repeatmodeler_run_count'),
                     },

      -flow_into  => {
                       '2' => ['run_repeatmodeler'],
                     },
    },


    {
      -logic_name => 'run_repeatmodeler',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
                       cmd => 'sleep 120; cd #repeatmodeler_run_dir# && '.$self->o('repeatmodeler_path').'/RepeatModeler -engine ncbi -numAddlRounds 1 -pa '.$self->o('num_cores').' -database #path_to_genomic_fasta#/repeatmodeler_db',
                     },
      -rc_name    => 'default_himem',
      -max_retry_count => 1,
      -flow_into  => {
                       '-1' => ['failed_repeatmodeler_jobs_50GB'],
                     },
      -analysis_capacity => 250,
    },
    
    {
      -logic_name => 'failed_repeatmodeler_jobs_50GB',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
                       cmd => 'sleep 120; cd #repeatmodeler_run_dir# && '.$self->o('repeatmodeler_path').'/RepeatModeler -engine ncbi -numAddlRounds 1 -pa '.$self->o('num_cores').' -database #path_to_genomic_fasta#/repeatmodeler_db',
                     },
      -rc_name    => 'default_himem_50',
      -max_retry_count => 1,
      -flow_into  => {
                       '-1' => ['failed_repeatmodeler_jobs_75GB'],
                     },
      -analysis_capacity => 250,
    },

   {
      -logic_name => 'failed_repeatmodeler_jobs_75GB',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
                       cmd => 'sleep 120; cd #repeatmodeler_run_dir# && '.$self->o('repeatmodeler_path').'/RepeatModeler -engine ncbi -numAddlRounds 1 -pa '.$self->o('num_cores').' -database #path_to_genomic_fasta#/repeatmodeler_db',
                     },
      -rc_name    => 'default_himem_75',
      -max_retry_count => 1,
      -analysis_capacity => 250,
    },


    {
      -logic_name => 'store_generated_libraries',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRepeatmodelerResults',
      -parameters => {
                       path_to_assembly_libs => catfile($self->o('base_repeat_dir'),'GCA'),
                       path_to_species_libs => catfile($self->o('base_repeat_dir'),'species'),
                       min_consensi_files    => $self->o('min_consensi_files'),
                       assembly_registry_db  => $self->o('assembly_registry_db'),
                     },

      -rc_name    => 'default',
      -flow_into => {
	      1 => ['update_run_count_table'],
	      2 => {'update_ftp' => { 'species_name' => '#species_name#', 'species_lib' => '#species_lib#'}
		   },
	   }
    },


    {
      -logic_name => 'update_ftp',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
	               cmd => "rsync -ahvW #species_lib# " . $self->o('public_ftp') . '#species_name#/', 
                     },
       -rc_name    => 'move_library_ftp',
    },
    {
      -logic_name => 'update_run_count_table',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
                       db_conn => $self->o('pipeline_db'),
                         sql => ['UPDATE run_records set run_count=#run_count# WHERE accession="#iid#"',
                                 'UPDATE run_records set status="complete" WHERE accession="#iid#"'],
                     },
      -max_retry_count => 0,
      -rc_name => 'default',
      -flow_into  => {
                       '1' => ['delete_output_dir'],
                     },

    },


    {
      -logic_name => 'delete_output_dir',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
                       cmd => 'perl '.catfile($ENV{'ENSEMBL_ROOT_DIR'},'ensembl-analysis/scripts','delete_big_dir.pl').' -dir #path_to_genomic_fasta#',
                     },
      -rc_name    => 'default',
    },

  ]
} # end pipeline analyses

sub resource_classes {
  my $self = shift;

  return {
    'default'       => { 'LSF' => $self->lsf_resource_builder('production',$self->default_options->{'default_mem'}) },
    'default_himem' => { 'LSF' => $self->lsf_resource_builder('production',$self->default_options->{'default_himem'},undef,undef,$self->o('num_cores')) },
    'default_himem_50' => { 'LSF' => $self->lsf_resource_builder('production',$self->default_options->{'default_himem_50'},undef,undef,$self->o('num_cores')) },
    'default_himem_75' => { 'LSF' => $self->lsf_resource_builder('production',$self->default_options->{'default_himem_75'},undef,undef,$self->o('num_cores')) },
    'move_library_ftp' => { 'LSF' => $self->lsf_resource_builder('datamover', $self->default_options->{'datamover_5'},undef,undef,$self->o('num_cores')) },
  }
} # end resource_classes

1;
