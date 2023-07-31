package Bio::EnsEMBL::Analysis::Hive::Config::RedRepeatPipeline_conf;

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
     user_r => $self->o('user_r'),
     password => $self->o('password'),
     user_w => $self->o('user_w'),
     'pipeline_db' => {
	-dbname => 'red_repeat_masking',
        -host   => $self->o('host'),
        -port   => $self->o('port'),
        -user   => $self->o('user_w'),
        -pass   => $self->o('password'),
	-driver => 'mysql',
     },

# Handy SQL for looking at failed or passed on red jobs:
# select count(*),substring(data, 12,15) as failed from job join analysis_data on analysis_data_id=replace(input_id,'_extended_data_id ','') 
# where analysis_id=(select analysis_id from analysis_base where logic_name='run_red') and status in ('FAILED','PASSED_ON') group by failed;


# Things to set
'farm_user_name'              => '' || $ENV{EHIVE_USER} || $ENV{USER},

'base_repeat_dir'             => '/hps/nobackup/flicek/ensembl/genebuild/custom_repeat_libraries/red/',
'red_path'                    => '/hps/software/users/ensembl/ensw/C8-MAR21-sandybridge/linuxbrew/bin/Red',

'assembly_registry_dbname'    => $self->o('reg_db'),,
'assembly_registry_db_server' => $self->o('host'),
'assembly_registry_db_port'   => $self->o('port'),

# Variables (with suggested values)
'min_contig_n50'              => 100000, # min contig n50 to process
'assembly_group'              => 'dtol', # assembly group to process
'clade_group'              => "'amphibians','aves','humans','mammalia','marsupials','primates','reptiles','rodentia','sharks','teleostei','vertebrates'", # set of clades to process
'min_scaffold_n50'            => 0, # min scaffold n50 to process
'min_total_length'            => 50000000, # min amount of total sequence
'sleep_length_hours'          => 24, # How long to sleep between checking for new assemblies
'run_count'                   => 1, # Number of runs of red per assembly
'num_cores'                   => 10, # Number of cores for BLAST phase of red

# Mostly constant stuff
'ncbi_base_ftp'               => 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all',
'base_output_path'            => catfile($self->o('base_repeat_dir'),'pipeline_output_dir'),
'red_library_path'            => catfile($self->o('base_repeat_dir'),'libraries'),
'default_mem'                 => '2900',
'default_himem'               => '20000',
'default_himem_50'            => '50000',
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
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSeedRedPipeline',
      -parameters => {
                       assembly_registry_db => $self->o('assembly_registry_db'),
                       base_repeat_dir      => $self->o('base_repeat_dir'),
                       min_contig_n50       => $self->o('min_contig_n50'),
                       assembly_group       => $self->o('assembly_group'),
		       clade_group          => $self->o('clade_group'),
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
	               '1' => ['generate_red_repeat_jobs'],
                     },
      -analysis_capacity => 20,
    },


    {
      -logic_name => 'generate_red_repeat_jobs',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenerateRepeatmodelerJobs',
      -parameters => {
                       run_count => $self->o('run_count'),
                     },

      -flow_into  => {
                       '2->A' => ['run_red'],
		       'A->1' => ['store_generated_libraries'],
                     },
    },


    {
      -logic_name => 'run_red',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
                       cmd => 'sleep 120; mkdir '.$self->o('red_library_path').'/#species#/#iid# -p; ' .'[ -f #path_to_genomic_fasta#/genomic.fna ] && mv #path_to_genomic_fasta#/genomic.fna #path_to_genomic_fasta#/#iid#_red.fa && '.$self->o('red_path').' -frm 2  -gnm #path_to_genomic_fasta# -msk '.$self->o('red_library_path').'/#species#/#iid# -rpt '.$self->o('red_library_path').'/#species#/#iid#',
                     },
      -rc_name    => 'default_himem',
      -max_retry_count => 1,
      -flow_into  => {
                       '-1' => ['failed_red_jobs_50GB'],
                     },
      -analysis_capacity => 250,
    },
    
    {
      -logic_name => 'failed_red_jobs_50GB',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
	      cmd => 'sleep 120; mkdir '.$self->o('red_library_path').'/#species#/#iid# -p; ' .'[ -f #path_to_genomic_fasta#/genomic.fna ] && mv #path_to_genomic_fasta#/genomic.fna #path_to_genomic_fasta#/#iid#_red.fa && '.$self->o('red_path').' -frm 2  -gnm #path_to_genomic_fasta# -msk '.$self->o('red_library_path').'/#species#/#iid# -rpt '.$self->o('red_library_path').'/#species#/#iid#',
	                           },
      -rc_name    => 'default_himem_50',
      -max_retry_count => 1,
      -analysis_capacity => 250,
    },


    {
      -logic_name => 'store_generated_libraries',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRedResults',
      -parameters => {
                       assembly_registry_db  => $self->o('assembly_registry_db'),
		       path_to_repeat_libs   => $self->o('red_library_path').'/#species#/#iid#',
                     },

      -rc_name    => 'default',
      -flow_into => {
	      '1' => ['update_run_count_table'],
       },
    },


    {
      -logic_name => 'update_run_count_table',
      -module => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
                       db_conn => $self->o('pipeline_db'),
                       sql => ['UPDATE run_records set run_count=#run_count# WHERE accession="#iid#"',
                                 'UPDATE run_records set status="complete" WHERE accession="#iid#"'],
		       run_count             => $self->o('run_count'),
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
  }
} # end resource_classes

1;
