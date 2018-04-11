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

'farm_user_name'        => '', # for output_db prefix
'user_r'                => '',
'user_w'                => '',
'password'              => '',

'base_repeat_dir'         => '',
'repeatmodeler_path'      => '',

'pipeline_name'         => '',
'pipe_db_server'        => '',
'port'                  => '',

'assembly_registry_dbname'    => 'do1_stable_id_space_assembly_registry',
'assembly_registry_db_server' => $ENV{GBS5},
'assembly_registry_db_port'   => $ENV{GBP5},

'min_contig_n50'          => 20000,
'min_scaffold_n50'        => 0,
'min_total_length'        => 500000000,
'sleep_length_hours'      => 1,

'base_output_path'        => catfile($self->o('base_repeat_dir'),'pipeline_output_dir'),
'min_consensi_files'      => 5,

'repeatmodeler_run_count' => 10,
'num_cores'               => 10,

'ncbi_base_ftp'           => 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all',

'base_output_path'        => catfile($self->o('base_repeat_dir'),'pipeline_output_dir'),
'binary_base'           => '/nfs/software/ensembl/RHEL7/linuxbrew/bin',
'default_mem'           => '2900',
'default_himem'         => '20000',

'user'                  => $self->o('user_w'),
'dna_db_name'           => '',

'assembly_registry_db' => {
              -dbname => $self->o('assembly_registry_dbname'),
              -host   => $self->o('assembly_registry_db_server'),
              -port   => $self->o('assembly_registry_db_port'),
              -user   => $self->o('user_r'),
              -driver => $self->o('hive_driver'),
            },

 } # end return
} # end default_options


sub pipeline_create_commands {
    my ($self) = @_;
    return [
    # inheriting database and hive tables' creation
	    @{$self->SUPER::pipeline_create_commands},
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
                      cmd => 'sleep 30; cd #path_to_genomic_fasta# && '.$self->o('repeatmodeler_path').'/BuildDatabase -name repeatmodeler_db -engine ncbi #path_to_genomic_fasta#/genomic.fna',
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
                       cmd => 'sleep 30; cd #repeatmodeler_run_dir# && '.$self->o('repeatmodeler_path').'/RepeatModeler -engine ncbi -pa '.$self->o('num_cores').' -database #path_to_genomic_fasta#/repeatmodeler_db',
                     },
      -rc_name    => 'default_himem',
      -max_retry_count => 1,
      -flow_into  => {
                       'ANYFAILURE' => ['failed_repeatmodeler_jobs'],
                     },
      -analysis_capacity => 250,
    },


    {
      -logic_name => 'failed_repeatmodeler_jobs',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -parameters => {
                     },
      -rc_name          => 'default',
      -can_be_empty  => 1,
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
      -flow_into  => {
                       '1' => ['delete_output_dir'],
                     },

    },


    {
      -logic_name => 'delete_output_dir',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
                       cmd => 'perl '.catfile($self->o('enscode_root_dir'),'ensembl-analysis/scripts','delete_big_dir.pl').' -dir #path_to_genomic_fasta#',
                     },
      -rc_name    => 'default',
    },

  ]
} # end pipeline analyses

sub resource_classes {
  my $self = shift;

  return {
    'default'       => { 'LSF' => $self->lsf_resource_builder('production-rh7',$self->default_options->{'default_mem'}) },
    'default_himem' => { 'LSF' => $self->lsf_resource_builder('production-rh7',$self->default_options->{'default_himem'},undef,undef,$self->o('num_cores')) },
  }
} # end resource_classes

1;
