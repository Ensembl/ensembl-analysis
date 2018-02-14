package HiveCesar_conf;

use warnings;
use strict;
use feature 'say';

use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

use Bio::EnsEMBL::ApiVersion qw/software_version/;

sub default_options {
    my ($self) = @_;
    return {
        # inherit other stuff from the base class
	    %{ $self->SUPER::default_options() },

'enscode_root_dir' => '/path/to/enscode/',
'cesar_path' => '/path/to/CESAR2.0/',
'user' => 'READ_USER',
'user_r' => 'READ_USER',
'user_w' => 'WRITE_USER',
'password' => 'WRITE_PASS',
'driver' => 'mysql',
'method_link_type' => 'LASTZ_NET', # no need to modify this
'clone_db_script_path' => $self->o('enscode_root_dir').'/ensembl-analysis/scripts/clone_database.ksh', # no need to modify this

# database details for the eHive pipe database
'server1' => '',
'port1' => '',
'pipeline_dbname' => 'USER_SPECIES_NAME2_cesar_pipe_RELEASE_ASSEMBLY', # this db will be created

# database details for the other databases
'server2' => '',
'port2' => '',
'compara_dbname' => 'USER_SPECIES_NAME2_lastz_RELEASE_ASSEMBLY', # compara db containing the lastz alignments
'source_dna_dbname' => 'USER_SPECIES_NAME1_core_RELEASE_ASSEMBLY', # core db containing the dna for the gene set to be projected
'source_transcript_dbname' => 'USER_SPECIES_NAME1_core_RELEASE_ASSEMBLY', # core db containing the gene set to be projected
'projection_dna_dbname' => 'USER_SPECIES_NAME2_core_RELEASE_ASSEMBLY', # core db containing the dna for the projected gene set
'projection_dbname' => 'USER_SPECIES_NAME2_cesar_RELEASE_ASSEMBLY', # core db which will be created to store the projected gene set

'output_path' => '/path/to/cesar_output',

'pipeline_name' => 'cesar',

        'pipeline_db' => {
          # connection parameters
            -dbname => $self->o('pipeline_dbname'),
            -host   => $self->o('server1'),
            -port   => $self->o('port1'),
            -user   => $self->o('user_w'),
            -pass   => $self->o('password'),
            -driver => $self->o('driver'),
        },

        'source_dna_db' => {
                            -dbname => $self->o('source_dna_dbname'),
                            -host   => $self->o('server2'),
                            -port   => $self->o('port2'),
                            -user   => $self->o('user'),
                          },

        'source_transcript_db' => {
                            -dbname => $self->o('source_transcript_dbname'),
                            -host   => $self->o('server2'),
                            -port   => $self->o('port2'),
                            -user   => $self->o('user'),
                          },

        'projection_dna_db' => {
                            -dbname => $self->o('projection_dna_dbname'),
                            -host   => $self->o('server2'),
                            -port   => $self->o('port2'),
                            -user   => $self->o('user'),
                          },

        'projection_db' => {
                            -dbname => $self->o('projection_dbname'),
                            -host   => $self->o('server2'),
                            -port   => $self->o('port2'),
                            -user   => $self->o('user_w'),
                            -pass   => $self->o('password'),
                          },

         'compara_db' => {
                            -dbname => $self->o('compara_dbname'),
                            -host   => $self->o('server2'),
                            -port   => $self->o('port2'),
                            -user   => $self->o('user'),
                          },

    };
  }

sub pipeline_create_commands {
    my ($self) = @_;
    return [
    # inheriting database and hive tables' creation
	    @{$self->SUPER::pipeline_create_commands},
    ];
  }


## See diagram for pipeline structure
sub pipeline_analyses {
    my ($self) = @_;

    return [

	    {
        -input_ids  => [ {} ],
        -logic_name => 'create_output_db',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
        -parameters => {
                         source_db => $self->o('projection_dna_db'),
                         target_db => $self->o('projection_db'),
                         create_type => 'clone',
                         script_path => $self->o('clone_db_script_path'),
                         user_r => $self->o('user_r'),
                         user_w => $self->o('user_w'),
                         pass_w => $self->o('password'),
                       },
        -rc_name    => 'default',
        -flow_into => { 1 => ['create_projection_input_ids'] },
      },

      {
        -logic_name => 'create_projection_input_ids',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
        -parameters => {
                         target_db => $self->o('source_transcript_db'),
                         iid_type => 'feature_id',
                         batch_size => 20,
                         feature_type => 'gene',
                         feature_restriction => 'projection',
                         allowed_biotypes => {'protein_coding' => 1},
                       },

        -flow_into => {
                        2 => ['cesar'],
                      },
         -rc_name    => 'default',
      },

      {
        -logic_name => 'cesar',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCesar',
        -parameters => {
                         'output_path' => $self->o('output_path')."/project_exons/",
                         'source_dna_db' => $self->o('source_dna_db'),
                         'target_dna_db' => $self->o('projection_dna_db'),
                         'source_db' => $self->o('source_transcript_db'),
                         'target_db' => $self->o('projection_db'),
                         'compara_db' => $self->o('compara_db'),
                         'method_link_type' => $self->o('method_link_type'),
                         'exon_region_padding' => 50,
                         'cesar_path' => $self->o('cesar_path'),
                         TRANSCRIPT_FILTER => {
                           OBJECT     => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
                           PARAMETERS => {
                             -coverage => 50,
                             -percent_id => 50,
                           },
                         },
                       },
        -rc_name    => 'default',
        -analysis_capacity => 300,
        -flow_into => {
                        -1 => ['cesar_himem'],
                      },
      },

      {
        -logic_name => 'cesar_himem',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCesar',
        -parameters => {
                         'output_path' => $self->o('output_path')."/project_exons/",
                         'source_dna_db' => $self->o('source_dna_db'),
                         'target_dna_db' => $self->o('projection_dna_db'),
                         'source_db' => $self->o('source_transcript_db'),
                         'target_db' => $self->o('projection_db'),
                         'compara_db' => $self->o('compara_db'),
                         'method_link_type' => $self->o('method_link_type'),
                         'exon_region_padding' => 50,
                         'cesar_path' => $self->o('cesar_path'),
                         'cesar_mem' => '20', # mem in GB to be used by cesar (parameter --max-memory)
                         TRANSCRIPT_FILTER => {
                           OBJECT     => 'Bio::EnsEMBL::Analysis::Tools::ExonerateTranscriptFilter',
                           PARAMETERS => {
                             -coverage => 50,
                             -percent_id => 50,
                           },
                         },
                       },
        -rc_name    => 'default_20GB',
        -analysis_capacity => 300,
      },
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
      'default' => { LSF => '-M1900 -R"select[mem>1900] rusage[mem=1900]"' },
      'default_20GB' => { LSF => '-M20000 -R"select[mem>20000] rusage[mem=20000]"' },
    }
  }

1;
