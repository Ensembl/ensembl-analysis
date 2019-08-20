=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

# reference_db - basically get a copy of the core
# Also get the gene models from relevant species rnaseq_db
# parse augustus output to gene model format and put in target_db (basically core plus new output)
# add in option to read in soft repeat masked genome from flat file as well as from db


package augustus_test_conf;

use strict;
use warnings;
use File::Spec::Functions;
use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_analysis_settings);
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');


sub default_options {
  my ($self) = @_;

  return {
    # inherit other stuff from the base class
    %{ $self->SUPER::default_options() },

    'enscode_root_dir'        => '/homes/jma/enscode', # path to the code checkout
    'clone_db_script_path'    => $self->o('enscode_root_dir').'/ensembl-analysis/scripts/clone_database.ksh',
    'output_path'             => '/hps/nobackup2/production/ensembl/jma',
    'softmasked_genome_file'  => '/hps/nobackup2/production/ensembl/genebuild/blastdb/human/toplevel.with_nonref_and_GRCh38_p7.no_duplicate.softmasked_dusted.fa',

    'reference_db' => {
      -driver   => 'mysql',
      -host     => 'ensembldb.ensembl.org',
      -user     => 'anonymous',
      -dbname   => 'homo_sapiens_98_38'
    },

#    'reference_db' => {
#      -driver   => 'mysql',
#      -host     => 'mysql-ens-genebuild-prod-1',
#      -user     => 'anonymous',
#      -dbname   => 'homo_sapiens_98_38'
#    },

    'pipeline_db' => {
      -driver   => 'mysql',
      -host     => 'mysql-ens-genebuild-prod-1',
      -port     => 4527,
      -user     => 'ensadmin',
      -password => 'ensembl',
      -dbname   => 'augustus_testing_pipe'
    },

    'target_db' => {
      -driver   => 'mysql',
      -host     => 'mysql-ens-genebuild-prod-2',
      -port     => 4528,
      -user     => 'ensadmin',
      -password => 'ensembl',
      -dbname   => 'augustus_testing_output'
    },

  };
}


sub pipeline_analyses {
  my ($self) = @_;
	return [
    {

      -logic_name => 'create_augustus_slices',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
                         target_db => $self->o('dna_db'),
                         iid_type => 'slices',
                      },
      -rc_name    => 'default',
      -flow_into => {
                        '2' => ['run_augustus'],
                    },
    },


    {
  	 -logic_name => 'run_augustus',
  	 -module => 'Bio::EnsEMBL::Analysis::RunnableDB::Finished::Augustus'
  	 -parameters => {
  		 target_db      => $self->o('reference_db'),
       logic_name     => 'augustus',
       module         => 'HiveAugustus',
      }
  	},

  ];
}

sub pipeline_wide_parameters {
  my ($self) = @_;
  return {
    %{ $self->SUPER::pipeline_wide_parameters() },  # inherit other stuff from the base class
    softmasked_genome_file  => $self->o('genome_file'),
    augustus_path           => $self->o('augustus_path'),
    SOFT_MASKED_REPEATS     => $self->o('repeat_masking_logic_names'),
  };
}

sub resource_classes {
  my $self = shift;

  return {
    'default' => { LSF => $self->lsf_resource_builder('production-rh7', 2000, [$self->default_options->{'pipe_db_server'}, $self->default_options->{'dna_db_server'}])},  
  }
}
1;
