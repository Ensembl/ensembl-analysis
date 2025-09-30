=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#Copyright [2016-2024] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 DESCRIPTION

This is the config file for all blast analysis. You should use it in your Hive configuration file to
specify the parameters of an analysis. You can either choose an existing config or you can create
a new one based on the default hash. 

=head1 SYNOPSIS

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_analysis_settings);
use parent ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

sub pipeline_analyses {
    my ($self) = @_;

    return [
      {
        -logic_name => 'run_uniprot_blast',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlastGenscanPep',
        -parameters => {
                         blast_db_path => $self->o('uniprot_blast_db_path'),
                         blast_exe_path => $self->o('uniprot_blast_exe_path'),
                         commandline_params => '-cpus 3 -hitdist 40',
                         repeat_masking_logic_names => ['repeatmasker_'.$self->o('repeatmasker_library')],
                         prediction_transcript_logic_names => ['genscan'],
                         iid_type => 'feature_id',
                         logic_name => 'uniprot',
                         module => 'HiveBlastGenscanPep',
                         %{get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic','BlastGenscanPep')},
                      },
        -flow_into => {
                        -1 => ['run_uniprot_blast_himem'],
                        -2 => ['run_uniprot_blast_long'],
                      },
        -rc_name    => 'blast',
      },
  ];
}

=head1 METHODS

  _master_config_settings: contains all possible parameters

=cut

package Bio::EnsEMBL::Analysis::Hive::Config::BlastStatic;

use strict;
use warnings;

use parent ('Bio::EnsEMBL::Analysis::Hive::Config::BaseStatic');

sub _master_config {
  my ($self, $key) = @_;

  my %config = (
      default => {
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
            -type => 'ncbi',
          }
      },

      BlastGenscanPep => {
        PARSER_PARAMS => {
                           -regex => '^\s*(\w+\W\d+)',
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

      BlastGenscanPep_non_vert => {
        PARSER_PARAMS => {
                           -regex => '^\s*([^\s]+)',
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

      BlastGenscanVertRNA => {
        BLAST_PARSER => 'Bio::EnsEMBL::Analysis::Tools::FilterBPlite',
        PARSER_PARAMS => {
                           -regex => '^\s*(\w+\W\d+)',
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

      BlastGenscanUnigene => {
        BLAST_PARSER => 'Bio::EnsEMBL::Analysis::Tools::FilterBPlite',
        PARSER_PARAMS => {
                           -regex => '^\s*(\w+\.\d+)',
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
        BlastRFam => {
          BLAST_PARSER => 'Bio::EnsEMBL::Analysis::Tools::FilterBPlite',
          PARSER_PARAMS => {
                             -regex => '(\w+)\.\w+',
                             -query_type => 'dna',
                             -database_type => 'dna',
                             -threshold => 0.01,
                           },
        },
        BlastmiRBase => {
          BLAST_PARSER => 'Bio::EnsEMBL::Analysis::Tools::BPliteWrapper',
          PARSER_PARAMS => {
                             -regex => '\S+\s+(\w+)',
                             -query_type => 'dna',
                             -database_type => 'dna',
                           },
          },

      BlastUniProtToGenome => {
        BLAST_PARSER => 'Bio::EnsEMBL::Analysis::Tools::FilterBlastGenome',
        PARSER_PARAMS => {
                           -regex => '^(\S+:\S+:\S+:\d+:\d+:1)',
                           -query_type => 'pep',
                           -database_type => 'dna',
                           -threshold_type => 'PVALUE',
                           -threshold => 0.00001,
                           -filter => 0,
                         },
        BLAST_FILTER => 'Bio::EnsEMBL::Analysis::Tools::FeatureFilterOnGenome',
        FILTER_PARAMS => {
                           -prune => 1,
                           -coverage => 3,
                         },
      },

  );
  return $config{$key};
}

1;
