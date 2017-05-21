=head1 LICENSE

Copyright [2017] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

=head1 DESCRIPTION

This is the config file for the RecoverRNASeqTissuesModels analysis. You should use it in your Hive configuration
file to specify the parameters of an analysis. You can either choose an existing config or you can create
a new one based on the default hash.


=head1 SYNOPSIS

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_analysis_settings);
use parent ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

sub pipeline_analyses {
    my ($self) = @_;

    return [
      {
        -logic_name => 'run_recover_tissue_models',
        -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::RecoverRNASeqTissuesModels',
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


package Bio::EnsEMBL::Analysis::Hive::Config::RecoverRNASeqTissuesModelsStatic;

use strict;
use warnings;

use parent ('Bio::EnsEMBL::Analysis::Hive::Config::BaseStatic');

sub _master_config {
  my ($self, $key) = @_;

  my %config = (
      default => {
        # set the verbosity level
        VERBOSE => 0,

        # Values between START and END shouldn't be changed
        # unless you know what you are doing
        ## START ##
        CODING_ONLY      => 0,
        IGNORE_STRAND    => 0,
        CHECK_FRAMESHIFT => 0,
        CHECK_STOP       => 0,

        CHECK_TWOWAY => 1,
        CHECK_SINGLE => 1,
        CHECK_ORPHAN => 1,
        ## END ##

        # If the value of the biotype is '', it's not changed
        # Otherwise it is concatenated to the original one
        # biotype of good genes
        DELETE_BIOTYPE => 'deleted',
        SINGLE_BIOTYPE => 'single',
        ORPHAN_BIOTYPE => 'orphan',

        # If GOOD_BIOTYPE is set, it will overwrite the biotypes
        # unless it is a "deleted model"
        GOOD_BIOTYPE => 'best',

        MERGED_SET_DATABASE => ['REFINED_DB'],
        MERGED_UNIPROT_DATABASE => 'BLAST_DB',

        # list of biotypes to be fetched from the databases
        # If the arrayref is empty, fetch all
        MERGED_BIOTYPES => {
          REFINE_MERGED_DB => [],
        },
        TISSUES_SET_DATABASE => ['REFINED_DB'],
        TISSUES_UNIPROT_DATABASE => 'BLAST_DB',

        # list of biotypes to be fetched from the databases
        # If the arrayref is empty, fetch all
        TISSUES_BIOTYPES => {
          REFINE_TISSUES_DB => [],
        },

        OUTPUT_DATABASE => 'RECOVER_RNASEQ_DB',

        INTRON_BAM_FILE => '',

        # It is supposed to be the number of tissues you have for your species.
        # If the recovered model has less than half of NUM_GROUPS of tissues to
        # support it, the biotype will contain 'low'
        NUM_GROUPS      => 1,

        # Modifies the merge model score
        # Useful to remove merged models with the last/first intron poorly supported
        # Set it higher to favor the merge model
        INTRON_ALLOWANCE      => 23,

        # When 2 tissue model have a deep 'coverage', use the longest translation if
        # score_less_deep > score_deeper*ABUNDANCE_THRESHOLD
        ABUNDANCE_THRESHOLD      => 0.65,

      }, # End default

  );
  return $config{$key};
}


1;
