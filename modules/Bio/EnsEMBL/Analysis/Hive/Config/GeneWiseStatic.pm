=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2024] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::Config::GeneWiseStatic

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

=head1 DESCRIPTION

This is the config file for all genewise analysis. You should use it in your Hive configuration file to
specify the parameters of an analysis. You can either choose an existing config or you can create
a new one based on the default hash.

=head1 METHODS

  _master_config_settings: contains all possible parameters

=cut

package Bio::EnsEMBL::Analysis::Hive::Config::GeneWiseStatic;

use strict;
use warnings;

use parent ('Bio::EnsEMBL::Analysis::Hive::Config::BaseStatic');

sub _master_config {
  my ($self, $key) = @_;

  my %config = (
      default => {
        PAF_LOGICNAMES => [], #an array of logic names found in the
        #protein align feature table,
        PAF_MIN_SCORE_THRESHOLD => 0, # 200 for Similarity
        PAF_UPPER_SCORE_THRESHOLD => undef,
        PAF_SOURCE_DB => '#genewise_db#',
        GENE_SOURCE_DB => '#genewise_db#',
        GENEWISE_PARAMETERS => {},
        MINIGENEWISE_PARAMETERS => {},
        MULTIMINIGENEWISE_PARAMETERS => {},
        BLASTMINIGENEWISE_PARAMETERS => {
                                         -fullseq => 1,
                                        },
        EXONERATE_PARAMETERS => {},  # for ExonerateForGenewise / TargettedExonerate
        #example exonerate parameters => {
        #                                 -options => '--model protein2genome
        #                                              --bestn 1
        #                                              --maxintron 700000'
        #                                 }
        FILTER_OBJECT => 'Bio::EnsEMBL::Analysis::Tools::Filter::BlastMiniGenewise',
        #path to object
        FILTER_PARAMS => {},
        BIOTYPES_TO_MASK => [], #empty means no masking
        #specified types will be masked,
        EXON_BASED_MASKING => 1,
        GENE_BASED_MASKING => 0,
        PRE_GENEWISE_MASK => 1,
        POST_GENEWISE_MASK => 1,
        REPEATMASKING => [],
        SOFTMASKING => 0,
        SEQFETCHER_OBJECT => 'Bio::EnsEMBL::Analysis::Tools::SeqFetcher::OBDAIndexSeqFetcher',
        SEQFETCHER_PARAMS => {
                              -db => '#seqfetcher_index#',
                             },
        USE_KILL_LIST => 1,
        LIMIT_TO_FEATURE_RANGE => undef,
        FEATURE_RANGE_PADDING => 0,
        WRITE_REJECTED => 0,
        REJECTED_BIOTYPE => 'rejected',
        OPTIMAL_LENGTH => 1000000,
      },
      targetted_genewise => {
         PAF_LOGICNAMES => ['bestpmatch'],
         PAF_SOURCE_DB => '#target_db#',
         GENE_SOURCE_DB => '#target_db#',
         OUTPUT_DB => '#target_db#',
         OUTPUT_BIOTYPE => '#biotype#',

         GENEWISE_PARAMETERS => {
                                 -endbias => 1,
                                 -matrix => 'BLOSUM80.bla',
                                 -gap => 20,
                                -extension => 8,
                                -splice_model => '#gtag#',
                                },
         MINIGENEWISE_PARAMETERS => {
                                     -terminal_padding => 20000,
                                     -exon_padding => 200,
                                     -minimum_intron => 1000,
                                    },
         MULTIMINIGENEWISE_PARAMETERS =>{
                                          -minimum_feature_length => 50,
                                        },
         FILTER_PARAMS => {
                           -max_exon_length => '20000',
                           -multi_exon_min_coverage => '40',
                           -single_exon_min_coverage => '80',
                           -max_intron_length => '#max_intron_length#',
                           -min_split_coverage => 95,
                           -max_low_complexity => 101,
                          },
         LIMIT_TO_FEATURE_RANGE => 1,
         FEATURE_RANGE_PADDING => 20000,
      },
      similarity => {
         PAF_LOGICNAMES => ['uniprot'],
         PAF_SOURCE_DB => '#target_db#',
         GENE_SOURCE_DB => '#target_db#',
         OUTPUT_DB => '#target_db#',
         OUTPUT_BIOTYPE => '#biotype#',

         GENEWISE_PARAMETERS => {
                                 -endbias => 1,
                                 -matrix => 'BLOSUM80.bla',
                                 -gap => 20,
                                -extension => 8,
                                -splice_model => '#gtag#',
                                },
         MINIGENEWISE_PARAMETERS => {
                                     -terminal_padding => 20000,
                                     -exon_padding => 200,
                                     -minimum_intron => 1000,
                                    },
         MULTIMINIGENEWISE_PARAMETERS =>{
                                          -minimum_feature_length => 50,
                                        },
         PAF_MIN_SCORE_THRESHOLD => 200,
         FILTER_PARAMS => {
                           -max_exon_length => '20000',
                           -multi_exon_min_coverage => '70',
                           -single_exon_min_coverage => '90',
                           -max_intron_length => '#max_intron_length#',
                           -min_split_coverage => '90',
                           -max_low_complexity => '60',
                          },
         REPEATMASKING => '#repeatmasking#',
      },
      targetted_exonerate => {
         PAF_LOGICNAMES => ['bestpmatch'],
         PAF_SOURCE_DB => '#target_db#',
         GENE_SOURCE_DB => '#target_db#',
         OUTPUT_DB => '#target_db#',
         OUTPUT_BIOTYPE => '#biotype#',
         EXONERATE_PARAMETERS => {
                                  -options => '--model protein2genome --bestn 1 --maxintron #max_intron_length#'
                                 },

         FILTER_PARAMS => {
                           -max_exon_length => '20000',
                           -multi_exon_min_coverage => '40',
                           -single_exon_min_coverage => '80',
                           -max_intron_length => '#max_intron_length#',
                           -min_split_coverage => 95,
                           -max_low_complexity => 101,
                          },
         LIMIT_TO_FEATURE_RANGE => 1,
         FEATURE_RANGE_PADDING => 20000,
         GENEWISE_PARAMETERS => {
                                 -endbias => 1,
                                 -matrix => 'BLOSUM80.bla',
                                 -gap => 20,
                                -extension => 8,
                                -splice_model => '#gtag#',
                                },
         MINIGENEWISE_PARAMETERS => {
                                     -terminal_padding => 20000,
                                     -exon_padding => 200,
                                     -minimum_intron => 1000,
                                    },
         MULTIMINIGENEWISE_PARAMETERS =>{
                                          -minimum_feature_length => 50,
                                        },
      },
  );
  return $config{$key};
}

1;
