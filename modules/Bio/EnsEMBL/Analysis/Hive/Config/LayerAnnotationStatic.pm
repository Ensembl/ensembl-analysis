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

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::Config::LayerAnnotationStatic

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
                         get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::GenebuilderStatic',$self->o('uniprot_set'),
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

This is the config file for all layer annotation analysis. You should use it in your Hive configuration file to
specify the parameters of an analysis. You can either choose an existing config or you can create
a new one based on the default hash.

=head1 METHODS

  _master_config_settings: contains all possible parameters

=cut

package Bio::EnsEMBL::Analysis::Hive::Config::LayerAnnotationStatic;

use strict;
use warnings;

use parent ('Bio::EnsEMBL::Analysis::Hive::Config::BaseStatic');

sub _master_config {
  my ($self, $key) = @_;

  my %config = (
    default => [],

# Note cdna represents cdnas/Iso-Seq
# self = pig specific proteins
# human/mouse/mammals = UniProt PE12 proteins aligned with GenBlast
# projection = protein coding GENCODE annotation from human mapped via a pairwise whole genome alignment
# noncanon/int/pseudo = alignments with multiple structural irregularities such as non-canoncial splice sites or frameshift introns (used in pseudogene calling)
# '_1' => [95, 90], -> >= 95 percent coverage, >= 90 percent identity of ORF translation to either realignment of original evidence or a UniProt PE12 protein
# '_2' => [90, 80], -> >= 95 percent coverage, >= 80 percent identity (and not in _1)
# '_3' => [90, 60], -> etc
# '_4' => [90, 40], -> etc
# '_5' => [80, 20], -> etc
# '_6' => [60, 20], -> etc
# '_7' => [0, 0], -> etc

    mammals_basic => [
             {
               ID         => 'LAYER1',
               BIOTYPES   => [
                               'IG_C_gene',
                               'IG_J_gene',
                               'IG_V_gene',
                               'IG_D_gene',
                               'TR_C_gene',
                               'TR_J_gene',
                               'TR_V_gene',
                               'TR_D_gene',
                               'seleno_self',
                             ],
               DISCARD    => 0,
             },


            {
              ID         => 'LAYER2',
              BIOTYPES   => [
                             'cdna2genome',
                             'edited',
                             'gw_gtag',
                             'gw_nogtag',
                             'gw_exo',
                             'cdna_1',
                             'cdna_2',
                             'cdna_3',
                             'rnaseq_merged_1',
                             'rnaseq_merged_2',
                             'rnaseq_merged_3',
                             'rnaseq_tissue_1',
                             'rnaseq_tissue_2',
                             'rnaseq_tissue_3',
                             'self_pe12_sp_1',
                             'self_pe12_tr_1',
                             'self_pe12_sp_2',
                             'self_pe12_tr_2',
                             'projection_1',
                             'projection_2',
                             'projection_3',
                            ],
              FILTER_AGAINST => ['LAYER1'],
              DISCARD    => 0,
            },


            {
              ID         => 'LAYER3',
              BIOTYPES   => [
                             'cdna_4',
                             'rnaseq_merged_4',
                             'rnaseq_tissue_4',
                             'human_pe12_sp_1',
                             'human_pe12_tr_1',
                             'mouse_pe12_sp_1',
                             'mouse_pe12_tr_1',
                             'human_pe12_tr_2',
                             'human_pe12_sp_2',
                             'mouse_pe12_sp_2',
                             'mouse_pe12_tr_2',
                             'genblast_rnaseq_top',
                             'projection_4',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER4',
              BIOTYPES   => [
                             'cdna_5',
                             'rnaseq_merged_5',
                             'rnaseq_tissue_5',
                             'mammals_pe12_sp_1',
                             'mammals_pe12_tr_1',
                             'mammals_pe12_sp_2',
                             'mammals_pe12_tr_2',
                             'self_pe3_sp_1',
                             'self_pe3_tr_1',
                             'genblast_rnaseq_high',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'],
              DISCARD    => 0,
            },


            {
              ID         => 'LAYER5',
              BIOTYPES   => [
                             'human_pe12_sp_3',
                             'human_pe12_tr_3',
                             'mouse_pe12_sp_3',
                             'mouse_pe12_tr_3',
                             'human_pe12_sp_4',
                             'human_pe12_tr_4',
                             'mouse_pe12_sp_4',
                             'mouse_pe12_tr_4',
                             'genblast_rnaseq_medium',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4'],
              DISCARD    => 0,
            },


            {
              ID         => 'LAYER6',
              BIOTYPES   => [
                             'mammals_pe12_sp_3',
                             'mammals_pe12_tr_3',
                             'mammals_pe12_sp_4',
                             'mammals_pe12_tr_4',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5'],
              DISCARD    => 0,
            },

 # This layer is mostly used to pick up potential pseudogenes later
             {
              ID         => 'LAYER7',
              BIOTYPES   => [
                             'cdna_6',
                             'rnaseq_merged_6',
                             'rnaseq_tissue_6',
                             'human_pe12_sp_int_1',
                             'human_pe12_tr_int_1',
                             'human_pe12_sp_int_2',
                             'human_pe12_tr_int_2',
                             'human_pe12_sp_int_3',
                             'human_pe12_tr_int_3',
                             'human_pe12_sp_int_4',
                             'human_pe12_tr_int_4',
                             'mouse_pe12_sp_int_1',
                             'mouse_pe12_tr_int_1',
                             'mouse_pe12_sp_int_2',
                             'mouse_pe12_tr_int_2',
                             'mouse_pe12_sp_int_3',
                             'mouse_pe12_tr_int_3',
                             'mouse_pe12_sp_int_4',
                             'mouse_pe12_tr_int_4',
                             'mammals_pe12_sp_int_1',
                             'mammals_pe12_tr_int_1',
                             'mammals_pe12_sp_int_2',
                             'mammals_pe12_tr_int_2',
                             'mammals_pe12_sp_int_3',
                             'mammals_pe12_tr_int_3',
                             'mammals_pe12_sp_int_4',
                             'mammals_pe12_tr_int_4',
                             'projection_1_noncanon',
                             'projection_2_noncanon',
                             'projection_3_noncanon',
                             'projection_4_noncanon',
                             'projection_1_pseudo',
                             'projection_2_pseudo',
                             'projection_3_pseudo',
                             'projection_4_pseudo',
                           ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6'],
              DISCARD    => 0,
            },

# This layer is just used for completeness when layer filtering lncRNAs
             {
              ID         => 'LAYER8',
              BIOTYPES   => [
                             'cdna_7',
                             'rnaseq_merged_7',
                             'rnaseq_tissue_7',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER7'],
              DISCARD    => 0,
            },


# This is the candidate lncRNA layer, anything here will not have had any significant evidence of
# protein coding overlap and also has no BLAST hit to the longest ORF
             {
              ID         => 'LAYER9',
              BIOTYPES   => [
                              'cdna',
                              'rnaseq_merged',
                              'rnaseq_tissue',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER7','LAYER8'],
              DISCARD    => 0,
            },

    ],

  );
  return $config{$key};
}

1;

