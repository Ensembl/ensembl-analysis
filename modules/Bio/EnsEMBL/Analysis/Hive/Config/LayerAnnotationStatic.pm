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

    primates_basic => [
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
                             'human_pe12_sp_2',
                             'genblast_rnaseq_top',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER4',
              BIOTYPES   => [
                             'human_pe12_tr_1',
                             'human_pe12_tr_2',
                             'cdna_5',
                             'rnaseq_merged_5',
                             'rnaseq_tissue_5',
                             'self_pe3_sp_1',
                             'self_pe3_tr_1',
                             'genblast_rnaseq_high',
                             'projection_3',
                             'projection_4',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'],
              DISCARD    => 0,
            },


            {
              ID         => 'LAYER5',
              BIOTYPES   => [
                             'genblast_rnaseq_medium',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4'],
              DISCARD    => 0,
            },


             {
              ID         => 'LAYER6',
              BIOTYPES   => [
                             'cdna_6',
                             'rnaseq_merged_6',
                             'rnaseq_tissue_6',
                             'projection_1_noncanon',
                             'projection_2_noncanon',
                             'projection_3_noncanon',
                             'projection_4_noncanon',
                             'projection_1_pseudo',
                             'projection_2_pseudo',
                             'projection_3_pseudo',
                             'projection_4_pseudo',
                           ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER7',
              BIOTYPES   => [
                             'cdna_7',
                             'rnaseq_merged_7',
                             'rnaseq_tissue_7',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER8',
              BIOTYPES   => [
                              'cdna',
                              'rnaseq_merged',
                              'rnaseq_tissue',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER7'],
              DISCARD    => 0,
            },

    ],


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


    fish_basic => [
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
                             'rnaseq_merged_1',
                             'rnaseq_merged_2',
                             'rnaseq_merged_3',
                             'rnaseq_merged_4',
                             'rnaseq_tissue_1',
                             'rnaseq_tissue_2',
                             'rnaseq_tissue_3',
                             'rnaseq_tissue_4',
                             'self_pe12_sp_1',
                             'self_pe12_tr_1',
                             'self_pe12_sp_2',
                             'self_pe12_tr_2',
                             'human_pe12_sp_1',
                             'human_pe12_sp_2',
                             'human_pe12_tr_1',
                             'human_pe12_tr_2',
			     'genblast_rnaseq_top',
			     'genblast_rnaseq_high',
                             'fish_pe12_sp_1',
                             'fish_pe12_tr_1',
                             'fish_pe12_sp_2',
                             'fish_pe12_tr_2',
                            ],
              FILTER_AGAINST => ['LAYER1'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER3',
              BIOTYPES   => [
                             'human_pe12_sp_3',
                             'human_pe12_tr_3',
                             'human_pe12_sp_4',
                             'human_pe12_tr_4',
                             'fish_pe12_sp_3',
                             'fish_pe12_tr_3',
                             'fish_pe12_sp_4',
                             'fish_pe12_tr_4',
                             'mammals_pe12_sp_1',
                             'mammals_pe12_tr_1',
                             'mammals_pe12_sp_2',
                             'mammals_pe12_tr_2',
                             'genblast_rnaseq_medium',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER4',
              BIOTYPES   => [
                             'human_pe12_sp_5',
                             'human_pe12_tr_5',
                             'vert_pe12_sp_1',
                             'vert_pe12_tr_1',
                             'vert_pe12_sp_2',
                             'vert_pe12_tr_2',
                             'vert_pe12_sp_3',
                             'vert_pe12_tr_3',
                             'fish_pe12_sp_5',
                             'fish_pe12_tr_5',
                             'mammals_pe12_sp_3',
                             'mammals_pe12_tr_3',
                             'mammals_pe12_sp_4',
                             'mammals_pe12_tr_4',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'],
              DISCARD    => 0,
            },


            {
              ID         => 'LAYER5',
              BIOTYPES   => [
		             'rnaseq_merged_5',
                             'rnaseq_tissue_5',
                             'rnaseq_merged_6',
                             'rnaseq_tissue_6',
                             'human_pe12_sp_6',
                             'human_pe12_tr_6',
                             'fish_pe12_sp_6',
                             'fish_pe12_tr_6',
                             'genblast_rnaseq_low',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER6',
              BIOTYPES   => [
                             'rnaseq_merged_7',
                             'rnaseq_tissue_7',
                             'mammals_pe12_sp_5',
                             'mammals_pe12_tr_5',
                             'mammals_pe12_sp_6',
                             'mammals_pe12_tr_6',
                             'mammals_pe12_sp_7',
                             'mammals_pe12_tr_7',
                             'human_pe12_sp_7',
                             'human_pe12_tr_7',
                             'vert_pe12_sp_4',
                             'vert_pe12_tr_4',
                             'vert_pe12_sp_5',
                             'vert_pe12_tr_5',
                             'vert_pe12_sp_6',
                             'vert_pe12_tr_6',
                             'vert_pe12_sp_7',
                             'vert_pe12_tr_7',
                             'fish_pe12_sp_7',
                             'fish_pe12_tr_7',
                             'genblast_rnaseq_weak',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER7',
              BIOTYPES   => [
                             'projection_1',
                             'projection_2',
                             'projection_3',
                             'projection_4',
                             'projection_5',
                             'projection_6',
                             'projection_7',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER8',
              BIOTYPES   => [
                             'fish_pe12_sp_int_1',
                             'fish_pe12_tr_int_1',
                             'mammals_pe12_sp_int_1',
                             'mammals_pe12_tr_int_1',
                             'vert_pe12_sp_int_1',
                             'vert_pe12_tr_int_1',
                             'human_pe12_sp_int_1',
                             'human_pe12_tr_int_1',
                             'fish_pe12_sp_int_2',
                             'fish_pe12_tr_int_2',
                             'mammals_pe12_sp_int_2',
                             'mammals_pe12_tr_int_2',
                             'vert_pe12_sp_int_2',
                             'vert_pe12_tr_int_2',
                             'human_pe12_sp_int_2',
                             'human_pe12_tr_int_2',
                             'fish_pe12_sp_int_3',
                             'fish_pe12_tr_int_3',
                             'mammals_pe12_sp_int_3',
                             'mammals_pe12_tr_int_3',
                             'vert_pe12_sp_int_3',
                             'vert_pe12_tr_int_3',
                             'human_pe12_sp_int_3',
                             'human_pe12_tr_int_3',
                             'fish_pe12_sp_int_4',
                             'fish_pe12_tr_int_4',
                             'mammals_pe12_sp_int_4',
                             'mammals_pe12_tr_int_4',
                             'vert_pe12_sp_int_4',
                             'vert_pe12_tr_int_4',
                             'human_pe12_sp_int_4',
                             'human_pe12_tr_int_4',
                             'fish_pe12_sp_int_5',
                             'fish_pe12_tr_int_5',
                             'mammals_pe12_sp_int_5',
                             'mammals_pe12_tr_int_5',
                             'vert_pe12_sp_int_5',
                             'vert_pe12_tr_int_5',
                             'human_pe12_sp_int_5',
                             'human_pe12_tr_int_5',
                             'fish_pe12_sp_int_6',
                             'fish_pe12_tr_int_6',
                             'mammals_pe12_sp_int_6',
                             'mammals_pe12_tr_int_6',
                             'vert_pe12_sp_int_6',
                             'vert_pe12_tr_int_6',
                             'human_pe12_sp_int_6',
                             'human_pe12_tr_int_6',
                             'fish_pe12_sp_int_7',
                             'fish_pe12_tr_int_7',
                             'mammals_pe12_sp_int_7',
                             'mammals_pe12_tr_int_7',
                             'vert_pe12_sp_int_7',
                             'vert_pe12_tr_int_7',
                             'human_pe12_sp_int_7',
                             'human_pe12_tr_int_7',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER7'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER9',
              BIOTYPES   => [
                              'rnaseq_merged',
                              'rnaseq_tissue',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER7','LAYER8'],
              DISCARD    => 0,
            },


    ],


    fish_complete => [
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
              FILTER_AGAINST => ['LAYER1'],
              BIOTYPES   => [
                             'cdna2genome',
                             'edited',
                             'gw_gtag',
                             'gw_nogtag',
                             'gw_exo',
                             'projection_95',
                             'projection_80',
                             'rnaseq_95',
                             'rnaseq_80',
                             'self_pe12_sp_95',
                             'self_pe12_tr_95',
                             'self_pe12_sp_80',
                             'self_pe12_tr_80',
                             'fish_pe12_sp_95',
                             'fish_pe12_tr_95',
                            ],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER3',
              BIOTYPES   => [
                             'fish_pe12_sp_80',
                             'fish_pe12_tr_80',
                             'human_pe12_sp_95',
                             'human_pe12_tr_95',
                             'mouse_pe12_sp_95',
                             'mouse_pe12_tr_95',
                             'self_pe3_sp_95',
                             'self_pe3_tr_95',
                             'vert_pe12_sp_95',
                             'vert_pe12_tr_95',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER4',
              BIOTYPES   => [
                             'human_pe12_sp_80',
                             'human_pe12_tr_80',
                             'mouse_pe12_sp_80',
                             'mouse_pe12_tr_80',
                             'vert_pe12_sp_80',
                             'vert_pe12_tr_80',
                             'mammals_pe12_sp_95',
                             'mammals_pe12_tr_95',
                             'mammals_pe12_sp_80',
                             'mammals_pe12_tr_80',
                             'projection_50',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'],
              DISCARD    => 0,
            },

    ],


    distant_vertebrate => [
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
                              'rnaseq_merged_1',
                              'rnaseq_merged_2',
                              'rnaseq_merged_3',
                              'rnaseq_tissue_1',
                              'rnaseq_tissue_2',
                              'rnaseq_tissue_3',
                              'genblast_select_1',
                              'genblast_select_2',
                              'vert_pe12_sp_1',
                              'vert_pe12_sp_2',
                              'projection_1',
                              'projection_2',
                              'self_pe12_sp_1',
                              'self_pe12_sp_2',
                              'self_pe12_tr_1',
                              'self_pe12_tr_2',
                           ],
              FILTER_AGAINST => ['LAYER1'],

              DISCARD    => 0,
            },


            {
              ID         => 'LAYER3',
              BIOTYPES   => [
                              'rnaseq_merged_4',
                              'rnaseq_tissue_4',
                              'genblast_select_3',
                              'vert_pe12_sp_3',
                              'vert_pe12_tr_1',
                              'vert_pe12_tr_2',
                              'projection_3',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2'],
              DISCARD    => 0,
            },


            {
              ID         => 'LAYER4',
              BIOTYPES   => [
                              'rnaseq_merged_5',
                              'rnaseq_tissue_5',
                              'genblast_select_4',
                              'vert_pe12_sp_4',
                              'vert_pe12_tr_3',
                              'projection_4',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'],
              DISCARD    => 0,
            },


            {
              ID         => 'LAYER5',
              BIOTYPES   => [
                              'rnaseq_merged_6',
                              'rnaseq_tissue_6',
                              'genblast_select_5',
                              'vert_pe12_sp_5',
                              'vert_pe12_tr_4',
                              'projection_5',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER6',
              BIOTYPES   => [
                              'genblast_select_6',
                              'vert_pe12_sp_6',
                              'projection_6',
                              'vert_pe12_tr_6',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER7',
              BIOTYPES   => [
                              'rnaseq_merged_7',
                              'rnaseq_tissue_7',
                              'genblast_select_7',
                              'vert_pe12_sp_7',
                              'projection_7',
                              'vert_pe12_tr_7',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6'],
              DISCARD    => 0,
            },

    ],


  reptiles_basic => [
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
                             'rnaseq_merged_1',
                             'rnaseq_merged_2',
                             'rnaseq_merged_3',
                             'rnaseq_merged_4',
                             'rnaseq_tissue_1',
                             'rnaseq_tissue_2',
                             'rnaseq_tissue_3',
                             'rnaseq_tissue_4',
                             'self_pe12_sp_1',
                             'self_pe12_tr_1',
                             'self_pe12_sp_2',
                             'self_pe12_tr_2',
                            ],
              FILTER_AGAINST => ['LAYER1'],
              DISCARD    => 0,
            },


           {
              ID         => 'LAYER3',
              BIOTYPES   => [
                             'projection_1',
                             'projection_2',
                             'rnaseq_merged_5',
                             'rnaseq_tissue_5',
                             'reptiles_pe12_sp_1',
                             'reptiles_pe12_tr_1',
                             'reptiles_pe12_tr_2',
                             'reptiles_pe12_sp_2',
                             'genblast_rnaseq_top',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2'],
              DISCARD    => 0,
           },

           {
              ID         => 'LAYER4',
              BIOTYPES   => [
                             'human_pe12_sp_1',
                             'human_pe12_tr_1',
                             'human_pe12_tr_2',
                             'human_pe12_sp_2',
                             'aves_pe12_sp_1',
                             'aves_pe12_tr_1',
                             'aves_pe12_tr_2',
                             'aves_pe12_sp_2',
                             'amphibians_sp_1',
                             'amphibians_tr_1',
                             'amphibians_sp_2',
                             'amphibians_tr_2',
                             'reptiles_sp_3',
                             'reptiles_tr_3',
                             'reptiles_sp_4',
                             'reptiles_tr_4',
                             'genblast_rnaseq_high',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'],
              DISCARD    => 0,
            },

           {
              ID         => 'LAYER5',
              BIOTYPES   => [
                             'projection_3',
                             'projection_4',
                             'mammals_pe12_sp_1',
                             'mammals_pe12_tr_1',
                             'mammals_pe12_sp_2',
                             'mammals_pe12_tr_2',
                             'self_pe3_sp_1',
                             'self_pe3_tr_1',
                             'genblast_rnaseq_medium',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4'],
              DISCARD    => 0,
            },


	   {
              ID         => 'LAYER6',
              BIOTYPES   => [
                             'human_pe12_sp_3',
                             'human_pe12_tr_3',
                             'human_pe12_sp_4',
                             'human_pe12_tr_4',
                             'aves_pe12_sp_3',
                             'aves_pe12_tr_3',
                             'aves_pe12_sp_4',
                             'aves_pe12_tr_4',
                             'amphibians_sp_3',
                             'amphibians_tr_3',
                             'amphibians_sp_4',
                             'amphibians_tr_4',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5'],
              DISCARD    => 0,
            },


	   {
              ID         => 'LAYER7',
              BIOTYPES   => [
                             'mammals_pe12_sp_3',
                             'mammals_pe12_tr_3',
                             'mammals_pe12_sp_4',
                             'mammals_pe12_tr_4',
                             'human_pe12_sp_5',
                             'human_pe12_tr_5',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6'],
              DISCARD    => 0,
            },


            {
              ID         => 'LAYER8',
              BIOTYPES   => [
                             'human_pe12_sp_6',
                             'human_pe12_tr_6',
                             'mammals_pe12_sp_5',
                             'mammals_pe12_tr_5',
                             'aves_pe12_sp_5',
                             'aves_pe12_tr_5',
                             'amphibians_sp_5',
                             'amphibians_tr_5',
                             'reptiles_sp_5',
                             'reptiles_tr_5',
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
                             'amphibians_pe12_sp_int_1',
                             'amphibians_pe12_tr_int_1',
                             'amphibians_pe12_sp_int_2',
                             'amphibians_pe12_tr_int_2',
                             'amphibians_pe12_sp_int_3',
                             'amphibians_pe12_tr_int_3',
                             'amphibians_pe12_sp_int_4',
                             'amphibians_pe12_tr_int_4',
                             'aves_pe12_sp_int_1',
                             'aves_pe12_tr_int_1',
                             'aves_pe12_sp_int_2',
                             'aves_pe12_tr_int_2',
                             'aves_pe12_sp_int_3',
                             'aves_pe12_tr_int_3',
                             'aves_pe12_sp_int_4',
                             'aves_pe12_tr_int_4',
                             'reptiles_pe12_sp_int_1',
                             'reptiles_pe12_tr_int_1',
                             'reptiles_pe12_sp_int_2',
                             'reptiles_pe12_tr_int_2',
                             'reptiles_pe12_sp_int_3',
                             'reptiles_pe12_tr_int_3',
                             'reptiles_pe12_sp_int_4',
                             'reptiles_pe12_tr_int_4',
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
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER7'],
              DISCARD    => 0,
            },

           {
              ID         => 'LAYER9',
              BIOTYPES   => [
                             'rnaseq_merged_7',
                             'rnaseq_tissue_7',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER7','LAYER8'],
              DISCARD    => 0,
            },

		  {
              ID         => 'LAYER10',
              BIOTYPES   => [
                              'rnaseq_merged',
                              'rnaseq_tissue',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER7','LAYER8','LAYER9'],
              DISCARD    => 0,
            },

    ],



  birds_basic => [
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
                             'cdna_4',
                             'rnaseq_merged_1',
                             'rnaseq_merged_2',
                             'rnaseq_merged_3',
                             'rnaseq_merged_4',
                             'rnaseq_tissue_1',
                             'rnaseq_tissue_2',
                             'rnaseq_tissue_3',
                             'rnaseq_tissue_4',
                             'self_pe12_sp_1',
                             'self_pe12_tr_1',
                             'self_pe12_sp_2',
                             'self_pe12_tr_2',
                            ],
              FILTER_AGAINST => ['LAYER1'],
              DISCARD    => 0,
            },


           {
              ID         => 'LAYER3',
              BIOTYPES   => [
                             'projection_1',
                             'projection_2',
                             'cdna_5',
                             'rnaseq_merged_5',
                             'rnaseq_tissue_5',
                             'aves_pe12_sp_1',
                             'aves_pe12_tr_1',
                             'aves_pe12_tr_2',
                             'aves_pe12_sp_2',
                             'genblast_rnaseq_top',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2'],
              DISCARD    => 0,
           },

           {
              ID         => 'LAYER4',
              BIOTYPES   => [
                             'human_pe12_sp_1',
                             'human_pe12_tr_1',
                             'human_pe12_tr_2',
                             'human_pe12_sp_2',
                             'aves_pe12_sp_3',
                             'aves_pe12_tr_3',
                             'aves_pe12_sp_4',
                             'aves_pe12_tr_4',
                             'reptiles_pe12_sp_1',
                             'reptiles_pe12_tr_1',
                             'reptiles_pe12_tr_2',
                             'reptiles_pe12_sp_2',
                             'genblast_rnaseq_high',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'],
              DISCARD    => 0,
            },

           {
              ID         => 'LAYER5',
              BIOTYPES   => [
                             'projection_3',
                             'projection_4',
                             'mammals_pe12_sp_1',
                             'mammals_pe12_tr_1',
                             'mammals_pe12_sp_2',
                             'mammals_pe12_tr_2',
                             'reptiles_sp_3',
                             'reptiles_tr_3',
                             'reptiles_sp_4',
                             'reptiles_tr_4',
                             'self_pe3_sp_1',
                             'self_pe3_tr_1',
                             'genblast_rnaseq_medium',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4'],
              DISCARD    => 0,
            },


	   {
              ID         => 'LAYER6',
              BIOTYPES   => [
                             'human_pe12_sp_3',
                             'human_pe12_tr_3',
                             'human_pe12_sp_4',
                             'human_pe12_tr_4',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5'],
              DISCARD    => 0,
            },


	   {
              ID         => 'LAYER7',
              BIOTYPES   => [
                             'mammals_pe12_sp_3',
                             'mammals_pe12_tr_3',
                             'mammals_pe12_sp_4',
                             'mammals_pe12_tr_4',
                             'human_pe12_sp_5',
                             'human_pe12_tr_5',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6'],
              DISCARD    => 0,
            },


            {
              ID         => 'LAYER8',
              BIOTYPES   => [
                             'human_pe12_sp_6',
                             'human_pe12_tr_6',
                             'mammals_pe12_sp_5',
                             'mammals_pe12_tr_5',
                             'aves_pe12_sp_5',
                             'aves_pe12_tr_5',
                             'reptiles_sp_5',
                             'reptiles_tr_5',
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
                             'aves_pe12_sp_int_1',
                             'aves_pe12_tr_int_1',
                             'aves_pe12_sp_int_2',
                             'aves_pe12_tr_int_2',
                             'aves_pe12_sp_int_3',
                             'aves_pe12_tr_int_3',
                             'aves_pe12_sp_int_4',
                             'aves_pe12_tr_int_4',
                             'reptiles_pe12_sp_int_1',
                             'reptiles_pe12_tr_int_1',
                             'reptiles_pe12_sp_int_2',
                             'reptiles_pe12_tr_int_2',
                             'reptiles_pe12_sp_int_3',
                             'reptiles_pe12_tr_int_3',
                             'reptiles_pe12_sp_int_4',
                             'reptiles_pe12_tr_int_4',
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
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER7'],
              DISCARD    => 0,
            },

           {
              ID         => 'LAYER9',
              BIOTYPES   => [
                             'cdna_7',
                             'rnaseq_merged_7',
                             'rnaseq_tissue_7',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER7','LAYER8'],
              DISCARD    => 0,
            },

		  {
              ID         => 'LAYER10',
              BIOTYPES   => [
                              'cdna',
                              'rnaseq_merged',
                              'rnaseq_tissue',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER7','LAYER8','LAYER9'],
              DISCARD    => 0,
            },

    ],


    hemiptera_basic => [
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
			     'genblast_rnaseq_top',
			     'genblast_rnaseq_high',
			     'genblast_rnaseq_medium',
                             'rnaseq_merged_1',
                             'rnaseq_merged_2',
                             'rnaseq_merged_3',
                             'rnaseq_merged_4',
                             'rnaseq_merged_5',
                             'rnaseq_tissue_1',
                             'rnaseq_tissue_2',
                             'rnaseq_tissue_3',
                             'rnaseq_tissue_4',
                             'rnaseq_tissue_5',
                             'self_pe12_sp_1',
                             'self_pe12_tr_1',
                             'self_pe12_sp_2',
                             'self_pe12_tr_2',
                             'genblast_rnaseq_top',
                             'genblast_rnaseq_high',
                            ],
              FILTER_AGAINST => ['LAYER1'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER3',
              BIOTYPES   => [
                             'rnaseq_merged_6',
                             'rnaseq_merged_7',
                             'rnaseq_tissue_6',
                             'rnaseq_tissue_7',
			     'genblast_rnaseq_low',
			     'genblast_rnaseq_weak',
                             'self_pe12_sp_3',
                             'self_pe12_tr_3',
                             'self_pe12_sp_4',
                             'self_pe12_tr_4',
                             'hemiptera_pe12_sp_1',
                             'pisum_pe12_sp_1',
                             'drosophila_sp_1',
                             'flies_sp_1',
                             'hemiptera_pe12_sp_2',
                             'pisum_pe12_sp_2',
                             'drosophila_sp_2',
                             'flies_sp_2',
                             'hemiptera_pe12_sp_3',
                             'pisum_pe12_sp_3',
                             'drosophila_sp_3',
                             'flies_sp_3',
                             'hemiptera_pe12_sp_4',
                             'pisum_pe12_sp_4',
                             'drosophila_sp_4',
                             'flies_sp_4',
                             'genblast_rnaseq_medium',
                             'genblast_rnaseq_low',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER4',
              BIOTYPES   => [
                             'hemiptera_pe12_tr_1',
                             'pisum_pe12_tr_1',
                             'drosophila_tr_1',
                             'flies_tr_1',
                             'hemiptera_pe12_tr_2',
                             'pisum_pe12_tr_2',
                             'drosophila_tr_2',
                             'flies_tr_2',
                             'hemiptera_pe12_tr_3',
                             'pisum_pe12_tr_3',
                             'drosophila_tr_3',
                             'flies_tr_3',
                             'hemiptera_pe12_tr_4',
                             'pisum_pe12_tr_4',
                             'drosophila_tr_4',
                             'flies_tr_4',
                             'genblast_rnaseq_weak',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'],
              DISCARD    => 0,
            },


            {
              ID         => 'LAYER5',
              BIOTYPES   => [
                              'rnaseq_merged',
                              'rnaseq_tissue',
			     'hemiptera_pe3_tr_1',
			     'hemiptera_pe3_tr_2',
			     'hemiptera_pe3_tr_3',
			     'hemiptera_pe3_tr_4',
			     'hemiptera_pe3_tr_5',
			     'hemiptera_pe3_tr_6',
			     'hemiptera_pe3_tr_7',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4'],
              DISCARD    => 0,
            },

	   {
              ID         => 'LAYER6',
              BIOTYPES   => [
			     'hemiptera_pe12_tr_int_1',
			     'hemiptera_pe12_tr_int_2',
			     'hemiptera_pe12_tr_int_3',
			     'hemiptera_pe12_tr_int_4',
			     'hemiptera_pe12_tr_int_5',
			     'hemiptera_pe12_tr_int_6',
			     'hemiptera_pe12_tr_int_7',
			     'hemiptera_pe3_tr_int_1',
			     'hemiptera_pe3_tr_int_2',
			     'hemiptera_pe3_tr_int_3',
			     'hemiptera_pe3_tr_int_4',
			     'hemiptera_pe3_tr_int_5',
			     'hemiptera_pe3_tr_int_6',
			     'hemiptera_pe3_tr_int_7',
			     'pisum_pe12_tr_int_1',
			     'pisum_pe12_tr_int_2',
			     'pisum_pe12_tr_int_3',
			     'pisum_pe12_tr_int_4',
			     'pisum_pe12_tr_int_5',
			     'pisum_pe12_tr_int_6',
			     'pisum_pe12_tr_int_7',
			     'self_pe12_tr_int_1',
			     'self_pe12_tr_int_2',
			     'self_pe12_tr_int_3',
			     'self_pe12_tr_int_4',
			     'self_pe12_tr_int_5',
			     'self_pe12_tr_int_6',
			     'self_pe12_tr_int_7'
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4'],
              DISCARD    => 0,
            },

    ],


    insects_basic => [
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
                             'rnaseq_merged_1',
                             'rnaseq_merged_2',
                             'rnaseq_merged_3',
                             'rnaseq_merged_4',
                             'rnaseq_merged_5',
                             'rnaseq_tissue_1',
                             'rnaseq_tissue_2',
                             'rnaseq_tissue_3',
                             'rnaseq_tissue_4',
                             'rnaseq_tissue_5',
                             'self_pe12_sp_1',
                             'self_pe12_tr_1',
                             'self_pe12_sp_2',
                             'self_pe12_tr_2',
                             'dicondylia_pe12_sp_1',
                             'dicondylia_pe12_tr_1',
                             'dicondylia_pe12_sp_2',
                             'dicondylia_pe12_tr_2',
                             'monocondylia_pe12_sp_1',
                             'monocondylia_pe12_tr_1',
                             'monocondylia_pe12_sp_2',
                             'monocondylia_pe12_tr_2',
                             'human_pe12_sp_1',
                             'human_pe12_tr_1',
                             'human_pe12_sp_2',
                             'human_pe12_tr_2',
                             'human_pe12_sp_3',
                             'human_pe12_tr_3',
                             'genblast_rnaseq_top',
                             'genblast_rnaseq_high',
                             ],
              FILTER_AGAINST => ['LAYER1'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER3',
              BIOTYPES   => [
                             'rnaseq_merged_6',
                             'rnaseq_merged_7',
                             'rnaseq_tissue_6',
                             'rnaseq_tissue_7',
                             'self_pe12_sp_3',
                             'self_pe12_tr_3',
                             'self_pe12_sp_4',
                             'self_pe12_tr_4',
                             'dicondylia_pe12_sp_3',
                             'dicondylia_pe12_tr_3',
                             'dicondylia_pe12_sp_4',
                             'dicondylia_pe12_tr_4',
                             'monocondylia_pe12_sp_3',
                             'monocondylia_pe12_tr_3',
                             'monocondylia_pe12_sp_4',
                             'monocondylia_pe12_tr_4',
                             'human_pe12_sp_4',
                             'human_pe12_tr_4',
                             'human_pe12_sp_5',
                             'human_pe12_tr_5',
                             'genblast_rnaseq_medium',
                             'genblast_rnaseq_low',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER4',
              BIOTYPES   => [
                             'self_pe3_sp_1',
                             'self_pe3_tr_1',
                             'self_pe3_sp_2',
                             'self_pe3_tr_2',
                             'dicondylia_pe3_sp_1',
                             'dicondylia_pe3_tr_1',
                             'dicondylia_pe3_sp_2',
                             'dicondylia_pe3_tr_2',
                             'monocondylia_pe3_sp_1',
                             'monocondylia_pe3_tr_1',
                             'monocondylia_pe3_sp_2',
                             'monocondylia_pe3_tr_2',
                             'human_pe12_sp_6',
                             'human_pe12_tr_6',
                             'genblast_rnaseq_weak',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER5',
              BIOTYPES   => [
                             'self_pe12_sp_5',
                             'self_pe12_tr_5',
                             'self_pe12_sp_6',
                             'self_pe12_tr_6',
                             'dicondylia_pe12_sp_5',
                             'dicondylia_pe12_tr_5',
                             'dicondylia_pe12_sp_6',
                             'dicondylia_pe12_tr_6',
                             'monocondylia_pe12_sp_5',
                             'monocondylia_pe12_tr_5',
                             'monocondylia_pe12_sp_6',
                             'monocondylia_pe12_tr_6',
                             'dicondylia_pe3_sp_3',
                             'dicondylia_pe3_tr_3',
                             'dicondylia_pe3_sp_4',
                             'dicondylia_pe3_tr_4',
                             'monocondylia_pe3_sp_3',
                             'monocondylia_pe3_tr_3',
                             'monocondylia_pe3_sp_4',
                             'monocondylia_pe3_tr_4',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4'],
              DISCARD    => 0,
            },


            {
              ID         => 'LAYER6',
              BIOTYPES   => [
                              'rnaseq_merged',
                              'rnaseq_tissue',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5'],
              DISCARD    => 0,
            },

    ],


    self_patch => [
            {
              ID         => 'LAYER1',
              BIOTYPES   => [
                             'self_pe12_sp_95',
                             'self_pe12_tr_95',
                            ],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER2',
              BIOTYPES   => [
                              'self_frag_pe12_sp_95',
                              'self_frag_pe12_tr_95',
                            ],
              FILTER_AGAINST => ['LAYER1'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER3',
              BIOTYPES   => [
                             'self_pe12_sp_80',
                             'self_pe12_tr_80',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER4',
              BIOTYPES   => [
                             'self_frag_pe12_sp_80',
                             'self_frag_pe12_tr_80',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'],
              DISCARD    => 0,
            },
    ],

  );
  return $config{$key};
}

1;

