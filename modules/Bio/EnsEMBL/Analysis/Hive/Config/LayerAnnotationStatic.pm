=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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
                             'realign_95',
                             'realign_80',
                             'rnaseq_merged_95',
                             'rnaseq_merged_80',
                             'self_pe12_sp_95',
                             'self_pe12_tr_95',
                             'self_pe12_sp_80',
                             'self_pe12_tr_80',
                             'human_pe12_sp_95',
                             'human_pe12_tr_95',
                             'primates_pe12_sp_95',
                             'primates_pe12_tr_95',
                             'mammals_pe12_sp_95',
                             'mammals_pe12_tr_95',
                            ],
              FILTER_AGAINST => ['LAYER1'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER3',
              BIOTYPES   => [
                             'rnaseq_tissue_95',
                             'human_pe12_sp_80',
                             'human_pe12_tr_80',
                             'primates_pe12_sp_80',
                             'primates_pe12_tr_80',
                             'mammals_pe12_sp_80',
                             'mammals_pe12_tr_80',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER4',
              BIOTYPES   => [
                             'rnaseq_tissue_80',
                             'realign_50',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER5',
              BIOTYPES   => [
                             'human_pe12_sp_50',
                             'human_pe12_tr_50',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4'],
              DISCARD    => 0,
            },

    ],

    rodents_basic => [

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
              BIOTYPES   => ['realign_95',
                             'rnaseq_95',
                             'rnaseq_80',
                             'self_pe12_sp_95',
                             'self_pe12_sp_80',
                             'mouse_pe12_sp_95',
                             'rodents_pe12_sp_95',
                             'human_pe12_sp_95',
                             'cdna2genome',
                             'edited',
                             'gw_gtag',
                             'gw_nogtag',
                             'gw_exo',
                            ],
              FILTER_AGAINST => ['LAYER1'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER3',
              BIOTYPES   => ['self_pe12_tr_95',
                             'mouse_pe12_tr_95',
                             'rodents_pe12_tr_95',
                             'human_pe12_tr_95',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2'],
              DISCARD    => 0,

            },

             {
              ID         => 'LAYER4',
              BIOTYPES   => ['realign_80',
                             'mouse_pe12_sp_80',
                             'rodents_pe12_sp_80',
                             'human_pe12_sp_80',
                             'mammals_pe12_sp_95',
                             'vert_pe12_sp_95',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER5',
              BIOTYPES   => ['self_pe12_tr_80',
                             'mouse_pe12_tr_80',
                             'rodents_pe12_tr_80',
                             'human_pe12_tr_80',
                             'mammals_pe12_tr_95',
                             'vert_pe12_tr_95',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER6',
              BIOTYPES   => ['rodents_pe3_sp_95',
                             'rodents_pe3_tr_95',
                             'mammals_pe12_sp_80',
                             'vert_pe12_sp_80',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER7',
              BIOTYPES   => ['realign_50',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6'],
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
                             'realign_1',
                             'realign_2',
                             'rnaseq_merged_1',
                             'rnaseq_merged_2',
                             'rnaseq_merged_3',
                             'rnaseq_tissue_1',
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
                             'realign_3',
                             'rnaseq_tissue_2',
                             'rnaseq_tissue_3',
                             'human_pe12_sp_1',
                             'human_pe12_tr_1',
                             'mouse_pe12_sp_1',
                             'mouse_pe12_tr_1',
                             'mammals_pe12_sp_1',
                             'mammals_pe12_tr_1',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER4',
              BIOTYPES   => [
                             'realign_4',
                             'self_pe3_sp_1',
                             'self_pe3_tr_1',
                             'human_pe12_tr_2',
                             'human_pe12_sp_2',
                             'mouse_pe12_sp_2',
                             'mouse_pe12_tr_2',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'],
              DISCARD    => 0,
            },


            {
              ID         => 'LAYER5',
              BIOTYPES   => [
                             'mammals_pe12_sp_2',
                             'mammals_pe12_tr_2',
                             'human_pe12_sp_3',
                             'human_pe12_tr_3',
                             'mouse_pe12_sp_3',
                             'mouse_pe12_tr_3',
                             'vert_pe12_sp_1',
                             'vert_pe12_tr_1',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER6',
              BIOTYPES   => [
                             'rnaseq_merged_4',
                             'mammals_pe12_sp_3',
                             'mammals_pe12_tr_3',
                             'vert_pe12_sp_3',
                             'vert_pe12_tr_3',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER7',
              BIOTYPES   => [
                             'human_pe12_sp_4',
                             'human_pe12_tr_4',
                             'mouse_pe12_sp_4',
                             'mouse_pe12_tr_4',
                             'rnaseq_merged_5',
                             'rnaseq_tissue_5',
                             'human_pe12_sp_5',
                             'mouse_pe12_sp_5',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER8',
              BIOTYPES   => [
                             'human_pe12_tr_5',
                             'mouse_pe12_tr_5',
                             'mammals_pe12_sp_4',
                             'mammals_pe12_tr_4',
                             'vert_pe12_sp_4',
                             'vert_pe12_tr_4',
                             'realign_5',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER7'],
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
                             'realign_1',
                             'realign_2',
                             'rnaseq_merged_1',
                             'rnaseq_merged_2',
                             'rnaseq_merged_3',
                             'rnaseq_tissue_1',
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
                             'realign_3',
                             'rnaseq_tissue_2',
                             'rnaseq_tissue_3',
                             'fish_pe12_sp_1',
                             'fish_pe12_tr_1',
                             'fish_pe12_sp_2',
                             'fish_pe12_tr_2',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER4',
              BIOTYPES   => [
                             'realign_4',
                             'fish_pe12_sp_3',
                             'fish_pe12_tr_3',
                             'human_pe12_sp_1',
                             'human_pe12_tr_1',
                             'mouse_pe12_sp_1',
                             'mouse_pe12_tr_1',
                             'vert_pe12_sp_1',
                             'vert_pe12_tr_1',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'],
              DISCARD    => 0,
            },


            {
              ID         => 'LAYER5',
              BIOTYPES   => [
                             'human_pe12_sp_2',
                             'human_pe12_tr_2',
                             'mouse_pe12_sp_2',
                             'mouse_pe12_tr_2',
                             'vert_pe12_sp_2',
                             'vert_pe12_tr_2',
                             'mammals_pe12_sp_1',
                             'mammals_pe12_tr_1',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER6',
              BIOTYPES   => [
                             'human_pe12_sp_3',
                             'human_pe12_tr_3',
                             'mouse_pe12_sp_3',
                             'mouse_pe12_tr_3',
                             'vert_pe12_sp_3',
                             'vert_pe12_tr_3',
                             'mammals_pe12_sp_2',
                             'mammals_pe12_tr_2',
                             'mammals_pe12_sp_3',
                             'mammals_pe12_tr_3',
                             'rnaseq_merged_4',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER7',
              BIOTYPES   => [
                             'fish_pe12_sp_4',
                             'fish_pe12_tr_4',
                             'human_pe12_sp_4',
                             'human_pe12_tr_4',
                             'mouse_pe12_sp_4',
                             'mouse_pe12_tr_4',
                             'rnaseq_merged_5',
                             'rnaseq_tissue_4',
                             'rnaseq_tissue_5',
                             'human_pe12_sp_5',
                             'mouse_pe12_sp_5',
                             'vert_pe12_sp_4',
                             'vert_pe12_tr_4',

                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER8',
              BIOTYPES   => [
                             'fish_pe12_sp_5',
                             'fish_pe12_tr_5',
                             'human_pe12_tr_5',
                             'mouse_pe12_tr_5',
                             'mammals_pe12_sp_4',
                             'mammals_pe12_tr_4',
                             'mammals_pe12_sp_5',
                             'mammals_pe12_tr_5',
                             'vert_pe12_sp_5',
                             'vert_pe12_tr_5',
                             'realign_5',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER7'],
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
                             'realign_95',
                             'realign_80',
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
                             'realign_50',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'],
              DISCARD    => 0,
            },

    ],


    distant_vertebrate => [
             {
              ID         => 'LAYER1',
              BIOTYPES   => [
                             'rnaseq_merged_1',
                             'rnaseq_merged_2',
                             'rnaseq_merged_90_70',
                             'rnaseq_merged_90_60',
                             'rnaseq_tissue_90_70',
                             'rnaseq_tissue_90_60',
                             'rnaseq_merged_80_50',
                             'rnaseq_tissue_80_50',
                             'rnaseq_tissue_1',
                             'rnaseq_tissue_2',
                            ],
              DISCARD    => 0,
            },


             {
              ID         => 'LAYER2',
              BIOTYPES   => [
                              'genblast_1',
                              'genblast_2',
                              'genblast_90_70',
                              'genblast_90_60',
                            ],
              FILTER_AGAINST => ['LAYER1'],

              DISCARD    => 0,
            },

            {
              ID         => 'LAYER3',
              BIOTYPES   => [
                              'rnaseq_merged_70_40',
                              'rnaseq_tissue_70_40',
                              'genblast_80_50'
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER4',
              BIOTYPES   => [
                              'genblast_70_30',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'],
              DISCARD    => 0,
            },


            {
              ID         => 'LAYER5',
              BIOTYPES   => [
                              'genblast_60_20'
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER6',
              BIOTYPES   => [
                             'rnaseq_merged_50_25',
                             'rnaseq_tissue_50_25',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER7',
              BIOTYPES   => [
                              'genblast_0_0',
                              'rnaseq_tissue_0_0',
                              'rnaseq_merged_0_0',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6'],
              DISCARD    => 0,
            },

    ],

    bird_basic => [
             {
              ID         => 'LAYER1',
              BIOTYPES   => [
                             'cdna2genome',
                             'edited',
                             'gw_gtag',
                             'gw_nogtag',
                             'gw_exo',
                             'realign_1',
                             'realign_2',
                             'rnaseq_merged_1',
                             'rnaseq_merged_2',
                             'rnaseq_merged_3',
                             'rnaseq_tissue_1',
                             'self_pe12_sp_1',
                             'self_pe12_tr_1',
                             'self_pe12_sp_2',
                             'self_pe12_tr_2',
                            ],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER2',
              BIOTYPES   => [
                             'rnaseq_tissue_2',
                             'rnaseq_tissue_3',
                             'bird_pe12_sp_1',
                             'bird_pe12_tr_1',
                             'bird_pe12_sp_2',
                             'bird_pe12_tr_2',
                            ],
              FILTER_AGAINST => ['LAYER1'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER3',
              BIOTYPES   => [
                             'bird_pe12_sp_3',
                             'bird_pe12_tr_3',
                             'human_pe12_sp_1',
                             'human_pe12_tr_1',
                             'vert_pe12_sp_1',
                             'vert_pe12_tr_1',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER4',
              BIOTYPES   => [
                             'human_pe12_sp_2',
                             'human_pe12_tr_2',
                             'vert_pe12_sp_2',
                             'vert_pe12_tr_2',
                             'mammals_pe12_sp_1',
                             'mammals_pe12_tr_1',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'],
              DISCARD    => 0,
            },


            {
              ID         => 'LAYER5',
              BIOTYPES   => [
                             'human_pe12_sp_3',
                             'human_pe12_tr_3',
                             'vert_pe12_sp_3',
                             'vert_pe12_tr_3',
                             'mammals_pe12_sp_2',
                             'mammals_pe12_tr_2',
                             'mammals_pe12_sp_3',
                             'mammals_pe12_tr_3',
                             'rnaseq_merged_4',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER6',
              BIOTYPES   => [
                             'bird_pe12_sp_4',
                             'bird_pe12_tr_4',
                             'human_pe12_sp_4',
                             'human_pe12_tr_4',
                             'rnaseq_merged_5',
                             'rnaseq_tissue_4',
                             'rnaseq_tissue_5',
                             'human_pe12_sp_5',
                             'vert_pe12_sp_4',
                             'vert_pe12_tr_4',

                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER7',
              BIOTYPES   => [
                             'bird_pe12_sp_5',
                             'bird_pe12_tr_5',
                             'human_pe12_tr_5',
                             'mammals_pe12_sp_4',
                             'mammals_pe12_tr_4',
                             'mammals_pe12_sp_5',
                             'mammals_pe12_tr_5',
                             'vert_pe12_sp_5',
                             'vert_pe12_tr_5',
                             'realign_5',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6'],
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
    cat => [ { ID         => 'LAYER1', BIOTYPES   => [ 'IG_C_gene', 'IG_J_gene', 'IG_V_gene', 'IG_D_gene', 'TR_C_gene', 'TR_J_gene', 'TR_V_gene', 'TR_D_gene', ], DISCARD    => 0, }, { ID         => 'LAYER2', BIOTYPES   => [ 'realign_1', 'realign_2', 'rnaseq_merged_1', 'rnaseq_merged_2', 'rnaseq_merged_3', 'rnaseq_tissue_1', 'self_pe12_sp_1', 'self_pe12_sp_2', ], FILTER_AGAINST => ['LAYER1'], DISCARD    => 0, }, { ID         => 'LAYER3', BIOTYPES   => [ 'rnaseq_tissue_2', 'rnaseq_tissue_3', 'human_pe12_sp_1', 'mouse_pe12_sp_1', 'mammals_pe12_sp_1', ], FILTER_AGAINST => ['LAYER1','LAYER2'], DISCARD    => 0, }, { ID         => 'LAYER4', BIOTYPES   => [ 'self_pe3_sp_1', 'self_pe3_tr_1', 'human_pe12_sp_2', 'mouse_pe12_sp_2', ], FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'], DISCARD    => 0, }, { ID         => 'LAYER5', BIOTYPES   => [ 'mammals_pe12_sp_2', 'human_pe12_sp_3', 'mouse_pe12_sp_3', 'vert_pe12_sp_1', ], FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4'], DISCARD    => 0, }, { ID         => 'LAYER6', BIOTYPES   => [ 'rnaseq_merged_4', 'mammals_pe12_sp_3', 'vert_pe12_sp_3', ], FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5'], DISCARD    => 0, }, { ID         => 'LAYER7', BIOTYPES   => [ 'human_pe12_sp_4', 'mouse_pe12_sp_4', 'rnaseq_merged_5', 'rnaseq_tissue_5', 'human_pe12_sp_5', 'mouse_pe12_sp_5', ], FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6'], DISCARD    => 0, }, { ID         => 'LAYER8', BIOTYPES   => [ 'human_pe12_tr_5', 'mouse_pe12_tr_5', 'mammals_pe12_sp_4', 'mammals_pe12_tr_4', 'vert_pe12_sp_4', 'vert_pe12_tr_4', 'realign_5', 'human_pe12_tr_2', 'self_pe12_tr_1', 'self_pe12_tr_2', 'human_pe12_tr_1', 'mouse_pe12_tr_1', 'mammals_pe12_tr_1', 'mouse_pe12_tr_2', 'mammals_pe12_tr_2', 'human_pe12_tr_3', 'mouse_pe12_tr_3', 'vert_pe12_tr_1', 'mammals_pe12_tr_3', 'vert_pe12_tr_3', 'human_pe12_tr_4', 'mouse_pe12_tr_4', ], FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER7'], DISCARD    => 0, }, ],
  );
  return $config{$key};
}

genebuild1 [{'BIOTYPES' => ['IG_C_gene','IG_J_gene','IG_V_gene','IG_D_gene','TR_C_gene','TR_J_gene','TR_V_gene','TR_D_gene'],'DISCARD' => 0,'ID' => 'LAYER1'},{'BIOTYPES' => ['realign_1','realign_2','rnaseq_merged_1','rnaseq_merged_2','rnaseq_merged_3','rnaseq_tissue_1','self_pe12_sp_1','self_pe12_tr_1','self_pe12_sp_2','self_pe12_tr_2'],'DISCARD' => 0,'FILTER_AGAINST' => ['LAYER1'],'ID' => 'LAYER2'},{'BIOTYPES' => ['rnaseq_tissue_2','rnaseq_tissue_3','human_pe12_sp_1','human_pe12_tr_1','mouse_pe12_sp_1','mouse_pe12_tr_1','mammals_pe12_sp_1','mammals_pe12_tr_1'],'DISCARD' => 0,'FILTER_AGAINST' => ['LAYER1','LAYER2'],'ID' => 'LAYER3'},{'BIOTYPES' => ['self_pe3_sp_1','self_pe3_tr_1','human_pe12_tr_2','human_pe12_sp_2','mouse_pe12_sp_2','mouse_pe12_tr_2'],'DISCARD' => 0,'FILTER_AGAINST' => ['LAYER1','LAYER2','LAYER3'],'ID' => 'LAYER4'},{'BIOTYPES' => ['mammals_pe12_sp_2','mammals_pe12_tr_2','human_pe12_sp_3','human_pe12_tr_3','mouse_pe12_sp_3','mouse_pe12_tr_3','vert_pe12_sp_1','vert_pe12_tr_1'],'DISCARD' => 0,'FILTER_AGAINST' => ['LAYER1','LAYER2','LAYER3','LAYER4'],'ID' => 'LAYER5'},{'BIOTYPES' => ['rnaseq_merged_4','mammals_pe12_sp_3','mammals_pe12_tr_3','vert_pe12_sp_3','vert_pe12_tr_3'],'DISCARD' => 0,'FILTER_AGAINST' => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5'],'ID' => 'LAYER6'},{'BIOTYPES' => ['human_pe12_sp_4','human_pe12_tr_4','mouse_pe12_sp_4','mouse_pe12_tr_4','rnaseq_merged_5','rnaseq_tissue_5','human_pe12_sp_5','mouse_pe12_sp_5'],'DISCARD' => 0,'FILTER_AGAINST' => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6'],'ID' => 'LAYER7'},{'BIOTYPES' => ['human_pe12_tr_5','mouse_pe12_tr_5','mammals_pe12_sp_4','mammals_pe12_tr_4','vert_pe12_sp_4','vert_pe12_tr_4','realign_5'],'DISCARD' => 0,'FILTER_AGAINST' => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6','LAYER7'],'ID' => 'LAYER8'}]

1;
count(*) | biotype           |
+----------+-------------------+
|     2576 | human_pe12_sp_1   |
|     4596 | human_pe12_sp_2   |
|      490 | human_pe12_sp_3   |
|     5667 | human_pe12_sp_4   |
|      445 | human_pe12_sp_5   |
|     2072 | human_pe12_sp_6   |
|       36 | human_pe12_sp_7   |
|     4672 | human_pe12_tr_1   |
|     9508 | human_pe12_tr_2   |
|     1880 | human_pe12_tr_3   |
|    12295 | human_pe12_tr_4   |
|     1513 | human_pe12_tr_5   |
|     3619 | human_pe12_tr_6   |
|       58 | human_pe12_tr_7   |
|     6505 | mammals_pe12_sp_1 |
|    11501 | mammals_pe12_sp_2 |
|     1282 | mammals_pe12_sp_3 |
|    12724 | mammals_pe12_sp_4 |
|     1288 | mammals_pe12_sp_5 |
|     4495 | mammals_pe12_sp_6 |
|       85 | mammals_pe12_sp_7 |
|    33263 | mammals_pe12_tr_1 |
|    53077 | mammals_pe12_tr_2 |
|     7228 | mammals_pe12_tr_3 |
|    69721 | mammals_pe12_tr_4 |
|     5797 | mammals_pe12_tr_5 |
|    20746 | mammals_pe12_tr_6 |
|      254 | mammals_pe12_tr_7 |
|     3353 | mouse_pe12_sp_1   |
|     7834 | mouse_pe12_sp_2   |
|      868 | mouse_pe12_sp_3   |
|     8694 | mouse_pe12_sp_4   |
|      731 | mouse_pe12_sp_5   |
|     3171 | mouse_pe12_sp_6   |
|       48 | mouse_pe12_sp_7   |
|     4920 | mouse_pe12_tr_1   |
|    12158 | mouse_pe12_tr_2   |
|     2018 | mouse_pe12_tr_3   |
|    17439 | mouse_pe12_tr_4   |
|     1615 | mouse_pe12_tr_5   |
|     6071 | mouse_pe12_tr_6   |
|      126 | mouse_pe12_tr_7   |
|      137 | self_pe12_sp_1    |
|       26 | self_pe12_sp_2    |
|        7 | self_pe12_sp_3    |
|       45 | self_pe12_sp_4    |
|        8 | self_pe12_sp_5    |
|       23 | self_pe12_sp_6    |
|      372 | self_pe12_tr_1    |
|       95 | self_pe12_tr_2    |
|       63 | self_pe12_tr_3    |
|      172 | self_pe12_tr_4    |
|       12 | self_pe12_tr_5    |
|       86 | self_pe12_tr_6    |
|        5 | self_pe12_tr_7    |
|       34 | self_pe3_sp_1     |
|        2 | self_pe3_sp_2     |
|        1 | self_pe3_sp_3     |
|        4 | self_pe3_sp_4     |
|        1 | self_pe3_sp_5     |
|        6 | self_pe3_sp_6     |
|     5422 | self_pe3_tr_1     |
|     1453 | self_pe3_tr_2     |
|      520 | self_pe3_tr_3     |
|     3766 | self_pe3_tr_4     |
|      368 | self_pe3_tr_5     |
|     1845 | self_pe3_tr_6     |
|       60 | self_pe3_tr_7     |
|      945 | vert_pe12_sp_1    |
|     3333 | vert_pe12_sp_2    |
|      498 | vert_pe12_sp_3    |
|     8282 | vert_pe12_sp_4    |
|      665 | vert_pe12_sp_5    |
|     3160 | vert_pe12_sp_6    |
|      123 | vert_pe12_sp_7    |
|     5599 | vert_pe12_tr_1    |
|    26035 | vert_pe12_tr_2    |
|     4127 | vert_pe12_tr_3    |
|    76690 | vert_pe12_tr_4    |
|     6548 | vert_pe12_tr_5    |
|    30028 | vert_pe12_tr_6    |
|      758 | vert_pe12_tr_7
| count(*) | biotype                |
+----------+------------------------+
|        9 | IG_C_gene              |
|       15 | IG_V_gene              |
|        2 | TR_C_gene              |
|        4 | TR_J_gene              |
|        7 | TR_V_gene              |
 count(*) | biotype         |
+----------+-----------------+
|     7433 | rnaseq_merged   |
|     4285 | rnaseq_merged_1 |
|      666 | rnaseq_merged_2 |
|      796 | rnaseq_merged_3 |
|     1109 | rnaseq_merged_4 |
|      588 | rnaseq_merged_5 |
|     1014 | rnaseq_merged_6 |
|      121 | rnaseq_merged_7 |
|     1433 | rnaseq_merged_8 |
|    41686 | rnaseq_tissue   |
|    96666 | rnaseq_tissue_1 |
|    10711 | rnaseq_tissue_2 |
|    14707 | rnaseq_tissue_3 |
|    24220 | rnaseq_tissue_4 |
|    15399 | rnaseq_tissue_5 |
|    32485 | rnaseq_tissue_6 |
|     2699 | rnaseq_tissue_7 |
|    66589 | rnaseq_tissue_8
count(*) | biotype    |
+----------+------------+
|       20 | realign_0  |
|     9016 | realign_50 |
|    24861 | realign_80 |
|    17496 | realign_95 |
