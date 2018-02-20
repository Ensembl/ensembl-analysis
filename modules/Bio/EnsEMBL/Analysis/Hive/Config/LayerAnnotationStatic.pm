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
                            ],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER2',
              BIOTYPES   => [
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
                             'human_pe12_sp_1',
                             'human_pe12_tr_1',
                             'mouse_pe12_sp_1',
                             'mouse_pe12_tr_1',
                             'mammals_pe12_sp_1',
                             'mammals_pe12_tr_1',
                            ],
              FILTER_AGAINST => ['LAYER1'],
              DISCARD    => 0,
            },


            {
              ID         => 'LAYER3',
              BIOTYPES   => [
                             'self_pe3_sp_1',
                             'self_pe3_tr_1',
                             'human_pe12_tr_2',
                             'human_pe12_sp_2',
                             'mouse_pe12_sp_2',
                             'mouse_pe12_tr_2',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2'],
              DISCARD    => 0,
            },

            {
              ID         => 'LAYER4',
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
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3'],
              DISCARD    => 0,
            },


            {
              ID         => 'LAYER5',
              BIOTYPES   => [
                             'rnaseq_merged_4',
                             'mammals_pe12_sp_3',
                             'mammals_pe12_tr_3',
                             'vert_pe12_sp_3',
                             'vert_pe12_tr_3',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER6',
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
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER7',
              BIOTYPES   => [
                             'human_pe12_tr_5',
                             'mouse_pe12_tr_5',
                             'mammals_pe12_sp_4',
                             'mammals_pe12_tr_4',
                             'vert_pe12_sp_4',
                             'vert_pe12_tr_4',
                             'realign_5',
                            ],
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6'],
              DISCARD    => 0,
            },

    ],

    fish_basic => [
             {
              ID         => 'LAYER1',
              BIOTYPES   => [
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
                             'fish_pe12_sp_3',
                             'fish_pe12_tr_3',
                             'human_pe12_sp_1',
                             'human_pe12_tr_1',
                             'mouse_pe12_sp_1',
                             'mouse_pe12_tr_1',
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
                             'mouse_pe12_sp_2',
                             'mouse_pe12_tr_2',
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
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER6',
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
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER7',
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
              FILTER_AGAINST => ['LAYER1','LAYER2','LAYER3','LAYER4','LAYER5','LAYER6'],
              DISCARD    => 0,
            },

    ],

    fish_complete => [
             {
              ID         => 'LAYER1',
              BIOTYPES   => [
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
              ID         => 'LAYER2',
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
              FILTER_AGAINST => ['LAYER1'],
              DISCARD    => 0,
            },

             {
              ID         => 'LAYER3',
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
              FILTER_AGAINST => ['LAYER1','LAYER2'],
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

