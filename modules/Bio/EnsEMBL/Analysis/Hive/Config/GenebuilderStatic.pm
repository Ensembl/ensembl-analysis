=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::Config::GenebuilderStatic

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

This is the config file for all genebuilder analysis. You should use it in your Hive configuration file to
specify the parameters of an analysis. You can either choose an existing config or you can create
a new one based on the default hash.

=head1 METHODS

  _master_config_settings: contains all possible parameters

=cut

package Bio::EnsEMBL::Analysis::Hive::Config::GenebuilderStatic;

use strict;
use warnings;

use parent ('Bio::EnsEMBL::Analysis::Hive::Config::BaseStatic');

sub _master_config {
  my ($self, $key) = @_;

  my %config = (
    default => [],
    primates_basic => [
                'IG_C_gene',
                'IG_J_gene',
		'IG_V_gene',
                'IG_D_gene',
                'TR_C_gene',
                'TR_J_gene',
                'TR_V_gene',
                'TR_D_gene',
                'realign_80',
                'realign_95',
                'realign_50',
                'self_pe12_sp_95',
                'self_pe12_sp_80',
                'self_pe12_tr_95',
                'self_pe12_tr_80',
                'rnaseq_merged_95',
                'rnaseq_merged_80',
                'rnaseq_tissue_95',
                'rnaseq_tissue_80',
                'primates_pe12_sp_80',
                'primates_pe12_tr_80',
                'primates_pe12_tr_95',
                'primates_pe12_sp_95',
                'human_pe12_sp_80',
                'human_pe12_sp_95',
                'human_pe12_tr_80',
                'human_pe12_tr_95',
                'mammals_pe12_sp_80',
                'mammals_pe12_sp_95',
                'mammals_pe12_tr_95',
                'mammals_pe12_tr_80',
    ],

    rodents_basic => [
                'IG_C_gene',
                'IG_J_gene',
		'IG_V_gene',
                'IG_D_gene',
                'TR_C_gene',
                'TR_J_gene',
                'TR_V_gene',
                'TR_D_gene',
                'realign_80',
                'realign_95',
                'realign_50',
                'self_pe12_sp_95',
                'self_pe12_sp_80',
                'self_pe12_tr_95',
                'self_pe12_tr_80',
                'self_pe3_tr_95',
                'self_pe3_tr_80',
                'rnaseq_95',
                'rnaseq_80',
                'mouse_pe12_sp_95',
                'mouse_pe12_sp_80',
                'mouse_pe12_tr_95',
                'mouse_pe12_tr_80',
                'vert_pe12_sp_95',
                'vert_pe12_sp_80',
                'vert_pe12_tr_80',
                'vert_pe12_tr_95',
                'rodents_pe12_sp_80',
                'rodents_pe12_tr_80',
                'rodents_pe12_tr_95',
                'rodents_pe12_sp_95',
                'human_pe12_sp_80',
                'human_pe12_sp_95',
                'human_pe12_tr_80',
                'human_pe12_tr_95',
                'rodents_pe3_sp_95',
                'rodents_pe3_tr_80',
                'rodents_pe3_tr_95',
                'rodents_pe3_sp_80',
                'rodents_pe45_sp_95',
                'rodents_pe45_tr_80',
                'rodents_pe45_tr_95',
                'rodents_pe45_sp_80',
                'mammals_pe12_sp_80',
                'mammals_pe12_sp_95',
                'mammals_pe12_tr_95',
                'mammals_pe12_tr_80',
    ],

    mammals_basic => [
                             'realign_95_95',
                             'realign_95_80',
                             'rnaseq_merged_95_95',
                             'rnaseq_merged_95_80',
                             'rnaseq_merged_90_80',
                             'self_pe12_sp_95_95',
                             'self_pe12_tr_95_95',
                             'self_pe12_sp_95_80',
                             'self_pe12_tr_95_80',
                             'human_pe12_sp_95_95',
                             'human_pe12_tr_95_95',
                             'human_pe12_tr_95_80',
                             'human_pe12_sp_95_80',
                             'mouse_pe12_sp_95_95',
                             'mouse_pe12_tr_95_95',
                             'mouse_pe12_sp_95_80',
                             'mouse_pe12_tr_95_80',
                             'mammals_pe12_sp_95_95',
                             'mammals_pe12_tr_95_95',
                             'mammals_pe12_sp_95_80',
                             'mammals_pe12_tr_95_80',
                             'rnaseq_tissue_95_95',
                             'rnaseq_tissue_95_80',
                             'rnaseq_tissue_90_80',
                             'human_pe12_sp_90_80',
                             'human_pe12_tr_90_80',
                             'mouse_pe12_sp_90_80',
                             'mouse_pe12_tr_90_80',
                             'self_pe3_sp_95_95',
                             'self_pe3_tr_95_95',
                             'mammals_pe12_sp_90_80',
                             'mammals_pe12_tr_90_80',
                             'vert_pe12_sp_95_95',
                             'vert_pe12_tr_95_95',
                             'rnaseq_merged_80_60',
                             'human_pe12_sp_80_60',
                             'human_pe12_tr_80_60',
                             'mouse_pe12_sp_80_60',
                             'mouse_pe12_tr_80_60',
                             'mammals_pe12_sp_80_60',
                             'mammals_pe12_tr_80_60',
                             'vert_pe12_sp_90_80',
                             'vert_pe12_tr_90_80',
                             'rnaseq_merged_70_60',
                             'rnaseq_tissue_70_60',
                             'human_pe12_sp_70_60',
                             'human_pe12_tr_70_60',
                             'mouse_pe12_sp_70_60',
                             'mouse_pe12_tr_70_60',
                             'vert_pe12_sp_80_60',
                             'vert_pe12_tr_80_60',
                             'realign_70_60',
    ],

    fish_basic => [
                             'realign_95_95',
                             'realign_95_80',
                             'rnaseq_merged_95_95',
                             'rnaseq_merged_95_80',
                             'rnaseq_merged_90_80',
                             'rnaseq_tissue_95_95',
                             'self_pe12_sp_95_95',
                             'self_pe12_tr_95_95',
                             'self_pe12_sp_95_80',
                             'self_pe12_tr_95_80',
                             'rnaseq_tissue_95_80',
                             'rnaseq_tissue_90_80',
                             'fish_pe12_sp_95_95',
                             'fish_pe12_tr_95_95',
                             'fish_pe12_sp_95_80',
                             'fish_pe12_tr_95_80',
                             'fish_pe12_sp_90_80',
                             'fish_pe12_tr_90_80',
                             'human_pe12_sp_95_95',
                             'human_pe12_tr_95_95',
                             'mouse_pe12_sp_95_95',
                             'mouse_pe12_tr_95_95',
                             'vert_pe12_sp_95_95',
                             'vert_pe12_tr_95_95',
                             'human_pe12_sp_95_80',
                             'human_pe12_tr_95_80',
                             'mouse_pe12_sp_95_80',
                             'mouse_pe12_tr_95_80',
                             'vert_pe12_sp_95_80',
                             'vert_pe12_tr_95_80',
                             'mammals_pe12_sp_95_95',
                             'mammals_pe12_tr_95_95',
                             'human_pe12_sp_90_80',
                             'human_pe12_tr_90_80',
                             'mouse_pe12_sp_90_80',
                             'mouse_pe12_tr_90_80',
                             'vert_pe12_sp_90_80',
                             'vert_pe12_tr_90_80',
                             'mammals_pe12_sp_95_80',
                             'mammals_pe12_tr_95_80',
                             'mammals_pe12_sp_90_80',
                             'mammals_pe12_tr_90_80',
                             'rnaseq_merged_80_60',
                             'fish_pe12_sp_80_60',
                             'fish_pe12_tr_80_60',
                             'human_pe12_sp_80_60',
                             'human_pe12_tr_80_60',
                             'mouse_pe12_sp_80_60',
                             'mouse_pe12_tr_80_60',
                             'rnaseq_merged_70_60',
                             'rnaseq_tissue_80_60',
                             'rnaseq_tissue_70_60',
                             'human_pe12_sp_70_60',
                             'mouse_pe12_sp_70_60',
                             'vert_pe12_sp_80_60',
                             'vert_pe12_tr_80_60',
                             'fish_pe12_sp_70_60',
                             'fish_pe12_tr_70_60',
                             'human_pe12_tr_70_60',
                             'mouse_pe12_tr_70_60',
                             'mammals_pe12_sp_80_60',
                             'mammals_pe12_tr_80_60',
                             'mammals_pe12_sp_70_60',
                             'mammals_pe12_tr_70_60',
                             'vert_pe12_sp_70_60',
                             'vert_pe12_tr_70_60',
                             'realign_70_60',
    ],

    fish_complete => [
                'realign_80',
                'realign_95',
                'realign_50',
                'self_pe12_sp_95',
                'self_pe12_sp_80',
                'self_pe12_tr_95',
                'self_pe12_tr_80',
                'self_pe3_tr_95',
                'self_pe3_tr_80',
                'rnaseq_95',
                'rnaseq_80',
                'fish_pe12_sp_80',
                'fish_pe12_sp_95',
                'fish_pe12_tr_80',
                'fish_pe12_tr_95',
                'human_pe12_sp_80',
                'human_pe12_sp_95',
                'human_pe12_tr_80',
                'human_pe12_tr_95',
                'mouse_pe12_sp_95',
                'mouse_pe12_sp_80',
                'mouse_pe12_tr_80',
                'mouse_pe12_tr_95',
                'mammals_pe12_sp_80',
                'mammals_pe12_sp_95',
                'mammals_pe12_tr_95',
                'mammals_pe12_tr_80',
                'vert_pe12_sp_95',
                'vert_pe12_sp_80',
                'vert_pe12_tr_80',
                'vert_pe12_tr_95',
    ],

    hagfish => [
                             'rnaseq_merged_95_95',
                             'rnaseq_merged_95_80',
                             'fish_pe12_sp_95_95',
                             'fish_pe12_tr_95_95',
                             'fish_pe12_sp_95_80',
                             'fish_pe12_tr_95_80',
                             'rnaseq_tissue_95_95',
                             'rnaseq_tissue_95_80',
                             'self_pe12_sp_95_95',
                             'self_pe12_tr_95_95',
                             'self_pe12_sp_95_80',
                             'self_pe12_tr_95_80',
                              'rnaseq_merged_90_70',
                              'rnaseq_merged_90_60',
                              'rnaseq_tissue_90_70',
                              'rnaseq_tissue_90_60',
                              'rnaseq_merged_80_50',
                              'rnaseq_tissue_80_50',
                              'fish_pe12_sp_90_70',
                              'fish_pe12_tr_90_70',
                              'fish_pe12_sp_90_60',
                              'fish_pe12_tr_90_60',
                             'human_pe12_sp_95_95',
                             'human_pe12_tr_95_95',
                             'human_pe12_sp_95_80',
                             'human_pe12_tr_95_80',
                             'mammals_pe12_sp_95_95',
                             'mammals_pe12_tr_95_95',
                             'mammals_pe12_sp_95_80',
                             'mammals_pe12_tr_95_80',
                             'vert_pe12_sp_95_95',
                             'vert_pe12_tr_95_95',
                             'vert_pe12_sp_95_80',
                             'vert_pe12_tr_95_80',
                            'fish_pe12_sp_80_50',
                              'fish_pe12_tr_80_50',
                              'human_pe12_sp_80_50',
                              'human_pe12_tr_80_50',
                             'mammals_pe12_sp_90_70',
                             'mammals_pe12_tr_90_70',
                             'mammals_pe12_sp_90_60',
                             'mammals_pe12_tr_90_60',
                             'vert_pe12_sp_90_70',
                             'vert_pe12_tr_90_70',
                             'vert_pe12_sp_90_60',
                             'vert_pe12_tr_90_60',
                             'vert_pe12_sp_80_50',
                             'vert_pe12_tr_80_50',
                             'mammals_pe12_sp_80_50',
                             'mammals_pe12_tr_80_50',
                              'rnaseq_merged_70_40',
                              'rnaseq_tissue_70_40',
                              'fish_pe12_sp_70_40',
                              'fish_pe12_tr_70_40',
                              'human_pe12_sp_70_40',
                              'human_pe12_tr_70_40',
                              'vert_pe12_sp_70_40',
                              'vert_pe12_tr_70_40',
                              'mammals_pe12_sp_70_40',
                              'mammals_pe12_tr_70_40',
                              'fish_pe12_sp_70_30',
                              'fish_pe12_tr_70_30',
                              'human_pe12_sp_70_30',
                              'human_pe12_tr_70_30',
                              'vert_pe12_sp_70_30',
                              'vert_pe12_tr_70_30',
                              'mammals_pe12_sp_70_30',
                              'mammals_pe12_tr_70_30',
                              'fish_pe12_sp_60_20',
                              'fish_pe12_tr_60_20',
                              'human_pe12_sp_60_20',
                              'human_pe12_tr_60_20',
                              'vert_pe12_sp_60_20',
                              'vert_pe12_tr_60_20',
                              'mammals_pe12_sp_60_20',
                              'mammals_pe12_tr_60_20',
                              'rnaseq_merged_50_25',
                              'rnaseq_tissue_50_25',
    ],

    self_patch => [
                'self_pe12_sp_95',
                'self_pe12_sp_80',
                'self_pe12_tr_95',
                'self_pe12_tr_80',
                'self_frag_pe12_sp_95',
                'self_frag_pe12_tr_95',
                'self_frag_pe12_sp_80',
                'self_frag_pe12_tr_80',
    ],
  );
  return $config{$key};
}

1;

