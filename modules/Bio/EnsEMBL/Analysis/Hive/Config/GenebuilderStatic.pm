=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

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
                             'seleno_self',
                             'self_pe12_seleno_1',
                             'self_pe12_seleno_2',
                             'self_pe12_seleno_3',
                             'self_pe12_seleno_4',
                             'cdna2genome',
                             'edited',
                             'gw_gtag',
                             'gw_nogtag',
                             'gw_exo',
                             'projection_1',
                             'projection_2',
                             'rnaseq_tissue_1',
                             'rnaseq_tissue_2',
	                     'self_pe12_sp_1',
                             'self_pe12_tr_1',
                             'self_pe12_sp_2',
                             'self_pe12_tr_2',
                             'human_pe12_sp_1',
                             'human_pe12_tr_1',
	                     'human_pe12_sp_2',
                             'human_pe12_tr_2',
	                     'primates_pe12_sp_1',
		             'primates_pe12_sp_2',
	                     'projection_3',
                             'projection_4',
	                     'projection_5',
                             'projection_6',
	                     'rnaseq_tissue_3',
	                     'rnaseq_tissue_4',
	                     'rnaseq_tissue',
	                     'genblast_rnaseq_top',
	                     'genblast_rnaseq_high',
	                     'genblast_rnaseq_medium',
                             'human_pe12_sp_6',
	                     'human_pe12_tr_6',
	                     'projection_1_noncanon',
                             'projection_2_noncanon',
                             'projection_3_noncanon',
                             'projection_4_noncanon',
                             'projection_1_pseudo',
                             'projection_2_pseudo',
                             'projection_3_pseudo',
                             'projection_4_pseudo',
    ],


    mammals_basic => [
                             'IG_C_gene',
                             'IG_J_gene',
                             'IG_V_gene',
                             'IG_D_gene',
                             'TR_C_gene',
                             'TR_J_gene',
                             'TR_V_gene',
                             'TR_D_gene',
                             'seleno_self',
                             'self_pe12_seleno_1',
                             'self_pe12_seleno_2',
                             'self_pe12_seleno_3',
                             'self_pe12_seleno_4',
	                     'cdna2genome',
                             'edited',
                             'gw_gtag',
                             'gw_nogtag',
	                     'gw_exo',
	                     'genblast_rnaseq_top',
	                     'genblast_rnaseq_high',
	                     'genblast_rnaseq_medium',
	                     'self_pe12_sp_1',
	                     'self_pe12_sp_2',
	                     'human_pe12_sp_1',
	                     'human_pe12_sp_2',
                             'human_pe12_sp_3',
                             'human_pe12_sp_4',
	                     'mouse_pe12_sp_1',
	                     'mouse_pe12_sp_2',
	                     'mouse_pe12_sp_3',
	                     'mouse_pe12_sp_4',
	                     'mammals_pe12_sp_1',
                             'mammmals_pe12_sp_2',
                             'mammals_pe12_sp_3',
                             'mammals_pe12_sp_4',
	                     'rnaseq_tissue_1',
                             'rnaseq_tissue_2',
	                     'rnaseq_tissue_3',
	                     'rnaseq_tissue_4',
	                     'rnaseq_tissue',
	                     'self_pe3_sp_1',
	                     'projection_1',
	                     'projection_2',
	                     'projection_3',
	                     'projection_4',
	                     'cdna_1',
	                     'cdna_2',
	                     'cdna_3',
	                     'cdna_4',
	                     'cdna',
	                     'projection_1_noncanon',
                             'projection_2_noncanon',
                             'projection_3_noncanon',
                             'projection_4_noncanon',
                             'projection_1_pseudo',
                             'projection_2_pseudo',
                             'projection_3_pseudo',
                             'projection_4_pseudo',
      ],

    reptiles_basic => [
                             'IG_C_gene',
                             'IG_J_gene',
                             'IG_V_gene',
                             'IG_D_gene',
                             'TR_C_gene',
                             'TR_J_gene',
                             'TR_V_gene',
                             'TR_D_gene',
                             'seleno_self',
                             'self_pe12_seleno_1',
                             'self_pe12_seleno_2',
                             'self_pe12_seleno_3',
                             'self_pe12_seleno_4',
                             'cdna2genome',
                             'edited',
                             'gw_gtag',
                             'gw_exo',
                             'rnaseq_tissue_1',
                             'rnaseq_tissue_2',
                             'rnaseq_tissue_3',
                             'rnaseq_tissue_4',
                             'self_pe12_sp_1',
                             'self_pe12_tr_1',
                             'self_pe12_sp_2',
                             'self_pe12_tr_2',
                             'projection_1',
                             'projection_2',
                             'reptiles_pe12_sp_1',
                             'reptiles_pe12_tr_1',
                             'reptiles_pe12_tr_2',
                             'reptiles_pe12_sp_2',
                             'genblast_rnaseq_top',
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
                             'projection_3',
                             'projection_4',
                             'mammals_pe12_sp_1',
                             'mammals_pe12_sp_2',
                             'self_pe3_sp_1',
                             'self_pe3_tr_1',
                             'genblast_rnaseq_medium',
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
                             'mammals_pe12_sp_3',
                             'mammals_pe12_sp_4',
                             'human_pe12_sp_5',
                             'human_pe12_sp_6',
                             'mammals_pe12_sp_5',
                             'aves_pe12_sp_5',
                             'amphibians_sp_5',
                             'reptiles_sp_5',
                             'projection_1_noncanon',
                             'projection_2_noncanon',
                             'projection_3_noncanon',
                             'projection_4_noncanon',
                             'projection_1_pseudo',
                             'projection_2_pseudo',
                             'projection_3_pseudo',
                             'projection_4_pseudo',
                             'rnaseq_tissue',
	                     'cdna_1',
                             'cdna_2',
                             'cdna_3',
                             'cdna_4',
                             'cdna',
      ],

    fish_basic => [
                             'IG_C_gene',
                             'IG_J_gene',
                             'IG_V_gene',
                             'IG_D_gene',
                             'TR_C_gene',
                             'TR_J_gene',
                             'TR_V_gene',
                             'TR_D_gene',
                             'seleno_self',
                             'self_pe12_seleno_1',
                             'self_pe12_seleno_2',
                             'self_pe12_seleno_3',
                             'self_pe12_seleno_4',
	                     'cdna2genome',
                             'edited',
                             'gw_gtag',
    	                     'gw_exo',
	                     'genblast_rnaseq_top',
	                     'genblast_rnaseq_high',
	                     'genblast_rnaseq_medium',
	                     'rnaseq_tissue_1',
                             'rnaseq_tissue_2',
                             'rnaseq_tissue_3',
                             'rnaseq_tissue_4',
	                     'rnaseq_tissue_5',
	                     'rnaseq_tissue_6',
	                     'rnaseq_tissue',
	                     'self_pe12_sp_1',
                             'self_pe12_tr_1',
                             'self_pe12_sp_2',
                             'self_pe12_tr_2',
	                     'self_pe12_sp_3',
	                     'self_pe12_tr_3',
	                     'self_pe12_sp_4',
	                     'self_pe12_tr_4',
	                     'self_pe12_sp_5',
	                     'self_pe12_tr_5',
                             'self_pe12_sp_6',
               	       	     'self_pe12_tr_6',	
	                     'self_pe3_sp_1',
	                     'self_pe3_sp_2',
	                     'self_pe3_sp_3',
	                     'self_pe3_sp_4',
	                     'fish_pe12_sp_1',
                             'fish_pe12_tr_1',
                             'fish_pe12_sp_2',
                             'fish_pe12_tr_2',
                             'fish_pe12_sp_3',
                             'fish_pe12_tr_3',
                             'fish_pe12_sp_4',
                             'fish_pe12_tr_4',
                             'fish_pe12_sp_5',
                             'human_pe12_sp_1',
                             'human_pe12_tr_1',
                             'human_pe12_sp_2',
                             'human_pe12_tr_2',
                             'human_pe12_sp_3',
                             'human_pe12_tr_3',
                             'human_pe12_sp_4',
                             'human_pe12_tr_4',
                             'human_pe12_sp_5',
	                     'mammals_pe12_sp_1',
                             'mammals_pe12_tr_1',
                             'mammals_pe12_sp_2',
                             'mammals_pe12_tr_2',
                             'mammals_pe12_sp_3',
                             'mammals_pe12_sp_4',
	                     'vert_pe12_sp_1',
                             'vert_pe12_tr_1',
                             'vert_pe12_sp_2',
                             'vert_pe12_tr_2',
                             'vert_pe12_sp_3',
	                     'vert_pe12_sp_4',
	                     'cdna_1',
                             'cdna_2',
                             'cdna_3',
                             'cdna_4',
	                     'cdna_5',
	                     'cdna_6',
                             'cdna',
    ],

    distant_vertebrate => [
                             'IG_C_gene',
                             'IG_J_gene',
                             'IG_V_gene',
                             'IG_D_gene',
                             'TR_C_gene',
                             'TR_J_gene',
                             'TR_V_gene',
                             'TR_D_gene',
                             'seleno_self',
                             'self_pe12_seleno_1',
                             'self_pe12_seleno_2',
                             'self_pe12_seleno_3',
                             'self_pe12_seleno_4',
                             'cdna2genome',
                             'edited',
                             'gw_gtag',
                             'gw_nogtag',
                             'gw_exo',
                             'rnaseq_tissue_1',
                             'rnaseq_tissue_2',
                             'rnaseq_tissue_3',
                             'vert_pe12_sp_1',
                             'vert_pe12_sp_2',
                             'self_pe12_sp_1',
                             'self_pe12_sp_2',
                             'self_pe12_tr_1',
                             'self_pe12_tr_2',
	                     'rnaseq_tissue_4',
               	       	     'rnaseq_tissue',
                             'vert_pe12_sp_3',
                             'vert_pe12_tr_1',
                             'vert_pe12_tr_2',
                             'vert_pe12_sp_4',
                             'vert_pe12_tr_3',
                             'vert_pe12_sp_5',
                             'vert_pe12_sp_6',
                             'cdna_1',
                             'cdna_2',
                             'cdna_3',
                             'cdna_4',
                             'cdna',
    ],

    aves_basic => [
                             'IG_C_gene',
                             'IG_J_gene',
                             'IG_V_gene',
                             'IG_D_gene',
                             'TR_C_gene',
                             'TR_J_gene',
                             'TR_V_gene',
                             'TR_D_gene',
                             'seleno_self',
	                     'cdna2genome',
                             'edited',
                             'gw_gtag',
                             'gw_exo',
                             'cdna_1',
                             'cdna_2',
                             'cdna_3',
                             'cdna_4',
	                     'rnaseq_tissue_1',
                             'rnaseq_tissue_2',
                             'rnaseq_tissue_3',
                             'rnaseq_tissue_4',
                             'self_pe12_sp_1',
                             'self_pe12_tr_1',
                             'self_pe12_sp_2',
                             'self_pe12_tr_2',
                             'projection_1',
                             'projection_2',
	                     'aves_pe12_sp_1',
                             'aves_pe12_tr_1',
                             'aves_pe12_tr_2',
                             'aves_pe12_sp_2',
                             'genblast_rnaseq_top',
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
                             'human_pe12_sp_3',
                             'human_pe12_tr_3',
                             'human_pe12_sp_4',
                             'human_pe12_tr_4',
                             'mammals_pe12_sp_3',
                             'mammals_pe12_tr_3',
                             'mammals_pe12_sp_4',
                             'mammals_pe12_tr_4',
                             'human_pe12_sp_5',
                             'human_pe12_sp_6',
                             'mammals_pe12_sp_5',
                             'aves_pe12_sp_5',
                             'reptiles_sp_5',
                             'projection_1_noncanon',
                             'projection_2_noncanon',
                             'projection_3_noncanon',
                             'projection_4_noncanon',
                             'projection_1_pseudo',
                             'projection_2_pseudo',
                             'projection_3_pseudo',
                             'projection_4_pseudo',
	                     'cdna',
	                     'rnaseq_tissue',
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
