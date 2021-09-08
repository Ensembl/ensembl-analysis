=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Analysis::Hive::Config::TranscriptSelection;

use strict;
use warnings;
use File::Spec::Functions;

use Bio::EnsEMBL::ApiVersion qw/software_version/;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(get_analysis_settings);
use base ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

sub default_options {
  my ($self) = @_;
  return {
    # inherit other stuff from the base class
    %{ $self->SUPER::default_options() },

######################################################
#
# Variable settings- You change these!!!
#
######################################################
########################
# Misc setup info
########################
    'dbowner' => '' || $ENV{EHIVE_USER} || $ENV{USER},
    'pipeline_name' => '' || $self->o('production_name') . '_' . $self->o('ensembl_release'),
    'user_r'                    => '',                                                          # read only db user
    'user'                      => '',                                                          # write db user
    'password'                  => '',                                                          # password for write db user
    'server_set'                => '',                                                          # What server set to user, e.g. set1
    'pipe_db_host'              => '',                                                          # host for pipe db
    'databases_host'            => '',                                                          # host for general output dbs
    'dna_db_host'               => '',                                                          # host for dna db
    'pipe_db_port'              => '',                                                          # port for pipeline host
    'databases_port'            => '',                                                          # port for general output db host
    'dna_db_port'               => '',                                                          # port for dna db host
    'release_number'            => '' || $self->o('ensembl_release'),
    'species_name'              => '',                                                          # e.g. mus_musculus
    'production_name'           => '',                                                          # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
    'uniprot_set'               => '',                                                          # e.g. mammals_basic, check UniProtCladeDownloadStatic.pm module in hive config dir for suitable set,
    'output_path'               => '',                                                          # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses
    'use_genome_flatfile'       => '1',                                                         # This will read sequence where possible from a dumped flatfile instead of the core db
    'skip_projection'           => '0',                                                         # Will skip projection process if 1
    'skip_rnaseq'               => '0',                                                         # Will skip rnaseq analyses if 1
    'skip_cleaning'             => '0',                                                         # Will skip the cleaning phase, will keep more genes/transcripts but some lower quality models may be kept

########################
# Pipe and ref db info
########################

    'cdna_db_host'   => $self->o('databases_host'),
    'cdna_db_port'   => $self->o('databases_port'),
    'cdna_db_user'   => $self->o('user'),
    'cdna_db_pass'   => $self->o('password'),

    'genblast_nr_db_host'   => $self->o('databases_host'),
    'genblast_nr_db_port'   => $self->o('databases_port'),
    'genblast_nr_db_user'   => $self->o('user'),
    'genblast_nr_db_pass'   => $self->o('password'),

    'genblast_rnaseq_support_nr_db_host'   => $self->o('databases_host'),
    'genblast_rnaseq_support_nr_db_port'   => $self->o('databases_port'),
    'genblast_rnaseq_support_nr_db_user'   => $self->o('user'),
    'genblast_rnaseq_support_nr_db_pass'   => $self->o('password'),

    'ig_tr_db_host'   => $self->o('databases_host'),
    'ig_tr_db_port'   => $self->o('databases_port'),
    'ig_tr_db_user'   => $self->o('user'),
    'ig_tr_db_pass'   => $self->o('password'),

    'best_targeted_db_host'   => $self->o('databases_host'),
    'best_targeted_db_port'   => $self->o('databases_port'),
    'best_targeted_db_user'   => $self->o('user'),
    'best_targeted_db_pass'   => $self->o('password'),

    'selected_projection_db_host'   => $self->o('databases_host'),
    'selected_projection_db_port'   => $self->o('databases_port'),
    'selected_projection_db_user'   => $self->o('user'),
    'selected_projection_db_pass'   => $self->o('password'),

    'long_read_final_db_host'   => $self->o('databases_host'),
    'long_read_final_db_port'   => $self->o('databases_port'),
    'long_read_final_db_user'   => $self->o('user'),
    'long_read_final_db_pass'   => $self->o('password'),

    'rnaseq_for_layer_db_host'   => $self->o('databases_host'),
    'rnaseq_for_layer_db_port'   => $self->o('databases_port'),
    'rnaseq_for_layer_db_user'   => $self->o('user'),
    'rnaseq_for_layer_db_pass'   => $self->o('password'),

    'rnaseq_for_layer_nr_db_host'   => $self->o('databases_host'),
    'rnaseq_for_layer_nr_db_port'   => $self->o('databases_port'),
    'rnaseq_for_layer_nr_db_user'   => $self->o('user'),
    'rnaseq_for_layer_nr_db_pass'   => $self->o('password'),

    # Layering is one of the most intesnive steps, so separating it off the main output server helps
    # Have also set module to use flatfile seq retrieval, so even if it's on the same server as the
    # core, the core should not be accessed
    'layering_db_host'   => $self->o('dna_db_host'),
    'layering_db_port'   => $self->o('dna_db_port'),
    'layering_db_user'   => $self->o('user'),
    'layering_db_pass'   => $self->o('password'),

    'utr_db_host'   => $self->o('databases_host'),
    'utr_db_port'   => $self->o('databases_port'),
    'utr_db_user'   => $self->o('user'),
    'utr_db_pass'   => $self->o('password'),

    'genebuilder_db_host'   => $self->o('databases_host'),
    'genebuilder_db_port'   => $self->o('databases_port'),
    'genebuilder_db_user'   => $self->o('user'),
    'genebuilder_db_pass'   => $self->o('password'),

    'pseudogene_db_host'   => $self->o('databases_host'),
    'pseudogene_db_port'   => $self->o('databases_port'),
    'pseudogene_db_user'   => $self->o('user'),
    'pseudogene_db_pass'   => $self->o('password'),

    'ncrna_db_host'   => $self->o('databases_host'),
    'ncrna_db_port'   => $self->o('databases_port'),
    ncrna_db_name     => $self->o('dbowner') . '_' . $self->o('production_name') . '_ncrna_' . $self->o('release_number'),
    'ncrna_db_user'   => $self->o('user'),
    'ncrna_db_pass'   => $self->o('password'),

    'final_geneset_db_host'   => $self->o('databases_host'),
    'final_geneset_db_port'   => $self->o('databases_port'),
    'final_geneset_db_user'   => $self->o('user'),
    'final_geneset_db_pass'   => $self->o('password'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'ensembl_release'      => $ENV{ENSEMBL_RELEASE},     # this is the current release version on staging to be able to get the correct database
    'production_db_host'   => 'mysql-ens-meta-prod-1',
    'production_db_port'   => '4483',

    databases_to_delete => ['layering_db', 'utr_db', 'genebuilder_db', 'pseudogene_db', 'final_geneset_db'],    #, 'projection_realign_db'

######################################################
#
# Mostly constant settings
#
######################################################

    genome_dumps => catdir( $self->o('output_path'), 'genome_dumps' ),
    # This one is used in replacement of the dna table in the core db, so where analyses override slice->seq. Has simple headers with just the seq_region name. Also used by bwa in the RNA-seq analyses. Not masked
    faidx_genome_file => catfile( $self->o('genome_dumps'), $self->o('species_name') . '_toplevel.fa' ),

    create_toplevel_dbs => [
      $self->o('genblast_nr_db'),
      $self->o('selected_projection_db'),
      $self->o('rnaseq_for_layer_db'),
    ],
    layering_input_gene_dbs => [
      $self->o('genblast_nr_db'),
      $self->o('genblast_rnaseq_support_nr_db'),
      $self->o('rnaseq_for_layer_nr_db'),
      $self->o('selected_projection_db'),
      $self->o('ig_tr_db'),
      $self->o('best_targeted_db'),
      $self->o('long_read_final_db'),
    ],

    split_intergenic_dbs => [
      $self->o('genblast_nr_db'),
      $self->o('genblast_rnaseq_support_nr_db'),
      $self->o('rnaseq_for_layer_nr_db'),
      $self->o('selected_projection_db'),
      $self->o('ig_tr_db'),
      $self->o('best_targeted_db'),
      $self->o('long_read_final_db'),
      $self->o('cdna_db'),
    ],

    utr_donor_dbs => [
      $self->o('cdna_db'),
      $self->o('rnaseq_for_layer_db'),
      $self->o('long_read_final_db'),
    ],

    utr_acceptor_dbs => [
      $self->o('layering_db'),
    ],

    'utr_biotype_priorities' => {
      'rnaseq' => 2,
      'cdna'   => 1,
    },

    'cleaning_blessed_biotypes' => {
      'pseudogene'           => 1,
      'processed_pseudogene' => 1,
      'IG_C_gene'            => 1,
      'IG_V_gene'            => 1,
      'TR_C_gene'            => 1,
      'TR_D_gene'            => 1,
      'TR_V_gene'            => 1,
      'lncRNA'               => 1,
    },

    ensembl_analysis_script => catdir( $self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts' ),
    remove_duplicates_script_path => catfile( $self->o('ensembl_analysis_script'), 'find_and_remove_duplicates.pl' ),
    remove_small_orf_script => catfile( $self->o('ensembl_analysis_script'), 'genebuild', 'remove_small_orf.pl' ),

    # cutoffs for removing small_orf genes
    'small_orf_cutoff' => '100',
    'intron_cutoff'    => '75',


########################
# Executable paths
########################

    'blast_type' => 'ncbi',    # It can be 'ncbi', 'wu', or 'legacy_ncbi'

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# No option below this mark should be modified
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#

    cdna_db_name => $self->o('dbowner') . '_' . $self->o('production_name') . '_cdna_' . $self->o('release_number'),
    genblast_nr_db_name => $self->o('dbowner').'_'.$self->o('production_name').'_genblast_nr_'.$self->o('release_number'),
    genblast_rnaseq_support_nr_db_name => $self->o('dbowner').'_'.$self->o('production_name').'_gb_rnaseq_nr_'.$self->o('release_number'),
    ig_tr_db_name => $self->o('dbowner').'_'.$self->o('production_name').'_igtr_'.$self->o('release_number'),
    best_targeted_db_name => $self->o('dbowner').'_'.$self->o('production_name').'_bt_'.$self->o('release_number'),
    long_read_final_db_name => $self->o('dbowner') . '_' . $self->o('production_name') . '_lrfinal_' . $self->o('release_number'),
    selected_projection_db_name => $self->o('dbowner') . '_' . $self->o('production_name') . '_sel_proj_' . $self->o('release_number'),
    rnaseq_for_layer_db_name => $self->o('dbowner') . '_' . $self->o('production_name') . '_rnalayer_nr_' . $self->o('release_number'),
    rnaseq_for_layer_nr_db_name => $self->o('dbowner').'_'.$self->o('production_name').'_rnalayer_nr_'.$self->o('release_number'),
    layering_db_name => $self->o('dbowner') . '_' . $self->o('production_name') . '_layer_' . $self->o('release_number'),
    utr_db_name => $self->o('dbowner') . '_' . $self->o('production_name') . '_utr_' . $self->o('release_number'),
    genebuilder_db_name => $self->o('dbowner') . '_' . $self->o('production_name') . '_gbuild_' . $self->o('release_number'),
    pseudogene_db_name => $self->o('dbowner') . '_' . $self->o('production_name') . '_pseudo_' . $self->o('release_number'),
    final_geneset_db_name => $self->o('dbowner') . '_' . $self->o('production_name') . '_final_' . $self->o('release_number'),

########################
# db info
########################

    'cdna_db' => {
      -dbname => $self->o('cdna_db_name'),
      -host   => $self->o('cdna_db_host'),
      -port   => $self->o('cdna_db_port'),
      -user   => $self->o('cdna_db_user'),
      -pass   => $self->o('cdna_db_pass'),
      -driver => $self->o('hive_driver'),
    },

    'genblast_nr_db' => {
      -dbname => $self->o('genblast_nr_db_name'),
      -host   => $self->o('genblast_nr_db_host'),
      -port   => $self->o('genblast_nr_db_port'),
      -user   => $self->o('genblast_nr_db_user'),
      -pass   => $self->o('genblast_nr_db_pass'),
      -driver => $self->o('hive_driver'),
    },

    'genblast_rnaseq_support_nr_db' => {
      -dbname => $self->o('genblast_rnaseq_support_nr_db_name'),
      -host   => $self->o('genblast_rnaseq_support_nr_db_host'),
      -port   => $self->o('genblast_rnaseq_support_nr_db_port'),
      -user   => $self->o('genblast_rnaseq_support_nr_db_user'),
      -pass   => $self->o('genblast_rnaseq_support_nr_db_pass'),
      -driver => $self->o('hive_driver'),
    },

    'ig_tr_db' => {
      -dbname => $self->o('ig_tr_db_name'),
      -host   => $self->o('ig_tr_db_host'),
      -port   => $self->o('ig_tr_db_port'),
      -user   => $self->o('ig_tr_db_user'),
      -pass   => $self->o('ig_tr_db_pass'),
      -driver => $self->o('hive_driver'),
    },

    'best_targeted_db' => {
      -dbname => $self->o('best_targeted_db_name'),
      -host   => $self->o('best_targeted_db_host'),
      -port   => $self->o('best_targeted_db_port'),
      -user   => $self->o('best_targeted_db_user'),
      -pass   => $self->o('best_targeted_db_pass'),
      -driver => $self->o('hive_driver'),
    },

    long_read_final_db => {
      -dbname => $self->o('long_read_final_db_name'),
      -host   => $self->o('long_read_final_db_host'),
      -port   => $self->o('long_read_final_db_port'),
      -user   => $self->o('long_read_final_db_user'),
      -pass   => $self->o('long_read_final_db_pass'),
      -driver => $self->o('hive_driver'),
    },

    'selected_projection_db' => {
      -dbname => $self->o('selected_projection_db_name'),
      -host   => $self->o('selected_projection_db_host'),
      -port   => $self->o('selected_projection_db_port'),
      -user   => $self->o('selected_projection_db_user'),
      -pass   => $self->o('selected_projection_db_pass'),
      -driver => $self->o('hive_driver'),
    },

    'rnaseq_for_layer_db' => {
      -dbname => $self->o('rnaseq_for_layer_db_name'),
      -host   => $self->o('rnaseq_for_layer_db_host'),
      -port   => $self->o('rnaseq_for_layer_db_port'),
      -user   => $self->o('rnaseq_for_layer_db_user'),
      -pass   => $self->o('rnaseq_for_layer_db_pass'),
      -driver => $self->o('hive_driver'),
    },

    'rnaseq_for_layer_nr_db' => {
      -dbname => $self->o('rnaseq_for_layer_nr_db_name'),
      -host   => $self->o('rnaseq_for_layer_nr_db_host'),
      -port   => $self->o('rnaseq_for_layer_nr_db_port'),
      -user   => $self->o('rnaseq_for_layer_nr_db_user'),
      -pass   => $self->o('rnaseq_for_layer_nr_db_pass'),
      -driver => $self->o('hive_driver'),
    },

    'layering_db' => {
      -dbname => $self->o('layering_db_name'),
      -host   => $self->o('layering_db_host'),
      -port   => $self->o('layering_db_port'),
      -user   => $self->o('layering_db_user'),
      -pass   => $self->o('layering_db_pass'),
      -driver => $self->o('hive_driver'),
    },

    'utr_db' => {
      -dbname => $self->o('utr_db_name'),
      -host   => $self->o('utr_db_host'),
      -port   => $self->o('utr_db_port'),
      -user   => $self->o('utr_db_user'),
      -pass   => $self->o('utr_db_pass'),
      -driver => $self->o('hive_driver'),
    },

    'genebuilder_db' => {
      -dbname => $self->o('genebuilder_db_name'),
      -host   => $self->o('genebuilder_db_host'),
      -port   => $self->o('genebuilder_db_port'),
      -user   => $self->o('genebuilder_db_user'),
      -pass   => $self->o('genebuilder_db_pass'),
      -driver => $self->o('hive_driver'),
    },

    'pseudogene_db' => {
      -dbname => $self->o('pseudogene_db_name'),
      -host   => $self->o('pseudogene_db_host'),
      -port   => $self->o('pseudogene_db_port'),
      -user   => $self->o('pseudogene_db_user'),
      -pass   => $self->o('pseudogene_db_pass'),
      -driver => $self->o('hive_driver'),
    },

    'ncrna_db' => {
      -dbname => $self->o('ncrna_db_name'),
      -host   => $self->o('ncrna_db_host'),
      -port   => $self->o('ncrna_db_port'),
      -user   => $self->o('ncrna_db_user'),
      -pass   => $self->o('ncrna_db_pass'),
      -driver => $self->o('hive_driver'),
    },

    'final_geneset_db' => {
      -dbname => $self->o('final_geneset_db_name'),
      -host   => $self->o('final_geneset_db_host'),
      -port   => $self->o('final_geneset_db_port'),
      -user   => $self->o('final_geneset_db_user'),
      -pass   => $self->o('final_geneset_db_pass'),
      -driver => $self->o('hive_driver'),
    },

  };
}

sub pipeline_create_commands {
  my ($self) = @_;
  return [
    # inheriting database and hive tables' creation
    @{ $self->SUPER::pipeline_create_commands },
  ];
}

sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{ $self->SUPER::pipeline_wide_parameters },
    skip_projection     => $self->o('skip_projection'),
    skip_rnaseq         => $self->o('skip_rnaseq'),
    use_genome_flatfile => $self->o('use_genome_flatfile'),
    genome_file         => $self->o('faidx_genome_file'),
    }
}

## See diagram for pipeline structure
sub pipeline_analyses {
  my ($self) = @_;

  return [

############################################################################
#
# Finalisation analyses
#
############################################################################
    {
      -logic_name => 'create_layering_output_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('dna_db'),
        target_db   => $self->o('layering_db'),
        create_type => 'clone',
      },
      -rc_name   => 'default',
      -input_ids  => [{}],
      -flow_into => {
        1 => ['create_utr_db'],
      },
    },

    {
      -logic_name => 'create_utr_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('dna_db'),
        target_db   => $self->o('utr_db'),
        create_type => 'clone',
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['create_genebuilder_db'],
      },
    },

    {
      -logic_name => 'create_genebuilder_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('dna_db'),
        target_db   => $self->o('genebuilder_db'),
        create_type => 'clone',
      },
      -rc_name   => 'default',
      -flow_into => {
        '1->A' => ['create_toplevel_slices'],
        'A->1' => ['layer_annotation_sanity_checks'],
      },
    },

    {
      -logic_name => 'create_toplevel_slices',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db             => $self->o('dna_db'),
        iid_type              => 'slice',
        coord_system_name     => 'toplevel',
        include_non_reference => 0,
        top_level             => 1,
        # These options will create only slices that have a gene on the slice in one of the feature dbs
        feature_constraint => 1,
        feature_type       => 'gene',
        feature_dbs        => $self->o('create_toplevel_dbs'),
      },
      -flow_into => {
        '2' => ['split_slices_on_intergenic'],
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'split_slices_on_intergenic',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveFindIntergenicRegions',
      -parameters => {
        dna_db         => $self->o('dna_db'),
        input_gene_dbs => $self->o('split_intergenic_dbs'),
        iid_type       => 'slice',
      },
      -batch_size    => 100,
      -hive_capacity => $self->hive_capacity_classes->{'hc_medium'},
      -rc_name       => '5GB',
      -flow_into     => {
        '2' => ['layer_annotation'],
      },
    },

    {
      -logic_name => 'layer_annotation',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLayerAnnotation',
      -parameters => {
        use_genome_flatfile => 1,
        genome_file         => $self->o('faidx_genome_file'),
        logic_name          => 'layer_annotation',
        module              => 'HiveLayerAnnotation',
        TARGETDB_REF        => $self->o('layering_db'),
        SOURCEDB_REFS       => $self->o('layering_input_gene_dbs'),
        # Filtering is using done at the exon-overlap level
        # When no FILTER exists in this file, this is the default behaviour
        # If you would like to filter in a different way, please specify filter
        #FILTER => 'Bio::EnsEMBL::Analysis::Tools::GenomeOverlapFilter',
        #FILTER => 'Bio::EnsEMBL::Analysis::Tools::AllExonOverlapFilter',
        FILTER => 'Bio::EnsEMBL::Analysis::Tools::CodingExonOverlapFilter',
        # ordered list of annotation layers. Genes from lower layers
        # are only retained if they do not "interfere" with genes from
        # higher layers. Genes in "Discard" layers are when assessing
        # interference, but are not written to the final database
        LAYERS => get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::LayerAnnotationStatic', $self->o('uniprot_set'), undef, 'ARRAY' ),
      },
      -rc_name   => '4GB',
      -flow_into => {
        '1->A' => ['run_utr_addition'],
        'A->1' => ['genebuilder'],
      },
      -batch_size => 50,
    },

    {
      -logic_name => 'run_utr_addition',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveUTRAddition',
      -parameters => {
        logic_name             => 'utr_addition',
        dna_db                 => $self->o('dna_db'),
        donor_dbs              => $self->o('utr_donor_dbs'),
        acceptor_dbs           => $self->o('utr_acceptor_dbs'),
        utr_biotype_priorities => $self->o('utr_biotype_priorities'),
        target_db              => $self->o('utr_db'),
        iid_type               => 'slice',
      },
      -batch_size    => 20,
      -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
      -rc_name       => '5GB',
      -flow_into     => {
        -1 => ['run_utr_addition_10GB'],
      },
    },

    {
      -logic_name => 'run_utr_addition_10GB',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveUTRAddition',
      -parameters => {
        logic_name             => 'utr_addition',
        dna_db                 => $self->o('dna_db'),
        donor_dbs              => $self->o('utr_donor_dbs'),
        acceptor_dbs           => $self->o('utr_acceptor_dbs'),
        utr_biotype_priorities => $self->o('utr_biotype_priorities'),
        target_db              => $self->o('utr_db'),
        iid_type               => 'slice',
      },
      -batch_size    => 20,
      -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
      -rc_name       => '10GB',
      -flow_into     => {
        -1 => ['run_utr_addition_30GB'],
      },
    },

    {
      -logic_name => 'run_utr_addition_30GB',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveUTRAddition',
      -parameters => {
        logic_name             => 'utr_addition',
        dna_db                 => $self->o('dna_db'),
        donor_dbs              => $self->o('utr_donor_dbs'),
        acceptor_dbs           => $self->o('utr_acceptor_dbs'),
        utr_biotype_priorities => $self->o('utr_biotype_priorities'),
        target_db              => $self->o('utr_db'),
        iid_type               => 'slice',
      },
      -batch_size    => 20,
      -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
      -rc_name       => '30GB',
      -flow_into     => {
        -1 => ['utr_memory_failover'],
      },
    },

    {
      -logic_name => 'utr_memory_failover',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveUTRAddition',
      -parameters => {
        logic_name             => 'utr_addition',
        dna_db                 => $self->o('dna_db'),
        donor_dbs              => $self->o('utr_donor_dbs'),
        acceptor_dbs           => $self->o('utr_acceptor_dbs'),
        utr_biotype_priorities => $self->o('utr_biotype_priorities'),
        target_db              => $self->o('utr_db'),
        iid_type               => 'slice',
        copy_only              => 1,
      },
      -batch_size    => 20,
      -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
      -rc_name       => '10GB',
    },

    {
      -logic_name => 'genebuilder',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGeneBuilder',
      -parameters => {
        source_db   => $self->o('utr_db'),
        target_db   => $self->o('genebuilder_db'),
        dna_db      => $self->o('dna_db'),
        logic_name  => 'ensembl',
        module      => 'HiveGeneBuilder',
        INPUT_GENES => {
          source_db => get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::GenebuilderStatic',
            $self->o('uniprot_set'), undef, 'ARRAY' ),
        },
        OUTPUT_BIOTYPE              => 'protein_coding',
        MAX_TRANSCRIPTS_PER_CLUSTER => 10,
        MIN_SHORT_INTRON_LEN        => 7,                  #introns shorter than this seem
                                                           #to be real frame shifts and shoudn't be ignored
        MAX_SHORT_INTRON_LEN        => 15,
        BLESSED_BIOTYPES            => {
          'ccds_gene' => 1,
          'IG_C_gene' => 1,
          'IG_J_gene' => 1,
          'IG_V_gene' => 1,
          'IG_D_gene' => 1,
          'TR_C_gene' => 1,
          'TR_J_gene' => 1,
          'TR_V_gene' => 1,
          'TR_D_gene' => 1,
        },
        #the biotypes of the best genes always to be kept
        MAX_EXON_LENGTH => 20000,
        #if the coding_only flag is set to 1, the transcript clustering into genes is done over coding exons only
        # the current standard way is to cluster only on coding exons
        CODING_ONLY => 1,
      },
      -rc_name       => '4GB',
      -hive_capacity => $self->hive_capacity_classes->{'hc_high'},
    },

    {
      -logic_name => 'layer_annotation_sanity_checks',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
      -parameters => {
        skip_rnaseq                => $self->o('skip_rnaseq'),
        skip_projection            => $self->o('skip_projection'),
        target_db                  => $self->o('layering_db'),
        sanity_check_type          => 'gene_db_checks',
        min_allowed_feature_counts => get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
          'gene_db_checks' )->{ $self->o('uniprot_set') }->{'layer'},
      },
      -rc_name   => '4GB',
      -flow_into => {
        1 => ['genebuilder_sanity_checks'],
      },
    },

    {
      -logic_name => 'genebuilder_sanity_checks',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
      -parameters => {
        skip_rnaseq                => $self->o('skip_rnaseq'),
        skip_projection            => $self->o('skip_projection'),
        target_db                  => $self->o('genebuilder_db'),
        sanity_check_type          => 'gene_db_checks',
        min_allowed_feature_counts => get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
          'gene_db_checks' )->{ $self->o('uniprot_set') }->{'genebuilder'},
      },
      -rc_name   => '4GB',
      -flow_into => {
        1 => ['restore_ig_tr_biotypes'],
      },
    },

    {
      -logic_name => 'restore_ig_tr_biotypes',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('genebuilder_db'),
        sql     => [
          'UPDATE gene JOIN transcript USING(gene_id) SET gene.biotype = transcript.biotype' .
            ' WHERE transcript.biotype LIKE "IG\_%" OR transcript.biotype LIKE "TR\_%"',
        ],
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['create_pseudogene_db'],
      },
    },

    {
      -logic_name => 'create_pseudogene_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('dna_db'),
        target_db   => $self->o('pseudogene_db'),
        create_type => 'clone',
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['update_lncrna_biotypes'],
      },
    },

    {
      -logic_name => 'update_lncrna_biotypes',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('genebuilder_db'),
        sql     => [
          'UPDATE transcript SET biotype="pre_lncRNA" WHERE biotype IN ("rnaseq_merged","rnaseq_tissue","cdna")',
          'UPDATE gene JOIN transcript USING(gene_id) SET gene.biotype="pre_lncRNA" WHERE transcript.biotype="pre_lncRNA"',
        ],
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['pseudogenes'],
      },
    },

    {
      -logic_name => 'pseudogenes',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivePseudogenes',
      -parameters => {
        single_multi_file => 1,
        output_path       => $self->o('output_path') . '/pseudogenes/',
        input_gene_db     => $self->o('genebuilder_db'),
        repeat_db         => $self->o('dna_db'),
        output_db         => $self->o('pseudogene_db'),
        dna_db            => $self->o('dna_db'),
        logic_name        => 'pseudogenes',
        module            => 'HivePseudogenes',
        %{ get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::PseudoGeneStatic', 'pseudogenes' ) },
      },
      -rc_name   => '30GB',
      -flow_into => {
        1 => ['remove_small_orf'],
      },
    },

    {
      -logic_name => 'remove_small_orf',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . $self->o('remove_small_orf_script') .
          ' -host ' . $self->o( 'pseudogene_db', '-host' ) .
          ' -port ' . $self->o( 'pseudogene_db', '-port' ) .
          ' -user_w ' . $self->o( 'pseudogene_db', '-user' ) .
          ' -pass ' . $self->o( 'pseudogene_db', '-pass' ) .
          ' -dbname ' . $self->o( 'pseudogene_db', '-dbname' ) .
          ' -dna_host ' . $self->o( 'dna_db', '-host' ) .
          ' -dna_port ' . $self->o( 'dna_db', '-port' ) .
          ' -user_r ' . $self->o( 'dna_db', '-user' ) .
          ' -dna_dbname ' . $self->o( 'dna_db', '-dbname' ) .
          ' -orf_cutoff ' . $self->o('small_orf_cutoff') .
          ' -intron_cutoff ' . $self->o('intron_cutoff'),
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['format_blast_db'],
      },
    },

    {
      -logic_name => 'format_blast_db',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'if [ "' . $self->o('blast_type') .
          '" = "ncbi" ];then makeblastdb -dbtype nucl -in ' .
          $self->o('output_path') . '/pseudogenes/all_multi_exon_genes.fasta;' .
          ' else xdformat -n ' . $self->o('output_path') . '/pseudogenes/all_multi_exon_genes.fasta;fi'
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['spliced_elsewhere'],
      },
    },

    {
      -logic_name => 'spliced_elsewhere',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSplicedElsewhere',
      -parameters => {
        multi_exon_db_path => $self->o('output_path') . '/pseudogenes/',
        input_gene_db      => $self->o('genebuilder_db'),
        repeat_db          => $self->o('dna_db'),
        output_db          => $self->o('pseudogene_db'),
        dna_db             => $self->o('dna_db'),
        logic_name         => 'spliced_elsewhere',
        module             => 'HiveSplicedElsewhere',
        %{ get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::PseudoGeneStatic', 'pseudogenes' ) },
      },
      -rc_name   => '5GB',
      -flow_into => {
        1 => ['create_final_geneset_db'],
      },
    },

    {
      -logic_name => 'create_final_geneset_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('pseudogene_db'),
        target_db   => $self->o('final_geneset_db'),
        create_type => 'copy',
      },
      -rc_name   => 'default',
      -flow_into => {
        '1' => ['filter_lncrnas'],
      },
    },

    {
      -logic_name => 'filter_lncrnas',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveFilterlncRNAs',
      -parameters => {
        input_gene_db => $self->o('final_geneset_db'),
        dna_db        => $self->o('dna_db'),
        logic_name    => 'filter_lncrnas',
        module        => 'HiveFilterlncRNAs',
      },
      -rc_name   => '3GB',
      -flow_into => {
        1 => ['change_biotype_for_weak_cds'],
      },
    },

    {
      -logic_name => 'change_biotype_for_weak_cds',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('final_geneset_db'),
        sql     => [
          'UPDATE transcript JOIN transcript_supporting_feature USING(transcript_id) JOIN protein_align_feature ON feature_id=protein_align_feature_id SET biotype="low_coverage" WHERE feature_type="protein_align_feature" AND hcoverage < 50 AND biotype not like "pseudo%"',
          'UPDATE gene JOIN transcript USING(gene_id) SET gene.biotype="low_coverage" WHERE transcript.biotype="low_coverage" AND gene_id NOT IN ('
            .'SELECT DISTINCT(gbio.gene_id) FROM ('
              .'SELECT count(*) AS biotype_cnt,g.gene_id FROM gene g,transcript t WHERE g.gene_id=t.gene_id AND g.gene_id IN ('
                .'SELECT DISTINCT(gene_id) FROM transcript WHERE biotype="low_coverage"'
              .') GROUP BY g.gene_id HAVING biotype_cnt > 1'
            .') AS gbio,transcript t WHERE t.gene_id=gbio.gene_id AND t.transcript_id IN (SELECT transcript_id FROM transcript WHERE biotype="low_coverage")'
          .')',
          'UPDATE gene SET biotype="low_coverage" WHERE gene_id IN (SELECT gene_id FROM (SELECT count(*) AS cnt,gene_id,biotype FROM (SELECT count(*),gene_id,biotype FROM transcript GROUP BY gene_id,biotype) AS numbioperg GROUP BY gene_id having cnt=1) AS geneswith1bio WHERE biotype="low_coverage")'
        ],
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['update_rnaseq_ise_logic_names'],
      },
    },

    {
      -logic_name => 'update_rnaseq_ise_logic_names',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -parameters => {
        db_conn => $self->o('final_geneset_db'),
        sql     => [
          'UPDATE analysis SET logic_name = REPLACE(logic_name, "_rnaseq_gene", "_rnaseq_ise")',
        ],
      },
      -rc_name   => 'default',
      -flow_into => {
        1 => ['run_cleaner'],
      },
    },

    {
      -logic_name => 'run_cleaner',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCleanGeneset',
      -parameters => {
        skip_analysis                         => $self->o('skip_cleaning'),
        input_db                              => $self->o('final_geneset_db'),
        dna_db                                => $self->o('dna_db'),
        output_path                           => $self->o('output_path') . '/clean_genes/',
        blessed_biotypes                      => $self->o('cleaning_blessed_biotypes'),
        flagged_redundancy_coverage_threshold => 95,
        general_redundancy_coverage_threshold => 95,
      },
      -rc_name   => '8GB',
      -flow_into => {
        '1' => ['delete_flagged_transcripts'],
      },
    },

    {
      -logic_name => 'delete_flagged_transcripts',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDeleteTranscripts',
      -parameters => {
        skip_analysis       => $self->o('skip_cleaning'),
        dbhost              => $self->o( 'final_geneset_db', '-host' ),
        dbname              => $self->o( 'final_geneset_db', '-dbname' ),
        dbuser              => $self->o('user'),
        dbpass              => $self->o('password'),
        dbport              => $self->o( 'final_geneset_db', '-port' ),
        transcript_ids_file => catfile( $self->o('output_path'), 'clean_genes', 'transcript_ids_to_remove.txt' ),
        delete_transcripts_path => catdir( $self->o('ensembl_analysis_script'), 'genebuild/' ),
        delete_genes_path       => catdir( $self->o('ensembl_analysis_script'), 'genebuild/' ),
        delete_transcripts_script_name => '/delete_transcripts.pl',
        delete_genes_script_name       => '/delete_genes.pl',
        output_path                    => catdir( $self->o('output_path'), 'clean_genes' ),
        output_file_name               => 'delete_transcripts.out',
      },
      -max_retry_count => 0,
      -flow_into       => {
        1 => ['transfer_ncrnas'],
      },
    },

    {
      -logic_name => 'transfer_ncrnas',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => 'perl ' . catfile( $self->o('ensembl_analysis_script'), 'genebuild', 'copy_genes.pl' ) .
          ' -sourcehost ' . $self->o( 'ncrna_db', '-host' ) .
          ' -sourceuser ' . $self->o('user_r') .
          ' -sourceport ' . $self->o( 'ncrna_db', '-port' ) .
          ' -sourcedbname ' . $self->o( 'ncrna_db', '-dbname' ) .
          ' -dnauser ' . $self->o('user_r') .
          ' -dnahost ' . $self->o( 'dna_db', '-host' ) .
          ' -dnaport ' . $self->o( 'dna_db', '-port' ) .
          ' -dnadbname ' . $self->o( 'dna_db', '-dbname' ) .
          ' -targetuser ' . $self->o('user') .
          ' -targetpass ' . $self->o('password') .
          ' -targethost ' . $self->o( 'final_geneset_db', '-host' ) .
          ' -targetport ' . $self->o( 'final_geneset_db', '-port' ) .
          ' -targetdbname ' . $self->o( 'final_geneset_db', '-dbname' ) .
          ' -all'
      },
      -rc_name   => 'default',
      -flow_into => {
        '1' => ['delete_duplicate_genes'],
      },
    },

    {
      -logic_name => 'delete_duplicate_genes',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -parameters => {
        cmd => "perl " . $self->o('remove_duplicates_script_path')
          . " -dbhost " . $self->o( 'final_geneset_db', '-host' )
          . " -dbuser " . $self->o('user')
          . " -dbpass " . $self->o('password')
          . " -dbname " . $self->o( 'final_geneset_db', '-dbname' )
          . " -dbport " . $self->o( 'final_geneset_db', '-port' )
          . " -dnadbhost " . $self->o( 'dna_db', '-host' )
          . " -dnadbuser " . $self->o('user_r')
          . " -dnadbname " . $self->o( 'dna_db', '-dbname' )
          . " -dnadbport " . $self->o( 'dna_db', '-port' ),
      },
      -max_retry_count => 0,
      -rc_name         => '2GB',
      -flow_into       => {
        '1' => ['final_db_sanity_checks'],
      },
    },

    {
      -logic_name => 'final_db_sanity_checks',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAnalysisSanityCheck',
      -parameters => {
        skip_rnaseq                => $self->o('skip_rnaseq'),
        skip_projection            => $self->o('skip_projection'),
        target_db                  => $self->o('final_geneset_db'),
        sanity_check_type          => 'gene_db_checks',
        min_allowed_feature_counts => get_analysis_settings( 'Bio::EnsEMBL::Analysis::Hive::Config::SanityChecksStatic',
          'gene_db_checks' )->{ $self->o('uniprot_set') }->{'final'},
      },
      -rc_name => '4GB',
    },

  ];
}

sub resource_classes {
  my $self = shift;

  return {
    '2GB'     => { LSF => $self->lsf_resource_builder( 'production', 2000 ) },
    '3GB'     => { LSF => $self->lsf_resource_builder( 'production', 3000 ) },
    '4GB'     => { LSF => $self->lsf_resource_builder( 'production', 4000 ) },
    '5GB'     => { LSF => $self->lsf_resource_builder( 'production', 5000 ) },
    '8GB'     => { LSF => $self->lsf_resource_builder( 'production', 8000 ) },
    '10GB'    => { LSF => $self->lsf_resource_builder( 'production', 10000 ) },
    '30GB'    => { LSF => $self->lsf_resource_builder( 'production', 30000 ) },
    'default' => { LSF => $self->lsf_resource_builder( 'production', 900 ) },
    }
}

1;
