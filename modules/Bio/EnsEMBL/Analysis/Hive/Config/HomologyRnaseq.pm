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
=cut

package Bio::EnsEMBL::Analysis::Hive::Config::HomologyRnaseq;

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
    'dbowner'                         => '' || $ENV{EHIVE_USER} || $ENV{USER},
    'pipeline_name'                   => '' || $self->o('production_name').'_'.$self->o('ensembl_release'),
    'user_r'                          => '', # read only db user
    'user'                            => '', # write db user
    'password'                        => '', # password for write db user
    'pipe_db_host'                    => '', # host for pipe db
    'databases_host'                  => '', # host for general output dbs
    'dna_db_host'                     => '', # host for dna db
    'pipe_db_port'                    => '', # port for pipeline host
    'databases_port'                  => '', # port for general output db host
    'dna_db_port'                     => '', # port for dna db host
    'release_number'                  => '' || $self->o('ensembl_release'),
    'species_name'                    => '', # e.g. mus_musculus
    'production_name'                 => '', # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
    dbname_accession          => '', # This is the assembly accession without [._] and all lower case, i.e gca001857705v1
    'uniprot_set'                     => '', # e.g. mammals_basic, check UniProtCladeDownloadStatic.pm module in hive config dir for suitable set,
    'use_genome_flatfile'             => '1',# This will read sequence where possible from a dumped flatfile instead of the core db
    'output_path'                     => '', # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses

########################
# Pipe and ref db info
########################

    'pipe_db_name'                  => $self->o('dbowner').'_'.$self->o('dbname_accession').'_pipe_'.$self->o('release_number'),
    'dna_db_name'                   => $self->o('dbowner').'_'.$self->o('dbname_accession').'_core_'.$self->o('release_number'),

    'genblast_db_host'             => $self->o('databases_host'),
    'genblast_db_port'             => $self->o('databases_port'),

    'genblast_rnaseq_support_db_host'    => $self->o('databases_host'),
    'genblast_rnaseq_support_db_port'    => $self->o('databases_port'),

    'rnaseq_refine_db_host'         => $self->o('databases_host'),
    'rnaseq_refine_db_port'         => $self->o('databases_port'),
    'rnaseq_refine_db_name'         => $self->o('dbowner').'_'.$self->o('dbname_accession').'_refine_'.$self->o('release_number'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    'ensembl_release'              => $ENV{ENSEMBL_RELEASE}, # this is the current release version on staging to be able to get the correct database

    databases_to_delete => ['genblast_rnaseq_support_db', 'genblast_rnaseq_support_nr_db'],

######################################################
#
# Mostly constant settings
#
######################################################
    'min_toplevel_slice_length'   => 250,
    genome_dumps                  => catdir($self->o('output_path'), 'genome_dumps'),
    faidx_genome_file             => catfile($self->o('genome_dumps'), $self->o('species_name').'_toplevel.fa'),
    # This is used for "messaging" other sub pipeline
    transcript_selection_url => undef,


# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# No option below this mark should be modified
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

########################
# db info
########################
    'genblast_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('dbname_accession').'_genblast_'.$self->o('release_number'),
      -host   => $self->o('genblast_db_host'),
      -port   => $self->o('genblast_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'genblast_rnaseq_support_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('dbname_accession').'_gb_rnaseq_'.$self->o('release_number'),
      -host   => $self->o('genblast_rnaseq_support_db_host'),
      -port   => $self->o('genblast_rnaseq_support_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'genblast_rnaseq_support_nr_db' => {
      -dbname => $self->o('dbowner').'_'.$self->o('dbname_accession').'_gb_rnaseq_nr_'.$self->o('release_number'),
      -host   => $self->o('genblast_rnaseq_support_db_host'),
      -port   => $self->o('genblast_rnaseq_support_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    'rnaseq_refine_db' => {
      -dbname => $self->o('rnaseq_refine_db_name'),
      -host   => $self->o('rnaseq_refine_db_host'),
      -port   => $self->o('rnaseq_refine_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },
  };
}

sub pipeline_create_commands {
  my ($self) = @_;
  return [
    # inheriting database and hive tables' creation
    @{$self->SUPER::pipeline_create_commands},
  ];
}

sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{$self->SUPER::pipeline_wide_parameters},
    use_genome_flatfile  => $self->o('use_genome_flatfile'),
    genome_file          => $self->o('faidx_genome_file'),
  }
}

## See diagram for pipeline structure
sub pipeline_analyses {
  my ($self) = @_;
  return [

###############################################################################
#
# Homology Pipeline with RNAseq support
#
###############################################################################

    {
      -logic_name => 'create_genblast_rnaseq_support_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('dna_db'),
        target_db   => $self->o('genblast_rnaseq_support_db'),
        create_type => 'clone',
      },
      -rc_name    => 'default',
      -input_ids  => [{}],
      -flow_into  => {
        '1->A' => ['create_genblast_rnaseq_slice_ids'],
        'A->1' => ['create_genblast_rnaseq_nr_db'],
      },
    },

    {
      -logic_name => 'create_genblast_rnaseq_slice_ids',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db             => $self->o('dna_db'),
        coord_system_name     => 'toplevel',
        iid_type              => 'slice',
        include_non_reference => 0,
        top_level             => 1,
        min_slice_length      => 0,
      },
        -rc_name    => 'default',
        -flow_into  => {
          '2' => ['genblast_rnaseq_support'],
      },
    },

    {
      -logic_name => 'genblast_rnaseq_support',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveHomologyRNASeqIntronsCheck',
      -parameters => {
        dna_db            => $self->o('dna_db'),
        source_db         => $self->o('genblast_db'),
        intron_db         => $self->o('rnaseq_refine_db'),
        target_db         => $self->o('genblast_rnaseq_support_db'),
        logic_name        => 'genblast_rnaseq_support',
        classify_by_count => 1,
        update_genes      => 0,
        module => 'HiveHomologyRNASeqIntronsCheck',
      },
      -rc_name    => '4GB',
      -flow_into  => {
        '-1' => ['genblast_rnaseq_support_himem'],
      },
     },

    {
      -logic_name => 'genblast_rnaseq_support_himem',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveHomologyRNASeqIntronsCheck',
      -parameters => {
        dna_db            => $self->o('dna_db'),
        source_db         => $self->o('genblast_db'),
        intron_db         => $self->o('rnaseq_refine_db'),
        target_db         => $self->o('genblast_rnaseq_support_db'),
        logic_name        => 'genblast_rnaseq_support',
        classify_by_count => 1,
        update_genes      => 0,
        module            => 'HiveHomologyRNASeqIntronsCheck',
      },
        -rc_name    => '8GB',
    },

    {
      -logic_name => 'create_genblast_rnaseq_nr_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('genblast_rnaseq_support_db'),
        target_db   => $self->o('genblast_rnaseq_support_nr_db'),
        create_type => 'copy',
      },
      -rc_name    => 'default',
      -flow_into  => {
        '1' => ['create_genblast_rnaseq_nr_slices'],
      },
    },

    {
      -logic_name => 'create_genblast_rnaseq_nr_slices',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db              => $self->o('dna_db'),
        coord_system_name      => 'toplevel',
        iid_type               => 'slice',
        slice_size             => 20000000,
        include_non_reference  => 0,
        top_level              => 1,
        min_slice_length       => $self->o('min_toplevel_slice_length'),
        batch_slice_ids        => 1,
        batch_target_size      => 20000000,
      },
      -rc_name    => '2GB',
      -flow_into  => {
        '2->A'  => ['remove_redundant_genblast_rnaseq_genes'],
        'A->1'  => ['notification_pipeline_is_done'],
      },
    },

    {
      -logic_name => 'remove_redundant_genblast_rnaseq_genes',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::RemoveRedundantGenes',
      -parameters => {
        target_db   => $self->o('genblast_rnaseq_support_nr_db'),
        target_type => 'biotype_priority',
        layers      => get_analysis_settings('Bio::EnsEMBL::Analysis::Hive::Config::LayerAnnotationStatic', $self->o('uniprot_set'), undef, 'ARRAY'),
      },
      -hive_capacity => $self->o('hc_normal'),
      -rc_name    => '5GB',
    },

    {
      -logic_name => 'notification_pipeline_is_done',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::MessagePipeline',
      -parameters => {
        tweak_script => catfile($self->o('enscode_root_dir'), 'ensembl-hive', 'scripts', 'tweak_pipeline.pl'),
        messages   => [
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'layer_annotation',
          param => 'SOURCEDB_REFS',
          data => $self->o('genblast_rnaseq_support_nr_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'create_toplevel_slices',
          param => 'feature_dbs',
          data => $self->o('genblast_rnaseq_support_nr_db'),
          update => 1,
        },
        {
          url => $self->o('transcript_selection_url'),
          logic_name => 'split_slices_on_intergenic',
          param => 'input_gene_dbs',
          data => $self->o('genblast_rnaseq_support_nr_db'),
          update => 1,
        }],
      },
      -rc_name    => 'default',
    },
  ];
}

sub resource_classes {
  my $self = shift;

  return {
    #inherit other stuff from the base class
     %{ $self->SUPER::resource_classes() },
  }
}
1;
