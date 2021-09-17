=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::Config::CleanUTRs

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::Config::CleanUTRs;

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
    dbowner => '' || $ENV{EHIVE_USER} || $ENV{USER},
    pipeline_name => '' || $self->o('production_name') . '_' . $self->o('ensembl_release'),
    release_number => '' || $self->o('ensembl_release'),
    production_name => '',     # usually the same as species name but currently needs to be a unique entry for the production db, used in all core-like db names
    desired_slice_length => 10000000,
    store_rejected => 0,

    ########################
    # Pipe and ref db info
    ########################
    # The following might not be known in advance, since the come from other pipelines
    # These values can be replaced in the analysis_base table if they're not known yet
    # If they are not needed (i.e. no projection or rnaseq) then leave them as is

    pipe_db_name => $self->o('dbowner').'_'.$self->o('production_name').'_clean_utr_pipe_'.$self->o('release_number'),
    dna_db_name  => $self->o('dbowner').'_'.$self->o('production_name').'_core_'.$self->o('release_number'),

    # This is used for the ensembl_production and the ncbi_taxonomy databases
    ensembl_release => $ENV{ENSEMBL_RELEASE},    # this is the current release version on staging to be able to get the correct database
    target_db_exists => 0,

    ######################################################
    #
    # Mostly constant settings
    #
    ######################################################
    gene_db_host => $self->o('databases_host'),
    gene_db_port => $self->o('databases_port'),
    gene_db_name => $self->o('dbowner').'_'.$self->o('production_name').'_gene_'.$self->o('release_number'),

    clean_utr_db_host => $self->o('databases_host'),
    clean_utr_db_port => $self->o('databases_port'),
    clean_utr_db_name => $self->o('dbowner').'_'.$self->o('production_name').'_clean_utr_'.$self->o('release_number'),

    ########################
    # db info
    ########################
    gene_db => {
      -dbname => $self->o('gene_db_name'),
      -host   => $self->o('gene_db_host'),
      -port   => $self->o('gene_db_port'),
      -user   => $self->o('user_r'),
      -pass   => $self->o('password_r'),
      -driver => $self->o('hive_driver'),
    },

    clean_utr_db => {
      -dbname => $self->o('clean_utr_db_name'),
      -host   => $self->o('clean_utr_db_host'),
      -port   => $self->o('clean_utr_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

    databases_to_delete => ['clean_utr_db'],
  };
}


## See diagram for pipeline structure
sub pipeline_analyses {
  my ($self) = @_;

  return [
    {
      -logic_name => 'create_cleanutr_db',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDatabase',
      -parameters => {
        source_db   => $self->o('gene_db'),
        target_db   => $self->o('clean_utr_db'),
        create_type => 'clone',
        skip_analysis => $self->o('target_db_exists'),
      },
      -rc_name   => 'default',
      -input_ids => [{}],
      -flow_into => {
        1 => ['create_toplevel_slices'],
      },
    },

    {
      -logic_name => 'create_toplevel_slices',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db             => $self->o('gene_db'),
        iid_type              => 'slice',
        coord_system_name     => 'toplevel',
        include_non_reference => 0,
        top_level             => 1,
      },
      -flow_into => {
        2 => ['split_slices_on_intergenic'],
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'split_slices_on_intergenic',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db             => $self->o('gene_db'),
        iid_type              => 'intergenic_slice',
        coord_system_name     => 'toplevel',
        include_non_reference => 0,
        top_level             => 1,
        desired_slice_length  => $self->o('desired_slice_length'),
      },
      -batch_size => 300,
      -flow_into => {
        2 => ['clean_utr'],
      },
      -rc_name => 'default',
    },

    {
      -logic_name => 'clean_utr',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::CleanUTRs',
      -parameters => {
        source_db => $self->o('gene_db'),
        target_db => $self->o('clean_utr_db'),
        dna_db    => $self->o('dna_db'),
        store_rejected => $self->o('store_rejected'),
      },
      -rc_name => '4GB',
      -max_retry_count => 1,
      -analysis_capacity => 400,
    },

  ];
}

sub resource_classes {
  my $self = shift;

  return {
    'default' => { LSF => $self->lsf_resource_builder( 'production', 2000 ) },
    '4GB' => { LSF => $self->lsf_resource_builder( 'production', 4000 ) },
    '8GB' => { LSF => $self->lsf_resource_builder( 'production', 8000 ) },
  }
}

1;

