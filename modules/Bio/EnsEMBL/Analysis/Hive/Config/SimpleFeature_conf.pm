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

Bio::EnsEMBL::Analysis::Hive::Config::SimpleFeature_conf

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::Config::SimpleFeature_conf;

use strict;
use warnings;

use File::Spec::Functions;

use base ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');

sub default_options {
  my ($self) = @_;
  return {

    # inherit other stuff from the base class
    %{ $self->SUPER::default_options() },

######################################################
#
# Variable settings - You change these!!!
#
######################################################
    dbowner             => '' || $ENV{EHIVE_USER} || $ENV{USER},
    pipeline_name       => '' || $self->o('production_name') . '_' . $self->o('ensembl_release'),
    assembly_name       => '',
    dbname_accession    => '', # This is the assembly accession without [._] and all lower case, i.e gca001857705v1
    user_r              => '',                                  # read only db user
    user                => '',                                  # write db user
    password            => '',                                  # password for write db user
    pipe_db_host        => '',                                  # host for pipe db
    pipe_db_port        => '',                                  # port for pipeline host
    dna_db_host         => '',                                  # host for dna db
    dna_db_port         => '',                                  # port for dna db host
    release_number      => '' || $self->o('ensembl_release'),
    species_name        => '',                                  # e.g. mus_musculus
    production_name     => '',
    output_path         => '',                                  # Lustre output dir. This will be the primary dir to house the assembly info and various things from analyses
    use_genome_flatfile => '1',                                 # This will read sequence where possible from a dumped flatfile instead of the core db

######################################################
#
# Mostly constant settings
#
######################################################

########################
# Pipe and ref db info
########################
    pipe_db_name => $self->o('dbowner') . '_' . $self->o('dbname_accession') . '_simplefeature_pipe_' . $self->o('release_number'),
    dna_db_name  => $self->o('dbowner') . '_' . $self->o('dbname_accession') . '_core_' . $self->o('release_number'),

    reference_db_name   => $self->o('dna_db_name'),
    reference_db_host   => $self->o('dna_db_host'),
    reference_db_port   => $self->o('dna_db_port'),

########################
# Misc setup info
########################
    # This is used for the ensembl_production and the ncbi_taxonomy databases
    ensembl_release      => $ENV{ENSEMBL_RELEASE},     # this is the current release version on staging to be able to get the correct database

    min_toplevel_slice_length => 250,

    genome_dumps => catdir( $self->o('output_path'), 'genome_dumps' ),
# This one is used in replacement of the dna table in the core db, so where analyses override slice->seq.
# Has simple headers with just the seq_region name. Also used by bwa in the RNA-seq analyses. Not masked
    faidx_genome_file => catfile( $self->o('genome_dumps'), $self->o('species_name') . '_toplevel.fa' ),

    ensembl_analysis_script => catdir($self->o('enscode_root_dir'), 'ensembl-analysis', 'scripts'),
    sequence_dump_script    => catfile($self->o('ensembl_analysis_script'), 'sequence_dump.pl'),


########################
# Executable paths
########################
    samtools_path => catfile($self->o('binary_base'), 'samtools'), #You may need to specify the full path to the samtools binary
    eponine_java_path => 'java',
    eponine_jar_path => '/hps/software/users/ensembl/ensw/swenv/env/prod-base/bin/eponine-scan.jar',
    cpg_path => catfile($self->o('binary_base'), 'cpg_lh'),
    trnascan_path => catfile($self->o('binary_base'), 'tRNAscan-SE'),

########################
# db info
########################
    reference_db => {
      -dbname => $self->o('reference_db_name'),
      -host   => $self->o('reference_db_host'),
      -port   => $self->o('reference_db_port'),
      -user   => $self->o('user'),
      -pass   => $self->o('password'),
      -driver => $self->o('hive_driver'),
    },

  };
}


sub pipeline_wide_parameters {
  my ($self) = @_;

  return {
    %{ $self->SUPER::pipeline_wide_parameters },
    use_genome_flatfile     => $self->o('use_genome_flatfile'),
    genome_file             => $self->o('faidx_genome_file'),
  }
}


sub pipeline_analyses {
  my ($self) = @_;

  return [
    {
      -logic_name => 'create_faidx_genome_file',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => 'default',
      -parameters => {
        cmd => 'if [ ! -s "'.$self->o('faidx_genome_file').'" ]; then mkdir -p '.$self->o('genome_dumps').'; perl '.$self->o('sequence_dump_script').' -dbhost '.$self->o('dna_db_host').' -dbuser '.$self->o('dna_db_user').' -dbport '.$self->o('dna_db_port').' -dbname '.$self->o('dna_db_name').' -coord_system_name '.$self->o('assembly_name').' -toplevel -onefile -header rnaseq -filename '.$self->o('faidx_genome_file').';fi',
      },
      -input_ids  => [{}],
      -flow_into => {
        1 => [ 'create_faidx'],
      },
    },

    {
      -logic_name => 'create_faidx',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
      -rc_name => '5GB',
      -parameters => {
        cmd => 'if [ ! -e "'.$self->o('faidx_genome_file').'.fai" ]; then '.$self->o('samtools_path').' faidx '.$self->o('faidx_genome_file').';fi',
      },
      -flow_into => {
        1 => [ 'create_10mb_slice_ids'],
      },
    },
    {
      -logic_name => 'create_10mb_slice_ids',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSubmitAnalysis',
      -parameters => {
        target_db             => $self->o('dna_db'),
        coord_system_name     => 'toplevel',
        iid_type              => 'slice',
        slice_size            => 10000000,
        include_non_reference => 0,
        top_level             => 1,
        min_slice_length      => $self->o('min_toplevel_slice_length'),
        batch_slice_ids       => 1,
        batch_target_size     => 10000000,
      },
      -rc_name   => '3GB',
      -flow_into => {
        2 => ['eponine'],
      },
    },

    {
      -logic_name => 'eponine',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveEponine',
      -parameters => {
        target_db => $self->o('reference_db'),
        logic_name => 'eponine',
        module => 'HiveEponine',
        eponine_path => $self->o('eponine_java_path'),
        commandline_params => '-epojar => '.$self->o('eponine_jar_path').', -threshold => 0.999',
      },
      -rc_name    => '3GB',
      -flow_into => {
        MAIN => ['cpg'],
        ANYFAILURE => ['cpg'],
      },
      -hive_capacity => $self->o('hc_normal'),
      -batch_size => 20,
    },

    {
      -logic_name => 'cpg',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveCPG',
      -parameters => {
        target_db => $self->o('reference_db'),
        logic_name => 'cpg',
        module => 'HiveCPG',
        cpg_path => $self->o('cpg_path'),
      },
      -rc_name    => '3GB',
      -flow_into => {
        MAIN => ['trnascan'],
        ANYFAILURE => ['trnascan'],
      },
      -hive_capacity => $self->o('hc_normal'),
      -batch_size => 20,
    },

    {
      -logic_name => 'trnascan',
      -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveTRNAScan',
      -parameters => {
        target_db => $self->o('reference_db'),
        logic_name => 'trnascan',
        module => 'HiveTRNAScan',
        trnascan_path => $self->o('trnascan_path'),
      },
      -rc_name    => '3GB',
      -hive_capacity => $self->o('hc_normal'),
      -batch_size => 20,
    },
  ];
}

sub resource_classes {
  my $self = shift;

  return {
    'default' => {
      LSF => $self->lsf_resource_builder( 'production', 1000 ),
      SLURM => $self->slurm_resource_builder( 'standard', 1000, '30:00'),
    },
    '3GB' => {
      LSF => $self->lsf_resource_builder( 'production', 3000 ),
      SLURM => $self->slurm_resource_builder( 'standard', 3000, '1:00:00'),
    },
    '5GB' => {
      LSF => $self->lsf_resource_builder( 'production', 5000 ),
      SLURM => $self->slurm_resource_builder( 'standard', 5000, '5'),
    },
  }
}

1;
