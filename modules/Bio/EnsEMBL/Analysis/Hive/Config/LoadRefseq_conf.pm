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

Bio::EnsEMBL::Analysis::Hive::Config::LoadRefseq_conf

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::Config::LoadRefseq_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Production::Pipeline::PipeConfig::Base_conf');

sub default_options {
  my ($self) = @_;

  return {
    %{$self->SUPER::default_options},
    species      => [],
    antispecies  => [],
    division     => [],
    run_all      => 0,
    meta_filters => {},
  };
}

sub pipeline_analyses {
  my ($self) = @_;

  return [
    {
      -logic_name      => 'DbFactory',
      -module          => 'Bio::EnsEMBL::Production::Pipeline::Common::DbFactory',
      -max_retry_count => 1,
      -parameters      => {
        species         => $self->o('species'),
        antispecies     => $self->o('antispecies'),
        division        => $self->o('division'),
        run_all         => $self->o('run_all'),
        meta_filters    => $self->o('meta_filters'),
        chromosome_flow => 0,
        regulation_flow => 0,
        variation_flow  => 0,
      },
      -flow_into       => {
        '2->A' => ['BackupTables'],
        'A->2' => ['AnnotateProteinFeatures'],
      },
      -meadow_type     => 'LOCAL',
    },
    {
      -logic_name      => 'check_synonyms_loaded',
      -module          => 'Bio::EnsEMBL::Production::Pipeline::Common::DbFactory',
      -max_retry_count => 1,
      -parameters      => {
        species         => $self->o('species'),
        antispecies     => $self->o('antispecies'),
        division        => $self->o('division'),
        run_all         => $self->o('run_all'),
        meta_filters    => $self->o('meta_filters'),
        chromosome_flow => 0,
        regulation_flow => 0,
        variation_flow  => 0,
      },
      -flow_into       => {
        '1->A' => ['BackupTables'],
        'A->1' => ['AnnotateProteinFeatures'],
      },
      -meadow_type     => 'LOCAL',
    },
    {
      -logic_name      => 'load_gff',
      -module          => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadRefSeqGFF3',
      -max_retry_count => 1,
      -parameters      => {
      },
      -flow_into       => {
        '2->A' => ['BackupTables'],
        'A->2' => ['AnnotateProteinFeatures'],
      },
      -meadow_type     => 'LOCAL',
    },
    {
      -logic_name      => 'DbFactory',
      -module          => 'Bio::EnsEMBL::Production::Pipeline::Common::DbFactory',
      -max_retry_count => 1,
      -parameters      => {
        species         => $self->o('species'),
        antispecies     => $self->o('antispecies'),
        division        => $self->o('division'),
        run_all         => $self->o('run_all'),
        meta_filters    => $self->o('meta_filters'),
        chromosome_flow => 0,
        regulation_flow => 0,
        variation_flow  => 0,
      },
      -flow_into       => {
        '2->A' => ['BackupTables'],
        'A->2' => ['AnnotateProteinFeatures'],
      },
      -meadow_type     => 'LOCAL',
    },
  ];
}


1;

