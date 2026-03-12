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

=head1 DESCRIPTION

Minimal eHive pipeline for smoke-testing the HiveRunNextflow runnable
against the ensembl-genes-nf example pipeline. Runs locally with no
containers or SLURM -- if this works, the eHive/Nextflow plumbing is good.

=head2 Typical init command

  init_pipeline.pl Bio::EnsEMBL::Analysis::Hive::Config::NextflowExample_conf \
    -pipeline_db -host=<host> -pipeline_db -port=<port>          \
    -pipeline_db -user=$GBUSER -pipeline_db -pass=$GBPASS        \
    -pipeline_db -dbname=${USER}_nextflow_example_pipe            \
    -nextflow_work_root   /path/to/nextflow_work                  \
    -nextflow_output_root /path/to/nextflow_output                \
    -nextflow_pipeline_dir /path/to/ensembl-genes-nf/pipelines/example \
    -example_input /any/existing/file

=cut

package Bio::EnsEMBL::Analysis::Hive::Config::NextflowExample_conf;

use strict;
use warnings;
use feature 'say';

use File::Spec::Functions qw(catdir);

use base ('Bio::EnsEMBL::Analysis::Hive::Config::NextflowPipelineBase_conf');


sub default_options {
    my ($self) = @_;
    return {
        %{ $self->SUPER::default_options() },

        # --- eHive DB credentials ---
        user        => $ENV{GBUSER},
        password    => $ENV{GBPASS},
        user_r      => $ENV{USER_R} || 'ensro',
        password_r  => undef,
        dna_db_name => '',

        # --- eHive pipeline identity ---
        pipeline_name          => 'nextflow_example',

        # --- Nextflow pipeline identity ---
        nextflow_pipeline_name => 'example',
        nextflow_pipeline_dir  => '/path/to/ensembl-genes-nf/pipelines/example',

        # Run locally -- the example pipeline defines no profiles
        nextflow_profile => undef,
        nextflow_resume_mode => 'attempt',

        # --- Example pipeline params ---
        # --input is required by the schema; any existing file satisfies it
        example_input => undef,
    };
}


sub pipeline_analyses {
    my ($self) = @_;

    return [

        {
            -logic_name => 'run_example',
            -module     => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRunNextflow',
            -parameters => {
                nextflow_binary        => $self->o('nextflow_binary'),
                nextflow_pipeline_dir  => $self->o('nextflow_pipeline_dir'),
                nextflow_pipeline_name => $self->o('nextflow_pipeline_name'),
                nextflow_work_root     => $self->o('nextflow_work_root'),
                nextflow_output_dir    => catdir(
                    $self->o('nextflow_output_root'),
                    $self->o('nextflow_pipeline_name'),
                ),
                nextflow_resume_mode   => $self->o('nextflow_resume_mode'),
                nextflow_profile       => $self->o('nextflow_profile'),
                nextflow_extra_flags   => $self->o('nextflow_extra_flags'),
                nextflow_params        => {
                    input => $self->o('example_input'),
                },
            },
            -rc_name     => '1GB',     # local run, minimal resources
            -input_ids   => [{}],
            -flow_into   => { 1 => ['pipeline_complete'] },
        },

        {
            -logic_name => 'pipeline_complete',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -rc_name    => 'default',
        },

    ];
}


1;
