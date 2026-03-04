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

Base configuration for eHive pipelines that orchestrate Nextflow pipelines
via HiveRunNextflow. Inheriting configs set nextflow_pipeline_dir,
nextflow_pipeline_name, nextflow_work_root, and nextflow_output_root, then
define their pipeline_analyses using HiveRunNextflow as the module.

=head2 Runtime directory layout

  ${nextflow_work_root}/
    ${nextflow_pipeline_name}/
      job_${hive_job_id}/

        # resume_mode = never or attempt:
        attempt_0/
          work/            <-- Nextflow workDir
          manifest.json    <-- run metadata written by HiveRunNextflow
        attempt_1/         <-- on eHive retry; attempt mode gets fresh dir
          work/
          manifest.json

        # resume_mode = job:
        work/              <-- shared workDir across all retries
        manifest.json

  ${nextflow_output_root}/
    ${nextflow_pipeline_name}/
      <pipeline-specific publishDir output>

=head2 Resource classes

Two launcher resource classes are added on top of the base classes:

  nextflow_launcher        2 GB, 7 days  -- typical pipeline
  nextflow_launcher_himem  8 GB, 7 days  -- if the NF head process needs more

The launcher job is lightweight (it just submits to SLURM and polls), but
it must outlast the entire Nextflow pipeline run.

=head1 METHODS

=head2 default_options

Returns the default options, extending HiveBaseConfig_conf with:

  nextflow_binary        - nextflow executable (default: 'nextflow')
  nextflow_resume_mode   - never | attempt | job (default: 'attempt')
  nextflow_work_root     - base path for per-job work directories
  nextflow_output_root   - base path for pipeline output
  nextflow_pipeline_dir  - path to directory containing main.nf
  nextflow_pipeline_name - identifier used in directory naming
  nextflow_profile       - Nextflow -profile string (default: 'slurm,singularity')
  nextflow_extra_flags   - arrayref of verbatim CLI flags

=cut

package Bio::EnsEMBL::Analysis::Hive::Config::NextflowPipelineBase_conf;

use strict;
use warnings;
use feature 'say';

use File::Spec::Functions qw(catdir);

use base ('Bio::EnsEMBL::Analysis::Hive::Config::HiveBaseConfig_conf');


sub default_options {
    my ($self) = @_;
    return {
        %{ $self->SUPER::default_options() },

        # --- Nextflow binary ---
        nextflow_binary       => 'nextflow',

        # --- Resume behaviour ---
        # never:   always run fresh; -resume never passed
        # attempt: -resume within same eHive retry attempt; fresh dir on retry
        # job:     -resume always; single workDir across all retries
        nextflow_resume_mode  => 'attempt',

        # --- Directory roots (must be set by the inheriting config) ---
        nextflow_work_root    => undef,   # base for Nextflow work directories
        nextflow_output_root  => undef,   # base for published pipeline output

        # --- Pipeline identity (must be set by the inheriting config) ---
        nextflow_pipeline_dir  => undef,  # path to the directory containing main.nf
        nextflow_pipeline_name => undef,  # used in directory naming

        # --- Execution profile ---
        # Common values: 'slurm,singularity'  'local,docker'  'slurm,conda'
        nextflow_profile      => 'slurm,singularity',

        # Extra Nextflow CLI flags passed verbatim, e.g. ['-ansi-log', 'false']
        nextflow_extra_flags  => [],
    };
}


sub resource_classes {
    my $self = shift;

    return {
        %{ $self->SUPER::resource_classes() },

        # Launcher jobs: small memory, very long wall-time to outlast the pipeline
        'nextflow_launcher' => {
            SLURM => $self->slurm_resource_builder(2000, '7-00:00:00'),
        },
        'nextflow_launcher_himem' => {
            SLURM => $self->slurm_resource_builder(8000, '7-00:00:00'),
        },
    };
}


1;
