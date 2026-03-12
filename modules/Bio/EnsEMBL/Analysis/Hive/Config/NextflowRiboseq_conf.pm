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

eHive pipeline configuration for running the ensembl-genes-nf riboseq
Nextflow pipeline via HiveRunNextflow.

=head2 Typical init command

  init_pipeline.pl Bio::EnsEMBL::Analysis::Hive::Config::NextflowRiboseq_conf \
    -pipeline_db -host=<host> -pipeline_db -port=<port>             \
    -pipeline_db -dbname=<user>_riboseq_pipe                         \
    -nextflow_work_root   /hps/nobackup/flicek/ensembl/genebuild/nextflow_work  \
    -nextflow_output_root /hps/nobackup/flicek/ensembl/genebuild/nextflow_output \
    -nextflow_pipeline_dir /path/to/ensembl-genes-nf/pipelines/riboseq          \
    -riboseq_sample_sheet /path/to/samples.csv                       \
    -riboseq_fasta        /path/to/genome.fa                         \
    -riboseq_star_index   /path/to/star_index                        \
    -riboseq_gtf          /path/to/annotations.gtf

=head2 Directory layout produced

  ${nextflow_work_root}/riboseq/job_${id}/attempt_0/work/
  ${nextflow_work_root}/riboseq/job_${id}/attempt_0/manifest.json

  ${nextflow_output_root}/riboseq/
    <riboseq pipeline publishDir output: collapsed_reads/, star_genome/, etc.>

=cut

package Bio::EnsEMBL::Analysis::Hive::Config::NextflowRiboseq_conf;

use strict;
use warnings;
use feature 'say';

use File::Spec::Functions qw(catdir);

use base ('Bio::EnsEMBL::Analysis::Hive::Config::NextflowPipelineBase_conf');


sub default_options {
    my ($self) = @_;
    return {
        %{ $self->SUPER::default_options() },

        # --- eHive DB credentials (pulled from standard genebuild env vars) ---
        user       => $ENV{GBUSER},
        password   => $ENV{GBPASS},
        user_r     => $ENV{USER_R} || 'ensro',
        password_r => undef,
        dna_db_name => '',   # no Ensembl core DB needed; stub to satisfy base class

        # --- eHive pipeline identity ---
        pipeline_name          => 'nextflow_riboseq',

        # --- Nextflow pipeline identity ---
        nextflow_pipeline_name => 'riboseq',
        nextflow_pipeline_dir  => '/path/to/ensembl-genes-nf/pipelines/riboseq',

        # --- Resume strategy ---
        # 'attempt' is the safe default: resume within an attempt, fresh on retry
        nextflow_resume_mode   => 'attempt',

        # Override base default of 'slurm,singularity' -- the riboseq pipeline's
        # slurm profile already enables Singularity; there is no separate profile for it
        nextflow_profile       => 'slurm',

        # --- Directory roots (must be overridden on command line or here) ---
        nextflow_work_root     => undef,
        nextflow_output_root   => undef,

        # --- Riboseq pipeline parameters (passed as --param to nextflow run) ---

        # Required
        riboseq_sample_sheet  => undef,   # CSV with Run and study_accession columns
        riboseq_fasta         => undef,   # Reference genome FASTA

        # Strongly recommended
        riboseq_star_index    => undef,   # STAR genome index directory
        riboseq_gtf           => undef,   # GTF annotation file

        # Optional: passed only when defined (undef => not included in command)
        riboseq_rrna_index              => undef,
        riboseq_chrom_sizes_file        => undef,
        riboseq_ribometric_annotation   => undef,

        # Optional with sensible defaults matching the Nextflow pipeline.
        # Booleans must be the strings 'true'/'false' -- nf-schema 2.x rejects 0/1.
        riboseq_fetch                   => 'false',
        riboseq_force_fetch             => 'false',
        riboseq_mismatches              => 3,
        riboseq_allow_introns           => 'true',
        riboseq_max_multimappers        => 10,
        riboseq_alignment_type          => 'EndToEnd',
    };
}


sub pipeline_analyses {
    my ($self) = @_;

    # Build the nextflow --param hashref.
    # Values that remain undef are skipped in HiveRunNextflow::_build_command.
    my $nf_params = {
        # Required
        sample_sheet     => $self->o('riboseq_sample_sheet'),
        fasta            => $self->o('riboseq_fasta'),

        # Strongly recommended
        star_index       => $self->o('riboseq_star_index'),
        gtf              => $self->o('riboseq_gtf'),

        # Optional: pass undef if not set; HiveRunNextflow will skip them
        rrna_index              => $self->o('riboseq_rrna_index'),
        chrom_sizes_file        => $self->o('riboseq_chrom_sizes_file'),
        ribometric_annotation   => $self->o('riboseq_ribometric_annotation'),

        # Flags with explicit defaults
        fetch            => $self->o('riboseq_fetch'),
        force_fetch      => $self->o('riboseq_force_fetch'),
        mismatches       => $self->o('riboseq_mismatches'),
        allow_introns    => $self->o('riboseq_allow_introns'),
        max_multimappers => $self->o('riboseq_max_multimappers'),
        alignment_type   => $self->o('riboseq_alignment_type'),
    };

    return [

        # 1. Sanity-check that the essential input files exist before doing anything
        {
            -logic_name  => 'check_inputs',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                riboseq_sample_sheet => $self->o('riboseq_sample_sheet'),
                riboseq_fasta        => $self->o('riboseq_fasta'),
                cmd => 'test -f #riboseq_sample_sheet# && test -f #riboseq_fasta#',
            },
            -rc_name     => '1GB',
            -input_ids   => [{}],
            -flow_into   => { 1 => ['run_riboseq'] },
        },

        # 2. Launch the Nextflow riboseq pipeline and wait for its exit status.
        #    The eHive worker just monitors the process; real work runs in SLURM
        #    jobs submitted by Nextflow.
        {
            -logic_name      => 'run_riboseq',
            -module          => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRunNextflow',
            -parameters      => {
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
                nextflow_params        => $nf_params,
            },
            -rc_name         => 'nextflow_launcher',
            -max_retry_count => 2,
            -flow_into       => { 1 => ['pipeline_complete'] },
        },

        # 3. Placeholder analysis to mark completion in the pipeline graph.
        #    Can be extended to send notifications, update a registry, etc.
        {
            -logic_name  => 'pipeline_complete',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                nextflow_output_dir => catdir(
                    $self->o('nextflow_output_root'),
                    $self->o('nextflow_pipeline_name'),
                ),
                cmd => 'echo "Riboseq pipeline complete. Results: #nextflow_output_dir#"',
            },
            -rc_name     => '1GB',
        },

    ];
}


1;
