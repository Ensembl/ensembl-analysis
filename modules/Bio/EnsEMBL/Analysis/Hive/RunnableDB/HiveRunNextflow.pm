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

Runnable that launches a Nextflow pipeline and waits for its exit status.
The eHive job is a lightweight launcher that submits work to Nextflow
(which in turn submits to SLURM/LSF), so the worker resource class should
be small but have a wall-time long enough to outlast the full pipeline run.

On failure the runnable parses the Nextflow log to extract the actual error
reason and includes it in the eHive job message, so you don't have to dig
through the raw log.

=head2 Resume modes

  never   - Always run fresh; -resume is never passed. A new work directory
            is created per eHive retry attempt.

  attempt - Pass -resume; work directory is stable within a single eHive
            retry attempt (indexed by retry_count). On eHive retry the
            retry_count increments, yielding a fresh work directory.
            This is the recommended default.

  job     - Pass -resume; work directory is stable across ALL retries of
            the same eHive job (retry_count not in the path).

=head2 Directory layout produced

  ${nextflow_work_root}/
    ${nextflow_pipeline_name}/
      job_${hive_job_id}/
        attempt_${retry_count}/     # resume=never or resume=attempt
          work/                     # Nextflow workDir
          nextflow.log              # Nextflow log (parsed on failure)
          manifest.json             # run metadata
        # OR for resume=job:
        work/
        nextflow.log
        manifest.json

=head1 PARAMETERS

  nextflow_pipeline_dir  - Path to directory containing main.nf (required)
  nextflow_pipeline_name - Used in work directory naming (required)
  nextflow_work_root     - Base directory for per-job work dirs (required)
  nextflow_output_dir    - Passed as --outdir to Nextflow (optional)
  nextflow_binary        - Path/name of nextflow binary (default: nextflow)
  nextflow_resume_mode   - never | attempt | job (default: attempt)
  nextflow_profile       - Nextflow -profile value, e.g. 'slurm,singularity'
  nextflow_params        - Hashref of --param => value pairs for the pipeline
  nextflow_extra_flags   - Arrayref of extra CLI flags passed verbatim

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRunNextflow;

use strict;
use warnings;
use feature 'say';

use File::Path qw(make_path);
use File::Spec::Functions qw(catdir catfile);
use POSIX qw(strftime);

use base ('Bio::EnsEMBL::Hive::Process');


sub param_defaults {
    return {
        nextflow_binary       => 'nextflow',
        nextflow_resume_mode  => 'attempt',   # never | attempt | job
        nextflow_profile      => undef,
        nextflow_params       => {},           # hashref --key => value
        nextflow_extra_flags  => [],           # extra flags passed verbatim
        nextflow_output_dir   => undef,        # passed as --outdir
    };
}


sub run {
    my ($self) = @_;

    my $job_id  = $self->input_job->dbID;
    my $attempt = $self->input_job->retry_count;

    my $work_dir = $self->_compute_work_dir($job_id, $attempt);
    (my $run_dir = $work_dir) =~ s{/work$}{};
    make_path($work_dir);   # also creates run_dir

    my $log_file = catfile($run_dir, 'nextflow.log');
    my ($cmd, $use_resume) = $self->_build_command($work_dir, $log_file);

    $self->_write_manifest($run_dir, $cmd, $use_resume, $attempt, 'running');

    $self->warning("HiveRunNextflow: launching pipeline");
    $self->warning("  workDir:     $work_dir");
    $self->warning("  log:         $log_file");
    $self->warning("  resume mode: " . $self->param('nextflow_resume_mode'));
    $self->warning("  command:     $cmd");

    my $exit_status = system($cmd);

    my $final_status = ($exit_status == 0) ? 'complete' : 'failed';
    $self->_write_manifest($run_dir, $cmd, $use_resume, $attempt, $final_status, $exit_status >> 8);

    if ($exit_status != 0) {
        my $nf_error = $self->_extract_nf_error($log_file);
        $self->throw(sprintf(
            "Nextflow exited with status %d.\n%sResume mode: %s\nLog: %s\nworkDir: %s",
            $exit_status >> 8,
            $nf_error ? "Nextflow error: $nf_error\n" : '',
            $self->param('nextflow_resume_mode'),
            $log_file,
            $work_dir,
        ));
    }
}


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

sub _compute_work_dir {
    my ($self, $job_id, $attempt) = @_;

    my $mode = $self->param('nextflow_resume_mode');
    my $base = catdir(
        $self->param_required('nextflow_work_root'),
        $self->param_required('nextflow_pipeline_name'),
        sprintf('job_%d', $job_id),
    );

    if ($mode eq 'job') {
        return catdir($base, 'work');
    }
    else {
        return catdir($base, sprintf('attempt_%d', $attempt), 'work');
    }
}


sub _build_command {
    my ($self, $work_dir, $log_file) = @_;

    my $mode = $self->param('nextflow_resume_mode');
    my $use_resume = ($mode eq 'never') ? 0 : 1;

    # -log is a global Nextflow option and must come before the subcommand
    my @parts = (
        $self->param('nextflow_binary'),
        '-log', $log_file,
        'run',
        $self->param_required('nextflow_pipeline_dir'),
        '-work-dir', $work_dir,
    );

    push @parts, '-resume' if $use_resume;

    if (my $profile = $self->param('nextflow_profile')) {
        push @parts, '-profile', $profile;
    }

    if (my $outdir = $self->param('nextflow_output_dir')) {
        push @parts, '--outdir', $outdir;
    }

    my $params = $self->param('nextflow_params') // {};
    for my $key (sort keys %$params) {
        my $val = $params->{$key};
        next unless defined $val;
        push @parts, "--${key}", $val;
    }

    push @parts, @{ $self->param('nextflow_extra_flags') };

    return (join(' ', @parts), $use_resume);
}


# Parse the Nextflow log and return the human-readable error reason.
# Looks for the patterns Nextflow uses to report pipeline-level failures.
sub _extract_nf_error {
    my ($self, $log_file) = @_;

    return '' unless -f $log_file;

    open my $fh, '<', $log_file or return '';
    my @lines = <$fh>;
    close $fh;
    chomp @lines;

    my @errors;

    for my $line (@lines) {
        # "Session aborted -- Cause: <reason>" is the most useful single line
        if ($line =~ /Session aborted\s*--\s*Cause:\s*(.+)/) {
            push @errors, $1;
        }
        # Named exception with message on same line
        elsif ($line =~ /(?:WorkflowScriptErrorException|AbortOperationException|SchemaValidationException): (.+)/) {
            push @errors, $1;
        }
        # Generic ERROR lines from the Nextflow launcher
        elsif ($line =~ /ERROR\s+nextflow\S+\s+-\s+(.+)/) {
            push @errors, $1;
        }
    }

    # Deduplicate while preserving order
    my %seen;
    my @unique = grep { !$seen{$_}++ } @errors;

    return join(' | ', @unique);
}


sub _write_manifest {
    my ($self, $run_dir, $cmd, $use_resume, $attempt, $status, $exit_code) = @_;

    make_path($run_dir);
    my $path = catfile($run_dir, 'manifest.json');

    my $escaped_cmd = $cmd;
    $escaped_cmd =~ s/\\/\\\\/g;
    $escaped_cmd =~ s/"/\\"/g;

    my $ts = strftime('%Y-%m-%dT%H:%M:%SZ', gmtime);

    my @lines = (
        '{',
        sprintf('  "status": "%s",',       $status),
        sprintf('  "timestamp": "%s",',    $ts),
        sprintf('  "hive_job_id": %d,',    $self->input_job->dbID),
        sprintf('  "retry_count": %d,',    $attempt),
        sprintf('  "resume_mode": "%s",',  $self->param('nextflow_resume_mode')),
        sprintf('  "resume_flag": %s,',    $use_resume ? 'true' : 'false'),
        sprintf('  "work_dir": "%s/work",', $run_dir),
        sprintf('  "pipeline_dir": "%s",', $self->param('nextflow_pipeline_dir')),
        sprintf('  "command": "%s"',       $escaped_cmd),
    );

    if (defined $exit_code) {
        $lines[-1] .= ',';
        push @lines, sprintf('  "exit_code": %d', $exit_code);
    }

    push @lines, '}';

    if (open my $fh, '>', $path) {
        say $fh join("\n", @lines);
        close $fh;
    }
    else {
        $self->warning("HiveRunNextflow: could not write manifest to $path: $!");
    }
}


1;
