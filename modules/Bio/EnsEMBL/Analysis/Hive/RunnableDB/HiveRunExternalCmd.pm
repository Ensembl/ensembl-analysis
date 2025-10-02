package Bio::EnsEMBL::Hive::RunnableDB::HiveRunExternalCmd;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::Process');
use Bio::EnsEMBL::Compara::Utils::RunCommand;

sub param_defaults {
    return {
        'cmd' => undef,   # command to run
    };
}

sub run {
    my $self = shift;

    my $cmd = $self->param_required('cmd');

    my $rc = Bio::EnsEMBL::Compara::Utils::RunCommand
                ->new_and_exec($cmd, { die_on_failure => 1 });

    # Save stdout into a hive param so it can be used in flow_into
    my $stdout = $rc->out;
    chomp $stdout;
    $self->param('stdout', $stdout);
    $self->param('stderr', $rc->err);
}

sub write_output {
    my $self = shift;

    # Flow stdout as #stdout#
    $self->dataflow_output_id({ stdout => $self->param('stdout') }, 1);
}

1;
