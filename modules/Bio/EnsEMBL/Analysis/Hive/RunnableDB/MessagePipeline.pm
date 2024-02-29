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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::MessagePipeline

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::MessagePipeline;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::Utils qw(stringify);
use parent 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB';


sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    tweak_script => 'tweak_pipeline.pl',
  }
}

sub fetch_input {
  my ($self) = @_;

  my $messages = $self->param_required('messages');
  my @messages_to_pass;
  if (ref($messages) eq 'HASH') {
    $self->throw('Malformed message, expects url, logic_name and data, not: '.join(' ', keys %$messages))
      unless (exists $messages->{url} and exists $messages->{logic_name} and exists $messages->{data});
    if (defined $messages->{url}) {
      push(@messages_to_pass, $messages);
    }
  }
  elsif (ref($messages) eq 'ARRAY') {
    foreach my $message (@$messages) {
      if (ref($message) eq 'HASH') {
        $self->throw('Malformed message, expects url, logic_name and data, not: '.join(' ', keys %$message))
          unless (exists $message->{url} and exists $message->{logic_name} and exists $message->{data});
        if (defined $message->{url}) {
          push(@messages_to_pass, $message);
        }
      }
      else {
        $self->throw("Malformed message, expects a hash with at least url, logic_name and data not: $message");
      }
    }
  }
  else {
    $self->throw("Malformed messages, expects a hash or array of hashes with at least url, logic_name and data not: $messages");
  }
  if (@messages_to_pass) {
    $self->param('messages_to_pass', \@messages_to_pass);
  }
  else {
    $self->complete_early('No messages to pass, check url of your messages if you expected some');
  }
}


sub run {
  my ($self) = @_;

  foreach my $message (@{$self->param('messages_to_pass')}) {
    my $data;
    my $base_cmd = $self->param('tweak_script').' --url '.$message->{url};
    foreach my $key ('reg_conf', 'reg_type', 'reg_alias') {
      if (exists $message->{$key} and $message->{$key}) {
        $base_cmd .= " --$key ".$message->{$key}
      }
    }
    $self->say_with_header($base_cmd);
    my $analysis_pattern = 'analysis['.$message->{logic_name}.'].param['.$message->{param}.']';
    if (exists $message->{update} and $message->{update}) {
      my $check_command = "$base_cmd --SHOW $analysis_pattern";
      $self->say_with_header($check_command);
      open(IO, "$check_command |") || $self->throw("Could not execute '$check_command'");
      while (my $line = <IO>) {
        if ($line =~ /Tweak.Error/) {
          $self->throw("Failed to process the command: '$line'");
        }
        elsif ($line =~ /Tweak.Show\s+\S+\s+::\t(.*)/) {
          if ($1 ne '(missing value)') {
            $data = eval($1);
            $self->say_with_header($data);
          }
          else {
            $self->say_with_header('Value does not exist as a parameter');
          }
        }
      }
      close(IO) || $self->throw("'$check_command' did not exited cleanly");
      if (ref($data) eq 'ARRAY') {
        push(@$data, $message->{data});
      }
      elsif(ref($data) eq 'HASH') {
        %$data = %{$message->{data}};
      }
    }
    else {
      $data = $message->{data};
    }
    $self->output(["$base_cmd --tweak '$analysis_pattern=".stringify($data)."'"]);
  }
}


sub write_output {
  my ($self) = @_;

  foreach my $cmd (@{$self->output}) {
    $self->say_with_header($cmd);
    my (undef, undef, undef, $stdout) = $self->run_system_command($cmd);
    if ($stdout =~ /Tweak.Error/) {
      $self->throw("Failed to execute '$cmd'");
    }
  }
}

1;
