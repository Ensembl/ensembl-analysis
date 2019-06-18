=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::JiraTicket

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::JiraTicket;

use strict;
use warnings;

use MIME::Base64 qw(encode_base64url);
use LWP::UserAgent;
use JSON;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: Default parameters:
               jira_base_url => 'https://www.ebi.ac.uk/panda/jira/rest/api/2',
               authentification_type => 'basic', # can be none, basic, cookie, oauth
               type => 'Bug',
               project_name => 'Ensembl Genebuild',
               project_key => 'ENSGENEBUI',
               assignee => 'any', # can be self, unassigned, any
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    jira_base_url => 'https://www.ebi.ac.uk/panda/jira/rest/api/2',
    authentification_type => 'basic', # can be none, basic, cookie, oauth
    type => 'Bug',
    project_name => 'Ensembl Genebuild',
    project_key => 'ENSGENEBUI',
    assignee => 'any', # can be self, unassigned, any
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: If tickect_name is not set, it completes early. Otherwise it checks
              the 'cmd' parameter and the authentification type. If the type is a
              sub task, it checks that 'parent' is present. It also makes sure that
              the 'jira_base_url' ends with "/"
 Returntype : None
 Exceptions : Throws if 'cmd' is not supported
              Throws if 'authentification_type' is not supported
              Throws if 'parent' is not set for a sub task

=cut

sub fetch_input {
  my ($self) = @_;

  if ($self->param_is_defined('ticket_name')) {
    my $ticket_name = $self->param('ticket_name');
    if ($ticket_name) {
      my $cmd = $self->param_required('cmd');
      $self->param_required('parent')
        if ($cmd eq 'create' and $self->param('type') eq 'Sub-task');
      $self->throw("command $cmd not supported") unless ($cmd eq 'create' or $cmd eq 'comment');
      my $auth_type = $self->param('authentification_type');
      if ($auth_type eq 'none') {

      }
      elsif ($auth_type eq 'basic') {
        if ($self->param_is_defined('base64')) {
          $self->throw('You authentification is empty') unless ($self->param('base64'));
        }
        else {
          $self->param('base64', encode_base64url($self->param_required('user').':'.$self->param_required('password')));
        }
      }
      elsif ($auth_type eq 'cookie') {
        $self->throw("'authentification_type' '$auth_type' is not available");
      }
      elsif ($auth_type eq 'oauth') {
        $self->throw("'authentification_type' '$auth_type' is not available");
      }
      else {
        $self->throw("'authentification_type' '$auth_type' is not supported");
      }
    }
    else {
      $self->complete_early('"ticket_name" is empty, not using JIRA');
    }
    $self->param('jira_base_url', $self->param('jira_base_url').'/') unless ($self->param('jira_base_url') =~ '/$');
  }
  else {
    $self->complete_early('"ticket_name" does not exist.');
  }
}


=head2 run

 Arg [1]    : None
 Description: It creates the JSON string and the REST URL to create, update or
              comment a JIRA ticket. If 'params' is set, it adds the hash to the
              hash to be changed to JSON.
 Returntype : None
 Exceptions : Throws if an HTTP request fails
              Throws if the number of returned element is not the expected number, (1)

=cut

sub run {
  my ($self) = @_;

  my $ua = LWP::UserAgent->new();
  my $json_coder = JSON->new;
  $ua->default_header('Content-Type' => 'application/json');
  if ($self->param_is_defined('base64')) {
    $ua->default_header('Authorization' => 'Basic '.$self->param('base64'));
  }
  $self->param('user_agent', $ua);
  my $cmd = $self->param('cmd');
  my @output;
  if ($cmd eq 'create') {
    my %json_hash = (
      fields => {
        project => {
          key => $self->param('project_key'),
        },
        issuetype => {
          name => $self->param('type'),
        },
        summary => $self->param('ticket_name'),
      },
    );
    if ($self->param('type') eq 'Sub-task') {
      my %json_search  = (
        jql => 'project="'.$self->param('project_name').'" AND summary~"'.$self->param('parent').'"',  
        fields => ['id'],
      );
      my $response = $ua->post($self->param('jira_base_url').'search', 'Content-Type' => 'application/json', 'Content' => $json_coder->encode(\%json_search));
      if ($response->is_success) {
        my $data = $json_coder->decode($response->decoded_content);
        if ($data->{total}) {
          if ($data->{total} == 1) {
            $json_hash{fields}->{parent} = { id => $data->{issues}->[0]->{id}};
          }
          else {
            $self->throw('You have more than one result for '.$response->decoded_content);
          }
        }
        else {
          $self->throw('Something went wrong for '.$response->decoded_content);
        }
      }
      else {
        $self->throw('Could not find parent ticket: '.$self->param('parent').' '.$response->status_line);
      }
    }
    if ($self->param('assignee') eq 'self') {
      my $response = $ua->get($self->param('jira_base_url').'myself');
      if ($response->is_success) {
        my $data = $json_coder->decode($response->decoded_content);
        $json_hash{fields}->{assignee}->{name} = $data->{name};
      }
      else {
        $self->throw($response->status_line);
      }
    }
    elsif ($self->param('assignee') eq 'unassigned') {
        $json_hash{fields}->{assignee}->{name} = '';
    }
    if ($self->param_is_defined('params')) {
      %{$json_hash{fields}} = (%{$json_hash{fields}}, %{$self->param('params')});
    }
    push(@output, [$self->param('jira_base_url').'issue/', $json_coder->encode(\%json_hash)]);
  }
  elsif ($cmd eq 'comment') {
    my %json_hash  = (
      jql => 'project="'.$self->param('project_name').'" AND summary~"'.$self->param('ticket_name').'"',  
      fields => ['id', 'status'],
    );
    if ($self->param('type') eq 'Sub-task') {
      push(@{$json_hash{fields}}, 'parent');
    }
    my $response = $ua->post($self->param('jira_base_url').'search', 'Content-Type' => 'application/json', 'Content' => $json_coder->encode(\%json_hash));
    if ($response->is_success) {
      my $data = $json_coder->decode($response->decoded_content);
      if ($data->{total}) {
        if ($data->{total} == 1) {
          if ($data->{issues}->[0]->{fields}->{status}->{name} eq 'Reported') {
            my %update = (
              transition => {
                id => 4,
              },
            );
            if (exists $data->{issues}->[0]->{fields}->{parent} and $data->{issues}->[0]->{fields}->{parent}->{fields}->{status}->{name} eq 'Reported') {
              push(@output, [$self->param('jira_base_url').'issue/'.$data->{issues}->[0]->{fields}->{parent}->{key}.'/transitions', $json_coder->encode(\%update)]);
            }
            push(@output, [$self->param('jira_base_url').'issue/'.$data->{issues}->[0]->{key}.'/transitions', $json_coder->encode(\%update)]);
          }
          my %comment = (
            body => $self->param('comment'),
          );
          push(@output, [$self->param('jira_base_url').'issue/'.$data->{issues}->[0]->{key}.'/comment/', $json_coder->encode(\%comment)]);
        }
        else {
          $self->throw('You have more than one result for '.$response->decoded_content);
        }
      }
      else {
        $self->throw('Something went wrong for '.$response->decoded_content);
      }
    }
    else {
      $self->throw('Could not query JIRA for '.$json_hash{jql}."\n".$response->status_line);
    }
  }
  $self->output(\@output);
}


=head2 write_output

 Arg [1]    : None
 Description: Send the POST request to create, modify,... a JIRA ticket
 Returntype : None
 Exceptions : Throws if the HTTP request failed

=cut

sub write_output {
  my ($self) = @_;

  my $ua = $self->param('user_agent');
  foreach my $query (@{$self->output}) {
    my $response = $ua->post($query->[0], 'Content-Type' => 'application/json', 'Content' => $query->[1]);
    if ($response->is_success) {
      $self->warning($response->status_line);
    }
    else {
      $self->throw('Could not process '.$query->[0]."\n".$query->[1]."\n".$response->status_line);
    }
  }
}

1;
