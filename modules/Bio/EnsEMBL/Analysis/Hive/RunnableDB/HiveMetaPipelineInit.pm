=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMetaPipelineInit

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMetaPipelineInit;

use strict;
use warnings;

use File::Spec::Functions qw(catfile);
use Bio::EnsEMBL::Hive::Utils qw(stringify);

use parent ('Bio::EnsEMBL::Hive::RunnableDB::SystemCmd');


=head2 param_defaults

 Arg [1]    : None
 Description: Set 'hive_init_script' to use 'enscode_root_dir'
               hive_init_script => catfile('#enscode_root_dir#', 'ensembl-hive', 'scripts', 'init_pipeline.pl'),
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    hive_init_script => catfile('#enscode_root_dir#', 'ensembl-hive', 'scripts', 'init_pipeline.pl'),
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Create the command for initialising the sub pipeline. You need to provide
              'ehive_url', 'meta_pipeline_db' and 'hive_config'. If you want to specify
              databases you need to add the name of the database(s) in 'databases'. The
              parameters will be created by using "${db_name}_name" where $db_name is the
              name of the databases in 'hive_config'. You should also have paramters in
              the config like "${db_name}_name", "${db_name}_port", "${db_name}_host",
              "${db_name}_user", "${db_name}_password".
 Returntype : None
 Exceptions : Throws if 'ehive_url', 'pipeline_name', 'meta_pipeline_db' or 'hive_config' is not set
              Throws if 'hive_init_script' does not exists

=cut

sub fetch_input {
  my $self = shift;

  my $url = $self->param_required('ehive_url');
  my $pipeline_db = $self->param_required('meta_pipeline_db');
  $self->throw('Your pipeline dbname should be all lowercase, not '.$pipeline_db->{'-dbname'})
    if ($pipeline_db->{'-dbname'} =~ /[[:upper:]]/);
  $self->throw('"hive_init_script" does not exist: '.$self->param('hive_init_script'))
    unless (-e $self->param('hive_init_script'));
  my $dna_db;
  if ($self->param_is_defined('dna_db')) {
    $dna_db = $self->param_required('dna_db');
  }
  else {
    $self->warning('No "dna_db" was defined');
  }
  my @cmd;
  push(@cmd, 'perl', $self->param('hive_init_script'), $self->param_required('hive_config'),
    '-host', $pipeline_db->{'-host'},
    '-port', $pipeline_db->{'-port'},
    '-user', $pipeline_db->{'-user'},
    '-password', $pipeline_db->{'-pass'},
    '-pipe_db_host', $pipeline_db->{'-host'},
    '-pipe_db_port', $pipeline_db->{'-port'},
    '-pipe_db_user', $pipeline_db->{'-user'},
    '-pipe_db_pass', $pipeline_db->{'-pass'},
    '-pipe_db_name', $pipeline_db->{'-dbname'},
    '-pipeline_name', $self->param_required('pipeline_name'));
  push(@cmd, '-enscode_root_dir', $self->param('enscode_root_dir'));
  if ($dna_db) {
    $self->warning('Your dna dbname has upper case character, it might cause problems, '.$dna_db->{'-dbname'})
      if ($dna_db->{'-dbname'} =~ /[[:upper:]]/);
    push(@cmd, '-dna_db_name', $dna_db->{'-dbname'});
    push(@cmd, '-dna_db_host', $dna_db->{'-host'});
    push(@cmd, '-dna_db_port', $dna_db->{'-port'});
    push(@cmd, '-dna_db_user', $dna_db->{'-user'}) if (exists $dna_db->{'-user'});
    push(@cmd, '-dna_db_pass', $dna_db->{'-pass'}) if (exists $dna_db->{'-pass'} and $dna_db->{'-pass'});
  }
  foreach my $db_title (@{$self->param('databases')}) {
    my $db = $self->param($db_title);
    $self->warning("Your $db_title dbname has upper case character, it might cause problems, ".$db->{'-dbname'})
      if ($db->{'-dbname'} =~ /[[:upper:]]/);
    push(@cmd, '-'.$db_title.'_name', $db->{'-dbname'});
    push(@cmd, '-'.$db_title.'_host', $db->{'-host'});
    push(@cmd, '-'.$db_title.'_port', $db->{'-port'});
    push(@cmd, '-'.$db_title.'_user', $db->{'-user'}) if (exists $db->{'-user'});
    push(@cmd, '-'.$db_title.'_pass', $db->{'-pass'}) if (exists $db->{'-pass'} and $db->{'-pass'});
  }
  if ($self->param_is_defined('extra_parameters')) {
    my $extra_parameters = $self->param('extra_parameters');
    if (ref($extra_parameters) eq 'HASH') {
      foreach my $key (keys %$extra_parameters) {
        if (ref($extra_parameters->{$key}) eq 'ARRAY') {
          # We need to make sure that an arrayref/hashref is correctly passed to the init script
          foreach my $arraydata (@{$extra_parameters->{$key}}) {
            my $value = ref($arraydata) ? stringify($arraydata) : $arraydata;
            push(@cmd, "-$key", $value);
          }
        }
        else {
          # We need to make sure that an arrayref/hashref is correctly passed to the init script
          my $value = ref($extra_parameters->{$key}) ? stringify($extra_parameters->{$key}) : $extra_parameters->{$key};
          push(@cmd, "-$key", $value) if (defined($value));
        }
      }
    }
  }
  if ($self->param_is_defined('commandline_params')) {
    push(@cmd, $self->param('commandline_params'));
  }
  $self->param('cmd', join(' ', @cmd));
}


1;
