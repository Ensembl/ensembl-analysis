#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDumpGenome;

use strict;
use warnings;
use feature 'say';

use File::Spec::Functions;
use File::Path qw(make_path);

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(execute_with_wait);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: Default parameters
               alt_as_scaffolds => 0, # Dump the haplotypes and patches as scaffold instead of chromosome if 1
               patch_only => 0, # Only dump that patches if set to 1
               lfs_stripe => 0, # The lfs stripe value, 0 disable it, -1 for all or > 0 for fine setting
               lfs_options => undef, # If they are more settings you want to give to  lfs setstripe
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    alt_as_scaffolds => 0,
    patch_only => 0,
    lfs_stripe => 0,
    lfs_options => undef,
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Check that the paremeters are set and create the output directory if needed. If 'lfs_stripe'
              is different from 0, it stripes the directory based on the 'lfs_stripe' value.
 Returntype : None
 Exceptions : Throws if 'target_db' is not set
              Throws if 'species_name' is not set
              Throws if 'output_path' is not set
              Throws if 'coord_system_name' is not set
              Throws if 'enscode_root_dir' is not set

=cut

sub fetch_input {
  my $self = shift;

  unless($self->param('target_db')) {
    $self->throw("target_db not passed into parameters hash. The core db to load the assembly info ".
                 "into must be passed in with write access. You need to pass in the connection hash with 'target_db'");
  }

  unless($self->param('species_name')) {
    $self->throw("species_name not passed into parameters hash. You need to specify what species you're working on with 'species_name'");
  }

  if ($self->param_is_defined('output_path')) {
    my $output_path = $self->param('output_path');
    if (!-e $output_path) {
      make_path($output_path);
    }
    if ($self->param('lfs_stripe')) {
      my $lfs_command = 'lfs setstripe -c '.$self->param('lfs_stripe');
      $lfs_command .= ' '.$self->param('lfs_options') if ($self->param('lfs_options'));
      execute_with_wait("$lfs_command $output_path");
    }
  }
  else {
    $self->throw("Output path not passed into parameters hash. You need to specify where the output dir will be with 'output_path'");
  }

  unless($self->param('coord_system_name')) {
    $self->throw("Coord system name was not passed in using the 'coord_system_version' parameter");
  }

  unless($self->param('enscode_root_dir')) {
    $self->throw("enscode_dir not passed into parameters hash. You need to specify where your code checkout is");
  }

  return 1;
}


=head2 run

 Arg [1]    : None
 Description: Create the command line for the dump script.
              It add the masking options if 'repeat_logic_names' has values in its arrayref
              Uses 'alt_as_scaffolds' to know how to dump hapoltypes and patches if needed
              Uses 'patch_only' to know which sequence to dump
 Returntype : None
 Exceptions : None

=cut

sub run {
  my $self = shift;

  my $db_info = $self->param('target_db');
  my $repeat_logic_names = $self->param('repeat_logic_names');
  my $repeat_string;
  if(scalar(@{$repeat_logic_names})) {
    $repeat_string .= ' -mask -softmask ';
    foreach my $repeat_logic_name (@{$repeat_logic_names}) {
      $repeat_string .= ' -mask_repeat '.$repeat_logic_name;
    }
  }

  my $cmd = 'perl '.catfile($self->param('enscode_root_dir'), 'ensembl-analysis', 'scripts', 'sequence_dump.pl').
            ' -dbuser '.$db_info->{-user}.
            ' -dbport '.$db_info->{-port}.
            ' -dbhost '.$db_info->{-host}.
            ' -dbname '.$db_info->{-dbname}.
            ' -coord_system_name '.$self->param('coord_system_name').
            ' -toplevel'.
            ' -onefile'.
            ' -nonref'.
            ' -filename '.catfile($self->param('output_path'), $self->param('species_name').'_softmasked_toplevel.fa');

  $cmd .= $repeat_string if ($repeat_string);
  $cmd .= ' -dbpass '.$db_info->{-pass} if ($db_info->{-pass});
  $cmd .= ' -alt_as_scaffolds ' if ($self->param('alt_as_scaffolds'));

  if($self->param('patch_only')) {
    $cmd .= " -patch_only ";
  }

  say "Running command:\n".$cmd;

  execute_with_wait($cmd);
}


=head2 write_output

 Arg [1]    : None
 Description: Nothing to do store in the database
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my $self = shift;

  return 1;
}

1;
