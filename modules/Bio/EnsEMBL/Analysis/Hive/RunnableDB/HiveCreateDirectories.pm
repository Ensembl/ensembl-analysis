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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDirectories

=head1 SYNOPSIS

To create two simple directories
{
  -logic_name => 'pipeline_directories',
  -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDirectories',
  -parameters => {
    paths => [
      'input_dir',
      'output_dir',
    ],
  }
  -rc_name => 'default',
}

To create a simple directory and one directory with striping and different permissions
{
  -logic_name => 'pipeline_dir',
  -module => 'Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDirectories',
  -parameters => {
    paths => [
      'input_dir',
      {
        path => 'output_dir',
        stripe => 1,
        mode => 0755,
      },
    ],
  }
  -rc_name => 'default',
}

=head1 DESCRIPTION

The module will create directories based on the parameters given.
This is most useful when you are using LFS and you need to have
some of the directories striped. But also if you need to share
directories between people member of the same UNIX group.
The list of path to create should be given using the 'paths'
parameter.
If you need to stripe directories or to set different permissions,
you need to provide a hash instead of a string

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCreateDirectories;

use strict;
use warnings;

use File::Spec::Functions qw(catdir splitdir);


use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_production_directory);

use parent qw(Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB);


=head2 fetch_input

 Arg [1]    : None
 Description: Checks that 'paths' has been set
 Returntype : None
 Exceptions : None

=cut

sub fetch_input {
  my ($self) = @_;

  $self->param_required('paths');
}


=head2 run

 Arg [1]    : None
 Description: Create the commands to run to create all needed directories
              It will recursively generate the commands needed if you need
              multiple directories created
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;

  my @paths_to_create;
  foreach my $path (@{$self->param('paths')}) {
    my $dir;
    my $do_stripe;
    my $mode;
    if (ref($path) eq 'HASH') {
      $dir = $path->{path};
      $do_stripe = $path->{stripe};
      $mode = $path->{mode};
    }
    else {
      $dir = $path;
    }
    my $full_path;
    foreach my $directory (splitdir($dir)) {
      next unless $directory;
      $full_path = catdir($full_path, $directory);
      next if (-d $full_path);
      push(@paths_to_create, [$full_path, $do_stripe, $mode]);
      $self->say_with_header("About to create '$full_path'"
        .($do_stripe ? ", stripe: '$do_stripe'" : '')
        .($mode ? ", permissions: '".sprintf("%04o", $mode)."'" : ""));
    }
  }
  $self->output(\@paths_to_create);
}


=head2 write_output

 Arg [1]    : None
 Description: Create the directories, changes the permissions if needed
              and stripe the directories if needed and using Lfs
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;

  foreach my $path (@{$self->output}) {
    create_production_directory(@$path);
  }
}

1;
