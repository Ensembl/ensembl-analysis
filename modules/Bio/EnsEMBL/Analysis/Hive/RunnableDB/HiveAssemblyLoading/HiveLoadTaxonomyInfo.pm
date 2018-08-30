=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadTaxonomyInfo

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadTaxonomyInfo;

use strict;
use warnings;

use File::Spec::Functions qw(catfile);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: Returns the defaults parameters
               ensembl_release => $ENV{ENSEMBL_RELEASE},
 Returntype : None
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    ensembl_release => $ENV{ENSEMBL_RELEASE},
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Prepare the command to run the populate_species_meta.pl script. 'enscode_root_dir', 'target_db'
              are required. You can also specify the 'production_db' and the 'taxonomy_db' if needed.
              The 'ensembl_release' defaults to $ENV{ENSEMBL_RELEASE}
 Returntype : None
 Exceptions : Throws if 'enscode_root_dir' is not set
              Throws if 'target_db' is not set

=cut

sub fetch_input {
  my $self = shift;

  if ($self->param('skip_analysis')) {
    $self->complete_early('Skipping the analysis');
  }
  my $taxonomy_script = catfile($self->param_required('enscode_root_dir'), 'ensembl-production', 'scripts', 'production_database', 'populate_species_meta.pl');
  my $target_db = $self->param_required('target_db');
  my $cmd = 'perl '.$taxonomy_script.
            ' -h '.$target_db->{'-host'}.
            ' -u '.$target_db->{'-user'}.
            ' -P '.$target_db->{'-port'}.
            ' -d '.$target_db->{'-dbname'}.
            ' -rel '.$self->param('ensembl_release');
  $cmd .= ' -p '.$target_db->{'-pass'} if (defined $target_db->{'-pass'});

  if ($self->param_is_defined('production_db')) {
    my $production_db = $self->param('production_db');
    $cmd .= ' -mh '.$production_db->{-host}.
            ' -mu '.$production_db->{-user}.
            ' -md '.$production_db->{-dbname}.
            ' -mP '.$production_db->{-port};
    $cmd .= ' -mp '.$production_db->{-pass} if (defined $production_db->{-pass});
  }
  if ($self->param_is_defined('taxonomy_db')) {
    my $taxonomy_db = $self->param('taxonomy_db');
    $cmd .= ' -th '.$taxonomy_db->{-host}.
            ' -tu '.$taxonomy_db->{-user}.
            ' -td '.$taxonomy_db->{-dbname}.
            ' -tP '.$taxonomy_db->{-port};
    $cmd .= ' -tp '.$taxonomy_db->{-pass} if (defined $taxonomy_db->{-pass});
  }

  $self->param('cmd', $cmd);
}


=head2 run

 Arg [1]    : None
 Description: Execute the command prepare in fetch_input
 Returntype : None
 Exceptions : Throws if the command failed

=cut

sub run {
  my $self = shift;

  if(system($self->param('cmd'))) {
    $self->throw("The load_taxonomy script returned a non-zero exit value. Commandline used:\n".$self->param('cmd'));
  }
}


=head2 write_output

 Arg [1]    : None
 Description: Do nothing, can probably be removed
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my $self = shift;

  return 1;
}

1;
