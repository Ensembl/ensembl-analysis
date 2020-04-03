package PseudoFactory;

use warnings;
use strict;
use feature 'say';

sub param_defaults {
  my ($self) = @_;

  return {
	  %{$self->SUPER::param_defaults},
	 };
}

sub run {
  my ($self) = @_;
  my $reg = 'Bio::EnsEMBL::Registry';
  if ($self->param_is_defined('registry_file')) {
    $reg->load_all($self->param('registry_file'));
  }

  my $all_dbas = $reg->get_all_DBAdaptors();
  my %dbs;

  if ( ! scalar(@$all_dbas) ) ) {
    $self->throw("No databases found in the registry");
  }

  foreach my $dba (@$all_dbas) {
    $$dbs{$dba->dbc->dbname}{$dba->species} = $dba;
    $dba->dbc->disconnect_if_idle();
  }

  $self->param( 'dbs', \%dbs );
}

sub write_output {
  my ($self) = @_;

  


}

1;
