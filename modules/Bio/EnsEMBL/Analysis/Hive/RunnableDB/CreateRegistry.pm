=head1 LICENSE

Copyright [2019-2020] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::CreateRegistry

=head1 SYNOPSIS


=head1 DESCRIPTION

Create registry file for use in different parts of the pipeline (LastZ, Loading analysis descriptions from Produciton db)

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::CreateRegistry;

use warnings;
use strict;
use feature 'say';

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

=head2 fetch_input

 Arg [1]    : None
 Description: Set db connections, the datbases must exist.
              If you provide the parameter source_dbs, it will use the hashes
              in this arrayref and the production_db to generate a registry file.
              Otherwise it needs:
                compara_db
                target_db
                production_db
                projection_source_db
 Returntype : None
 Exceptions : Throws if source_dbs is set but empty

=cut

sub fetch_input {
  my ($self) = @_;

  $self->param_required('registry_file');
  if ($self->param_is_defined('source_dbs')) {
    my @source_dbs;
    foreach my $hashref (@{$self->param('source_dbs')}) {
      push(@source_dbs, $self->hrdb_get_dba($hashref));
    }
    $self->throw("source_dbs is empty") unless (@source_dbs);
    $self->param('databases', \@source_dbs);
  }
  else {
    my $compara_dba = $self->hrdb_get_dba($self->param_required('compara_db'),undef,'Compara');
    $self->hrdb_set_con($compara_dba,'compara_db');

    my $target_dba = $self->hrdb_get_dba($self->param_required('target_db'));
    $self->hrdb_set_con($target_dba,'target_db');

    my $projection_source_dba = $self->hrdb_get_dba($self->param_required('projection_source_db'));
    $self->hrdb_set_con($projection_source_dba,'projection_source_db');
  }

  my $production_dba = $self->hrdb_get_dba($self->param_required('production_db'),undef,'Production');
  $self->hrdb_set_con($production_dba,'production_db');

}

=head2 run

 Arg [1]    : None
 Description: Create registry file, Databases.pm
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;

  open(OUT,">".$self->param('registry_file')) || $self->throw('Could not open '.$self->param('registry_file'));
  if ($self->param_is_defined('databases')) {
    say OUT 'package Reg;';
    say OUT 'use warnings;';
    say OUT 'use strict;';
    say OUT 'use Bio::EnsEMBL::DBSQL::DBAdaptor;';
    say OUT 'use Bio::EnsEMBL::Production::DBSQL::DBAdaptor;';
    say OUT '{';
    foreach my $db (@{$self->param('databases')}) {
      # The group is important for datachecks
      my ($group) = $db->dbc->dbname =~ /(otherfeatures|rnaseq|cdna)/;
      # Fancy way of saying if $group is not defined set it to $db->group
      $group //= $db->group;
      say OUT ref($db), '->new(';
      foreach my $key (qw(host port dbname user pass)) {
        if ($db->dbc->$key) {
          say OUT "-$key => '", $db->dbc->$key, "',";
        }
      }
      say OUT "-species => '", $db->get_MetaContainer->get_production_name, "',";
      say OUT "-group => '", $group, "',";
      say OUT ');';
    }
    my $db = $self->hrdb_get_con('production_db');
    say OUT ref($db), '->new(';
    foreach my $key (qw(host port dbname user pass)) {
      if ($db->dbc->$key) {
        say OUT "-$key => '", $db->dbc->$key, "',";
      }
    }
    say OUT "-species => 'multi',";
    say OUT "-group => 'production',";
    say OUT ');';
    say OUT '}';
  }
  else {
    my $compara_dba =  $self->hrdb_get_con('compara_db');
    my $target_dba = $self->hrdb_get_con('target_db');
    my $projection_source_dba = $self->hrdb_get_con('projection_source_db');
    my $production_dba = $self->hrdb_get_con('production_db');

    my $source_production_name = $projection_source_dba->get_MetaContainer->get_production_name();
    my $target_production_name = $target_dba->get_MetaContainer->get_production_name();

    my $compara_dbname = $compara_dba->dbc->dbname;
    my $compara_host = $compara_dba->dbc->host;
    my $compara_port = $compara_dba->dbc->port;
    my $compara_user = $compara_dba->dbc->user;
    my $compara_pass = $compara_dba->dbc->pass;

    my $target_dbname = $target_dba->dbc->dbname;
    my $target_host = $target_dba->dbc->host;
    my $target_port = $target_dba->dbc->port;
    my $target_user = $target_dba->dbc->user;
    my $target_pass = $target_dba->dbc->pass;

    my $source_dbname = $projection_source_dba->dbc->dbname;
    my $source_host = $projection_source_dba->dbc->host;
    my $source_port = $projection_source_dba->dbc->port;
    my $source_user = $projection_source_dba->dbc->user;
    my $source_pass = $projection_source_dba->dbc->pass;

    my $production_dbname = $production_dba->dbc->dbname;
    my $production_host = $production_dba->dbc->host;
    my $production_port = $production_dba->dbc->port;
    my $production_user = $production_dba->dbc->user;
    my $production_pass = $production_dba->dbc->pass;

    say OUT <<DATABASES
  use warnings;
  use strict;

  use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;
  use Bio::EnsEMBL::Production::DBSQL::DBAdaptor;

  # Target species
  Bio::EnsEMBL::Registry->load_registry_from_url('mysql://$target_user:$target_pass\@$target_host:$target_port/$target_dbname?group=core&species=$target_production_name');

  # Reference species
  Bio::EnsEMBL::Registry->load_registry_from_url('mysql://$source_user:$source_pass\@$source_host:$source_port/$source_dbname?group=core&species=$source_production_name');

  Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -host    => '$source_host',
    -user    => '$source_user',
    -port    => $source_port,
    -species => '$source_production_name',
    -group   => 'core',
    -dbname  => '$source_dbname',
    -pass    => '$source_pass',
 );

  # Compara master database:
  Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->new(
    -host    => '$compara_host',
    -user    => '$compara_user',
    -port    => $compara_port,
    -species => 'compara_master',
    -dbname  => '$compara_dbname',
    -pass    => '$compara_pass',
  );

  # Production db
  Bio::EnsEMBL::Production::DBSQL::DBAdaptor->new(
    -host    => '$production_host',
    -port    => $production_port,
    -user    => '$production_user',
    -dbname  => '$production_dbname',
    -species => 'multi',
    -group   => 'production',
  );

  1;
DATABASES
  ;
  }
  close(OUT) || $self->throw('Could not close '.$self->param('registry_file'));

}

=head2 write_output

 Arg [1]    : None
 Description: Nothing to do here
 Returntype : None
 Exceptions : None

=cut

sub write_output {

}

1;
