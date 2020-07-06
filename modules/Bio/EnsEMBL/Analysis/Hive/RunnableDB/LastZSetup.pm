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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadmRNAs

=head1 SYNOPSIS


=head1 DESCRIPTION

Module to load mRNA into a customised table in the Hive database

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::LastZSetup;

use strict;
use warnings;
use feature 'say';

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my ($self) = @_;

  say "Creating database adaptors";
  my $compara_dba = $self->hrdb_get_dba($self->param_required('compara_db'));
  $self->hrdb_set_con($compara_dba,'compara_db');

  my $target_dba = $self->hrdb_get_dba($self->param_required('target_db'));
  $self->hrdb_set_con($target_dba,'target_db');

  my $projection_source_dba = $self->hrdb_get_dba($self->param_required('projection_source_db'));
  $self->hrdb_set_con($projection_source_dba,'projection_source_db');

  my $pipeline_dba = $self->hrdb_get_dba($self->param_required('pipeline_db'));
  $self->hrdb_set_con($pipeline_dba,'pipeline_db');

  say "...finished creating database adaptors";

}


sub run {
  my ($self) = @_;

  my $compara_dba =  $self->hrdb_get_con('compara_db');
  my $target_dba = $self->hrdb_get_con('target_db');
  my $projection_source_dba = $self->hrdb_get_con('projection_source_db');
  my $pipeline_dba = $self->hrdb_get_con('pipeline_db');
  my $source_assembly_name = $projection_source_dba->get_MetaContainer->single_value_by_key('assembly.default');
  my $target_assembly_name = $target_dba->get_MetaContainer->single_value_by_key('assembly.default');
  my $source_production_name = $projection_source_dba->get_MetaContainer->get_production_name();
  my $target_production_name = $target_dba->get_MetaContainer->get_production_name();
  my $databases_conf_path = $self->param_required('registry_path');

  unless($databases_conf_path) {
    $self->throw($databases_conf_path." does not exist");
  }

  my $update_genome_db_cmd = 'perl '.$self->param_required('compara_genome_db_update_path').
                             ' --reg_conf '.$databases_conf_path.
                             ' --compara '.$self->param_required('compara_db_url').
                             ' --species "'.$target_production_name.'" --force';

  say "Running update script";
  if(system($update_genome_db_cmd)) {
     $self->throw("The update script exited with a non-zero exit code. Commandline used:\n".$update_genome_db_cmd);
  }
  say "... finished running update script";

  say "Fetching genome db ids from the compara db";
  my $sql = "SELECT genome_db_id FROM genome_db WHERE assembly=? and name=?";
  my $sth = $compara_dba->dbc->prepare($sql);
  $sth->bind_param(1,$source_assembly_name);
  $sth->bind_param(2,$source_production_name);
  $sth->execute();
  my $source_genome_db_id = $sth->fetchrow_array();

  $sth = $compara_dba->dbc->prepare($sql);
  $sth->bind_param(1,$target_assembly_name);
  $sth->bind_param(2,$target_production_name);
  $sth->execute();
  my $target_genome_db_id = $sth->fetchrow_array();

  unless($source_genome_db_id && $target_genome_db_id && $source_genome_db_id != $target_genome_db_id) {
    $self->throw("Something went wrong with fetching the genome db ids from the compara db");
  }

  say "...fetched the following genome db ids: ".$target_genome_db_id.",".$source_genome_db_id;

  my $insert_mlss_cmd = 'perl '.$self->param_required('compara_mlss_script_path').
                             ' --method_link_type  LASTZ_NET '.
                             ' --genome_db_id '.$target_genome_db_id.','.$source_genome_db_id.
                             ' --source "ensembl" '.
                             ' --reg_conf '.$self->param_required('compara_mlss_reg_conf_path').
                             ' --compara '.$self->param('compara_db_url').
                             ' --force 1';
# Note this won't work if it thinks the mlss already exists. So for this there is some sql later
# on to insert the tag directly
#                             ' --ref_species '.$source_production_name.
  say "Running mlss script";
  if(system($insert_mlss_cmd)) {
    $self->throw("The mlss script exited with a non-zero exit code. Commandline used:\n".$insert_mlss_cmd);
  }
  say "...finished running mlss script";

  say "Fetching mlss id from the db";
  my $ss_id_sql = "SELECT max(ss1.species_set_id) from species_set ss1 join species_set ss2 using(species_set_id) where ss1.species_set_id=ss2.species_set_id and ss1.genome_db_id=? and ss2.genome_db_id=?";
  $sth = $compara_dba->dbc->prepare($ss_id_sql);
  $sth->bind_param(1,$source_genome_db_id);
  $sth->bind_param(2,$target_genome_db_id);
  $sth->execute();
  my $ss_id = $sth->fetchrow_array();

  my $mlss_id_sql = "SELECT method_link_species_set_id FROM method_link_species_set WHERE species_set_id=?";
  $sth = $compara_dba->dbc->prepare($mlss_id_sql);
  $sth->bind_param(1,$ss_id);
  $sth->execute();
  my $mlss_id = $sth->fetchrow_array();

  unless($mlss_id) {
    $self->throw("Failed to fetch an mlss id. SQL used:\n".$mlss_id_sql);
  }
  say "...fetched the following mlss id: ".$mlss_id;

  $self->param('mlss_id',$mlss_id);

#  say "Inserting tag info";
#  my $set_reference_mlss_tag_sql = "INSERT INTO method_link_species_set_tag VALUES (?, 'reference_species', ?)";
#  my $pipe_sth = $pipeline_dba->dbc->prepare($set_reference_mlss_tag_sql);
#  $pipe_sth = $pipeline_dba->dbc->prepare($set_reference_mlss_tag_sql);
#  $pipe_sth->bind_param(1,$mlss_id);
#  $pipe_sth->bind_param(2,$source_production_name);
#  $pipe_sth->execute();
#  say "...finished inserting tag info";

  say "Updating release numbers for any previous assemblies of the target species";
  my $update_old_assemblies_sql = "UPDATE genome_db SET first_release=77 WHERE genome_db_id != ? AND name = ?";
  $sth = $compara_dba->dbc->prepare($update_old_assemblies_sql);
  $sth->bind_param(1,$target_genome_db_id);
  $sth->bind_param(2,$target_production_name);
  $sth->execute();

  $update_old_assemblies_sql = "UPDATE genome_db SET last_release=77 WHERE genome_db_id != ? AND name = ?";
  $sth = $compara_dba->dbc->prepare($update_old_assemblies_sql);
  $sth->bind_param(1,$target_genome_db_id);
  $sth->bind_param(2,$target_production_name);
  $sth->execute();
  say "...finished updating release numbers";
}


sub write_output {
  my ($self) = @_;
}


sub post_cleanup {
}

1;
