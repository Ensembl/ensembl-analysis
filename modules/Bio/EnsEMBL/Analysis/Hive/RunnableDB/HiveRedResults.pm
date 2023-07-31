=head1 LICENSE

Copyright [2018-2020] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRedResults

=head1 SYNOPSIS


=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveRedResults;

use strict;
use warnings;
use feature 'say';
use File::Spec::Functions;
use DateTime;
use Bio::EnsEMBL::Hive::Utils qw(destringify);
use Bio::EnsEMBL::IO::Parser::Fasta;
use Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my ($self) = @_;

  my $path_to_genomic = $self->param_required('path_to_genomic_fasta');
  my $assembly_registry_dba = new Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor(%{$self->param_required('assembly_registry_db')});
  $self->hrdb_set_con($assembly_registry_dba,'assembly_registry_db');
}


sub run {
  my ($self) = @_;
  my $full_path = catfile($self->param_required('path_to_genomic_fasta'),'output');
  my $gca = $self->param_required('iid');
  my ($chain,$version) = $self->split_gca($gca);
  my $assembly_registry_dba = $self->hrdb_get_con('assembly_registry_db');
  my $species_name = $assembly_registry_dba->fetch_species_name_by_gca($gca);
  $species_name =~ s/\s+/\_/g;
  $species_name = lc($species_name);

  unless($species_name) {
    $self->throw("Could not find species name for the GCA in the assembly registry. GCA used: ".$gca);
  }


  my $output_file = catfile($self->param_required('path_to_repeat_libs'),$gca."_red.msk");
  unless (-e $output_file) {
     $self->throw("No repeat library file exist for assembly $gca. Do check the previous analysis output. Path:\n".$output_file);
  }

  $self->update_registry($gca,$assembly_registry_dba);
}


sub write_output {
  my ($self) = @_;
}


sub update_registry{
  #Update registry to show that library has been generated
  my ($self,$gca,$assembly_registry_dba) = @_;
  my ($chain,$version) = $self->split_gca($gca);
  my $dt = DateTime->now();
  my $date = $dt->ymd;
  $assembly_registry_dba = $self->hrdb_get_con('assembly_registry_db');
  my $library_name = $gca . '_red_msk';
  my $sql = "update repeat_library_status set library_status =? where assembly_accession = ? and library_source = ?";
  my $sth = $assembly_registry_dba->dbc->prepare($sql);
  $sth->bind_param(1,'completed');
  $sth->bind_param(2,$gca);
  $sth->bind_param(3,'red');
  unless($sth->execute()){
    throw("Could not update repeatmodeler status for assembly with accession ".$gca);
  }
  $sql = "update repeat_library_status set date_completed =? where assembly_accession = ? and library_source = ?";
  $sth = $assembly_registry_dba->dbc->prepare($sql);
  $sth->bind_param(1,$date);
  $sth->bind_param(2,$gca);
  $sth->bind_param(3,'red');
  unless($sth->execute()){
    throw("Could not update completion date for assembly with accession ".$gca);
  }
  $sql = "update repeat_library_status set library_name =? where assembly_accession = ? and library_source = ?";
  $sth = $assembly_registry_dba->dbc->prepare($sql);
  $sth->bind_param(1,$library_name);
  $sth->bind_param(2,$gca);
  $sth->bind_param(3,'red');
  unless($sth->execute()){
    throw("Could not update completion date for assembly with accession ".$gca);
  }
}

sub split_gca {
  my ($self,$chain_version) = @_;
  unless($chain_version =~ /^(GCA\_\d{9})\.(\d+)$/) {
    $self->throw("Could not parse versioned GCA. GCA used: ".$chain_version);
  }

  my $chain = $1;
  my $version = $2;

  return($chain,$version);
}

1;
