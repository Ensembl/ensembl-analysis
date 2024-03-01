# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

=head1 NAME 

HiveLoadPDBProteinFeatures.pm

=head1 DESCRIPTION

This module inserts protein features into an Ensembl core database based on the PDB-UniProt mappings found in the EMBL-EBI PDB SIFTS data and the UniProt-ENSP mappings found in the GIFTS database in order to make the link between PDB and ENSP having a PDB entry as a protein feature for a given ENSP protein.

It also populates the "pdb_ens" table in the GIFTS database with similar data.

=head1 OPTIONS

-ftp_path       FTP path where the PDB chain Uniprot file is located (including file name).

-output_path    Output path where output and log files will be written.

-core_dbhost    Core database host name.

-core_dbport    Core database port.

-core_dbname    Core database name.

-core_dbuser    Core database username to connect as.

-core_dbpass    Core database password to use.

-cs_version     Coordinate system version.

-rest_server    GIFTS rest server to fetch the perfect matches data from. 

=head1 EXAMPLE USAGE

standaloneJob.pl Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadPDBProteinFeatures -ftp_path https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz -output_path OUTPUT_PATH -core_dbhost genebuild3 -core_dbport 4500 -core_dbname carlos_homo_sapiens_core_89_test -core_dbuser *** -core_dbpass *** -cs_version GRCh38 -rest_server https://www.ebi.ac.uk/gifts/api/

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadPDBProteinFeatures;

use strict;
use warnings;

# Bio::DB::HTS::Faidx used in Bio::EnsEMBL::GIFTS::DB needs Perl 5.14.2
use 5.014002;
use Bio::EnsEMBL::Analysis::Tools::Utilities;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');
use Bio::EnsEMBL::Analysis::Runnable::MakePDBProteinFeatures;

use Net::FTP;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use File::Basename;
use File::Find;
use List::Util qw(sum);

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(run_command);

sub param_defaults {
    my ($self) = @_;
 
    return {
      %{$self->SUPER::param_defaults},
      ftp_path => undef,
      output_path => undef,
      core_dbhost => undef,
      core_dbport => undef,
      core_dbname => undef,
      core_dbuser => undef,
      core_dbpass => undef,
      cs_version => undef,
      species => undef,
      rest_server => undef,
    }
}

sub fetch_input {
  my $self = shift;

  $self->param_required('ftp_path');
  $self->param_required('output_path');
  $self->param_required('core_dbhost');
  $self->param_required('core_dbport');
  $self->param_required('core_dbname');
  $self->param_required('core_dbuser');
  $self->param_required('core_dbpass');
  $self->param_required('cs_version');
  $self->param_required('species');
  $self->param_required('rest_server');

  #add / at the end of the paths if it cannot be found to avoid possible errors
  if (!($self->param('output_path') =~ /\/$/)) {
    $self->param('output_path',$self->param('output_path')."/");
  }

  # create output dir if it does not exist
  if (not -e $self->param('output_path')) {
    run_command("mkdir -p ".$self->param('output_path'),"Create output path.");
  }
  
  # connect to the core database
  my $core_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                   '-no_cache' => 1,
                   '-host'     => $self->param('core_dbhost'),
                   '-port'     => $self->param('core_dbport'),
                   '-user'     => $self->param('core_dbuser'),
                   '-pass' => $self->param('core_dbpass'),
                   '-dbname' => $self->param('core_dbname'),
  ) or die('Failed to connect to the core database.');

  $self->hrdb_set_con($core_dba,"core");

  # download the PDB file from the FTP
  my $pdb_filepath = $self->download_pdb_file($self->param('ftp_path'),
                                              $self->param('output_path'));

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::MakePDBProteinFeatures->new(
    -analysis => new Bio::EnsEMBL::Analysis(-logic_name => 'sifts_import',
                                            -db => 'sifts',
                                            #-db_version => , it will be populated in the MakePDBProteinFeatures module when parsing the file
                                            -db_file => $pdb_filepath,
                                            -display_label => 'SIFTS import',
                                            -displayable => '1',
                                            -description => 'Protein features based on the PDB-UniProt mappings found in the EMBL-EBI PDB SIFTS data and the UniProt-ENSP mappings found in the GIFTS database.'),
    -core_dba => $self->hrdb_get_con("core"),
    -pdb_filepath => $pdb_filepath,
    -species => $self->param('species'),
    -cs_version => $self->param('cs_version'),
    -rest_server => $self->param('rest_server')
    );
  $self->runnable($runnable);

  return 1;
}

sub run {
  my $self = shift;

  foreach my $runnable (@{$self->runnable}) {
    $runnable->run();
    $self->output($runnable->output);
  }
}

sub write_output {
  my $self = shift;

  # insert the Ensembl-PDB links into the protein_feature table in the core database
  # and add its associated xrefs
  $self->insert_protein_features();

  return 1;
}

sub download_pdb_file() {
# download the SIFTS PDB chain file from the ftp_path (including file name)
# to the local directory 'local_dir' and return the file path of the downloaded file

  my ($self,$ftp_path,$local_dir) = @_;

  my $wget_verbose = "-nv";

  if (system("wget ".$wget_verbose." -nH -P ".$local_dir." ".$ftp_path)) {
    $self->throw("Could not download SIFTS PDB chain file from ".$ftp_path." to ".$local_dir.". Please, check that both paths are valid.");
  }
  else {
    print($ftp_path." file was downloaded successfully.\n");
  }

  if (system("gunzip $local_dir/*.gz")) {
    $self->throw("Could not gunzip the .gz files.");
  } else {
  	print("File gunzipped successfully.\n");
  }
  
  my $filename_without_gz = $ftp_path;
  $filename_without_gz =~ s/.*\///;
  $filename_without_gz =~ s/\.gz$//;
  return $local_dir.$filename_without_gz;
}

sub insert_protein_features() {
# insert the protein features and their associated xrefs into the database 'core_dba'

  my $self = shift;

  my $core_dba = $self->hrdb_get_con("core");
  my $pfa = $core_dba->get_ProteinFeatureAdaptor();

  foreach my $pf_hashref (@{$self->output}) {  
    my ($translation_id) = keys %$pf_hashref;
    $pfa->store($pf_hashref->{$translation_id},$translation_id);
    
    # This is not done because we only want to insert "real" IDs in the dbprimary_acc column.
    # Fake "abcd.A" PDB.chain IDs are not linkable to any external DB at the moment.
    #$self->insert_protein_features_xrefs($core_dba,$pf_hashref->{$translation_id},$translation_id);
  }
}

sub insert_protein_features_xrefs {
# Inserts a PDB protein feature translation xref
  my ($self,$db_adaptor,$pf,$translation_id) = @_;

  my $dbe_adaptor = $db_adaptor->get_DBEntryAdaptor();
  my ($pdb_acc,$pdb_chain) = split(/\./,$pf->hseqname());
  my $pf_description = $pf->hdescription();
  my @pf_description_array = split(' ',$pf_description);
  
  my $pdb_xref = new Bio::EnsEMBL::DBEntry(
                                             -adaptor => $dbe_adaptor,
                                             -primary_id => $pdb_acc.".".$pdb_chain,
                                             -version => 0,
                                             -dbname  => 'PDB',
                                             -release => 1,
                                             -display_id => $pdb_acc.".".$pdb_chain,
                                             -description => $pf_description,
                                             -priority => 5,
                                             -db_display_name => 'PDB',
                                             -info_type => 'DEPENDENT',
                                             -type => 'MISC'
                                           );
  $pdb_xref->status('XREF');
  $pdb_xref->analysis($pf->analysis());
  $dbe_adaptor->store($pdb_xref,$translation_id,'Translation');
}

1;
