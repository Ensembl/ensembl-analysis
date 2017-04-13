# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2017] EMBL-European Bioinformatics Institute
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

=head1 OPTIONS

-ftp_path       FTP path where the PDB chain Uniprot file is located (including file name).

-output_path    Output path where output and log files will be written.

-core_dbhost    Core database host name

-core_dbport    Core database port (default 3306)

-core_dbname    Core database name

-core_dbuser    Core database username to connect as

-core_dbpass    Core database password to use

-cs_version     Coordinate system version.

=head1 EXAMPLE USAGE

standaloneJob.pl Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadPDBProteinFeatures -ftp_path ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_uniprot.tsv.gz -output_path OUTPUT_PATH -core_dbhost genebuild3 -core_dbport 4500 -core_dbname carlos_homo_sapiens_core_89_test -core_dbuser *** -core_dbpass *** -cs_version GRCh38 

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadPDBProteinFeatures;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Tools::Utilities;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

use Net::FTP;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use File::Basename;
use File::Find;
use List::Util qw(sum);

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(run_command);
use Bio::EnsEMBL::GIFTS::DB qw(get_default_gifts_dbc fetch_latest_uniprot_enst_perfect_matches);

sub param_defaults {
    return {
      ftp_path => undef,
      output_path => undef,
      core_dbhost => undef,
      core_dbport => undef,
      core_dbname => undef,
      core_dbuser => undef,
      core_dbpass => undef,
      cs_version => undef,
    }
}

sub fetch_input {
  my $self = shift;



  return 1;
}

sub run {

  my $self = shift;

  $self->param_required('ftp_path');
  $self->param_required('output_path');
  $self->param_required('core_dbhost');
  $self->param_required('core_dbport');
  $self->param_required('core_dbname');
  $self->param_required('core_dbuser');
  $self->param_required('core_dbpass');
  $self->param_required('cs_version');
  
  # connect to the GIFTS database
  my $gifts_dbc = get_default_gifts_dbc();
  
  # connect to the core database
  my $core_dba = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
                   '-no_cache' => 1,
                   '-host'     => $self->param('core_dbhost'),
                   '-port'     => $self->param('core_dbport'),
                   '-user'     => $self->param('core_dbuser'),
                   '-pass' => $self->param('core_dbpass'),
                   '-dbname' => $self->param('core_dbname'),
  ) or die('Failed to connect to the core database.');
  
  #add / at the end of the paths if it cannot be found to avoid possible errors
  if (!($self->param('output_path') =~ /\/$/)) {
    $self->param('output_path',$self->param('output_path')."/");
  }

  # create output dir if it does not exist
  if (not -e $self->param('output_path')) {
    run_command("mkdir -p ".$self->param('output_path'),"Create output path.");
  }

  my $pdb_filepath = $self->download_pdb_file($self->param('ftp_path'),
                                              $self->param('output_path'));

  my @pdb_info = parse_pdb_file($pdb_filepath);
  my %perfect_matches = fetch_latest_uniprot_enst_perfect_matches($gifts_dbc,"Homo sapiens","GRCh38");
  insert_protein_features($core_dba,\%perfect_matches,\@pdb_info);
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

sub parse_pdb_file() {
# Parse the CSV PDB file containing the following 9 columns:
# PDB CHAIN SP_PRIMARY RES_BEG RES_END PDB_BEG PDB_END SP_BEG SP_END
# Lines starting with '#' or 'PDB' are ignored.
# Return array of hashes with SP_PRIMARY, PDB, CHAIN, RES_BEG, RES_END, SP_BEG and SP_END as keys.

  my $pdb_filepath = shift;

  open(my $pdb_fh,'<',$pdb_filepath) or die "Cannot open: $pdb_filepath";
  my @pdb_info = ();
  
  while (my $line = <$pdb_fh>) {

    next if $line =~ /^#/;
    next if $line =~ /^PDB/; 

    my ($pdb,$chain,$sp_primary,$res_beg,$res_end,undef,undef,$sp_beg,$sp_end) = split(/\s+/,$line);

    push(@pdb_info,{'PDB' => $pdb,
                    'CHAIN' => $chain,
                    'SP_PRIMARY' => $sp_primary,
                    'RES_BEG' => $res_beg,
                    'RES_END' => $res_end,
                    'SP_BEG' => $sp_beg,
                    'SP_END' => $sp_end
                   });
  }
  return @pdb_info;
}

sub insert_protein_features() {
# create and insert the protein features into the database 'dba'
# linking the ENSP proteins and the PDB entries
# from 'ref_perfect_matches' and 'pdb_info'

  my ($dba,$ref_perfect_matches,$pdb_info) = @_;

  # get list of transcript stable IDs from the keys of the perfect matches hash
  my @t_sids = keys %{$ref_perfect_matches};

  # loop through all pdb lines and find the corresponding ENSTs (if any)
  # in the perfect matches hash,
  # fetch their proteins and add their protein features
  my $ta = $dba->get_TranscriptAdaptor();
  my $pfa = $dba->get_ProteinFeatureAdaptor();
  my $analysis = new Bio::EnsEMBL::Analysis(-logic_name => 'sifts_import',
                                            -display_label => 'SIFTS import',
                                            -displayable => '1',
                                            -description => 'Protein features based on the PDB-UniProt mappings found in the EMBL-EBI PDB SIFTS data and the UniProt-ENSP mappings found in the GIFTS database.',
  );
  
  foreach my $pdb_line (@{$pdb_info}) {
    my $pdb_uniprot = $$pdb_line{'SP_PRIMARY'};    

    if (defined($$ref_perfect_matches{$pdb_uniprot})) {
      my @ensts = @{$$ref_perfect_matches{$pdb_uniprot}};
      if (scalar(@ensts) > 0) {
        foreach my $enst (@ensts) {
          my $t = $ta->fetch_by_stable_id($enst);
      
          my $translation = $t->translation();
          my $translation_sid = $translation->stable_id();

          my $pf = Bio::EnsEMBL::ProteinFeature->new(
                  -start    => $$pdb_line{'SP_BEG'},
                  -end      => $$pdb_line{'SP_END'},
                  -hseqname => $$pdb_line{'PDB'},
                  -hstart   => $$pdb_line{'RES_BEG'},
                  -hend     => $$pdb_line{'RES_END'},
                  -analysis => $analysis,
                  -hdescription => "Chain ".$$pdb_line{'CHAIN'}.". Via UniProt protein ".$$pdb_line{'SP_PRIMARY'}." isoform exact match to Ensembl protein $translation_sid",
               );
          $pfa->store($pf,$translation->dbID());
        } # foreach my enst
      } # if scalar
    } # if ensts
  } # foreach my pdb_line
}

sub write_output {
  my $self = shift;

  return 1;
}

1;
