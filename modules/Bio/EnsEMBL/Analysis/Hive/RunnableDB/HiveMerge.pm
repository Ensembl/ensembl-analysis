# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

HiveMerge.pm

=head1 DESCRIPTION

This module runs merge.pl on eHive.

=head1 OPTIONS

-dbhost         database host name

-dbport         database port (default 3306)

-dbname         database name

-dbuser         database username to connect as

-dbpass         database password to use

=head1 EXAMPLE USAGE

standaloneJob.pl Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMerge -dbhost genebuildX -dbuser USER_W -dbpass PASS_W -dbname DBNAME -dbport DBPORT

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMerge;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Tools::Utilities;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

use Net::FTP;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use File::Basename;
use File::Find;
use List::Util qw(sum);

sub param_defaults {
    return {
      ensembl_analysis_base => '',
      host_secondary => '',
      user_secondary => '',
      password_secondary => '',
      database_secondary => '',
      host_primary => '',
      user_primary => '',
      password_primary => '',
      database_primary => '',
      host_dna => '',
      port_dna => '',
      user_dna => '',
      database_dna => '',
      host_ccds => '',
      user_ccds => '',
      password_ccds => '',
      database_ccds => '',
      host_output => '',
      user_output => '',
      password_output => '',
      database_output => '',
      secondary_include => '',
      secondary_exclude => '',
      primary_include => '',
      primary_exclude => '',
      secondary_tag => '',
      primary_tag => '',
      primary_gene_xref => '',
      primary_transcript_xref => '',
      primary_translation_xref => '',
      njobs => '',
      job => '',
      file => '',
    }
}

sub fetch_input {
  my $self = shift;

  return 1;
}

sub run {

  my $self = shift;
  
  if (not $self->param('host_secondary') or not $self->param('host_primary')) {
    throw("Parameters missing");
  }

  # add / at the end of the paths if it cannot be found to avoid possible errors
  #if (!($self->param('output_path') =~ /\/$/)) {
  #  $self->param('output_path',$self->param('output_path')."/");
  #}

  # create output dir if it does not exist
  #if (not -e $self->param('output_path')) {
  #  run_command("mkdir -p ".$self->param('output_path'),"Create output path.");
  #}

  run_command("perl ".$self->param('ensembl_analysis_base')."/scripts/Merge/merge.pl".
              " --host_secondary=".$self->param('host_secondary').
              " --user_secondary=".$self->param('user_secondary').
              " --password_secondary=".$self->param('password_secondary').
              " --database_secondary=".$self->param('database_secondary').
              " --host_primary=".$self->param('host_primary').
              " --user_primary=".$self->param('user_primary').
              " --password_primary=".$self->param('password_primary').
              " --database_primary=".$self->param('database_primary').
              " --host_dna=".$self->param('host_dna').
              " --port_dna=".$self->param('port_dna').
              " --user_dna=".$self->param('user_dna').
              " --database_dna=".$self->param('database_dna').
              " --host_ccds=".$self->param('host_ccds').
              " --user_ccds=".$self->param('user_ccds').
              " --password_ccds=".$self->param('password_ccds').
              " --database_ccds=".$self->param('database_ccds').
              " --host_output=".$self->param('host_output').
              " --user_output=".$self->param('user_output').
              " --password_output=".$self->param('password_output').
              " --database_output=".$self->param('database_output').
              " --secondary_include=".$self->param('secondary_include').
              " --secondary_exclude=".$self->param('secondary_exclude').
              " --primary_include=".$self->param('primary_include').
              " --primary_exclude=".$self->param('primary_exclude').
              " --secondary_tag=".$self->param('secondary_tag').
              " --primary_tag=".$self->param('primary_tag').
              " --primary_gene_xref='".$self->param('primary_gene_xref')."'".
              " --primary_transcript_xref='".$self->param('primary_transcript_xref')."'".
              " --primary_translation_xref='".$self->param('primary_translation_xref')."'".
              " --njobs=".$self->param('njobs').
              " --job=".$self->param('job').
              " --file=".$self->param('file').
              " &> ".$self->param('file').".merge-run.out"
              );

  return 1;
}

sub write_output {
  my $self = shift;

  return 1;
}

1;
