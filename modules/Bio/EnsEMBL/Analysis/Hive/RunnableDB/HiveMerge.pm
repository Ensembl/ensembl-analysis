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

=head1 NAME 

HiveMerge.pm

=head1 DESCRIPTION

This module runs merge.pl on eHive. Please see merge.pl for further details about the parameters.

=head1 OPTIONS

-dbhost         database host name

-dbport         database port

-dbname         database name

-dbuser         database username to connect as

-dbpass         database password to use

=head1 EXAMPLE USAGE

standaloneJob.pl Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMerge
                               -ensembl_analysis_base $ENSCODE/ensembl-analysis
                               -host_secondary ENSEMBLHOST
                               -port_secondary ENSEMBLPORT
                               -user_secondary READ_ONLY_USER
                               -password_secondary READ_ONLY_PASS
                               -database_secondary ENSEMBLDBNAME
                               -host_primary VEGADBHOST
                               -port_primary VEGADBPORT
                               -user_primary READ_ONLY_USER
                               -password_primary READ_ONLY_PASS
                               -database_primary VEGADBNAME
                               -host_output MERGEDBHOST
                               -port_output MERGEDBPORT
                               -user_output WRITE_USER
                               -password_output WRITE_PASS
                               -database_output MERGEDBNAME
                               
                               # Tagging:  Will be used as suffix for logic names ("_tag") and for
                               # source.  With the default settings, merged genes and transcripts will
                               # get the source "secondary_primary".
                               
                               -secondary_tag ensembl
                               -primary_tag havana
                               
                               # Xrefs:  The format is a comma-separated list of
                               # "db_name,db_display_name,type"
                               
                               -primary_gene_xref "OTTG,Havana gene,ALT_GENE"
                               -primary_transcript_xref "OTTT,Havana transcript,ALT_TRANS"
                               -primary_translation_xref "OTTP,Havana translation,MISC"
                               
                               # as the chunks (and a job per chunk) are created in the step before,
                               # these parameters would define how many jobs per chunk we want, just 1 as we don't want chunks of chunks
                               # and we cannot use the LSF job index on the ehive to create chunks of chunks here anyway
                               -njobs => 1,
                               -job => 1,
                               
                               -file => OUTPUT/vega_genes_for_merge.ids_chunk_2.txt, this parameter would normally come from 'chunk_genes_for_merge' output, see FileFactory.pm
=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveMerge;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(run_command);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

use Net::FTP;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use File::Basename;
use File::Find;
use List::Util qw(sum);

sub param_defaults {
    return {
      ensembl_analysis_base => '',
      host_secondary => '',
      port_secondary => '',
      user_secondary => '',
      password_secondary => '',
      database_secondary => '',
      host_primary => '',
      port_primary => '',
      user_primary => '',
      password_primary => '',
      database_primary => '',
      host_dna => '',
      port_dna => '',
      user_dna => '',
      database_dna => '',
      host_output => '',
      port_output => '',
      user_output => '',
      password_output => '',
      database_output => '',
      secondary_include => '',
      secondary_exclude => '',
      primary_include => '',
      primary_exclude => '',
      secondary_tag => '',
      primary_tag => '',
      primary_xref => '',
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
              " --port_secondary=".$self->param('port_secondary').
              " --user_secondary=".$self->param('user_secondary').
              " --password_secondary=".$self->param('password_secondary').
              " --database_secondary=".$self->param('database_secondary').
              " --host_primary=".$self->param('host_primary').
              " --port_primary=".$self->param('port_primary').
              " --user_primary=".$self->param('user_primary').
              " --password_primary=".$self->param('password_primary').
              " --database_primary=".$self->param('database_primary').
              " --host_dna=".$self->param('host_dna').
              " --port_dna=".$self->param('port_dna').
              " --user_dna=".$self->param('user_dna').
              " --database_dna=".$self->param('database_dna').
              " --host_output=".$self->param('host_output').
              " --port_output=".$self->param('port_output').
              " --user_output=".$self->param('user_output').
              " --password_output=".$self->param('password_output').
              " --database_output=".$self->param('database_output').
              " --secondary_include=".$self->param('secondary_include').
              " --secondary_exclude=".$self->param('secondary_exclude').
              " --primary_include=".$self->param('primary_include').
              " --primary_exclude=".$self->param('primary_exclude').
              " --secondary_tag=".$self->param('secondary_tag').
              " --primary_tag=".$self->param('primary_tag').
              " --primary_xref='".$self->param('primary_xref')."'".
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
