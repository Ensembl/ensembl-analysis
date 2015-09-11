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

HiveDeleteTranscripts.pm

=head1 DESCRIPTION

This module deletes transcripts by biotype from a given database. It also deletes the empty genes if there is any after deleting the transcripts and a report is sent by email if required.

=head1 OPTIONS

-dbhost         database host name

-dbport         database port (default 3306)

-dbname         database name

-dbuser         database username to connect as

-dbpass         database password to use

=head1 EXAMPLE USAGE

standaloneJob.pl Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDeleteTranscripts -dbhost genebuildX -dbuser USER_W -dbpass PASS_W -dbname DBNAME -dbport DBPORT -biotype TRANSCRIPT_BIOTYPE_TO_DELETE -delete_transcripts_path PATH_TO_SCRIPT_TO_DELETE_TRANSCRIPTS -script_name DELETE_TRANSCRIPTS_SCRIPT_NAME -output_path OUTPUT_PATH -output_file_name DELETE_TRANSCRIPTS_SCRIPT_OUTPUT_FILE_NAME

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveDeleteTranscripts;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use Bio::EnsEMBL::Analysis::Tools::Utilities;
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
      dbhost => undef,
      dbname => undef,
      dbuser => undef,
      dbpass => undef,
      dbport => 3306,
      biotype => undef,
      delete_transcripts_path => '$ENSCODE/ensembl-analysis/scripts/genebuild/',
      delete_genes_path => '$ENSCODE/ensembl-analysis/scripts/genebuild/',
      delete_transcripts_script_name => 'delete_transcripts.pl',
      delete_genes_script_name => 'delete_genes.pl',
      output_path => undef,
      output_file_name => 'delete_transcripts_'.time().'.out',
      email => 'report_to@ebi.ac.uk',
      from => 'your_email@ebi.ac.uk'
    }
}

sub fetch_input {
  my $self = shift;

  return 1;
}

sub run {

  my $self = shift;
  
  $self->param_required('biotype');
  $self->param_required('dbhost');
  $self->param_required('dbname');
  $self->param_required('dbuser');
  $self->param_required('dbpass');
  $self->param_required('output_path');

  #add / at the end of the paths if it cannot be found to avoid possible errors
  if (!($self->param('output_path') =~ /\/$/)) {
    $self->param('output_path',$self->param('output_path')."/");
  }

  # create output dir if it does not exist
  if (not -e $self->param('output_path')) {
    run_command("mkdir -p ".$self->param('output_path'),"Create output path.");
  }

  my $transcript_ids_file = $self->param('biotype')."_transcript_ids.txt";
  my $gene_ids_file = $self->param('biotype')."_gene_ids.txt";
  my $sql_get_biotype = 'select transcript_id from transcript where biotype='."'".$self->param('biotype')."'";
  my $sql_count_biotype = 'select count(transcript_id) from transcript where biotype='."'".$self->param('biotype')."'";

  # get the transcript identifiers for the given transcript biotype and write them into a file
  run_command("mysql -NB -u".$self->param('dbuser')
                      ." -h".$self->param('dbhost')
                      ." -p".$self->param('dbpass')
                      ." -D".$self->param('dbname')
                      ." -P".$self->param('dbport')
                      ." -e".'"'.$sql_get_biotype.'"'
                      ." > ".$self->param('output_path').$transcript_ids_file);

  # delete the transcripts
  run_command("perl ".$self->param('delete_transcripts_path')
                     .$self->param('delete_transcripts_script_name')
                     ." -dbhost ".$self->param('dbhost')
                     ." -dbuser ".$self->param('dbuser')
                     ." -dbpass ".$self->param('dbpass')
                     ." -dbname ".$self->param('dbname')
                     ." -dbport ".$self->param('dbport')
                     ." ".$self->param('output_path').$transcript_ids_file
                     ." > ".$self->param('output_path').$self->param('output_file_name'));
  
  # check that there is no transcript of the given biotype after the deletion
  run_command("mysql -NB -u".$self->param('dbuser')
                      ." -h".$self->param('dbhost')
                      ." -p".$self->param('dbpass')
                      ." -D".$self->param('dbname')
                      ." -P".$self->param('dbport')
                      ." -e".'"'.$sql_count_biotype.'"'
             ,"Checking deleted transcripts..."
             ,0); # expected transcript count

  # parse the deletion output file to find the empty genes
  run_command("grep empty ".$self->param('output_path').$self->param('output_file_name').' | awk '."'".'{print $5}'."'"." | cut -d')' -f1 > ".$self->param('output_path').$gene_ids_file,"Generating file containing empty gene IDs...");
  
  # delete the genes which are empty after the transcripts deletion
  run_command("perl ".$self->param('delete_genes_path')
                     .$self->param('delete_genes_script_name')
                     ." -dbhost ".$self->param('dbhost')
                     ." -dbuser ".$self->param('dbuser')
                     ." -dbpass ".$self->param('dbpass')
                     ." -dbname ".$self->param('dbname')
                     ." -dbport ".$self->param('dbport')
                     ." -idfile ".$self->param('output_path').$gene_ids_file
                     ." 2>> ".$self->param('output_path').$self->param('output_file_name'));
  
  # check that the right number of genes has been deleted
  my $num_genes_to_delete = run_command("grep empty ".$self->param('output_path').$self->param('output_file_name')." | wc -l","Counting number of empty genes to delete...");
  run_command("grep Deleted ".$self->param('output_path').$self->param('output_file_name')." | grep -v transcript | wc -l","Checking that number of genes to delete and number of deleted genes match...",$num_genes_to_delete);
  
  # send report
  if ($self->param('email')) {
  	
  	# 'cat' at the end prevents the command from failing as if there was an error due to grep not finding any pattern
  	my $broken_genes = run_command("grep -v Deleted ".$self->param('output_path').$self->param('output_file_name')." | grep -v Gene | grep -v crash | cat","Preparing list of broken genes for email...");
  	
  	my $empty_genes = run_command("grep empty ".$self->param('output_path').$self->param('output_file_name')." | cat","Preparing list of empty genes for email...");
  	
  	my $body = <<'END_BODY';
Hi,
 
Please find below both the genes which are broken and the genes which have been deleted because they were empty after deleting the corresponding transcripts.

Regards,
The Ensembl Genebuild team's automatic pipeline powered by eHive

----- broken and empty genes below -----

END_BODY

    $body .= $broken_genes;
    $body .= $empty_genes;

    send_email($self->param('email'),$self->param('from'),$self->param('biotype')." transcripts deleted report",$body);
  }
  
  return 1;
}

sub write_output {
  my $self = shift;

  return 1;
}

1;
