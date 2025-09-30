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

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(run_command send_email);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

use Net::FTP;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
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
      email => undef, # email to send the report to
      from => 'HiveDeleteTranscripts@ebi.ac.uk',
      fix_broken_genes => 0
    }
}

sub fetch_input {
  my $self = shift;

  if($self->param('skip_analysis')) {
    $self->complete_early('Skip analysis flag is enabled, so no cleaning will occur');
  }

  return 1;
}

sub run {

  my $self = shift;

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

  my $gene_ids_file;
  if($self->param('transcript_ids_file')) {
    my $transcript_ids_file = $self->param('transcript_ids_file');
    $gene_ids_file = 'gene_ids_to_check.txt';
    # delete the transcripts
    run_command("perl ".$self->param('delete_transcripts_path')
                       .$self->param('delete_transcripts_script_name')
                       ." -dbhost ".$self->param('dbhost')
                       ." -dbuser ".$self->param('dbuser')
                       ." -dbpass ".$self->param('dbpass')
                       ." -dbname ".$self->param('dbname')
                       ." -dbport ".$self->param('dbport')
                       ." ".$transcript_ids_file
                       ." > ".$self->param('output_path').$self->param('output_file_name'));

  } else {
    $self->param_required('biotype');
    my $transcript_ids_file = $self->param('biotype')."_transcript_ids.txt";
    $gene_ids_file = $self->param('biotype')."_gene_ids.txt";
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

    } # End else

    # parse the deletion output file to find the empty genes
    my $cmd = "grep empty ".$self->param('output_path').$self->param('output_file_name');
    my @empty_gene_id_list = `$cmd`;
    open(OUT,">".$self->param('output_path')."/".$gene_ids_file);
    foreach my $line (@empty_gene_id_list) {
      $line =~ /Gene.+\(id = (\d+)\)/;
      say OUT $1;
    }
    close OUT;

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

    # prepare report

    # 'cat' at the end prevents the command from failing as if there was an error due to grep not finding any pattern
    my $broken_genes = run_command("grep -v Deleted ".$self->param('output_path').$self->param('output_file_name')." | grep -v Gene | grep -v crash | cat","Preparing list of broken genes for email...");

  my $empty_genes = run_command("grep empty ".$self->param('output_path').$self->param('output_file_name')." | cat","Preparing list of empty genes for email...");

  my $body;

  if ($self->param('fix_broken_genes')) {
    $body = <<'END_BODY';
Hi,

Please find below both the genes which were broken (and have been fixed) and the genes which have been deleted because they were empty after deleting the corresponding transcripts.

Regards,
The Ensembl Genebuild team's automatic pipeline powered by eHive

----- broken and empty genes below -----

END_BODY
  } else {
    $body = <<'END_BODY';
Hi,
 
Please find below both the genes which are broken and the genes which have been deleted because they were empty after deleting the corresponding transcripts.

Regards,
The Ensembl Genebuild team's automatic pipeline powered by eHive

----- broken and empty genes below -----

END_BODY
  }

  $body .= $broken_genes;
  $body .= $empty_genes;

  # fix broken genes
  if ($self->param('fix_broken_genes')) {

    my $current_gene_id;
    my $current_transcript_id;
    my @lines = split /\n/,$broken_genes;
    foreach my $line (@lines) {

      if ($line =~ /^WARNING:(\s+)(\S+)(\s+)(\S+)(\s+)(\S+)(\s+)(\S+)(\s+)(\S+)(\s+)(\S+)\),/) {
        $current_gene_id = $12;

        run_command("mysql -NB -u".$self->param('dbuser')
                            ." -h".$self->param('dbhost')
                            ." -p".$self->param('dbpass')
                            ." -D".$self->param('dbname')
                            ." -P".$self->param('dbport')
                            ." -e".'"'."insert into gene (biotype,analysis_id,seq_region_id,seq_region_start,seq_region_end,seq_region_strand,display_xref_id,source,description,is_current,canonical_transcript_id,stable_id,version,created_date,modified_date) (select biotype,analysis_id,seq_region_id,seq_region_start,seq_region_end,seq_region_strand,display_xref_id,source,description,is_current,canonical_transcript_id,stable_id,version,created_date,modified_date from gene where gene_id=$current_gene_id);".'"',
                            "Inserting copy of gene row with gene id $current_gene_id");
      } elsif ($line =~ /(\s+)(\S+)(\s+)(\S+)(\s+)(\S+)(\s+)(\S+)\)/) {
        if ($4 eq '(id') {
          $current_transcript_id = $8;

          run_command("mysql -NB -u".$self->param('dbuser')
                              ." -h".$self->param('dbhost')
                              ." -p".$self->param('dbpass')
                              ." -D".$self->param('dbname')
                              ." -P".$self->param('dbport')
                              ." -e".'"'."update transcript set gene_id=(select max(gene_id) from gene) where transcript_id in ($current_transcript_id);".'"',
                              "Updating gene id for transcript $current_transcript_id");
        }
      }
    }
    
    # set new gene coordinates
    run_command("mysql -NB -u".$self->param('dbuser')
                        ." -h".$self->param('dbhost')
                        ." -p".$self->param('dbpass')
                        ." -D".$self->param('dbname')
                        ." -P".$self->param('dbport')
                        ." -e".'"'."update gene g,(select t.gene_id,min(t.seq_region_start) as minim,max(t.seq_region_end) as maxim from transcript t,gene g where g.gene_id=t.gene_id group by t.gene_id) as c set g.seq_region_start=c.minim,g.seq_region_end=c.maxim where c.gene_id=g.gene_id;".'"',
                              "Updating gene coordinates...");
  }

  # send report
  if ($self->param('email')) {
    send_email($self->param('email'),$self->param('from'),$self->param('biotype')." transcripts deleted report",$body);
  }
  
  return 1;
}

sub write_output {
  my $self = shift;

  return 1;
}

1;
