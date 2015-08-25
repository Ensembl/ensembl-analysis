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

HiveCopyGenes.pm

=head1 DESCRIPTION

This module copies the genes whose identifiers are listed in a given file from a source database into an output database.

=head1 OPTIONS

-dbhost         database host name

-dbport         database port (default 3306)

-dbname         database name

-dbuser         database username to connect as

-dbpass         database password to use

=head1 EXAMPLE USAGE

standaloneJob.pl Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCopyGenes -dbhost genebuildX -dbuser USER_W -dbpass PASS_W -dbname DBNAME -dbport DBPORT -biotype TRANSCRIPT_BIOTYPE_TO_DELETE -delete_transcripts_path PATH_TO_SCRIPT_TO_DELETE_TRANSCRIPTS -script_name DELETE_TRANSCRIPTS_SCRIPT_NAME -output_path OUTPUT_PATH -output_file_name DELETE_TRANSCRIPTS_SCRIPT_OUTPUT_FILE_NAME

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCopyGenes;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(run_command);
use Bio::EnsEMBL::Utils::Exception qw(warning throw);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

use Net::FTP;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use File::Basename;
use File::Find;
use List::Util qw(sum);

sub param_defaults {
    return {
      copy_genes_path => '$ENSCODE/ensembl-analysis/scripts/genebuild/',
      copy_genes_script_name => 'copy_genes.pl',

      # copy_genes.pl script parameters
      logic => '', # to set the genes and transcripts analysis (logic names) (optional)
      sourcehost => '',
      sourceuser => '',
      sourceport => '3306',
      sourcepass => '',
      sourcedbname => '',
      outhost => '',
      outuser => '',
      outpass => '',
      outdbname => '',
      outport => '3306',
      dnahost => '',
      dnadbname => '',
      dnauser => '',
      dnaport => '3306',
      file => '' #$SCR9/gene_ids_to_copy.txt
    }
}

sub fetch_input {
  my $self = shift;

  return 1;
}

sub run {

  my $self = shift;
  
  if (not $self->param('sourcehost') or not $self->param('sourceuser') or not $self->param('sourcedbname') or not $self->param('outhost') or not $self->param('outuser')
      or not $self->param('outdbname') ) {
    throw("Parameters missing");
  }

#  my $num_genes_to_copy = int(run_command("wc -l ".$self->param('file'),"Counting number of genes to copy..."));
#  print "es $num_genes_to_copy\n";
#  my $num_genes_source_db = int(run_command("mysql -NB -u".$self->param('sourceuser')
#                                                                ." -h".$self->param('sourcehost')
#                                                                #." -p".$self->param('sourcepass')
#                                                                ." -D".$self->param('sourcedbname')
#                                                                ." -P".$self->param('sourceport')
#                                                                ." -e".'"'."select count(*) from gene".'"'
#                                                                ,"Counting number of genes on the source database before copying the genes..."));
#
#  my @tables_to_count = ('gene','transcript','exon','translation');
#  if ($num_genes_to_copy != $num_genes_source_db) {
#  	# if not all the genes are copied, we will only check if the gene number is right. Otherwise we would need to full-fetch all genes to count their transcripts, etc. 
#  	@tables_to_count = ('gene');
#  }
#
#  my %source_table_counts_before = ();
#  my %output_table_counts_before = ();
#  foreach my $table_to_count (@tables_to_count) {
#    my $sql_count = "select count(*) from $table_to_count";
#    $source_table_counts_before{$table_to_count} = run_command("mysql -NB -u".$self->param('sourceuser')
#                                                                ." -h".$self->param('sourcehost')
#                                                                #." -p".$self->param('sourcepass')
#                                                                ." -D".$self->param('sourcedbname')
#                                                                ." -P".$self->param('sourceport')
#                                                                ." -e".'"'.$sql_count.'"'
#                                                                ,"Counting $table_to_count rows on the source database before copying the genes...");
#    $output_table_counts_before{$table_to_count} = run_command("mysql -NB -u".$self->param('outuser')
#                                                                ." -h".$self->param('outhost')
#                                                                ." -p".$self->param('outpass')
#                                                                ." -D".$self->param('outdbname')
#                                                                ." -P".$self->param('outport')
#                                                                ." -e".'"'.$sql_count.'"'
#                                                                ,"Counting $table_to_count rows on the output database before copying the genes...");
#  }
  
  my $command = "perl ".$self->param('copy_genes_path')
                       .$self->param('copy_genes_script_name')
                       ." -sourcehost ".$self->param('sourcehost')
                       ." -sourceport ".$self->param('sourceport')
                       ." -sourcedbname ".$self->param('sourcedbname')
                       ." -outuser ".$self->param('outuser')
                       ." -outpass ".$self->param('outpass')
                       ." -outdbname ".$self->param('outdbname')
                       ." -outport ".$self->param('outport')
                       ." -outhost ".$self->param('outhost')
                       ." -dnahost ".$self->param('dnahost')
                       ." -dnadbname ".$self->param('dnadbname')
                       ." -dnauser ".$self->param('dnauser')
                       ." -dnaport ".$self->param('dnaport')
                       ." -file ".$self->param('file');

  if ($self->param('logic')) {
  	$command .= " -logic ".$self->param('logic');
  }
  
  run_command($command,"Copying genes...");
  
#  foreach my $table_to_count (@tables_to_count) {
#    my $sql_count = "select count(*) from $table_to_count";
#    run_command("mysql -NB -u".$self->param('sourceuser')
#                        ." -h".$self->param('sourcehost')
#                        #." -p".$self->param('sourcepass')
#                        ." -D".$self->param('sourcedbname')
#                        ." -P".$self->param('sourceport')
#                        ." -e".'"'.$sql_count.'"'
#                        ,"Counting $table_to_count rows on the source database after copying the genes..."
#                        ,$source_table_counts_before{$table_to_count});
#                        
#    if ($num_genes_to_copy == $num_genes_source_db) {
#      run_command("mysql -NB -u".$self->param('outuser')
#                          ." -h".$self->param('outhost')
#                          ." -p".$self->param('outpass')
#                          ." -D".$self->param('outdbname')
#                          ." -P".$self->param('outport')
#                          ." -e".'"'.$sql_count.'"'
#                          ,"Counting $table_to_count rows on the output database after copying the genes..."
#                          ,$source_table_counts_before{$table_to_count}+$output_table_counts_before{$table_to_count});
#    } else {
#      run_command("mysql -NB -u".$self->param('outuser')
#                          ." -h".$self->param('outhost')
#                          ." -p".$self->param('outpass')
#                          ." -D".$self->param('outdbname')
#                          ." -P".$self->param('outport')
#                          ." -e".'"'.$sql_count.'"'
#                          ,"Counting $table_to_count rows on the output database after copying the genes..."
#                          ,$output_table_counts_before{$table_to_count}+$num_genes_to_copy);
#    }
#  }
  return 1;
}

sub write_output {
  my $self = shift;

  return 1;
}

1;
