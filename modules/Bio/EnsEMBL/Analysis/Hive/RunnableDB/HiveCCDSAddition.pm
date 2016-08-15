# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016] EMBL-European Bioinformatics Institute
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

HiveCCDSAddition.pm

=head1 DESCRIPTION

This module:

- emails the missing CCDS
- copies the missing CCDS from the CCDS database to the specified target database
- adds the CCDS transcripts as supporting features
- sets the CCDS analyses

=head1 OPTIONS

-copy_genes_path            Directory where the script to copy the genes is located.
-copy_genes_script_name     Filename of the script to copy the genes.
-add_ccds_path              Directory where the script to add the CCDS as supporting features is located.
-add_ccds_script_name       Filename of the script to add the CCDS as supporting features.
-email                      Email address to send the missing CCDS report to.
-from                       Email address where the recipient in 'email' should reply to in case of questions.
-ccds_comparison_output_dir Directory where the output files from the CCDS comparison step are located.
-ccds_filename_prefix       Prefix used for the output files from the CCDS comparison step.
-ccds_dbname                CCDS database name.
-ccds_host                  CCDS database host.
-ccds_user                  CCDS database user name.
-dna_dbname                 DNA database name.
-dna_host                   DNA database host.
-dna_user                   DNA database user name.
-merge_dbname               Merge database name.
-merge_host                 Merge database host.
-merge_user                 Merge database user name.
-merge_pass                 Merge database user pass.
-assembly_path              Assembly path.
-logic_name                 Logic name to be assigned to the CCDS genes and transcripts copied into the merge database.

=head1 EXAMPLE USAGE

standaloneJob.pl Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCCDSAddition -copy_genes_path $ENSCODE/ensembl-analysis/scripts/genebuild/ -copy_genes_script_name copy_genes.pl -add_ccds_path $ENSCODE/ensembl-personal/genebuilders/scripts/ -add_ccds_script_name add_ccds_support.pl -email report_to@ebi.ac.uk -from me@ebi.ac.uk -ccds_comparison_output_dir OUTPUT_DIR -ccds_filename_prefix CCDS_PREFIX -ccds_dbname CCDSDBNAME -ccds_host CCDSHOST -ccds_user CCDSUSER -dna_dbname DNADBNAME -dna_host DNAHOST -dna_user DNAUSER -merge_dbname MERGEDBNAME -merge_host MERGEHOST -merge_user MERGEUSER -merge_pass MERGEPASS -assembly_path ASSPATH -logic_name ensembl

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveCCDSAddition;

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
      copy_genes_path => undef,
      copy_genes_script_name => undef,
      add_ccds_path => undef,
      add_ccds_script_name => undef,
      email => 'report_to@ebi.ac.uk',
      from => 'your_email@ebi.ac.uk',
      ccds_comparison_output_dir => undef,
      ccds_filename_prefix => undef,
      ccds_dbname => undef,
      ccds_host => undef,
      ccds_user => undef,
      dna_dbname => undef,
      dna_host => undef,
      dna_user => undef,
      merge_dbname => undef,
      merge_host => undef,
      merge_user => undef,
      merge_pass => undef,
      assembly_path => undef,
      logic_name => '', # logic name to be assigned to the missing CCDS genes and their transcripts once they have been copied into the merge db
    }
}

sub fetch_input {
  my $self = shift;

  return 1;
}

sub run {

  my $self = shift;
  
  $self->param_required('copy_genes_path');
  $self->param_required('copy_genes_script_name');
  $self->param_required('add_ccds_path');
  $self->param_required('add_ccds_script_name');
  $self->param_required('email');
  $self->param_required('from');
  $self->param_required('ccds_comparison_output_dir');
  $self->param_required('ccds_filename_prefix');
  $self->param_required('ccds_dbname');
  $self->param_required('ccds_host');
  $self->param_required('ccds_user');
  $self->param_required('dna_dbname');
  $self->param_required('dna_host');
  $self->param_required('dna_user');
  $self->param_required('merge_dbname');
  $self->param_required('merge_host');
  $self->param_required('merge_user');
  $self->param_required('merge_pass');
  $self->param_required('assembly_path');

  #add / at the end of the paths if it cannot be found to avoid possible errors
  if (!($self->param('ccds_comparison_output_dir') =~ /\/$/)) {
    $self->param('ccds_comparison_output_dir',$self->param('ccds_comparison_output_dir')."/");
  }

  # throw exception if ccds_comparison_output_dir does not exist
  if (not -e $self->param('ccds_comparison_output_dir')) {
    run_command("mkdir -p ".$self->param('output_path'),"Create output path.");
  }

  # send report of missing CCDS
  my $missing_ccds_file = emailMissingCCDS($self->param('email'),
                                           $self->param('from'),
                                           $self->param('ccds_comparison_output_dir'),
                                           $self->param('ccds_filename_prefix'),
                                           $self->param('ccds_host'),
                                           $self->param('ccds_dbname'),
                                           $self->param('ccds_user'));

  # copy missing CCDS from CCDS db into merge db
  run_command("perl ".$self->param('copy_genes_path').$self->param('copy_genes_script_name').
              " -sourcedbname ".$self->param('ccds_dbname').
              " -sourcehost ".$self->param('ccds_host').
              " -sourceuser ".$self->param('ccds_user').
              " -dnadbname ".$self->param('dna_dbname').
              " -dnahost ".$self->param('dna_host').
              " -dnauser ".$self->param('dna_user').
              " -outdbname ".$self->param('merge_dbname').
              " -outhost ".$self->param('merge_host').
              " -outuser ".$self->param('merge_user').
              " -outpass ".$self->param('merge_pass').
              " -logic ".$self->param('logic_name').
              " -stable_id".
              " -file $missing_ccds_file".
              " -verbose",
              "Copying CCDS missing genes into merge database...");

  # add CCDS models as supporting features
  run_command("grep -r 'MATCH:' ".$self->param('ccds_comparison_output_dir').$self->param('ccds_filename_prefix')."* | ".
              "awk '{print ".'$5,$8'."}' | ". 
              "sort -u ".
              "> ".$self->param('ccds_comparison_output_dir')."/ccds_id_vs_merge_transcript_stableids.list",
              "Getting the list of CCDS matches with the corresponding transcript stable IDs in preparation for adding the CCDS models as supporting features...");

  run_command("perl ".$self->param('add_ccds_path').$self->param('add_ccds_script_name').
              " -dbname ".$self->param('merge_dbname').
              " -host ".$self->param('merge_host').
              " -port 3306".
              " -user ".$self->param('merge_user').
              " -pass ".$self->param('merge_pass').
              " -path ".$self->param('assembly_path').
              " -file ".$self->param('ccds_comparison_output_dir')."/ccds_id_vs_merge_transcript_stableids.list".
              ">& ".$self->param('ccds_comparison_output_dir')."/transcript_stable_ids_for_ccds_genes.out");

  run_command("tail -1 ".$self->param('ccds_comparison_output_dir')."/transcript_stable_ids_for_ccds_genes.out | awk '{print ".'$NF'."}' ","Checking that all the stable IDs have been processed...",0);

  return 1;
}

                    
sub emailMissingCCDS {
# It parses the output produced by the analysis 'ccds_comparison' in the merge pipeline
# to get the list of missing CCDS models, which is emailed to 'email'.
# It returns the full file path to the file containing the list of missing CCDS.
  my ($email,
      $from,
      $output_dir,
      $ccds_filename_prefix,
      $ccds_host,
      $ccds_name,
      $ccds_user) = @_;

  my $missing_ccds_file = "$output_dir/missing_ccds_models.stable_ids";

  run_command("grep -r 'MATCH:' $output_dir$ccds_filename_prefix* | awk '{print ".'$5'."}' | sort -u > $output_dir/unique_ccds_in_merge_db.list",
              "Getting the list of CCDS matches in preparation for emailing the missing CCDS...");

  run_command("mysql -h $ccds_host -D $ccds_name -u $ccds_user -e 'select stable_id from gene' -BN | sort > $output_dir/ccds_models.list",
              "Getting the list of CCDS stable IDs in preparation for emailing the missing CCDS...");

  run_command("comm -23 $output_dir/ccds_models.list $output_dir/unique_ccds_in_merge_db.list > $missing_ccds_file");
  my $missing_ccds = run_command("cat $missing_ccds_file");

  my $body = <<'END_BODY';
Hi,
 
The merge for the next release has been run and this is the list of missing CCDS models that I got this time.

Regards,
The Ensembl Genebuild team's automatic pipeline powered by eHive

----- missing CCDS models -----

END_BODY

  $body .= $missing_ccds;

  send_email($email,$from,"Missing CCDS models after the merge",$body);

  return $missing_ccds_file;
}

sub write_output {
  my $self = shift;

  return 1;
}

1;
