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

HiveVegaPreparation.pm

=head1 DESCRIPTION

This module automates the preparation steps on the vega database which need to be done before running the Ensembl-HAVANA merge:

Tag read-through transcripts
Methionine to stop codon check
Set ncrna_host attribute
Truncate transcript_supporting_feature table
Add gene attribute to the GAGE cluster (human only)

=head1 OPTIONS

-output_path    Output path where output and log files will be written.

-dbhost         vega database host name

-dbport         vega database port (default 3306)

-dbname         vega database name

-dbuser         vega database username to connect as

-dbpass         vega database password to use

-dnadbhost      dna database host name

-dnadbport      dna database port

-dnadbname      dna database name

-check_vega_met_stop_dir   path to the directory containing check_vega_met_stop.pl

-skip           Number of steps to be skipped. 0 by default.
                For example: if skip=2, the script will begin by running
                the step 3.

-only           Run only the specified step.
                For example, it can be set like "-only 4" to run only
                step number 4.  The 'skip' parameter will be ignored if
                set.

=head1 EXAMPLE USAGE

standaloneJob.pl Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveVegaPreparation -output_path $SCR9/vega_prep.out -dbhost genebuildX -dbname vega_81 -dbuser *** -dbpass *** -dbport 3306 -dnadbhost genebuildX -dnadbport 3306 -dnadbname ensembl_79 -check_vega_met_stop_dir $ENSCODE/ensembl-analysis/scripts/Merge

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveVegaPreparation;

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
      output_path => undef,
      dbhost => undef,
      dbname => undef,
      dbuser => undef,
      dbpass => undef,
      dbport => 3306,
      dnadbhost => undef,
      dnadbport => 3306,
      dnadbname => undef,
      check_vega_met_stop_dir => undef,
      skip => 0,
      only => 0,
    }
}

sub fetch_input {
  my $self = shift;

  return 1;
}

sub run {

  my $self = shift;

  $self->param_required('output_path');
  $self->param_required('check_vega_met_stop_dir');
  $self->param_required('dbhost');
  $self->param_required('dbname');
  $self->param_required('dbuser');
  $self->param_required('dnadbhost');
  $self->param_required('dnadbname');

  #add / at the end of the paths if it cannot be found to avoid possible errors
  if (!($self->param('output_path') =~ /\/$/)) {
    $self->param('output_path',$self->param('output_path')."/");
  }

  if (!($self->param('check_vega_met_stop_dir') =~ /\/$/)) {
    $self->param('check_vega_met_stop_dir',$self->param('check_vega_met_stop_dir')."/");
  }

  # if parameter 'only' is set, make sure that we 'skip' anything else
  if ($self->param('only') > 0) {
    $self->param('skip',99999); # hopefully, we wont have 99999 steps in a genebuild
  }

  # create output dir if it does not exist
  if (not -e $self->param('output_path')) {
    run_command("mkdir -p ".$self->param('output_path'),"Create output path.",0);
  }

  # print the number of skipped steps
  print("The first step will be skipped.\n") if ($self->param('skip') == 1 and $self->param('only') == 0);
  print("The ".$self->param('skip')." first steps will be skipped.\n") if ($self->param('skip') > 1 and $self->param('only') == 0);

  # print the number of the only step executed
  print("Only the step number ".$self->param('only')." will be executed.\n") if ($self->param('only') > 0);

  $self->readthrough_transcripts_tagged($self->param('dbhost'),
                                 $self->param('dbport'),
                                 $self->param('dbuser'),
                                 $self->param('dbpass'),
                                 $self->param('dbname')) if ($self->param('skip') < 1 or $self->param('only') == 1);

  $self->methionine_to_stop_codon($self->param('dbhost'),
                           $self->param('dbport'),
                           $self->param('dbuser'),
                           $self->param('dbpass'),
                           $self->param('dbname'),
                           $self->param('dnadbhost'),
                           $self->param('dnadbport'),
                           $self->param('dnadbname'),
                           $self->param('check_vega_met_stop_dir'),
                           $self->param('output_path')) if ($self->param('skip') < 2 or $self->param('only') == 2);

  $self->set_ncrna_host_gene_attribute($self->param('dbhost'),
                                $self->param('dbport'),
                                $self->param('dbuser'),
                                $self->param('dbpass'),
                                $self->param('dbname')) if ($self->param('skip') < 3 or $self->param('only') == 3);

  $self->truncate_tsf_table($self->param('dbhost'),
                     $self->param('dbport'),
                     $self->param('dbuser'),
                     $self->param('dbpass'),
                     $self->param('dbname')) if ($self->param('skip') < 4 or $self->param('only') == 4);

  $self->add_attribute_to_GAGE_cluster($self->param('dbhost'),
                                $self->param('dbport'),
                                $self->param('dbuser'),
                                $self->param('dbpass'),
                                $self->param('dbname')) if ($self->param('skip') < 5 or $self->param('only') == 5);
  return 1;
}

sub write_output {
  my $self = shift;

  return 1;
}

sub readthrough_transcripts_tagged {
  my ($self, $dbhost,$dbport,$dbuser,$dbpass,$dbname) = @_;

  my $readthrough_transcripts_sql = "update transcript_attrib ta, transcript t set ta.attrib_type_id = (select attrib_type_id from attrib_type where code = 'readthrough_tra') where ta.value like '%read%through%' and ta.transcript_id = t.transcript_id and ta.value not like '%_stop_codon_rt%';";

  run_command("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e\"$readthrough_transcripts_sql\"",
              "\n===== STEP 1: Readthrough transcripts tagged =====\n");

  my $num_readthrough_transcripts_sql = "select count(*) from transcript_attrib ta,transcript t, attrib_type at where code = 'readthrough_tra' and (ta.value like '%read%through%') and ta.transcript_id = t.transcript_id;";

  my $num_updated = run_command("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e\"$num_readthrough_transcripts_sql\"");

  print("\nThere are ".int($num_updated)." readthrough transcripts tagged now.\n");
}

sub methionine_to_stop_codon {
  my ($self, $dbhost,$dbport,$dbuser,$dbpass,$dbname,$dnadbhost,$dnadbport,$dnadbname,$check_vega_met_stop_dir,$output_path) = @_;

  my $ids_file = "$output_path/havana_coding_transcript_ids.txt";
  my $full_ids_file = "$output_path/full_length_havana_coding_transcript_ids.txt";

  my $dump_ids_sql = "select distinct(tr.transcript_id) from translation tr, transcript t, transcript_attrib ta, gene g where t.gene_id = g.gene_id and t.transcript_id = ta.transcript_id and tr.transcript_id = t.transcript_id;";
  run_command("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e\"$dump_ids_sql\" > $ids_file",
              "\n===== STEP 2: Methionine to stop codon =====\n");

  my $num_ids = int(run_command("wc -l $ids_file | cut -f1 -d' '"));

  run_command("perl $check_vega_met_stop_dir/check_vega_met_stop.pl -dbname $dbname -host $dbhost -port $dbport -dna_dbname $dnadbname -dna_host $dnadbhost -dna_port $dnadbport < $ids_file > $full_ids_file");

  my $num_wrong_ids = int(
     run_command(
       "awk '\$3 != \"ok\" { nok++ } END { print nok }' $full_ids_file")
  );
  my $num_correct_ids = $num_ids - $num_wrong_ids;

  # Check that $num_wrong_ids is similar to the previous release one (if previous release exists).
  $self->warning("Methionine to stop codon check: Please check manually that the number of wrong transcripts $num_wrong_ids is similar to the one in previous release (if exists). Some of the previous releases numbers are: e71: 29260; e70: 29381 e69: 56594; e68: 53703");
}

sub set_ncrna_host_gene_attribute {
  my ($self, $dbhost,$dbport,$dbuser,$dbpass,$dbname) = @_;

  my $ncrna_host_sql = "update gene_attrib ga, gene g set ga.attrib_type_id = (select attrib_type_id from attrib_type where code = 'ncrna_host') where ga.gene_id = g.gene_id and (value like 'transcribed%' or value = 'ncrna_host');";
  my $check_sql_1 = "select count(*) from gene_attrib where (value like 'transcribed%' or value = 'ncrna_host');";
  my $check_sql_2 = "select count(*) from gene g, gene_attrib ga where ga.attrib_type_id = (select attrib_type_id from attrib_type where code = 'ncrna_host') and ga.gene_id = g.gene_id;";

  run_command("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e\"$ncrna_host_sql\"",
              "\n===== STEP 3: Set ncrna_host gene attribute =====\n");

  run_command("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e\"$check_sql_2\"",
              "",
              run_command("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e\"$check_sql_1\"",
                          "\nChecking that the number of updated attrib_type_id matches the number of values that were found.\n")
             );
}

sub truncate_tsf_table {
  my ($self, $dbhost,$dbport,$dbuser,$dbpass,$dbname) = @_;

  my $tsf_bak_name = "tsf_bak_".time();

  my $num_tsf = run_command("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e'select count(*) from transcript_supporting_feature;'",
                            "\n===== STEP 4: Truncate transcript_supporting_feature table =====\n");

  run_command("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e'rename table transcript_supporting_feature to $tsf_bak_name; create table transcript_supporting_feature like $tsf_bak_name;'");

  run_command("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e'select count(*) from $tsf_bak_name;'","",$num_tsf);

  run_command("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e'select count(*) from transcript_supporting_feature;'","",0);

}

sub add_attribute_to_GAGE_cluster {
  my ($self, $dbhost,$dbport,$dbuser,$dbpass,$dbname) = @_;


  my $num_att = run_command("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e\"select count(*) from transcript,gene_attrib ga where value = 'gene_cluster_GAGE' and transcript.gene_id=ga.gene_id;\"",
                            "\n===== STEP 5: Add attribute to GAGE cluster =====\n");

  my $insert_sql = <<END;
insert ignore into transcript_attrib (transcript_id,attrib_type_id,value)
select t.transcript_id,(select attrib_type_id from attrib_type where code = 'gene_cluster'),'gene_cluster_GAGE'
from transcript t, gene_attrib ga, attrib_type at
where t.gene_id = ga.gene_id and ga.attrib_type_id = at.attrib_type_id and ga.value = 'gene_cluster_GAGE';
END

  run_command("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e\"$insert_sql\"");

  my $check_sql = "select count(*) from transcript_attrib where attrib_type_id=(select attrib_type_id from attrib_type where code = 'gene_cluster');";
  run_command("mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e\"$check_sql\"",
              "",
              $num_att);
}

1;
