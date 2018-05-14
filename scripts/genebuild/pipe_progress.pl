use warnings;
use strict;
use feature 'say';

use Getopt::Long qw(:config no_ignore_case);
use File::Spec::Functions;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor;
use Net::FTP;
use Data::Dumper;


my $db_pattern;
my $pipeline_type = 'main';
my $db_hash;
GetOptions('db_pattern:s'    => \$db_pattern,
           'pipeline_type:s' => \$pipeline_type,);

my $db_port = 4527;
my $matched_db_list = "";
for(my $db_number = 1; $db_number <= 7; $db_number++) {
  my $query_db = "mysql -uensro -hmysql-ens-genebuild-prod-".$db_number.".ebi.ac.uk -P".$db_port." -NB -e \"show databases like '".$db_pattern."%'\"";
  my $matched_db_list .= `$query_db`;
  my @matched_dbs = split(/\n/,$matched_db_list);
  foreach my $matched_db (@matched_dbs) {
    say $matched_db.":";
    assess_db($matched_db,$db_number,$db_port,$pipeline_type);
  }
  $db_port++;
}

exit;

sub assess_db {
  my ($db_name,$db_number,$db_port,$pipeline_type) = @_;

  my $query_base = "mysql -uensro -hmysql-ens-genebuild-prod-".$db_number.".ebi.ac.uk -P".$db_port." ".$db_name." -NB -e ";
  my @analyses = @{checkpoint_analyses($pipeline_type)};
  my $analyses_count = scalar(@analyses);
  unless($analyses_count) {
    say "No checkpoint reached"
  }
  for(my $i=0; $i < $analyses_count; $i++) {
    my $analysis = $analyses[$i];
    my ($logic_name,$message) = split(/\:/,$analysis);
    my $analysis_query = $query_base."\"select status from job join analysis_base using(analysis_id) where logic_name='".$logic_name."'\"";
    my $status = `$analysis_query`;
    chomp($status);
    if($status eq "DONE") {
      my $checkpoint_index = $analyses_count - $i;
      say "Checkpoint ".$checkpoint_index." of ".$analyses_count.": ".$message." (".$logic_name.")";
      last;
    }
  }

#  my $job_status_query = $query_base."\"select count(*),status from job where status in ('FAILED','READY','RUN') group by status order by status\"";
  my $job_status_query = $query_base."\"select count(*),status from job group by status order by status\"";
  my $jobs = `$job_status_query`;

  say "Status of jobs:";
  print $jobs;
  if($jobs =~ /FAILED/) {
    say "Warning: found failed jobs in the pipeline db!!!!!!!!!";
    say "Failed analyses:";
    my $failed_analyses_query = $query_base."\"select logic_name,count(*) from job join analysis_base using(analysis_id) where ".
                                "status='FAILED' group by logic_name\"";
    my $failed_analyses = `$failed_analyses_query`;
    print $failed_analyses;
  }

  my $last_completed_job_query = $query_base."\"select when_completed,logic_name from job join analysis_base using(analysis_id) where ".
                                 "when_completed=(select max(when_completed) from job) limit 1\"";
  my $last_completed = `$last_completed_job_query`;
  print "Last job completed: ".$last_completed;

  my $time_diff_query = $query_base."\"select time_to_sec(timediff(now(), max(when_completed) )) / 3600 from job\"";
  my $time_diff = `$time_diff_query`;
  chomp($time_diff);
  if($time_diff >= 5) {
    say "Warning: ".$time_diff." hours have passed since a job completed!!!!!!!!!!";
  }

  my $highest_active_analysis_query = $query_base."\"select distinct(logic_name) from job join analysis_base using(analysis_id) where ".
                                      "analysis_id=(select max(analysis_id) from job where status in ('RUN'))\"";
  my $highest_active_analysis = `$highest_active_analysis_query`;
  chomp $highest_active_analysis;
  my $highest_aa_breakdown_query = $query_base."\"select count(*),status from job join analysis_base using(analysis_id) where ".
                                   "logic_name='".$highest_active_analysis."' group by status\"";
  my $highest_aa = `$highest_aa_breakdown_query`;
  my @highest_aa_array = split('\n',$highest_aa);
  my $highest_aa_breakdown = "";
  foreach my $row (@highest_aa_array) {
    my ($count,$status) = split('\t',$row);
    $highest_aa_breakdown .= " ".$count." ".$status.",";
  }
  chop $highest_aa_breakdown;
  say "Highest analysis id active: ".$highest_active_analysis.$highest_aa_breakdown."\n";
}

sub checkpoint_analyses {
  my ($pipeline_type) = @_;

  if($pipeline_type eq 'main') {
    return([
             'core_healthchecks:core healthchecks complete',
             'clean_unused_analyses:optimise complete',
             'set_meta_coords:Genes copied to core',
             'final_db_sanity_checks:finished making final db',
             'create_final_geneset_db:finished pseudogenes',
             'create_pseudogene_db:starting pseudogenes',
             'layer_annotation_sanity_checks:finished layering/genebuilder phase',
             'create_layering_output_db:entering layering/genebuilder phase',
             'copy_rnaseq_blast_db:preparing the rnaseq db',
             'cluster_ig_tr_genes:projection complete',
             'create_projection_coding_input_ids:projection started',
             'create_projection_coding_db:ncRNAs complete',
             'fan_ncrna:starting ncRNAs',
             'ig_tr_sanity_checks:IG/TR complete',
             'create_ig_tr_db:cdnas complete',
             'create_cdna_db:starting cdnas',
             'classify_genblast_models:genblast complete',
             'download_uniprot_files:about to start genblast',
             'genome_prep_sanity_checks:genome prep complete',
             'load_toplevel_sequences:toplevel loaded',
             'create_core_db:empty core db created',
          ]);
  }

}
