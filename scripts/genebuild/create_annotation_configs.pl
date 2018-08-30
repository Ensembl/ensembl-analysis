#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

use warnings;
use strict;
use feature 'say';

use Getopt::Long qw(:config no_ignore_case);
use File::Spec::Functions;
use File::Basename;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor;
use Net::FTP;
use Data::Dumper;


my $config_file;
my $config_only = 0;

my $base_guihive = 'http://guihive.ebi.ac.uk:8080';
my $ftphost = "ftp.ncbi.nlm.nih.gov";
my $ftpuser = "anonymous";
my $ftppassword = "";
my $assembly_registry_host = $ENV{GBS5};
my $assembly_registry_port = $ENV{GBP5};


GetOptions('config_file:s' => \$config_file,
           'config_only!'  => \$config_only,
           'assembly_registry_host:s' => \$assembly_registry_host,
           'assembly_registry_port:s' => \$assembly_registry_port,
          );

unless(-e $config_file) {
  die "Could not find the config file. Path used:\n".$config_file;
}

my $assembly_registry = new Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor(
  -host    => $assembly_registry_host,
  -port    => $assembly_registry_port,
  -user    => 'ensro',
  -dbname  => 'do1_stable_id_space_assembly_registry');

my $ncbi_taxonomy = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => 4240,
  -user    => 'ensro',
  -host    => 'mysql-ensembl-mirror',
  -dbname  => 'ncbi_taxonomy');

my $general_hash = {};
my $assembly_accessions;

open(IN,$config_file) || die("Could not open $config_file");
while(<IN>) {
    my $line = $_;

    $line =~ s/\s//g;
    if($line =~ /^(.+)\=(.+)/) {
      my $key = $1;
      my $value = $2;
      say "Found key/value pair: ".$key." => ".$value;

      # Note that in the ini file the key for write user is user_w for clarity, but in the actual
      # hive config it's user (since hive expects a user key). This is just a substitution to account
      # for this
      if($key eq 'user_w') {
        $key = 'user';
      }

      $general_hash->{$key} = $value;
    }elsif($line eq "\n") {
	# Skip
    }else {
      say "Line format not recognised. Skipping line:\n".$line;
    }
}
close IN || die("Could not close $config_file");

$assembly_accessions = $general_hash->{'assembly_accessions'};
$assembly_accessions =~ /\[(.+)\]/;
my $accession_string = $1;
$accession_string =~ s/ //g;
my @accession_array = split(',',$accession_string);
unless(scalar(@accession_array)) {
  die "Issue parsing assembly_accessions line. Format expected:\n".
      "assembly_accessions=[GCA_000952055.2,GCA_000164805.2,GCA_001604975.1]";
}

unless($general_hash->{'output_path'}) {
  die "Could not find an output path setting in the config. Expected setting".
      "output_path=/path/to/output/dir/";
}

my $output_path = $general_hash->{'output_path'};
unless(-e $output_path) {
  system("mkdir -p ".$general_hash->{'output_path'});
}

my $clade = $general_hash->{'clade'};
unless($clade) {
  die "No clade selected. Need a clade to specify parameters. Format expected:\n".
      "clade=primates";
}

my $clade_hash = clade_settings($clade);
foreach my $key (keys(%{$clade_hash})) {
  $general_hash->{$key} = $clade_hash->{$key};
}

my $ftp_base_dir = '/genomes/all/';

unless($config_only) {
  my $unique_file_name = basename($config_file).".cmds";
  open(LOOP_CMD,">".$general_hash->{'output_path'}."/".$unique_file_name) || die("Could not open $unique_file_name");
}

foreach my $accession (@accession_array) {
  my $assembly_hash = {};

  say "Processing accession: ".$accession;
  unless($accession =~ /GCA_([\d]{3})([\d]{3})([\d]{3})\.\d+/) {
    die "Found an assembly accession that did not match the regex. Offending accession: ".$accession;
  }


  # Get stable id prefix
  my $stable_id_prefix = $assembly_registry->fetch_stable_id_prefix_by_gca($accession);
  say "Fetched the following stable id prefix for ".$accession.": ".$stable_id_prefix;
  $assembly_hash->{'stable_id_prefix'} = $stable_id_prefix;

  # Get stable id start
  my $stable_id_start = $assembly_registry->fetch_stable_id_start_by_gca($accession);
  say "Fetched the following stable id start for ".$accession.": ".$stable_id_start;
  $assembly_hash->{'stable_id_start'} = $stable_id_start;

  my $assembly_ftp_path = $ftp_base_dir.'GCA/'.$1.'/'.$2.'/'.$3.'/';
  my $full_assembly_path;
  my $assembly_name;

  my $ftp = Net::FTP->new($ftphost) or die "Can't open $ftphost\n";
  $ftp->login($ftpuser, $ftppassword) or die "Can't log $ftpuser in\n";
  $ftp->cwd($assembly_ftp_path);
  my @ftp_dir_contents = $ftp->ls;
  foreach my $entry (@ftp_dir_contents) {
    if($entry =~ /^$accession\_(.+)$/) {
      $full_assembly_path = $assembly_ftp_path."/".$&."/";
      $assembly_name = $1;
    }
  }

  unless($full_assembly_path && $assembly_name) {
    die "Issue finding ftp path for the following GCA: $accession $assembly_ftp_path";
  }

  say "Setting assembly name for ".$accession." to ".$assembly_name;

  $assembly_hash->{'assembly_accession'} = $accession;
  $assembly_hash->{'assembly_name'} = $assembly_name;

  parse_assembly_report($ftp,$general_hash,$assembly_hash,$accession,$assembly_name,$full_assembly_path,$output_path);

  # Add in the species url by uppercasing the first letter of the species name
  my $species_url = $assembly_hash->{'species_name'};
  $species_url = ucfirst($species_url);
  $assembly_hash->{'species_url'} = $species_url;

  # Get repeatmodeler library path if one exists
  my $repeatmodeler_file = catfile($ENV{REPEATMODELER_DIR},'species',$assembly_hash->{'species_name'},$assembly_hash->{'species_name'}.'.repeatmodeler.fa');
  if(-e $repeatmodeler_file) {
    say "Found the following repeatmodeler file for the species:\n".$repeatmodeler_file;
    $assembly_hash->{'repeatmodeler_library'} = $repeatmodeler_file;
  } else {
    say "Did not find an repeatmodeler species library for ".$assembly_hash->{'species_name'}." on path:\n".$assembly_hash->{'species_name'};
  }

  create_config($assembly_hash);

  unless($config_only) {
    chdir($assembly_hash->{'output_path'});
    my $cmd = "init_pipeline.pl Genome_annotation_conf.pm -hive_force_init 1";
    my $result = `$cmd`;
    unless($result =~ /beekeeper.+\-sync/) {
      die "Failed to run init_pipeline for ".$assembly_hash->{'species_name'}."\nCommandline used:\n".$cmd;
    }

    my $sync_command = $&;
    my $return = system($sync_command);
    if($return) {
      die "Failed to sync the pipeline for ".$assembly_hash->{'species_name'}."\nCommandline used:\n".$cmd;
    }

    my $loop_command = $sync_command;
    $loop_command =~ s/sync/loop \-sleep 0.3/;
    my ($ehive_url) = $sync_command =~ /url\s+(\S+)/;
    my ($driver, $user, $password, $host, $port, $dbname) = $ehive_url =~ /^(\w+)[:\/]+(\w*):(\w+)@([^:]+):(\d+)\/(\w+)$/;
    if ($password) {
      $password = '&passwd=xxxxx';
    }
    say LOOP_CMD "#GuiHive: $base_guihive/?driver=$driver&username=$user&host=$host&port=$port&dbname=$dbname$password";
    say LOOP_CMD "export EHIVE_URL=$ehive_url";
    say LOOP_CMD $loop_command;
  }
}

unless($config_only) {
  close(LOOP_CMD) || die("Could not close the cmd file");
}

exit;

sub parse_assembly_report {
  my ($ftp,$general_hash,$assembly_hash,$accession,$assembly_name,$full_assembly_path,$output_path) = @_;

  $ftp->cwd();
  $ftp->cwd($full_assembly_path);

  my $report_file_name = $accession."_".$assembly_name."_assembly_report.txt";
  my $report_file_content;
  my $report_file_handle;
  my $taxon_id;
  my $species_name;
  my $refseq_accession;
  my $assembly_level;
  my $wgs_id;
  open($report_file_handle, '>', \$report_file_content) || die("could not open $report_file_name");

  unless($ftp->get($report_file_name, $report_file_handle)) {
    die "Failed to retrieve the assembly report file: ", $ftp->message;
  }

  my @report_file_content = split("\n",$report_file_content);
  foreach my $line (@report_file_content) {
    unless($line =~ /^\# /) {
      next;
    }

    if($line =~ /Taxid\:\s*(\d+)/) {
      $taxon_id = $1;
    } elsif($line =~ /Assembly level\:\s*([a-zA-Z]+)/) {
      $assembly_level = $1;
    } elsif($line =~ /RefSeq assembly accession\:\s*(GCF\_\d{9}\.\d+)/){
      $refseq_accession = $1;
    } elsif($line =~ /WGS project\:\s*([A-Z]{4}\d{2})/) {
      $wgs_id = $1;
    }

  }

  unless($taxon_id && $assembly_level) {
    die "Failed to fully parse the assembly report file";
  }

  if(exists($general_hash->{'load_toplevel_only'}) && $general_hash->{'load_toplevel_only'} == 0) {
    unless($wgs_id) {
      die "Need the wgs id as the load_toplevel_only flag was set to 0. Failed to find the id in the report file";
    }
  }

  my $sth = $ncbi_taxonomy->dbc->prepare("SELECT name from ncbi_taxa_name where taxon_id=? and name_class='scientific name'");
  $sth->bind_param(1,$taxon_id);
  $sth->execute;
  $species_name = $sth->fetchrow_array();
  $species_name =~ s/\s+/\_/g;
  $species_name = lc($species_name);
  unless($species_name) {
    die "Was not able to retrieve the species name from the NCBI taxonomy db using the taxon id. Taxon id used: ".$taxon_id;
  }

  $assembly_hash->{'taxon_id'} = $taxon_id;

  if($refseq_accession) {
    $assembly_hash->{'assembly_refseq_accession'} = $refseq_accession;
  } else {
    say "Found no RefSeq accession for this assembly";
  }

  $assembly_hash->{'assembly_level'} = $assembly_level;
  $assembly_hash->{'wgs_id'} = $wgs_id;
  $assembly_hash->{'species_name'} = $species_name;
  $assembly_hash->{'production_name'} = $species_name;
  @{$assembly_hash}{keys(%$general_hash)} = values(%$general_hash);
  $assembly_hash->{'output_path'} .= "/".$species_name."/".$accession."/";
}


sub create_config {
  my ($assembly_hash) = @_;

  foreach my $key (keys(%{$assembly_hash})) {
    say "ASSEMBLY: ". $key." => ".$assembly_hash->{$key};
  }

  my $output_path = $assembly_hash->{'output_path'};
  unless(-e $output_path) {
   system('mkdir -p '.$output_path);
  }

  my $config_string = "";
  my $past_default_options = 0;
  open(CONFIG,$ENV{ENSCODE}."/ensembl-analysis/modules/Bio/EnsEMBL/Analysis/Hive/Config/Genome_annotation_conf.pm") || die("Could not open the config file");
  while(my $line = <CONFIG>) {
    if($line =~ /sub pipeline_create_commands/) {
      $past_default_options = 1;
    } elsif($past_default_options == 0) {
      if($line =~ /SUPER\:\:default_options/) {
        if($assembly_hash->{'dbowner'}) {
          $line .= "'dbowner' => '".$assembly_hash->{'dbowner'}."',";
        }
      }
      if($line =~ /\'([^\']+)\'\s*\=\>\s*('[^\']*\')/) {
        my $conf_key = $1;
        my $conf_val = $2;
        if($assembly_hash->{$conf_key}) {
          my $sub_val = "'".$assembly_hash->{$conf_key}."'";
          $line =~ s/$conf_val/$sub_val/;
        }
      }
    }

    $config_string .= $line;
  }
  close(CONFIG) || die("CouLd not close the config file");

  $config_string =~ s/package Bio::EnsEMBL::Analysis::Hive::Config::([^;]+;)/package $1/;
  open(OUT_CONFIG,">".$output_path."/Genome_annotation_conf.pm") || die("Could not open the config file for writing");
  print OUT_CONFIG $config_string;
  close OUT_CONFIG || die("Could not close the config file for writing");
}

sub clade_settings {
  my ($clade) = @_;
  my $clade_settings = {
    'primates' => {
      'repbase_library'    => 'primates',
      'repbase_logic_name' => 'primates',
      'uniprot_set'        => 'primates_basic',
    },

    'rodents' => {
      'repbase_library'    => 'rodents',
      'repbase_logic_name' => 'rodents',
      'uniprot_set'        => 'rodents_basic',
    },

    'mammals' => {
      'repbase_library'    => 'mammals',
      'repbase_logic_name' => 'mammals',
      'uniprot_set'        => 'mammals_basic',
    },

    'fish_teleost' => {
      'repbase_library'     => 'Teleostei',
      'repbase_logic_name'  => 'teleost',
      'uniprot_set'         => 'fish_basic',
      'ig_tr_fasta_file'    => 'fish_ig_tr.fa',
      'masking_timer_long'  => '6h',
      'masking_timer_short' => '3h',
    },

    'distant_vertebrate' => {
      'repbase_library'    => 'vertebrates',
      'repbase_logic_name' => 'vertebrates',
      'uniprot_set'        => 'distant_vertebrate',
      'masking_timer_long'  => '6h',
      'masking_timer_short' => '3h',
    },

  };

  unless($clade_settings->{$clade}) {
    die "Could not find the clade specified in the clade_settings hash. Clade specified: ".$clade;
  }

  return($clade_settings->{$clade});
}

