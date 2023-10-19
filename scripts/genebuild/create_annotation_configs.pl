#!/usr/bin/env perl

# Copyright [2017-2022] EMBL-European Bioinformatics Institute
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
use File::Spec::Functions qw(catfile splitdir catdir updir devnull);
use File::Basename;
use File::Copy;
use Bio::EnsEMBL::Utils::Exception qw (warning throw);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor;
use Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor;
use Net::FTP;
use Cwd qw(realpath chdir getcwd);
use Data::Dumper;
use DateTime;
use JSON;

my $config_file;
my $config_only = 0;
my $custom_load = 0;
my $early_load = 0;

my $total_running_workers_max = 2000;
my $base_guihive = 'http://guihive.ebi.ac.uk:8080';
my $ftphost = "ftp.ncbi.nlm.nih.gov";
my $ftpuser = "anonymous";
my $ftppassword = "";
my $assembly_registry_host = $ENV{GBS1};
my $assembly_registry_port = $ENV{GBP1};
my $force_init = 0;
my $check_for_transcriptomic = 0;
my $selected_db;
my $current_genebuild = 0;

### change to to 1 for non-verts
my $is_non_vert = 0;


GetOptions('config_file:s' => \$config_file,
           'config_only!'  => \$config_only,
           'custom_load!'  => \$custom_load,
           'early_load!'   => \$early_load,
           'assembly_registry_host:s' => \$assembly_registry_host,
           'assembly_registry_port:s' => \$assembly_registry_port,
           'force_init!' => \$force_init,
           'is_non_vert!' => \$is_non_vert,
           'check_for_transcriptomic!' => \$check_for_transcriptomic,
           'current_genebuild!' => \$current_genebuild,
          );

unless(-e $config_file) {
  throw("Could not find the config file. Path used:\n".$config_file);
}

my $general_hash = {};

if ($is_non_vert == 1) {
  $general_hash->{'replace_repbase_with_red_to_mask'} = '1';
  $general_hash->{'skip_projection'} = '1';
  $general_hash->{'is_non_vert'} = '1';
  $general_hash->{'protein_blast_db_file'} = 'PE12';
  $general_hash->{'protein_entry_loc_file'} = 'entry_loc';
}

$selected_db = "gb_assembly_registry";

my $taxonomy_adaptor = new Bio::EnsEMBL::Taxonomy::DBSQL::TaxonomyDBAdaptor(
  -host    => 'mysql-ens-meta-prod-1',
  -port    => 4483,
  -user    => 'ensro',
  -dbname  => 'ncbi_taxonomy');

my $ncbi_taxonomy = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => 4483,
  -user    => 'ensro',
  -host    => 'mysql-ens-meta-prod-1',
  -dbname  => 'ncbi_taxonomy');

open(IN,$config_file) || throw("Could not open $config_file");
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
        $general_hash->{$key} = $value;
      }
      if($key eq 'password') {
        $general_hash->{$key} = $value;
      }
      #Ignore clade settings from .ini file if set
      if($key eq 'clade' && !$custom_load && !$early_load) {
       # $key = 'user';
       say "Warning: Ignoring clade value from .ini file";
      }
      else{
        $general_hash->{$key} = $value;
      }
    }elsif($line eq "\n") {
	# Skip
    }else {
      say "Line format not recognised. Skipping line:\n".$line;
    }
}
close IN || throw("Could not close $config_file");

# If replace_repbase_with_red_to_mask is true then force first_choice_repeat to be red logic name
if (exists $general_hash->{'replace_repbase_with_red_to_mask'} and $general_hash->{'replace_repbase_with_red_to_mask'} and !exists $general_hash->{first_choice_repeat}) {
  $general_hash->{first_choice_repeat} = $general_hash->{red_logic_name} || 'repeatdetector';
}

 my $assembly_registry = new Bio::EnsEMBL::Analysis::Hive::DBSQL::AssemblyRegistryAdaptor(
  -host    => $assembly_registry_host,
  -port    => $assembly_registry_port,
  -user    => $general_hash->{'user'},
  -pass    => $general_hash->{'password'},
  -dbname  => $selected_db
  );

#Adding registry details to hash for populating main config
$general_hash->{registry_host} = $assembly_registry_host;
$general_hash->{registry_port} = $assembly_registry_port;
$general_hash->{registry_db} = $assembly_registry->{_dbc}->{_dbname};

unless($general_hash->{'output_path'}) {
  throw("Could not find an output path setting in the config. Expected setting".
        "output_path=/path/to/output/dir/");
}

my $output_path = $general_hash->{'output_path'};
unless(-e $output_path) {
  system("mkdir -p ".$general_hash->{'output_path'});
}

my $fh;
unless($config_only) {
  my $unique_file_name = basename($config_file).".cmds";
  open($fh,">".$general_hash->{'output_path'}."/".$unique_file_name) || throw("Could not open $unique_file_name");
}

$total_running_workers_max = $general_hash->{total_running_workers_max} if (exists $general_hash->{total_running_workers_max});
my @path = splitdir(realpath($0));
while (my ($dir) = pop @path) {
  if ($dir eq 'ensembl-analysis') {
    last;
  }
}
my $enscode_directory = catdir(@path);
my $hive_directory;
if (-d catdir($enscode_directory, 'ensembl-hive')) {
  $hive_directory = catdir($enscode_directory, 'ensembl-hive');
}
elsif (-d catdir($enscode_directory, updir, 'ensembl-hive')) {
  $hive_directory = catdir($enscode_directory, updir, 'ensembl-hive');
}
if ($hive_directory) {
  print "Your hive directory is $hive_directory\n";
  my $json = '';
  open(JH, catfile($hive_directory, 'hive_config.json')) || throw('Could not open '.catfile($hive_directory, 'hive_config.json'));
  while(<JH>) {
    chomp;
    $json .= $_;
  }
  close(JH) || throw('Could not open '.catfile($hive_directory, 'hive_config.json'));
  my $json_decoder = JSON->new->relaxed;
  my $json_data = $json_decoder->decode($json);
  if (exists $json_data->{Meadow}->{LSF}->{TotalRunningWorkersMax}) {
    $json_data->{Meadow}->{LSF}->{TotalRunningWorkersMax} = $total_running_workers_max;
  }
  if (exists $json_data->{Meadow}->{LSF}->{EBI}->{TotalRunningWorkersMax}) {
    $json_data->{Meadow}->{LSF}->{EBI}->{TotalRunningWorkersMax} = $total_running_workers_max;
  }
  open(JH, '>'.catfile($hive_directory, 'hive_config.json')) || throw('Could not open '.catfile($hive_directory, 'hive_config.json'));
  $json_decoder->pretty(1);
  print JH $json_decoder->encode($json_data);
  close(JH) || throw('Could not open '.catfile($hive_directory, 'hive_config.json'));
}
else {
  print "Cannot change the TotalRunningWorkersMax in you hive_config.json\n";
}

assign_server_info($general_hash);

# If someone is trying to custom load an assembly, then we need to do things quite differently
# Also since the keys basically need to be manually entered, this will only run on one species at a time
if($custom_load) {
  custom_load_data($general_hash,$ncbi_taxonomy);
  create_config($general_hash);
  unless($config_only) {
    init_pipeline($general_hash,$hive_directory,$force_init,$fh);
    close($fh) || throw("Could not close the cmd file");
  }
  exit;
}

# If we get here, it means standard loading using assembly report files for any GCAs provided
my $ftp_base_dir = '/genomes/all/';
my $assembly_accessions = $general_hash->{'assembly_accessions'};
$assembly_accessions =~ /\[(.+)\]/;
my $accession_string = $1;
$accession_string =~ s/ //g;
my @accession_array = split(',',$accession_string);
unless(scalar(@accession_array)) {
  throw("Issue parsing assembly_accessions line. Format expected:\n".
        "assembly_accessions=[GCA_000952055.2,GCA_000164805.2,GCA_001604975.1]");
}

foreach my $accession (@accession_array) {
  my $assembly_hash = {};

  say "Processing accession: ".$accession;
  unless($accession =~ /GCA_([\d]{3})([\d]{3})([\d]{3})\.\d+/) {
    throw("Found an assembly accession that did not match the regex. Offending accession: ".$accession);
  }
  
  #Check genebuild status of assembly
  if ($current_genebuild == 1) {
  } else {
    check_annotation_status($accession);
  }

  # Get stable id prefix
  my $stable_id_prefix;
  if ($general_hash->{'stable_id_prefix'}){
    $stable_id_prefix = $general_hash->{'stable_id_prefix'};
  }
  else {
    $stable_id_prefix = $assembly_registry->fetch_stable_id_prefix_by_gca($accession);
  }
  say "Fetched the following stable id prefix for ".$accession.": ".$stable_id_prefix;
  $assembly_hash->{'stable_id_prefix'} = $stable_id_prefix;

  #Get clade
  my $clade;
  if ($general_hash->{'clade'}) {
    $clade = $general_hash->{'clade'};
  }
  else {
    $clade = $assembly_registry->fetch_clade_by_gca($accession);
  }
  say "Fetched the following clade for ".$accession.": ".$clade;
  #      #Note: this is to assign repeat library settings for clades that do not have defined settings yet
  if (($clade eq 'amphibians') || ($clade eq 'sharks') || ($clade eq 'vertebrates')){
    $clade = 'distant_vertebrate';
  }
  $assembly_hash->{'clade'} = $clade;
  #set a db for validating models from transcriptomic data
  if ($clade eq 'mammals'){
      $general_hash->{'protein_blast_db_file'} = 'uniprot_mammalia_sp';
  }
  else{
      $general_hash->{'protein_blast_db_file'} = 'uniprot_vertebrata_sp';
  }
  
  # Get stable id start
  my $stable_id_start;
  if (exists ($general_hash->{'stable_id_start'}) && $general_hash->{'stable_id_start'} >=0) {
    $stable_id_start = $general_hash->{'stable_id_start'};
  }
  else {
    $stable_id_start = $assembly_registry->fetch_stable_id_start_by_gca($accession);
  }
 unless (defined($stable_id_start)) {
   throw ("Could not find stable id start");
 }
  say "Fetched the following stable id start for ".$accession.": ".$stable_id_start;
  $assembly_hash->{'stable_id_start'} = $stable_id_start;

  my $assembly_ftp_path = $ftp_base_dir.'GCA/'.$1.'/'.$2.'/'.$3.'/';
  my $full_assembly_path;
  my $assembly_name;

  my $ftp = Net::FTP->new($ftphost) or throw("Can't open $ftphost");
  $ftp->login($ftpuser, $ftppassword) or throw("Can't log $ftpuser in");
  $ftp->cwd($assembly_ftp_path);
  my @ftp_dir_contents = $ftp->ls;
  foreach my $entry (@ftp_dir_contents) {
    if($entry =~ /^$accession\_(.+)$/) {
      $full_assembly_path = $assembly_ftp_path."/".$&."/";
      $assembly_name = $1;
    }
  }

  unless($full_assembly_path && $assembly_name) {
    throw("Issue finding ftp path for the following GCA: $accession $assembly_ftp_path");
  }

  say "Setting assembly name for ".$accession." to ".$assembly_name;

  $assembly_hash->{'assembly_accession'} = $accession;
  $assembly_hash->{'assembly_name'} = $assembly_name;

  #get settings per clade
  my $clade_hash = clade_settings($clade);
  foreach my $key (keys(%{$clade_hash})) {
    $general_hash->{$key} = $clade_hash->{$key};
  }

  parse_assembly_report($ftp,$general_hash,$assembly_hash,$accession,$assembly_name,$full_assembly_path,$output_path);

  # Add in the species url by uppercasing the first letter of the species name
  my $species_url = $assembly_hash->{'species_name'};
  $species_url = ucfirst($species_url).'_'.$accession;
  $assembly_hash->{'species_url'} = $species_url;
  
  #Set up taxonomy adaptor to fetch rank of species using NCBI taxonomy database
  my $node_adaptor = $taxonomy_adaptor->get_TaxonomyNodeAdaptor();
  my $taxon_node = $node_adaptor->fetch_by_taxon_id($assembly_hash->{'taxon_id'});

 # Check for repeatmodeler library at species level
  if ($taxon_node->rank() eq 'subspecies' || $taxon_node->rank() eq 'no rank'){
    my $parent_node = $node_adaptor->fetch_by_taxon_id($taxon_node->parent_id());
    my $parent_name = $parent_node->names()->{'scientific name'}->[0];
    $parent_name =~ s/\s+/\_/g;
    $parent_name = lc($parent_name);
   #  Species is at subspecies level so get repeatmodeler library path if one exists at the species level
    my $repeatmodeler_file = catfile($ENV{REPEATMODELER_DIR},'species',$parent_name,$parent_name.'.repeatmodeler.fa');
    if(-e $repeatmodeler_file) {
      say "Found the following repeatmodeler file for the species:\n".$repeatmodeler_file;
      $assembly_hash->{'repeatmodeler_library'} = $repeatmodeler_file;
    } else {
      say "Did not find a repeatmodeler species library for ".$assembly_hash->{'species_name'}." on path:\n".$assembly_hash->{'species_name'};
    }
  }
  elsif ($taxon_node->rank() eq 'species'){
  #   Get repeatmodeler library path if one exists
    my $repeatmodeler_file = catfile($ENV{REPEATMODELER_DIR},'species',$assembly_hash->{'species_name'},$assembly_hash->{'species_name'}.'.repeatmodeler.fa');
    if(-e $repeatmodeler_file) {
      say "Found the following repeatmodeler file for the species:\n".$repeatmodeler_file;
      $assembly_hash->{'repeatmodeler_library'} = $repeatmodeler_file;
    } else {
      say "Did not find a repeatmodeler species library for ".$assembly_hash->{'species_name'}." on path:\n".$assembly_hash->{'species_name'};
    }
  }
  else{
    say "Rank is ",$taxon_node->rank();
  }

  create_config($assembly_hash);
  copy_general_module();
  check_compara_release_version();

  unless($config_only) {
    init_pipeline($assembly_hash,$hive_directory,$force_init,$fh);
  }
}

unless($config_only) {
  close($fh) || throw("Could not close the cmd file");
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
  my $genus_taxon_id;
  my $species_taxon_id;
  my $species_name;
  my $refseq_accession;
  my $assembly_level;
  my $wgs_id;
  my $modifier;
  open($report_file_handle, '>', \$report_file_content) || throw("could not open $report_file_name");

  unless($ftp->get($report_file_name, $report_file_handle)) {
    throw("Failed to retrieve the assembly report file: ", $ftp->message);
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
    throw("Failed to fully parse the assembly report file");
  }

  if(exists($general_hash->{'load_toplevel_only'}) && $general_hash->{'load_toplevel_only'} == 0) {
    if ($wgs_id) {
      $assembly_hash->{'wgs_id'} = $wgs_id;
    }
    else {
      throw("Need the wgs id as the load_toplevel_only flag was set to 0. Failed to find the id in the report file");
    }
  }

  my $sth = $ncbi_taxonomy->dbc->prepare("SELECT name from ncbi_taxa_name where taxon_id=? and name_class='scientific name'");
  $sth->bind_param(1,$taxon_id);
  $sth->execute;
  $species_name = $sth->fetchrow_array();
  $species_name =~ s/\s+/\_/g;
  $species_name = lc($species_name);
  unless($species_name) {
    throw("Was not able to retrieve the species name from the NCBI taxonomy db using the taxon id. Taxon id used: ".$taxon_id);
  }

  $assembly_hash->{'taxon_id'} = $taxon_id;

 #use taxon id to retrieve genus taxon id
 #genus taxon_id will be used to download genus level rnaseq data
  my $node_adaptor = $taxonomy_adaptor->get_TaxonomyNodeAdaptor();
  my $taxon_node = $node_adaptor->fetch_by_taxon_id($taxon_id);
  foreach my $ancestor ( @{ $node_adaptor->fetch_ancestors($taxon_node)}){
    #store genus level taxon id
     if ($ancestor->rank eq 'genus'){
        $genus_taxon_id = $ancestor->taxon_id;
        $assembly_hash->{'genus_taxon_id'} = $genus_taxon_id;
     }
    #store taxonomy at species level if species is a subspecies
     elsif ($ancestor->rank eq 'species'){
        $species_taxon_id = $ancestor->taxon_id;
        $assembly_hash->{'species_taxon_id'} = $species_taxon_id;
     }
     else{}
  }

  if($refseq_accession) {
    $assembly_hash->{'assembly_refseq_accession'} = $refseq_accession;
  } else {
    say "Found no RefSeq accession for this assembly";
  }

# create production name - must be unique, no more than trinomial, and include GCA
# Add dbname_accession which replace production_name in the dbname and is just the GCA MySQL friendly
  my $underscore_count = $species_name =~ tr/_//;
  $species_name =~ /(^[^_]+_[^_]+)_*.*/;
  my $binomial_name = $1;
  my $accession_append = lc($accession);
  $accession_append =~ s/\./v/g;
  $accession_append =~ s/_//g;
  $assembly_hash->{'production_name'} = $binomial_name.'_'.$accession_append;
  $assembly_hash->{'dbname_accession'} = $accession_append;

  $assembly_hash->{'strain'} = $species_name;
  $assembly_hash->{'assembly_level'} = $assembly_level;
  $assembly_hash->{'species_name'} = $species_name;
  @{$assembly_hash}{keys(%$general_hash)} = values(%$general_hash);
  $assembly_hash->{'output_path'} .= "/".$species_name."/".$accession."/";
}

=pod

=head1 Description of method

This method updates the registry database with the timestamp of when the annotation started. 
It also updates the registry with the status of the annotation as well as the user who started it.

=cut

sub update_annotation_status{
  my $accession = shift;
  my $dt   = DateTime->now;   # Stores current date and time as datetime object
  my $date = $dt->ymd;
  my $assembly_id = $assembly_registry->fetch_assembly_id_by_gca($accession);
  my ($sql,$sth);
  if ($current_genebuild == 1){
    say "Updating genebuild status to overwrite";
    $sql = "update genebuild_status set is_current = ? where assembly_accession = ?";
    $sth = $assembly_registry->dbc->prepare($sql);
    $sth->bind_param(1,0);
    $sth->bind_param(2,$accession);
    unless($sth->execute()){
      throw("Could not update annoation status for assembly with accession ".$accession);
    }
  }
  $sql = "insert into genebuild_status(assembly_accession,progress_status,date_started,genebuilder,assembly_id,is_current,annotation_source) values(?,?,?,?,?,?,?)";
  $sth = $assembly_registry->dbc->prepare($sql);
  $sth->bind_param(1,$accession);
  $sth->bind_param(2,'in progress');
  $sth->bind_param(3,$date);
  $sth->bind_param(4,$ENV{EHIVE_USER} || $ENV{USER});
  $sth->bind_param(5,$assembly_id);
  $sth->bind_param(6,1);
  $sth->bind_param(7,'pending');
  say "Accession being worked on is $accession";
  unless($sth->execute()){
   throw("Could not update annotation status for assembly with accession ".$accession);
  }
}

=pod

=head1 Description of method

This method checks if there is an existing genebuild entry for the assembly.  
If yes, genebuilder must decide whether to continue annotation or not.
If genebuilder decides to continue, then rerun with option: -current_genebuild 1
This would automatically make this new genebuild the current annotation for tracking purposes

=cut

sub check_annotation_status{
  my $accession = shift;
  my @status = $assembly_registry->fetch_genebuild_status_by_gca($accession);
  if (@status){
    if ($status[2]){
      throw("A genebuild pipeline has already been run for this assembly. "."$accession.\nGenebuild status: $status[0]\nDate started: $status[1]\nDate completed: $status[2]\nGenebuilder: $status[3],\nAnnotation source: $status[4]"."\nTo proceed with this genebuild, re-run script with option: -current_genebuild 1"  );
    }
    else{
      throw("A genebuild pipeline is currently in progress for this assembly. "."$accession\nGenebuild status: $status[0]\nDate started: $status[1]\nDate completed: Pending\nGenebuilder: $status[3]\nAnnotation source: $status[4]"."\nTo proceed with this genebuild, re-run script with option: -current_genebuild 1"  );
    }
  }
}

sub create_config {
  my ($assembly_hash) = @_;

  $assembly_hash->{'registry_path'} = catfile($assembly_hash->{'output_path'},"Databases.pm");

  foreach my $key (keys(%{$assembly_hash})) {
    say "ASSEMBLY: ". $key." => ".$assembly_hash->{$key};
  }

  my $output_path = $assembly_hash->{'output_path'};
  unless(-e $output_path) {
   system('mkdir -p '.$output_path);
  }

  my $config_string = "";
  my $past_default_options = 0;
  open(CONFIG,$ENV{ENSCODE}."/ensembl-analysis/modules/Bio/EnsEMBL/Analysis/Hive/Config/Genome_annotation_conf.pm") || throw("Could not open the config file");
  while(my $line = <CONFIG>) {
    if($line =~ /sub pipeline_create_commands/) {
      $past_default_options = 1;
    } elsif($past_default_options == 0) {
      if($line =~ /SUPER\:\:default_options/) {
        if($assembly_hash->{'dbowner'}) {
          $line .= "'dbowner' => '".$assembly_hash->{'dbowner'}."',";
        }
      }
      if($line =~ /'?([^' ]+)'?\s*=>\s*('?[^']*'?)/) {
        my $conf_key = $1;
        my $conf_val = $2;
        if(exists $assembly_hash->{$conf_key}) {
          print "REPLACING ".$conf_key." with ".$assembly_hash->{$conf_key}."\n";
          my $sub_val = "'".$assembly_hash->{$conf_key}."'";
          $line =~ s/$conf_val/$sub_val/;
        }
      }
    }

    $config_string .= $line;
  }
  close(CONFIG) || throw("CouLd not close the config file");

  $config_string =~ s/package Bio::EnsEMBL::Analysis::Hive::Config::([^;]+;)/package $1/;
  open(OUT_CONFIG,">".$output_path."/Genome_annotation_conf.pm") || throw("Could not open the config file for writing");
  print OUT_CONFIG $config_string;
  close OUT_CONFIG || throw("Could not close the config file for writing");
}

sub clade_settings {
  my ($clade) = @_;
  my $clade_settings = {

    'primates' => {
      'repbase_library'    => 'primates',
      'repbase_logic_name' => 'primates',
      'uniprot_set'        => 'primates_basic',
      'sanity_set'         => 'primates_basic',
      'projection_source_production_name' => 'homo_sapiens',
      'projection_source_db_name' => current_projection_source_db('homo_sapiens'),
    },

    'rodentia' => {
      'repbase_library'    => 'rodentia',
      'repbase_logic_name' => 'rodentia',
      'uniprot_set'        => 'mammals_basic',
      'sanity_set'         => 'rodentia_basic',      
      'projection_source_production_name' => 'mus_musculus',
      'projection_source_db_name' => current_projection_source_db('mus_musculus'),
    },

    'mammalia' => {
      'repbase_library'    => 'mammals',
      'repbase_logic_name' => 'mammals',
      'uniprot_set'        => 'mammals_basic',
      'sanity_set'         => 'mammals_basic',      
      'ig_tr_fasta_file'    => 'imgt_mammals_ig_tr.fa',
      'projection_source_production_name' => 'homo_sapiens',
      'projection_source_db_name' => current_projection_source_db('homo_sapiens'),
    },

   'marsupials' => {
      'repbase_library'    => 'mammals',
      'repbase_logic_name' => 'mammals',
      'uniprot_set'        => 'mammals_basic',
      'sanity_set'         => 'mammals_basic',      
      'ig_tr_fasta_file'    => 'multispecies_ig_tr.fa',
      'projection_source_production_name' => 'homo_sapiens',
      'projection_source_db_name' => current_projection_source_db('homo_sapiens'),
    },

    'aves' => {
      'repbase_library'    => 'Birds',
      'repbase_logic_name' => 'aves',
      'uniprot_set'        => 'birds_basic',
      'sanity_set'         => 'birds_basic',      
      'projection_source_production_name' => 'homo_sapiens',
      'projection_source_db_name' => current_projection_source_db('homo_sapiens'),
    },

    'reptiles' => {
      'repbase_library'    => 'vertebrates',
      'repbase_logic_name' => 'vertebrates',
      'uniprot_set'        => 'reptiles_basic',
      'sanity_set'         => 'reptiles_basic',      
      'masking_timer_long'  => '6h',
      'masking_timer_short' => '3h',
      'projection_source_production_name' => 'homo_sapiens',
      'projection_source_db_name' => current_projection_source_db('homo_sapiens'),
    },

    'teleostei' => {
      'repbase_library'     => 'Teleostei',
      'repbase_logic_name'  => 'teleost',
      'uniprot_set'         => 'fish_basic',
      'sanity_set'          => 'fish_basic',      
      'ig_tr_fasta_file'    => 'fish_ig_tr.fa',
      'masking_timer_long'  => '6h',
      'masking_timer_short' => '3h',
      'skip_projection'    => 1,
      'skip_lastz'         => 1,
      # need a default projection source db set
      'projection_source_production_name' => 'danio_rerio',
      'projection_source_db_name' => current_projection_source_db('danio_rerio'),
      # need a different value for creating repeatmasker slices
      repeatmasker_slice_size => 500000,
    },

    'distant_vertebrate' => {
      'repbase_library'    => 'vertebrates',
      'repbase_logic_name' => 'vertebrates',
      'uniprot_set'        => 'distant_vertebrate',
      'sanity_set'         => 'distant_vertebrate',      
      'masking_timer_long'  => '6h',
      'masking_timer_short' => '3h',
      'projection_source_production_name' => 'homo_sapiens',
      'projection_source_db_name' => current_projection_source_db('homo_sapiens'),
      # need a different value for creating repeatmasker slices and batching slices to avoid long-running jobs
      repeatmasker_slice_size => 300000,
      batch_target_size => 300000,
    },


    'insects' => {
      'repbase_library'    => 'insecta',
      'repbase_logic_name' => 'insects',
      'uniprot_set'        => 'insects_basic',
      'sanity_set'         => 'insects_basic',      
      'projection_source_production_name' => 'homo_sapiens',
      'projection_source_db_name' => current_projection_source_db('homo_sapiens'),
    },

    'lepidoptera' => {
      'repbase_library'    => 'lepidoptera',
      'repbase_logic_name' => 'lepidoptera',
      'uniprot_set'        => 'lepidoptera_basic',
      'sanity_set'         => 'lepidoptera_basic',      
      'protein_blast_db'   => '/hps/nobackup/flicek/ensembl/genebuild/blastdb/proteomes/HMLEP',
      'protein_blast_index'=> '/hps/nobackup/flicek/ensembl/genebuild/blastdb/proteomes/HMLEP_index',
      'skip_projection'    => 1,
      'skip_lastz'         => 1,
      # need a default projection source db set - for now use human and projection is skipped, will consider updating to use a butterfly annotation
      'projection_source_production_name' => 'homo_sapiens',
      'projection_source_db_name' => current_projection_source_db('homo_sapiens'),
    },

    'hymenoptera' => {
      'repbase_library'    => 'hymenoptera',
      'repbase_logic_name' => 'hymenoptera',
      'uniprot_set'        => 'hymenoptera_basic',
      'sanity_set'         => 'hymenoptera_basic',      
      'skip_projection'    => 1,
      'skip_lastz'         => 1,
      # need a default projection source db set - for now use human and projection is skipped, will consider updating to use a hymenoptera annotation
      'projection_source_production_name' => 'homo_sapiens',
      'projection_source_db_name' => current_projection_source_db('homo_sapiens'),
    },

    'atroparvus' => {
      'repbase_library'    => 'insecta',
      'repbase_logic_name' => 'insects',
      'uniprot_set'        => 'insects_basic',
      'sanity_set'         => 'insects_basic',      
      'protein_blast_db'   => '/hps/nobackup/flicek/ensembl/genebuild/blastdb/proteomes/5_01_21-A_atroparvus/Combined_A.atroparvus_n356016',
      'protein_blast_index'=> '/hps/nobackup/flicek/ensembl/genebuild/blastdb/proteomes/5_01_21-A_atroparvus/Combined_A.atroparvus_n356016_index',
      'skip_projection'    => 1,
      'skip_lastz'         => 1,
      'projection_source_production_name' => 'homo_sapiens',
      'projection_source_db_name' => current_projection_source_db('homo_sapiens'),
    },

   'perniciosus' => {
      'repbase_library'    => 'insecta',
      'repbase_logic_name' => 'insects',
      'uniprot_set'        => 'insects_basic',
      'sanity_set'         => 'insects_basic',      
      'protein_blast_db'   => '/hps/nobackup/flicek/ensembl/genebuild/blastdb/proteomes/P_perniciosus/Phlebotomus_Uniprot_Ensembl_DB',
      'protein_blast_index'=> '/hps/nobackup/flicek/ensembl/genebuild/blastdb/proteomes/P_perniciosus/Phlebotomus_Uniprot_Ensembl_DB_index',
      'skip_projection'    => 1,
      'skip_lastz'         => 1,
      'projection_source_production_name' => 'homo_sapiens',
      'projection_source_db_name' => current_projection_source_db('homo_sapiens'),
    },

    'non_vertebrates' => {
      'repbase_library'    => 'non_vertebrates',
      'repbase_logic_name' => 'non_vertebrates',
      'uniprot_set'        => 'non_vertebrates_basic',
      'sanity_set'         => 'non_vertebrates_basic',      
      'projection_source_production_name' => 'homo_sapiens',
      'projection_source_db_name' => current_projection_source_db('homo_sapiens'),
    },

    'fungi' => {
      'repbase_library'    => 'fungi',
      'repbase_logic_name' => 'fungi',
      'uniprot_set'        => 'fungi_basic',
      'sanity_set'         => 'fungi_basic',      
      'projection_source_production_name' => 'homo_sapiens',
      'projection_source_db_name' => current_projection_source_db('homo_sapiens'),
    },

    'metazoa' => {
      'repbase_library'    => 'metazoa',
      'repbase_logic_name' => 'metazoa',
      'uniprot_set'        => 'metazoa_basic',
      'sanity_set'         => 'metazoa_basic',      
      'projection_source_production_name' => 'homo_sapiens',
      'projection_source_db_name' => current_projection_source_db('homo_sapiens'),
    },

    'plants'  => {
      'repbase_library'    => 'plants',
      'repbase_logic_name' => 'plants',
      'uniprot_set'        => 'plants_basic',
      'sanity_set'         => 'plants_basic',      
      'projection_source_production_name' => 'homo_sapiens',
      'projection_source_db_name' => current_projection_source_db('homo_sapiens'),
    },

    'protists'  => {
      'repbase_library'    => 'protists',
      'repbase_logic_name' => 'protists',
      'uniprot_set'        => 'protists_basic',
      'sanity_set'         => 'protists_basic',      
      'projection_source_production_name' => 'homo_sapiens',
      'projection_source_db_name' => current_projection_source_db('homo_sapiens'),
    },

  };

  unless($clade_settings->{$clade}) {
    throw("Could not find the clade specified in the clade_settings hash. Clade specified: ".$clade);
  }

  return($clade_settings->{$clade});
}


sub custom_load_data {
  my ($general_hash,$ncbi_taxonomy) = @_;

  my $required_keys = {"custom_toplevel_file_path"    => 1,
                       "clade"                        => 1,
                       "assembly_name"                => 1,
                       "taxon_id"                     => 1,
                       "custom_protein_blastdb"       => 1,
                       "custom_protein_blastdb_index" => 1};

  my $expected_keys = {"load_toplevel_only"        => 2,
                       "use_repeatmodeler_to_mask" => 0,
                       "rnaseq_summary_file"       => '',
                       "skip_post_repeat_analyses" => 1,
                       "skip_projection"           => 1,
                       "stable_id_prefix"          => "CUSTOM",
                       "species_division"          => "CUSTOM",
                       "release_number"            => 777,
                       "genebuilder_id"            => 0};

  my $critical_failure = 0;
  foreach my $required_key (keys(%$required_keys)) {
    unless($general_hash->{$required_key}) {
      warning("Could not find a value for the following required key in the ini file: ".$required_key);
      $critical_failure++;
    }
  }

  if($critical_failure) {
    throw("A critical key was missing a value");
  }

  # Get clade
  my $clade;
  if ($general_hash->{'clade'}) {
    $clade = $general_hash->{'clade'};
  }
  else {
    $clade = $assembly_registry->fetch_clade_by_taxon_id($general_hash->{'taxon_id'});
  }

  my $clade_hash = clade_settings($clade);
  foreach my $key (keys(%{$clade_hash})) {
    $general_hash->{$key} = $clade_hash->{$key};
  }

  # Get the species name
  my $taxon_id = $general_hash->{'taxon_id'};
  my $species_name;
  my $sth = $ncbi_taxonomy->dbc->prepare("SELECT name from ncbi_taxa_name where taxon_id=? and name_class='scientific name'");
  $sth->bind_param(1,$taxon_id);
  $sth->execute;
  $species_name = $sth->fetchrow_array();
  $species_name =~ s/\s+/\_/g;
  $species_name = lc($species_name);
  unless($species_name) {
    throw("Was not able to retrieve the species name from the NCBI taxonomy db using the taxon id. Taxon id used: ".$taxon_id);
  }
  $general_hash->{'species_name'} = $species_name;


  foreach my $expected_key (keys(%$expected_keys)) {
    unless(exists $general_hash->{$expected_key}) {
      warning("Could not find a value for the following expected key in the ini file: ".$expected_key."\nUsing default value: ".$expected_keys->{$expected_key});
      $general_hash->{$expected_key} = $expected_keys->{$expected_key};
    }
  }

  #Set up taxonomy adaptor to fetch rank of species using NCBI taxonomy database
  my $node_adaptor = $taxonomy_adaptor->get_TaxonomyNodeAdaptor();
  my $taxon_node = $node_adaptor->fetch_by_taxon_id($general_hash->{'taxon_id'});

  #Check for repeatmodeler library at species level
  if ($taxon_node->rank() eq 'subspecies' || $taxon_node->rank() eq 'no rank'){
    my $parent_node = $node_adaptor->fetch_by_taxon_id($taxon_node->parent_id());
    my $parent_name = $parent_node->names()->{'scientific name'}->[0];
    $parent_name =~ s/\s+/\_/g;
    $parent_name = lc($parent_name);
    # Species is at subspecies level so get repeatmodeler library path if one exists at the species level. Note that this will be overwritten with a custom path if one has been provided
    my $repeatmodeler_file = catfile($ENV{REPEATMODELER_DIR},'species',$parent_name,$parent_name.'.repeatmodeler.fa');
    if(-e $repeatmodeler_file) {
      say "Found the following repeatmodeler file for the species:\n".$repeatmodeler_file;
      $general_hash->{'repeatmodeler_library'} = $repeatmodeler_file;
    } else {
      say "Did not find a repeatmodeler species library for ".$general_hash->{'species_name'}." on path:\n".$general_hash->{'species_name'};
    }
  }
  elsif ($taxon_node->rank() eq 'species'){
   #  Get repeatmodeler library path if one exists. Note that this will be overwritten with a custom path if one has been provided
    my $repeatmodeler_file = catfile($ENV{REPEATMODELER_DIR},'species',$general_hash->{'species_name'},$general_hash->{'species_name'}.'.repeatmodeler.fa');
    if(-e $repeatmodeler_file) {
      say "Found the following repeatmodeler file for the species:\n".$repeatmodeler_file;
      $general_hash->{'repeatmodeler_library'} = $repeatmodeler_file;
    } else {
      say "Did not find a repeatmodeler species library for ".$general_hash->{'species_name'}." on path:\n".$general_hash->{'species_name'};
    }
  }
  else{
    say "Rank is ",$taxon_node->rank();
  }

}

sub init_pipeline {
    my ($assembly_hash,$hive_directory,$force_init,$fh) = @_;

    chdir($assembly_hash->{'output_path'});
    my $cmd;
    if ($hive_directory) {
      $cmd = 'perl '.catfile($hive_directory, 'scripts', 'init_pipeline.pl');
    }
    else {
      print "WARNING using init_pipeline.pl from your PATH, may not be the same as your PERL5LIB\n";
      $cmd = 'init_pipeline.pl';
    }
    if ($force_init) {
      $cmd .= ' Genome_annotation_conf.pm -hive_force_init 1';
    }
    else {
      $cmd .= ' Genome_annotation_conf.pm';
    }
    my $result = `$cmd`;
    unless($result =~ /beekeeper.+\-sync/) {
      throw("Failed to run init_pipeline for ".$assembly_hash->{'species_name'}."\nCommandline used:\n".$cmd);
    }
    unless ($custom_load){
	update_annotation_status($assembly_hash->{'assembly_accession'});
    }
    my $sync_command = $&;
    if ($hive_directory) {
      $sync_command = 'perl '.catdir($hive_directory, 'scripts').catfile('','').$sync_command; # The crazy catfile in the middle is to get the path separator
    }
    my $return = system($sync_command);
    if($return) {
      throw("Failed to sync the pipeline for ".$assembly_hash->{'species_name'}."\nCommandline used:\n".$sync_command);
    }


    if($check_for_transcriptomic) {
      my $run_command = $sync_command;
      $run_command =~ s/\-sync/\-run/;
      my $return = system($run_command);
      if($return) {
        throw("Failed to run a loop the pipeline for ".$assembly_hash->{'species_name'}."\nCommandline used:\n".$run_command);
      }
    }


    my $loop_command = $sync_command;
    $loop_command =~ s/\-sync/\-loop \-sleep 0.3/;
    my ($ehive_url) = $sync_command =~ /url\s+(\S+)/;
    my ($driver, $user, $password, $host, $port, $dbname) = $ehive_url =~ /^"?(\w+)[:\/]+(\w*):(\w+).([^:]+):(\d+)\/(\w+)"?$/;
    if ($password) {
      $password = '&passwd=xxxxx';
    }

    say $fh "#GuiHive: $base_guihive/?driver=$driver&username=$user&host=$host&port=$port&dbname=$dbname$password";
    say $fh "export EHIVE_URL=$ehive_url";
    say $fh $loop_command;
}


=head2 assign_server_info

 Arg [1]    : Hashref, the hash with all the information about the config
 Description: Assign servers connection details based on 'server_set' if it is present.
              If the set does not exists, it uses the values for "set1"
              If 'server_set' is not defined, it looks for the following keys which should
              all be defined:
                databases_host
                databases_port
                pipe_db_host
                pipe_db_port
                dna_db_host
                dna_db_port
 Returntype : None
 Exceptions : Throws if 'server_set' is not set and none of the connection details are set.

=cut

sub assign_server_info {
  my ($general_hash) = @_;

  my $servers = {
    set1 => {
              pipe_db_host => "mysql-ens-genebuild-prod-4",
              pipe_db_port   => 4530,
              databases_host  => "mysql-ens-genebuild-prod-3",
              databases_port    => 4529,
              dna_db_host  => "mysql-ens-genebuild-prod-2",
              dna_db_port    => 4528,
            },

    set2 => {
              pipe_db_host => "mysql-ens-genebuild-prod-7",
              pipe_db_port   => 4533,
              databases_host  => "mysql-ens-genebuild-prod-5",
              databases_port    => 4531,
              dna_db_host  => "mysql-ens-genebuild-prod-6",
              dna_db_port    => 4532,
            },
  };

  my $server_set = $general_hash->{'server_set'};
  if ($server_set) {
    if (! exists $servers->{$server_set}) {
      warning("Could not find an associated server set entry in the HiveBaseConfig for ".$server_set.". Will default to set1");
      $server_set = 'set1';
    }
    $general_hash->{databases_host} = $servers->{$server_set}->{'databases_host'};
    $general_hash->{databases_port} = $servers->{$server_set}->{'databases_port'};
    $general_hash->{pipe_db_host} = $servers->{$server_set}->{'pipe_db_host'};
    $general_hash->{pipe_db_port} = $servers->{$server_set}->{'pipe_db_port'};
    $general_hash->{dna_db_host} = $servers->{$server_set}->{'dna_db_host'};
    $general_hash->{dna_db_port} = $servers->{$server_set}->{'dna_db_port'};

  }
  else {
    throw("You are missing connection details for at least one of them: databases_host, databases_port, pipe_db_host, pipe_db_port, dna_db_host, dna_db_port")
      unless (
        exists $general_hash->{databases_host} and
        exists $general_hash->{databases_port} and
        exists $general_hash->{pipe_db_host} and
        exists $general_hash->{pipe_db_port} and
        exists $general_hash->{dna_db_host} and
        exists $general_hash->{dna_db_port}
      );
  }
}

=head2 copy_general_module

 Arg [1]    : None
 Description: Find the path to 'ensembl-analysis' in @INC, then copy Bio/EnsEMBL/Analysis/Config/General.pm.example
              to Bio/EnsEMBL/Analysis/Config/General.pm to be able to run the projection pipeline and Compara parts
 Returntype : None
 Exceptions : Throws if it cannot find 'ensembl-analysis'

=cut

sub copy_general_module {
  my ($analysis_path) = grep {/ensembl-analysis/} @INC;
  throw("I cannot find ensembl-analysis in your PERL5LIB") unless($analysis_path);
  my @dirs = splitdir($analysis_path);
  while (my $dir = pop(@dirs)) {
    if ($dir eq 'ensembl-analysis') {
      my $file = catfile(@dirs, $dir, 'modules', 'Bio', 'EnsEMBL', 'Analysis', 'Config', 'General.pm');
      if (! -e $file) {
        copy("$file.example", $file) or warning("Could not copy $file.example to $file");
      }
      return;
    }
  }
}


=head2 check_compara_release_version

 Arg [1]    : None
 Description: Checks if the branch in ensembl-compara is set to the correct version
              in order to run the projection part successfully.
              If the branch is not the wanted version, it tries to switch to the wanted
              version
 Returntype : None
 Exceptions : Throws if it fails to move to the ensembl-compara directory
              Throws if it fails to move back to the current directory
              Throws if it fails to run 'git branch'
              Throws if it fails to close the git command pipe
              Throws if it fails to run 'git fetch'
              Throws if it fails to run 'git checkout <wanted version>'

=cut

sub check_compara_release_version {
  my $wanted_branch = 'release/98';
  foreach my $dir (@INC) {
    if ($dir =~ /ensembl-compara/) {
      my $current_dir = getcwd();
      chdir($dir) or throw("Could not go to '$dir'");
      open(CH, 'git branch |') or throw('Could not run "git branch"');
      while(<CH>) {
        my ($current, $branch) = $_ =~ /^(.)\s+(\S+)/;
        if ($current eq '*') {
          if ($branch ne $wanted_branch) {
            if (system('git fetch > '.devnull().' 2>&1')) {
              throw('Could not fetch branches with git');
            }
            elsif (system("git checkout $wanted_branch > ".devnull().' 2>&1')) {
              throw("Could not checkout '$wanted_branch'");
            }
            else {
              warning("Switched branch '$branch' to '$wanted_branch' in '$dir'");
            }
          }
          last;
        }
      }
      close(CH) or throw('Could not close the git command pipe');
      chdir($current_dir) or throw("Could not move back to '$current_dir'");
      last;
    }
  }
}

=head2 

 Arg [1]    : Species name, e.g. mus_musculus
 Description: Looks for the most recent core on the mirror server for the given 
              species - to be used as reference for projection.
 Returntype : String
 Exceptions : Warning if no core exists on mirror for given species

=cut
sub current_projection_source_db{
  my ($species) = @_;
  my $current_db;

  my @out= `mysql-ens-mirror-1 -NB -e "SHOW DATABASES LIKE '${species}_core%'"`;
  if(!@out){
    warning("No core database available on mirror for species +$species+ - will use latest homo_sapiens core.");
    @out= `mysql-ens-mirror-1 -NB -e "SHOW DATABASES LIKE 'homo_sapiens_core%'"`;
  }
  chomp($out[-1]); #the 3 most recent versions of a core will be available on mirror
  return $out[-1];
}
