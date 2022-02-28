#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveSetupAssemblyPipeline;

use strict;
use warnings;
use feature 'say';


use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;

  my $output_path = $self->param('output_path');
  my $ftp_link_file  = $self->param('ftp_link_file');
  my $genebuilder_id = $self->param('genebuilder_id');
  my $farm_user_name = $self->param('farm_user_name');
  my $db_write_user = $self->param('user_w');
  my $db_write_pass = $self->param('pass_w');
  my $reference_db_server = $self->param('reference_db_server');
  my $reference_db_port = $self->param('reference_db_port');
  my $create_type = $self->param('create_type');
  my $pipe_db_server = $self->param('pipe_db_server');
  my $production_db = $self->param('production_db');
  my $taxonomy_db = $self->param('taxonomy_db');
  my $enscode_root_dir = $self->param('enscode_root_dir');
  my $primary_assembly_dir_name = $self->param('primary_assembly_dir_name');
  my $contigs_source = $self->param('contigs_source');
  # Add production db and tax db here once modules done

  unless(-e $ftp_link_file) {
   $self->throw("The ftp link file does not exist on the path provided. Path used:\n".$ftp_link_file);
  }

  return 1;

}

sub run {
  my $self = shift;

  my $ftp_link_file  = $self->param('ftp_link_file');
  $self->parse_ftp_link_file($ftp_link_file);

  return 1;
}

sub write_output {
  my $self = shift;

  my $output_array = $self->param('output_array');
  unless(scalar(@{$output_array})) {
    $self->throw('The output array was empty, meaning no species were processed into output ids');
  }

  foreach my $output_hash (@{$output_array}) {
    $self->dataflow_output_id($output_hash,1);
    $self->dataflow_output_id($output_hash,2);
    $self->dataflow_output_id($output_hash,3);
  }

  return 1;
}


sub parse_ftp_link_file {
  my ($self) = @_;

  my $output_path = $self->param('output_path');
  my $ftp_link_file  = $self->param('ftp_link_file');
  my $genebuilder_id = $self->param('genebuilder_id');
  my $farm_user_name = $self->param('farm_user_name');
  my $db_write_user = $self->param('user_w');
  my $db_write_pass = $self->param('pass_w');
  my $reference_db_server = $self->param('reference_db_server');
  my $reference_db_port = $self->param('reference_db_port');
  my $create_type = $self->param('create_type');
  my $production_db = $self->param('production_db');
  my $taxonomy_db = $self->param('taxonomy_db');
  my $pipe_db_server = $self->param('pipe_db_server');
  my $enscode_root_dir = $self->param('enscode_root_dir');
  my $primary_assembly_dir_name = $self->param('primary_assembly_dir_name');
  my $contigs_source = $self->param('contigs_source');
  my $chromosomes_present = 0;
  my $taxon_id;
  my $species_name;
  my $wgs_id;
  my $assembly_level;
  my $output_array = [];

  open(REPORT_FILE,">".$output_path."/setup_report_file.txt");
  say REPORT_FILE "#############################################################";
  say REPORT_FILE "# GENERAL SETTINGS                                          #";
  say REPORT_FILE "#############################################################";
  say REPORT_FILE "Master output path: ".$output_path;
  say REPORT_FILE "FTP link file: ".$ftp_link_file;
  say REPORT_FILE "Genebuilder id: ".$genebuilder_id;
  say REPORT_FILE "Farm user name: ".$farm_user_name;
  say REPORT_FILE "DB write user: ".$db_write_user;
  say REPORT_FILE "DB write pass: ".$db_write_pass;
  say REPORT_FILE "Reference DB server: ".$reference_db_server;
  say REPORT_FILE "Reference DB port: ".$reference_db_port;
  say REPORT_FILE "Creation type for reference DB: ".$create_type;
  say REPORT_FILE "Pipeline DB server: ".$pipe_db_server;
  say REPORT_FILE "ENSCODE root dir: ".$enscode_root_dir;
  say REPORT_FILE "Primary assembly dir name: ".$primary_assembly_dir_name;
  say REPORT_FILE "Contigs source: ".$contigs_source."\n";

  open(FTP_LINK_FILE,$ftp_link_file);
  while(<FTP_LINK_FILE>) {
    my $link = $_;
    chomp $link;

    # Chuck all the generic stuff into the output hash (suspect there is an easier way of doing this)
    my $output_hash = {
                        'output_path' => $output_path,
                        'ftp_link_file' => $ftp_link_file,
                        'genebuilder_id' => $genebuilder_id,
                        'farm_user_name' => $farm_user_name,
                        'user_w' => $db_write_user,
                        'pass_w' => $db_write_pass,
                        'reference_db_server' => $reference_db_server,
                        'reference_db_port' => $reference_db_port,
                        'create_type' => $create_type,
                        'production_db' => $production_db,
                        'taxonomy_db' => $taxonomy_db,
                        'pipe_db_server' => $pipe_db_server,
                        'enscode_root_dir' => $enscode_root_dir,
                        'primary_assembly_dir_name' => $primary_assembly_dir_name,
                        'contigs_source' => $contigs_source,
                      };

    unless($link =~ /^ftp/) {
      say "Skipping the line that did not start with ftp";
      next;
    }

    unless($link =~ /\/genomes\/genbank\/[^\/]+\/([^\/]+)\/all_assembly_versions\/([^\/]+)\//) {
      $self->throw("Failed to parse the following line:\n".$link."\n\nExpected a link to the main dir for the species. For example:\n".
                   "ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Propithecus_coquereli/all_assembly_versions/GCA_000956105.1_Pcoq_1.0/");
    }

    $species_name = $1;
    $species_name = lc($species_name);

    my $gca_and_assembly_version = $2;
    unless($gca_and_assembly_version =~ /^(GCA\_\d+\.\d+)\_(.+)/) {
      $self->throw("Couldn't parse the assembly version out of the ftp link. Tried to parse off end of:\n".$gca_and_assembly_version);
    }

    my $gca = $1;

    my $assembly_version = $2;
    $output_hash->{'assembly_version'} = $assembly_version;

    my $report_species_name = uc($species_name);
    say REPORT_FILE "-------------------------------------------------------------";
    say REPORT_FILE "- ".$report_species_name." SETTINGS";
    say REPORT_FILE "-------------------------------------------------------------";
    say REPORT_FILE "Species dir path: ".$output_path."/".$species_name;
    say REPORT_FILE "Assembly name (coord_system_version): ".$assembly_version;
    say REPORT_FILE "GCA accession: ".$gca;

    my $cmd = "wget ".$link."/".$gca_and_assembly_version."_assembly_report.txt";
    my $return = system($cmd);

    if($return) {
      $self->throw("wget for assembly report failed. Commandline used:\n".$cmd);
    }

    open(ASSEMBLY_REPORT,$gca_and_assembly_version."_assembly_report.txt");
    my @assembly_report = <ASSEMBLY_REPORT>;
    close ASSEMBLY_REPORT;

    my %wgs_code = ();
    foreach my $assembly_line (@assembly_report) {
      chomp $assembly_line;

      if($assembly_line =~ /^#/) {
        if($assembly_line =~ /Taxid\:[^\d]+(\d+)/) {
          $taxon_id = $1;
        } elsif($assembly_line =~ /Assembly level: ([a-zA-Z]+)/) {
          $assembly_level = $1;
        }

        next;
      }

      my @assembly_columns = split("\t",$assembly_line);
      my $accession = $assembly_columns[4];
      unless($accession && ($accession =~ /(^[A-Z]{4})/ || $accession =~ /^gb\|([A-Z]{4})/)) {
        next;
      }

      my $code = $1;
      if(exists($wgs_code{$code})) {
        $wgs_code{$code}++;
      } else {
        $wgs_code{$code} = 1;
      }
    }

    unless($taxon_id) {
      $self->throw("Failed to find and parse 'Taxid' line from report file. File used:\n".$gca_and_assembly_version."_assembly_report.txt");
    }
    say REPORT_FILE "Taxon id: ".$taxon_id;

    unless($assembly_level) {
      $self->throw("Failed to find and parse 'Assembly level' line from report file. File used:\n".$gca_and_assembly_version."_assembly_report.txt");
    }

    $assembly_level = lc($assembly_level);
    unless($assembly_level eq 'scaffold' || $assembly_level eq 'chromosome') {
      $self->throw("Parsed assembly level from report file but it was not 'scaffold' or 'chromosome'. Level found:\n".$assembly_level);
    }
    say REPORT_FILE "Assembly level: ".$assembly_level;

    if($assembly_level eq 'chromosome') {
      $chromosomes_present = 1;
    }

    # This section deals with either getting 0, 1 or many potential wgs codes. Getting 1 is ideal
    my @code_types = keys(%wgs_code);
    if(scalar(@code_types == 0)) {
      $self->throw("Failed to parse any potential 4 letter wgs code from the accessions in ".$gca_and_assembly_version."_assembly_report.txt");
    }

    if(scalar(@code_types > 1)) {
      $self->throw("Failed because of multiple potential 4 letter wgs code from the accessions in ".$gca_and_assembly_version."_assembly_report.txt, ".
            " codes parsed:\n".@code_types);
    }

    $wgs_id = $code_types[0];
    say REPORT_FILE "wgs id: ".$wgs_id;

    my $full_ftp_path = $link."/".$gca_and_assembly_version."_assembly_structure";

    my $reference_db_name = $farm_user_name."_".$species_name."_ref";
    say REPORT_FILE "Reference DB name: ".$reference_db_name."\n";
    my $target_db = { -dbname => $reference_db_name,
                      -user => $db_write_user,
                      -pass => $db_write_pass,
                      -host => $reference_db_server,
                      -port => $reference_db_port,
                    };

    $output_hash->{'gca'} = $gca;
    $output_hash->{'gca_and_assembly_version'} = $gca_and_assembly_version;
    $output_hash->{'taxon_id'} = $taxon_id;
    $output_hash->{'assembly_level'} = $assembly_level;
    $output_hash->{'chromosomes_present'} = $chromosomes_present;
    $output_hash->{'species_name'} = $species_name;
    $output_hash->{'full_ftp_path'} = $full_ftp_path;
    $output_hash->{'wgs_id'} = $wgs_id;
    $output_hash->{'target_db'} = $target_db;
    $output_hash->{'coord_system_version'} = $assembly_version;

    push(@{$output_array},$output_hash);
  }
  close FTP_LINK_FILE;
  close REPORT_FILE;

  $self->param('output_array',$output_array);
}

1;
