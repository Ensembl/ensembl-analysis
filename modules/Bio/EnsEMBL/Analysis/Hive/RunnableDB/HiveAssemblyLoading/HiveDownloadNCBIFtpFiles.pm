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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveDownloadNCBIFtpFiles;

use strict;
use warnings;
use feature 'say';

use File::Spec::Functions qw(splitdir catdir);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;
  return 1;
  unless($self->param('full_ftp_path') && $self->param('output_path')) {
    $self->throw("Must pass in the following parameters:\nfull_ftp_path e.g /genomes/genbank/vertebrate_mammalian/Canis_lupus/".
                 "all_assembly_versions/GCA_000002285.2_CanFam3.1/GCA_000002285.2_CanFam3.1_assembly_structure/\n".
                 "output_path e.g /path/to/work/dir\n".
                 "primary_assembly_dir_name e.g. Primary_Assembly (usually this on NCBI ftp)");
  }

}

sub run {
  my $self = shift;

  my $ftp_species_path = $self->param('full_ftp_path')."/";
  say "Using the following as the ftp path:\n".$ftp_species_path;
  my $local_path = $self->param('output_path');
  $self->download_ftp_dir($ftp_species_path,$local_path);
  $self->download_non_nuclear_ftp_dir($self->param('full_ftp_path'),$self->param('output_path'));
  $self->unzip($local_path);

  say "Finished downloading NCBI ftp structure and files";
  return 1;
}

sub download_ftp_dir {
  my ($self,$ftp_path,$local_dir) = @_;

  my $wget_verbose = "-nv";

  # Test
  #$ftp_path = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/151/905/GCA_000151905.3_gorGor4/';
  #$ftp_path = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/006/655/GCA_900006655.1_GSMRT3/';

  my @ftp_path_array = split('/',$ftp_path);
  my $top_dir = pop(@ftp_path_array);

  my $cmd = "wget --no-cache --spider ".$ftp_path."/".$top_dir."_assembly_structure/ 2>&1 |";
  say "Checking for single or multilevel assembly with:\n".$cmd;

  my $single_level = 0;

  # If there is a failure to find the assembly_structure dir then it is a single level assembly (assuming the path is correct)
  open CMDOUT, $cmd;
  while (my $output = <CMDOUT>) {
    if($output =~ /No such directory/) {
      $single_level = 1;
    }
  }
  close CMDOUT;

  system('mkdir '.$local_dir."/Primary_Assembly");
  if($single_level) {
    $cmd = "wget ".$wget_verbose." -nH -P ".$local_dir." ".$ftp_path."/".$top_dir."_assembly_report.txt -O ".$local_dir."/Primary_Assembly/assembly_report.txt";
    if (system($cmd)) {
      $self->throw("Could not download *_assembly_report.txt file to ".$local_dir.". Please, check that both paths are valid.\n".
                   "Commandline used:\n".$cmd);
    }
  } else {
    # This is some magic to get the correct number for --cut-dir
    my @dirs = splitdir($ftp_path."/".$top_dir."_assembly_structure/Primary_Assembly");
    my $numDirs = scalar(@dirs)-2;
    $cmd = "wget --no-proxy ".$wget_verbose." -r -nH --cut-dirs=".$numDirs." --reject *.rm.out.gz -P ".$local_dir."/Primary_Assembly/ ".$ftp_path."/".$top_dir."_assembly_structure/Primary_Assembly";
    my $return = system($cmd);
    if($return) {
      $self->throw("Could not download the AGP, FASTA and info files. Commandline used:\n".$cmd);
    }

    # check if no file was downloaded
    my $num_files = int(`find $local_dir -type f | wc -l`);
    if ($num_files == 0) {
      $self->throw("No file was downloaded from ".$ftp_path." to ".$local_dir.". Please, check that both paths are valid.");
    } else {
      say "$num_files files were downloaded";
    }

    # get the name of the directory where the assembly report file is, replace any // with / (except for ftp://) in the path to stop issues with the pop
    $ftp_path =~ s/([^\:])\/\//$1\//;

    my @ftp_path_array = split('/',$ftp_path);
    $cmd = "wget ".$wget_verbose." -nH -P ".$local_dir." ".$ftp_path."/".$top_dir."_assembly_report.txt -O ".$local_dir."/Primary_Assembly/assembly_report.txt";
    if (system($cmd)) {
      $self->throw("Could not download *_assembly_report.txt file to ".$local_dir.". Please, check that both paths are valid.\n".
                   "Commandline used:\n".$cmd);
    }
    else {
      say "Assembly report file was downloaded";
    }
  } # end else for single_level
}

sub download_non_nuclear_ftp_dir {
  my ($self,$ftp_path,$local_dir) = @_;
  my $non_nuclear = "/non-nuclear/assembled_chromosomes/AGP/";
  my $full_ftp_path = $ftp_path.$non_nuclear;
  my $full_local_dir_path = $local_dir."/non_nuclear_agp/";

  my $cmd = "wget -P ".$full_local_dir_path." ".$full_ftp_path."*.agp.gz";
  if (system($cmd)) {
    $self->warning("No non-nuclear AGP files found. Continuing anyway\n".
                   "Commandline used:\n".$cmd);
    return;
  } else {
    say "Downloaded non-nuclear AGP files (can be used to identify mito later in contigs.fa)";
  }

  $self->unzip($full_local_dir_path);
  $self->concat_non_nuclear($full_local_dir_path);
}


sub unzip {
  my ($self,$path) = @_;
  say "Unzipping the compressed files...";
  $self->throw("gunzip operation failed. Please check your error log file.") if (system("gunzip -r $path") == 1);
  say "Unzipping finished!";
}

sub concat_non_nuclear {
  my ($self,$path) = @_;
  my $cmd = "cat ".$path."/*.agp > ".$path."/concat_non_nuclear.agp";
  if (system($cmd)) {
    $self->throw("Problem concatenation the non nuclear AGP files into concat_non_nuclear.agp\n".
                 "Commandline used:\n".$cmd);
  }
}


sub write_output {
  my $self = shift;

  return 1;
}

1;
