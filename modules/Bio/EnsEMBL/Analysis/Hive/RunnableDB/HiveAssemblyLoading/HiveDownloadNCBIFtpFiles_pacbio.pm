#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2017] EMBL-European Bioinformatics Institute
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
  unless($self->param('full_ftp_path') && $self->param('output_path') && $self->param('primary_assembly_dir_name')) {
    $self->throw("Must pass in the following parameters:\nfull_ftp_path e.g /genomes/genbank/vertebrate_mammalian/Canis_lupus/".
                 "all_assembly_versions/GCA_000002285.2_CanFam3.1/GCA_000002285.2_CanFam3.1_assembly_structure/\n".
                 "output_path e.g /path/to/work/dir\n".
                 "primary_assembly_dir_name e.g. Primary_Assembly (usually this on NCBI ftp)");
  }

}

sub run {
  my $self = shift;

  my $ftp_species_path = $self->param('full_ftp_path')."/".$self->param('primary_assembly_dir_name');
  my $local_path = catdir($self->param('output_path'), $self->param('primary_assembly_dir_name'));
  $self->download_ftp_dir($ftp_species_path,$local_path);
  $self->download_non_nuclear_ftp_dir($self->param('full_ftp_path'),$self->param('output_path'));
  $self->unzip($local_path);

  say "Finished downloading NCBI ftp structure and files";
  return 1;
}

sub download_ftp_dir {
  my ($self,$ftp_path,$local_dir) = @_;

  my $wget_verbose = "-nv";

  # This is some magic to get the correct number for --cut-dir
  my @dirs = splitdir($ftp_path);
  my $numDirs = scalar(@dirs)-2;

  my $cmd = "wget --no-proxy ".$wget_verbose." -r -nH --cut-dirs=".$numDirs." --reject *.rm.out.gz -P ".$local_dir." ".$ftp_path;
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
 # my $ass_report_dir = pop(@ftp_path_array);
 # $ass_report_dir = pop(@ftp_path_array);
 # $ass_report_dir = pop(@ftp_path_array);

 # $cmd = "wget ".$wget_verbose." -nH -P ".$local_dir." ".$ftp_path."/../../".$ass_report_dir."_assembly_report.txt -O ".$local_dir."/assembly_report.txt";
  $cmd = "wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/006/655/GCA_900006655.1_GSMRT3/GCA_900006655.1_GSMRT3_assembly_report.txt"
  if (system($cmd)) {
    $self->throw("Could not download *_assembly_report.txt file to ".$local_dir.". Please, check that both paths are valid.\n".
                 "Commandline used:\n".$cmd);
  }
  else {
    say "Assembly report file was downloaded";
  }
}

sub download_non_nuclear_ftp_dir {
  my ($self,$ftp_path,$local_dir) = @_;
  my $non_nuclear = "/non-nuclear/assembled_chromosomes/AGP/";
  my $non_nuclear_test = "/non-nuclear/";
  my $full_ftp_path = $ftp_path.$non_nuclear;
  my $full_local_dir_path = $local_dir."/non_nuclear_agp/";
  my $non_nuclear_http_check = $ftp_path.$non_nuclear_test ;


  use LWP::Simple qw($ua head);
  $ua->timeout(100);
 
  if (head($non_nuclear_http_check)) {
    # my $cmd = "wget -P ".$full_local_dir_path." ".$full_ftp_path."*.agp.gz";
    $non_nuclear_http_check =~ s/ftp\:\/\/ftp.ncbi.nlm.nih.gov\//ftp.ncbi.nlm.nih.gov\:\:/g;
    my $cmd = "rsync $non_nuclear_http_check/*/AGP/*agp.gz  -P $full_local_dir_path"; 
	

    if (system($cmd)) {
      	$self->throw("No non-nuclear AGP files found, but there is a non-nuclear dir.  Will not continue. Check me! $non_nuclear_http_check \n".
                   "Commandline used:\n".$cmd);                   
      	return;
    } 
    else {
      	say "Downloaded non-nuclear AGP files (can be used to identify mito later in contigs.fa)";
    }
  	$self->unzip($full_local_dir_path);
  	$self->concat_non_nuclear($full_local_dir_path);
  } 
  
  else {
    	print "non-nuclear dir/ftp does not exist or timeout, so nothing to download\n";;
  	}
}


sub unzip {
  my ($self,$path) = @_;
  say "Unzipping the compressed files...$path";
  $self->throw("gunzip operation failed. Please check your error log file.") if (system("gunzip -r $path") == 1);
  say "Unzipping finished!";
}

sub concat_non_nuclear {
  my ($self,$path) = @_;
  my $cmd = "cat ".$path."/*.agp > ".$path."/concat_non_nuclear.agp";
  my $filename = $path."/concat_non_nuclear.agp"; 
  if (!-e $filename) {
  	print "DEBUG:::The file does not exist!";
  	if (system($cmd)) {
    	$self->throw("Problem concatenation the non nuclear AGP files into concat_non_nuclear.agp\n".
                "Commandline used:\n".$cmd);
  		}
  else{
  	print $filename.": File exists, cannot concat";
  }	
	}
}


sub write_output {
  my $self = shift;

  return 1;
}

1;
