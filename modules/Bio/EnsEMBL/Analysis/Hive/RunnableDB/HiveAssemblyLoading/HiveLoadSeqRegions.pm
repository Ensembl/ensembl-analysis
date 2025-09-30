#!/usr/bin/env perl

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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadSeqRegions;

use strict;
use warnings;
use feature 'say';

use File::Find;
use File::Spec::Functions qw(catfile catdir);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    contig_file => 'contigs.fa',
    scaffold_contig => 'scaf_all.agp',
    chromosome_contig => 'comp_all.agp',
    chromosome_scaffold => 'chr_all.agp',
    replace_ambiguous_bases => 0,
    ignore_ambiguous_bases => 0,
  }
}


sub fetch_input {
  my $self = shift;

  unless($self->param('target_db')) {
    $self->throw("target_db flag not passed into parameters hash. The target db to load the assembly info ".
                 "into must be passed in with write access");
  }

  unless($self->param('output_path')) {
    $self->throw("output_path flag not passed into parameters hash. This should be the path to the working directory ".
                 "that you downloaded the ftp files to earlier in the pipeline");
  }

  unless($self->param('enscode_root_dir')) {
    $self->throw("enscode_dir flag not passed into parameters hash. You need to specify where your code checkout is");
  }

  if ($self->param_is_defined('assembly_name') and !$self->param_is_defined('coord_system_version')) {
    $self->param('coord_system_version', $self->param('assembly_name'));
  }
  unless($self->param('coord_system_version')) {
    $self->throw("coord_system_version flag not passed into parameters hash. You need to specify the assembly version e.g. GRCh38");
  }

  return 1;

}

sub run {
  my $self = shift;

  say "Loading seq regions into reference db";
  my $target_db = $self->param('target_db');
  my $path_to_files = $self->param('output_path');
  if ($self->param_is_defined('primary_assembly_dir_name')) {
    $path_to_files = catdir($self->param('output_path'), $self->param_required('primary_assembly_dir_name'));
    if(-e $path_to_files."/AGP/") {
      $self->concat_files_for_loading($path_to_files);
    }
  }
  my $enscode_dir = $self->param('enscode_root_dir');
  my $coord_system_version = $self->param('coord_system_version');

  $self->load_seq_regions($target_db,$path_to_files,$enscode_dir,$coord_system_version);

  say "Finished downloading contig files";
  return 1;
}

sub write_output {
  my $self = shift;

  return 1;
}

sub concat_files_for_loading {
  my ($self,$path_to_files) = @_;

  # remove possibly (previous executions) created agp files to avoid inappropriate concatenation
  `rm $path_to_files/AGP/*_all.agp`;
  `rm $path_to_files/FASTA/*_all.agp.fastawc`;

  # Copy the agp files to their own subdirs to separate for concat
  system(
<<COMMAND
    for i in `ls $path_to_files/AGP/*.agp | grep -v '.comp.agp' | grep -v '.scaf.agp'`; do
      cat \$i >> $path_to_files/AGP/chr_all.agp
    done
    for i in `ls $path_to_files/AGP/*.scaf.agp`; do
      cat \$i >> $path_to_files/AGP/scaf_all.agp
    done
    for i in `ls $path_to_files/AGP/chr*.comp.agp`; do
      cat \$i >> $path_to_files/AGP/comp_all.agp
    done
COMMAND
        );

  # Also need to take out the headers for checking loading later
  system(
<<COMMAND
    for i in `ls $path_to_files/FASTA/*.fna | grep -v '.scaf.fna'`; do
      grep '>' \$i >> $path_to_files/FASTA/chr_all.agp.fastawc
    done
    for i in `ls $path_to_files/FASTA/*.scaf.fna`; do
      grep '>' \$i >> $path_to_files/FASTA/scaf_all.agp.fastawc
    done
COMMAND
        );

}

sub load_seq_regions {
  my ($self,$target_db,$path_to_files,$enscode_dir,$coord_system_version) = @_;

  my $dbhost = $target_db->{'-host'};
  my $dbport = $target_db->{'-port'};
  my $dbuser = $target_db->{'-user'};
  my $dbpass = $target_db->{'-pass'};
  my $dbname = $target_db->{'-dbname'};

  unless($dbhost && $dbuser && $dbpass && $dbname) {
    $self->throw("Database connections info not fully present. Must pass in using the target_db flag");
  }

  # Check if chromosomes exist
  my $chromo_present = 1;
  my $single_level = 0;
  my $rank = 3;
  my $middle_dir = $self->param_is_defined('primary_assembly_dir_name') ? 'AGP' : '';
  my $file_chr = catfile($path_to_files, $middle_dir, $self->param('chromosome_scaffold'));

  unless(-e $path_to_files."/AGP/") {
    $self->warning("Path to AGP dir does not exist. Assuming single level assembly. Will set rank to 2 for contig. ".
                   "Path checked:\n".$path_to_files."/AGP/");
    $rank = 2;
    $single_level = 1;
  } elsif(-e $file_chr) {
    say "Chromosome agp file found. Chromosome level will be of rank 1, scaffold level will be of rank 2 and contig level will be of rank 3";
  }
  elsif (-e catfile($path_to_files, $middle_dir, $self->param('chromosome_contig'))) {
    $file_chr = catfile($path_to_files, $middle_dir, $self->param('chromosome_contig'));
    say "Chromosome agp file found. Chromosome level will be of rank 1, scaffold level will be of rank 2 and contig level will be of rank 3";
  }
  else {
    say "Chromosome agp file not found. Scaffold level will be of rank 1 and contig level will be of rank 2";
    $rank = 2;
    $chromo_present = 0;
  }

  $middle_dir = $self->param_is_defined('primary_assembly_dir_name') ? 'contigs' : '';
  my $contigs_file_path = catfile($path_to_files, $middle_dir, $self->param('contig_file'));
  unless (-e $contigs_file_path) {
    $self->throw("Could not locate contigs file. Expected to find here:\n".$contigs_file_path);
  }

  my $base_cmd = 'perl '.catfile($enscode_dir, 'ensembl-analysis', 'scripts', 'assembly_loading', 'load_seq_region.pl').
            ' -dbhost '.$dbhost.
            ' -dbuser '.$dbuser.
            ' -dbpass '.$dbpass.
            ' -dbport '.$dbport.
            ' -dbname '.$dbname.
            ' -coord_system_version '.$coord_system_version.
            ' -default_version'.
            ' -noverbose';
  my $base_sql = "mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname";
  # use the load seq_regions script to load the contigs
  say "\nLoading the contigs...";
  say "Processing file:\n".$contigs_file_path;

  my $cmd = $base_cmd.
            " -coord_system_name contig".
            " -rank ".$rank.
            " -fasta_file ".$contigs_file_path.
            " -sequence_level";

  if($self->param('replace_ambiguous_bases')) {
    $cmd .= " -replace_ambiguous_bases ";
  } elsif($self->param('ignore_ambiguous_bases')) {
    $cmd .= " -ignore_ambiguous_bases ";
  }

  my $result = system($cmd);
  if($result) {
    $self->throw("The load_seq_regions script returned a non-zero exit code when loading the contigs. Commandline used:\n".$cmd);
  }

  # check contigs
  my $num_contigs = int(`$base_sql -NB -e'select count(*) from seq_region'`);
  my $num_dna = int(`$base_sql -NB -e'select count(*) from dna'`);

  if ($num_contigs != $num_dna) {
    $self->throw("The number of 'seq_region' table loaded contigs ($num_contigs) differs from the number of 'dna' table rows ($num_dna)");
  } else {
    say "The number of 'seq_region' table loaded contigs (".$num_contigs.") is the same as the number of 'dna' table rows ($num_dna). Great!";
  }

  # ought to update contig version to NULL
  `$base_sql -e'update coord_system set version=NULL where name = "contig"'`;

  # check that it was updated
  my $contig_version = `$base_sql -e'select version from coord_system where name = "contig"'`;

  if ($contig_version =~ /NULL/) {
    say "Column 'version' in 'coord_system' table set to NULL for contig level as required";
  } else {
    $self->throw("Column 'version' in 'coord_system' table could not be set to NULL for contig level");
  }

  if($single_level) {
    # At this point we want to stop for single level assemblies
    return;
  }

  $rank--;

  # use the load seq_regions script to load the scaffolds
  # for each existing scaffold agp file (unplaced, placed, unlocalized)
  say "\nLoading the scaffolds...";

  $middle_dir = $self->param_is_defined('primary_assembly_dir_name') ? 'AGP' : '';
  $cmd = $base_cmd.
         " -coord_system_name scaffold".
         " -rank ".$rank.
         " -agp_file ".catfile($path_to_files, $middle_dir, $self->param('scaffold_contig'));

  $result = system($cmd);
  if($result) {
    $self->throw("The load_seq_regions script returned a non-zero exit code when loading the contigs. Commandline used:\n".$cmd);
  }

  $rank--;

  # If we have chromosomes
  if ($rank > 0) {
    say "Loading the chromosomes...";
    my $chromo_mapping_file;

    $cmd = $base_cmd.
           " -coord_system_name chromosome".
           " -rank ".$rank.
           " -agp_file ".$file_chr;

  $result = system($cmd);
  if($result) {
    $self->throw("The load_seq_regions script returned a non-zero exit code when loading the contigs. Commandline used:\n".$cmd);
  }


    # maybe check for exceptions in the output file here

    # maybe check fasta lengths here
    # maybe check warnings about non-ATGCN (RYKMSWBDHV) bases
  }



}
1;
