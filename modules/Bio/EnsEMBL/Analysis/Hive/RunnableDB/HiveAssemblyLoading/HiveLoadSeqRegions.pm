#!/usr/bin/env perl

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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadSeqRegions;

use strict;
use warnings;
use feature 'say';

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;

  unless ( $self->param('target_db') ) {
    $self->throw( "target_db flag not passed into parameters hash. The target db to load the assembly info " .
                  "into must be passed in with write access" );
  }

  unless ( $self->param('primary_assembly_dir_name') ) {
    $self->throw("primary_assembly_dir_name flag not passed into parameters hash. This is usually Primary_Assembly ");
  }

  unless ( $self->param('output_path') ) {
    $self->throw( "output_path flag not passed into parameters hash. This should be the path to the working directory " .
                  "that you downloaded the ftp files to earlier in the pipeline" );
  }

  unless ( $self->param('enscode_root_dir') ) {
    $self->throw("enscode_dir flag not passed into parameters hash. You need to specify where your code checkout is");
  }

  unless ( $self->param('coord_system_version') ) {
    $self->throw("coord_system_version flag not passed into parameters hash. You need to specify the assembly version e.g. GRCh38");
  }

  return 1;

} ## end sub fetch_input

sub run {
  my $self = shift;

  say "Loading seq regions into reference db";

  $self->concat_files_for_loading();

  my $assembly_loading_style = $self->determine_assembly_style();
  $self->load_seq_regions($assembly_loading_style);

  say "Finished loading seq regions and AGP info";

  return 1;
}

sub write_output {
  my $self = shift;

  return 1;
}

sub concat_files_for_loading {
  my ($self) = @_;

  my $path_to_files = $self->param('output_path') . "/" . $self->param('species_name') . "/" . $self->param('primary_assembly_dir_name');

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
    for i in `ls $path_to_files/AGP/*.comp.agp`; do
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

} ## end sub concat_files_for_loading

sub determine_assembly_style {
  my ($self) = @_;

  my $path_to_files = $self->param('output_path') . "/" . $self->param('species_name') . "/" . $self->param('primary_assembly_dir_name');
  my $coord_system_version = $self->param('coord_system_version');
  my $assembly_level       = $self->param('assembly_level');

  my $assembly_loading_style;
  my $chr_all_pres  = 0;
  my $comp_all_pres = 0;
  my $scaf_all_pres = 0;

  if ( -e $path_to_files . "/AGP/chr_all.agp" ) {
    $chr_all_pres = 1;
  }

  if ( -e $path_to_files . "/AGP/scaf_all.agp" ) {
    $scaf_all_pres = 1;
  }

  if ( -e $path_to_files . "/AGP/comp_all.agp" ) {
    $comp_all_pres = 1;
  }

  # This one should represent things like salmon, where there are 'contigs' that are effectively scaffolds
  # and only a mapping from chromosome to component
  if ( !$chr_all_pres && $comp_all_pres && $assembly_level eq 'chromosome' ) {
    $assembly_loading_style = "pseudocontig";
  }

  # This if for the traditional chromosome, scaffold, contig model
  elsif ( $chr_all_pres && $scaf_all_pres && $assembly_level eq 'chromosome' ) {
    $assembly_loading_style = "traditional";
  }

  # This is for the traditional scaffold, contig model
  elsif ( !$chr_all_pres && $scaf_all_pres && $assembly_level eq 'scaffold' ) {
    $assembly_loading_style = "traditional";
  }

  else {
    $self->throw( "The file structure of this assembly is not covered in the code" .
                 "\nCoordSystemVersion: " . $coord_system_version . "\nAssembly level: " . $assembly_level . "\nChomosome AGP present: " .
                 $chr_all_pres . "\nScaffold AGP present: " . $scaf_all_pres . "\nComponent AGP present: " . $comp_all_pres );
  }

  return $assembly_loading_style;
} ## end sub determine_assembly_style

sub load_seq_regions {
  my ( $self, $assembly_loading_style ) = @_;

  my $target_db     = $self->param('target_db');
  my $path_to_files = $self->param('output_path') . "/" . $self->param('species_name') . "/" . $self->param('primary_assembly_dir_name');
  my $enscode_dir   = $self->param('enscode_root_dir');
  my $coord_system_version = $self->param('coord_system_version');
  my $assembly_level       = $self->param('assembly_level');

  my $dbhost = $target_db->{'-host'};
  my $dbport = $target_db->{'-port'};
  my $dbuser = $target_db->{'-user'};
  my $dbpass = $target_db->{'-pass'};
  my $dbname = $target_db->{'-dbname'};

  unless ( $dbhost && $dbuser && $dbpass && $dbname ) {
    $self->throw("Database connections info not fully present. Must pass in using the target_db flag");
  }

  my $contig_rank     = 0;
  my $scaffold_rank   = 0;
  my $chromosome_rank = 0;
  my $chromosome_agp;
  my $scaffold_agp;
  if ( $assembly_loading_style eq 'pseudocontig' && $assembly_level eq 'chromosome' ) {
    $contig_rank     = 2;
    $chromosome_rank = 1;
    $chromosome_agp  = $path_to_files . "/AGP/comp_all.agp";
    say "Pseudocontig chromosome level AGP structure found. Will load pseudocontigs as rank 2 and chromosomes as rank 1";
  }
  elsif ( $assembly_loading_style eq 'traditional' ) {
    if ( $assembly_level eq 'chromosome' ) {
      $contig_rank     = 3;
      $scaffold_rank   = 2;
      $chromosome_rank = 1;
      $scaffold_agp    = $path_to_files . "/AGP/scaf_all.agp";
      $chromosome_agp  = $path_to_files . "/AGP/chr_all.agp";
      say "Traditional chromosome level AGP structure found. Will load contigs as rank 3, scaffolds as rank 2 and chromosomes as rank 1";
    }
    elsif ( $assembly_level eq 'scaffold' ) {
      $contig_rank   = 2;
      $scaffold_rank = 1;
      $scaffold_agp  = $path_to_files . "/AGP/scaf_all.agp";
      say "Traditional scaffold level AGP structure found. Will load contigs as rank 2 and scaffolds as rank 1";
    }
  }

  my $contigs_file_path = $path_to_files . "/contigs/contigs.fa";
  unless ( -e $contigs_file_path ) {
    $self->throw( "Could not locate contigs file. Expected to find here:\n" . $contigs_file_path );
  }

  # use the load seq_regions script to load the contigs
  say "\nLoading the contigs...";
  say "Processing file:\n" . $contigs_file_path;

  my $cmd =
    "perl " . $enscode_dir . "/ensembl-pipeline/scripts/load_seq_region.pl" .
    " -dbhost " . $dbhost . " -dbuser " . $dbuser . " -dbpass " . $dbpass . " -dbport " . $dbport . " -dbname " . $dbname .
    " -coord_system_name contig" . " -coord_system_version " . $coord_system_version . " -rank " . $contig_rank . " -default_version" .
    " -fasta_file " . $contigs_file_path . " -sequence_level" . " -noverbose" . " > " . $path_to_files . "/load_seq_region_contigs.out";

  my $result = system($cmd);
  if ($result) {
    $self->throw( "The load_seq_regions script returned a non-zero exit code when loading the contigs. Commandline used:\n" . $cmd );
  }

  say "The output  was written to:\n" . $path_to_files . "/load_seq_region_contigs.out";

  # check contigs
  my $num_contigs = int(`mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e'select count(*) from seq_region'`);
  my $num_dna     = int(`mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e'select count(*) from dna'`);

  if ( $num_contigs != $num_dna ) {
    $self->throw("The number of 'seq_region' table loaded contigs ($num_contigs) differs from the number of 'dna' table rows ($num_dna)");
  }
  else {
    say "The number of 'seq_region' table loaded contigs (" .
      $num_contigs . ") is the same as the number of 'dna' table rows ($num_dna). Great!";
  }

  # use the load seq_regions script to load the scaffolds
  # for each existing scaffold agp file (unplaced, placed, unlocalized)
  if ( $assembly_loading_style eq 'traditional' ) {
    say "\nLoading the scaffolds...";

    $cmd =
      "perl " . $enscode_dir . "/ensembl-pipeline/scripts/load_seq_region.pl" .
      " -dbhost " . $dbhost . " -dbuser " . $dbuser . " -dbpass " . $dbpass . " -dbport " . $dbport . " -dbname " .
      $dbname . " -coord_system_name scaffold" . " -coord_system_version " . $coord_system_version . " -rank " . $scaffold_rank .
      " -default_version" . " -agp_file " . $scaffold_agp . " -noverbose" . " > " . $path_to_files . "/load_seq_region_scaffolds.out";

    $result = system($cmd);
    if ($result) {
      $self->throw( "The load_seq_regions script returned a non-zero exit code when loading the contigs. Commandline used:\n" . $cmd );
    }

    say "The output  was written to:\n" . $path_to_files . "/load_seq_region_scaffolds.out\n";
  }

  # If we have chromosomes
  if ( $assembly_level eq 'chromosome' ) {
    say "Loading the chromosomes...";

    $cmd =
      "perl " . $enscode_dir . "/ensembl-pipeline/scripts/load_seq_region.pl" .
      " -dbhost " . $dbhost . " -dbuser " . $dbuser . " -dbpass " . $dbpass . " -dbport " . $dbport . " -dbname " .
      $dbname . " -coord_system_name chromosome" . " -coord_system_version " . $coord_system_version . " -rank " . $chromosome_rank .
      " -default_version" . " -agp_file " . $chromosome_agp . " -noverbose" . " > " . $path_to_files . "/load_seq_region_chromosomes.out";

    $result = system($cmd);
    if ($result) {
      $self->throw( "The load_seq_regions script returned a non-zero exit code when loading the contigs. Commandline used:\n" . $cmd );
    }

    say "The output  was written to:\n" . $path_to_files . "/load_seq_region_scaffolds.out";

    # maybe check for exceptions in the output file here

    # maybe check fasta lengths here
    # maybe check warnings about non-ATGCN (RYKMSWBDHV) bases
  }

  # ought to update contig version to NULL
  `mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -e'update coord_system set version=NULL where name = "contig"'`;

  # check that it was updated
  my $contig_version =
    `mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -e'select version from coord_system where name = "contig"'`;

  if ( $contig_version =~ /NULL/ ) {
    say "Column 'version' in 'coord_system' table set to NULL for contig level as required";
  }
  else {
    $self->throw("Column 'version' in 'coord_system' table could not be set to NULL for contig level");
  }

} ## end sub load_seq_regions
1;
