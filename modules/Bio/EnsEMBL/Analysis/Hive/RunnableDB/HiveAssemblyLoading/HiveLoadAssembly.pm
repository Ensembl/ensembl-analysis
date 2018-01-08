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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadAssembly;

use strict;
use warnings;
use feature 'say';

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;

  unless($self->param('target_db')) {
    $self->throw("target_db flag not passed into parameters hash. The target db to load the assembly info ".
                 "into must be passed in with write access");
  }

  unless($self->param('primary_assembly_dir_name')) {
    $self->throw("primary_assembly_dir_name flag not passed into parameters hash. This is usually Primary_Assembly ");
  }

  unless($self->param('output_path')) {
    $self->throw("output_path flag not passed into parameters hash. This should be the path to the working directory ".
                 "that you downloaded the ftp files to earlier in the pipeline");
  }

  unless($self->param('enscode_root_dir')) {
    $self->throw("enscode_dir flag not passed into parameters hash. You need to specify where your code checkout is");
  }

  return 1;
}

sub run {
  my $self = shift;

  say "Loading seq regions into reference db";
  my $target_db = $self->param('target_db');
  my $path_to_files = $self->param('output_path')."/".$self->param('primary_assembly_dir_name')."/AGP/";
  my $enscode_dir = $self->param('enscode_root_dir');

  $self->load_assembly($target_db,$path_to_files,$enscode_dir);

  say "Finished downloading contig files";
  return 1;
}

sub write_output {
  my $self = shift;

  return 1;
}

sub load_assembly {
  my ($self,$target_db,$path_to_files,$enscode_dir) = @_;

  my $dbhost = $target_db->{'-host'};
  my $dbport = $target_db->{'-port'};
  my $dbuser = $target_db->{'-user'};
  my $dbpass = $target_db->{'-pass'};
  my $dbname = $target_db->{'-dbname'};

  unless(-e $enscode_dir."/ensembl-pipeline/scripts/load_agp.pl") {
    $self->throw("Could not find the load_agp script in enscode dir. Path checked:\n".$enscode_dir."/ensembl-pipeline/scripts/load_agp.pl");
  }

  my $chromo_present = 0;

  my $assembled;
  my $component;
  if(-e $path_to_files."/scaf_all.agp") {
    say "Loading scaffold to contig mappings from:\n".$path_to_files."/scaf_all.agp";
    $assembled = "scaffold";
    $component = "contig";

    my $assembly_count_before = int(`mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e 'select count(*) from assembly'`);

    my $cmd = "perl ".$enscode_dir."/ensembl-pipeline/scripts/load_agp.pl".
              " -dbhost ".$dbhost.
              " -dbuser ".$dbuser.
              " -dbpass ".$dbpass.
              " -dbport ".$dbport.
              " -dbname ".$dbname.
              " -assembled_name ".$assembled.
              " -component_name ".$component.
              " -agp_file ".$path_to_files."/scaf_all.agp".
              "  > ".$path_to_files."/load_assembly_scaffold_contig.out";
    my $return = system($cmd);

    if($return) {
      $self->throw("The load agp script returned a non-zero exit code. Commandline used:\n".$cmd);
    }

    say "Loading of scaffold to contig mappings completed. Output written to:\n".$path_to_files."/load_assembly_scaffold_contig.out\n";

    my $assembly_count_after = int(`mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e 'select count(*) from assembly'`);
    my $assembly_count_diff = $assembly_count_after - $assembly_count_before;
    my $assembly_count_file = int(`grep -v "#" $path_to_files/scaf_all.agp | grep -v "\\WN\\W" | grep -c -v "\\WU\\W"`);
    unless($assembly_count_file == $assembly_count_diff) {
      $self->throw("Difference in count between scaff_all.agp and entries loaded into the database:\n".
            "scaff_all.agp: ".$assembly_count_file."\ndatabase: ".$assembly_count_diff);
    }

    # Check output file for any issues
    my $num_warnings = int(`grep "WARN" $path_to_files/scaf_all.agp | wc -l`);
    if ($num_warnings > 0) {
      $self->warning("There are some warning messages in the output file that should be checked (grep 'WARN'):\n".$path_to_files."/scaf_all.agp");
    }

    $num_warnings = int(`grep "You are already using" $path_to_files/scaf_all.agp | wc -l`);
    if ($num_warnings > 0) {
      $self->warning("There are some warning messages in the output file that should be checked (grep 'You are already using'):\n".$path_to_files."/scaf_all.agp");
    }

    my $exception = int(`grep "EXCEPTION" $path_to_files/scaf_all.agp | wc -l`);
    if($exception) {
      $self->throw("Exception found in output file:\n".$path_to_files."/scaf_all.agp");
    }

  } else {
    $self->throw("No scaffold to contig agp file exists, expected the following file:\n".$path_to_files."/scaf_all.agp");
  }

  # Need to rewrite the following at some point to just be a call to a single subroutine with the file name, assembled and component
  # vars passed in. At the moment it's needlessly duplicated code
  if(-e $path_to_files."/comp_all.agp") {
    say "Loading chromosome to contig mappings from:\n".$path_to_files."/comp_all.agp";
    $assembled = "chromosome";
    $component = "contig";
    $chromo_present++;

    my $assembly_count_before = int(`mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e 'select count(*) from assembly'`);

    my $cmd = "perl ".$enscode_dir."/ensembl-pipeline/scripts/load_agp.pl".
              " -dbhost ".$dbhost.
              " -dbuser ".$dbuser.
              " -dbpass ".$dbpass.
              " -dbport ".$dbport.
              " -dbname ".$dbname.
              " -assembled_name ".$assembled.
              " -component_name ".$component.
              " -agp_file ".$path_to_files."/comp_all.agp".
              "  > ".$path_to_files."/load_assembly_chromosome_contig.out";
    my $return = system($cmd);

    if($return) {
      $self->throw("The load agp script returned a non-zero exit code. Commandline used:\n".$cmd);
    }

    say "Loading of chromosome to contig mappings completed. Output written to:\n".$path_to_files."/load_assembly_chromosome_contig.out\n";

    my $assembly_count_after = int(`mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e 'select count(*) from assembly'`);
    my $assembly_count_diff = $assembly_count_after - $assembly_count_before;
    my $assembly_count_file = int(`grep -v "#" $path_to_files/comp_all.agp | grep -v "\\WN\\W" | grep -c -v "\\WU\\W"`);
    unless($assembly_count_file == $assembly_count_diff) {
      $self->throw("Difference in count between comp_all.agp and entries loaded into the database:\n".
            "comp_all.agp: ".$assembly_count_file."\ndatabase: ".$assembly_count_diff);
    }

    # Check output file for any issues
    my $num_warnings = int(`grep "WARN" $path_to_files/comp_all.agp | wc -l`);
    if ($num_warnings > 0) {
      $self->warning("There are some warning messages in the output file that should be checked (grep 'WARN'):\n".$path_to_files."/comp_all.agp");
    }

    $num_warnings = int(`grep "You are already using" $path_to_files/comp_all.agp | wc -l`);
    if ($num_warnings > 0) {
      $self->warning("There are some warning messages in the output file that should be checked (grep 'You are already using'):\n".$path_to_files."/comp_all.agp");
    }

    my $exception = int(`grep "EXCEPTION" $path_to_files/comp_all.agp | wc -l`);
    if($exception) {
      $self->throw("Exception found in output file:\n".$path_to_files."/comp_all.agp");
    }

  }

  if(-e $path_to_files."/chr_all.agp") {
    say "Loading chromosome to scaffold mappings from:\n".$path_to_files."chr_all.agp";
    $assembled = "chromosome";
    $component = "scaffold";
    $chromo_present++;

    my $assembly_count_before = int(`mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e 'select count(*) from assembly'`);

    my $cmd = "perl ".$enscode_dir."/ensembl-pipeline/scripts/load_agp.pl".
              " -dbhost ".$dbhost.
              " -dbuser ".$dbuser.
              " -dbpass ".$dbpass.
              " -dbport ".$dbport.
              " -dbname ".$dbname.
              " -assembled_name ".$assembled.
              " -component_name ".$component.
              " -agp_file ".$path_to_files."/chr_all.agp".
              "  > ".$path_to_files."/load_assembly_chromosome_scaffold.out";
    my $return = system($cmd);

    if($return) {
      $self->throw("The load agp script returned a non-zero exit code. Commandline used:\n".$cmd);
    }

    say "Loading of chromosome to scaffold mappings completed. Output written to:\n".$path_to_files."/load_assembly_chromosome_scaffold.out\n";

    my $assembly_count_after = int(`mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e 'select count(*) from assembly'`);
    my $assembly_count_diff = $assembly_count_after - $assembly_count_before;
    my $assembly_count_file = int(`grep -v "#" $path_to_files/chr_all.agp | grep -v "\\WN\\W" | grep -c -v "\\WU\\W"`);
    unless($assembly_count_file == $assembly_count_diff) {
      $self->throw("Difference in count between chr_all.agp and entries loaded into the database:\n".
            "chr_all.agp: ".$assembly_count_file."\ndatabase: ".$assembly_count_diff);
    }

    # Check output file for any issues
    my $num_warnings = int(`grep "WARN" $path_to_files/chr_all.agp | wc -l`);
    if ($num_warnings > 0) {
      $self->warning("There are some warning messages in the output file that should be checked (grep 'WARN'):\n".$path_to_files."/chr_all.agp");
    }

    $num_warnings = int(`grep "You are already using" $path_to_files/chr_all.agp | wc -l`);
    if ($num_warnings > 0) {
      $self->warning("There are some warning messages in the output file that should be checked (grep 'You are already using'):\n".$path_to_files."/chr_all.agp");
    }

    my $exception = int(`grep "EXCEPTION" $path_to_files/chr_all.agp | wc -l`);
    if($exception) {
      $self->throw("Exception found in output file:\n".$path_to_files."/chr_all.agp");
    }

  }

  # Have to ask about this scenario and whether it warrants a warning or a throw. It may not
  # actually be needed. But if it is then I need to make it more useful
  if($chromo_present) {
    unless($chromo_present == 2) {
      $self->warning("Chromosomes present, but you only have one mapping file instead of two");
    }
  }

  # check meta table new rows for assembly.mapping
  my $num_meta_rows = int(`mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e'select count(*) from meta where meta_key = "assembly.mapping"'`);
  if ($chromo_present) {
    if ($num_meta_rows == 3) {
      say "3 new rows found in 'meta' table for 'assembly.mapping' meta_key. Great!";
    } else {
      $self->warning($num_meta_rows." new rows found in 'meta' table for 'assembly.mapping' meta_key. 3 expected.");
    }
  } else { # if there is no chromosome level
    if ($num_meta_rows == 1) {
      say "1 new row found in 'meta' table for 'assembly.mapping' meta_key. Great!";
    } else {
      $self->throw($num_meta_rows." new rows found in 'meta' table for 'assembly.mapping' meta_key. 1 expected");
    }
  }

}
1;
