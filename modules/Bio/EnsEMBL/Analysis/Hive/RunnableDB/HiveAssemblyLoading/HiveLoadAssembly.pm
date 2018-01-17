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

use File::Spec::Functions qw(catfile catdir);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveLoadSeqRegions');

sub fetch_input {
  my $self = shift;

  unless($self->param('target_db')) {
    $self->throw("target_db flag not passed into parameters hash. The target db to load the assembly info ".
                 "into must be passed in with write access");
  }

  $self->throw("Could not find the load_agp script in enscode dir. Path checked:\n". catfile($self->param('enscode_root_dir'), 'ensembl-pipeline', 'scripts', 'load_agp.pl'))
    unless(-e catfile($self->param_required('enscode_root_dir'), 'ensembl-pipeline', 'scripts', 'load_agp.pl'));

  $self->param_required('output_path');
  return 1;
}

sub run {
  my $self = shift;

  my $path_to_files = $self->param('output_path');
  if ($self->param_is_defined('primary_assembly_dir_name')) {
    $path_to_files = catdir($self->param('output_path'), $self->param_required('primary_assembly_dir_name'), 'AGP');
  }

  unless(-e $path_to_files) {
    $self->warning("The AGP directory was not found. Assuming the assembly is single level. Path checked:\n".$path_to_files);
    return;
  }

  say "Loading seq regions into reference db";
  my $target_db = $self->param('target_db');
  my $dbhost = $target_db->{'-host'};
  my $dbport = $target_db->{'-port'};
  my $dbuser = $target_db->{'-user'};
  my $dbpass = $target_db->{'-pass'};
  my $dbname = $target_db->{'-dbname'};

  my $chromo_present = 0;

  my %assembled = (
    $self->param('scaffold_contig') => 'scaffold',
    $self->param('chromosome_contig') => 'chromosome',
    $self->param('chromosome_scaffold')  => 'chromosome',
  );
  my %component = (
    $self->param('scaffold_contig') => 'contig',
    $self->param('chromosome_contig') => 'contig',
    $self->param('chromosome_scaffold')  => 'scaffold',
  );
  my $sql_count = "mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e 'SELECT COUNT(*) FROM assembly'";
  my $base_cmd = 'perl '.catfile($self->param('enscode_root_dir'), 'ensembl-pipeline', 'scripts', 'load_agp.pl').
                ' -dbhost '.$dbhost.
                ' -dbuser '.$dbuser.
                ' -dbpass '.$dbpass.
                ' -dbport '.$dbport.
                ' -dbname '.$dbname;
  foreach my $filename (keys %assembled) {
    my $file = catfile($path_to_files, $filename);
    if (-e $file) {
      ++$chromo_present if ($assembled{$filename} eq 'chromosome');
      say 'Loading '.$assembled{$filename}.' to '.$component{$filename}." mappings from:\n$file";

      my $assembly_count_before = int(`$sql_count`);

      my $cmd = $base_cmd.' -assembled_name '.$assembled{$filename}.
                ' -component_name '.$component{$filename}.
                ' -agp_file '.$file;

      my $num_warnings = 0;
      my $assembly_count_file = 0;
      open(CMD, $cmd.' 2>&1 | ') || $self->throw("Could not execute $cmd");
      while (<CMD>) {
        if (/WARN/) {
          ++$num_warnings;
        }
        elsif (/You are already using/) {
          ++$num_warnings;
        }
        elsif (/EXCEPTION/) {
          $self->throw('Something went wrong!');
        }
      }
      close(CMD) || $self->throw("$cmd failed");

#      my $assembly_count_after = int(`$sql_count`);
#      my $assembly_count_diff = $assembly_count_after - $assembly_count_before;
#      if ($assembly_count_file != $assembly_count_diff) {
#        $self->throw("Difference in count between $filename and entries loaded into the database:\n".
#              "$filename: ".$assembly_count_file."\ndatabase: ".$assembly_count_diff);
#      }
#
#      if ($num_warnings > 0) {
#        $self->warning("There are some warning messages in the output file that should be checked (grep 'You are already using'):\n$file");
#      }
    }
  }

  # Have to ask about this scenario and whether it warrants a warning or a throw. It may not
  # actually be needed. But if it is then I need to make it more useful
  if($chromo_present) {
    if ($chromo_present != 2) {
      $self->warning("Chromosomes present, but you only have $chromo_present mapping file instead of two");
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
  say "Finished loading the AGP files";
  return 1;
}

sub write_output {
  my $self = shift;

  return 1;
}

1;
