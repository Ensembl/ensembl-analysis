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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveSetAndCheckToplevel;

use strict;
use warnings;
use feature 'say';

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;

  unless ( $self->param('target_db') ) {
    $self->throw( "target_db flag not passed into parameters hash. The core db to load the assembly info " .
                  "into must be passed in with write access" );
  }

  unless ( $self->param('enscode_root_dir') ) {
    $self->throw( "enscode_root_dir flag not passed into parameters hash. You need to specify where your code checkout is" );
  }

  return 1;
}

sub run {
  my $self = shift;

  say "Loading seq regions into reference db";
  my $target_db                 = $self->param('target_db');
  my $enscode_dir               = $self->param('enscode_root_dir');
  my $primary_assembly_dir_name = $self->param('primary_assembly_dir_name');
  my $path_to_files             = $self->param('output_path') . "/" . $self->param('species_name') . "/" . $primary_assembly_dir_name;

  $self->set_toplevel( $target_db, $enscode_dir );
  $self->check_toplevel( $target_db, $enscode_dir, $path_to_files );
  say "Finished downloading contig files";
  return 1;
}

sub write_output {
  my $self = shift;

  return 1;
}

sub set_toplevel {
  my ( $self, $target_db, $enscode_dir ) = @_;

  my $dbhost = $target_db->{'-host'};
  my $dbport = $target_db->{'-port'};
  my $dbuser = $target_db->{'-user'};
  my $dbpass = $target_db->{'-pass'};
  my $dbname = $target_db->{'-dbname'};

  my $cmd =
    "perl " . $enscode_dir . "/ensembl-pipeline/scripts/set_toplevel.pl " .
    " -dbhost " . $dbhost . " -dbuser " . $dbuser . " -dbpass " . $dbpass . " -dbport " . $dbport . " -dbname " . $dbname;

  my $return = system($cmd);
  if ($return) {
    $self->throw( "The set_toplevel script returned a non-zero exit code. Commandline used:\n" . $cmd );
  }

  my $num_toplevel = int(
      `mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e'select count(*) from seq_region_attrib where attrib_type_id = 6'` );
  if ( $num_toplevel > 0 ) {
    say "Found " . $num_toplevel . " toplevel (non-redundant) seq region attributes found in seq_region_attrib table (code 6)";
  }
  else {
    $self->throw(
              "No toplevel (non-redundant) seq region attributes found in seq_region_attrib table (code 6). Commandline used:\n" . $cmd );
  }

} ## end sub set_toplevel

sub check_toplevel {
  my ( $self, $target_db, $enscode_dir, $path_to_files ) = @_;
  my $check_dump_path = $path_to_files . "/check_and_dump_toplevel";

  my $dbhost = $target_db->{'-host'};
  my $dbport = $target_db->{'-port'};
  my $dbuser = $target_db->{'-user'};
  my $dbpass = $target_db->{'-pass'};
  my $dbname = $target_db->{'-dbname'};

  if ( -d $check_dump_path ) {
    `rm -r $check_dump_path`;
  }

  `mkdir $check_dump_path`;

  my $cmd = "ls " . $path_to_files . "/FASTA/*.fna | grep -v '\\.placed.' | xargs cat >> " . $check_dump_path . "/downloaded_toplevel.fa";
  my $return = system($cmd);
  if ($return) {
    $self->throw( "Concatenation of downloaded fna sequences filed. Commandline used:\n" . $cmd );
  }

  my $num_toplevel_downloaded = int(`grep -c ">" $check_dump_path/downloaded_toplevel.fa`);
  if ( $num_toplevel_downloaded == 0 ) {
    $self->throw( "Found no headers in the concatenated downloaded files. Commandline used:\n" .
                  "grep -c \">\" " . $check_dump_path . "/downloaded_toplevel.fa" );
  }

  say "\nDumping top level sequence...";
  $cmd =
    "perl " . $enscode_dir . "/ensembl-analysis/scripts/sequence_dump.pl" .
    " -dbhost " . $dbhost . " -dbuser " . $dbuser . " -dbpass " . $dbpass . " -dbport " . $dbport . " -dbname " .
    $dbname . " -toplevel " . " -onefile " . " -output_dir " . $check_dump_path . " > " . $check_dump_path . "/sequence_dump.out";

  $return = system($cmd);
  if ($return) {
    $self->throw( "Dump of toplevel sequences failed. Commandline used:\n" . $cmd );
  }

  say "Finished dumping toplevel sequences. Output written to:\n" . $check_dump_path . "/sequence_dump.out";

  my $num_toplevel_dumped = int(`grep -c ">" $check_dump_path/toplevel.fa`);
  if ( $num_toplevel_dumped == 0 ) {
    $self->throw(
             "Found no headers in the dumped toplevel file. Commandline used:\n" . "grep -c \">\" " . $check_dump_path . "/toplevel.fa" );
  }

  if ( $num_toplevel_dumped != $num_toplevel_downloaded ) {
    $self->throw( "The number of dumped top level sequences (" .
              $num_toplevel_dumped . ") does not match the number of downloaded top level sequences (" . $num_toplevel_downloaded . ")" );
  }
  else {
    say "The number of dumped top level sequences (" .
      $num_toplevel_dumped . ") matches the number of downloaded top level sequences (" . $num_toplevel_downloaded . "). Great!";
  }

  say "Combining dumped top level file and downloaded top level file into " . $check_dump_path . "/combined_toplevel.fa...\n";

  $cmd = "cat " .
    $check_dump_path . "/toplevel.fa " . $check_dump_path . "/downloaded_toplevel.fa > " . $check_dump_path . "/combined_toplevel.fa";
  $return = system($cmd);
  if ($return) {
    $self->throw( "Concatenation of dumped and downloaded toplevel sequence files failed. Commandline used:\n" . $cmd );
  }

  my $num_toplevel_combined = int(`grep -c ">" $check_dump_path/combined_toplevel.fa`);

  if ( $num_toplevel_combined != ( $num_toplevel_dumped + $num_toplevel_downloaded ) ) {
    $self->throw( "The number of combined top level sequences (" .
                  $num_toplevel_combined . ") does not match the sum of dumped top level " . "and downloaded top level sequences (" .
                  ( $num_toplevel_dumped + $num_toplevel_downloaded ) . ")" );
  }
  else {
    say "The number of combined top level sequences (" .
      $num_toplevel_combined . ") matched the sum of dumped top level " . "and downloaded top level sequences (" .
      ( $num_toplevel_dumped + $num_toplevel_downloaded ) . "). Great!";
  }

  # running fastanrdb to check that both sets of sequences are equivalent
  say "\nRunning fastanrdb...\n";

  $cmd    = "fastanrdb " . $check_dump_path . "/combined_toplevel.fa > " . $check_dump_path . "/fastanrdb.out";
  $return = system($cmd);

  if ($return) {
    $self->throw( "fastanrdb died unexpectedly. Commandline used:\n" . $cmd );
  }

  say "The output file " . $check_dump_path . "/fastanrdb.out was written";

  my $num_toplevel_fastanrdb = int(`grep -c ">" $check_dump_path/fastanrdb.out`);

  if ( $num_toplevel_fastanrdb != $num_toplevel_dumped ) {
    $self->throw( "The number of sequences in file " . $check_dump_path . "/fastanrdb.out is " .
                  $num_toplevel_fastanrdb . ", but it should be " . $num_toplevel_dumped . ". It means that we have duplicates.\n" );

    # maybe add more checks to inform about the sequences that does not match
    # cat $SCRA/fastanrdb.log | perl -ne 'next unless ($a,$b) = /^>(\w+\.\w+)\s*(\w+\.\w+)/; if ($a eq $b) {print "$a\n"}' | wc -l
    # to look at cases where NCBI and Ref DB agree
    # bsub -o $SCRA/fastanrdb_minus_r.log fastanrdb temp.fa -r
    # to compare sequences in reverse complement too

  }
  else {
    say "The number of sequences in file " .
      $check_dump_path . "/fastanrdb.out is " . $num_toplevel_fastanrdb . ", which is good since it " .
      "is the expected number, but there might be odd cases where duplicates can be present. More checks will be run";
  }

  # check that the number of identifiers in headers matches the expected number of top level sequences*2
  my $num_ids = int(`grep '>' $check_dump_path/fastanrdb.out | wc -w`);
  if ( $num_ids != $num_toplevel_dumped*2 ) {
    $self->throw( "The number of identifiers in headers in " .
                  $check_dump_path . "/fastanrdb.out (" . $num_ids . ") does not match the number of toplevel sequences * 2 (" .
                  ( $num_toplevel_dumped*2 ) . ")" );
  }
  else {
    say "The number of identifiers in headers in " .
      $check_dump_path . "/fastanrdb.out (" . $num_ids . ") matched the number of toplevel sequences * 2 (" .
      ( $num_toplevel_dumped*2 ) . "). Great!";
  }

  # check that there is no header with one identifier only
  my $num_non_matching_seq = int(`grep '>' $check_dump_path/fastanrdb.out | grep -c -v '\\w*[ ]\\w*'`);
  if ( $num_non_matching_seq > 0 ) {
    $self->throw( "There are " . $num_non_matching_seq .
                  " sequences with one identifier only in " . $check_dump_path . "/fastanrdb.out\nfastanrdb check FAILED!" );
  }
  else {
    say "There is no one-identifier-only header in " . $check_dump_path .
      "/fastanrdb.out\nThe dumped top level sequence and the downloaded top " . "level sequence are equivalent.\nfastanrdb check passed!";
  }

} ## end sub check_toplevel

1;
