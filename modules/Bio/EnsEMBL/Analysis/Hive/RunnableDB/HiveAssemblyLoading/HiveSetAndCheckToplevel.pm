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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveSetAndCheckToplevel;

use strict;
use warnings;
use feature 'say';
use File::Spec::Functions;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;

  unless($self->param('target_db')) {
    $self->throw("target_db flag not passed into parameters hash. The core db to load the assembly info ".
                 "into must be passed in with write access");
  }

  unless($self->param('enscode_root_dir')) {
    $self->throw("enscode_root_dir flag not passed into parameters hash. You need to specify where your code checkout is");
  }

  return 1;
}

sub run {
  my $self = shift;

  say "Setting toplevel flag and checking";
  my $target_db = $self->param('target_db');
  my $enscode_dir = $self->param('enscode_root_dir');
  my $primary_assembly_dir_name = $self->param('primary_assembly_dir_name');
  my $path_to_files = catdir($self->param('output_path'), $primary_assembly_dir_name);

  $self->set_toplevel($target_db,$enscode_dir);
  $self->check_toplevel($target_db,$enscode_dir,$path_to_files)
    if ($self->param_is_defined('primary_assembly_dir_name'));
  say "Finished checking contig files";
  return 1;
}

sub write_output {
  my $self = shift;

  return 1;
}

sub set_toplevel {
  my ($self,$target_db,$enscode_dir) = @_;

  my $dbhost = $target_db->{'-host'};
  my $dbport = $target_db->{'-port'};
  my $dbuser = $target_db->{'-user'};
  my $dbpass = $target_db->{'-pass'};
  my $dbname = $target_db->{'-dbname'};

  my $cmd = 'perl '.catfile($enscode_dir, 'ensembl-analysis', 'scripts', 'assembly_Loading', 'set_toplevel.pl').
            " -dbhost ".$dbhost.
            " -dbuser ".$dbuser.
            " -dbpass ".$dbpass.
            " -dbport ".$dbport.
            " -dbname ".$dbname;

  my $return = system($cmd);
  if($return) {
    $self->throw("The set_toplevel script returned a non-zero exit code. Commandline used:\n".$cmd);
  }

  my $num_toplevel = int(`mysql -h$dbhost -P$dbport -u$dbuser -p$dbpass -D$dbname -NB -e'select count(*) from seq_region_attrib where attrib_type_id = 6'`);
  if ($num_toplevel > 0) {
    say "Found ".$num_toplevel." toplevel (non-redundant) seq region attributes found in seq_region_attrib table (code 6)";
  } else {
    $self->throw("No toplevel (non-redundant) seq region attributes found in seq_region_attrib table (code 6). Commandline used:\n".$cmd);
  }

}

sub check_toplevel {
  my ($self,$target_db,$enscode_dir,$path_to_files) = @_;
  my $check_dump_path = $path_to_files."/check_and_dump_toplevel";

  my $dbhost = $target_db->{'-host'};
  my $dbport = $target_db->{'-port'};
  my $dbuser = $target_db->{'-user'};
  my $dbpass = $target_db->{'-pass'};
  my $dbname = $target_db->{'-dbname'};

  if(-d $check_dump_path) {
    `rm -r $check_dump_path`;
  }

 `mkdir $check_dump_path`;
  my $cmd = "";
  my $return;
  unless(-e $path_to_files."/FASTA/") {
    $self->warning("The path to the FASTA dir does not exist. Assuming assembly is single level");
    $cmd = "cp ".$path_to_files."/contigs/contigs.fa ".$check_dump_path."/downloaded_toplevel.fa";
    $return = system($cmd);
    if($return) {
      $self->throw("Copy of contig sequences filed. Commandline used:\n".$cmd);
    }
  } else {
    $cmd = "ls ".$path_to_files."/FASTA/*.fna | grep -v '\\.placed.' | xargs cat >> ".$check_dump_path."/downloaded_toplevel.fa";
    $return = system($cmd);
    if($return) {
      $self->throw("Concatenation of downloaded fna sequences filed. Commandline used:\n".$cmd);
    }
  }

  my $num_toplevel_downloaded = int(`grep -c ">" $check_dump_path/downloaded_toplevel.fa`);
  if($num_toplevel_downloaded == 0) {
    $self->throw("Found no headers in the concatenated downloaded files. Commandline used:\n".
                 "grep -c \">\" ".$check_dump_path."/downloaded_toplevel.fa");
  }

  say "\nDumping top level sequence...";
  $cmd = "perl ".$enscode_dir."/ensembl-analysis/scripts/sequence_dump.pl".
         " -dbhost ".$dbhost.
         " -dbuser ".$dbuser.
         " -dbpass ".$dbpass.
         " -dbport ".$dbport.
         " -dbname ".$dbname.
         " -toplevel ".
         " -onefile ".
         " -output_dir ".$check_dump_path.
         " > ".$check_dump_path."/sequence_dump.out";

  $return = system($cmd);
  if($return) {
    $self->throw("Dump of toplevel sequences failed. Commandline used:\n".$cmd);
  }

  say "Finished dumping toplevel sequences. Output written to:\n".$check_dump_path."/sequence_dump.out";

  my $num_toplevel_dumped = int(`grep -c ">" $check_dump_path/toplevel.fa`);
  if($num_toplevel_dumped == 0) {
    $self->throw("Found no headers in the dumped toplevel file. Commandline used:\n".
                 "grep -c \">\" ".$check_dump_path."/toplevel.fa");
  }

  if ($num_toplevel_dumped != $num_toplevel_downloaded) {
    $self->throw("The number of dumped top level sequences (".$num_toplevel_dumped.") does not match the number of downloaded top level sequences (".$num_toplevel_downloaded.")");
  } else {
    say "The number of dumped top level sequences (".$num_toplevel_dumped.") matches the number of downloaded top level sequences (".$num_toplevel_downloaded."). Great!";
  }

  say "Combining dumped top level file and downloaded top level file into ".$check_dump_path."/combined_toplevel.fa...\n";

  $cmd = "cat ".$check_dump_path."/toplevel.fa ".$check_dump_path."/downloaded_toplevel.fa > ".$check_dump_path."/combined_toplevel.fa";
  $return = system($cmd);
  if($return) {
    $self->throw("Concatenation of dumped and downloaded toplevel sequence files failed. Commandline used:\n".$cmd);
  }

  my $fastanrdb_out = catfile($check_dump_path, 'fastanrdb.out');
  $cmd = "fastanrdb ".$check_dump_path."/combined_toplevel.fa > ".$fastanrdb_out;
  $return = system($cmd);

  if($return) {
    $self->throw("fastanrdb died unexpectedly. Commandline used:\n".$cmd);
  }

  say "The output file $fastanrdb_out was written";

  # check that the number of identifiers in headers matches the expected number of top level sequences*2
  my $num_toplevel_fastanrdb = 0;
  open(FH, $fastanrdb_out) || $self->throw("Could not open $fastanrdb_out");
  while (<FH>) {
    if (/>/) {
      my @ids = split(' ', $_);
      if (@ids == 2) {
        ++$num_toplevel_fastanrdb;
        $self->throw("Something is wrong in $fastanrdb_out: $_") unless ($ids[0] =~ /:$ids[1]:/);
      }
      elsif (scalar(@ids)%2 == 0) {
        my %uniq;
        foreach my $id (@ids) {
          $id =~ s/>?[^:]+:[^:]+:([^:]+):\S+/$1/;
          $uniq{$id}++;
        }
        if (grep {$_ != 2} values %uniq) {
          $self->throw("Something is wrong in $fastanrdb_out: $_");
        }
        else {
          $num_toplevel_fastanrdb += keys %uniq;
          $self->warning("You have duplicated sequences in the assembly $_");
        }
      }
      else {
        $self->throw("Something is wrong in $fastanrdb_out: $_");
      }
    }
  }
  close(FH) || $self->throw("Could not close $fastanrdb_out");
  if ($num_toplevel_fastanrdb != $num_toplevel_dumped) {
    $self->throw("The number of sequences in file $fastanrdb_out is $num_toplevel_fastanrdb, but it should be ".
                 $num_toplevel_dumped.". It means that we have duplicates.\n");

    # maybe add more checks to inform about the sequences that does not match
    # cat $SCRA/fastanrdb.log | perl -ne 'next unless ($a,$b) = /^>(\w+\.\w+)\s*(\w+\.\w+)/; if ($a eq $b) {print "$a\n"}' | wc -l
    # to look at cases where NCBI and Ref DB agree
    # bsub -o $SCRA/fastanrdb_minus_r.log fastanrdb temp.fa -r
    # to compare sequences in reverse complement too

  } else {
    say "The number of sequences in file ".$check_dump_path."/fastanrdb.out is ".$num_toplevel_fastanrdb.", which is good since it ".
        "is the expected number, but there might be odd cases where duplicates can be present.";
  }

}

1;
