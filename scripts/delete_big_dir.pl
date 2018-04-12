#!/usr/bin/env perl

# Copyright [2018] EMBL-European Bioinformatics Institute
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

=head1 NAME

  delete_big_dir.pl

=head1 DESCRIPTION

This script take a single argument, the full path to a dir to remove and then removes all
files and subdirs before finally removing the dir itself. This is designed to remove very
complex subdir structures or dirs with a very large number of files in them. Perl is much
faster at this kind of task than rm or rsync to an empty dir

=cut

use Cwd;
use File::Spec;
use warnings;
use strict;
use feature 'say';
use Getopt::Long;

my $full_master_dir_path;
GetOptions('dir:s' => \$full_master_dir_path);

unless($full_master_dir_path) {
  die "No agruments entered. You need to pass in the name of the dir in the current directory to delete";
}

$full_master_dir_path = File::Spec->rel2abs($full_master_dir_path);

unless(-d $full_master_dir_path) {
  die "The argument you entered is not a dir. Argument entered: ".$full_master_dir_path;
}

say "The full path for the dir to be deleted is:\n".$full_master_dir_path;

say "Getting subdir list...";
my @subdirs = `lfs find $full_master_dir_path -type d`;

say "Found ".(scalar(@subdirs) - 1)." subdirs";

for (my $i=5; $i>0; $i--) {
  say "Beginning file deletion in ".$i."...";
  sleep(1);
}

print "\n";

foreach my $dir (@subdirs) {
  chomp $dir;
  say "Removing files from:\n".$dir;
  foreach my $file (<$dir/*>) {
    unless($file =~ /^$full_master_dir_path/) {
      die "Potential issue with file path, path didn't match to the master dir path. Path found:\n".$file;
    }
    unlink($file);
  }
}

say "\nFinished removing files. Now removing empty dirs...";
my $result = system('rm -r '.$full_master_dir_path);
if($result) {
  die "Could not remove the master dir, something potentially went wrong with the deletion!";
}

exit;

