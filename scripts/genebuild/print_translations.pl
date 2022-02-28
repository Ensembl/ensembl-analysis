#!/usr/bin/env perl

# Copyright [2017-2022] EMBL-European Bioinformatics Institute
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

 print_translations.pl

=head1 DESCRIPTION

  A script to dump translations in fasta format (with the stable id as a header) based on a list of input ids.
  The script will check each id to see if it's an ENSP, ENST or CCDS stable id and act appropriately. If the
  it is none of the above it will assume the id is a translation db id. The script will write to STD out by
  default, but if an output file is provided then it will write to that

=head1 OPTIONS
  -host
  -port
  -user
  -dbname
  -dnahost
  -dnaport
  -dnauser
  -dnadbname
  -id_file
  -output_file

=head1 EXAMPLE

  perl print_translations.pl \
      --host=host --user=ensro \
      --dbname=core_db \
      --port=XXX
      --id_file=translation_ids_to_dump
      --output_file=outfile.fa

=cut

use warnings;
use strict;
use feature 'say';
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long qw(:config no_ignore_case);

my $dbname = '';
my $user   = 'ensro';
my $host   = '';
my $port;
my $pass;
my $id_file = '';
my $output_file = '';

my $dnadbname = '';
my $dnauser   = 'ensro';
my $dnahost   = '';
my $dnaport;
my $dnapass;

my $result = GetOptions ("user|dbuser|u=s"      => \$user,
                         "host|dbhost|h=s"      => \$host,
                         "port|dbport|P=i"      => \$port,
                         "dbname|db|D=s"    => \$dbname,
                         "dbpass|pass|p=s"  => \$pass,
                         "dnauser=s"   => \$dnauser,
                         "dnahost=s"   => \$dnahost,
                         "dnaport=i"   => \$dnaport,
                         "dnadbname=s" => \$dnadbname,
                         "dnadbpass=s"  => \$dnapass,
                         "id_file=s"   => \$id_file,
                         "output_file=s" => \$output_file);

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $port,
  -user    => $user,
  -host    => $host,
  -dbname  => $dbname,
  -pass    => $pass);

my $dnadb;
if ($dnadbname && $dnauser  && $dnahost && $dnaport) {
  $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $dnaport,
  -user    => $dnauser,
  -host    => $dnahost,
  -dbname  => $dnadbname,
  -pass   => $dnapass);
  $db->dnadb($dnadb);
}

my $translation_adaptor = $db->get_TranslationAdaptor();

unless(-e $id_file) {
  throw("Input file does not exist. Path: ".$id_file);
}

my @translation_ids = ();
open(IN,$id_file);
@translation_ids = <IN>;
close IN;

if($output_file) {
  open(OUT,">".$output_file);
}

foreach my $translation_id (@translation_ids) {
  chomp $translation_id;
  my $translation;
  if ( $translation_id =~ /^ENS.*P\d+/){
    $translation = $translation_adaptor->fetch_by_stable_id($translation_id); 
  } elsif ( $translation_id =~ /^ENS.*T\d+/){
    my $transcript = $db->get_TranscriptAdaptor->fetch_by_stable_id($translation_id);
    $translation = $transcript->translation;
  } elsif ( $translation_id =~ /^CCDS\d+\.\d+/){
    my $transcript = $db->get_TranscriptAdaptor->fetch_by_stable_id($translation_id);
    $translation = $transcript->translation;
  } else {
    $translation = $translation_adaptor->fetch_by_dbID($translation_id);
  }

  if (!$translation && $translation_id !~ /\.\d+/) {
    # try to add a version on
    $translation_id .= ".1";
    $translation = $translation_adaptor->fetch_by_stable_id($translation_id);
  }

  if (!$translation) {
    throw("Could not find translation for $translation_id");
  }

  # format into 60-base lines
  my $translation_seq = $translation->seq;
  $translation_seq =~ s/(\w{60})/$1\n/g;

  if($output_file) {
    say OUT ">".$translation->stable_id().".".$translation->version()."\n".$translation_seq;
  } else {
    say ">".$translation->stable_id().".".$translation->version()."\n".$translation_seq;
  }
}

if($output_file) {
  close OUT;
}

exit;
