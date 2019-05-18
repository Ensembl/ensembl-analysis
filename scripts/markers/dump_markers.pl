#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

# Dumps markers from the marker table into a file, in the format expected by
# the EPCR Runnable

=head1 NAME

ensembl-analysis/scripts/markers/dump_markers.pl

=head1 SYNOPSIS

a script to dump the marker table out into a format for epcr

=head1 DESCRIPTION

This script takes the entries from the marker table and dumps them out
in the format expected by the runnable Bio::EnsEMBL::Analysis::Runnable::EPCR
The file should be be pointed at by the dbfile column of the analysis table
for the analysis which is to run the EPCR

=head1 OPTIONS
    -dbhost    host name for database (gets put as host= in locator)

    -dbport    For RDBs, what port to connect to (port= in locator)

    -dbname    For RDBs, what name to connect to (dbname= in locator)

    -dbuser    For RDBs, what username to connect as (user= in locator)

    -dbpass    For RDBs, what password to use (pass= in locator)

    -outfile   the path to the file to dump the data into

    -help      prints out the perl docs

=head1 EXAMPLES

perl dump_markers -dbhost myhost -dbuser myuser -dbpass mypass -dbname
  mydatabase -dbport 3306 -outfile /path/to/epcr.marker.file

=cut

use warnings ;
use strict;
use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

my $host   = '';
my $user   = 'ensro';
my $pass   = undef;
my $dbname = '';
my $port = 3306;
my $outfile;
my $help;


$| = 1;

GetOptions(
            'dbhost|host|h:s'   => \$host,
            'dbuser|user|u:s'   => \$user,
            'dbname|db|D:s' => \$dbname,
            'dbport|port|P:n'   => \$port,
            'dbpass|pass|p:s'   => \$pass,
            'outfile:s' => \$outfile,
            'help!' => \$help,
           ) or($help = 1);

if ($help) {
    exec('perldoc', $0);
}

if(!$host || !$dbname || !$user){
  throw("Need -dbhost $host -dbuser $user and -dbname $dbname to run ".
        " use -help for docs");
}



# Open database
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host   => $host,
  -user   => $user,
  -pass   => $pass,
  -port   => $port,
  -dbname => $dbname,
);


my $sts = $db->get_MarkerAdaptor->fetch_all;
die "No markers in database" unless @{$sts};
#die scalar(@{$sts}) . " markers in database";


dump_sts_file($outfile, $sts);

sub dump_sts_file {
  my ($dest, $sts) = @_;
  
  open DEST, "> $dest";
  foreach my $m (@{$sts}) {
    unless (ref $m && $m->isa("Bio::EnsEMBL::Map::Marker")) {
      die "Object not a Bio::EnsEMBL::Map::Marker: [$m]";
    }
    
    next unless length($m->left_primer) > 0;
    next unless length($m->right_primer) > 0;
    
    my ($min_dist, $max_dist);
    if ($m->min_primer_dist == 0) {
      $min_dist = 80;
    } else {
      $min_dist = $m->min_primer_dist;
    }
    
    if ($m->max_primer_dist == 0) {
      $max_dist = 600;
    } else {
      $max_dist = $m->max_primer_dist;
    }
    
    my $dist; 
    if ($min_dist == $max_dist) {
      $dist = $min_dist;
    } else {
      $dist = join("-", $min_dist, $max_dist);
    }
    
    print DEST join("\t",
                    $m->dbID,
                    $m->left_primer,
                    $m->right_primer,
                    $dist,
                    ), "\n";
  }
  close DEST;
}

