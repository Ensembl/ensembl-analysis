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

use strict;
use warnings;
$| = 1;

use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::UnmappedObject;
use Bio::EnsEMBL::DBSQL::UnmappedObjectAdaptor;
use Bio::EnsEMBL::Map::DBSQL::MarkerAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $host = '';
my $user = 'enrso';
my $dbname = '';
my $port = '';
my $max_duplicates = 5;
my @duplicates;
my @unmapped;
my $pass = '';
my $ln = 'Marker';
my $dry_run = 0;

GetOptions(
	   'host|dbhost|h:s'      => \$host,
	   'user|dbuser|u:s'      => \$user,
	   'dbname|db|D:s'           => \$dbname,
	   'port|dbport|P:n'      => \$port,
	   'logic_name:s'       => \$ln,
	   'pass|dbpass|p:s'      => \$pass,
	   'max_duplicates:n'   => \$max_duplicates,
           'dry_run!'           => \$dry_run,
	  );

unless ($host && $user && $dbname && $port){
  die("unmapped_markers.pl - removes duplicated and unaligned markers from the database and populates the unmapped feature table
-host    $host
-user    $user
-dbname  $dbname
-port    $port
-logic_name $ln
-pass    $pass
-max_duplicates $max_duplicates
");
}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor
  (
   -host   => $host,
   -user   => $user,
   -port   => $port,
   -dbname => $dbname,
   -pass   => $pass,
  );

die("Cannot open connection to $dbname:$host.$port for  $user\n") unless $db;

if ($dry_run) {
  print STDERR "This is a dry run.\n";
}

print "You will be removing markers which map to the genome more than $max_duplicates times.\n";

my $aa = $db->get_AnalysisAdaptor;
my $analysis = $aa->fetch_by_logic_name($ln);
die("Cannnot find analysis with logic name $ln \n") unless $analysis;
my $ma = $db->get_MarkerAdaptor;
my $umma = $db->get_UnmappedObjectAdaptor;

# fetch Markers
my @markers = @{$ma->fetch_all};
die("No markers found\n") unless  scalar(@markers) > 0 ;
print STDERR "Found " . scalar(@markers) . " markers, checking / setting marker_feature map weights....\n";
foreach my $marker ( @markers ) {
  my @features = @{$marker->get_all_MarkerFeatures} ;
  my @map_locations = @{$marker->get_all_MapLocations} ;
  if ( (scalar@features > $max_duplicates) && scalar(@map_locations) ==  0 ) {
    push @duplicates,$marker;
  } elsif ( scalar(@features == 0 ) && scalar(@map_locations) ==  0 ) {
    push @unmapped,$marker;
  }

  # update map weights
  foreach my $feat (@features) {
    if ( scalar(@features) != $feat->map_weight ) {
      sql("UPDATE marker_feature SET map_weight = " . scalar(@features) . " WHERE marker_feature_id = " . $feat->dbID . ";",$db);
    }
  }
}

print STDERR "Found " . scalar(@duplicates) . " duplicates and " . scalar(@unmapped) . " unmapped markers\n";

if ($dry_run) {
  print STDERR "Dry run finished. Unmapped entries not created.  No markers or marker_features have been deleted from the DB.\n";
}

if (!$dry_run) {
  print STDERR "Creating unmapped entries...\n";

  my @unmapped_objects = @{make_entries("duplicates",\@duplicates)};
  push @unmapped_objects,@{make_entries("unmapped",\@unmapped)};
  foreach my $umo ( @unmapped_objects ) {
    $umma->store($umo);
  }
  print STDERR "Deleting duplicated and unmapped markers...\n";

  # delete the markers that dont map 
  foreach my $marker ( @duplicates ) {
    sql("DELETE from marker where marker_id = " . $marker->dbID . ";",$db);
  }
  foreach my $marker ( @unmapped ) {
    sql("DELETE from marker where marker_id = " . $marker->dbID . ";",$db);
  }

  # tidy up associated tables
  # marker feature
  sql("DELETE marker_feature from marker_feature LEFT JOIN marker
  ON marker_feature.marker_id =  marker.marker_id 
  WHERE marker.marker_id is null ",$db);
  # marker feature where the marker hits lots of times but is also in the marker_map_locations table
  # just delete the multiple hitting markers
  sql("DELETE marker_feature from marker_feature 
  WHERE map_weight > $max_duplicates",$db);
  # marker synonym
  sql("DELETE marker_synonym from marker_synonym LEFT JOIN marker
  ON marker_synonym.marker_id =  marker.marker_id 
  WHERE marker.marker_id is null",$db);

  print STDERR "Done\n";
}

exit;

sub make_entries{
  my ($type,$array_ref) = @_;
  my @umos;
  foreach my $marker ( @$array_ref ){
    # one entry per synonym
    foreach my $synonym ( @{$marker->get_all_MarkerSynonyms} ) {
      my $umo;
      # dont add synonyms that are just numbers
      next unless $synonym->name =~ /\D/;
      if ( $type eq 'duplicates' ) {
	$umo = Bio::EnsEMBL::UnmappedObject->new 
	  (
	   -type => 'Marker',
	   -analysis => $analysis,
	   -identifier => $synonym->name,
	   -summary => "Marker matches multiple times",
	   -full_desc => "Marker aligns to the genome > $max_duplicates times"
	  );
      } elsif ( $type eq 'unmapped' ) {
	$umo = Bio::EnsEMBL::UnmappedObject->new 
	  (
	   -type => 'Marker',
	   -analysis => $analysis,
	   -identifier => $synonym->name,
	   -summary => "Marker does not align",
	   -full_desc => "Unable to align to the genome"
	  );
      }
      push @umos,$umo;
    }
  }
  return \@umos;
}

sub sql {
  my ($query,$db) = @_;
  my $sth = $db->dbc->prepare($query);
  $sth->execute();
  return;
}
