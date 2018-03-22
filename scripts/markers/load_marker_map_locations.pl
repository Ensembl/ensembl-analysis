#!/usr/bin/env perl

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

=head1 NAME

ensembl-analysis/scripts/markers/

=head1 SYNOPSIS

loads marker map locations into marker_map_location table

=head1 DESCRIPTION

this file loads the marker maps into the appropriate tables by parsing files
generally from UniSTS the file format is expected to be

numeric id  name chromosome position 

with any other lines (header info) started by a #

=head1 OPTIONS

    -dbhost    host name for database (gets put as host= in locator)

    -dbport    For RDBs, what port to connect to (port= in locator)

    -dbname    For RDBs, what name to connect to (dbname= in locator)

    -dbuser    For RDBs, what username to connect as (user= in locator)

    -dbpass    For RDBs, what password to use (pass= in locator)
 
    -map_file  File containing map locations

    -map_name  Name of map

    -create    Insert the marker set name in the map table
 
    -write     writes data into the marker_map_location table
 
    -help      prints out the perl docs

=head1 EXAMPLES

perl load_marker_map_locations.pl -dbhost myhost -dbuser myuser -dbpass 
 mypass -dbname mydatabase -dbport 3306 -map_file marker_map_location_file.txt
-map_name MAPNAME

=cut

use warnings ;
use strict;
use Getopt::Long qw(:config no_ignore_case);

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

my $host;
my $port=3306;
my $dbname;
my $user;
my $pass;
my $map_name;
my $map_file;
my $help;
my $write;
my $create = 0;

GetOptions( 
            'dbhost|host|h:s'      => \$host,
            'dbport|port|P:n'      => \$port,
            'dbname|db|D:s'      => \$dbname,
            'dbuser|user|u:s'      => \$user,
            'dbpass|pass|p:s'      => \$pass,
            'map_name:s'    => \$map_name,
            'map_file:s'    => \$map_file,
            'write'         => \$write,
            'create!'         => \$create,
            'help!' => \$help,
	     ) or ($help = 1);


if ($help) {
  exec('perldoc', $0);
}

if (!$host || !$dbname || !$user) {
  throw("Need -dbhost $host -dbuser $user and -dbname $dbname to run ".
        " use -help for docs");
}


my %data;

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $host,
					    '-user'   => $user,
					    '-dbname' => $dbname,
					    '-pass'   => $pass,
					    '-port'   => $port,
					   );

if ($create) {
    my $create_map_cmd = 'INSERT INTO map (map_name) VALUES ("'.$map_name.'")';
    $db->dbc->do($create_map_cmd);
}
my $map_id = &get_map_id($map_name);
throw("could not get id for $map_name") if not $map_id;

open(FH, $map_file) or throw("Couldn't open ".$map_file);

while (<FH>) {
  /\#/ and next;

  #print $_,"\n";
  my @l = split /\t/;

  if (not $l[0]) {
    print "no UniSTS id for $l[1]; skipping\n";
    next;
  }
  
  my $lod;

  if ($l[4] && $l[4] =~ /^P?\>?\d+\.?\d*/){
    my @lod_val = $l[4]=~ /^P?\>?(\d+\.?\d*)/;
    $lod = $lod_val[0];
   # print "LOD  ",$lod_val[0],"\n";
  }else{
    $lod = "NULL";
   # print "LOD  IS NULL",,"\n";
  }

  my $en =  {
             id   => $l[0],
             name => $l[1],
             chr  => $l[2],
             pos  => $l[3],
             lod  => $lod,
            };
  
  push @{$data{$l[0]}}, $en;
}

close(FH) or throw("Couldn't close ".$map_file);

my $marker_id_sql = "select marker_id from marker_synonym where source = \"UniSTS_NUM\" and name = ?";
my $marker_id_sth = $db->dbc->prepare($marker_id_sql);

my $syn_check_sql = "select marker_synonym_id from marker_synonym where marker_id = ? and name = ?";
my $syn_check_sth = $db->dbc->prepare($syn_check_sql);

my $syn_insert_sql = "insert into marker_synonym(marker_id, source, name) values(?,?,?)";
my $syn_insert_sth = $db->dbc->prepare($syn_insert_sql);

my $check_map_loc_sql = "select marker_id from marker_map_location where marker_id = ? and  map_id = ? and chromosome_name = ? and marker_synonym_id = ? and  position = ? and lod_score = ?";
my $check_map_loc_sth = $db->dbc->prepare($check_map_loc_sql);


my $map_loc_sql = "insert into marker_map_location(marker_id, map_id, chromosome_name, marker_synonym_id, position, lod_score ) values (?,?,?,?,?,?)";
my $map_loc_sth = $db->dbc->prepare($map_loc_sql);

my $map_loc_sql_null = "insert into marker_map_location(marker_id, map_id, chromosome_name, marker_synonym_id, position ) values (?,?,?,?,?)";
my $map_loc_sth_null = $db->dbc->prepare($map_loc_sql_null);


my $synonym_sql = "select marker_synonym_id from marker_synonym where marker_id = ?";
my $syn_sth     = $db->dbc->prepare($synonym_sql);

if($write){
  foreach my $id (keys %data) {
    my $en = $data{$id};
    
    my ($name, $chr, $pos, $lod) = ($en->[0]->{name}, $en->[0]->{chr}, 
                              $en->[0]->{pos}, $en->[0]->{lod});
    $marker_id_sth->execute($id);
    my ($marker_id) = $marker_id_sth->fetchrow;
    if(!$marker_id){
      throw("Marker ".$id." ".$name." has failed to produce an internal id ".
            "cannot continue");
    }
    $syn_check_sth->execute($marker_id, $name);
    my ($marker_synonym_id) = $syn_check_sth->fetchrow;
    if(!$marker_synonym_id){
      print "No marker synonym found for ".$name." inserting one\n";
      $syn_insert_sth->execute($marker_id, $map_name, $name);
      $marker_synonym_id = $syn_insert_sth->{mysql_insertid};
    }
    
    $check_map_loc_sth->execute($marker_id, $map_id, $chr, $marker_synonym_id, $pos, $lod);
    my ($entry_check) = $check_map_loc_sth->fetchrow;
    
    if (!$entry_check){
      if ($lod eq "NULL"){
        $map_loc_sth_null->execute($marker_id, $map_id, $chr, $marker_synonym_id, $pos);
      }else{
        #print "ENTRY: ",$marker_id,"  ", $map_id,"  ", $chr,"  ", $marker_synonym_id,"  ", $pos,"  ", $lod,"\n";
        $map_loc_sth->execute($marker_id, $map_id, $chr, $marker_synonym_id, $pos, $lod);
      }
    }else{
      print "Entry already exists ",$marker_id,"  ", $map_id,"  ", $chr,"  ", $marker_synonym_id,"  ", $pos,"  ", $lod,"\n";
    }
  }
}

sub get_map_id {
  my $mname = shift;
  
  my $map_query = "select map_id from map where map_name = '$mname'";
  my $map_sth = $db->dbc->prepare($map_query);
  
  $map_sth->execute;
  
  my $map_ref = $map_sth->fetchrow_hashref;

  if ($map_ref) {
    return $map_ref->{map_id};
  } else {
    return 0;
  }
  $map_sth->finish;
}



