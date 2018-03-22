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

ensembl-analysis/scripts/markers/load_markers.pl

=head1 SYNOPSIS

this script loads markers in the appropriate tables

=head1 DESCRIPTION

this script loads the marker and the marker_synonym tables with data from
a species markers. The file format expected is that from UniSTS and should
go like this

display_id  left primer  right primer  distance  name - accession  species

note this leads to two possibly three synonyms being stored. The first
is the display id which should be numeric and will be assigned the source
source_NUM to differentiate from the symbols

=head1 OPTIONS

    -dbhost    host name for database (gets put as host= in locator)

    -dbport    For RDBs, what port to connect to (port= in locator)

    -dbname    For RDBs, what name to connect to (dbname= in locator)

    -dbuser    For RDBs, what username to connect as (user= in locator)

    -dbpass    For RDBs, what password to use (pass= in locator)

    -file      name of file to parse

    -source    name of database markers have come from

    -write     write data into the marker and marker_synonym tables

    -help      prints out the perl docs

=head1 EXAMPLES

perl load_markers.pl -dbhost myhost -dbuser myuser -dbpass 
 mypass -dbname mydatabase -dbport 3306 -file my/marker/file -source UniSTS


=cut

use warnings ;
use strict;
use Getopt::Long qw(:config no_ignore_case);

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

my $host;
my $port=3306;
my $name;
my $user;
my $pass;
my $file;
my $source;
my $help;
my $write = 1;

GetOptions(
            'dbhost|host|h:s'      => \$host,
            'dbport|port|P:n'      => \$port,
            'dbname|db|D:s'      => \$name,
            'dbuser|user|u:s'      => \$user,
            'dbpass|pass|p:s'      => \$pass,
            'file:s'        => \$file,
            'source:s'      => \$source,
            'help!'         => \$help,
            'write!'        => \$write,
            ) or ($help = 1);

if ($help) {
  exec('perldoc', $0);
}

if(!$host || !$name || !$user){
  throw("Need -dbhost $host -dbuser $user and -dbname $name to run ".
        " use -help for docs");
}


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					    '-host'   => $host,
					    '-user'   => $user,
					    '-dbname' => $name,
					    '-pass'   => $pass,
					    '-port'   => $port,
                                            );

my $query          = "insert into marker(left_primer, right_primer, min_primer_dist, ".
    "max_primer_dist, priority, type) values (?,?,?,?,1,'est')";
my $sth            = $db->dbc->prepare($query);
my $syn_sql        = "insert into marker_synonym (marker_id, source, name) values(?, ?, ?)";
my $syn_sth        = $db->dbc->prepare($syn_sql);
my $disp_query     = "update marker set display_marker_synonym_id=? where marker_id=?";
my $disp_sth       = $db->dbc->prepare($disp_query);
my $display_id_sql = "select marker_synonym_id from marker_synonym ".
    "where name = ? and source = ? and marker_id = ?";
my $disp_id_sth    = $db->dbc->prepare($display_id_sql);

open FP, $file or throw("Failed opening".$file);

#check length of names to avoid problems
my $maxlen      = 0;
my $fieldlength = 0;
my $lenght_sth  = $db->dbc->prepare("describe marker_synonym name;");
$lenght_sth->execute();
my @row = $lenght_sth->fetchrow_array();
$fieldlength = $row[1];
$fieldlength =~ /\w+\(([\d]+)\)/;
$fieldlength = $1 or warn "\nCould not verify lenght of name field in table marker_synonym".
    "\nplease check manually!\n\n";

if($fieldlength){
  while (<FP>) {
    my @line = split /\t/, $_;
    if(scalar(@line) != 8){
      print "Somethig wrong: ",join (" ",@line),"\n"; 
    }
    my ($display_id, $lprim, $rprim, $dist, $name,
        $junk, $acc, $species) = @line;
    if($name){
      my @names = split /\;/, $name;
      foreach my $name(@names){
        if(length($name) > $maxlen){
          $maxlen = length($name);
        }
      }
    }
  }
  if($maxlen > $fieldlength){
    print STDERR "\nERROR: Length of name field in table marker_synonym has length ".
        $fieldlength.", longest entry in file has length ".$maxlen."\n".
	  "Suggested correction:\n".
	  "mysql -u$user -p$pass -h$host -P$port -D$name -e".
	  "'alter table marker_synonym modify name varchar(".($maxlen+1).");'\n\n";
    exit 0;
  }
}


seek(FP, 0, 0);

MARKER: while (<FP>) {
  my ($display_id, $lprim, $rprim, $dist, $name, 
      $junk, $acc, $species) = split /\t/, $_;
  my %names;
  my @dists = split /-/, $dist;
  #insert marker
  if (scalar(@dists)>1){
    $sth->execute($lprim, $rprim, $dists[0], $dists[1]) if($write);
  }else{
    $sth->execute($lprim, $rprim, $dist, $dist) if($write);
  }
  my $mid = $sth->{mysql_insertid};
  my $dis_source = $source."_NUM";
  #insert marker_synonym
  $syn_sth->execute($mid, $dis_source, $display_id) if($write);
  my @names = split /\;/, $name if($name);
  my $display_name;
  my $check_source;
  foreach my $id(@names){
    if($id && $id ne '-'){
      $display_name = $id if(!$display_name);
      #print "Display name is ".$display_name."\n";
      $names{$id} = 1;
    }
  }
  $check_source = $source;
  my @accs = split  /\;/, $acc if($acc);
  foreach my $id(@accs){
    if($id && $id ne '-'){
      $display_name = $id if(!$display_name);
      #print "Display name is ".$display_name."\n";
      $names{$id} = 1;
    }
  }
  foreach my $id(keys(%names)){
    $syn_sth->execute($mid, $source, $id) if($write);
  }
  if(!$display_name){
    $display_name = $display_id;
    $check_source =  $source."_NUM";
  }
  #$disp_id_sth->execute($display_name, $check_source) if($write);
  $disp_id_sth->execute($display_name, $check_source, $mid) if($write);
  my ($dis_syn_id) = $disp_id_sth->fetchrow;
  if(!$dis_syn_id && $write){
    print "SQL = ".$display_id_sql." ".$display_name." ".$source." didnt ".
        "get an id\n";
    exit;
  }
#  if($dis_syn_id == $mid && $write){
#    print "Display synonym id and marker id are both the same:  $dis_syn_id:$mid\n";
#    exit;
#  }
  $disp_sth->execute($dis_syn_id, $mid) if($write);
  
}

close(FP) or throw("Failed to close ".$file);

