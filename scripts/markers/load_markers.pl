#!/usr/local/ensembl/bin/perl -w

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

    -file   name of file to parse

    -source name of database markers have come from

    -help      prints out the perl docs

=head1 EXAMPLES

perl load_markers.pl -dbhost myhost -dbuser myuser -dbpass 
 mypass -dbname mydatabase -dbport 3306 -file my/marker/file -source UniSTS


=cut



use strict;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

my $host;
my $port;
my $name;
my $user;
my $pass;
my $file;
my $source;
my $help;
my $write = 1;

&GetOptions( 
            'dbhost:s'      => \$host,
            'dbport:n'      => \$port,
            'dbname:s'      => \$name,
            'dbuser:s'      => \$user,
            'dbpass:s'      => \$pass,
            'file:s'        => \$file,
            'source:s'      => \$source,
            'help!' => \$help,
            'write!' => \$write,
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

my $query = "insert into marker(left_primer, right_primer, min_primer_dist, max_primer_dist, priority, type) values (?,?,?,?,1,'est')";
my $sth = $db->dbc->prepare($query);
my $syn_sql = "insert into marker_synonym (marker_id, source, name) ".
  "values(?, ?, ?)";
my $syn_sth = $db->dbc->prepare($syn_sql);
my $disp_query = "update marker set display_marker_synonym_id=? where marker_id=?";
my $disp_sth = $db->dbc->prepare($disp_query);
my $display_id_sql = "select marker_synonym_id from marker_synonym ".
  "where name = ? and source = ?";
my $disp_id_sth = $db->dbc->prepare($display_id_sql);
  
open FP, $file or throw("Failed opening".$file);

MARKER:

while (<FP>) {
  my ($display_id, $lprim, $rprim, $dist, $name, 
      $junk, $acc, $species) = split /\t/, $_;
  my %names;
  $sth->execute($lprim, $rprim, $dist, $dist) if($write);
  my $mid = $sth->{mysql_insertid};
  my $dis_source = $source."_NUM";
  $syn_sth->execute($mid, $dis_source, $display_id) if($write);
  my @names = split /\;/, $name if($name);
  my $display_name;
  foreach my $id(@names){
    if($id && $id ne '-'){
      $display_name = $id if(!$display_name);
      #print "Display name is ".$display_name."\n";
      $names{$id} = 1;
    }
  } 
  my $check_source = $source;
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
  $disp_id_sth->execute($display_name, $check_source) if($write);
  my ($dis_syn_id) = $disp_id_sth->fetchrow;
  if(!$dis_syn_id && $write){
    print "SQL = ".$display_id_sql." ".$display_name." ".$source." didnt ".
      "get an id\n";
    exit;
  }
  if($dis_syn_id == $mid && $write){
    print "Dis syn id and mid seem be both be $dis_syn_id:$mid\n";
    exit;
  }
  $disp_sth->execute($dis_syn_id, $mid) if($write);
  
}


close(FP) or throw("Failed to close ".$file);
