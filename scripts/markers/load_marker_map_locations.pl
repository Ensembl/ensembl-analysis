#!/usr/local/ensembl/bin/perl -w

=head1 NAME

ensembl-analysis/scripts/markers/

=head1 SYNOPSIS

loads marker map locations into marker_map_location table

=head1 DESCRIPTION

this file loads the marker maps into the appropriate tables by parsing files
generally from UniSTS the file format is expected to be

numeric id  name chromosome position 

with any other lines marked by a #

=head1 OPTIONS

    -dbhost    host name for database (gets put as host= in locator)

    -dbport    For RDBs, what port to connect to (port= in locator)

    -dbname    For RDBs, what name to connect to (dbname= in locator)

    -dbuser    For RDBs, what username to connect as (user= in locator)

    -dbpass    For RDBs, what password to use (pass= in locator)

    -help      prints out the perl docs
=head1 EXAMPLES

perl load_marker_map_locations.pl -dbhost myhost -dbuser myuser -dbpass 
 mypass -dbname mydatabase -dbport 3306

=cut


use strict;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

my $host;
my $port;
my $dbname;
my $user;
my $pass;
my $help;

&GetOptions( 
            'dbhost:s'      => \$host,
            'dbport:n'      => \$port,
            'dbname:s'      => \$dbname,
            'dbuser:s'      => \$user,
            'dbpass:s'      => \$pass,
            'help!' => \$help,
	     ) or ($help = 1);


if ($help) {
    exec('perldoc', $0);
}

if(!$host || !$dbname || !$dbuser){
  throw("Need -dbhost $host -dbuser $dbuser and -dbname $dbname to run ".
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

my $map_name = shift;


my $map_id = &get_map_id($map_name);
die "could not get id for $map_name" if not $map_id;


while (<>) {
  /\#/ and next;

  my @l = split /\t/;

  if (not $l[0]) {
    print "no UniSTS id for $l[1]; skipping\n";
    next;
  }
  
  my $en =  {
             id => $l[0],
             name => $l[1],
             chr => $l[2],
             pos => $l[3],
            };
  
  push @{$data{$l[0]}}, $en;
}


my $marker_id_sql = "select marker_id from marker_synonym where name = ?";
my $marker_id_sth = $db->dbc->prepare($marker_id_sql);
my $syn_check_sql = "select marker_synonym_id from marker_synonym where ".
  "marker_id = ? and name = ?";
my $syn_check_sth = $db->dbc->prepare($syn_check_sql);
my $syn_insert_sql = "insert into marker_synonym(marker_id, source, name) values(?,?,?)";
my $syn_insert_sth = $db->dbc->prepare($syn_insert_sql);
my $map_loc_sql = "insert into marker_map_location(marker_id, map_id, chromosome_name, marker_synonym_id, position ) values (?,?,?,?,?)";
my $map_loc_sth = $db->dbc->prepare($map_loc_sql);
my $synonym_sql = "select marker_synonym_id from marker_synonym where marker_id = ?";
my $syn_sth = $db->dbc->prepare($synonym_sql);
foreach my $id (keys %data) {
  my $en = $data{$id};

  my ($name, $chr, $pos) = ($en->[0]->{name}, $en->[0]->{chr}, 
                               $en->[0]->{pos});
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
  $map_loc_sth->execute($marker_id, $map_id, $chr, $marker_synonym_id, $pos);
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



