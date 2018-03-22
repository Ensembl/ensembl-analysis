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

  import_transcript_support_levels

=head1 DESCRIPTION

  import_transcript_support_levels reads in a two-columns file and uses it to attach transcript attibrutes 
  
=head1 OPTIONS

  -host/dbhost     host name for database (gets put as host= in locator)
  -port/dbport     For RDBs, what port to connect to (port= in locator)
  -dbname          For RDBs, what name to connect to (dbname= in locator)
  -user/dbuser     For RDBs, what username to connect as (dbuser= in locator)
  -pass/dbpass     For RDBs, what password to use (dbpass= in locator)

  -dnahost/dnadbhost     host name for dna database (gets put as host= in locator)
  -dnaport/dnadbport     For RDBs, what port to connect to (port= in locator)
  -dnadbname             For RDBs, what name to connect to (dbname= in locator)
  -dnauser/dnadbuser     For RDBs, what username to connect as (dbuser= in locator)

  -path            Name of the assembly
  -file            Input file to read
  -verbose         Verbose options adds print statements 

=cut

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long qw(:config no_ignore_case);
use Carp;


# this ewill help when debgugging:
$| = 1;

# # #
# Set up some variables...
# # #
my $host   = '';
my $user   = 'ensro';
my $pass   = '';
my $port   = 3306;
my $dbname = '';
my $dnahost   = '';
my $dnauser   = 'ensro';
my $dnaport   = 3306;
my $dnadbname;

my $coord_system_name = 'toplevel';
my $coord_system_version;
my $code = 'TSL';

my $write; # boolean
my $verbose; # boolean
my $file;


# # #
# These are the options that can be used on the commandline
# # #
GetOptions(
  'host|h|dbhost:s'        => \$host,
  'user|u|dbuser:s'        => \$user,
  'dbname|db|D:s'          => \$dbname,
  'pass|dbpass|p:s'        => \$pass,
  'port|dbport|P:n'        => \$port,

  'dnahost|dnadbhost:s'    => \$dnahost,
  'dnauser|dnadbuser:s'    => \$dnauser,
  'dnadbname:s'            => \$dnadbname,
  'dnaport|dnadbport:n'    => \$dnaport,

  'path|cs_version:s'      => \$coord_system_version,
  'file:s'                 => \$file,
  'write'                  => \$write,
  'verbose'                => \$verbose,
);


if (!$file) {
  throw("Please define -file");
  } else {
  print STDERR "Reading from file $file\n\n";
}


# # #
# Connect to databases
# # #
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host   => $host,
  -user   => $user,
  -port   => $port,
  -pass   => $pass,
  -dbname => $dbname
);
if ($dnadbname) {
  my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => $dnahost,
                                                 -user => $dnauser,
                                                 -port => $dnaport,
                                                 -dbname => $dnadbname,
                                                 );
  $db->dnadb($dnadb);
}
my $sa  = $db->get_SliceAdaptor();
my $aa  = $db->get_AttributeAdaptor();


# delete old attribs
if ( $write ) {
  print STDERR " Deleting old attributes...\n" if $verbose;
  delete_old_attrib($db, $code);
}


# read in file containing new attributes
my %support_levels;
open(INFILE, "<$file") or die ("Can't read $file $! \n");
while(my $line = <INFILE>){
  chomp $line;
  my @fields = split(/\t/, $line);
  my ($transcript_id, $version) = split(/\./, $fields[0]);
  my $support_level = $fields[1];
  $support_levels{$transcript_id}{'version'} = $version;
  $support_levels{$transcript_id}{'support_level'} = $support_level;
}
close INFILE;
print STDERR "Fetched ".(scalar(keys %support_levels))." new attributes\n";


# # #
# Fetch the sequences we are interested in - all or subset
# # #
my @slices = @{ $sa->fetch_all( $coord_system_name, $coord_system_version, 1, undef ) };
print STDERR "Got ".( scalar( @slices) )." slices\n" if $verbose;



# # #
# Now loop through each slices
# and then each gene on the slice
# # #
my $stable_id_in_file = 0;
foreach my $slice ( @slices ) {
  print STDERR "Doing slice ".$slice->seq_region_name."\n" if $verbose;
  my $gene_cnt = 0;
  my $transc_cnt = 0;
  my $transc_uptodate = 0;
  my $transc_no_data = 0;
  my $transc_updated = 0;

  # now look for new candidates
  foreach my $gene ( @{$slice->get_all_Genes} ) {
    #print STDERR "Gene ".$gene->stable_id."\n";
    $gene_cnt++;
    foreach my $transcript ( @{$gene->get_all_Transcripts} ) {
      $transc_cnt++;
      
      if (exists $support_levels{$transcript->stable_id}) {
        # oh good, found this stable ID
        if ( $transcript->version == $support_levels{$transcript->stable_id}{'version'} ) {
          # up to date
          $transc_uptodate++;
          $stable_id_in_file++;
          if ($write) {
            store_attrib($aa, $transcript, $code, $support_levels{$transcript->stable_id}{'support_level'});
            print STDERR "  writing ".$transcript->stable_id." version ".$transcript->version." TSL ".$support_levels{$transcript->stable_id}{'support_level'}."\n";
          }
        } else {
          # annotation likely changed since last release
          $transc_updated++;
          $stable_id_in_file++;
          warning("Transcript annotation mismatch ".$transcript->stable_id. " version in Ensembl=".$transcript->version." vs version in file=".$support_levels{$transcript->stable_id}{'version'});
          if ($write) {
            store_attrib($aa, $transcript, $code, $support_levels{$transcript->stable_id}{'support_level'}." (assigned to previous version ".$support_levels{$transcript->stable_id}{'version'}.")");
            print STDERR "  writing ".$transcript->stable_id." version ".$transcript->version." TSL ".$support_levels{$transcript->stable_id}{'support_level'}."\n";
          }
        }
      } else {
        # this is likely a new transcript that wasn't annotated last release
        $transc_no_data++;
        warning("No data in file for ".$transcript->stable_id);
      }
    }
  }
  print  "Slice ".$slice->seq_region_name." has genes $gene_cnt with $transc_cnt transcripts. There are $transc_uptodate with current attributes, $transc_updated transcripts with updated annotation and $transc_no_data with no attributes\n";
}
print "Matched stable_ids for ".$stable_id_in_file." of ".(scalar(keys %support_levels))." transcripts in file\n";
print STDERR "DONE!\n\nNow grep for:\n^Slice and ^Matched\n\n";



sub delete_old_attrib {
  my ($db, $code) = @_;

  my $sql = q{
    DELETE ta
    FROM transcript_attrib ta, attrib_type att 
    WHERE att.attrib_type_id = ta.attrib_type_id 
    AND att.code = ? };

  my $sth = $db->dbc->prepare($sql);
  $sth->bind_param( 1, $code);
  $sth->execute();
  $sth->finish();

  return;
}

sub store_attrib {
  my ($aa, $transcript, $code, $value) = @_;

  my ($attrib_type_id, $newcode, $name, $description ) = $aa->fetch_by_code($code);
  if (!$attrib_type_id) {
    throw("Unable to fetch attrib_type with code $code");
  }
  my $attrib = Bio::EnsEMBL::Attribute->new(
    -NAME        => $name,
    -CODE        => $code,
    -VALUE       => $value,
    -DESCRIPTION => $description
  );
  $aa->store_on_Transcript($transcript, [$attrib]);

  return;
}





