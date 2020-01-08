#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

  import_transcript_appris_results

=head1 DESCRIPTION

  import_transcript_appris_results reads in a two-columns file and uses it to attach transcript attibrutes 
  
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
use DateTime;
use Carp;
use feature 'say';

# this ewill help when debgugging:
$| = 1;

# # #
# Set up some variables...
# # #
my $host   = '';
my $user   = '';
my $pass   = '';
my $port   = 3306;
my $dbname = '';
my $dnahost   = '';
my $dnauser   = 'ensro';
my $dnaport   = 3306;
my $dnadbname;

my $coord_system_name = 'toplevel';
my $coord_system_version = '';

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
  'infile:s'               => \$file,
  'write'                  => \$write,
  'verbose'                => \$verbose,
);

my $time = DateTime->today()->dmy;

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



# # # 
# hard code the mapping
# # # 
my $label2code = {
                   'PRINCIPAL:1'                 => [ 'appris', 'principal1' ],
                   'PRINCIPAL:2'                 => [ 'appris', 'principal2' ],
                   'PRINCIPAL:3'                 => [ 'appris', 'principal3' ],
                   'PRINCIPAL:4'                 => [ 'appris', 'principal4' ],
                   'PRINCIPAL:5'                 => [ 'appris', 'principal5' ],
                   'ALTERNATIVE:1'               => [ 'appris', 'alternative1' ],
                   'ALTERNATIVE:2'               => [ 'appris', 'alternative2' ],
                 };


# delete old attribs
if ( $write ) {
  print STDERR " Deleting old attributes...\n" if $verbose;
  foreach my $code ( values %{$label2code} ) {
    delete_old_attrib($db, $code->[0]);
  }
}


# read in file containing new attributes
my %appris_results;
my %appris_transcripts;# be lazy
open(INFILE, "<$file") or die ("Can't read $file $! \n");
while(my $line = <INFILE>){
  chomp $line;
  my ($gene_id, $transcript_id, $label) = split(/[ \t]+/, $line);
  $appris_transcripts{$transcript_id} = 1;

# "appris_principal", transcript(s) expected to code for the main functional isoform based on a range of protein features (APPRIS pipeline, Nucleic Acids Res. 2013 Jan;41(Database issue):D110-7).
# "appris_candidate", where there is no 'appris_principal' variant(s) the main functional isoform will be translated from one of the 'appris_candidate' genes. APPRIS selects one of the `appris_candidates' to be the most probable principal isoform, these have one of three labels:
 ## "appris_candidate_ccds", the `appris_candidate' transcript that has an unique CCDS.
 ## "appris_candidate_longest_ccds", the `appris_candidate' transcripts where there are several CCDS, in this case APPRIS labels the longest CCDS.
 ## "appris_candidate_longest_seq", where there is no 'appris_candidate_ccds' or ` appris_candidate_longest_ccds' variant, the longest protein of the 'appris_candidate' variants is selected as the primary variant.
  if (exists $label2code->{$label}) {
    $appris_results{$gene_id}{$transcript_id}{'attrib_type_code'} = $label2code->{$label}->[0];
    $appris_results{$gene_id}{$transcript_id}{'attrib_value'} = $label2code->{$label}->[1];
  } else {
    throw("Label $label not recognised in label2code hash");
  }
}
close INFILE;
print STDERR "Fetched ".(scalar(keys %appris_results))." genes and ".(scalar(keys %appris_transcripts))." transcripts\n";


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
  my $transc_no_data = 0;

  # now look for new candidates
  foreach my $gene ( @{$slice->get_all_Genes} ) {
    #print STDERR "Gene ".$gene->stable_id."\n";
    $gene_cnt++;
    foreach my $transcript ( @{$gene->get_all_Transcripts} ) {
      $transc_cnt++;
      
      if (exists $appris_results{$gene->stable_id}{$transcript->stable_id}) {
        # oh good, found this stable ID
        $stable_id_in_file++;
        if ($write) {
          store_attrib($aa, $transcript, $appris_results{$gene->stable_id}{$transcript->stable_id}{'attrib_type_code'},$appris_results{$gene->stable_id}{$transcript->stable_id}{'attrib_value'});
          print STDERR "  writing gene ".$gene->stable_id." transcript". $transcript->stable_id." APPRIS ".$appris_results{$gene->stable_id}{$transcript->stable_id}{'attrib_type_code'}."\n";
        }
      } else {
        # this is likely a new transcript that wasn't annotated last release
        $transc_no_data++;
        #warning("No data in file for ".$transcript->stable_id);
      }
    }
  }
  print  "Slice ".$slice->seq_region_name." has genes $gene_cnt with $transc_cnt transcripts. There are transcripts $transc_no_data with no attributes\n";
}
print "Matched stable_ids for ".$stable_id_in_file." of ".(scalar(keys %appris_results))." transcripts in file\n";
if ($stable_id_in_file != scalar(keys %appris_transcripts)) {
  throw("Not all transcripts found in database");
}
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


