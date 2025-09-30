#!/usr/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;


# Connection to the target DB
my $host;
my $port = '3306';
my $user;
my $pass;
my $dbname;
my $ucsc_host = 'genome-mysql.cse.ucsc.edu';
my $ucsc_port = '3306';
my $ucsc_user = 'genome';
my $ucsc_pass;
my $ucsc_dbname;
my $write = 0;
my $external_db_id;
my $external_dbname = 'UCSC';
my $query = 'SELECT chrom, name FROM ucscToINSDC';
my $help = 0;

&GetOptions (
            'h|host|dbhost=s'   => \$host,
            'P|port|dbport=i'   => \$port,
            'u|user|dbuser=s'   => \$user,
            'p|pass|dbpass=s'   => \$pass,
            'd|dbname=s'        => \$dbname,
            'ucsc_host=s'       => \$ucsc_host,
            'ucsc_port=i'       => \$ucsc_port,
            'ucsc_user=s'       => \$ucsc_user,
            'ucsc_pass=s'       => \$ucsc_pass,
            'ucsc_dbname=s'     => \$ucsc_dbname,
            'external_db_id=i'  => \$external_db_id,
            'external_dbname=s' => \$external_dbname,
            'write!'            => \$write,
            'help!'             => \$help,
        );

if ($help) {
  usage();
  exit(0);
}
if (!$host or !$user or !$dbname or !$ucsc_dbname) {
  usage();
  exit(1);
}

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
            '-host'   => $host,
            '-port'   => $port,
            '-user'   => $user,
            '-pass'   => $pass,
            '-dbname' => $dbname,
);

my $ucsc_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
            '-host'   => $ucsc_host,
            '-port'   => $ucsc_port,
            '-user'   => $ucsc_user,
            '-pass'   => $ucsc_pass,
            '-dbname' => $ucsc_dbname,
);

if (!$external_db_id) {
  $external_db_id = $db->get_DBEntryAdaptor->get_external_db_id($external_dbname);
}
my @synonyms;
my $slice_adaptor = $db->get_SliceAdaptor;
my $sth = $ucsc_db->dbc->prepare($query);
$sth->execute();
while (my $row = $sth->fetchrow_arrayref) {
  my $slice = $slice_adaptor->fetch_by_region(undef, $row->[0]);
  if (!$slice) {
    $slice = $slice_adaptor->fetch_by_region(undef, $row->[1]);
    if ($slice) {
      push(@synonyms, Bio::EnsEMBL::SeqRegionSynonym->new(-synonym => $row->[0],
                                                          -seq_region_id => $slice->get_seq_region_id,
                                                          -external_db_id => $external_db_id)
      );
    }
    else {
      die("Could not find region $row->[0]/$row->[1] in your database");
    }
  }
}

if ($write) {
  my $synonym_adaptor = $db->get_SeqRegionSynonymAdaptor;
  foreach my $synonym (@synonyms) {
    $synonym_adaptor->store($synonym);
  }
  my $meta_adaptor = $db->get_MetaContainerAdaptor;
  my $ucsc_syn = $meta_adaptor->single_value_by_key('assembly.ucsc_alias');
  if ($ucsc_syn) {
    if ($ucsc_syn ne $ucsc_dbname) {
      print "Updating 'assembly.ucsc_alias' from $ucsc_syn to $ucsc_dbname\n";
      $meta_adaptor->update_key_value('assembly.ucsc_alias', $ucsc_dbname);
    }
  }
  else {
    $meta_adaptor->store_key_value('assembly.ucsc_alias', $ucsc_dbname);
  }
}

sub usage {
  print <<EOF

 load_UCSC_synonyms.pl -host <host> -user <user_rw> -pass <pass> -port <port> -dbname <dbname> -ucsc_dbname <ucsc_dbname> [-write]
    [-ucsc_host <host>] [-ucsc_user <user_rw>] [-ucsc_pass <pass>] [-ucsc_port <port>] [-external_db_id <ext db id>] [-external_dbname <ext dbname>]
    [-help]
  h|host|dbhost    Name of the host of your Ensembl database
  P|port|dbport    Port of your Ensembl database, default '3306'
  u|user|dbuser    User with write access
  p|pass|dbpass    Password
  d|dbname         Name of the Ensembl database
  ucsc_host        UCSC server, default 'genome-mysql.cse.ucsc.edu'
  ucsc_port        Port, default '3306'
  ucsc_user        User, default 'genome'
  ucsc_pass        Password
  ucsc_dbname      Name of the UCSC database, it corresponds to the UCSC alias and will be added if not present
  external_db_id   Id of the external database, if not set it will query the Ensembl database using 'external_dbname'
  external_dbname  Name of the database to assign the synonyms, default 'UCSC'
  write            Write the synonyms and the 'assembly.ucsc_alias' meta key in the Ensembl database
  help             This message

 The script will throw if a region from the UCSC database cannot be found in the Ensembl database as this should no happen!

EOF
}
