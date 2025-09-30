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

# script to calculate map_weights in a database that has markers
# and marker_features. Recreates the marker_feature table with weights set

use warnings ;
use strict;
use DBI;

use Getopt::Long qw(:config no_ignore_case);

my ( $host, $user, $pass, $port, $dbname );
my $verbose = 0;
$port = 3306;

GetOptions( "host|dbhost|h=s", \$host,
	    "user|dbuser|u=s", \$user,
	    "pass|dbpass|p=s", \$pass,
	    "port|dbport|P=i", \$port,
	    "dbname|db|D=s", \$dbname,
	    "verbose", \$verbose
	  );

if( !$host ) {
  usage();
}



my $dsn = "DBI:mysql:host=$host;dbname=$dbname";
if( $port ) {
  $dsn .= ";port=$port";
}

my $db = DBI->connect( $dsn, $user, $pass );


$db->do( "
  CREATE TABLE tmp_m_weight
  SELECT marker_id, count(*) as count 
  FROM marker_feature
  GROUP BY marker_id
" );

$db->do( "
  CREATE TABLE new_marker_feature
  SELECT mf.marker_feature_id, mf.marker_id, mf.seq_region_id, mf.seq_region_start,
         mf.seq_region_end, mf.analysis_id, tmw.count
  FROM   marker_feature mf, tmp_m_weight tmw
  WHERE  mf.marker_id = tmw.marker_id
" );

$db->do( "delete from marker_feature" );
$db->do( "insert into marker_feature select * from new_marker_feature" );
$db->do( "drop table tmp_m_weight" );
$db->do( "drop table new_marker_feature" );

sub usage {
  print <<EOF;
    
Usage: perl map_weight.pl [options]
   -user username for a write enabled user
   -host hostname
   -port portnumber
   -pass password
   -dbname database name where the markers and the features are

EOF

  exit;
}
