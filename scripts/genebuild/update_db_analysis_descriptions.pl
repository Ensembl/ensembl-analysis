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

use warnings;
use strict;
use feature 'say';

use HTTP::Tiny;
use Time::HiRes qw/sleep/;
use JSON qw/decode_json/;
use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Data::Dumper;

my ($help, $dbname, $host, $port);

GetOptions(
	      'help|h'   => \$help,
	      'dbname=s' => \$dbname,
	      'host=s'   => \$host,
	      'port=s'   => \$port,
);

die &helptext if ( $help );

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
      -port    => $port,
      -user    => 'ensro',
      -host    => $host,
      -dbname  => $dbname);

open(OUT, '>', "./".$dbname."_update_ana_desc.sql");
print OUT "USE ".$dbname.";\n";

my $sth_logic = $db->dbc->prepare("select logic_name, analysis_id from analysis");
$sth_logic->execute;
my $http = HTTP::Tiny->new();
while (my @analysis_data = $sth_logic->fetchrow) {
  my $logic_name = $analysis_data[0];
  my $analysis_id = $analysis_data[1];

  my $server = 'http://production-services.ensembl.org';
  my $ext = '/production_db/api/analysisdescription/';
  my $response = $http->request('GET', $server.$ext.$logic_name, {
		     headers => {
		         'Content-type' => 'application/json',
				},
		 });
  if ($response->{success}){
    my $hash_ref = decode_json($response->{content});
    my %hash = %$hash_ref;

    local $Data::Dumper::Terse = 1;
    local $Data::Dumper::Indent = 0;
    my $web_data = Dumper($hash{'web_data'}->{data});
    if ($web_data eq 'undef') {
      $web_data = "NULL";
    }
    else {
      $web_data = '"'.$web_data.'"';
    }
    my $desc = $hash{'description'};
    $desc =~ s/\'/\\\'/g;;

    say "Creating SQL command for the analysis description table for logic_name ".$logic_name;
    my $insert = "INSERT INTO analysis_description (analysis_id, description, display_label, displayable, web_data) VALUES ($analysis_id, '$desc', '$hash{'display_label'}', $hash{'displayable'}, $web_data);";
    print OUT $insert."\n";

  }
  else{
    say "The logic_name ".$logic_name." does not exist in the Production database, should it be in the database ".$dbname."?";
  }

}

sub helptext {
  my $msg = <<HELPEND;

Create the file, dbname_update_ana_desc.sql - SQL commands to update the analysis descrition table in your db. 
The table will be populated with information form the Production database. If the analysis logic name does not exist in the Production db, then it will need to be added.

Usage: perl update_db_analysis_descriptions.pl -dbname db_name -host host -port port

To update your db: writable_host_connection_info db_name < dbname_update_ana_desc.sql

HELPEND
  return $msg;
}
