#!/usr/env perl
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

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;

my ($dbname, $dbhost, $dbport, $dbuser, $dbpass, $mirnas_list);

GetOptions( 'dbhost|host|h:s'        => \$dbhost,
            'dbport|port|p:n'        => \$dbport,
            'dbname|d:s'        => \$dbname,
            'dbuser|user|u:s'        => \$dbuser,
            'password|e:s' => \$dbpass,
            'rnaproducts|r:s' => \$mirnas_list);

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
	-DBNAME => $dbname,
  	-HOST => $dbhost,
  	-PORT => $dbport,
  	-USER => $dbuser,
    -PASS => $dbpass,
	-DRIVER => 'mysql',
);

my $command = "mysql -h $dbhost -P $dbport -u $dbuser -p$dbpass -D $dbname -e \"insert ignore into rnaproduct_type (code, name, description) values ('miRNA', 'Mature miRNAs', 'MiRBase annotated miRNA products mapped to identified stem-loops');
alter table rnaproduct change created_date created_date datetime DEFAULT CURRENT_TIMESTAMP;\"";

system($command);

my $rpa = $db->get_RNAProductAdaptor();

open(FH, $mirnas_list) or die "Could not open provided list of mirnas";

my @line;

while(<FH>){
  @line = split(/\t/,$_);
  my @temp = split(/:/, $line[0]);
  my $precursor_dbid = $temp[0];

  my $rnaproduct = Bio::EnsEMBL::MicroRNA->new(
        -SEQ_START => $line[1],
        -SEQ_END   => $line[2],
        -ARM => $line[3]
  );

  $rpa->store($rnaproduct, $precursor_dbid);

}

close(FH);
