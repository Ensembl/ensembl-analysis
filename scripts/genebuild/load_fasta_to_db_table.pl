#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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
use Getopt::Long qw(:config no_ignore_case);

use Bio::EnsEMBL::Hive::DBSQL::DBAdaptor;
use Bio::EnsEMBL::IO::Parser::Fasta;
use parent ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');
use Data::Dumper;

my $dbname = '';
my $dbhost = '';
my $dbport = '';
my $dbuser = '';
my $dbpass = '';

my $fasta_file = "";
my $sequence_table_name;
my $create_table = 1;
my $force_uniq = 0;

my $result = GetOptions ("user|dbuser|u=s"      => \$dbuser,
                         "host|dbhost|h=s"      => \$dbhost,
                         "port|dbport|P=i"      => \$dbport,
                         "dbname|db|D=s"    => \$dbname,
                         "dbpass|pass|p=s"  => \$dbpass,
                         "fasta_file=s"   => \$fasta_file,
                         "sequence_table_name=s" => \$sequence_table_name,
                         "create_table!" => \$create_table,
                         "force_uniq_hitnames!" => \$force_uniq);


# Now connect to the pipe db
my $url = 'mysql://'.$dbuser.':'.$dbpass.'@'.$dbhost.':'.$dbport.'/'.$dbname;
my $db = Bio::EnsEMBL::Hive::DBSQL::DBAdaptor->new(-url => $url);

if($create_table) {
  if($force_uniq){
    my $sth_create_table = $db->dbc->prepare('CREATE table '.$sequence_table_name.' ('.
					     'accession varchar(255),'.
					     'biotype varchar(255),'.
					     'hit_name varchar(255),'.
					     'source varchar(255),'.
					     'seq text NOT NULL,'.
					     'PRIMARY KEY (accession))');
    $sth_create_table->execute();
  }
  else{
    my $sth_create_table = $db->dbc->prepare('CREATE table '.$sequence_table_name.' ('.
                                             'accession int NOT NULL AUTO_INCREMENT,'.
                                             'biotype varchar(255),'.
                                             'hit_name varchar(255),'.
                                             'source varchar(255),'.
                                             'seq text NOT NULL,'.
                                             'PRIMARY KEY (accession))');
    $sth_create_table->execute();
  }
}


my $table_adaptor = $db->get_NakedTableAdaptor();
$table_adaptor->table_name($sequence_table_name);

my $parser = Bio::EnsEMBL::IO::Parser::Fasta->open($fasta_file);
my $header;
my $seq;

my %uniq_accessions;

while($parser->next()) {
  $header = $parser->getHeader();
  $seq = $parser->getSequence();
  my ($biotype,$hit_name,$source) = split('\|',$header);

  if($force_uniq){
    my $accession = $hit_name;
    if(exists $uniq_accessions{$accession}){
      $uniq_accessions{$accession} += 1;
      $accession = $accession ."_". ($uniq_accessions{$accession});
    }
    else{
      $uniq_accessions{$accession} = 1;
    }

    my $db_row = [{
                    'accession'  => $accession,
                    'biotype'    => $biotype,
                    'hit_name'   => $hit_name,
                    'source'     => $source,
                    'seq'        => $seq,
                   }];
    $table_adaptor->store($db_row);
  }

  else{
    my $db_row = [{
                    'biotype'    => $biotype,
                    'hit_name'   => $hit_name,
    		    'source'     => $source,
                    'seq'        => $seq,
                 }];
    $table_adaptor->store($db_row);
  }
}
exit;
