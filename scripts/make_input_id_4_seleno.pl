#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

use warnings;
use strict;

use Getopt::Long qw(:config no_ignore_case);


use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;

my $host;
my $port=3306;
my $user;
my $pass;
my $dbname;
my $input_id_type;
my $analysis_id;

GetOptions(
            'host|dbhost|h=s'               => \$host,
            'port|dbport|P=s'               => \$port,
            'user|dbuser|u=s'               => \$user,
            'pass|dbpass|p=s'               => \$pass,
            'dbname|db|D=s'             => \$dbname, 
            'input_id_type=s'      => \$input_id_type,
            'analysis_id=s'        => \$analysis_id,
            );


$|=1;

my $db = new Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor(-dbname => $dbname,
                                                      -user => $user,
                                                      -pass => $pass,
                                                      -port => $port,
                                                      -host => $host);

my $q= 'select seq_region.name, coord_system.version, gene.seq_region_start, gene.seq_region_end, gene.seq_region_strand from gene, seq_region, coord_system where seq_region.seq_region_id = gene.seq_region_id and seq_region.coord_system_id = coord_system.coord_system_id order by gene.seq_region_id, gene.seq_region_start, gene.seq_region_end, gene.seq_region_strand';

my $sth = $db->prepare($q) || $db->throw("can't prepare: $q");

my $slice_start = 100000000000000000000000000000000;
my $slice_end = 0;
my $slice_chr = 0;
my $slice_strand = 0;

#while(<>){
#  chomp;
#  my $res = $sth->execute($_) || $db->throw("can't execute: $q");
   my $res = $sth->execute() || $db->throw("can't execute: $q");
 
  while( my ($chr, $version, $start, $end, $strand) = $sth->fetchrow_array) {
    next if ($chr eq 'MT');
    if ($slice_chr == 0){
      $slice_chr = $chr;
      $slice_start = $start;
      $slice_end = $end;
      $slice_strand = $strand;
    }
    if ($slice_chr eq $chr && $start < $slice_end && $end >$slice_start && $slice_strand == $strand){
      if ($start < $slice_start){ $slice_start = $start;}
      if ($end > $slice_end){ $slice_end = $end;}
    }else{
      print "insert into input_id_analysis values ('toplevel:".$version.":$chr:$start:$end:1','".$input_id_type."',".$analysis_id.", now(), '', '', '0');\n";
      $slice_chr = $chr;
      $slice_start = $start;
      $slice_end = $end;
      $slice_strand = $strand;   
    }
  }
#}

