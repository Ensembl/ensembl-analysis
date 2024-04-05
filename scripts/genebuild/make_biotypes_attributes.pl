#!/usr/bin/env perl

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

use warnings ;
use strict;
use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::Utilities;
use Bio::EnsEMBL::Utils::Exception;
use Bio::EnsEMBL::Analysis::Config::Databases;

my $host;
my $port=3306;
my $dbname;
my $user;
my $pass;
my $config_dbname;
my $seq_region;
my $coord_system = 'toplevel';
my $transcript = 1;
my $gene = 1;
my $write;

GetOptions( 
            'dbhost|host|h:s'      => \$host,
            'dbport|port|P:n'      => \$port,
            'dbname|db|D:s'    => \$dbname,
            'dbuser|user|u:s'    => \$user,
            'dbpass|pass|p:s'      => \$pass,
            'config_dbname:s' => \$config_dbname,
            'seq_region_name:s' => \$seq_region,
            'coord_system:s' => \$coord_system,
            'write!' => \$write,
           );


my $db;

if($config_dbname){
  $db = get_db_adaptor_by_string($config_dbname);
}elsif($dbname && $host){
  $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                            -host   => $host,
                                            -user   => $user,
                                            -port   => $port,
                                            -dbname => $dbname,
                                            -pass => $pass,
                                          );
}else{
  throw("Need to pass either -dbhost $host and -dbname $dbname or ".
        "-config_dbname $config_dbname for the script to work");
}

print "Have database ".$db."\n";
my $sa = $db->get_SliceAdaptor;

my $slices;
if($seq_region){
  my $slice = $sa->fetch_by_region($coord_system, $seq_region);
  push(@$slices, $slice);
}else{
  $slices = $sa->fetch_all($coord_system);
}
print "Have ".@$slices." slices\n";
my %attributes;
my $aa = $db->get_AttributeAdaptor;

SLICE:while(my $slice = shift(@$slices)){
  print "Looking at ".$slice->name."\n";
  my $genes = $slice->get_all_Genes;
  print "Have ".@$genes." genes\n";
 GENE:while(my $gene = shift(@$genes)){
    if(!$attributes{$gene->biotype}){
      my $attrib = create_Attribute($gene->biotype);
      $attributes{$gene->biotype} = $attrib;
    }
    my $attrib = $attributes{$gene->biotype};
    $gene->add_Attributes($attrib);
    if($write){
      $aa->store_on_Gene($gene->dbID, [$attrib]);
    }
    foreach my $transcript(@{$gene->get_all_Transcripts}){
      if(!$attributes{$transcript->biotype}){
        my $attrib = create_Attribute($transcript->biotype);
        $attributes{$transcript->biotype} = $attrib;
      }
      my $attrib = $attributes{$transcript->biotype};
      $transcript->add_Attributes($attrib);
      if($write){
        $aa->store_on_Transcript($transcript->dbID, [$attrib]);
      }
    }
  }
}


sub create_Attribute{
  my ($reason, $adaptor) = @_;
  #print "Creating attribute for ".$reason."\n";
  my $sql = "select attrib_type_id from attrib_type where code = ?";
  my $sth = $aa->dbc->prepare($sql);
  my $length = length($reason);
  my $code;
  if($length >= 16){
  CODE:foreach (my $offset = 15;$offset > 0; $offset--){
      #print "Looking at ".$offset." compared to ".$length."\n";
      my $real_offset = $offset - $length;
      $code = substr($reason, 0, $real_offset);
      throw("Still have a ".$code." which is too long from ".$reason." ".
            "using ".$real_offset) if(length($code) >= 16);
      #print "Executing sql\n";
      $sth->execute($code);
      #print "RUNNING ".$sql." with ".$code."\n";
      my ($id) = $sth->fetchrow;
      #print "HAVE found ".$id." id\n" if($id);
      last CODE if(!$id);
    }
  }else{
    $code = $reason;
  }
  my $attrib = Bio::EnsEMBL::Attribute->new(
                                            -code => $code,
                                            -name => $reason,
                                            -description => $reason,
                                            -value => $reason,
                                           );
  return $attrib;
}
