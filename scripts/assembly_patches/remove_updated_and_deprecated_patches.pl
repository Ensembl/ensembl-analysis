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

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Getopt::Long;

use strict;
use warnings ;

my $pass;
my $patchtype_file = "./data/patch_type";
my $dbname;
my $host;
my $user;
my $port;
my $central_coord_system = 'supercontig';
my $toplevel_coord_system = 'chromosome';
my $sqlfile = 'delete_patch.sql';

&GetOptions(
            'pass=s'         => \$pass,
            'patchtype_file=s' => \$patchtype_file,
            'host=s'         => \$host,
            'dbname=s'       => \$dbname,
            'user=s'         => \$user,
            'port=n'         => \$port,
            'central_cs=s'   => \$central_coord_system,
            'sqlfile=s'   => \$sqlfile,
           );

open(TYPE,"<".$patchtype_file)         || die "Could not open file $patchtype_file";
open(SQL,">$sqlfile") || die "Could not open $sqlfile for writing\n";
#connect to the database

my $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    '-host'    => $host,
    '-port'    => $port,
    '-user'    => $user,
    '-pass'    => $pass,
    '-dbname'  => $dbname,
    '-species' => "load"
    );

#get the full list of existing synonyms/accessions
# This should be more robust than SQL queries but I hate to have these "PATCH" regex
my %existing_synonyms;
my $ae_adaptor = $dba->get_AssemblyExceptionFeatureAdaptor;
foreach my $assembly_exception (@{$ae_adaptor->fetch_all}) {
    if ($assembly_exception->type =~ /^PATCH\w+$/) {
        my $projection = $assembly_exception->slice->project($central_coord_system);
        foreach my $segment (@$projection) {
            my $slice = $segment->to_Slice;
            if ($slice->seq_region_name !~ /\w{2}\d+\.\d+/) {
                print 'Synonym ', $slice->get_all_synonyms('INSDC')->[0]->name, "\n";
                $existing_synonyms{$slice->get_all_synonyms('INSDC')->[0]->name} = $slice->seq_region_name;
            }
        }
    }
}

my $get_synonym_sth = $dba->dbc->prepare("select synonym from seq_region_synonym, seq_region, coord_system where seq_region. name = ? and seq_region.coord_system_id=coord_system.coord_system_id and coord_system.name = ? and seq_region.seq_region_id=seq_region_synonym.seq_region_id")
  || die "Could not prepare to get synonym";

if ($get_synonym_sth == 0) {
  throw("Could not prepare to get synonym");
}


#for the new files
while (<TYPE>) {
  chomp;
  next if (/^#/);
  my ($alt_scaf_name,$alt_scaf_acc,$type) = split(/\t/,$_);
  if (!$alt_scaf_name || !$alt_scaf_acc || !$type) {
    throw("Unable to find name, accession or type");
  }

  #check to see if this is an updated patch
  my $synonym;
  $get_synonym_sth->execute($alt_scaf_name, $central_coord_system) || die "problem executing get synonym";
  $get_synonym_sth->bind_columns(\$synonym) || die "problem binding";
  $get_synonym_sth->fetch();

  if(defined($synonym)){
    #if patch still exists (even if modified) record by setting to 0
    #modified patches are dealt with here already
    $existing_synonyms{$synonym} = 0;
    #unchanged
    if($synonym eq $alt_scaf_acc){
      print $alt_scaf_name." exists and is unchanged\n";
    }
    #modified
    else{
      print $alt_scaf_name." has been modified acc ".$synonym." to ".$alt_scaf_acc."\n";
      remove_patch($alt_scaf_name);
    }
  }
  #new
  else{
    print $alt_scaf_name." is a new patch\n";
  }

}
#removed
foreach my $acc (keys %existing_synonyms){
  if($existing_synonyms{$acc}){
    print $existing_synonyms{$acc}, " with acc $acc has been removed from the new set\n";
    remove_patch($existing_synonyms{$acc});
  }
}


sub remove_patch{

  my $alt_scaf_name = shift;
  #get the patch seq_region_ids (chrom and scaffold)
  my($scaf_id, $chrom_id);
  my $slice_adaptor = $dba->get_SliceAdaptor;
  my $slice = $slice_adaptor->fetch_by_region($central_coord_system, $alt_scaf_name);
  $scaf_id = $slice->get_seq_region_id;

  my $projection = $slice->project($toplevel_coord_system);
  foreach my $segment (@$projection) {
      my $projected_slice = $segment->to_Slice();
      if ($projected_slice->assembly_exception_type !~ '^PATCH') {
          warning("$alt_scaf_name which is on ".$projected_slice->seq_region_name.' is of type '.$projected_slice->assembly_exception_type."\n  You nay have a problem!");

      }
      else {
          $chrom_id = $projected_slice->get_seq_region_id;
      }
  }
  throw("Could not find the $toplevel_coord_system seq_region_id for $alt_scaf_name") unless ($chrom_id);

  #check for components only used in patch
  my $scaf_comp;
  my $get_region_components_sth = $dba->dbc->prepare("select distinct cmp_seq_region_id from assembly where asm_seq_region_id = ?") || die "Couldn't get components";
  $get_region_components_sth->execute($scaf_id) || die "problem executing get components";
  $get_region_components_sth->bind_columns(\$scaf_comp) || die "problem binding";
  COMP: while($get_region_components_sth->fetch()){
    #check each component
    my $asm_id;
    my $get_asm_ids_sth = $dba->dbc->prepare("select distinct asm_seq_region_id from assembly where cmp_seq_region_id = ?") || die "Couldn't get asm ids";
    $get_asm_ids_sth->execute($scaf_comp) || die "problem executing get components";
    $get_asm_ids_sth->bind_columns(\$asm_id) || die "problem binding";
    while($get_asm_ids_sth->fetch()){
      if($asm_id != $chrom_id and $asm_id != $scaf_id){
        #print "Used outside of the patch\n";
        next COMP;
      }
    }
    #component only used in patch
    #print $scaf_comp." only used in patch\n";
    print SQL "delete from dna where seq_region_id = ".$scaf_comp.";\n";
    print SQL "delete from seq_region_attrib where seq_region_id = ".$scaf_comp.";\n";
    print SQL "delete from seq_region where seq_region_id = ".$scaf_comp.";\n";
  }

  print SQL "delete from seq_region_synonym where seq_region_id = ".$scaf_id.";\n";

  print SQL "delete from seq_region_attrib where seq_region_id = ".$chrom_id.";\n";

  print SQL "delete from assembly_exception where seq_region_id = ".$chrom_id.";\n";

  print SQL "delete from assembly where asm_seq_region_id = ".$chrom_id.";\n";
  print SQL "delete from assembly where asm_seq_region_id = ".$scaf_id.";\n";

  print SQL "delete from seq_region where seq_region_id = ".$chrom_id.";\n";
  print SQL "delete from seq_region where seq_region_id = ".$scaf_id.";\n";

  print SQL "delete from marker_feature where seq_region_id = ".$chrom_id.";\n";
}

