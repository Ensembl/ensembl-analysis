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

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor;
use Bio::EnsEMBL::AssemblyExceptionFeature;

$| = 1;

my $dbname       = '';
my $dbhost    = '';
my $dbuser    = '';
my $dbport    = '';
my $mapping_analysis_id = '';
my @patch_types = ('PATCH_FIX','PATCH_NOVEL');
my $patch_cmp_in_mapping;
my $ref_cmp_in_mapping;
my $external_db_id = 50633; # GRC_primary_assembly

&GetOptions( 'dbhost:s'                 => \$dbhost,
             'dbuser:s'                 => \$dbuser,
             'dbname:s'                 => \$dbname,
             'dbport:n'                 => \$dbport,
             'mapping_analysis_id:n'    => \$mapping_analysis_id,
             'patch_cmp_in_mapping'     => \$patch_cmp_in_mapping,
             'ref_cmp_in_mapping'       => \$ref_cmp_in_mapping);


if($patch_cmp_in_mapping){
  print "Looking for ref chromosomes as asm and patch chromosomes as cmp in mapping.\n";
}
if($ref_cmp_in_mapping){
  print "Looking for patch chromosomes as asm and ref chromosomes as cmp in mapping.\n";
}
if($patch_cmp_in_mapping and $ref_cmp_in_mapping){
  throw("Ref or patch must be asm and other cmp, not both.\n");
}
if(!$ref_cmp_in_mapping and !$patch_cmp_in_mapping){
  throw("Specify -patch_cmp_in_mapping or -ref_cmp_in_mapping to indicate if patch or ref is the cmp in the mapping between the two.\n");
}
#database with mapping between coord systems
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $dbhost,
                                             -user   => $dbuser,
                                             -port   => $dbport,
                                             -dbname => $dbname );

if(!$mapping_analysis_id){
  my $analysis_adaptor = $db->get_AnalysisAdaptor;
  my $analysis = $analysis_adaptor->fetch_by_logic_name('alt_seq_mapping');
  $mapping_analysis_id = $analysis->dbID;
  throw ("Need to specify analysis id mapping will be stored with using -mapping_analysis_id.\n") unless ($mapping_analysis_id);
}

my $sa = $db->get_SliceAdaptor();

#get the patches
print "Getting patches...\n";
my $asm_exc_adaptor = $db->get_AssemblyExceptionFeatureAdaptor();
my @exceptions = @{$asm_exc_adaptor->fetch_all()};
my @patches;
EXC: foreach my $exc (@exceptions){
  foreach my $type (@patch_types){
    if($exc->type() =~ m/$type/){
      push(@patches, $exc);
      next EXC;
    }
  }
}
#Assuming that AssemblyExceptionFeatureAdaptor's fetch_all will always return two 
#entries for each patch and that those two entries are adjacent
my $num_patches = scalar(@patches)/2;
print "Have ".$num_patches." patches.\n";

open(SQL,">mapping.sql") || die "Could not open mapping.sql for writing\n";

#for each patch
for (my $i = 0;$i < $num_patches;$i++){

  #get the two slices
  my $ref_slice;
  my $patch_slice;
  for(my $j = 0;$j < 2; $j++){
    my $exc = pop(@patches);
    #if this is the ref version
    if($exc->type =~ m/REF/){
      #alt is only the patch slice
      $patch_slice = $exc->alternate_slice();
    }
    else{
      #alt is replaced region of ref
      $ref_slice = $exc->alternate_slice();
    }    
  }
  if(!($patch_slice and $ref_slice)){
    throw("Something is wrong, the patch and ref slices were not set correctly.\n");
  }

  my $patch_chrom_id = $sa->get_seq_region_id($patch_slice);
  
  #cmp or asm in mapping?
  my $get_mapping_sth;
  my $get_primary_sth;#primary assembly
  my($seq_region_id, $start, $end, $primary_id, $pri_start, $pri_end, $ori);

  if($patch_cmp_in_mapping){
    $get_mapping_sth = $db->dbc->prepare("select cmp_seq_region_id, cmp_start, cmp_end, asm_seq_region_id, asm_start, asm_end, ori from assembly, seq_region, coord_system where cmp_seq_region_id = ? and asm_seq_region_id = seq_region_id and seq_region.coord_system_id=coord_system.coord_system_id and coord_system.name='chromosome'")
  || die "Couldn't get mapping";
  }
  elsif($ref_cmp_in_mapping){
    $get_mapping_sth = $db->dbc->prepare("select asm_seq_region_id, asm_start, asm_end, cmp_seq_region_id, cmp_start, cmp_end, ori from assembly, seq_region, coord_system where asm_seq_region_id = ? and cmp_seq_region_id = seq_region_id and seq_region.coord_system_id=coord_system.coord_system_id and coord_system.name='chromosome'")
  || die "Couldn't get mapping";
  }
  #get mapping
  $get_mapping_sth->execute($patch_chrom_id) || die "problem executing get mapping";
  $get_mapping_sth->bind_columns(\$seq_region_id, \$start, \$end, \$primary_id, \$pri_start, \$pri_end, \$ori) || die "problem binding";

  my $mapping_counter = 0;

  MAPPING: while($get_mapping_sth->fetch()){
    $mapping_counter++;

    #print $seq_region_id." ".$start." ".$end."\n";
    #is mapping outwith the patch region of the alt chromosome?
    if($end < $patch_slice->start){
      #print "Before\n";
      next MAPPING;
    }
    if($start > $patch_slice->end){
      #print "After\n";
      next MAPPING;
    }
    #does mapping extend beyond patch?
    if($end >= $patch_slice->start and $start < $patch_slice->start){
      #print "Start boundary\n";
      my $diff = $patch_slice->start - $start;
 
     #which way round is the primary relative to the patch?
      if($ori == 1){
        $pri_start = $pri_start + $diff;
      }
      else{
        $pri_end = $pri_end - $diff;
      }
      $start = $patch_slice->start;

    }
    if($start <= $patch_slice->end and $end > $patch_slice->end){
      #print "End boundary\n";
      my $diff = $end - $patch_slice->end;

      #which way round is patch
      if($ori == 1){
        $pri_end = $pri_end - $diff;
      }
      else{
        $pri_start = $pri_start + $diff;
      }
      $end = $patch_slice->end;

    }
    my $length = $end - ($start - 1);
    my $pri_length = $pri_end - ($pri_start - 1);

    if($length != $pri_length){
      throw('Lengths should be the same!');
    }

    my $hit_name = $ref_slice->seq_region_name;

    print SQL "insert into dna_align_feature (analysis_id, seq_region_id, seq_region_start, seq_region_end, seq_region_strand, hit_start, hit_end, hit_strand, hit_name, cigar_line,external_db_id)
values(".$mapping_analysis_id.", ".$patch_chrom_id.", ".$start.", ".$end.", 1, ".$pri_start.", ".$pri_end.", ".$ori.", '".$hit_name."', '".$length."M', $external_db_id);\n";

  }
  print $patch_slice->name." ".$mapping_counter." mapping entries processed\n";
}

print "Completed\n";
