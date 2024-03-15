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

my $host         = '';
my $user         = '';
my $port         = '';
my $dbname       = '';
my @patch_types = ('PATCH_FIX','PATCH_NOVEL');

&GetOptions( 'dbhost:s'       => \$host,
             'dbuser:s'       => \$user,
             'dbname:s'       => \$dbname,
             'dbport:n'       => \$port);

my $mapdb = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $host,
                                             -user   => $user,
                                             -port   => $port,
                                             -dbname => $dbname );



my $map_sa = $mapdb->get_SliceAdaptor();

#get the patches
print "Getting patches...\n";
my $asm_exc_adaptor = $mapdb->get_AssemblyExceptionFeatureAdaptor();
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

my $total_patch = 0;
my $total_patch_mapped = 0;
my $total_segs = 0;

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

  my @patch_segs = @{$patch_slice->project('chromosome')};

  print $patch_slice->name()." has ".scalar(@patch_segs)." segs mapping to ".$ref_slice->name()."\n";

  my $length = 0;

  foreach my $p_seg (@patch_segs){
    my $patch_seg_slice = $p_seg->to_Slice;
    $length = $length + $patch_seg_slice->length;
    print "Proj seg: ".$patch_seg_slice->name."\n";
  }
  my $perc = ($length/$patch_slice->length)*100;
  print $perc."% of patch bases mapped from ".$patch_slice->name()."\n"; 
  $total_patch = $total_patch + $patch_slice->length;
  $total_patch_mapped = $total_patch_mapped + $length;
  $total_segs = $total_segs + scalar(@patch_segs);
}
my $total_perc = ($total_patch_mapped/$total_patch)*100;
print $total_perc."% of patch bases mapped in ".$total_segs." segments.\n";


print "Completed\n";
