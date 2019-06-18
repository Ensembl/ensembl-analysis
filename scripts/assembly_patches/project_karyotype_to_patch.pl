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
my $coord_sys_ver = '';
my @patch_types = ('PATCH_FIX','PATCH_NOVEL');
my $sql_file = '';

&GetOptions( 'mapdbhost:s'       => \$host,
             'mapdbuser:s'       => \$user,
             'mapdbname:s'       => \$dbname,
             'mapdbport:n'       => \$port,
             'sql_file:s'        => \$sql_file,
             'coord_sys_ver:s'   => \$coord_sys_ver);

open (SQL, ">", $sql_file);

#database with mapping between coord systems
my $mapdb = new Bio::EnsEMBL::DBSQL::DBAdaptor( -host   => $host,
                                             -user   => $user,
                                             -port   => $port,
                                             -dbname => $dbname );

my $map_sa = $mapdb->get_SliceAdaptor();
my $ka = $mapdb->get_KaryotypeBandAdaptor();

#get the patches
print "Getting patches...\n";
my $asm_exc_adaptor = $mapdb->get_AssemblyExceptionFeatureAdaptor();
my @exceptions = @{$asm_exc_adaptor->fetch_all()};

print scalar(@exceptions)." assembly exceptions were returned\n";

my @patches;
EXC: foreach my $exc (@exceptions){
  foreach my $type (@patch_types){
    if($exc->type() eq $type){
      push(@patches, $exc);
      next EXC;
    }
  }
}

print "Have ".scalar(@patches)." patches.\n";
#for each patch
foreach my $patch_exc (@patches){
  #get the slice
  my $patch_slice = $patch_exc->feature_Slice;
  print "\n\n".$patch_slice->name."\n";
  my $ref_slice = $patch_exc->alternate_slice;
  print $ref_slice->name."\n";
  my $coord_system_name = $patch_slice->seq_region_name;
  
  #get and sort the bands
  my @ref_bands = @{$ka->fetch_all_by_Slice($ref_slice)};
  print scalar(@ref_bands)." bands returned\n";
  @ref_bands = sort {$a->start cmp $b->start} @ref_bands;

  my $curr_patch_max_coord;
  my $band_count = 1;

  foreach my $band (@ref_bands){

    my $patch_start;
    my $patch_end;
    #check where band is on ref
    my $band_ref_start = ($ref_slice->start-1)+$band->start;
    my $band_ref_end = ($ref_slice->start-1)+$band->end;
    print $band->name."\n";
    print "ref start: ".$band_ref_start." ref end: ".$band_ref_end."\n";

    #project and work out extent of projection
    my $band_max;
    my $band_min;
    my $first_seg = 1;
    my @band_proj_segs = @{$band->project($coord_system_name, $coord_sys_ver)};
    print scalar(@band_proj_segs)." proj segs returned\n";
    foreach my $ps (@band_proj_segs){
      my $ps_slice = $ps->to_Slice;
      if($first_seg){
        $band_min = $ps_slice->start;
        $band_max = $ps_slice->end;
        $first_seg = 0;
      }
      else{
        if($ps_slice->start < $band_min){
          $band_min = $ps_slice->start;
        }
        if($ps_slice->end > $band_max){
          $band_max = $ps_slice->end;
        }
      }
      print "seg ".$ps->to_Slice->name."\n" if $patch_slice->name =~ /MANN/;
    }
    print "min:$band_min max:$band_max\n";

    #work out patch chrom start and end
    #if it's the first band
    if($band_count == 1){
      print "FIRST BAND\n";
      $patch_start = $band_ref_start;
      $patch_end = $band_max;
      $curr_patch_max_coord = $patch_end;
    }
    else{
      print "MID BAND\n";
      $patch_start = $curr_patch_max_coord+1;
      $curr_patch_max_coord++;
      $patch_end = $band_max;
      $curr_patch_max_coord = $patch_end;
    }
    #if it is the last band
    if($band_count == scalar(@ref_bands)){
      print "LAST BAND\n";
      #print "executing last band code\n";
      $patch_end = $band_ref_end + ($patch_slice->end - $ref_slice->end);
    }
    $band_count++;
    print "patch start: ".$patch_start." patch end: ".$patch_end."\n";
    
    if ($patch_start > $patch_end) {
      print "Karyotype band has been mapped via the reverse orientation.\n" .
            "Since the bands are not a stranded feature we will flip the orientation:\n\n";
      my $tmp = $patch_start;
      $patch_end = $patch_start;
      $patch_start = $tmp;

      print "patch start: ".$patch_start." patch end: ".$patch_end."\n";
    }

    print SQL "insert into karyotype(seq_region_id, seq_region_start, seq_region_end, band, stain) values(".$patch_slice->get_seq_region_id.", ".$patch_start.", ".$patch_end.", '".$band->name."', '".$band->stain."');\n";
  }
}

print "Completed\n";
