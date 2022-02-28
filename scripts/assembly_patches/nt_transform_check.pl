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

# This script works around the three levels of the assembly checking that
# it is possible to map the whole way around the triangle. The script deals with
# the branching that occurs when the patches are in place (i.e. two mappings
# at one level of the assembly).

use strict;
use warnings;

use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Feature;

my $host;
my $user;
my $pass;
my $dbname;
my $port;
my $central_coordsystem;
my $chromosome_name;
my $path;
my $verbose = 0;

$| = 1;

&GetOptions(
  'host:s'   => \$host,
  'user:s'   => \$user,
  'dbname:s' => \$dbname,
  'port:n'   => \$port,
  'pass:s'   => \$pass,
  'path:s'   => \$path,
  'central_coordsystem:s' => \$central_coordsystem,
  'chromosome_name:s' => \$chromosome_name,
  'verbose:s' => \$verbose,
);


# Open database
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -host   => $host,
  -user   => $user,
  -pass   => $pass,
  -port   => $port,
  -dbname => $dbname,
  -path   => $path,
);

my $sa = $db->get_SliceAdaptor;


my @chrslices = $sa->fetch_all('chromosome', undef, 1);
if ($chromosome_name) {
  @chrslices = ();
  push(@chrslices,$sa->fetch_by_region('chromosome',$chromosome_name));
}

print "Done fetch\n";

my @chr_name;
my @multi_match_contigs;
my @no_match_contigs;
my @part_match_contigs;

foreach my $chrslice (@chrslices) {
  print "Doing chr " . $chrslice->name . "\n";

  # this is for only checking a specific chr
   #next unless ($chrslice->name eq "chromosome:GRCh37:HG417_PATCH:1:342635:1");
   next if ($chrslice->name =~ /$path:Y/);
   next if ($chrslice->name =~ /$path:MT/);
  #next if ($chrslice->name ne "chromosome:GRCh37:HSCHR5_1_CTG1:68505250:70125573:1");#"chromosome:GRCh37:5:1:180915260:1");

  push @chr_name, $chrslice->name;
  my @multi_match_c;
  my @no_match_c;
  my @part_match_c;

  my $ctgsegs = $chrslice->project('contig');

  foreach my $ctgseg (@$ctgsegs) {
    my $ctgslice = $ctgseg->[2];
    print "Slice: " . $ctgslice->name . "\n" if ($verbose);
    my $f = new Bio::EnsEMBL::Feature( -start  => 1,
                                       -end    => $ctgslice->length,
                                       -strand => 1,
                                       -slice  => $ctgslice); 

    my $matched = 0;
    my $part_match = 0;
    my @feats = check_feature($f);
    print "Feats: ".scalar(@feats)."\n" if ($verbose);
    if(scalar(@feats)>1){
      print "Multiple features returned\n" if ($verbose);
    }
    foreach my $f (@feats){
      print "Returned feature start: ".$f->start." end: ".$f->end." strand: ".$f->strand." ".$f->slice->name."\n" if ($verbose);
      my $f_slice = $f->feature_Slice;
      if($f_slice->name eq $ctgslice->name){
        print "Contig ".$ctgslice->name." matched\n" if ($verbose);
        $matched++;
      }#check if partial match
      elsif(($f_slice->seq_region_name eq $ctgslice->seq_region_name) and
            ($f_slice->strand == $ctgslice->strand) and
            ($f_slice->start <= $ctgslice->end) and
            ($f_slice->end >= $ctgslice->start) and
            ($f_slice->coord_system_name eq $ctgslice->coord_system_name)){
        print "Partial match\n" if ($verbose);
        print $f_slice->name."\n" if ($verbose);
        $part_match++;
        push @part_match_c, $f_slice;
      }
    }
    if(!$matched){
      print "NOTHING fully matched contig ".$ctgslice->name."\n" if ($verbose);
      push @no_match_c, $ctgslice;
    }
    else{
      print $matched." contigs matched ".$ctgslice->name."\n" if ($verbose);
      if($matched > 1){
        print "There were multiple matches to the contig\n" if ($verbose);
        push @multi_match_c, $ctgslice;
      }
    }
    
    print "\n" if ($verbose);
  }
  push @multi_match_contigs, [@multi_match_c];
  push @no_match_contigs, [@no_match_c];
  push @part_match_contigs, [@part_match_c];
}

#print summary
my $i = 0;
foreach my $chr_n (@chr_name) {
  print "Chromosome " . $chr_n . "\n";
  print scalar(@{$multi_match_contigs[$i]})." contigs with multiple matches.\n";
  foreach my $c (@{$multi_match_contigs[$i]}){
    print $c->name."\n";
  }
  print scalar(@{$part_match_contigs[$i]})." partial matches to contigs (not num contigs with partial matches).\n";
  foreach my $fs (@{$part_match_contigs[$i]}){
    print $fs->name." was a partial match to a contig\n";
  }

  print scalar(@{$no_match_contigs[$i]})." contigs with no matches.\n";
  foreach my $c (@{$no_match_contigs[$i]}){
    print "PROBLEM: ".$c->name." - this shouldn't happen\n";
  }
  if(scalar(@{$part_match_contigs[$i]})==0 and scalar(@{$multi_match_contigs[$i]})==0){
    print $chr_n." MAPPED ONLY TO SELF\n\n";
  }
  $i++;
}

sub check_feature{
  my $f = shift;
  print "Check start: ".$f->start." end: ".$f->end." strand: ".$f->strand." ".$f->slice->name."\n" if ($verbose);
  my @features;
  my $coord_sys = $f->coord_system_name;
  my $next_sys = '';
  if($coord_sys eq 'contig'){
    $next_sys = $central_coordsystem;
  }
  elsif($coord_sys eq $central_coordsystem){
    $next_sys = 'chromosome';
  }
  elsif($coord_sys eq 'chromosome'){
    $next_sys = 'contig';
  }
  else{
    print "Unrecognised coord system\n" if ($verbose);
  }

  #first try to transform to next coord_system
  my $next_f = $f->transform($next_sys);
  #if that didn't work, it could be that there are multiple items 
  #on that coord system that the feature may project to
  #in that case we want to investigate them all
  if(!defined($next_f)){
    my $projection = $f->project($next_sys);
    if(@$projection > 1){
      print "Multiple options on coord system $next_sys\n" if ($verbose);
      foreach my $proj_seg(@$projection){
        #print "PROJ_SEG:".$proj_seg->from_start." ".$proj_seg->from_end." ".$proj_seg->to_Slice->name."\n";
        if ($proj_seg->to_Slice->start() > $proj_seg->to_Slice->end()+1) {
          print "Projection segment ignored because start (".$proj_seg->to_Slice->start().") is greater than end+1 (".$proj_seg->to_Slice->end()."+1). ".$proj_seg->to_Slice->name."\n" if ($verbose);
        } else {
          my $proj_f = make_feature($proj_seg->to_Slice);

          if ($next_sys ne 'contig') {
            push @features, check_feature($proj_f);
          } else {
            push @features, $proj_f;
          }
        }
      }
    }
    elsif(@$projection == 1){
      print "Feature that didn't transform only had one projection - not expected\n" if ($verbose);
    }
  }elsif($next_sys ne 'contig'){
    push @features, check_feature($next_f);
  }
  else{
    push @features, $next_f;
  }
  return @features;
}

#lifted from the end of Feature.pm transform
sub make_feature{
  my $p_slice = shift;
  my $slice_adaptor = $db->get_SliceAdaptor;
  my $slice = $slice_adaptor->fetch_by_region($p_slice->coord_system()->name(),
             $p_slice->seq_region_name(),
             undef, #start
             undef, #end
             1, #strand
             $p_slice->coord_system()->version);

  my $new_feature = new Bio::EnsEMBL::Feature( -start  => $p_slice->start(),
                                               -end    => $p_slice->end(),
                                               -strand => $p_slice->strand(),
                                               -slice  => $slice);
  return $new_feature;
}
