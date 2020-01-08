#!/usr/bin/env perl

# Analyze the alignment of the ditag features
# Identify & remove ditags from areas where they align in both directions
# (on both strands). Like this:
#
#       v
#       ----->>
#       ------>>
#__________________
#
# <<-----
#  <<----
#
#or this:
#
#       v
#       ----->>
#       ------>>
# <<-----
#  <<----
#__________________
#

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


use warnings ;
use strict;
use Bio::EnsEMBL::Analysis::Tools::Default;
use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
$| = 0;

my $dbhost         = "";
my $dbname         = "";
my $dbport         = 3306;
my $dbuser         = "ensro";
#use only specific types:
my $ditagtype      = undef; #"ZZ13";
my $ditag_analysis = "FANTOM_GSC_PET";

my (@dfs, @check_this, @ditagids);
my (%ditag_pos);
my $slicename;
my $verbose        = 0;
my @slices;
my $slice;

GetOptions(
	   'chromosome:s'      => \$slicename,
	   'verbose!'          => \$verbose,
           'dbhost|host|h:s'          => \$dbhost,
           'dbname|db|D:s'          => \$dbname,
           'ditag_analysis:s'  => \$ditag_analysis,
	  );

#connect
my $dbconn     = dbconnect($dbhost, $dbport, $dbname, $dbuser);
my $df_adaptor = $dbconn->get_DitagFeatureAdaptor;

#get regions to look at
if($slicename){
    push(@slices, $dbconn->get_SliceAdaptor->fetch_by_region("toplevel", $slicename));
}
else{
  foreach my $slice (@{get_all_slices($dbconn, "toplevel")}){
    push(@slices, $slice);
  }
}

foreach my $slice (@slices){

  print STDERR "checking region ".$slice->seq_region_name."\n";
  #fetch paired (ditags)
  my $dfs = $df_adaptor->fetch_pairs_by_Slice($slice, $ditagtype, $ditag_analysis);

  if(!$dfs or !(scalar @$dfs)){
    print STDERR "no ditags found.\n";
    next
  }

  #cluster
  print STDERR "Have ".(scalar @$dfs)." ditags.\n" if $verbose;
  my $clusters = cluster_things($dfs);
  print STDERR "Maincluster: ".scalar @$clusters."\n" if $verbose;
  #collect ditag_features in clusters
  foreach my $cluster (@$clusters){
    my @collected_features = ();
    foreach my $ditag_pair (@{$cluster->{'features'}}){
      my $ditag_id = $ditag_pair->{'ditag'};
      my $pair_id  = $ditag_pair->{'pair_id'};
      my @ditag_features = @{$df_adaptor->fetch_all_by_ditagID($ditag_id, $pair_id)};
      push(@collected_features, @ditag_features);
    }
    #sub-cluster ditag_features
    my $clustered_features = cluster_features(\@collected_features);
    print STDERR "Subcluster: ".scalar @$clustered_features."\n" if $verbose;
    foreach my $subcluster (@$clustered_features){
      #looking at indiv. ditag-features
      check_positions($subcluster);
    }

  }

}


sub check_positions {
  my ($df_cluster) = @_;

  my $minus_start = 0;
  my $plus_start  = 0;
  my $minus_end   = 0;
  my $plus_end    = 0;
  my @ditagids    = ();
  my $minstart    = 999999999;
  my $maxend      = 0;

  print STDERR "Checking Cluster: ".$df_cluster->{'start'}."-".$df_cluster->{'end'}."\n";
  foreach my $df (@{$df_cluster->{'features'}}){
    #check whether this ditag_feature is start or end of the ditag
    my $localstart = $df->start;
    my $localend   = $df->end;
    if($df->start < $df->end){
      $localstart = $df->start;
      $localend   = $df->end;
    }
    else{
      print STDERR "VERDREHT!\n";
      $localstart = $df->end;
      $localend   = $df->start;
    }
    my ($totalstart, $totalend, $strand) = $df->get_ditag_location;
    my ($ditagpair) = $df_adaptor->fetch_all_by_ditagID($df->ditag_id, $df->ditag_pair_id, $df->analysis->dbID);
    if(scalar @$ditagpair != 2){
      warn "funny number for ".$df->ditag_id.", ".$df->ditag_pair_id.", ".$df->analysis->dbID.", ".
	   scalar @$ditagpair."\n";
    }
    push(@ditagids, $ditagpair->[0]->dbID);
    push(@ditagids, $ditagpair->[1]->dbID);

    if($totalstart < $minstart){ $minstart = $totalstart }
    if($totalend   > $maxend){   $maxend = $totalend }

    if($localstart == $totalstart){
      if($strand == -1){
	$minus_start++;
	print "ms." if $verbose;
      }
      else{
	$plus_start++;
	print "ps." if $verbose;
      }
    }
    if($localend == $totalend){
      if($strand == -1){
	$minus_end++;
	print "me." if $verbose;
      }
      else{
	$plus_end++;
	print "pe." if $verbose;
      }
    }

    print "\t".$df->ditag_id.":  ".$df->start." - ".$df->end.", ".$df->strand." (".$totalstart."-".$totalend.")\n"
      if $verbose;
  }

  if((($minus_start > 4) && ($plus_end > 4)) || (($minus_end > 4) && ($plus_start > 4)) ||
     (($minus_start > 4) && ($minus_end > 4)) || (($plus_start > 4) && ($plus_end > 4))){
    print STDERR "evil pos: ".$minstart."-".$maxend." [".$df_cluster->{'start'}."-".$df_cluster->{'end'}."] ".(scalar @ditagids)."\n";
    print STDERR "UPDATE ditag_feature SET analysis_id=1001 WHERE ditag_feature_id IN(".(join(", ", @ditagids)).");\n";

  }
  print "\n";

}

