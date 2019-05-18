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
=head1 NAME

  fsk::Default

=head1 SYNOPSIS

  use Bio::EnsEMBL::Analysis::Tools::Default;
  Default->dbconnect($dbhost, $dbport, $dbname, $dbuser, $dbpass, $species);

=head1 DESCRIPTION

  various functions used regulary in other scripts.

=head1 CONTACT

  fsk@sanger.ac.uk


=cut

package Bio::EnsEMBL::Analysis::Tools::Default;

use warnings ;
use Exporter;
use vars qw(@ISA @EXPORT);
our @ISA    = ("Exporter");
our @EXPORT = qw( dbconnect get_all_slices cluster_features cluster_things);
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my $verbosity = 0;

=head2 dbconnect

 Arg[1]      : dbhost
 Arg[2]      : (optional) dbport
 Arg[3]      : dbname
 Arg[4]      : dbuser
 Arg[5]      : (optional) dbpass
 Arg[6]      : (optional) species
 Arg[7]      : (optional) dnadbname
 Arg[8]      : (optional) dnadbhost
 Arg[8]      : (optional) dnadbport
 Description : create connection to ensembl database
               use env variables if none given
 Returntype  : Bio::EnsEMBL::DBSQL::DBAdaptor

=cut

sub dbconnect {
  my ($dbhost, $dbport, $dbname, $dbuser, $dbpass, $species, $dnadbname, $dnadbhost, $dnadbport) = @_;

  my $db;
  #try to use env vars or default vals
  if(!$dbhost)   { $dbhost=$ENV{"DBHOST"} }
  if(!$dbport)   { $dbport=$ENV{"DBPORT"} || 3306    }
  if(!$dbuser)   { $dbuser=$ENV{"DBUSER"} || "ensro" }
  if(!$dbpass)   { $dbpass=$ENV{"DBPASS"} || ""      }
  if(!$dbname)   { $dbpass=$ENV{"DBNAME"} || die "\nno db name given.\n" }

  $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
					   -host    => $dbhost,
					   -user    => $dbuser,
					   -pass    => $dbpass,
					   -port    => $dbport,
					   -dbname  => $dbname,
					   #-species => $species,
					  ) or die("cant connect to $dbname");
  if($dnadbname && $dnadbhost){
    my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor (
						    -host   => $dnadbhost,
						    -dbname => $dnadbname,
						    -port   => $dnadbport,
						    -user   => $dbuser,
						    -pass   => $dbpass,
						   ) or die "can t connect to database $dnadbname.";
    $db->dnadb($dnadb);
  }
  print "connected to $dbname.\n" if($verbosity);

  return $db;
}




=head2 get_all_slices

 Arg[1]    : dbconnection
 Arg[2]    : coordssytem
 Arg[3]    : (optinal) coordsystem version
 Arg[4]    : (optinal) flag to remove NT slices
 Function  : fetch slices (chromosomes) from the database, sort them
 Returntype: array ref of Bio::EnsEMBL::Slice

=cut

sub get_all_slices {
  my ($dbObj, $coordsys, $version, $remove_NTs) = @_;

  my @slices = ();
  my $slice_adaptor = $dbObj->get_SliceAdaptor;

  foreach my $chr (@{$slice_adaptor->fetch_all($coordsys, $version, undef, 1)}){
    if($remove_NTs and (($chr->seq_region_name =~ /_NT_/) or ($chr->seq_region_name =~/_RANDOM_/))){ next; }
    #if($chr->seq_region_name =~ /MT/){ next; }
    print $chr->seq_region_name.", ".$chr->start." - ".$chr->end."\n" if($verbosity);
    push @slices, $chr;
  }
  @slices = sort sortbychrnum @slices;

  return \@slices;
}





sub sortbychrnum {

  my @awords = split /_/,$a->seq_region_name;
  my @bwords = split /_/,$b->seq_region_name;

  my $anum = $awords[0];
  my $bnum = $bwords[0];

  $anum =~ s/chr//;
  $bnum =~ s/chr//;

  if ($anum !~ /^[0-9]*$/) {
    if ($bnum !~ /^[0-9]*$/) {
      return $anum cmp $bnum;
    } else {
      return 1;
    }
  }
  if ($bnum !~ /^[0-9]*$/) {
    return -1;
  }

  if ($anum <=> $bnum) {
    return $anum <=> $bnum;
  } else {
    if ($#awords == 0) {
      return -1;
    } elsif ($#bwords == 0) {
      return 1;
    } else {
	return $awords[1] cmp $bwords[1];
      }
  }
}



=head2 cluster_features

 Arg[1]    : array ref of features with start & end functions
 Function  : cluster features into overlapping groups
 Returntype: array ref of cluster objects with $start, $end and @features

 # UPDATED CODE TO USE ABSOLUTE COORDINATES!

=cut

sub cluster_features {
  my ($features, $verbose, $use_strand) = @_;

  my @clusters      = ();
  my @finalclusters = ();

 FEATURE:
  foreach my $feature (@$features) {
    print $feature->seq_region_start."-".$feature->seq_region_end."\n" if $verbose;
  CLUSTER:
    foreach my $cluster (@clusters) {
      next CLUSTER if($use_strand and ($feature->strand ne $cluster->{strand}));
      if (($feature->seq_region_end > $cluster->{start}) and ($feature->seq_region_start < $cluster->{end})) {
	#add to existing cluster
	print "adding to ".$cluster->{start}." - ".$cluster->{end}."\n" if $verbose;
	if($feature->seq_region_start < $cluster->{start}){
	  $cluster->{start} = $feature->seq_region_start;
	}
	if($feature->seq_region_end > $cluster->{end}){
	  $cluster->{end} = $feature->seq_region_end;
	}
	push (@{$cluster->{features}}, $feature);
	$cluster->{count} += 1;
	next FEATURE;
      }
    }
    #create new cluster
    print "creating new cluster.\n" if $verbose;
    my %newcluster;
    $newcluster{start}    = $feature->seq_region_start;
    $newcluster{end}      = $feature->seq_region_end;
    $newcluster{strand}   = $feature->strand;
    $newcluster{features} = [$feature];
    $newcluster{count} = 1;
    push(@clusters, \%newcluster);
  }
  print "Have ".(scalar @clusters)." clusters. Cluster 1: ".$clusters[0]->{start}."-".$clusters[0]->{end}."\n" if $verbose;


  return(\@clusters);
  #rest is not needed anymore!?

  #join overlapping clusters
  @clusters = sort {$a->{start} <=> $b->{start}} @clusters;

  for(my $i = 0; $i < (scalar @clusters); $i++) {
    #print $i.": ".$clusters[$i]->{start}."-".$clusters[$i]->{end}."\n";
    if(exists($clusters[$i+1]) && ($clusters[$i]->{end} > $clusters[$i+1]->{start})){
      $clusters[$i]->{end} = $clusters[$i+1]->{end};
      push (@{$clusters[$i]->{features}}, @{$clusters[$i+1]->{features}});
      $clusters[$i]->{count} += $clusters[$i+1]->{count};
      push(@finalclusters, $clusters[$i]);
      $i++;
      print "joining cluster with next: ".$clusters[$i]->{end}.">".$clusters[$i+1]->{start}."\n"; # if $verbose;
      next;
    }
    push(@finalclusters, $clusters[$i]);
  }

  @clusters = ();

  return(\@finalclusters);
}


=head2 cluster_things

 Arg[1]    : array ref of features with start & end hash entries
 Function  : cluster features into overlapping groups
 Returntype: array ref of cluster objects with $start, $end and @features

=cut

sub cluster_things {
  my ($features, $dont_save_features) = @_;

  my @clusters      = ();
  my @finalclusters = ();

 FEATURE:
  foreach my $feature (sort  { $a->{start} == $b->{start} ?  ( $b->{end} <=> $a->{end} ) : ( $a->{start} <=> $b->{start} ) } @$features) {
    #print $feature->{start}."-".$feature->{end}."\n";
  CLUSTER:
    foreach my $cluster (@clusters) {
      if (($feature->{end} > $cluster->{start} and $feature->{start} < $cluster->{end})) {
	#add to existing cluster
	#print "adding to ".$cluster->{start}." - ".$cluster->{end}."\n";
	if($feature->{start} < $cluster->{start}){
	  $cluster->{start} = $feature->{start};
	}
	if($feature->{end} > $cluster->{end}){
	  $cluster->{end} = $feature->{end};
	}
	unless($dont_save_features){
	  push (@{$cluster->{features}}, $feature);
	}
	next FEATURE;
      }
    }

    #create new cluster
    #print "creating new cluster.\n";
    my %newcluster;
    $newcluster{start}    = $feature->{start};
    $newcluster{end}      = $feature->{end};
    unless($dont_save_features){
      $newcluster{features} = [$feature];
    }
    push(@clusters, \%newcluster);
  }
  #print "Have ".(scalar @clusters)." clusters. Cluster 1: ".$clusters[0]->{start}."-".$clusters[0]->{end}."\n";

  #join overlapping clusters
  @clusters = sort {$a->{start} <=> $b->{start}} @clusters;

  for(my $i = 0; $i < (scalar @clusters); $i++) {
    #print $i.": ".$clusters[$i]->{start}."-".$clusters[$i]->{end}." (".
    #  ($clusters[$i]->{end} - $clusters[$i]->{start}).")\n";
    if(exists($clusters[$i+1]) and ($clusters[$i]->{end} > $clusters[$i+1]->{start})){
      $clusters[$i]->{end} = $clusters[$i+1]->{end};
      unless($dont_save_features){
	push (@{$clusters[$i]->{features}}, @{$clusters[$i+1]->{features}});
      }
      push(@finalclusters, $clusters[$i]);
      $i++;
      #print "joining with next: ".$clusters[$i]->{end}.">".$clusters[$i+1]->{start}."\n";
      next;
    }
    push(@finalclusters, $clusters[$i]);
  }

  return(\@finalclusters);
}



1;
