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
use Bio::EnsEMBL::Mapper;

$| = 1;

my ($dbname, 
    $dbuser,
    $dbhost,
    $dbport,
    $verbose,
    $agp_outfile, $agp_outfh,
    $gene_outfile, $gene_outfh,
    $log_outfile, $log_outfh,
    @agp_files,
    @gene_files);

$dbport = 3306;
GetOptions('agp=s@' => \@agp_files,
            'genes=s@' => \@gene_files,
            'outagp=s' => \$agp_outfile,
            'outgenes=s' => \$gene_outfile,
            'outlog=s' => \$log_outfile,
            'verbose' => \$verbose,
            'dbname|db|D=s' => \$dbname,
            'dbhost|host|h=s' => \$dbhost,
            'dbuser|user|u=s' => \$dbuser,
            'dbport|port|P=s' => \$dbport);

$| = 1;

my $DB = Bio::EnsEMBL::DBSQL::DBAdaptor->new(-dbname => $dbname,
                                             -host   => $dbhost,
                                             -user   => $dbuser,
                                             -port   => $dbport);
my $FAKE_SCAFFOLD_PREFIX = 'DUMMYSCAFFOLD_';
my $GENE_SCAFFOLD_PREFIX = 'GeneScaffold_';
my $fake_scaffold_count = 1;
my $GENE_SCAFFOLD_PADDING = 100;

if (defined $agp_outfile) {
  open $agp_outfh, ">$agp_outfile" 
     or die "Could not open $agp_outfile for writing\n";
} else {
  $agp_outfh = \*STDOUT;
}
if (defined $gene_outfile) {
  open $gene_outfh, ">$gene_outfile" 
      or die "Could not open $gene_outfile for writing\n";
} else {
  $gene_outfh = \*STDOUT;
}
if (defined $log_outfile) {
  open $log_outfh, ">$log_outfile" 
      or die "Could not open $log_outfile for writing\n";
} else {
  $log_outfh = \*STDERR;
}


$verbose and print STDERR "Indexing AGP files...\n";
my ($gs_index, $s_index) = &index_gene_scaffold_files(@agp_files);

$verbose and print STDERR "Indexing Annotation files...\n";
&index_gene_annotation_files(@gene_files);

my @clusters;
if (@ARGV) {
  @clusters = ([@ARGV]);
} else {
  $verbose and print STDERR "Clustering gene scaffolds...\n";
  @clusters = &single_linkage_cluster($gs_index, $s_index);
}

my $out_gene_scaffold_count = 1;

for(my $cl_cnt=0; $cl_cnt < @clusters; $cl_cnt++) {
  my $cluster = $clusters[$cl_cnt];

  $verbose and print STDERR "Doing $cl_cnt Cluster @$cluster\n";
  print $log_outfh "CLUSTER: ", join("|", @$cluster), "\n";

  my $gene_scaffolds = {};

  &fetch_gene_scaffold_entries($gene_scaffolds,   @$cluster);  
  &fetch_gene_annotation_entries($gene_scaffolds, @$cluster); 

  # transcript filtering could have left gene scaffold fragments
  # that have no exons; remove these
  #&remove_unused_components_from_gene_scaffolds($gene_scaffolds);

  # add in the gaps we have "plugged" as fake components
  &identify_plugged_gaps($gene_scaffolds);

  # we need the slices and their seq-level component structure for 
  # assessing scaffold splits
  my ($target_slices, $target_slice_maps) = 
      &fetch_slices_and_component_maps($DB, $gene_scaffolds);

  # make list of "non trivial" gene scaffolds, which are those
  # that involve more than one component

  my (%simple_gene_scaffolds, %complex_gene_scaffolds);
  foreach my $gs_id (keys %$gene_scaffolds) {
    if (&is_simple_gene_scaffold($gs_id, $gene_scaffolds)) {
      $simple_gene_scaffolds{$gs_id} = $gene_scaffolds->{$gs_id};
    } else {
      $complex_gene_scaffolds{$gs_id} = $gene_scaffolds->{$gs_id};
    }
  }

  my $chains = &make_gene_scaffold_chains(\%complex_gene_scaffolds);

  my $extended_chains = 
      &extend_gene_scaffold_components($chains,
                                       \%complex_gene_scaffolds,
                                       \%simple_gene_scaffolds,
                                       $target_slices,
                                       $target_slice_maps);

  # make a map to/from Gene scaffold coords and component coords
  my $orig_map = Bio::EnsEMBL::Mapper->new('scaffold',
                                           'genescaffold');
    
  my $new_map     = Bio::EnsEMBL::Mapper->new('scaffold',
                                              'genescaffold');
  
  foreach my $gs_id (keys %$gene_scaffolds) {
    my $gs = $gene_scaffolds->{$gs_id};
    foreach my $seg (@{$gs->{components}}) {
      $orig_map->add_map_coordinates($seg->to->id,
                                     $seg->to->start,
                                     $seg->to->end,
                                     $seg->to->strand,
                                     $gs_id,
                                     $seg->from->start,
                                     $seg->from->end,
                                     );
    }
  }

  foreach my $ch (@$extended_chains) {
    my @gs_ids = @{$ch->{members}};

    my $new_gs_name = $GENE_SCAFFOLD_PREFIX . $out_gene_scaffold_count++;
    
    print $agp_outfh "##-AGP for $new_gs_name [@gs_ids]\n";

    my @coords   = @{$ch->{merged_coords}};
    
    my $last_end = 0;
    my $line_count = 1;
    for(my $i=0; $i < @coords; $i++) {
      my $c = $coords[$i];
      
      my $gs_start = $last_end + 1;
      $last_end    = $gs_start + $c->length - 1;
      
      $new_map->add_map_coordinates($c->id,
                                    $c->start,
                                    $c->end,
                                    $c->strand,
                                    $new_gs_name,
                                    $gs_start,
                                    $last_end,
                                    );
      
      if ($c->id =~ /^$FAKE_SCAFFOLD_PREFIX/) {
        printf($agp_outfh "%s\t%d\t%d\t%d\tN\t%d\n", 
               $new_gs_name, 
               $gs_start,
               $last_end,
               $line_count++,
               $last_end - $gs_start + 1);       
      } else {
        printf($agp_outfh "%s\t%d\t%d\t%d\tW\t%s\t%d\t%d\t%s\n", 
               $new_gs_name, 
               $gs_start, 
               $last_end,
               $line_count++,
               $c->id,
               $c->start,
               $c->end,
               $c->strand > 0 ? "+" : "-");        
      }
      
      unless ($i == @coords - 1) {
        $gs_start = $last_end + 1;
        $last_end = $gs_start + $GENE_SCAFFOLD_PADDING - 1;
        
        printf($agp_outfh "%s\t%d\t%d\t%d\tN\t%d\n", 
               $new_gs_name, 
               $gs_start,
               $last_end,
               $line_count++,
               $last_end - $gs_start + 1);
        
      }
    }
    
    printf $gene_outfh "# GENE for $new_gs_name [@gs_ids]\n";

    foreach my $gs_id (@gs_ids) {
      foreach my $line (@{$gene_scaffolds->{$gs_id}->{annotation}}) {
        if (ref($line) ne "ARRAY") {
          print $gene_outfh $line;
          next;
        } 
        
        my ($loc_in_orig) = $orig_map->map_coordinates($line->[0],
                                                       $line->[3],
                                                       $line->[4],
                                                       $line->[5],
                                                       'genescaffold');
        
        my ($loc_in_new) = $new_map->map_coordinates($loc_in_orig->id,
                                                     $loc_in_orig->start,
                                                     $loc_in_orig->end,
                                                     $loc_in_orig->strand,
                                                     'scaffold');
        
        $line->[0] = $new_gs_name;
        $line->[3] = $loc_in_new->start;
        $line->[4] = $loc_in_new->end;
        $line->[5] = $loc_in_new->strand;
        print $gene_outfh join("\t", @$line);
      }
    }
  }
 
  # we need to deal with the simple gene scaffolds which refer
  # to a component that has ended up in this gene scaffold
  foreach my $gs_id (keys %simple_gene_scaffolds) {

    my @cmps = @{$simple_gene_scaffolds{$gs_id}->{components}};

    my %result_gs;
    
    foreach my $cmp (@cmps) {
      my ($loc_in_orig) = $orig_map->map_coordinates($gs_id,
                                                     $cmp->from->start,
                                                     $cmp->from->end,
                                                     1,
                                                     'genescaffold');
      
      my ($loc_in_new) = $new_map->map_coordinates($loc_in_orig->id,
                                                   $loc_in_orig->start,
                                                   $loc_in_orig->end,
                                                   1,
                                                   'scaffold');
      
      if ($loc_in_new->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
        $result_gs{$loc_in_new->id}++;
      }
    }
    
    if (scalar(keys %result_gs) == 1) {
      my ($other_gs) = keys %result_gs;

      foreach my $line (@{$simple_gene_scaffolds{$gs_id}->{annotation}}) {
        if (ref($line) ne "ARRAY") {
          print $gene_outfh $line;
          next;
        } 
        
        my ($loc_in_orig) = $orig_map->map_coordinates($line->[0],
                                                       $line->[3],
                                                       $line->[4],
                                                       $line->[5],
                                                       'genescaffold');
        
        my ($loc_in_new) = $new_map->map_coordinates($loc_in_orig->id,
                                                     $loc_in_orig->start,
                                                     $loc_in_orig->end,
                                                     $loc_in_orig->strand,
                                                     'scaffold');
        
        $line->[0] = $loc_in_new->id;
        $line->[3] = $loc_in_new->start;
        $line->[4] = $loc_in_new->end;
        $line->[5] = $loc_in_new->strand;
        print $gene_outfh join("\t", @$line);
      }
      
      delete $simple_gene_scaffolds{$gs_id};

      print $log_outfh "WRITE: Accommodated annotation for $gs_id into $other_gs\n";

    } elsif (scalar(keys %result_gs)) {
      delete $simple_gene_scaffolds{$gs_id};

      print($log_outfh "WRITE: Simple gs $gs_id is split between >1 complex ones so rejecting : " 
            . join(",", keys %result_gs) 
            ."\n");    
    }

  }

  # for all of the remaining simple gene scaffolds, the annotation 
  # can be written directly without need for agp

  foreach my $gs_id (keys %simple_gene_scaffolds) {
    printf $gene_outfh "# GENE output for %s\n", $gs_id; 
      
    foreach my $line (@{$simple_gene_scaffolds{$gs_id}->{annotation}}) {
      if (ref($line) ne "ARRAY") {
        print $gene_outfh $line;
        next;
      } 
      
      my ($mapped_reg) = $orig_map->map_coordinates($line->[0],
                                                    $line->[3],
                                                    $line->[4],
                                                    1,
                                                    'genescaffold');
      $line->[0] = $mapped_reg->id;
      $line->[3] = $mapped_reg->start;
      $line->[4] = $mapped_reg->end; 
      
      if ($simple_gene_scaffolds{$gs_id}->{components}->[0]->to->strand < 0) {
        $line->[5] *= -1;
      }
      print $gene_outfh join("\t", @$line);
    }

    print $log_outfh "WRITE: Accommodated annotation for $gs_id without gene scaffold\n";
  }
}



#################################################################
# single_linkage_cluster
#
#   Single-linkage cluster the gene_scaffolds based on common 
#   scaffold
#################################################################
sub single_linkage_cluster {
  my ($gs_index, $s_index) = @_;
  
  my %clusters;
  
  foreach my $gs_id (keys %{$gs_index}) {
    # find out which scaffold are implicated in this gene scaf
    # find other gene scafs that share these scaffolds
    # form a cluster, and make all gene scaf ids, and their cluster 
    # members, 
    #  point to the same cluster

    my @other_gs_ids;
    
    foreach my $s_id (keys %{$gs_index->{$gs_id}}) {
      # look for other gene scaffold that are involved in this scaffold
      foreach my $other_gs_id (keys %{$s_index->{$s_id}}) {
        push @other_gs_ids, $other_gs_id if $other_gs_id ne $gs_id;
      }
    }
    
    my (%new_cluster);
    $new_cluster{$gs_id} = 1;

    # merge with existing clusters
    foreach my $other_gs_id (@other_gs_ids) {
      if (exists $clusters{$other_gs_id}) {
        foreach my $k (keys %{$clusters{$other_gs_id}}) {
          $new_cluster{$k} = 1;
        }
      }
      $new_cluster{$other_gs_id} = 1;
    }
    
    # make all members of the revised cluster point to it
    foreach my $gs_mem (keys %new_cluster) {
      $clusters{$gs_mem} = \%new_cluster;
    }
  }
  
  # finally, return a list of distinct clusters
  my %unique_clusters;
  foreach my $clust (values %clusters) {
    if (not exists $unique_clusters{$clust}) {
      my @members = keys %$clust;
      $unique_clusters{$clust} = \@members;
    }
  }

  my @clusts = values %unique_clusters;
  @clusts = sort { $a->[0] cmp $b->[0] } @clusts;

  return @clusts;
}



#################################################################
# identify_plugged_gaps
#
#   looks for exons that lie on gene_scaffold gaps
#################################################################
sub identify_plugged_gaps {
  my ($gene_scaffolds) = @_;

  my %plugged_gaps;

  foreach my $gs_id (keys %{$gene_scaffolds}) {
    my @lines;
    foreach my $line (@{$gene_scaffolds->{$gs_id}->{annotation}}) {
      if (ref($line) eq "ARRAY") {
        push @lines, $line;
      }
    }

    my @comps = sort {
      $a->from->start <=> $b->from->start;
    } @{$gene_scaffolds->{$gs_id}->{components}};

    $plugged_gaps{$gs_id} = [];

    my @gaps;

    foreach my $l (@lines) {
      my ($type, $start, $end) = ($l->[2], $l->[3], $l->[4]);
      
      next if $type ne 'Exon';

      # does this exon lie in a gap? 
      my $in_gap = 1;
      foreach my $seg (@comps) {
        if ($start >= $seg->from->start and $end <= $seg->from->end) {
          $in_gap = 0;
          last;
        }
      }
      
      if ($in_gap) {
        push @gaps, {
          start => $start, 
          end   => $end,
        };
      }
    }

    foreach my $gap (sort { $a->{start} <=> $b->{start} } @gaps) {
      if (@{$plugged_gaps{$gs_id}} and 
          $plugged_gaps{$gs_id}->[-1]->{end} >= $gap->{start}) {
        if ($gap->{end} > $plugged_gaps{$gs_id}->[-1]->{end}) {
          $plugged_gaps{$gs_id}->[-1]->{end} = $gap->{end};
        }
      } else {
        push @{$plugged_gaps{$gs_id}}, $gap;
      }
    }
  }
 
  foreach my $gs_id (keys %plugged_gaps) {

    foreach my $reg (@{$plugged_gaps{$gs_id}}) {

      my $fake_name = &get_unique_fake_scaffold_id;
      my $fake_start = 1;
      my $fake_end = $reg->{end} - $reg->{start} + 1;

      my $fake_q = Bio::EnsEMBL::Mapper::Coordinate->new($gs_id, 
                                                         $reg->{start}, 
                                                         $reg->{end},
                                                         1);
      my $fake_t = Bio::EnsEMBL::Mapper::Coordinate->new($fake_name, 
                                                         $fake_start, 
                                                         $fake_end,
                                                         1);
      my $pair   = Bio::EnsEMBL::Mapper::Pair->new($fake_q, 
                                                   $fake_t);
      
      push @{$gene_scaffolds->{$gs_id}->{components}}, $pair;
    }

    $gene_scaffolds->{$gs_id}->{components} = 
        [
         sort { 
           $a->from->start <=> $b->from->start; 
         } @{$gene_scaffolds->{$gs_id}->{components}} 
        ];
  }
}


#################################################################
# make_gene_scaffold_chains
#
#   returns a list of chains, each element of which is a
#   list of AGP references that can be glued together, in
#   order, to make a larger AGP. Each of the given AGPs
#   occurs in exactly one of the returned chains
#################################################################
sub make_gene_scaffold_chains {
  my ($gene_scaffolds) = @_; 

  my $go_further = 1;

  $go_further = 0 
      if scalar(keys %{$gene_scaffolds}) <= 1;

  my (@all_chains, %sub_clusters, @shared_sids);

  if ($go_further) {
    # form chains based on component composition

    foreach my $gs_id (keys %{$gene_scaffolds}) {
      foreach my $comp (@{$gene_scaffolds->{$gs_id}->{components}}) {
        push @{$sub_clusters{$comp->to->id}->{$gs_id}}, $comp;
      }
    }

    foreach my $s_id (keys %sub_clusters) {
      if (scalar(keys %{$sub_clusters{$s_id}}) > 1) {
        push @shared_sids, $s_id;
      }
    }
    
    if (scalar(@shared_sids) > 1) {
      print $log_outfh "CHAINING: cluster has more than one shared scaffold: @shared_sids\n";
    } elsif (not @shared_sids) {
      print $log_outfh "CHAINING: odd; cluster had no shared scaffolds\n";
      $go_further = 0;
    }
  }

  if ($go_further) {
    # check that the shared component is self-consistent in each 
    # implicated gene scaffold

    foreach my $shared_sid (@shared_sids) {
      my (@extents, @consistent_extents, @chains);

    GS:
      foreach my $gs_id (keys %{$sub_clusters{$shared_sid}}) {

        my @comps = @{$sub_clusters{$shared_sid}->{$gs_id}};
        
        my ($tstart, $tend, $strand);
        
        foreach my $comp (@comps) {
          if (not defined $strand) {
            $tstart = $comp->to->start;
            $tend   = $comp->to->end;
            $strand = $comp->to->strand;
          } else {
            # need to check for consistency
            if ($comp->to->strand == $strand) {
              if ($strand > 0 and
                  $comp->to->start > $tend) {
                $tend = $comp->to->end;
              } elsif ($strand < 0 and
                       $comp->to->end < $tstart) {
                $tstart = $comp->to->start;
              } else {
                next GS;
              }
            }
          }
        }
        
        push @extents, {
          gs_id   => $gs_id,
          qid     => $gene_scaffolds->{$gs_id}->{ref_name},
          qstart  => $gene_scaffolds->{$gs_id}->{ref_start},
          qend    => $gene_scaffolds->{$gs_id}->{ref_end},
          tstart  => $tstart,
          tend    => $tend,
          tstrand => $strand,
        };
      }
    
      # generate all possible selections of consistent chains from the list
    
      if (scalar(@extents) > 1) {
        my @sub_lists = _generate_sub_lists(@extents);
        
        foreach my $sl (@sub_lists) {
          my @these_extents = @$sl;
          
          # check if this is a consistent sub-list;
          my $consistent = 1;
          for(my $i=1; $i < @these_extents; $i++) {
            my $this = $these_extents[$i];
            my $last = $these_extents[$i-1];
            
            if ($this->{qid} eq $last->{qid} and
                $this->{qstart} > $last->{qend} and
                $this->{tstrand} == $last->{tstrand} and
                (($this->{tstrand} > 0 and $this->{tstart} > $last->{tend}) or
                 ($this->{tstrand} < 0 and $this->{tend} < $last->{tstart}))) {
              # consistent
            } else {
              $consistent = 0;
              last;
            }
          }
          if ($consistent) {
            push @consistent_extents, $sl;
          }
        }
      }
    
      @consistent_extents = sort { scalar(@$b) <=> scalar(@$b) } @consistent_extents;
      
      my %gs_in_chains;
      
      while (@consistent_extents) {
        my $first = shift @consistent_extents;
        my @gs_ids = map { $_->{gs_id} } @$first;
        
        if (not grep { exists($gs_in_chains{$_}) } @gs_ids) {
          push @chains, \@gs_ids;
          map { $gs_in_chains{$_} = 1 } @gs_ids;
        }
      }
      
      foreach my $ch (@chains) {
        print $log_outfh "CHAINING: formed chain " . join(" ", @$ch) . "\n";
        push @all_chains, $ch;
      }
    }
  }
  
  # attempt to merge consistent chains together
  @all_chains = sort { scalar(@$b) <=> scalar(@$a) } @all_chains;

  if (scalar(@all_chains) > 1) {
    print $log_outfh "CHAINING: Attempting to consolidate ", scalar(@all_chains), " chains\n";
  }

  my @merged_chains;

  foreach my $c (@all_chains) {
    my $merged = 0;

    OC: foreach my $o_c (@merged_chains) {      
      my ($i, $j);
      for($i=0; $i < @$o_c; $i++) {
        if ($o_c->[$i] eq $c->[0]) {
          my $mismatch = 0;
          for($j=0; $j < @$c and $i < @$o_c ; $j++, $i++) {
            if ($c->[$j] ne $o_c->[$i]) {
              $mismatch = 1;
            }
          }
          if (not $mismatch) {
            for(my $k=$j; $k < @$c; $k++) {
              push @$o_c, $c->[$k];
            }
            $merged = 1;
            print $log_outfh "CHAINING: merging ", join(":", @$c), " to ", join(":", @$o_c), "\n";
            last OC;
          } else {
            last;
          }
        } 
      }
      for($i=scalar(@$o_c)-1; $i>=0; $i--) {
        if ($o_c->[$i] eq $c->[-1]) {
          my $mismatch = 0;
          for($j=scalar(@$c)-1; $j>=0 and $i>=0; $j--, $i--) {
            if ($c->[$j] ne $o_c->[$i]) {
              $mismatch = 1;
            }
          }
          if (not $mismatch) {
            for(my $k=$j; $k >= 0; $k--) {
              unshift @$o_c, $c->[$k];
            }
            $merged = 1;
            print $log_outfh "CHAINING: merging ", join(":", @$c), " to ", join(":", @$o_c), "\n";
            last OC;
          } else {
            last;
          }
        }
      }
    }
    if (not $merged) {
      push @merged_chains, $c;
    }
  }

  # any remaining chains that share elements are removed
  my %chains_by_el;
  foreach my $c (@merged_chains) {
    map { push @{$chains_by_el{$_}}, $c } @$c;
  }

  my %delete;
  foreach my $el (keys %chains_by_el) {
    if (@{$chains_by_el{$el}} > 1) {
      foreach my $cref (@{$chains_by_el{$el}}) {
        map { $delete{$_} = 1 } @$cref;
        @merged_chains = grep { $_ ne $cref } @merged_chains;
      }
      print $log_outfh "CHAINING: $el occurs inconsistently in more than one chain; removing\n";
    }
  }

  map { delete $chains_by_el{$_} } keys %delete;

  # finally, we have to make singleton chains for all 
  # gene-scaffolds not in chains
  
  foreach my $gs_id (keys %$gene_scaffolds) {
    if (not exists $chains_by_el{$gs_id}) {
      push @merged_chains, [$gs_id];
    }
  }

  foreach my $c (@merged_chains) {
    print $log_outfh "CHAINING: final chain " . join(":", @$c), "\n";
  }
  
  return \@merged_chains;
}

sub _generate_sub_lists {
  my @list = @_;

  if (scalar(@list) == 2) {
    return ([$list[0], $list[1]],[$list[1],$list[0]]);
  } else {
    my ($first, @rest) = @list;

    my @s_lists = &_generate_sub_lists(@rest);

    if (scalar(@s_lists) > 50_000) {
      print "!!!! Bailing out of combinatorial explosion !!!!!\n";
      return @s_lists;
    }

    my @new_s_lists;
    foreach my $sl (@s_lists) {
      my @this_list = @$sl;
      my $size = scalar(@this_list);

      push @new_s_lists, [$first, @this_list];
      for(my $i = 0; $i < $size - 1; $i++) {
        my @new_list = (@this_list[0..$i], $first, @this_list[$i+1..$size-1]);

        push @new_s_lists, \@new_list;
      }

      push @new_s_lists, [@this_list, $first];
    }
    push @s_lists, @new_s_lists;

    return @s_lists;
  }
}

#################################################################
# sort_gene_scaffolds_by_common_component
#
#################################################################
sub sort_gene_scaffolds_by_common_component {
  my ($gene_scaffolds, $s_id, $gs_ids) = @_;

  # these gene scaffolds can be oriented if the extents
  # of the components that they have in common are in
  # orientation agreement and are non-overlapping
  my @gs_ids = @$gs_ids;

  my @summaries;
  foreach my $gs_id (@$gs_ids) {
    my $extents = 
        &get_gene_scaffold_component_extents($gene_scaffolds->{$gs_id});

    if (exists($extents->{$s_id}->{ori}->{1}) and 
        exists($extents->{$s_id}->{ori}->{-1})) {
      # this AGP has fragments from a single scaffold in
      # both orientation. We may be use a consensus 
      # orientation, but for now, give up. 
      
      printf($log_outfh "SORT: Cannot order [@gs_ids] on $s_id due to internal ori disagreement in %s\n", 
             $gs_id);
      return 0;
    }

    my ($ori) = keys %{$extents->{$s_id}->{ori}};

    push @summaries, {
      gs_id => $gs_id,
      start => $extents->{$s_id}->{start},
      end   => $extents->{$s_id}->{end},
      ori   => $ori,
    };
  }

  for(my $i=0; $i < @summaries; $i++) {
    for(my $j=0; $j < $i; $j++) {
      if ($summaries[$i]->{start} <= $summaries[$j]->{end} and 
          $summaries[$i]->{end}   >= $summaries[$j]->{start}) {
        # overlap
        printf($log_outfh "SORT: Cannot order [@gs_ids] on $s_id due to overlap in %s and %s\n", 
               $summaries[$i]->{gs_id}, $summaries[$j]->{gs_id});
        return 0;
      }
      if ($summaries[$i]->{ori} != $summaries[$j]->{ori}) {
        printf($log_outfh "SORT: Cannot order [@gs_ids] on $s_id due to ori disagreement between %s and %s\n", 
               $summaries[$i]->{gs_id}, $summaries[$j]->{gs_id});
        return 0;
      }
    }
  }

  if ($summaries[0]->{ori} < 0) {
    @summaries = sort { $b->{start} <=> $a->{start} } @summaries;
  } else {
    @summaries = sort { $a->{start} <=> $b->{start} } @summaries;
  }

  @$gs_ids = map { $_->{gs_id} } @summaries;
  return 1;
}

#################################################################
# extend_gene_scaffold_components
#
#################################################################
sub extend_gene_scaffold_components {
  my ($chains,
      $gene_scaffolds, 
      $simple_gene_scaffolds,
      $target_slices, 
      $target_slice_maps) = @_;

  
  # the annotation for the simple gene scaffolds can be
  # placed on the appropriate part of a a complex gene-scaffold
  # (or directly on the scaffold), but we must ensure that
  # do no split a scaffold in the middle of one of these
  # regions
  my %non_break_regions;
  foreach my $gs_id (keys %$simple_gene_scaffolds) {
    my ($min_start, $max_end, $id);
    foreach my $cmp (@{$simple_gene_scaffolds->{$gs_id}->{components}}) {
      $id = $cmp->to->id if not defined $id;
      $min_start = $cmp->to->start 
          if not defined $min_start or $cmp->to->start < $min_start;
      $max_end   = $cmp->to->end
          if not defined $max_end or $cmp->to->end > $max_end;
    }
    push @{$non_break_regions{$id}}, {
      start => $min_start,
      end   => $max_end,
    };
  }
  foreach my $sid (keys %non_break_regions) {
    my @segs = @{$non_break_regions{$sid}};
    @segs = sort { $a->{start} <=> $b->{start} } @segs;
    
    my @nov_segs;
    foreach my $seg (@segs) {
      if (not @nov_segs or $nov_segs[-1]->{end} < $seg->{start} - 1) {
        push @nov_segs, $seg;
      } elsif ($nov_segs[-1]->{end} < $seg->{end}) {
        $nov_segs[-1]->{end} = $seg->{end};
      }
    }
    
    $non_break_regions{$sid} = \@nov_segs;
  }
  
  
  my (%coords_by_scaffold);
  
  #
  # index-by-scaffold of all the agp fragments in this group
  #
  foreach my $chain (@$chains) {
    foreach my $gs_id (@$chain) {
      my $gs = $gene_scaffolds->{$gs_id};
      foreach my $pair (@{$gs->{components}}) {
        push @{$coords_by_scaffold{$pair->to->id}}, $pair->to;
      }
    }
  }
  #
  # join components together where possible
  #
  my @merged_chains;

  foreach my $chain (@$chains) {
    my @merged_coords;

    foreach my $gs_id (@$chain) {
    
      foreach my $pair (@{$gene_scaffolds->{$gs_id}->{components}}) {
        my $this_c = $pair->to;

        my $merged = 0;
        
        if (@merged_coords and
            $merged_coords[-1]->strand == $this_c->strand and
            $merged_coords[-1]->id eq $this_c->id) {
          
          if ($this_c->strand == 1) {
            if ($this_c->start > $merged_coords[-1]->end) {
              my $can_merge = 1;
              foreach my $o_c (@{$coords_by_scaffold{$this_c->id}}) {
                if ($o_c->start > $merged_coords[-1]->end and
                    $o_c->end   < $this_c->start) {
                  $can_merge = 0;
                  last;
                }
              }
              if ($can_merge) {
                $merged_coords[-1]->end($this_c->end);
                $merged = 1;
              }
            }
          } else {
            if ($this_c->end < $merged_coords[-1]->start) {
              my $can_merge = 1;
              foreach my $o_c (@{$coords_by_scaffold{$this_c->id}}) {
                if ($o_c->start > $this_c->end and
                    $o_c->end < $merged_coords[-1]->start) {
                  $can_merge = 0;
                  last;
                }
              }
              if ($can_merge) {
                $merged_coords[-1]->start($this_c->start);
                $merged = 1;
              }
            }
          }
        }
        
        if (not $merged) {
          push @merged_coords, Bio::EnsEMBL::Mapper::Coordinate
              ->new($this_c->id,
                    $this_c->start,
                    $this_c->end,
                    $this_c->strand);
        }
      }
    }

    push @merged_chains, {
      members => $chain,
      merged_coords => \@merged_coords,
    };
  }


  # now extend scaffold components so that for any scaffold involved
  # in a gene scaffold, every bp of it is accounted for somewhere

  %coords_by_scaffold = ();
  foreach my $ch (@merged_chains) {
    foreach my $c (@{$ch->{merged_coords}}) {
      if ($c->id !~ /^$FAKE_SCAFFOLD_PREFIX/) {
        push @{$coords_by_scaffold{$c->id}}, $c;
      }
    }
  }

  foreach my $sid (keys %coords_by_scaffold) {       
    my @sorted_els = sort {
      $a->start <=> $b->start;
    } @{$coords_by_scaffold{$sid}};

    my $tsl = $target_slices->{$sid};
    my $map = $target_slice_maps->{$sid};
    
    my (@contig_breaks, @scaffold_gaps, @non_break_regions);

    if (exists $non_break_regions{$sid}) {
      @non_break_regions = @{$non_break_regions{$sid}};
      @non_break_regions = sort { $a->{start} <=> $b->{start} } @non_break_regions;
    }
    
    my @bits = $map->map_coordinates($sid,
                                     1,
                                     $tsl->length,
                                     1,
                                     'scaffold');
    my $last_end = 0;
    for (my $i=0; $i < @bits; $i++) {
      my $bit = $bits[$i];
      
      my $this_start = $last_end + 1;
      my $this_end   = $this_start + $bit->length - 1;
      
      if ($bit->isa("Bio::EnsEMBL::Mapper::Gap")) {
        push @scaffold_gaps, $bit;
      } elsif ($i > 0 and 
               $bits[$i-1]->isa("Bio::EnsEMBL::Mapper::Coordinate")) {
        
        push @contig_breaks, Bio::EnsEMBL::Mapper::Gap->new($last_end,
                                                            $this_start);
      }
      
      $last_end = $bit->end;
    }

    
    $sorted_els[0]->start(1);
    for(my $i=1; $i<@sorted_els; $i++) {
      my $prev_el = $sorted_els[$i-1];
      my $this_el = $sorted_els[$i];
      
      my (@relevant_gaps, @relevant_breaks);

      foreach my $break (@scaffold_gaps) {
        if ($break->start <= $this_el->start and
            $break->end >= $prev_el->end) {
          # check that it does not overlap non-interruptable regions
          my $bad_break = 0;
          foreach my $reg (@non_break_regions) {
            if ($break->start >= $reg->{start} and
                $break->end   <= $reg->{end}) {
              $bad_break = 1; last;
            }
          }
          push @relevant_gaps, $break
              if not $bad_break;
        }
      }

      foreach my $break (@contig_breaks) {
        if ($break->start < $this_el->start and
            $break->end > $prev_el->end) {
          my $bad_break = 0;
          foreach my $reg (@non_break_regions) {
            if ($break->start >= $reg->{start} and
                $break->end   <= $reg->{end}) {
              $bad_break = 1; last;
            }
          }
          push @relevant_breaks, $break
              if not $bad_break;
        }
      }
      
      my $split_reg_start = $prev_el->end + 1;
      my $split_reg_end   = $this_el->start - 1;
      
      if (@relevant_gaps) {
        my $mid = int(scalar(@relevant_gaps) / 2);
        $split_reg_start = $relevant_gaps[$mid]->start;
        $split_reg_start = $prev_el->end + 1 if $split_reg_start < $prev_el->end + 1;
        
        $split_reg_end = $relevant_gaps[$mid]->end;
        $split_reg_end = $this_el->start - 1 if $split_reg_end > $this_el->start - 1;
        
        printf($log_outfh "EXTEND: Split in scaffold gap (%s/%d-%d -> %d-%d)\n", 
               $sid, 
               $prev_el->end, 
               $this_el->start,
               $split_reg_start,
               $split_reg_end);
      } elsif (@relevant_breaks) { 
        # relevant breaks are always a precise positiom
        my $mid = int(scalar(@relevant_breaks) / 2);
        
        $split_reg_start = $relevant_breaks[$mid]->end;
        $split_reg_end = $relevant_breaks[$mid]->start;
        
        printf($log_outfh "EXTEND: Split between adjacent contigs (%s/%d-%d -> %d-%d)\n",
               $sid, 
               $prev_el->end, 
               $this_el->start,
               $split_reg_start,
               $split_reg_end);
      } else {
        printf($log_outfh "EXTEND: Split within contig region (%s/%d-%d)\n",
               $sid, 
               $prev_el->end, 
               $this_el->start,
               );
        # make sure that we do not break in a non-break region
        if (@non_break_regions) {
          for(my $i=1; $i < @non_break_regions; $i++) {
            my $okay_reg_start = $non_break_regions[$i-1]->{end} + 1;
            my $okay_reg_end   = $non_break_regions[$i]->{start} - 1;

            if ($okay_reg_start >= $split_reg_start and
                $okay_reg_end  <= $split_reg_end) {
              $split_reg_start = $okay_reg_start;
              $split_reg_end   = $okay_reg_end;
              last;
            }
          }
        }
      }      
      
      my $new_prev_end = int((($split_reg_start + $split_reg_end)/2));
      my $new_this_start = $new_prev_end + 1;
            
      $prev_el->end($new_prev_end);
      $this_el->start($new_this_start);
    }     

    $sorted_els[-1]->end($tsl->length);
  } 

  return \@merged_chains;
}



#################################################################
# fetch_slices_and_component_maps
#
#################################################################
sub fetch_slices_and_component_maps {
  my ($db, $gene_scaffolds) = @_;

  my (%slices, %slice_maps);

  foreach my $gs_id (keys %$gene_scaffolds) {
    foreach my $comp (@{$gene_scaffolds->{$gs_id}->{components}}) {
      my $sid = $comp->to->id; 

      next if $sid =~ /^$FAKE_SCAFFOLD_PREFIX/;
      
      if (not exists($slices{$sid})) {
        $slices{$sid} = $db->get_SliceAdaptor->fetch_by_region('toplevel', 
                                                               $sid);
        
        # get components here
        my @seq_lev = @{$slices{$sid}->project('seqlevel')};

        my $map = Bio::EnsEMBL::Mapper->new('seqlevel',
                                            'scaffold');

        foreach my $comp (@seq_lev) {
          $map->add_map_coordinates($comp->to_Slice->seq_region_name,
                                    $comp->to_Slice->start,
                                    $comp->to_Slice->end,
                                    $comp->to_Slice->strand,
                                    $sid,
                                    $comp->from_start,
                                    $comp->from_end);

        }
        
        $slice_maps{$sid} = $map;
      }
    }
  }

  return (\%slices, \%slice_maps);
}



#################################################################
# is_simple_gene_scaffold
#
#################################################################
sub is_simple_gene_scaffold {
  my ($gene_scaf_id, $all_gene_scaffs) = @_;
  
  # a simple gene scaffold is one in which 
  # (a) all components come from the same scaffold
  # (b) all are in the same orientation
  # (c) the order is consistent

  my $gene_scaf = $all_gene_scaffs->{$gene_scaf_id};

  for(my $i=1; $i < @{$gene_scaf->{components}}; $i++) {
    my $this = $gene_scaf->{components}->[$i];
    my $prev = $gene_scaf->{components}->[$i-1];

    if ($this->to->id ne $prev->to->id or
        $this->to->strand != $prev->to->strand or
        ($this->to->strand > 0 and $this->to->start <= $prev->to->end) or
        ($this->to->strand < 0 and $this->to->end >= $prev->to->start)) {
      return 0;
    }
  }

  return 1; 
}

#################################################################
# get_unique_fake_scaffold_id
#################################################################
sub get_unique_fake_scaffold_id {
  
  my $fake_name = $FAKE_SCAFFOLD_PREFIX . $fake_scaffold_count;
  $fake_scaffold_count++;
  
  return $fake_name;
}


#####################################
#
# file index and read functions
#
#####################################

my (@gs_fhs, 
    @gene_fhs,
    %gs_file_index,
    %gene_file_index);

sub index_gene_scaffold_files {
  my @files = @_;

  my (%gs_index, %s_index, $last_gs_id);

  for(my $i=0; $i < @files; $i++) {
    my $file = $files[$i];

    open my $fh, $file or die "Could not open '$file' for reading\n";
    $gs_fhs[$i] = $fh;

    for(my $curpos = tell $fh; $_ = <$fh>; $curpos = tell $fh) {
      my ($gs_id, $s_id);
      
      if (/^\#.+AGP.+scaffold\s+(\S+)\s+/) {
        $gs_id = $1;
      } elsif ($_ !~ /^\#/) {
        my @l = split(/\t/, $_);
        
        if ($l[4] ne 'N') {
          $gs_id = $l[0];
          $s_id  = $l[5];        
        }
      }

      if (defined $gs_id) {
        if (not exists $gs_file_index{$gs_id}) {
          $gs_file_index{$gs_id} = {
            fh_index => $i,
            from     => $curpos,
            to       => $curpos,
          };
        } else {        
          if (defined $last_gs_id and $gs_id ne $last_gs_id) {
            die "Error in gene scaffold files(s) $gs_id appears in more than one block\n";
          }
          $gs_file_index{$gs_id}->{to} = $curpos;
        }

        if (defined $s_id) {
          # maintain 2 other indices for the single-linkage clustering

          $s_index{$s_id}->{$gs_id} = 1;
          $gs_index{$gs_id}->{$s_id} = 1;
        }

        $last_gs_id = $gs_id;
      }     
    }
  }

  return (\%gs_index, \%s_index);
}


sub index_gene_annotation_files {
  my @files = @_;

  my ($last_gs_id);

  for(my $i=0; $i < @files; $i++) {
    my $file = $files[$i];
    
    open my $fh, $file or die "Could not open '$file' for reading\n";
    $gene_fhs[$i] = $fh;

    for(my $curpos = tell $fh; $_ = <$fh>; $curpos = tell $fh) {
      my $gs_id;

      if (/^\#\s+Gene\s+report\s+for\s+(\S+)/) {
        $gs_id = $1;
      } elsif ($_ !~ /^\#/) {
        my @l = split(/\t/, $_);

        $gs_id = $l[0];
      }
       
      if (defined $gs_id) {
        if (not exists $gene_file_index{$gs_id}) {
          $gene_file_index{$gs_id} = {
            fh_index => $i,
            from     => $curpos,
            to       => $curpos,
          };
        } else {
          if (defined $last_gs_id and $gs_id ne $last_gs_id) {
            die "Error in gene annotation(s) files: $gs_id appears in more than one block\n";
          }
          $gene_file_index{$gs_id}->{to} = $curpos;
        }

        $last_gs_id = $gs_id;
      }
    }
  }
}


sub fetch_gene_scaffold_entries {
  my ($gene_scaffolds, @gs_ids) = @_;

  foreach my $gs_id (@gs_ids) {

    my $index = $gs_file_index{$gs_id};
    my $fh = $gs_fhs[$index->{fh_index}];

    seek $fh, $index->{from}, 0;

    my %scaffolds_seen;

    for(my $curpos = tell $fh; $_ = <$fh>; $curpos = tell $fh) {

      last if $curpos > $index->{to};

      /^\#.+AGP.+scaffold\s+(\S+)\s+.+region\=(\S+)\/(\d+)\-(\d+)/ and do {
        if ($1 ne $gs_id) {
          die "Index lookup error: Looking for $gs_id, found $1\n";
        }

        $gene_scaffolds->{$gs_id}->{id}        = $gs_id;
        $gene_scaffolds->{$gs_id}->{ref_name}  = $2;
        $gene_scaffolds->{$gs_id}->{ref_start} = $3;
        $gene_scaffolds->{$gs_id}->{ref_end}   = $4;

        next;
      };


      /^(\S+)\s+(\d+)\s+(\d+)\s+\S+\s+\S+\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)/ and do {
        my ($this_gs_id, 
            $gs_start, 
            $gs_end, 
            $scaf_id, 
            $scaf_start, 
            $scaf_end, 
            $scaf_ori) = ($1, $2, $3, $4, $5, $6, $7);

        $scaffolds_seen{$scaf_id} = 1;

        if ($this_gs_id ne $gs_id) {
          die "Index lookup error: Looking for $gs_id, found $this_gs_id\n";
        }

        $scaf_ori = ($scaf_ori eq '+') ? 1 : -1;
        
        my ($level) = ($gs_id =~ /\-(\d+)$/); 
        $gene_scaffolds->{$gs_id}->{level}      = $level;
        
        ###
        # finally, record the AGP line itself
        ###
        my $asm_unit = Bio::EnsEMBL::Mapper::Coordinate->new($gs_id, 
                                                             $gs_start, 
                                                             $gs_end,
                                                             1);
        my $cmp_unit = Bio::EnsEMBL::Mapper::Coordinate->new($scaf_id, 
                                                             $scaf_start, 
                                                             $scaf_end,
                                                             $scaf_ori);
        my $pair = Bio::EnsEMBL::Mapper::Pair->new($asm_unit, 
                                                   $cmp_unit); 

        push @{$gene_scaffolds->{$gs_id}->{components}}, $pair;

      }
    }    
  }
}



sub fetch_gene_annotation_entries {
  my ($gene_scaffolds, @gs_ids) = @_;

  foreach my $gs_id (@gs_ids) {
    my $index = $gene_file_index{$gs_id};
    my $fh = $gene_fhs[$index->{fh_index}];

    seek $fh, $index->{from}, 0;
    for(my $curpos = tell $fh; $_ = <$fh>; $curpos = tell $fh) {

      last if $curpos > $index->{to};
      
      if (/^\#\#\-ATTRIBUTE/) {
        push @{$gene_scaffolds->{$gs_id}->{annotation}}, $_;
      } elsif ($_ !~ /^\#/) {
        my @l = split(/\t/, $_);

        push @{$gene_scaffolds->{$gs_id}->{annotation}}, \@l;
      }
    }
  }
}
