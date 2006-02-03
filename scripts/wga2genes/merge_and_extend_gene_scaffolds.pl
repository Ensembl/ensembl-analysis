#!/usr/local/ensembl/bin/perl

use strict;
use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Mapper;


my ($dbname, 
    $dbuser,
    $dbhost,
    $dbport,
    $verbose,
    $agp_outfile, $agp_outfh,
    $gff_outfile, $gff_outfh,
    $log_outfile, $log_outfh,
    $reject_contig_splits,
    $filter_low_level,
    @agp_files,
    @gff_files);

&GetOptions('agp=s@' => \@agp_files,
            'gff=s@' => \@gff_files,
            'outagp=s' => \$agp_outfile,
            'outgff=s' => \$gff_outfile,
            'outlog=s' => \$log_outfile,
            'verbose' => \$verbose,
            'contig_split_reject' => \$reject_contig_splits,
            'low_level_reject'     => \$filter_low_level,
            'dbname=s' => \$dbname,
            'dbhost=s' => \$dbhost,
            'dbuser=s' => \$dbuser,
            'dbport=s' => \$dbport);

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
if (defined $gff_outfile) {
  open $gff_outfh, ">$gff_outfile" 
      or die "Could not open $gff_outfile for writing\n";
} else {
  $gff_outfh = \*STDOUT;
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
&index_gene_annotation_files(@gff_files);

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
  &remove_unused_components_from_gene_scaffolds($gene_scaffolds);

  # add in the gaps we have "plugged" as fake components
  &identify_plugged_gaps($gene_scaffolds);

  # 
  foreach my $gs_id (keys %$gene_scaffolds) {
    if (not &prune_gap_ends_from_gene_scaffold($gene_scaffolds->{$gs_id})) {
      delete $gene_scaffolds->{$gs_id};
    }
  }

  # gene scaffolds that share a scaffold but with different 
  # orientation are reoriented (if poss). 
  &reorient_gene_scaffolds($gene_scaffolds) if scalar(@$cluster) > 1;

  # if low-level gene scaffolds conflict with a higher one, it
  # is rejected
  &filter_low_level_gene_scaffolds($gene_scaffolds) 
      if $filter_low_level and scalar(@$cluster) > 1;

  # we need the slices and their seq-level component structure for 
  # assessing scaffold splits
  my ($target_slices, $target_slice_maps) = 
      &fetch_slices_and_component_maps($DB, $gene_scaffolds);

  # make list of "non trivial" gene scaffolds, which are those
  # that involve more than one component

  my (%simple_gene_scaffolds, %complex_gene_scaffolds);
  foreach my $gs_id (keys %$gene_scaffolds) {
    if (&is_simple_gene_scaffold($gene_scaffolds->{$gs_id})) {
      $simple_gene_scaffolds{$gs_id} = $gene_scaffolds->{$gs_id};
    } else {
      $complex_gene_scaffolds{$gs_id} = $gene_scaffolds->{$gs_id};
    }
  }

  my $extended_chains;
  while (not defined $extended_chains) {
    # cluster the members into sub-clusters based on a single,
    # shared scaffold
    my $chains = &make_gene_scaffold_chains(\%complex_gene_scaffolds);

    # extend components so that the whole of each scaffold is used 
    # somewhere
    $extended_chains = 
        &extend_gene_scaffold_components(\%complex_gene_scaffolds,
                                         $chains,
                                         $target_slices,
                                         $target_slice_maps,
                                         $reject_contig_splits);
  }

  # make a map to/from Gene scaffold coords and component coords
  my $orig_map = Bio::EnsEMBL::Mapper->new('scaffold',
                                           'genescaffold');
  
  foreach my $gs_id (keys %$gene_scaffolds) {
    my $gs = $gene_scaffolds->{$gs_id};
    foreach my $seg (@{$gs->{components}}) {
      $orig_map->add_map_coordinates($seg->to->id,
                                     $seg->to->start,
                                     $seg->to->end,
                                     $gs->{flipped} ? $seg->ori * -1 : $seg->ori,
                                     $gs_id,
                                     $seg->from->start,
                                     $seg->from->end,
                                     );
    }
  }

  foreach my $chain (@$extended_chains) {
    my $new_map     = Bio::EnsEMBL::Mapper->new('scaffold',
                                                'genescaffold');
    
    my $new_gs_name = $GENE_SCAFFOLD_PREFIX . $out_gene_scaffold_count++;
    print $agp_outfh "##-AGP for $new_gs_name [@{$chain->{members}}]\n";

    my @unit_objs   = @{$chain->{unit_objs}};
    
    my $last_end = 0;
    for(my $i=0; $i < @unit_objs; $i++) {
      my $unit_obj = $unit_objs[$i];
      my $unit_len = $unit_obj->{unit}->end - $unit_obj->{unit}->start + 1;
      
      my $gs_start = $last_end + 1;
      $last_end    = $gs_start + $unit_len - 1;
      
      $new_map->add_map_coordinates($unit_obj->{unit}->id,
                                    $unit_obj->{unit}->start,
                                    $unit_obj->{unit}->end,
                                    $unit_obj->{ori},
                                    $new_gs_name,
                                    $gs_start,
                                    $last_end,
                                    );
      
      if ($unit_obj->{unit}->id =~ /^$FAKE_SCAFFOLD_PREFIX/) {
        printf($agp_outfh "%s\t%d\t%d\tN\t%d\n", 
               $new_gs_name, 
               $gs_start,
               $last_end,
               $last_end - $gs_start + 1);       
      } else {
        printf($agp_outfh "%s\t%d\t%d\tW\t%s\t%d\t%d\t%s\n", 
               $new_gs_name, 
               $gs_start, 
               $last_end,
               $unit_obj->{unit}->id,
               $unit_obj->{unit}->start,
               $unit_obj->{unit}->end,
               $unit_obj->{ori} > 0 ? "+" : "-");        
      }
      
      if ($i < @unit_objs - 1) {
        $gs_start = $last_end + 1;
        $last_end = $gs_start + $GENE_SCAFFOLD_PADDING - 1;
        
        printf($agp_outfh "%s\t%d\t%d\tN\t%d\n", 
               $new_gs_name, 
               $gs_start,
               $last_end,
               $last_end - $gs_start + 1);
        
      }
    }
    
    printf $gff_outfh "# GFF for $new_gs_name [@{$chain->{members}}]\n";
    
    foreach my $member (@{$chain->{members}}) {
      
      foreach my $line (@{$gene_scaffolds->{$member}->{annotation}}) {
        if (ref($line) ne "ARRAY") {
          print $gff_outfh $line;
          next;
        } 
        
        my ($loc_in_orig) = $orig_map->map_coordinates($line->[0],
                                                       $line->[3],
                                                       $line->[4],
                                                       1,
                                                       'genescaffold');
        
        my ($loc_in_new) = $new_map->map_coordinates($loc_in_orig->id,
                                                     $loc_in_orig->start,
                                                     $loc_in_orig->end,
                                                     1,
                                                     'scaffold');
        
        $line->[0] = $new_gs_name;
        $line->[3] = $loc_in_new->start;
        $line->[4] = $loc_in_new->end;
        print $gff_outfh join("\t", @$line);
      }
    }
    
    # we need to deal with the simple gene scaffolds which refer
    # to a component that has ended up in this gene scaffold
    foreach my $gs_id (keys %simple_gene_scaffolds) {
      my @cmps = @{$simple_gene_scaffolds{$gs_id}->{components}};
      my $num_found = 0;
      
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
          $num_found++;
        }          
      }
      
      if ($num_found == scalar(@cmps)) {
        foreach my $line (@{$simple_gene_scaffolds{$gs_id}->{annotation}}) {
          if (ref($line) ne "ARRAY") {
            print $gff_outfh $line;
            next;
          } 
          
          my ($loc_in_orig) = $orig_map->map_coordinates($line->[0],
                                                         $line->[3],
                                                         $line->[4],
                                                         1,
                                                         'genescaffold');
          
          my ($loc_in_new) = $new_map->map_coordinates($loc_in_orig->id,
                                                       $loc_in_orig->start,
                                                       $loc_in_orig->end,
                                                       1,
                                                       'scaffold');
          
          $line->[0] = $new_gs_name;
          $line->[3] = $loc_in_new->start;
          $line->[4] = $loc_in_new->end;
          print $gff_outfh join("\t", @$line);
        }
        
        delete $simple_gene_scaffolds{$gs_id};

      } elsif ($num_found) {
        print $log_outfh "WRITE: Simple gs $gs_id is split between >1 complex ones\n";
      }
    }
  }

  # for all of the remaining simple gene scaffolds, the annotation 
  # can be written directly without need for agp

  foreach my $gs_id (keys %simple_gene_scaffolds) {
    printf $gff_outfh "# GFF for %s\n", $gs_id; 
      
    foreach my $line (@{$simple_gene_scaffolds{$gs_id}->{annotation}}) {
      if (ref($line) ne "ARRAY") {
        print $gff_outfh $line;
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
      
      if ($simple_gene_scaffolds{$gs_id}->{components}->[0]->ori < 0) {
        $line->[5] *= -1;
      }
      print $gff_outfh join("\t", @$line);
    }
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

      my $fake_q = Bio::EnsEMBL::Mapper::Unit->new($gs_id, 
                                                   $reg->{start}, 
                                                   $reg->{end});
      my $fake_t = Bio::EnsEMBL::Mapper::Unit->new($fake_name, 
                                                   $fake_start, 
                                                   $fake_end);
      my $pair   = Bio::EnsEMBL::Mapper::Pair->new($fake_q, 
                                                   $fake_t, 
                                                   1);
      
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
# remove_unused_components_from_gene_scaffolds
#
#   removes gene scaffold components that have no exons
#################################################################
sub remove_unused_components_from_gene_scaffolds {
  my ($gene_scaffolds) = @_;

  foreach my $gs_id (keys %$gene_scaffolds) {
    my @lines;
    foreach my $line (@{$gene_scaffolds->{$gs_id}->{annotation}}) {
      if (ref($line) eq "ARRAY" and
          $line->[2] eq 'Exon') {
        push @lines, $line;
      }
    }
    @lines = sort {$a->[3] <=> $b->[3]} @lines;

    my @comps = sort {
      $a->from->start <=> $b->from->start;
    } @{$gene_scaffolds->{$gs_id}->{components}};

    my %dropped_comps;

    my $j=0;
    for(my$i=0; $i < @comps; $i++) {
      my $comp = $comps[$i];

      my $retain = 0;
      for(; $j < @lines; $j++) {
        if ($lines[$j]->[3] > $comp->from->end) {
          last;
        } elsif ($lines[$j]->[4] >= $comp->from->start and 
                 $lines[$j]->[3] <= $comp->from->end) {
          # overlap; this comp has something on it
          $retain = 1;
          last;
        } 
      }

      if (not $retain) {
        printf($log_outfh "UNUSED: removed %s/%d-%d from $gs_id because had no exons\n",
               $comp->to->id, 
               $comp->to->start, 
               $comp->to->end);
        $dropped_comps{$comp} = 1;
      } 
    }

    my @retained_comps;
    foreach my $comp (@{$gene_scaffolds->{$gs_id}->{components}}) {
      if (not exists $dropped_comps{$comp}) {
        push @retained_comps, $comp;
      }
    }
 
    $gene_scaffolds->{$gs_id}->{components} = \@retained_comps;
  }
}

#################################################################
# reorient_gene_scaffolds
#
#   Attempt is made to reorient gene scaffolds that share a common
#   scaffold component but in a different orientation
#################################################################
sub reorient_gene_scaffolds {
  my $gene_scaffolds = shift;

  my (%s_index, %gs_index);
  foreach my $gs_id (keys %{$gene_scaffolds}) {    
    my $extents = 
        &get_gene_scaffold_component_extents($gene_scaffolds->{$gs_id});
    foreach my $s_id (keys %$extents) {
      $s_index{$s_id}->{$gs_id} = 1;
      $gs_index{$gs_id}->{$s_id} = 1;
    }
  }

  my @s_ids = grep { scalar(keys %{$s_index{$_}}) > 1 } keys %s_index;

  # for each scaffold x:
  #   does x occur in different orientations in different scaffolds? 
  #   are the different orientations completely consistent within 
  #   each single gene scaffold? 
  #   If yes to both, choose a scaffold to flip. If poss, the other 
  #   components in the gene-scaffold should not be involved elsewhere.
  #   To break ties, choose the gs with the fewest components to flip
  
  SHARED_SCAF: foreach my $s_id (@s_ids) {
    my (%extents_by_s_id);
    
    foreach my $gs_id (keys %{$s_index{$s_id}}) {
      my $extents = 
          &get_gene_scaffold_component_extents($gene_scaffolds->{$gs_id});
      foreach my $o_s_id (keys %$extents) {
        $extents_by_s_id{$o_s_id}->{$gs_id} = $extents->{$o_s_id}; 
      }
    }
   
    my @gs_ids = keys %{$extents_by_s_id{$s_id}};

    my (@forward, @reverse);
    foreach my $gs_id (@gs_ids) {
      if (exists($extents_by_s_id{$s_id}->{$gs_id}->{ori}->{1})) {
        push @forward, $gs_id;
      }
      if (exists($extents_by_s_id{$s_id}->{$gs_id}->{ori}->{-1})) {
        push @reverse, $gs_id;
      }
    }

    if (@forward and @reverse) {
      # this component appears in different orientations

      my %bad;
      foreach my $f (@forward) {
        foreach my $r (@reverse) {
          if ($f eq $r) {
            # forward and reverse component found in a single gene 
            $bad{$f} = 1;
            last;
          }
        }
      }

      @forward = grep { not exists $bad{$_} } @forward;
      @reverse = grep { not exists $bad{$_} } @reverse;
      
      foreach my $list (\@forward, \@reverse) {
        my @gs_ids = @$list;
        foreach my $gs_id (@gs_ids) {
          # need to check that none of the other components in this gene
          # scaffold exist in other gene_scaffold. If not, we can flip it
          # with impunity
          my $flip_would_have_consequences = 0;
          foreach my $o_s_id (keys %{$gs_index{$gs_id}}) {
            # is this scaffold used elsewhere?             
            next if $o_s_id eq $s_id;
            if (scalar(keys %{$s_index{$o_s_id}}) > 1) {
              $flip_would_have_consequences = 1;
              last;
            }
          }
          if ($flip_would_have_consequences) {
            # cannot therefore flip ANY in forward
            @$list = ();
            last;
          }
        }
      }

      # anything remaining in forward and reverse can be flipped
      # flip the side with the fewest total components
      if (@forward and not @reverse) {
        map { &flip_gene_scaffold($gene_scaffolds->{$_}) } @forward;
      } elsif (@reverse and not @forward) {
        map { &flip_gene_scaffold($gene_scaffolds->{$_}) } @reverse;
      } elsif (@forward and @reverse) {
        my (%total_components); 
        foreach my $obj ({ key => 'forward', list => \@forward }, 
                         { key => 'reverse', list => \@reverse }) {
          
          foreach my $el (@{$obj->{list}}) {
            $total_components{$obj->{key}} += 
                scalar(@{$gene_scaffolds->{$el}->{components}});
          }
        }
        
        if ($total_components{forward} > $total_components{reverse}) {
          map { &flip_gene_scaffold($gene_scaffolds->{$_}) } @reverse;
        } else {
          map { &flip_gene_scaffold($gene_scaffolds->{$_}) } @forward;
        }
      }
    }
  }
}



#################################################################
# filter_low_level_gene_scaffolds
#
#   this function filters out GSs that are inconsistent with GSs 
#   at a higher ("more confident") level. AGPs that are internally 
#   inconsistent (either within themselves, or with respect to 
#   other AGPs at the same level) are retained by this filter
#################################################################
sub filter_low_level_gene_scaffolds {
  my ($gene_scaffolds) = @_; 

  my %gene_scaffolds_by_level;
  foreach my $gs_id (keys %$gene_scaffolds) {
    my $gs = $gene_scaffolds->{$gs_id};
    $gene_scaffolds_by_level{$gs->{level}}->{$gs_id} = $gs;
  }

  my @rejected_gs_ids;

  my @levels = sort { $a <=> $b } keys %gene_scaffolds_by_level;
  for(my $i=0; $i<@levels; $i++) {
    my $this_level = $levels[$i];
    my @gs_at_this_level = values %{$gene_scaffolds_by_level{$this_level}};

    foreach my $gs (@gs_at_this_level) {
      # need to check for consistency against all lower-level 
      # gene_scaffolds

      my $gs_id = $gs->{id};
      my $extents = &get_gene_scaffold_component_extents($gs);

      my $inconsistent = 0;
      PREV_LEVELS: for(my $j=0; $j<$i; $j++) {
        my @gs_at_lower_level = 
            values %{$gene_scaffolds_by_level{$levels[$j]}};

        BETTER_AGPS: foreach my $better_gs (@gs_at_lower_level) {
          # inconsistency is defined as overlap in the extents
          # in the extents of any shared scaffolds
          my $better_gs_id = $better_gs->{id};
          my $better_extents = 
              &get_gene_scaffold_component_extents($better_gs);

          foreach my $s_id (keys %{$extents}) {
            if (exists $better_extents->{$s_id}) {
              my $extent = $extents->{$s_id};
              my $better_extent = $better_extents->{$s_id};

              if ($better_extent->{start} <= $extent->{end} and
                  $better_extent->{end}   >= $extent->{start}) {
                # overlap
                $inconsistent = 1;
                printf($log_outfh "FILTER: Rejecting %s due to overlap with %s in $s_id\n", 
                      $gs_id, 
                      $better_gs_id);
                last PREV_LEVELS;
              } elsif ((exists($better_extent->{ori}->{-1}) and
                        not exists($better_extent->{ori}->{1}) and
                        exists($extent->{ori}->{1})) or 
                       (exists($better_extent->{ori}->{1}) and
                        not exists($better_extent->{ori}->{-1}) and
                        exists($extent->{ori}->{-1}))) {
                $inconsistent = 1;
                printf($log_outfh "FILTER: Rejecting %s due to strand difference %s in $s_id\n", 
                       $gs_id, 
                       $better_gs_id);
                last PREV_LEVELS;
              }
            }
          }
        }
      }

      if ($inconsistent) {
        delete $gene_scaffolds_by_level{$this_level}->{$gs_id};
        delete $gene_scaffolds->{$gs_id};
      }
    }
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

  if (scalar(keys %{$gene_scaffolds}) == 1) {
    my ($gs_id) = keys %{$gene_scaffolds};
    return [[$gs_id]];
  }


  # form chains based on component composition
  my (@chains, @new_chains);

  my %sub_clusters;
  foreach my $gs_id (keys %{$gene_scaffolds}) {
    foreach my $comp (@{$gene_scaffolds->{$gs_id}->{components}}) {
      $sub_clusters{$comp->to->id}->{$gs_id} = 1;
    }
  }


  foreach my $s_id (keys %sub_clusters) {
    my $linked_gs_ids = [keys %{$sub_clusters{$s_id}}];
    
    if (@$linked_gs_ids > 1) {
      my $success = 
          &sort_gene_scaffolds_by_common_component($gene_scaffolds, 
                                                   $s_id, 
                                                   $linked_gs_ids);
      if ($success) {
        push @chains, $linked_gs_ids;
      }
    }    
  }

  @chains = sort { scalar(@$b) <=> scalar(@$a) } @chains;

  # Remove chains that:
  # - subsumed by another chain...
  foreach my $chain (@chains) {
    my $found = 0;
    OCHAINS: foreach my $o_chain (@new_chains) {

      OCHAIN: for(my $i=0; $i < @$o_chain; $i++) {
        if ($o_chain->[$i] eq $chain->[0]) {
          for(my $j=0; $j < @$chain; $j++) {
            if ($chain->[$j] ne $o_chain->[$i+$j]) {
              last OCHAIN;
            }
          }
          $found = 1;
          last OCHAINS;
        }
      }
    }
    if (not $found) {
      push @new_chains, $chain;
    }
  }
  @chains = @new_chains;

  # Indentify confliciting sub-chains. Sub-chains conflict when they
  # share common members, and the members are not situated head-to-tail

  @new_chains = ();
  my (@conflicting_pairs);
  for(my $i=0; $i < @chains; $i++) {

    my %chain_els;
    map { $chain_els{$_} = 1 } @{$chains[$i]};
    
    for(my $j=$i+1; $j < @chains; $j++) {
      
      my @in_common;
      foreach my $el (@{$chains[$j]}) {
        push @in_common, $el if exists($chain_els{$el});
      }
      if (@in_common > 1) {
        push @conflicting_pairs, [$i,$j];
      } elsif (@in_common == 1) {
        if ($chains[$i]->[0] eq $in_common[0]) {
          if ($chains[$j]->[-1] ne $in_common[0]) { 
            push @conflicting_pairs, [$i,$j];
          }
        } elsif ($chains[$i]->[-1] eq $in_common[0]) {
          if ($chains[$j]->[0] ne $in_common[0]) {
            push @conflicting_pairs, [$i,$j];
          }
        } else {
          push @conflicting_pairs, [$i,$j];
        }
      }
    }
  }

  #foreach my $p (@conflicting_pairs) {
  #  printf($log_outfh "CHAIN : conflicting_pair %s : %s\n", 
  #        join(" ", @{$chains[$p->[0]]}), 
  #        join(" ", @{$chains[$p->[1]]}));  
  #}

  # need to resolve all conflicts by removing chains
  my %indices_to_remove;

  while (@conflicting_pairs) {
    my %occurrences;
    foreach my $pair (@conflicting_pairs) {
      $occurrences{$pair->[0]}++;
      $occurrences{$pair->[1]}++;    
    }

    my (@candidates, $max_occur);
    foreach my $idx (keys %occurrences) {
      if (not defined $max_occur or
          $occurrences{$idx} > $max_occur) {
        $max_occur = $occurrences{$idx};
        @candidates = ($idx);
      } elsif ($occurrences{$idx} == $max_occur) {
        push @candidates, $idx;
      }      
    }

    my @to_remove;
    if (scalar(@candidates) == 1) {
      @to_remove = @candidates;
    } else {
      # if one of the chains is lower-ranked than the others, remove it,
      # else remove all of them
      my $lowest_rank;
      
      foreach my $idx (@candidates) {
        my $lowest_rank_of_chain;
        foreach my $mem (@{$chains[$idx]}) {
          my $rank = $gene_scaffolds->{$mem}->{level};
          # low level gene scaffolds are high rank!
          if (not defined $lowest_rank_of_chain or
              $rank > $lowest_rank_of_chain) {
            $lowest_rank_of_chain = $rank;
          }
        }
        
        if (not defined $lowest_rank or
            $lowest_rank_of_chain > $lowest_rank) {
          @to_remove = ($idx);
          $lowest_rank = $lowest_rank_of_chain;
        } elsif ($lowest_rank_of_chain == $lowest_rank) {
          push @to_remove, $idx;
        }
      }
    }

    map { $indices_to_remove{$_} = 1 } @to_remove;
    my @new_conflicting_pairs;
    foreach my $pair (@conflicting_pairs) {
      if (not exists($indices_to_remove{$pair->[0]}) and
          not exists($indices_to_remove{$pair->[1]})) {
        push @new_conflicting_pairs, $pair;
      }
    }
    @conflicting_pairs = @new_conflicting_pairs;
  }

  for(my $i=0; $i < @chains; $i++) {
    if (not exists($indices_to_remove{$i})) {
      push @new_chains, $chains[$i];
    }
  }
  @chains = @new_chains;

  # now merge chains that match head-to-tail
  @new_chains = ();
  for(my $i=0; $i < @chains; $i++) {
    my $sl = $chains[$i];
    my $merged_in = 0;
    for(my $j=$i+1; $j < @chains; $j++) {
      my $o_sl = $chains[$j];

      if ($sl->[0] eq $o_sl->[-1]) {
        for(my $k=1; $k < @$sl; $k++) {
          push @$o_sl, $sl->[$k];
        }
        $merged_in = 1;
        last;
      } elsif ($sl->[-1] eq $o_sl->[0]) {
        for(my $k=@$sl-2; $k >= 0; $k--) {
          unshift @$o_sl, $sl->[$k];
        }
        $merged_in = 1;
        last;
      }
    }

    if (not $merged_in) {
      push @new_chains, $sl;
    }
  }
  @chains = @new_chains;

  # add back in all of the singletons
  my %in_list;
  foreach my $chain (@chains) {
    map { $in_list{$_} = 1 } @$chain;
  }
  foreach my $gs_id (keys %{$gene_scaffolds}) {
    if (not $in_list{$gs_id}) {
      push @chains, [$gs_id];
    }
  }


  return \@chains;
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
  my ($gene_scaffolds, 
      $chains,
      $target_slices, 
      $target_slice_maps,
      $reject_contig_splits) = @_;

  my (%units_by_scaffold);
  
  #
  # index-by-scaffold of all the agp fragments in this group
  #
  foreach my $chain (@$chains) {
    foreach my $gs_id (@$chain) {
      foreach my $unit_pair (@{$gene_scaffolds->{$gs_id}->{components}}) {
        my $unit = $unit_pair->to;
        my $ori = $unit_pair->ori;
        push @{$units_by_scaffold{$unit->id}}, {
          unit => $unit,
          ori  => $ori,
        };
      }
    }
  }

  #
  # join components together where possible
  #
  my @merged_chains;

  foreach my $chain (@$chains) {
    my (@units, @merged_units);

    foreach my $gs_id (@$chain) {
      my @components = @{$gene_scaffolds->{$gs_id}->{components}};
      foreach my $comp (@components) {
        push @units, {
          unit => Bio::EnsEMBL::Mapper::Unit->new($comp->to->id,
                                                  $comp->to->start,
                                                  $comp->to->end),
          ori  => $comp->ori,
          gs_ids => [$gs_id],
        }; 
      }
    }

    foreach my $uo (@units) {

      my $merged = 0;

      if (@merged_units and
          $merged_units[-1]->{ori} == $uo->{ori} and
          $merged_units[-1]->{unit}->id eq $uo->{unit}->id) {

        if ($uo->{ori} == 1) {
          if ($uo->{unit}->start > $merged_units[-1]->{unit}->end) {
            my $can_merge = 1;
            foreach my $o_uo (@{$units_by_scaffold{$uo->{unit}->id}}) {
              if ($o_uo->{unit}->start > $merged_units[-1]->{unit}->end and
                  $o_uo->{unit}->end   < $uo->{unit}->start) {
                $can_merge = 0;
                last;
              }
            }
            if ($can_merge) {
              $merged_units[-1]->{unit}->end($uo->{unit}->end);
              if ($merged_units[-1]->{gs_ids}->[-1] ne 
                  $uo->{gs_ids}->[0]) {
                push @{$merged_units[-1]->{gs_ids}}, $uo->{gs_ids}->[0];
              }
              $merged = 1;
            }
          }
        } else {
          if ($uo->{unit}->end < $merged_units[-1]->{unit}->start) {
            my $can_merge = 1;
            foreach my $o_uo (@{$units_by_scaffold{$uo->{unit}->id}}) {
              if ($o_uo->{unit}->start > $uo->{unit}->end and
                  $o_uo->{unit}->end < $merged_units[-1]->{unit}->start) {
                $can_merge = 0;
                last;
              }
            }
            if ($can_merge) {
              $merged_units[-1]->{unit}->start($uo->{unit}->start);
              if ($merged_units[-1]->{gs_ids}->[0] ne $uo->{gs_ids}->[0]) {
                unshift @{$merged_units[-1]->{gs_ids}}, $uo->{gs_ids}->[0];
              }
              $merged = 1;
            }
          }
        }
      }

      if (not $merged) {
        push @merged_units, $uo;
      } 
    }
    push @merged_chains, {
      members     => $chain,
      unit_objs   => \@merged_units,
    };
  }

  #
  # finally, extend scaffold components so that for any scaffold involved
  # in a gene scaffold, every bp of it is accounted for somewhere
  #
  %units_by_scaffold = ();
  foreach my $msc (@merged_chains) {
    foreach my $unit_obj (@{$msc->{unit_objs}}) {
      if ($unit_obj->{unit}->id !~ /^$FAKE_SCAFFOLD_PREFIX/) {
        push @{$units_by_scaffold{$unit_obj->{unit}->id}}, $unit_obj;        
      }
    }
  }

  my (%contig_split_scaffolds);

  foreach my $sid (keys %units_by_scaffold) {       
    my @sorted_els = sort {
      $a->{unit}->start <=> $b->{unit}->start;
    } @{$units_by_scaffold{$sid}};

    my $tsl = $target_slices->{$sid};
    my $map = $target_slice_maps->{$sid};
    my @contig_breaks;
    foreach my $comp ($map->map_coordinates($sid,
                                            1,
                                            $tsl->length,
                                            1,
                                            'scaffold')) {
      if ($comp->isa("Bio::EnsEMBL::Mapper::Gap")) {
        push @contig_breaks, $comp;
      }
    }

    $sorted_els[0]->{unit}->start(1);
    for(my $i=1; $i<@sorted_els; $i++) {
      my $prev_el = $sorted_els[$i-1];
      my $this_el = $sorted_els[$i];

      my @between_breaks;
      foreach my $break (@contig_breaks) {
        if ($break->end   > $prev_el->{unit}->end and
            $break->start < $this_el->{unit}->start) {
          push @between_breaks, $break;
        }
      }

      my $split_reg_start = $prev_el->{unit}->end;
      my $split_reg_end   = $this_el->{unit}->start;

      if (@between_breaks) {
        my $mid = int(scalar(@between_breaks) / 2);
        $split_reg_start = $between_breaks[$mid]->start;
        $split_reg_end = $between_breaks[$mid]->end;

        printf($log_outfh "EXTEND: Split in gap (%s/%d-%d -> %d-%d)\n", 
               $sid, 
               $prev_el->{unit}->end, 
               $this_el->{unit}->start,
               $split_reg_start,
               $split_reg_end);
      } else {

        push @{$contig_split_scaffolds{$sid}}, {
          left_gs  => $prev_el->{gs_ids}->[-1],
          right_gs => $this_el->{gs_ids}->[0], 
        };

        printf($log_outfh "EXTEND: Split in contig region (%s/%d-%d) between %s and %s\n", 
               $sid, 
               $prev_el->{unit}->end, 
               $this_el->{unit}->start,
               $prev_el->{gs_ids}->[-1],
               $this_el->{gs_ids}->[0]
               );
      }      

      my $new_prev_end = int((($split_reg_start + $split_reg_end)/2));
      my $new_this_start = $new_prev_end + 1;

      $prev_el->{unit}->end($new_prev_end);
      $this_el->{unit}->start($new_this_start);
    }     
    $sorted_els[-1]->{unit}->end($tsl->length);
  }

  # if specified, junk gene scaffolds that involve split scaffolds

  if ($reject_contig_splits) {
    my ($to_reject, %records_per_sid);

    foreach my $sid (keys %contig_split_scaffolds) {

      my (%bad_gs_ids);

      foreach my $split (@{$contig_split_scaffolds{$sid}}) {
        $bad_gs_ids{$split->{left_gs}}++;
        $bad_gs_ids{$split->{right_gs}}++;
      }

      # sort the gene scaffolds by (a) number of split contigs they 
      # have, followed by (b) coverage of the split contig within
      # the gene scaffold
      my @records;
      foreach my $gs_id (keys %bad_gs_ids) {
        my $rank = $bad_gs_ids{$gs_id};
        my $coverage = 0;
        
        my $gs = $gene_scaffolds->{$gs_id};
        foreach my $comp (@{$gs->{components}}) {
          if ($comp->to->id eq $sid) {
            $coverage += $comp->to->end - $comp->to->start + 1;
          }
        }
        push @records, {
          s_id  => $sid,
          gs_id => $gs_id,
          rank => $rank,
          coverage => $coverage,
        };
      }
      @records = sort {
        $b->{rank} <=> $a->{rank} or $a->{coverage} <=> $b->{coverage}
      } @records;
      $records_per_sid{$sid} = \@records;
    }
    # choose sid with lowest coverage to reject
    foreach my $sid (keys %records_per_sid) {
      my $rej_rec = $records_per_sid{$sid}->[0];
      if (not defined $to_reject or 
          $rej_rec->{coverage} < $to_reject->{coverage}) {
        $to_reject = $rej_rec;
      }
    }
    # finally, massage the gene scaffold, removing the offending scaffold
    if (defined $to_reject) {
      my $success =  
          &remove_scaffold_from_gene_scaffold($gene_scaffolds->{$to_reject->{gs_id}},
                                              $to_reject->{s_id},
                                              0.50);
      if (not $success) {
        printf $log_outfh "EXTEND: Rejecting due to contig split: %s\n", $to_reject->{gs_id};
        delete $gene_scaffolds->{$to_reject->{gs_id}};
      } else {
        printf $log_outfh "EXTEND: Removed %s from %s due to contig split: %s\n", $to_reject->{s_id}, $to_reject->{gs_id};
      }
      # the following will force chains to be recalculated
      return undef;
    }
  }

  return \@merged_chains;
}


#################################################################
# remove_scaffold_from_gene_scaffold
#
#################################################################
sub remove_scaffold_from_gene_scaffold {
  my ($gs, $sid) = @_;

  # replace all $sid components with "fake" components
  my @comps = @{$gs->{components}};

  foreach my $comp (@{$gs->{components}}) {
    if ($comp->to->id eq $sid) {
      my $new_id = &get_unique_fake_scaffold_id;
      my $new_to = Bio::EnsEMBL::Mapper::Unit->new($new_id,
                                                   $comp->to->start,
                                                   $comp->to->end);
      $comp->to($new_to);
    } 
  }

  return &prune_gap_ends_from_gene_scaffold($gs);
}

  

#################################################################
# prune_gap_ends_from_gene_scaffold_and_transcripts
#
#################################################################
sub prune_gap_ends_from_gene_scaffold {
  my ($gs) = @_;

  my @comps = @{$gs->{components}};

  # remove terminal fake components
  while(@comps and $comps[0]->to->id =~ /^$FAKE_SCAFFOLD_PREFIX/) {
    shift @comps;
  }
  while(@comps and $comps[-1]->to->id =~ /^$FAKE_SCAFFOLD_PREFIX/) {
    pop @comps;
  }
  
  # the following will lead to the rejection of the gene scaffold
  # if no components remained
  return 0 if not @comps;
  $gs->{components} = \@comps;

  # build a mini-map
  my $mini_map = Bio::EnsEMBL::Mapper->new('scaffold',
                                           'genescaffold');
  foreach my $comp (@comps) {
    $mini_map->add_map_coordinates($comp->to->id,
                                   $comp->to->start,
                                   $comp->to->end,
                                   $comp->ori,
                                   "thegenescaffold",
                                   $comp->from->start,
                                   $comp->from->end,
                                   );

  }

  my (%trans, %exon_trans);
  # sort the annotation lines into transcripts and exons
  foreach my $l (@{$gs->{annotation}}) {
    if (ref($l) eq "ARRAY") {
      if ($l->[2] eq 'Exon') {
        my ($tid) = ($l->[8] =~ /transcript=(\S+);/);
        my ($eid) = ($l->[8] =~ /exon=(\S+);/);
        $exon_trans{$eid} = $tid;

        push @{$trans{$tid}->{exons}}, {
          eid => $eid,
          start => $l->[3],
          end  => $l->[4],
          length => $l->[4] - $l->[3] + 1,
          line => $l,
        };
      } elsif ($l->[2] eq 'Supporting') {
        my ($eid) = ($l->[6] =~ /exon=(\S+);/);
        my ($hstart) = ($l->[6] =~ /hstart=(\d+)/);
        my ($hend) = ($l->[6] =~ /hend=(\d+)/);
        my $tid = $exon_trans{$eid};

        foreach my $eobj (@{$trans{$tid}->{exons}}) {
          if ($eobj->{eid} eq $eid) {
            push @{$eobj->{support}}, {
              hstart => $hstart,
              hend   => $hend,
              hlength => $hend - $hstart + 1,
              line => $l,
            };
            last;
          }
        }
      }
    } else {
      my ($tid) = ($l =~ /transcript=(\S+)/);
      my ($code) = ($l =~ /code=(\S+)/);
      my ($val) = ($l =~ /value=(\S+)/);
      
      if ($code eq 'HitCoverage') {
        $trans{$tid}->{hit_coverage} = $val;
      }

      push @{$trans{$tid}->{attributes}}, $l;      
    }
  }

  my @kept_trans;
  TRAN:foreach my $tid (keys %trans) {
    my @eobjs = sort {
      $a->{start} <=> $b->{start}
    } @{$trans{$tid}->{exons}};

    my $total_hres   = 0;
    my $removed_hres = 0;

    foreach my $e (@eobjs) {
      foreach my $sup (@{$e->{support}}) {
        $total_hres += $sup->{hlength};
      }
    }

    my $projected_hit_bps = 
        $total_hres / ($trans{$tid}->{hit_coverage} / 100);

    # remove terminal exons
    while(@eobjs) {
      my ($c) = $mini_map->map_coordinates('thegenescaffold',
                                           $eobjs[0]->{start},
                                           $eobjs[0]->{end},
                                           1,
                                           'genescaffold');
      if ($c->isa("Bio::EnsEMBL::Mapper::Gap")) {
        my $removed = shift @eobjs;
        foreach my $sup (@{$removed->{support}}) {
          $removed_hres += $sup->{hlength};
        }
      } else {
        last;
      }
    }
    while(@eobjs) {
      my ($c) = $mini_map->map_coordinates('thegenescaffold',
                                           $eobjs[-1]->{start},
                                           $eobjs[-1]->{end},
                                           1,
                                           'genescaffold');
      if ($c->isa("Bio::EnsEMBL::Mapper::Gap")) {
        my $removed = pop @eobjs;
        foreach my $sup (@{$removed->{support}}) {
          $removed_hres += $sup->{hlength};
        }
      } else {
        last;
      }
    }
    
    # if the target coverage has dropped below 50%, reject
    my $target_cov = ($total_hres - $removed_hres) / $projected_hit_bps;
    if ($target_cov < 0.50) {
      printf $log_outfh "REMOVAL: less than 50 hit coverage; chucking $tid\n";
      next TRAN;
    } else {
      my $cov_string = sprintf("%.2f", $target_cov * 100);
      my @new_attrs;
      foreach my $attr_line (@{$trans{$tid}->{attributes}}) {
        if ($attr_line =~ /code=HitCoverage/) {
          $attr_line =~ s/value=\S+/value=$cov_string/;
        }
        push @new_attrs, $attr_line;
      }
      $trans{$tid}->{attributes} = \@new_attrs;
    }

    # finally, check the non-gap coverage
    my $total_cov = 0;
    my $gap_cov   = 0;

    foreach my $eobj (@eobjs) {
      my ($c) = $mini_map->map_coordinates('thegenescaffold',
                                           $eobjs[0]->{start},
                                           $eobjs[0]->{end},
                                           1,
                                           'genescaffold');
      if ($c->isa("Bio::EnsEMBL::Mapper::Gap")) {
        # none of the remaining exons should map to gaps!        
        printf $log_outfh "REMOVAL: retained exon maps to real gap; chucking $tid\n";
        next TRAN;
      } else {
        if ($c->id =~ /^$FAKE_SCAFFOLD_PREFIX/) {
          $gap_cov += $eobj->{end} - $eobj->{start} + 1;
        }
        $total_cov += $eobj->{end} - $eobj->{start} + 1;
      }
    }

    my $prop_non_gap = ($total_cov - $gap_cov) / $total_cov;

    if ($prop_non_gap < 0.5) {
      next TRAN;
    } else {
      my $cov_string = sprintf("%.2f", $prop_non_gap * 100);
      my @new_attrs;
      foreach my $attr_line (@{$trans{$tid}->{attributes}}) {
        if ($attr_line =~ /code=PropNonGap/) {
          $attr_line =~ s/value=\S+/value=$cov_string/;
        }
        push @new_attrs, $attr_line;
      }
      $trans{$tid}->{attributes} = \@new_attrs;
    }

    $trans{$tid}->{exons} = \@eobjs;

    push @kept_trans, $tid;
  }

  if (@kept_trans) {
    my @new_annot_lines;
    foreach my $tid (@kept_trans) {
      push @new_annot_lines, @{$trans{$tid}->{attributes}};
      foreach my $eobj (@{$trans{$tid}->{exons}}) {
        push @new_annot_lines, $eobj->{line};
        foreach my $sup (@{$eobj->{support}}) {
          push @new_annot_lines, $sup->{line};
        }
      }
    }
    $gs->{annotation} = \@new_annot_lines;
    return 1;
  } else {
    printf $log_outfh "REMOVAL: no transcripts remaining for %s\n", $gs->{id};
    return 0;
  }
}


#################################################################
# flip_gene_scaffold
#
#################################################################
sub flip_gene_scaffold {
  my ($gene_scaffold) = @_;

  $gene_scaffold->{components} = 
      [reverse @{$gene_scaffold->{components}}];
  foreach my $comp (@{$gene_scaffold->{components}}) {
    $comp->ori( $comp->ori() * -1 );
  }

  foreach my $line (@{$gene_scaffold->{annotation}}) {
    if (ref($line) eq "ARRAY") {
      $line->[5] *= -1;
    }
  }
  $gene_scaffold->{flipped} = 1;
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
# get_gene_scaffold_component_extents
#
#################################################################
sub get_gene_scaffold_component_extents {
  my (@gs_list) = @_;

  my %extents;

  foreach my $gs (@gs_list) {
    foreach my $comp (@{$gs->{components}}) {
      if (not exists $extents{$comp->to->id}) {
        $extents{$comp->to->id} = {
          start => $comp->to->start,
          end   => $comp->to->end,
        };
      } else {
        if ($comp->to->start < $extents{$comp->to->id}->{start}) {
          $extents{$comp->to->id}->{start} = $comp->to->start;
        }
        if ($comp->to->end > $extents{$comp->to->id}->{end}) {
          $extents{$comp->to->id}->{end} = $comp->to->end;
        }
      }
      $extents{$comp->to->id}->{ori}->{$comp->ori} = 1;
    }
  }
  
  return \%extents;
}


#################################################################
# is_simple_gene_scaffold
#
#################################################################
sub is_simple_gene_scaffold {
  my $gene_scaf = shift;
  
  # a simple gene scaffold is one in which 
  # (a) all components come from the same scaffold
  # (b) all are in the same orientation
  # (c) the order is consistent
   
  for(my $i=1; $i < @{$gene_scaf->{components}}; $i++) {
    my $this = $gene_scaf->{components}->[$i];
    my $prev = $gene_scaf->{components}->[$i-1];

    if ($this->to->id ne $prev->to->id or
        $this->ori != $prev->ori or
        ($this->ori > 0 and $this->to->start <= $prev->to->end) or
        ($this->ori < 0 and $this->to->end >= $prev->to->start)) {
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

  my (%gs_index, %s_index);

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
        
        if ($l[3] eq 'W') {
          $gs_id = $l[0];
          $s_id  = $l[4];        
        }
      }

      if (defined $gs_id) {
        if (not exists $gs_file_index{$gs_id}) {
          $gs_file_index{$gs_id} = {
            fh_index => $i,
            from     => $curpos,
          };
        } 
        $gs_file_index{$gs_id}->{to} = $curpos;

        if (defined $s_id) {
          # maintain 2 other indices for the single-linkage clustering

          $s_index{$s_id}->{$gs_id} = 1;
          $gs_index{$gs_id}->{$s_id} = 1;
        }
      }     
    }
  }

  return (\%gs_index, \%s_index);
}


sub index_gene_annotation_files {
  my @files = @_;

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
          };
        } 
        $gene_file_index{$gs_id}->{to} = $curpos;
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


      /^(\S+)\s+(\d+)\s+(\d+)\s+\S+\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)/ and do {
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
        my $asm_unit = Bio::EnsEMBL::Mapper::Unit->new($gs_id, 
                                                       $gs_start, 
                                                       $gs_end);
        my $cmp_unit = Bio::EnsEMBL::Mapper::Unit->new($scaf_id, 
                                                       $scaf_start, 
                                                       $scaf_end);
        my $pair = Bio::EnsEMBL::Mapper::Pair->new($asm_unit, 
                                                   $cmp_unit, 
                                                   $scaf_ori);
        push @{$gene_scaffolds->{$gs_id}->{components}}, $pair;

      }
    }

    $gene_scaffolds->{$gs_id}->{unique_elements} = 
        scalar(keys %scaffolds_seen);
    $gene_scaffolds->{$gs_id}->{flipped} = 0;
    
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
