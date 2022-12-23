use warnings;
use strict;
use feature 'say';

use POSIX;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Getopt::Long qw(:config no_ignore_case);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name);

#use Statistics::Descriptive;
#use PointEstimation;

my $coord_system = 'toplevel';
#my $dbname = 'caenorhabditis_elegans_core_101_269';
my $dbname = 'homo_sapiens_core_104_38';


my $user = 'ensro';
my $host   = 'mysql-ens-genebuild-prod-1';
my $port   = 4527;
my $pass;
#my $genome_index = "/hps/nobackup/flicek/ensembl/genebuild/fergal/hprc/chm13y1/GCA_000001405.28/homo_sapiens_reheadered_toplevel.fa.mmi";
my $genome_index = "/hps/nobackup/flicek/ensembl/genebuild/fergal/hprc/chm13y1/GCA_009914755.4/homo_sapiens_reheadered_toplevel.fa.mmi";
my $options = GetOptions ("user|dbuser|u=s"      => \$user,
                          "host|dbhost|h=s"      => \$host,
                          "port|dbport|P=i"      => \$port,
                          "dbname|db|D=s"    => \$dbname,
                          "dbpass|pass|p=s" => \$pass,
                          "genome|g=s" => \$genome_index);



my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
  -port    => $port,
  -user    => $user,
  -host    => $host,
  -dbname  => $dbname,
  -pass    => $pass);

my $slice_adaptor = $db->get_SliceAdaptor();
my $sequence_adaptor = $db->get_SequenceAdaptor();
my $slices = $slice_adaptor->fetch_all('toplevel');
my $gene_adaptor = $db->get_GeneAdaptor();
my $meta_adaptor = $db->get_MetaContainer();
my $production_name = $meta_adaptor->get_production_name;


# Change this to cluster genes based on the 100kb window principle, then get region bioundaries, then make anchors
foreach my $slice (@$slices) {
  # Just for testing, this only looks at chromosome 4
  next unless($slice->seq_region_name eq 'X');

  my $original_anchor_positions = {};

  my $region_seq = $slice->seq();
  my $genes = $slice->get_all_Genes();
  say "Processing: ".$slice->seq_region_name." (".scalar(@$genes)." genes)";

  my $gene_clusters = cluster_slice_genes($genes,$slice);
  print_clusters($gene_clusters);
  my $meta_clusters = {};
  set_cluster_regions($gene_clusters,$meta_clusters,$slice);
  set_region_anchors($meta_clusters,$original_anchor_positions,$region_seq,$slice);
  my $fasta_record = generate_fasta_record($meta_clusters);

  my $fasta_record_path = write_input_file($fasta_record);
  map_anchors($meta_clusters,$fasta_record_path,$genome_index,$original_anchor_positions);

  path_finder($meta_clusters);

  check_control_accuracy($meta_clusters);

  path_selector($meta_clusters,$original_anchor_positions);

  refine_path_boundaries($meta_clusters,$original_anchor_positions);

#  assess_conflict($meta_clusters);

#  assess_regions($meta_clusters);
}


sub check_control_accuracy {
  my ($meta_clusters) = @_;

  # This subroutine is specifically for checking against a control, so GRCh38 against itself or if we implement some method of
  # rearranging the GRCh38 assembly and tracking the re-arrangements, it could also be adapted for that

  foreach my $cluster_id (keys(%$meta_clusters)) {
    my $paths = $meta_clusters->{$cluster_id}->{'paths'};
    my $top_path = ${$paths}[0];
    my ($region_start,$region_end) = @{$meta_clusters->{$cluster_id}->{'region_boundaries'}};
    my $start_hit = ${$top_path}[0];
    my $end_hit = ${$top_path}[-1];

#    use Data:Dumper;
#    say "Dumper 1: ".Dumper();
#    if(${$end) {

#    }
    say "For cluster ".$cluster_id." the coverage of the original region with the top path is ".overlap_length($region_start,$region_end,${$start_hit}[1],${$end_hit}[2]);
  }
}


sub refine_path_boundaries {
  my ($meta_clusters,$original_anchor_positions) = @_;

  foreach my $cluster_id (keys(%$meta_clusters)) {
    say "Calculating region boundaries for selected paths in cluster ".$cluster_id, join(':', '', @{$meta_clusters->{$cluster_id}->{region_boundaries}}, '');

    my $gene_cluster = $meta_clusters->{$cluster_id}->{'gene_cluster'};
    my $anchors = $meta_clusters->{$cluster_id}->{'anchors'};
    my $selected_paths = $meta_clusters->{$cluster_id}->{'selected_paths'};
    my $selected_path_genes = $meta_clusters->{$cluster_id}->{'selected_path_genes'};
    for my $i (0..scalar(@$selected_paths)) {
      my $path = ${$selected_paths}[$i];
      my $bounded_genes = ${$selected_path_genes}[$i];
      unless(defined($bounded_genes) and scalar(@$bounded_genes)) {
        next;
      }

      my $min_gene_start;
      my $max_gene_end;
      foreach my $gene (@$bounded_genes) {
        unless($min_gene_start) {
          $min_gene_start = $gene->seq_region_start();
          $max_gene_end = $gene->seq_region_end();
          next;
        }
        if($gene->seq_region_start() < $min_gene_start) {
          $min_gene_start = $gene->seq_region_start();
        }

        if($gene->seq_region_end() > $max_gene_end) {
          $max_gene_end = $gene->seq_region_end();
        }
      } # End foreach my $gene


      # At this point we know the boundaries of the genes, now we want to find the closest anchor ids to the boundaries
      # Specifically we want to find the closest start/end if the boundary overlaps the anchor or the closest end/start
      # if the anchor does not overlap the gene. Because of the way the paths are built and the genes are considered
      # bounded, there will always be a part of anchors that go past the genes on either end
      my $source_region_start;
      my $source_region_end;
      my $target_region_start;
      my $target_region_end;
      foreach my $hit (@$path) {
#        use Data::Dumper;
#        say "DUMPER: ".Dumper($hit);
        my $hit_start = ${$hit}[1];
        my $hit_end = ${$hit}[2];
        my $hit_strand = ${$hit}[3];
        my $hit_length = $hit_end - $hit_start + 1;
        my $anchor_id = ${$hit}[4];
        my $anchor = $original_anchor_positions->{$anchor_id};
        my $anchor_start = ${$anchor}[1];
        my $anchor_end = ${$anchor}[2];
        my $anchor_length = $anchor_end - $anchor_start + 1;

        # Beacuse the hits might not be full length we want to also factor this into calculating the target region boundaries
        my $length_adjust = $anchor_length - $hit_length;
        if($length_adjust < 0) {
          $length_adjust = 0;
        }

        # NOTE!!!!!!!!!!! Need to consider the effect of anchors being on the opposite strand. At the moment the path
        #                 is sorted based on the hit starts, so need to double check the logic still works. I suspect
        #                 it might because really everything is decided based on the reference, which is always on the
        #                 forward strand, but should double check
        # NOTE2!!!!!!!!!! From tests there are sometimes negative values for the target region coords, so I assume this
        #                 means it does matter. Because even though the anchor hits are sorted from 5' to 3', if the
        #                 region is on the negative strand, it will still mean that things are in the wrong orientation

        # This is the start overlap cast
        if($anchor_start < $min_gene_start and $anchor_end >= $min_gene_start) {
          $source_region_start = $anchor_start;
          if($hit_strand eq '+') {
            $target_region_start = $hit_start - $length_adjust;
          } else {
            $target_region_end = $hit_end + $length_adjust;
          }
        } elsif($anchor_end < $min_gene_start and (!$source_region_start or $anchor_end > $source_region_start)) {
          # This is sort of confusing, this anchor does not overlap the min start, so we want to select it if the
          # end is less than the min start (so it's 5' to the genes) and either we haven't already selected an
          # anchor yet, or the end of the anchor closer to the min start than the current region start
          $source_region_start = $anchor_start;
          if($hit_strand eq '+') {
            $target_region_start = $hit_start - $length_adjust;
          } else {
            $target_region_end = $hit_end + $length_adjust;
          }
#          $target_region_start = $hit_start - $length_adjust;
        }

        # Do essentially the same for the 3' end
        if($anchor_start <= $max_gene_end and $anchor_end > $max_gene_end) {
          $source_region_end = $anchor_end;
          if($hit_strand eq '+') {
            $target_region_end = $hit_end + $length_adjust;
          } else {
            $target_region_start = $hit_start - $length_adjust;
          }
        } elsif($anchor_start > $max_gene_end and (!$source_region_end or $anchor_start < $source_region_end)) {
          # This is sort of confusing, this anchor does not overlap the min start, so we want to select it if the
          # end is less than the min start (so it's 5' to the genes) and either we haven't already selected an
          # anchor yet, or the end of the anchor closer to the min start than the current region start
          $source_region_end = $anchor_end;
          if($hit_strand eq '+') {
            $target_region_end = $hit_end + $length_adjust;
          } else {
            $target_region_start = $hit_start - $length_adjust;
          }
#          $target_region_end = $hit_end + $length_adjust;
        }
      } # End foreach my $hit

      say "  Min/max gene coords: ".$min_gene_start."/".$max_gene_end.", length: ".($max_gene_end - $min_gene_start + 1).", bounded gene count: ".scalar(@$bounded_genes);
      say "  Source region coords: ".$source_region_start."/".$source_region_end.", length: ".($source_region_end - $source_region_start + 1);
      say "  Target region coords: ".$target_region_start."/".$target_region_end.", length: ".($target_region_end - $target_region_start + 1);
    } # End for my $i (0..scalar
  } # End foreach my $cluster_id

}

sub path_selector {
  my ($meta_clusters,$original_anchor_positions) = @_;

  # This goes through the paths for each cluster and decides what set of paths to keep in order to get the most coverage of the original genes

  # This subroutine is specifically for checking against a control, so GRCh38 against itself or if we implement some method of
  # rearranging the GRCh38 assembly and tracking the re-arrangements, it could also be adapted for that
  my $total_gene_count = 0;
  my $recovered_gene_count = 0;
  foreach my $cluster_id (keys(%$meta_clusters)) {
    say "Calculating bounded genes for paths in cluster ".$cluster_id, join(':', '', @{$meta_clusters->{$cluster_id}->{region_boundaries}}, '');

    my $gene_cluster = $meta_clusters->{$cluster_id}->{'gene_cluster'};
    my $anchors = $meta_clusters->{$cluster_id}->{'anchors'};
    my $paths = $meta_clusters->{$cluster_id}->{'paths'};
    my $found_genes_by_db_id = {};
    my $selected_paths = [];
    my $bounded_path_genes = [];
    foreach my $path (@$paths) {
      # We want to get the start and end anchors in the path and then understand what genes they cover in the original cluster
      # We then mark these genes as found. After that we go on to subsequent paths and look for paths covering missing genes
      # only, if there's a found gene in the path we just skip it. Because all paths are calculated, even if there was a sub
      # cluster that had a duplicated copy of a found gene on a particular path with missing genes there should still be other
      # sub paths not containing the found gene. There could be issues if the gene fell awkwardly in terms of anchors

      my $start_hit = ${$path}[0];
      my $end_hit = ${$path}[-1];
      my $start_anchor_id = ${$start_hit}[-1];
      my $end_anchor_id = ${$end_hit}[-1];
      # If the start anchor id is less than the end anchor id we must have a set on the opposite strand. since we are just
      # going to use these to calculate the genes falling within the boundaries in the original gene cluster we can swap
      if($start_anchor_id > $end_anchor_id) {
        my $tmp = $start_anchor_id;
        $start_anchor_id = $end_anchor_id;
        $end_anchor_id = $tmp;
      }

      my $bounded_genes = calculate_bounded_genes($gene_cluster,$original_anchor_positions,$start_anchor_id,$end_anchor_id);
      unless(scalar(@$bounded_genes)) {
        next;
      }

      my $conflict = 0;
      foreach my $gene (@$bounded_genes) {
        my $gene_db_id = $gene->dbID();
        if($found_genes_by_db_id->{$gene_db_id}) {
          $conflict = 1;
          last;
        }
      }

      if($conflict) {
        next;
      }

      # No genes in the bounded set have been list as found elsewhere at this point
      foreach my $gene (@$bounded_genes) {
        my $gene_db_id = $gene->dbID();
        $found_genes_by_db_id->{$gene_db_id} = 1;
      }
      push(@$selected_paths,$path);
      push(@$bounded_path_genes,$bounded_genes);
    } # foreach my $path

    my $recovered_gene_location_count = scalar(keys(%$found_genes_by_db_id));
    my $original_gene_count = scalar(@$gene_cluster);
    my $perc_recovered = sprintf("%.2f",($recovered_gene_location_count/$original_gene_count * 100));
    say "  Recovered ".$recovered_gene_location_count." expected gene locations based on the paths chosen";
    say "  Number of genes in original cluster is ".$original_gene_count;
    say "  Missing gene count is ".($original_gene_count - $recovered_gene_location_count);
    say "  Percent of genes covered with mapped anchor boundaries for cluster ".$cluster_id." is ".$perc_recovered."%";

    $total_gene_count += $original_gene_count;
    $recovered_gene_count += $recovered_gene_location_count;
    $meta_clusters->{$cluster_id}->{'selected_paths'} = $selected_paths;
    $meta_clusters->{$cluster_id}->{'selected_path_genes'} = $bounded_path_genes;
  } # foreach my $cluster_id

  my $recovery_percent = sprintf("%.2f",($recovered_gene_count/$total_gene_count * 100));
  say "Final genes covered by mapped anchor boundaries: ".$recovered_gene_count."/".$total_gene_count."(".$recovery_percent."%)";

}


sub calculate_bounded_genes {
  my ($gene_cluster,$original_anchor_positions,$start_anchor_id,$end_anchor_id) = @_;

  my $bounded_genes = [];
  my $start_anchor = $original_anchor_positions->{$start_anchor_id};
  my $end_anchor = $original_anchor_positions->{$end_anchor_id};

  my $boundary_start = ${$start_anchor}[1];
  my $boundary_end = ${$end_anchor}[2];

#  say "  Boundary to evaluate: ".$boundary_start.":".$boundary_end;
#  say "  Genes in original cluster: ".scalar(@$gene_cluster);
  foreach my $gene (@$gene_cluster) {
    if($gene->seq_region_start() >= $boundary_start and $gene->seq_region_end() <= $boundary_end) {
      push(@$bounded_genes,$gene);
    } else {
#      say "  gene falls outside of anchor boundaries: ".$gene->seq_region_start().":".$gene->seq_region_end();
    }
  }

#  say "  Found ".scalar(@$bounded_genes)." within boundary anchors for current path";
  return($bounded_genes);
}


sub overlap_length {
# return the length of the overlap between featureA and featureB
  my ($region_start,$region_end,$path_start,$path_end) = @_;

  my $min_end = $region_end;
  if ($path_end < $min_end) {
    $min_end = $path_end;
  }

  my $max_start = $region_start;
  if ($path_start > $max_start) {
    $max_start = $path_start;
  }

  my $overlap = $min_end-$max_start+1;
  my $overlap_percent = sprintf("%.3f",($overlap/($region_end-$region_start+1) * 100));
  return $overlap_percent;
}


sub path_finder {
  my ($meta_clusters) = @_;

  # The number of anchors that can be skipped to continue building the path
  my $max_skips = 3;

  # The max distance to search between the end of the current anchor and the start of the next member
  # of the path. This should always be a number greater than the anchor spacing in set_region_anchors
  my $max_dist = 5000;

  # The expected distance is based on the anchor spacing, it should be the same value
  my $expected_dist = 500;

  foreach my $cluster_id (keys(%$meta_clusters)) {
    my $anchors = $meta_clusters->{$cluster_id}->{'anchors'};

    say "Processing cluster to find best paths for cluster ".$cluster_id." with ".scalar(@$anchors)." anchors";

    # What we want to do here is take each anchor in order of its id and then calculte the path for each hit
    # under the constraints above

    # First need to understand what the min/max anchor ids are and then also convert the anchor into into a hash with the
    # anchor ids as the keys. Anchor ids are

    # If the array is ordered by id then it should be a case of going element by element and assessing the maximal path
    # achievable given the constraints

    my $paths = [];
    my @sorted_anchors = sort { $a->[0] <=> $b->[0] } @{$anchors};
    for(my $i=0; $i<scalar(@sorted_anchors)-1; $i++) {
      my $anchor_i = $sorted_anchors[$i];
      my $anchor_i_id = ${$anchor_i}[0];
      my $anchor_i_hits = ${$anchor_i}[-1];
#      say "Processing paths for anchor ".$anchor_i_id.", ".scalar(@$anchor_i_hits)." starting points";
      # Foreach hit we want to see how far the path can go before we find another hit within the constraints
      # The assumption is that since the anchors are sequential we can just move through them until we find a
      # hit that's the expected next step on the path or the constraints are failed and we stop
      foreach my $hit_i (@{$anchor_i_hits}) {
        my $path_hit_i = [$hit_i];
        my $hit_i_region = ${$hit_i}[0];
        my $hit_i_start = ${$hit_i}[1];
        my $hit_i_end = ${$hit_i}[2];
        my $hit_i_strand = ${$hit_i}[3];
#        say "  Starting hit: ".$hit_i_region.":".$hit_i_start.":".$hit_i_end.":".$hit_i_strand;

        my $path_end_hit = $hit_i;
        my $path_end_id = $anchor_i_id;
        for(my $j=$i+1; $j<scalar(@sorted_anchors); $j++) {
          # Here we want to continue from j, the next anchor as far as possible through the remainder
          # of the anchors
          my $anchor_j = $sorted_anchors[$j];
          my $anchor_j_id = ${$anchor_j}[0];
          my $anchor_j_hits = ${$anchor_j}[-1];

#          say "  Next anchor, ".$anchor_j_id.", has ".scalar(@$anchor_j_hits)." hits";
          # If this anchor is more than the max allowed skips away from the anchor currently at the end of
          # the path then we end the path and move on
          if($anchor_j_id - $path_end_id > $max_skips) {
   #         say "  Number of skipped anchors exceeds the limit, so ending the current path";
            last;
          }

          # For the hits, there could be multiple hits within the allowed distance, so we will have to loop through all
          # hits and then select the one that is closest to the expected distance
          my $path_candidate_hits_j = [];
          foreach my $hit_j (@{$anchor_j_hits}) {
            my $hit_j_region = ${$hit_j}[0];
     	      my $hit_j_start = ${$hit_j}[1];
            my $hit_j_end = ${$hit_j}[2];
            my $hit_j_strand = ${$hit_j}[3];

#            say "  Considering hit: ".$hit_j_region.":".$hit_j_start.":".$hit_j_end.":".$hit_j_strand;

            # First check is on region
            unless($hit_i_region eq $hit_j_region) {
  #            say "  Skipping due to different region";
              next;
            }

            # Second check is on strand, if there is an inversion it would be better to have a separate path for
            # the inverted region as opposed to trying to run the path across the entire region
            unless($hit_i_strand eq $hit_j_strand) {
 #             say "  Skipping due to different strand";
              next;
            }

            # Third check is on the distance and overlap, we want hits within the distance limit and also
            # not overlapping with the current path end. As the path is stranded, we need to check whether
            # the path is going in the correct direction
            if($hit_j_strand eq '+') {
              my $current_path_end = ${$path_end_hit}[2];
              # If the start of the hit is less than or equal to the end of the current path then either
              # the hit overlaps the path end, or the hit is before the path end and thus in the wrong place
              # If the dist is greater than the max dist then it also fails
              # So in either case just move on to the next hit
              my $dist = $hit_j_start - $current_path_end;
              if($dist <= 0 or $dist > $max_dist) {
#                say "  Skipping because of the distance, dist calculated as ".$dist.", must be > 0 and < ".$max_dist;
                next;
              }
              push(@$hit_j,$dist);
            } else {
              # In the reverse orientation use the start of the last path hit as the end coord of the path
              my $current_path_end = ${$path_end_hit}[1];
               my $dist = $current_path_end - $hit_j_end;
              if($dist <= 0 or $dist > $max_dist) {
 #               say "  Skipping because of the distance, dist calculated as ".$dist.", must be > 0 and < ".$max_dist;
                next;
              }
              push(@$hit_j,$dist);
            }

            # All checks are completed and it passes so this hit becomes the current path end, move onto next anchor
#            say "  Adding hit to the list of candidates for path extension";
            push(@$path_candidate_hits_j,$hit_j);
          } # foreach my $hit_j

          # If we're here there are path candidates, so we pick the one closest to the expected distance
          my $selected_candidate_hit;
          my $closest_dist = 9999;
          foreach my $candidate_hit (@$path_candidate_hits_j) {
            my $hit_dist = pop(@$candidate_hit);
            my $dist_difference = abs($hit_dist - $expected_dist);
            if($dist_difference < $closest_dist) {
              $selected_candidate_hit = $candidate_hit;
              $closest_dist = $dist_difference;
            }
          }

          if($selected_candidate_hit) {
#            say "  Adding selected hit to end of path";
            push(@$path_hit_i,$selected_candidate_hit);
            $path_end_hit = $selected_candidate_hit;
            $path_end_id = $anchor_j_id;
          }
        } # for my $j=$i+1
#        say "  Path length: ".scalar(@$path_hit_i);
        push(@$paths,$path_hit_i);
      } # foreach my $hit_i

    } # for my $i=0

    # At this stage we have all the paths, so select the longest path, then see what anchors it should cover. Then try to
    # find any remaining anchors in non-overlapping paths. Ideally we'd want to chain paths together based on proximity,
    # though given that some assemblies will be scaffold level there can't be too much emphasis on chaining

    # Put in some code that removes subpaths. Stringify the path and then remove anything that's a subset of another path
    say "  Have ".scalar(@$paths)." pre filtering";
    my $filtered_paths = [];
#    my $all_paths_string = stringify_paths($paths);

#    say "  All paths string:\n".$all_paths_string;

# I've commented the below out because later on we want to consider all subpaths to find out if there are subpaths
# that contain genes not present in the top path, while also not containing any genes that are already covered by
# the chosen paths. If you filter them out here then you may miss subpaths that have just the missing genes
# Another approach to this filtering would have been to just look at the region overlap and filter
# out paths that fully contained within another path in terms of overlap. That would be simpler, not sure if there
# would be any ramifications in terms of anchors though, this way is probably safe in terms of ensuring that all
# anchor patterns across a region are considered
#    foreach my $path (@$paths) {
#      my $path_string = stringify_paths([$path]);
 #     say "  Path string:\n".$path_string;
#      my @matches = $all_paths_string =~ /$path_string/g;
#      if(scalar(@matches) == 1) {
#        push(@$filtered_paths,$path);
#      }
#    }

    $filtered_paths = $paths;
    say "  Have ".scalar(@$filtered_paths)." post filtering";
    my @sorted_paths = sort { @$b <=> @$a } @{$filtered_paths};


    # Now sort paths based on the seq region start of the hit, at this point the path will in reverse order in terms of the
    # hits if it's on the opposite strand. Sorting it here is easier to deal with later in terms of the code
    foreach my $path (@sorted_paths) {
      my @sorted_path_hits = sort {$a->[1] <=> $b->[1]} (@{$path});
      $path = [@sorted_path_hits];
    }

    # Maybe just concat all path strings and anything that matches more than once is removed
    my ($region_start,$region_end) = @{$meta_clusters->{$cluster_id}->{'region_boundaries'}};
    say "  Original region boundaries: ".$region_start.":".$region_end;
    say "  Top path length: ".print_path_summary($sorted_paths[0]);
    if(scalar(@sorted_paths) > 1) {
      say "  2nd path length: ".print_path_summary($sorted_paths[1]);
    }

    $meta_clusters->{$cluster_id}->{'paths'} = [@sorted_paths];
#    say "  2nd path length: ".scalar(@{$sorted_paths[1]});
#    say "  3rd path length: ".scalar(@{$sorted_paths[2]});

  } # foreach my $cluster_id
}


sub print_path_summary {
  my ($path) = @_;

  my $start_hit = ${$path}[0];
  my $end_hit = ${$path}[-1];

  my $summary_string;
#  if(${$start_hit}[3] eq '+') {
    $summary_string = scalar(@$path)." ".${$start_hit}[0].":".${$start_hit}[1].":".${$end_hit}[2].":".${$start_hit}[3];
#  } else {
#    $summary_string = scalar(@$path)." ".${$start_hit}[0].":".${$end_hit}[1].":".${$start_hit}[2].":".${$start_hit}[3];
#  }

  return($summary_string);
}


sub stringify_paths {
  my ($paths) = @_;

  my @path_string_array = ();
  foreach my $path (@$paths) {
    my $path_string = "";
    foreach my $hit (@$path) {
      $path_string .= ":".${$hit}[0].":".${$hit}[1].":".${$hit}[2].":".${$hit}[3];
    }
    push(@path_string_array,$path_string);
  }

  # Just in terms of cleaning up the final string, there were some issues in terms of special characters
  # then messing up the regexes later. These substitions are overkill, I'm fairly sure it was just the +
  # causing problems, but the below cleans the strings to be a little more readable too if they need to
  # be looked at
  my $final_string = join("\n", @path_string_array);
  $final_string =~ s/\-/0/g;
  $final_string =~ s/\+/1/g;
  $final_string =~ s/\:/\_/g;
  $final_string =~ s/\_+/\_/g;
  return($final_string);
}


sub cluster_slice_genes {
  my ($genes) = @_;

  my $gene_clusters = [];
  my @sorted_genes = sort {$a->start <=> $b->start} (@{$genes});
  my $window_length = 100000;
  my $current_cluster_start = 0;
  my $current_cluster_end = 0;
  my $current_cluster_genes = [];
  my $last_gene;
  foreach my $gene (@sorted_genes) {
    unless($current_cluster_start) {
      $current_cluster_start = $gene->seq_region_start();
      $last_gene = $gene;
    }

    # If the start of the current gene is > the allowed window length, then the gene is outside the window and a new cluster should be started
    if(($gene->seq_region_start - $current_cluster_start + 1) > $window_length) {
      # At this point we know the gene is outside the window length, but if the current last gene of the cluster hangs over the edge of the window then
      # we also want to include any genes that are contained within the boundary of the overlapping gene, since we'll be aligning the region they reside
      # in anyway

      # If this gene is not wholly overlapping the last gene in the current cluster the start a new cluster, else it is so add to the current cluster
      if($last_gene->seq_region_end() < $gene->seq_region_end()) {
        push(@$gene_clusters,$current_cluster_genes);
        $current_cluster_genes = [];
        push($current_cluster_genes,$gene);
        $current_cluster_start = $gene->seq_region_start();
        $last_gene = $gene;
      } else {
        # Note last gene is not updated here to ensure the last gene remains the last overlapping gene of the current window
        push(@$current_cluster_genes,$gene);
      }
    } else {
      push(@$current_cluster_genes,$gene);
      $last_gene = $gene;
    }
  }

  if(scalar(@$current_cluster_genes)) {
    push(@$gene_clusters,$current_cluster_genes);
  }

  return($gene_clusters);
}


sub map_anchors {
  my ($meta_clusters,$input_file,$genome_index,$original_anchor_positions) = @_;

  my $paf_file = $input_file.".paf";
  my $minimap2_command = "minimap2 --cs --secondary=yes -x map-ont -N 20 ".$genome_index." ".$input_file." > ".$paf_file;
  system($minimap2_command);

  unless(-e $paf_file) {
    throw("Could not find paf file from running minimap on anchor seqs");
  }


  my $coverage_cutoff = 0.80;
  my $identity_cutoff = 0.98;

  my $perfect_mappings = 0;
  my $partial_mappings = 0;
  my $off_target_mappings = 0;
  my $unique_off_target_mappings = 0;
  my $unique_off_target_mappings_tracker = {};
  my $hits_by_anchor_id = {};

  open(IN,$paf_file);
  while(<IN>) {
    my $line = $_;
    chomp($line);
    my @eles = split("\t",$line);
    my $anchor_name = $eles[0];
    my $source_length = $eles[1];
    my $source_hit_start = $eles[2]+1;
    my $source_hit_end = $eles[3];
    my $target_strand = $eles[4];
    my $target_genomic_name = $eles[5];
    my $target_genomic_length = $eles[6];
    my $target_genomic_start = $eles[7]+1;
    my $target_genomic_end = $eles[8];
    my $alignment_identities = $eles[9];
    my $alignment_length = $eles[10];
    my $source_hit_midpoint = ceil($source_length/2);
    my $source_hit_length = $source_hit_end - $source_hit_start + 1;

    my $cov = ($source_hit_end-$source_hit_start)/$source_length;
    my $pid = $alignment_identities/($source_hit_end-$source_hit_start);
    unless($cov >= $coverage_cutoff and $pid >= $identity_cutoff) {
      say "Warning: hit for anchor ".$anchor_name." failed coverage/identity cut-off, so will not keep: ".$cov.":".$pid;
      next;
    }

    my $orignal_location_details = $original_anchor_positions->{$anchor_name};
    my ($original_seq_region_name,$original_seq_region_start,$original_seq_region_end) = @{$orignal_location_details};
    say "Anchor: ".$anchor_name." ".$original_seq_region_start.":".$original_seq_region_end;
    if($original_seq_region_name eq $target_genomic_name and $target_genomic_start == $original_seq_region_start and $original_seq_region_end == $target_genomic_end) {
      $perfect_mappings++;
    } elsif($original_seq_region_name eq $target_genomic_name and $target_genomic_start >= $original_seq_region_start and $target_genomic_end <= $original_seq_region_end) {
      $partial_mappings++;
    } else {
      unless($unique_off_target_mappings_tracker->{$anchor_name}) {
        $unique_off_target_mappings++;
        $unique_off_target_mappings_tracker->{$anchor_name} = 1;
      }
      $off_target_mappings++;
    }
    say $line;

    unless($hits_by_anchor_id->{$anchor_name}) {
      $hits_by_anchor_id->{$anchor_name} = [];
    }

    my $result = [$target_genomic_name,$target_genomic_start,$target_genomic_end,$target_strand,$anchor_name];
    push(@{$hits_by_anchor_id->{$anchor_name}},$result);
  }
  close IN;

  my $off_target_percent = ($unique_off_target_mappings/scalar(keys(%$original_anchor_positions))) * 100;
  $off_target_percent = sprintf("%.2f",$off_target_percent);

  say "Found ".$perfect_mappings." perfectly mapped anchors out of a total of ".scalar(keys(%$original_anchor_positions)." anchors");
  say "Found ".$partial_mappings." partially mapped anchors";
  say "Found ".$off_target_mappings." off target mappings";
  say "Found ".$unique_off_target_mappings." anchors with off target mappings, ".$off_target_percent."% of total anchor count";

  foreach my $cluster_id (keys(%$meta_clusters)) {
    my $anchors = $meta_clusters->{$cluster_id}->{'anchors'};
    my $missing_anchors = [];
    foreach my $anchor (@$anchors) {
      my $anchor_id = ${$anchor}[0];
      my $results = $hits_by_anchor_id->{$anchor_id};
      push(@$anchor,$results);
      unless(defined($results) and scalar(@$results)) {
        push(@$missing_anchors,$anchor_id);
        say "Warning: anchor with id ".$anchor_id." was not mapped in target";
        say "Seq:\n".${$anchor}[4];
      }
    } # end foreach my $anchor
    $meta_clusters->{$cluster_id}->{'missing_anchors'} = $missing_anchors;
  } # end foreach my $cluster_id
}



sub set_region_anchors {
  my ($meta_clusters,$original_anchor_positions,$region_seq,$slice) = @_;

  # At 1000 anchor size, 3.42% of anchors have at least one off-target hit
  # At 2000 anchor size, 3.09% of anchors have at least one off-target hit
  # At 5000 anchor size, 2.70% of anchors have at least one off-target hit and one anchor is not found anymore
  # It's hard to tell in this control scenario whether reducing the percent having an off-target is better
  # than having more anchors in general. More anchors means a finer grained evaluation of the mappability, less
  # off targets means more signal and a simpler set of paths to evaluate
  my $anchor_size = 500;
  my $anchor_spacing = 250;
  my $anchor_increment = $anchor_size+$anchor_spacing;
  my $max_n_ratio = 0.05;
  my $anchor_id = 1;
  foreach my $cluster_id (keys(%{$meta_clusters})) {
    my ($region_start,$region_end) = @{$meta_clusters->{$cluster_id}->{'region_boundaries'}};
    say "Found RSE: ".$region_start.":".$region_end;

    # Anchors each will have [anchor_id,region_start,region_end,seq]
    my $anchors = [];
    # Starting at the region start, add an anchor of the size specified if the anchor start <= the end of the region and the end of the
    # anchor is < the slice end and the anchor contains < 5 percent Ns. Remove leading/lagging Ns. Then jump the spacing dist and repeat
    # substr is 0-index based, hence the $region_start-1 and then $i+1 to get the $anchor_start correctly
    for(my $i=$region_start-1; $i < $region_end; $i += $anchor_increment) {
      # In this case the anchor is off the slice, so just finish
      if($i+$anchor_size > $slice->length()) {
        last;
      }

      my $seq = substr($region_seq,$i,$anchor_size);
      my $n_count = $seq =~ tr/N//;
      my $n_ratio = $n_count/$anchor_size;
      if($n_ratio > $max_n_ratio) {
        say "Warning: skipping anchor because it fails the N ratio test";
        next;
      }

      $seq =~ s/^N+//;
      $seq =~ s/N+$//;
      if (length($seq) < $anchor_size) {
        say 'smaller anchor ', $anchor_id, ' ', $anchor_size, ' > ', length($seq);
      }
#      say "Anchor ".$region_start." ".$region_end." ".$anchor_id.":\n".$seq;
      my $anchor_start = $i+1;
      my $anchor_end = $i+$anchor_size;
      my $anchor = [$anchor_id,$slice->seq_region_name(),$anchor_start,$anchor_end,$seq];
      push(@$anchors,$anchor);
      $original_anchor_positions->{$anchor_id} = [$slice->seq_region_name(),$anchor_start,$anchor_end];
      $anchor_id++;
    }
    $meta_clusters->{$cluster_id}->{'anchors'} = $anchors;
  }
}


sub generate_fasta_record {
  my ($meta_clusters) = @_;

  my $fasta_record = "";
  foreach my $cluster_id (keys(%{$meta_clusters})) {
    my $anchors = $meta_clusters->{$cluster_id}->{'anchors'};
    foreach my $anchor (@$anchors) {
       my $anchor_id = ${$anchor}[0];
       my $seq = ${$anchor}[-1];
       $fasta_record .= ">".$anchor_id."\n".$seq."\n";
    }
  }
  return($fasta_record);
}

sub set_cluster_regions {
  my ($gene_clusters,$meta_clusters,$slice) = @_;

  my $minimum_region_length = 100000;
  my $minimum_flanking_region = 5000;
  my $cluster_id = 1;
  foreach my $gene_cluster (@$gene_clusters) {
    my $region_start;
    my $region_end;
    foreach my $gene (@$gene_cluster) {
      unless($region_start) {
        $region_start = $gene->seq_region_start();
        $region_end = $gene->seq_region_end();
        next;
      }

      if($gene->seq_region_start() < $region_start) {
        $region_start = $gene->seq_region_start();
      }

      if($gene->seq_region_end() > $region_end) {
        $region_end = $gene->seq_region_end();
      }
    }

#    my $region_start = ${$gene_cluster}[0]->seq_region_start();
#    my $region_end = ${$gene_cluster}[-1]->seq_region_end();
    my $initial_length = $region_end - $region_start + 1;
    if($initial_length < $minimum_region_length) {
      my $diff = ceil(($minimum_region_length - $initial_length)/2);
      $region_start -= $diff + $minimum_flanking_region;
      $region_end += $diff + $minimum_flanking_region;
    } else {
      $region_start -= $minimum_flanking_region;
      $region_end += $minimum_flanking_region;
    }

    if($region_start < 1) {
      $region_start = 1;
    }

    if($region_end > $slice->length()) {
      $region_end = $slice->length();
    }

    print_clusters([$gene_cluster],1);
    say "Adjusted regions: ".$region_start.":".$region_end.":".($region_end-$region_start+1);

    $meta_clusters->{$cluster_id}->{'gene_cluster'} = $gene_cluster;
    $meta_clusters->{$cluster_id}->{'region_boundaries'} = [$region_start,$region_end];
    $cluster_id++;
  }

#  my $record = "";
#  my $current_end = $start_offset;
#  my $end_offset = 5000000;
#  for(my $i=1; $i<=$num_regions && $current_end < ($slice->length()-$end_offset); $i++) {
#    my $current_start = ($i * $interval_length) + $start_offset;
#    my $current_end = $current_start + $region_length;
#    say "Region: ".$slice->seq_region_name.":".$current_start.":".$current_end;
#    my $header = ">".$header_label."_".($i-1);
#    my $seq = $slice->subseq($current_start,$current_end);
#    if($seq =~ /N/) {
#      $num_regions++;
#      next;
#    }
#    $record .= $header."\n".$seq."\n";
#  }
#  return($record);
}



sub print_clusters {
  my ($clusters,$skip_gene_string) = @_;

  my $total_cluster_genes = 0;
  foreach my $cluster (@$clusters) {
    $total_cluster_genes += scalar(@$cluster);
    my $cluster_region = ${$cluster}[0]->seq_region_name();
    my $cluster_start = ${$cluster}[0]->seq_region_start();
    my $cluster_end = ${$cluster}[-1]->seq_region_end();
    my $cluster_length = $cluster_end -$cluster_start + 1;
    say "Cluster: ".$cluster_region.": ".scalar(@$cluster)." genes, ".$cluster_start.":".$cluster_end.":".$cluster_length;
    my $gene_string = "";
    foreach my $gene (@$cluster) {
      my $string = feature_string($gene);
      $gene_string .= $string.",";
    }
    $gene_string =~ s/\,$//;
    unless($skip_gene_string) {
      say $gene_string;
    }
  }

  say "Total genes in all clusters: ".$total_cluster_genes;
}


sub feature_string {
  my ($feature) = @_;

  my $string = "(".$feature->seq_region_start().":".$feature->seq_region_end().":".$feature->strand().")";
  return($string);
}


sub write_input_file {
  my ($fasta_records) = @_;

  my $output_file = create_file_name();
  open(OUT,">".$output_file);
  say OUT $fasta_records;
  close OUT;

  return($output_file);
}

