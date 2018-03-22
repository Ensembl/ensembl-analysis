# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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


package Bio::EnsEMBL::Analysis::RunnableDB::Finished::ClusterDepthFilter;

use warnings ;
use strict;
use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Pipeline::Tools::MM_Taxonomy;

use Data::Dumper;
#use Scalar::Util 'weaken';

use base 'Bio::EnsEMBL::Analysis::RunnableDB::Finished::DepthFilter';

my $DEBUG                = 0;
my $SANITY_CHECK_STRANDS = 0;

# NB: not a method
sub overlap {

    # check if feature f1 overlaps feature f2

    my ( $f1_start, $f1_end, $f2_start, $f2_end ) = @_;
    
    return ($f1_end >= $f2_start and $f1_start <= $f2_end);
}

# NB: not a method
sub is_consistent {

	# check if the given hit can join the cluster

	my $hit     = shift;
	my $cluster = shift;

	# check if any of the new hit's exons encroach on any
	# of the cluster's introns

	for my $af ( @{ $hit->{features} } ) {
		for my $intron ( @{ $cluster->{introns} } ) {
			if ( overlap( $af->start, $af->end, $intron->[0], $intron->[1] ) ) {
				return 0;
			}
		}
	}

	# check if any of the clusters's exons encroach on any
	# of the new hit's introns

	my $start = $cluster->{start};

	for my $cluster_intron ( @{ $cluster->{introns} } ) {

		my $end = $cluster_intron->[0]-1;

		for my $intron ( @{ $hit->{introns} } ) {
			if ( overlap( $start, $end, $intron->[0], $intron->[1] ) ) {
				return 0;
			}
		}

		$start = $cluster_intron->[1]+1;
	}

	my $end = $cluster->{end};

	for my $intron ( @{ $hit->{introns} } ) {
		if ( overlap( $start, $end, $intron->[0], $intron->[1] ) ) {
			return 0;
		}
	}

	return 1;
}

# NB: not a method
sub find_new_introns {

	# identify and return any introns from set2 that are not in set1

	my ( $set1, $set2 ) = @_;

	my %h1 = map { $_->[0] . "-" . $_->[1] => $_ } @$set1;
	my %h2 = map { $_->[0] . "-" . $_->[1] => $_ } @$set2;

	my @new_keys = grep { !( exists $h1{$_} ) } keys %h2;

	return [ map { $h2{$_} } @new_keys ];
}

# NB: not a method
sub hit_str {
	my $hit = shift;
	my $s   = '';
	map { $s .= $_->start . "-" . $_->end . "," } @{ $hit->{features} };
	$s .= " (".$hit->{taxon_id}.")";
	return $s;
}

# NB: not a method
sub cluster_str {
	my $cluster = shift;
	my $s       = $cluster->{start} . "-";
	map { $s .= $_->[0] . "," . $_->[1] . "-" } @{ $cluster->{introns} };
	$s .= $cluster->{end};
	return $s;
}

# NB: not a method
sub meta_cluster_str {
	my $meta = shift;
	my $s = sprintf "%d-%d (%d)", $meta->{start}, $meta->{end}, $meta->{strand};
	return $s;
}

sub group_features_by_name {
	my ( $self, $orig_features, $percentid_cutoff, $no_filter ) = @_;

	my %grouped_by_name = ();

	for my $af (@$orig_features) {
		my ( $score, $percentid ) = ( $af->score(), $af->percent_id() );
		if ( $percentid < $percentid_cutoff ) {
			next;
		}

		my $hit = $grouped_by_name{ $af->hseqname() } ||= {};

		push @{ $hit->{features} }, $af;

		if ( $no_filter && $self->get_hit_description( $af->hseqname() ) ) {
			$hit->{taxon_id} =
			  $self->get_hit_description( $af->hseqname() )->taxon_id;
		}
	}

	return \%grouped_by_name;
}

sub break_discontinuities {

	my ( $self, $hits_by_name, $no_filter, $slice ) = @_;

	my @hits = ();

	for my $hit ( values %$hits_by_name ) {

            print STDERR "ClusterDepthFilter: Processing hit: ", $hit->{features}->[0]->hseqname,
            "\n"
                if $DEBUG;

            my @remaining_features = sort { $a->hstrand <=> $b->hstrand 
                                                        ||
                                            $a->start   <=> $b->start
                                                        ||
                                            $a->hstart  <=> $b->hstart } @{ $hit->{features} };

            while (@remaining_features) {

		my @features = @remaining_features;
                @remaining_features = ();

		my $last = shift @features;

		my $new_hit = {};
		$new_hit->{features}      = [$last];
		$new_hit->{tot_score}     = $last->score;
		$new_hit->{avg_percentid} = $last->percent_id;
		$new_hit->{strand}        = $last->hstrand;
		$new_hit->{taxon_id}      = $hit->{taxon_id} if $no_filter;

                print STDERR "ClusterDepthFilter: starting new hit: hstrand: ",
                $last->hstrand, " af: ", $last->hstart, "-", $last->hend,
                " (", $last->start, "-", $last->end, ")\n"
                    if $DEBUG;

              FEAT: for my $af ( @features ) {

			if (    $last->hstrand == $af->hstrand
                             && $last->hstart  == $af->hstart
                             && $last->hend    == $af->hend
                             && $last->start   == $af->start
                             && $last->end     == $af->end      ) {

				print STDERR "ClusterDepthFilter: Found duplicate align "
				  . "feature for hit " . $af->hseqname
				  . " [coords: " . $af->hstart . "-" . $af->hend
                                  . " (" . $af->start . "-" . $af->end . ")"
				  . ", slice: "
				  . $slice->name . "]\n";

				next FEAT;
			}

			my $delta =
			  		  $af->hstrand == $af->strand
			  		? ( $af->hstart - $last->hend ) 
			  		: ( $last->hstart - $af->hend );

                        if ( $af->hstrand == $last->hstrand && $af->start <= $last->end ) {

                                # found an overlap - bump for reprocessing
				print STDERR "ClusterDepthFilter: found overlap: hstrand: ",
				  $af->hstrand, " last: ", $last->hstart, "-", $last->hend,
				  " (", $last->start, "-", $last->end, ") af: ", $af->hstart,
				  "-", $af->hend, " (", $af->start, "-", $af->end, ")\n"
				  if $DEBUG;

                                push @remaining_features, $af;
                                next FEAT;
                        }
			elsif ( $last->hstrand != $af->hstrand || $delta != 1 ) {

				# found a discontinuity - it might belong to stuff on remaining features

				print STDERR "ClusterDepthFilter: discontinuity: hstrand: ",
				  $af->hstrand, " last: ", $last->hstart, "-", $last->hend,
				  " (", $last->start, "-", $last->end, ") af: ", $af->hstart,
				  "-", $af->hend, " (", $af->start, "-", $af->end, ")\n"
				  if $DEBUG;

                                push @remaining_features, $af;
                                next FEAT;
			}
			else {
				print STDERR "ClusterDepthFilter: continuing hit: hstrand: ",
				  $af->hstrand, " last: ", $last->hstart, "-", $last->hend,
				  " (", $last->start, "-", $last->end, ") af: ", $af->hstart,
				  "-", $af->hend, " (", $af->start, "-", $af->end, ")\n"
				  if $DEBUG;

				# continue the same hit

				push @{ $new_hit->{features} }, $af;
				$new_hit->{tot_score}     += $af->score;
				$new_hit->{avg_percentid} += $af->percent_id;
			}

			$last = $af;
		}

		# add the hit to the list

		$new_hit->{avg_percentid} /= scalar( @{ $new_hit->{features} } );
		push @hits, $new_hit;

                if (@remaining_features) {
                    print STDERR "ClusterDepthFilter: rescanning\n" if $DEBUG;
                }

            } # @remaining_features 
	}

	return \@hits;
}

sub cluster_hits {

	my ( $self, $hits ) = @_;

	# sort the hits so they are arranged sequentially against the
	# genomic sequence, with the longest sequence first for hits
	# with identical starts (this allows us to find all overlaps in
	# one pass)

	my @hits = sort {
		     ( $a->{features}->[0]->start <=> $b->{features}->[0]->start )
		  || ( $b->{features}->[-1]->end <=> $a->{features}->[-1]->end )
	} @$hits;

	# try to cluster hits according to overlaps and introns

	my @clusters = ();

	for my $hit (@hits) {

		my @afs = @{ $hit->{features} };

		$hit->{start} = $afs[0]->start;
		$hit->{end}   = $afs[-1]->end;

		print STDERR "ClusterDepthFilter: looking at hit: ", hit_str($hit), "\n"
		  if $DEBUG;

		my @introns = ();

		if ( @afs > 1 ) {    # we have at least 1 intron

			# identify all the introns

			my $start = $afs[0]->end;
			for my $af ( @afs[ 1 .. $#afs - 1 ] ) {
				my $end = $af->start;
				push @introns, [ $start+1, $end-1 ];
				$start = $af->end;
			}
			my $end = $afs[-1]->start;
			push @introns, [ $start+1, $end-1 ];
		}

		$hit->{introns} = [@introns];

		# check if we can add this hit to any existing cluster

		my $added = 0;

		for my $cluster (@clusters) {
			if (   $hit->{strand} == $cluster->{strand}
				&& $hit->{start} >= $cluster->{start}
				&& $hit->{end} <= $cluster->{end}
				&& is_consistent( $hit, $cluster ) )
			{

				# this hit is subsumed by the cluster and its
				# introns are consistent, so we can add it

				print STDERR "ClusterDepthFilter: adding hit to cluster: ",
				  cluster_str($cluster), "\n"
				  if $DEBUG;

				push @{ $cluster->{hits} }, $hit;

				$added = 1;

				last;
			}
			else {
				print STDERR "ClusterDepthFilter: can't add hit to cluster: ",
				  cluster_str($cluster), "\n"
				  if $DEBUG;
			}
		}

		# if we didn't find a matching cluster, start a new one
		if ( !$added ) {

			my $cluster = {
				start   => $hit->{start},
				end     => $hit->{end},
				introns => [ map { [@$_] } @introns ],
				hits    => [$hit],
				strand  => $hit->{strand}
			};
			push @clusters, $cluster;

			print STDERR "ClusterDepthFilter: creating new cluster: ",
			  cluster_str($cluster), "\n"
			  if $DEBUG;
		}
	}

	return \@clusters;
}

sub compute_vulgar_strings {
	my $hits = shift;

	for my $hit (@$hits) {
		
		my @afs = @{ $hit->{features} };
		
		my $vs = '';
		
		my $last;
		
		for my $af (@afs) {
			if ($last) {
				$vs .= '5 0 2 ';
				$vs .= 'I 0 '.($af->start - $last - 4).' ';
				$vs .= '3 0 2 ';
			}
			$vs .= 'M '.$af->length.' '.$af->length.' ';
			$last = $af->end;
		}
		
		$hit->{vulgar_string} = $vs;
		
		#print STDERR hit_str($hit), "\n\n";
		#print STDERR $vs, "\n\n";
	}
}

sub is_consistent_by_vulgar_string {
	my ($hit, $cluster) = @_;
	
	return 1 unless ($hit->{vulgar_string} =~ /I/ || $cluster->{vulgar_string} =~ /I/);
	
	my $vs = $hit->{vulgar_string};
	$vs =~ s/^M \d+ \d+ 5 0 2 //;
	$vs =~ s/ 3 0 2 M \d+ \d+ $//;
	
	if ($cluster->{vulgar_string} =~ /$vs/) {
		if ($vs ne $hit->{vulgar_string}) {
			print STDERR "orig   : ", $hit->{vulgar_string}, "\n";
			print STDERR "trimmed: ", $vs, "\n";
		}
		print STDERR "MATCH:\n";
		print STDERR $vs, "\n";
		print STDERR $cluster->{vulgar_string}, "\n";
		return 1;	
	}
	else {
		print STDERR "NO MATCH:\n";
		print STDERR $vs, "\n";
		print STDERR $cluster->{vulgar_string}, "\n";
		return 0;
	}
}

sub cluster_hits_by_strings {
	
	my $self = shift;
	
	my $hits = shift;
	
	compute_vulgar_strings($hits);
	
	my @hits = sort {
		     ( $a->{features}->[0]->start <=> $b->{features}->[0]->start )
		  || ( $b->{features}->[-1]->end <=> $a->{features}->[-1]->end )
	} @$hits;
	
	my @clusters = ();
	
	for my $hit (@hits) {
		
		my $added = 0;
		
		for my $cluster (@clusters) {
			
			if (   $hit->{strand} == $cluster->{strand}
				&& $hit->{start} >= $cluster->{start}
				&& $hit->{end} <= $cluster->{end}
				&& is_consistent_by_vulgar_string( $hit, $cluster ) )
			{
				push @{ $cluster->{hits} }, $hit;

				$added = 1;

				last;
			}
		}
		
		# if we didn't find a matching cluster, start a new one
		if ( !$added ) {

			my $cluster = {
				start   => $hit->{start},
				end     => $hit->{end},
				vulgar_string => $hit->{vulgar_string},
				hits    => [$hit],
				strand  => $hit->{strand}
			};
			push @clusters, $cluster;

			print STDERR "ClusterDepthFilter: creating new cluster: ",
			  cluster_str($cluster), "\n"
			  if $DEBUG;
		}
	}
	
	print "Created ",scalar(@clusters), " clusters using vulgar strings\n";
	
	return @clusters;
}

sub cluster_clusters {

	# cluster the clusters by overlaps

	my ( $self, $clusters ) = @_;

	my @meta_clusters = ();

	for my $cluster (@$clusters) {

		my $added = 0;

		for my $meta (@meta_clusters) {
			if (
				overlap(
					$cluster->{start}, $cluster->{end},
					$meta->{start},    $meta->{end}
				)
				&& $cluster->{strand} == $meta->{strand}
			  )
			{
				push @{ $meta->{clusters} }, $cluster;
				$added = 1;
				last;
			}
		}

		if ( !$added ) {
			my $new_meta = {
				start      => $cluster->{start},
				end        => $cluster->{end},
				strand     => $cluster->{strand},
				clusters   => [$cluster],
				hits_added => 0
			};

			push @meta_clusters, $new_meta;
		}
	}

	return \@meta_clusters;
}

sub build_summary_features {
	my ( $self, $meta_clusters ) = @_;

	my @summary_features = ();

	my $tot_filtered = 0;

	for my $meta_cluster (@$meta_clusters) {

		for my $cluster ( @{ $meta_cluster->{clusters} } ) {

			# create new 'summary' features representing the maximal
			# coverage of the remaining features (if we have any)

			# establish the start, end, introns list and score of
			# the summary feature

			my ( $start, $end, $score, $introns );

			my $filtered = 0;

			for my $hit ( @{ $cluster->{hits} } ) {

				# ignore hits that have already been added
				if ( $hit->{features}->[0]->analysis == $self->analysis ) {
					print STDERR
					  "ClusterDepthFilter: ignoring previously added feature\n"
					  if $DEBUG;
					next;
				}

				if ( defined $start ) {

					# check if this hit extends the summary feature

					my $s = $hit->{features}->[0]->start;
					$start = $s if $s < $start;

					my $e = $hit->{features}->[-1]->end;
					$end = $e if $e > $end;

					my $ms = $hit->{tot_score};
					$score = $ms if $ms > $score;

					my $new_introns =
					  find_new_introns( $introns, $hit->{introns} );
					push @$introns, map { [@$_] } @$new_introns;
					$introns = [ sort { $a->[0] <=> $b->[0] } @$introns ];
				}
				else {

					# initialise the summary feature data

					$start   = $hit->{features}->[0]->start;
					$end     = $hit->{features}->[-1]->end;
					$score   = $hit->{tot_score};
					$introns = $hit->{introns};
				}

				$filtered++;
			}

			if ( defined $start ) {

				# create our summary features with the maximal set of
				# introns

				for my $intron (@$introns) {
					my $summary = Bio::EnsEMBL::SimpleFeature->new(
						-start  => $start,
						-end    => $intron->[0],
						-strand => $cluster->{strand},
						-score  => $score,
						-display_label =>
						  sprintf( "summary feature from %d redundant hits",
							$filtered ),
						-analysis => $self->analysis()
					);

					push @summary_features, $summary;

					$start = $intron->[1];
				}

				my $summary = Bio::EnsEMBL::SimpleFeature->new(
					-start         => $start,
					-end           => $end,
					-strand        => $cluster->{strand},
					-score         => $score,
					-display_label =>
					  sprintf( "summary feature from %d redundant hits",
						$filtered ),
					-analysis => $self->analysis()
				);
				push @summary_features, $summary;
			}

			$tot_filtered += $filtered;
		}
	}

	print STDERR
"ClusterDepthFilter: built summary features from $tot_filtered filtered hits\n";

	return \@summary_features;
}

sub filter_features {

	my ( $self, $meta_clusters, $hits_to_keep, $no_filter, $matching_taxa,
		$max_hits_per_meta_cluster )
	  = @_;

	my @retained_features = ();
	my @summary_features  = ();

	my $retained = 0;

	for my $meta_cluster (@$meta_clusters) {

		my @clusters = @{ $meta_cluster->{clusters} };

		for my $cluster (@clusters) {

			# sort the hits in this cluster by tot_score, breaking ties
			# according to avg_percentid

			$cluster->{hits} = [
				sort {
					     ( $b->{tot_score} <=> $a->{tot_score} )
					  || ( $b->{avg_percentid} <=> $a->{avg_percentid} )
				  } @{ $cluster->{hits} }
			];

			# and add the features of the $hits_to_keep best hits to
			# the list of filtered features

			my $added = 0;

			print STDERR "ClusterDepthFilter: adding features of cluster: ",
			  $cluster->{start}, "-", $cluster->{end}, "\n"
			  if $DEBUG;

			while ( $added < $hits_to_keep && @{ $cluster->{hits} } ) {

				my $hit_to_keep = shift @{ $cluster->{hits} };

				$added++;

				print STDERR "ClusterDepthFilter: adding features of hit: ",
				  $hit_to_keep->{start}, "-", $hit_to_keep->{end}, "\n"
				  if $DEBUG;

				for my $af ( @{ $hit_to_keep->{features} } ) {
					print STDERR "ClusterDepthFilter: adding feature: ",
					  $af->start, "-", $af->end, "\n"
					  if $DEBUG;

					$af->analysis( $self->analysis );
					$af->dbID(0);
					$af->{adaptor} = undef;
					push @retained_features, $af;
				}
			}

			$meta_cluster->{hits_added} += $added;
		}

		print STDERR "ClusterDepthFilter: meta-cluster hits added: ".
			$meta_cluster->{hits_added}."\n" if $DEBUG;

		if (   $no_filter
			&& $max_hits_per_meta_cluster
			&& ( $meta_cluster->{hits_added} < $max_hits_per_meta_cluster ) )
		{

			# find candidate hits to retain from hits from matching taxa

			my @candidates = ();

			for my $cluster (@clusters) {
				for my $hit ( @{ $cluster->{hits} } ) {
					print STDERR "ClusterDepthFilter: candidate hit taxon_id: ".hit_str($hit)." ".$hit->{taxon_id}."\n" if $DEBUG;
					if ( grep { /^$hit->{taxon_id}$/ } @$matching_taxa ) {
						push @candidates, $hit;
						print STDERR "ClusterDepthFilter: candidate hit: ".hit_str($hit)."\n" if $DEBUG;
					}
				}
			}

			# sort these by score & percent_id

			@candidates = sort {
				     ( $b->{tot_score} <=> $a->{tot_score} )
				  || ( $b->{avg_percentid} <=> $a->{avg_percentid} )
			} @candidates;

			# and add as many as possible

			while ($meta_cluster->{hits_added} < $max_hits_per_meta_cluster
				&& @candidates )
			{

				my $hit_to_keep = shift @candidates;

				$meta_cluster->{hits_added}++;
				$retained++;

				for my $af ( @{ $hit_to_keep->{features} } ) {
					print STDERR "ClusterDepthFilter: adding extra feature: ",
					  $af->start, "-", $af->end, "\n"
					  if $DEBUG;

					$af->analysis( $self->analysis );
					$af->dbID(0);
					$af->{adaptor} = undef;
					push @retained_features, $af;
				}
			}
		}
	}

	print STDERR
	  "ClusterDepthFilter: $retained features retained from matching taxa\n";

	my $summary_features = $self->build_summary_features($meta_clusters);

	return ( \@retained_features, $summary_features );
}

sub _sneaky_new {
	my ($class) = @_;
  	my $self = bless {},$class;
  	return $self; 
}

sub depth_filter {

	my $self = shift;

	my ( $orig_features, $slice, $max_coverage, $percentid_cutoff,
		$hits_to_keep, $no_filter, $max_hits_per_meta_cluster, $hit_db )
	  = @_;

	print STDERR "ClusterDepthFilter: HitsToKeep=$hits_to_keep\n";
	print STDERR "ClusterDepthFilter: PercentIdCutoff=$percentid_cutoff\n";
	print STDERR "ClusterDepthFilter: NoFilter=$no_filter\n" if $no_filter;
	print STDERR
	  "ClusterDepthFilter: MaxHitsPerMetaCluster=$max_hits_per_meta_cluster\n";
	print STDERR "ClusterDepthFilter: "
	  . scalar(@$orig_features)
	  . " features before filtering\n";

	# features from these taxa should be retained

	my @matching_taxa = ();

	# (though only fill in this list if we'll actually use it)

	if ($no_filter) {
		@matching_taxa = @{ $self->get_taxon_id_child($no_filter) };
		push @matching_taxa, $no_filter;
	}

	my $grouped_by_name =
	  $self->group_features_by_name( $orig_features, $percentid_cutoff,
		$no_filter );

	print STDERR "ClusterDepthFilter: "
	  . scalar( keys %$grouped_by_name )
	  . " unique hitnames\n";

	# first break apart any hits that contain discontinuous hits

	my $hits =
	  $self->break_discontinuities( $grouped_by_name, $no_filter, $slice );

	if ($SANITY_CHECK_STRANDS) {
		for my $hit (@$hits) {
			for my $af ( @{ $hit->{features} } ) {
				if ( $af->hstrand != $hit->{strand} ) {
                                        my $af_name  = $af->hseqname;
					die "mismatching strands in hit group '$af_name'\n";
				}
			}
		}
	}

	print STDERR "ClusterDepthFilter: "
	  . scalar(@$hits)
	  . " hits after processing discontinuities\n";

	# identify clusters of hits

	my $clusters = $self->cluster_hits($hits);
	
	#my $clusters2 = $self->cluster_hits_by_strings($hits);

	print STDERR "ClusterDepthFilter: grouped "
	  . scalar(@$hits)
	  . " hits into "
	  . scalar(@$clusters)
	  . " clusters";

	print "\n\n";

	#die;

	if ($SANITY_CHECK_STRANDS) {
		for my $cluster ( @$clusters ) {
			for my $hit ( @{ $cluster->{hits} } ) {
				for my $af ( @{ $hit->{features} } ) {
					if ( $af->hstrand != $cluster->{strand} ) {
						die "mismatching strands in clusters\n";
					}
				}
			}
		}
	}

	print STDERR sprintf( " (average of %.2f hits per cluster)",
		( scalar(@$hits) / scalar(@$clusters) ) )
	  if @$clusters > 0;

	print STDERR "\n";

	map { print STDERR cluster_str($_), "\n" } @$clusters if $DEBUG;

	# compute the meta-clusters

	my $meta_clusters = $self->cluster_clusters($clusters);

	print STDERR "ClusterDepthFilter: grouped "
	  . scalar(@$clusters)
	  . " clusters into "
	  . scalar(@$meta_clusters)
	  . " meta clusters\n";

	if ($SANITY_CHECK_STRANDS) {
		for my $meta (@$meta_clusters) {
			for my $cluster ( @{ $meta->{clusters} } ) {
				for my $hit ( @{ $cluster->{hits} } ) {
					for my $af ( @{ $hit->{features} } ) {
						if ( $af->hstrand != $meta->{strand} ) {
							die "mismatching strands in meta clusters";
						}
					}
				}
			}
		}
	}

	# now actually identify the features we want to include

	my ( $filtered_features, $summary_features ) =
	  $self->filter_features( $meta_clusters, $hits_to_keep, $no_filter,
		\@matching_taxa, $max_hits_per_meta_cluster );

	print STDERR "ClusterDepthFilter: "
	  . scalar(@$summary_features)
	  . " summary features\n";

	print STDERR "ClusterDepthFilter: "
	  . scalar(@$filtered_features)
	  . " features after filtering\n";

	return ( $filtered_features, $summary_features );
}

1;

=head1 NAME - Bio::EnsEMBL::Analysis::RunnableDB::Finished::ClusterDepthFilter

=head2 AUTHOR

Graham Ritchie B<email> gr5@sanger.ac.uk
