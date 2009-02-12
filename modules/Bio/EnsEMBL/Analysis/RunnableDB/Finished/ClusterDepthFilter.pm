### Bio::EnsEMBL::Analysis::RunnableDB::ClusterDepthFilter

package Bio::EnsEMBL::Analysis::RunnableDB::Finished::ClusterDepthFilter;

use strict;
use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);

use base 'Bio::EnsEMBL::Analysis::RunnableDB::Finished::DepthFilter';

sub overlap {

	# check if feature f1 overlaps feature f2

	my ( $f1_start, $f1_end, $f2_start, $f2_end ) = @_;

	return ( ( $f1_start > $f2_start  && $f1_start < $f2_end )
		  || ( $f1_end > $f2_start    && $f1_end < $f2_end )
		  || ( $f1_start <= $f2_start && $f1_end >= $f2_end ) );
}

sub is_consistent {

	# check if the given node can join the cluster

	my $node    = shift;
	my $cluster = shift;

	# check if any of the new node's exons encroach on any
	# of the cluster's introns

	for my $af ( @{ $node->{features} } ) {
		for my $intron ( @{ $cluster->{introns} } ) {
			if ( overlap( $af->start, $af->end, $intron->[0], $intron->[1] ) ) {
				return 0;
			}
		}
	}

	# check if any of the clusters's exons encroach on any
	# of the new node's introns

	my $start = $cluster->{start};

	for my $cluster_intron ( @{ $cluster->{introns} } ) {

		my $end = $cluster_intron->[0];

		for my $intron ( @{ $node->{introns} } ) {
			if ( overlap( $start, $end, $intron->[0], $intron->[1] ) ) {
				return 0;
			}
		}

		$start = $cluster_intron->[1];
	}

	my $end = $cluster->{end};

	for my $intron ( @{ $node->{introns} } ) {
		if ( overlap( $start, $end, $intron->[0], $intron->[1] ) ) {
			return 0;
		}
	}

	return 1;
}

sub find_new_introns {

	# identify and return any introns from set2 that are not in set1

	my ( $set1, $set2 ) = @_;

	my %h1 = map { $_->[0] . "-" . $_->[1] => $_ } @$set1;
	my %h2 = map { $_->[0] . "-" . $_->[1] => $_ } @$set2;

	my @new_keys = grep { !( exists $h1{$_} ) } keys %h2;

	return [ map { $h2{$_} } @new_keys ];
}

sub node_str {
	my $node = shift;
	my $s = '';
	map {$s .= $_->start."-".$_->end.","} @{ $node->{features} };
	return $s;
}

sub cluster_str {
	my $cluster = shift;
	my $s = $cluster->{start}."-";
	map {$s .= $_->[0].",".$_->[1]."-"} @{ $cluster->{introns} };
	$s .= $cluster->{end};
	return $s;
}

sub depth_filter {
	
	my $DEBUG = 0;
	
	my $self = shift;
	
	my ( $orig_features, $slice, $max_coverage, $percentid_cutoff, $nodes_to_keep ) = @_;

	print STDERR "ClusterDepthFilter: NodesToKeep=$nodes_to_keep\n";
	print STDERR "ClusterDepthFilter: PercentIdCutoff=$percentid_cutoff\n";
	print STDERR "ClusterDepthFilter: "
	  . scalar(@$orig_features)
	  . " features before filtering\n";

	my %grouped_byname = ();

	for my $af (@$orig_features) {
		my ( $score, $percentid ) = ( $af->score(), $af->percent_id() );
		if ( $percentid < $percentid_cutoff ) {
			next;
		}

		my $node = $grouped_byname{ $af->hseqname() } ||= {};		
		push @{ $node->{features} }, $af;
	}

	print STDERR "ClusterDepthFilter: "
	  . scalar( keys %grouped_byname )
	  . " unique hitnames\n";

	# identify clusters of nodes

	my @nodes = ();	

	# first break apart any nodes that contain discontinuous hits

	for my $node (values %grouped_byname) {
		
		my @features = sort { $a->start <=> $b->start } @{ $node->{features} };
		
		my $last = $features[0];
			
		my $new_node = {};
		$new_node->{features} = [$last];
		$new_node->{tot_score} = $last->score;
		$new_node->{avg_percentid} = $last->percent_id;
		
		for my $af (@features[1 .. $#features]) {
			
			if ($last->hstart == $af->hstart && $last->hend == $af->hend) {
				throw("Found duplicate align feature for hit ".
					  $af->hseqname." (coords: ".$af->hstart."-".$af->hend.
					  ", slice: ".$slice->name.")");
			}
			
			my $delta = $af->hstrand == -1 ?
							($last->hstart - $af->hend) :
							($af->hstart - $last->hend);
			
			if ($delta != 1) {
				
				# found a discontinuity
				
				print STDERR "ClusterDepthFilter: creating new node: last->hend: ",
					$last->hend," af->hstart: ", $af->hstart, "\n" if $DEBUG;
				
				# we've finished the current node, so add it to the list 
				
				$new_node->{avg_percentid} /= scalar( @{ $new_node->{features} });
				push @nodes, $new_node;
				
				# and start a new node
				
				$new_node = {};
				$new_node->{features} = [$af];
				$new_node->{tot_score} = $af->score;
				$new_node->{avg_percentid} = $af->percent_id;
			}
			else {
				print STDERR "ClusterDepthFilter: continuing node: last->hend: ",
					$last->hend," af->hstart: ", $af->hstart, "\n" if $DEBUG;
				
				# continue the same node
				
				push @{ $new_node->{features} }, $af;
				$new_node->{tot_score} += $af->score;
				$new_node->{avg_percentid} += $af->percent_id;
			}
			
			$last = $af;
		}
		
		# add the final node to the list
		
		$new_node->{avg_percentid} /= scalar( @{ $new_node->{features} } );
		push @nodes, $new_node;
	}

	print STDERR "ClusterDepthFilter: "
	  . scalar( @nodes )
	  . " nodes after processing discontinuities\n";

	# sort the nodes so they are arranged sequentially against the
	# genomic sequence, with the longest sequence first for nodes 
	# with identical starts (this allows us to find all overlaps in 
	# one pass)

	@nodes = sort { 
			($a->{features}->[0]->start <=> $b->{features}->[0]->start) ||
	  		($b->{features}->[-1]->end <=> $a->{features}->[-1]->end) 
		} @nodes;
	  	
	# try to cluster nodes according to overlaps and introns

	my @clusters = ();

	for my $node (@nodes) {

		my @afs = @{ $node->{features} };

		$node->{start} = $afs[0]->start;
		$node->{end}   = $afs[-1]->end;
		
		print STDERR "ClusterDepthFilter: looking at node: ",
			node_str($node),"\n" if $DEBUG;

		my @introns = ();

		if ( @afs > 1 ) {    # we have at least 1 intron

			# identify all the introns

			my $start = $afs[0]->end;
			for my $af ( @afs[ 1 .. $#afs - 1 ] ) {
				my $end = $af->start;
				push @introns, [ $start, $end ];
				$start = $af->end;
			}
			my $end = $afs[-1]->start;
			push @introns, [ $start, $end ];
		}

		$node->{introns} = [@introns];

		# check if we can add this node to any existing cluster

		my $added = 0;

		for my $cluster (@clusters) {
			if ( $node->{start} >= $cluster->{start} &&
			     $node->{end} <= $cluster->{end} &&
			     is_consistent($node, $cluster) ) {                                                                       

				# this node is subsumed by the cluster and its 
				# introns are consistent, so we can add it

				print STDERR "ClusterDepthFilter: adding node to cluster: ",
					cluster_str($cluster),"\n" if $DEBUG;

				push @{ $cluster->{nodes} }, $node;

				$added = 1;

				last;
			}
			elsif ( $node->{start} <= $cluster->{start} &&
				    $node->{end} >= $cluster->{end} &&
				    is_consistent($node, $cluster) ) {
				die "Error: this should never happen - investigate";
			}
			else {
				print STDERR "ClusterDepthFilter: can't add node to cluster: ",
					cluster_str($cluster),"\n" if $DEBUG;
			}
		}

		# if we didn't find a matching cluster, start a new one
		if ( !$added ) {
			
			my $cluster = {
				start   => $node->{start},
				end     => $node->{end},
				introns => [ map { [@$_] } @introns ],
				nodes   => [$node]
			};
			push @clusters, $cluster;
			
			print STDERR "ClusterDepthFilter: creating new cluster: ",
					cluster_str($cluster),"\n" if $DEBUG;
		}
	}

	print STDERR "ClusterDepthFilter: grouped "
	  . scalar(@nodes)
	  . " nodes into "
	  . scalar(@clusters)
	  . " clusters";
	  
	print STDERR sprintf(" (average of %.2f nodes per cluster)",
	  ( scalar(@nodes) / scalar(@clusters) )) if @clusters > 0;
	
	print STDERR "\n";

	map {print STDERR cluster_str($_),"\n"} @clusters if $DEBUG;

	# now actually identify the features we want to include

	my @filtered_features = ();

	my @summary_features = ();

	for my $cluster (@clusters) {

		# sort the nodes in this cluster by tot_score, breaking ties 
		# according to avg_percentid

		my @sorted =
		  sort {
			     ( $b->{tot_score} <=> $a->{tot_score} )
			  || ( $b->{avg_percentid} <=> $a->{avg_percentid} )
		  } @{ $cluster->{nodes} };

		# and add the features of the $nodes_to_keep best hits to
		# the list of filtered features

		my $limit = $#sorted < $nodes_to_keep ? $#sorted : $nodes_to_keep - 1;

		for my $node_to_keep ( @sorted[ 0 .. $limit ] ) {
			
			print STDERR "ClusterDepthFilter: adding features of ",$limit,
				"/",$#sorted," nodes of cluster: ",
				$cluster->{start},"-",$cluster->{end},"\n" if $DEBUG;
			
			for my $af ( @{ $node_to_keep->{features} } ) {
				
				print STDERR "ClusterDepthFilter: adding feature: ",
					$af->start,"-",$af->end,"\n" if $DEBUG;
				
				$af->analysis( $self->analysis );
				$af->dbID(0);
				$af->{adaptor} = undef;
				push @filtered_features, $af;
			}
		}

		if ( $limit < $#sorted ) {

			# create new 'summary' features representing the maximal
			# coverage of the remaining features (if we have any)

			# establish the start, end, introns list and score of
			# the summary feature

			my $first = $sorted[ $limit + 1 ];

			my $start   = $first->{features}->[0]->start;
			my $end     = $first->{features}->[-1]->end;
			my $score   = $first->{tot_score};
			my $introns = $first->{introns};

			for my $node ( @sorted[ ( $limit + 2 ) .. $#sorted ] ) {
				
				my $s = $node->{features}->[0]->start;
				$start = $s if $s < $start;

				my $e = $node->{features}->[-1]->end;
				$end = $e if $e > $end;

				my $ms = $node->{tot_score};
				$score = $ms if $ms > $score;

				my $new_introns =
				  find_new_introns( $introns, $node->{introns} );
				push @$introns, map { [@$_] } @$new_introns;
				$introns = [ sort { $a->[0] <=> $b->[0] } @$introns ];
			}

			# create our summary features with the maximal set of
			# introns

			for my $intron (@$introns) {

				my $summary = Bio::EnsEMBL::SimpleFeature->new(
					-start         => $start,
					-end           => $intron->[0],
					-strand        => 0,              # we mix both strands here
					-score         => $score,
					-display_label =>
					  sprintf( "summary feature from %d redundant nodes",
						$#sorted - $limit ),
					-analysis => $self->analysis()
				);
				push @summary_features, $summary;

				$start = $intron->[1];
			}

			my $summary = Bio::EnsEMBL::SimpleFeature->new(
				-start         => $start,
				-end           => $end,
				-strand        => 0,        # we mix both strands here
				-score         => $score,
				-display_label =>
				  sprintf( "summary feature from %d redundant nodes",
					$#sorted - $limit ),
				-analysis => $self->analysis()
			);
			push @summary_features, $summary;
		}
	}

	print STDERR "ClusterDepthFilter: "
	  . scalar(@summary_features)
	  . " summary features\n";

	print STDERR "ClusterDepthFilter: "
	  . scalar(@filtered_features)
	  . " features after filtering\n";

	return ( \@filtered_features, \@summary_features );
}

1;

=head1 NAME - Bio::EnsEMBL::Analysis::RunnableDB::Finished::ClusterDepthFilter

=head2 AUTHOR
James Gilbert B<email> jgrg@sanger.ac.uk    - original implementation

=head2 AUTHOR
Mustapha Larbaoui B<email> ml6@sanger.ac.uk - new pipeline

=head2 AUTHOR
Leo Gordon B<email> lg4@sanger.ac.uk        - porting

=head2 AUTHOR
Graham Ritchie B<email> gr5@sanger.ac.uk    - new clustering filter algorithm
