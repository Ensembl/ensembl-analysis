### Bio::EnsEMBL::Analysis::RunnableDB::DepthFilter

package Bio::EnsEMBL::Analysis::RunnableDB::Finished::DepthFilter;

use strict;
use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::SimpleFeature;

use base 'Bio::EnsEMBL::Analysis::RunnableDB::Finished';

sub fetch_input {
	my ($self) = @_;
	my $slice =
	  $self->fetch_sequence( $self->input_id, $self->db,
		$ANALYSIS_REPEAT_MASKING, $SOFT_MASKING );
	$self->query($slice);

}

sub run {
	my ($self) = @_;

    my %params = %{ $self->parameters_hash() };
    my $max_coverage     = $params{max_coverage} || 10;
    my $percentid_cutoff = $params{percentid_cutoff} || 0.0;
	my $orig_analysis_name = $params{ori_analysis};
	my $hit_db = $params{hit_db};
	my $mode = $params{mode};

	# Get the blast db version from the raw analysis and save it
	my $analysis_adaptor = $self->db->get_AnalysisAdaptor();
	my $ana = $analysis_adaptor->fetch_by_logic_name($orig_analysis_name);
	my $stateinfocontainer = $self->db->get_StateInfoContainer;
	my $db_version = $stateinfocontainer->fetch_db_version($self->input_id,$ana);
	$self->db_version_searched($db_version);

    my $slice = $self->query();

	my $orig_features = $self->get_original_features($slice,$orig_analysis_name,$hit_db,$mode);

	my ( $filtered_features, $saturated_zones) =
        $self->depth_filter($orig_features, $slice, $max_coverage, $percentid_cutoff);

	$self->output($filtered_features, $saturated_zones);
}

sub get_original_features {
	my ($self,$slice,$analysis,$hit_db,$mode) = @_;
	my $hit_db_features = [];
    my $prot_feat_a = $self->db->get_ProteinAlignFeatureAdaptor;
    my $dna_feat_a  = $self->db->get_DnaAlignFeatureAdaptor;
    my $hit_desc_a = $self->db->get_HitDescriptionAdaptor;

    my $orig_features = $dna_feat_a->fetch_all_by_Slice( $slice, $analysis );
	$orig_features = $prot_feat_a->fetch_all_by_Slice( $slice, $analysis ) if (!(@$orig_features));

	# The block of code below discard old sequence versions from the original set
	# of features (should add "mode => single" in analysis parameters)
	if($mode && $mode eq "single") {
		print STDERR "DepthFilter: $mode mode\n";
		my $single_hash;
		my $old_seq;
		foreach my $feature (@$orig_features) {
			my $hit_name	= $feature->hseqname;
			my ($acc,$ver)	= $hit_name =~ /(\w+-?\w+)\.(\d+)/;
			my $key			= $acc;
			if($single_hash->{$key}) {
				my $stored_feature = $single_hash->{$key}->[0];
				my ($sacc,$sver)	= $stored_feature->hseqname =~ /(\w+-?\w+)\.(\d+)/;
				if($ver > $sver) {
					$single_hash->{$key} = [$feature];
					#print STDERR "DepthFilter: drop ".$stored_feature->hseqname." (".$feature->hseqname.")\n";
					push @$old_seq, $stored_feature;
				} elsif($ver == $sver) {
					push @{$single_hash->{$key}}, $feature;
				} else {
					#print STDERR "DepthFilter: drop ".$feature->hseqname." (".$stored_feature->hseqname.")\n";
					push @$old_seq, $feature;
				}
			} else {
				$single_hash->{$key} = [$feature];
			}
		}
		print STDERR "DepthFilter: drop ".scalar(@$old_seq)." old sequences\n" if @$old_seq;
		$orig_features = [];
		map ( push(@$orig_features,@$_) , values %$single_hash);
	}

    if($hit_db){
    	print STDERR "DepthFilter: hit db is $hit_db\n";
    	my $hit_hash = {map {$_->hseqname, undef} @$orig_features};
	    $hit_desc_a->fetch_HitDescriptions_into_hash($hit_hash);
	    foreach my $feat (@$orig_features) {
	        if (my $desc = $hit_hash->{$feat->hseqname}) {
	            push @$hit_db_features, $feat if($desc->db_name eq $hit_db);
	        }
	    }
	    print STDERR "DepthFilter: use ".scalar(@$hit_db_features)." features out of ".scalar(@$orig_features)."\n";
    } else {
    	return $orig_features;
    }

    return $hit_db_features
}


sub depth_filter {
	my ($self, $orig_features, $slice, $max_coverage, $percentid_cutoff) = @_;

	print STDERR "DepthFilter: MaxCoverage=$max_coverage\n";
	print STDERR "DepthFilter: PercentIdCutoff=$percentid_cutoff\n";
	print STDERR "DepthFilter: ".scalar(@$orig_features)." features before filtering\n";

    my %grouped_byname = ();

    for my $af (@$orig_features) {
        my ($score, $percentid) = ($af->score(), $af->percent_id());
        if($percentid < $percentid_cutoff) {
            next;
        }

        my $node = $grouped_byname{$af->hseqname()} ||= {};
        if(%$node) { # nonempty
            $node->{max_score} = $score if $score>$node->{max_score};

            $node->{max_percentid} = $percentid if $percentid>$node->{max_percentid};
        } else {
            $node->{max_score} = $score;
            $node->{max_percentid} = $percentid;
        }
        push @{$node->{features}}, $af;
    }

	print STDERR "DepthFilter: ".scalar(keys %grouped_byname)." unique hitnames\n";

    my @bisorted =
        sort { ($b->{max_score} <=> $a->{max_score})
            || ($b->{max_percentid} <=> $a->{max_percentid}) }
        values %grouped_byname;

    my @coverage_map = ();
    my @filtered_features = ();

    for my $node (@bisorted) {
        my $keep_node = 0;
        for my $af (sort {$a->start() <=> $b->start()} @{$node->{features}}) {
            for my $position ($af->start()..$af->end()) {
                my $depth = $coverage_map[$position] ||= 0;
                if($depth < $max_coverage) {
                    $keep_node = 1;
                }
                $coverage_map[$position]++;
            }
        }
        if($keep_node) {
            for my $af (@{$node->{features}}) {
                $af->analysis( $self->analysis );
                $af->dbID(0);
                $af->{adaptor} = undef;
                push @filtered_features, $af;
            }
        }
    }

	print STDERR "DepthFilter: ".scalar(@filtered_features)." features after filtering\n";

    my @saturated_zones = ();
    my $zone_start = undef;
    my $zone_score = 0;
    my $slice_length = $slice->length();
    for(my $i=1; $i<=$slice_length; $i++) {
        my $n = $coverage_map[$i] || 0;
        if ($zone_start) {
            $zone_score += $n;
            if ($n < $max_coverage) {
                my $new_zone = Bio::EnsEMBL::SimpleFeature->new(
                    -start  => $zone_start,
                    -end    => $i - 1,
                    -strand => 0,       # we mix both strands here
                    -score  => $zone_score,
                    -display_label => sprintf("avg depth = %.2f", $zone_score/($i-$zone_start)),
                    -analysis => $self->analysis(),
                );
                push(@saturated_zones, $new_zone);
                $zone_start = undef;
                $zone_score = 0;
            }
        }
        elsif ($n >= $max_coverage) {
            $zone_start = $i;
            $zone_score = $n;
        }
    }
    if ($zone_start) { # Are saturated up to end of contig:
        my $new_zone = Bio::EnsEMBL::SimpleFeature->new(
            -start  => $zone_start,
            -end    => $slice_length,
            -strand => 0,       # we mix both strands here
            -score  => $zone_score,
            -display_label => sprintf("avg depth = %.2f", $zone_score/($slice_length-$zone_start+1)),
            -analysis => $self->analysis(),
        );
        push(@saturated_zones, $new_zone);
    }

	print STDERR "DepthFilter: ".scalar(@saturated_zones)." saturated zones found\n";

    return (\@filtered_features, \@saturated_zones);
}

=head2 db_version_searched

    Title   :  db_version_searched
               [ distinguished from Runnable::*::get_db_version() ]
    Useage  :  $self->db_version_searched('version string')
               $obj->db_version_searched()
    Function:  Get/Set a blast database version that was searched
               The actual look up is done in Runnable::Finished::Blast
               This is just a holding place for the string in this
               module
    Returns :  String or undef
    Args    :  String
    Caller  :  $self::run()
               Job::run_module()

=cut

sub db_version_searched {

	my ( $self, $arg ) = @_;
	$self->{'_db_version_searched'} = $arg if $arg;
	return $self->{'_db_version_searched'};

}

1;

=head1 NAME - Bio::EnsEMBL::Analysis::RunnableDB::Finished::DepthFilter

=head2 AUTHOR
James Gilbert B<email> jgrg@sanger.ac.uk    - original implementation

=head2 AUTHOR
Mustapha Larbaoui B<email> ml6@sanger.ac.uk - new pipeline

=head2 AUTHOR
Leo Gordon B<email> lg4@sanger.ac.uk        - porting

