### Bio::EnsEMBL::Analysis::RunnableDB::DepthFilter

package Bio::EnsEMBL::Analysis::RunnableDB::Finished::DepthFilter;

use strict;
use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::SimpleFeature;

use base 'Bio::EnsEMBL::Analysis::RunnableDB::Finished';

sub new {
	my ( $new, @args ) = @_;
	my $self = $new->SUPER::new(@args);

	return $self;
}

sub fetch_input {
	my ($self) = @_;
	my $slice =
	  $self->fetch_sequence( $self->input_id, $self->db,
		$ANALYSIS_REPEAT_MASKING, $SOFT_MASKING );
	$self->query($slice);

}

sub run {
	my ($self) = @_;
	my $orig_analysis_name = $self->analysis->logic_name;
	$orig_analysis_name =~ s/df_//;

    my %params = %{ $self->parameters_hash() };

    my $max_coverage     = $params{max_coverage} || 10;
    my $percentid_cutoff = $params{percentid_cutoff} || 0.0;

	my $dna_feat_a  = $self->db->get_DnaAlignFeatureAdaptor;
    my $prot_feat_a = $self->db->get_ProteinAlignFeatureAdaptor;

    my $slice = $self->query();

        # Try the DNA features first:
	my $orig_features =
	  $dna_feat_a->fetch_all_by_Slice( $slice, $orig_analysis_name );

        # If we haven't found anything, try Protein features:
	if (!(@$orig_features))
	{
		$orig_features =
		  $prot_feat_a->fetch_all_by_Slice( $slice, $orig_analysis_name );
	}

	my ( $filtered_features, $saturated_zones) =
        $self->depth_filter($orig_features, $slice, $max_coverage, $percentid_cutoff); 

	$self->output($filtered_features);

    # TODO:  and store the @$saturated_zones as simple features
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
        for my $af (sort {$a->start() <=> $b->start()} @{$node->{features}}) {
            my $keep_hit = 0;
            for my $position ($af->start()..$af->end()) {
                my $depth = $coverage_map[$position] ||= 0;
                if($depth < $max_coverage) {
                    $keep_hit = 1;
                }
                $coverage_map[$position]++;
            }
            if($keep_hit) {
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
                    -strand => $slice->strand(),
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
            -strand => $slice->strand(),
            -score  => $zone_score,
            -display_label => sprintf("avg depth = %.2f", $zone_score/($slice_length-$zone_start+1)),
            -analysis => $self->analysis(),
        );
        push(@saturated_zones, $new_zone);
    }

	print STDERR "DepthFilter: ".scalar(@saturated_zones)." saturated zones found\n";

    return (\@filtered_features, \@saturated_zones);
}

1;

=head1 NAME - Bio::EnsEMBL::Analysis::RunnableDB::Finished::DepthFilter

=head2 AUTHOR
James Gilbert B<email> jgrg@sanger.ac.uk    - original implementation

=head2 AUTHOR
Mustapha Larbaoui B<email> ml6@sanger.ac.uk - new pipeline

=head2 AUTHOR
Leo Gordon B<email> lg4@sanger.ac.uk        - porting

