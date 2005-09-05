### Bio::EnsEMBL::Analysis::RunnableDB::DepthFilter

package Bio::EnsEMBL::Analysis::RunnableDB::Finished::DepthFilter;

use strict;
use Bio::EnsEMBL::Analysis::Config::General;

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
	my $orig_analysis = $self->analysis->logic_name;
	$orig_analysis =~ s/df_//;

	my $dna_feat_a  = $self->db->get_DnaAlignFeatureAdaptor;
    my $prot_feat_a = $self->db->get_ProteinAlignFeatureAdaptor;

        # Try the DNA features first:
	my $orig_features =
	  $dna_feat_a->fetch_all_by_Slice( $self->query, $orig_analysis );

        # If we haven't found anything, try Protein features:
	if (!(@$orig_features))
	{
		$orig_features =
		  $prot_feat_a->fetch_all_by_Slice( $self->query, $orig_analysis );
	}
	$self->output( $self->depth_filter($orig_features) );

}

sub depth_filter {
	my ($self, $orig_features, $max_coverage) = (@_, 10);

    my %grouped_byname = ();

    for my $af (@$orig_features) {
        my $node = $grouped_byname{$af->hseqname()} ||= {};

        my ($score, $percentid) = ($af->score(), $af->percent_id());
        if(%$node) { # nonempty
            $node->{max_score} = $score if $score>$node->{max_score};

            $node->{max_percentid} = $percentid if $percentid>$node->{max_percentid};
        } else {
            $node->{max_score} = $score;
            $node->{max_percentid} = $percentid;
        }
        push @{$node->{features}}, $af;
    }

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

    return \@filtered_features;
}

1;

=head1 NAME - Bio::EnsEMBL::Analysis::RunnableDB::Finished::DepthFilter

=head2 AUTHOR
James Gilbert B<email> jgrg@sanger.ac.uk    - original implementation

=head2 AUTHOR
Mustapha Larbaoui B<email> ml6@sanger.ac.uk - new pipeline

=head2 AUTHOR
Leo Gordon B<email> lg4@sanger.ac.uk        - porting

