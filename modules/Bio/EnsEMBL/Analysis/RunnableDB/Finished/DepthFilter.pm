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
	my $analysis = $self->analysis->logic_name;
	$analysis =~ s/df_//;
	my $prot_feat_a = $self->db->get_ProteinAlignFeatureAdaptor;
	my $dna_feat_a  = $self->db->get_DnaAlignFeatureAdaptor;
	my $db_features =
	  $dna_feat_a->fetch_all_by_Slice( $self->query, $analysis );
	if (!(@$db_features))
	{
		$db_features =
		  $prot_feat_a->fetch_all_by_Slice( $self->query, $analysis );
	}
	$self->output( $self->depth_filter($db_features) );

}

sub depth_filter {
	my ( $self, $db_features ) = @_;
	my @filtered_features;
	if (@$db_features) {

		my $feat = $db_features->[0];
		$feat->analysis( $self->analysis );
		$feat->{'dbID'} = 0;
		push @filtered_features, $feat;

	}
	return \@filtered_features;
}



1;

=head1 NAME - Bio::EnsEMBL::Analysis::RunnableDB::Finished::DepthFilter

=head1 AUTHOR

Mustapha Larbaoui B<email> ml6@sanger.ac.uk
