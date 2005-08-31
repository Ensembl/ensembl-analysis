### Bio::EnsEMBL::Analysis::RunnableDB::Finished

package Bio::EnsEMBL::Analysis::RunnableDB::Finished;

use strict;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Pfetch;
use base 'Bio::EnsEMBL::Analysis::RunnableDB';

sub new {
	my ( $new, @args ) = @_;
	my $self = $new->SUPER::new(@args);

	return $self;
}

sub write_output {
	my ($self) = @_;
	foreach my $runnable ( @{ $self->runnable } ) {
		my $db_version = $runnable->get_db_version
		  if $runnable->can('get_db_version');
		$self->db_version_searched($db_version);    # make sure we set this here

		if ( my @output = @{ $runnable->output } ) {
			my $dbobj      = $self->db;
			my $seqfetcher =
			  Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Pfetch->new;
			my %ids = map { $_->hseqname, 1 } @output;
			my @ids_keys = keys(%ids);
			$seqfetcher->write_descriptions( $dbobj, \@ids_keys );
		}

	}
	### Remove features already in db from each runnable
	if(my $adaptor = $self->get_adaptor){ 
	# if no adaptor is provided it means that there is no output feature
		
		## create a hashtable of key='feature signature' and value=1
		my $db_features =
		  $adaptor->fetch_all_by_Slice( $self->query,
			$self->analysis->logic_name );
		print STDERR "Found ", scalar(@$db_features), " features already in db\n";
		my %db_feat = map { $self->get_feature_key($_), 1 } @$db_features;
		## remove duplicated features
		my $tmp_array = [];
		my $ff        = $self->feature_factory;
		my $output = $self->output;
		for ( my $i = 0 ; $i < scalar(@$output) ; $i++ ) {
			my $f = $output->[$i];
			$f->analysis( $self->analysis );
			$f->slice( $self->query ) if ( !$f->slice );
			$ff->validate($f);
			my $f_key = $self->get_feature_key($f);
			if ( !$db_feat{ $f_key } ) {
				push( @$tmp_array, $f );
			}
			else {
				print "Feature ", $f_key , " already in table\n";
			}
		}
		$self->{'output'} = $tmp_array;
	}

	$self->SUPER::write_output();
	return;
}

sub get_feature_key {
	my ( $self, $feat ) = @_;
	throw(
"Must pass Bio::EnsEMBL::Analysis::RunnableDB::Finished::get_feature_key a Bio::EnsEMBL::BaseAlignFeature"
		  . "not a "
		  . $feat )
	  unless ( $feat->isa('Bio::EnsEMBL::BaseAlignFeature') );
	return join( ':',
		$feat->display_id, $feat->seq_region_name, $feat->cigar_string,
		$feat->start, $feat->end, $feat->hstart, $feat->hend );
}

sub get_adaptor {
	my ($self) = @_;
	my $feat;
	if(!($feat = $self->output->[0])){ return 0; }
	if ( $feat->isa('Bio::EnsEMBL::DnaPepAlignFeature') ) {
		return $self->db->get_ProteinAlignFeatureAdaptor;
	}
	elsif ( $feat->isa('Bio::EnsEMBL::DnaDnaAlignFeature') ) {
		return $self->db->get_DnaAlignFeatureAdaptor;
	}
	else {
		throw(
'feature must be a Bio::EnsEMBL::DnaPepAlignFeature or a Bio::EnsEMBL::DnaDnaAlignFeature.'
		);
	}
}

1;

=head1 NAME - Bio::EnsEMBL::Analysis::RunnableDB::Finished

=head1 AUTHOR

Mustapha Larbaoui B<email> ml6@sanger.ac.uk