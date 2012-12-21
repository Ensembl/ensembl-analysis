### Bio::EnsEMBL::Analysis::RunnableDB::Finished::Exonerate

# $Source: /tmp/ENSCOPY-ENSEMBL-ANALYSIS/modules/Bio/EnsEMBL/Analysis/RunnableDB/Finished/Exonerate.pm,v $
# $Version: $
package Bio::EnsEMBL::Analysis::RunnableDB::Finished::Exonerate;

use warnings ;
use strict;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis::Runnable::Finished::Exonerate;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Pfetch;
use Bio::EnsEMBL::Analysis::Config::General;

use base 'Bio::EnsEMBL::Analysis::RunnableDB::Finished';

############################################################

sub fetch_input {

	my ($self) = @_;
	my $slice =
	  $self->fetch_sequence( $self->input_id, $self->db,
		$ANALYSIS_REPEAT_MASKING, 1 );
	$self->query($slice);

	my $maskedslice =
	     $slice->get_repeatmasked_seq( $ANALYSIS_REPEAT_MASKING, $SOFT_MASKING )
	  or throw("Unable to fetch contig");
	my $maskedseq = $maskedslice->seq();

	if ( scalar( $maskedseq =~ s/([CATG])/$1/g ) > 3 ) {
		$self->input_is_void(0);
	}
	else {
		$self->input_is_void(1);
		warning("Need at least 3 nucleotides");
	}

	# Incremental updating of the embl blast db analysis
	# The embl blast dbs are made up of release files embl_*
	# and update files emnew_*. This block of code makes
	# sure that the analysis is only run against new version of either
	# of these files.

	my @files = split(",", $self->analysis->db_file);
	my @patches;
	if($files[-1] =~ /^embl_/){
		my $search_only_patch = 0;
		my $sic = $self->db->get_StateInfoContainer;
		my $db_version_saved = $sic->fetch_db_version($self->input_id, $self->analysis);
		my $db_version_current = $self->analysis->db_version;
		if($db_version_saved) {
			# split the embl blast db version "12-Mar-06 (85)" to
			# patch version "12-Mar-06" and release version "85"
			my ($patch_sv,$release_sv) = $db_version_saved =~ /^(\S+)\s+\((\d+)\)$/;
			my ($patch_cv,$release_cv) = $db_version_current =~ /^(\S+)\s+\((\d+)\)$/;
			if($release_sv eq $release_cv){
				$search_only_patch = 1;
				print STDOUT "blast db files [ @files ] version $release_sv already searched\n";
				# Just to make sure that nothing is going wrong with the incremental updating...
				throw("Problem with the embl blast db incremental updating, saved and current version identical !\n
				   saved [$db_version_saved] = current [$db_version_current]\n") unless($patch_sv ne $patch_cv)
			}
		}
	    foreach my $file (@files) {
	    	my $patch_file = $file;
	    	$patch_file =~ s/^embl_/emnew_/g;
	    	$search_only_patch ? $file = $patch_file : push @patches,$patch_file;
	    }
	}
	$self->analysis->db_file(join(",",@files,@patches));


	my %parameters = %{ $self->parameters_hash };

	my $runnable = Bio::EnsEMBL::Analysis::Runnable::Finished::Exonerate->new(
		-analysis => $self->analysis,
		-program  => $self->analysis->program,
		-query_db => $self->analysis->db_file,
		-target   => $slice,
		%parameters

	);

	$self->runnable($runnable);
	return 1;
}


sub db_version_searched {

	my ( $self, $arg ) = @_;
	$self->{'_db_version_searched'} = $arg if $arg;
	return $self->{'_db_version_searched'};

}

1;

=head1 NAME - Bio::EnsEMBL::Analysis::RunnableDB::Finished::Exonerate

=head1 AUTHOR

Mustapha Larbaoui B<email> ml6@sanger.ac.uk
