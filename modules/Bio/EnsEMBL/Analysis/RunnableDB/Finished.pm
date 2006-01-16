### Bio::EnsEMBL::Analysis::RunnableDB::Finished

package Bio::EnsEMBL::Analysis::RunnableDB::Finished;

use strict;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Pfetch;
use base 'Bio::EnsEMBL::Analysis::RunnableDB';

sub write_output {
	my ($self) = @_;
	my $outputs = $self->Bio::EnsEMBL::Analysis::RunnableDB::Finished::output;
    ### code db_version_searched method may be duplicated in several modules
    ### How does it get written into the input_id_analysis table?
    foreach my $runnable ( @{ $self->runnable } ) {
	    my $db_version = $runnable->get_db_version
	      if $runnable->can('get_db_version');
	    $self->db_version_searched($db_version);    # make sure we set this here
    }
    # We expect an array of arrays from output()
    foreach my $output (@$outputs) {
        next unless @$output;   # No feature output        
        my $feat = $output->[0];
        
        # The type of adaptor used to store the data depends
        # upon the type of the first element of @$output
        my $adaptor = $self->get_adaptor($feat);
        
	    # Remove AlignFeatures already in db from each output
        if (ref($feat) =~ /AlignFeature$/) {
            $self->remove_stored_AlignFeatures($adaptor, $output);       
            $self->write_descriptions($output);
        }
	    my $analysis = $self->analysis;
    	my $slice    = $self->query;
	    my $ff       = $self->feature_factory;
        # Store features in the database
        foreach my $feature (@$output) {
        	$feature->analysis($analysis);
        	$feature->slice($slice) if (!$feature->slice);
        	$ff->validate($feature);
            eval { $adaptor->store($feature); };
            if ($@) {
                throw("RunnableDB:store failed, failed to write " . $feature
                      . " to the database: $@");
            }
        }
    }
}

sub output{
  my ($self, @output) = @_;
  if(!$self->{'output'}){
    $self->{'output'} = [];
  }
  if(scalar(@output)){
    push(@{$self->{'output'}}, @output);
  }
  if(ref($self->{'output'}->[0]) ne 'ARRAY')
  {
  	return [ $self->{'output'} ] ;
  }
  return $self->{'output'};
}

sub replace_output {
  my ($self, @output) = @_;
  if(scalar(@output)){
  	$self->{'output'} = [];
    push(@{$self->{'output'}}, @output);
  }
  return $self->{'output'};
}

sub remove_stored_AlignFeatures {
    my ($self, $adaptor, $output) = @_;
    ## create a hashtable of key='feature signature' and value=1
    my $db_features =
      $adaptor->fetch_all_by_Slice($self->query, $self->analysis->logic_name);
    print STDERR "Found ", scalar(@$db_features), " features already in db\n";
    my %db_feat = map { $self->get_feature_key($_), 1 } @$db_features;
    ## remove duplicated features
    for (my $i = 0 ; $i < @$output ;) {
    	my $feature = $output->[$i];
        $feature->slice($self->query) if (!$feature->slice);
        my $f_key = $self->get_feature_key($feature);
        if ($db_feat{$f_key}) {
            print "Feature ", $f_key, " already in table\n";
            splice(@$output, $i, 1);
        }
        else {
            $i++;
        }
    }
}

sub write_descriptions {
    my( $self, $output ) = @_;
    my $dbobj = $self->db;
    my $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Pfetch->new;
    my $ids = [map $_->hseqname, @$output];
    $seqfetcher->write_descriptions( $dbobj, $ids );

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
	my ($self, $feat) = @_;

	if ( $feat->isa('Bio::EnsEMBL::DnaPepAlignFeature') ) {
		return $self->db->get_ProteinAlignFeatureAdaptor;
	}
	elsif ( $feat->isa('Bio::EnsEMBL::DnaDnaAlignFeature') ) {
		return $self->db->get_DnaAlignFeatureAdaptor;
	}
    elsif ( $feat->isa('Bio::EnsEMBL::SimpleFeature') ) {
        return $self->db->get_SimpleFeatureAdaptor;
    }
	else {
		throw('unexpected feature type: '. ref($feat));
	}
}

1;

=head1 NAME - Bio::EnsEMBL::Analysis::RunnableDB::Finished

=head1 AUTHOR

Mustapha Larbaoui B<email> ml6@sanger.ac.uk
