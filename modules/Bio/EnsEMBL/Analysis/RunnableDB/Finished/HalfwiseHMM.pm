# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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


=pod

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Finished::HalfwiseHMM

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Analysis::RunnableDB::Finished::HalfwiseHMM->new(
					     -db     => $db,
					     -input_id  => $id
                                             );
    $obj->fetch_input
    $obj->run

    my @newfeatures = $obj->output;


=head1 DESCRIPTION

runs HalfwiseHMM runnable and converts it output into genes which can be stored in an ensembl database

=head1 CONTACT

lec@sanger.ac.uk
refactored by Sindhu K. Pillai sp1@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Analysis::RunnableDB::Finished::HalfwiseHMM;

use warnings ;
use Bio::EnsEMBL::Analysis::RunnableDB::Finished;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::Runnable::Finished::HalfwiseHMM;
use Bio::EnsEMBL::Analysis::Tools::BlastDBTracking;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::Similarity qw (
  GB_SIMILARITY_DATABASES
);

use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::Finished);

=head2  new

    Arg      : all those inherited from RunnableDB
    Function : Make a new HalfwiseHMM object defining the above variables
    Exception: thrown if no input id is provided, this is found in RunnableDB
    Caller   :
    Example  : $runnable = Bio::EnsEMBL::Analysis::RunnableDB::HalfwiseHMM new->(-db => $db
										 -INPUT_ID => $id
										 -ANALYSIS => $analysis);

=cut

sub new {

	my ( $new, @args ) = @_;
	my $self = $new->SUPER::new(@args);

	# db, input_id, seqfetcher, and analysis objects are all set in
	# in superclass constructor (RunnableDB.pm)
	my ( $type, $threshold ) = rearrange( [qw(TYPE THRESHOLD)], @args );
	$self->{' '} = [];    #create key to an array of feature pairs
	return $self;

}

sub type {

	my ( $self, $type ) = @_;
	if ( defined($type) ) {
		$self->{_type} = $type;
	}
	return $self->{_type};

}

sub threshold {

	my ( $self, $threshold ) = @_;
	if ( defined($threshold) ) {
		$self->{_threshold} = $threshold;
	}
	return $self->{_threshold};

}

=head2  fetch_input

    Arg      : none
    Function : fetches the repeatmasked sequence and the uniprot features for the slice being run on and creates the HalfwiseHMM Runnable
    Exception: throws if no input_id has been provided
    Caller   :
    Example  :

=cut

sub fetch_input {

	my ($self) = @_;
	my %ests;
	my @estseqs;
	throw("No input id") unless defined( $self->input_id );
	my $sliceid = $self->input_id;
	my $sa      = $self->db->get_SliceAdaptor();
	my $slice   = $sa->fetch_by_name($sliceid);
	$slice->{'seq'} = $slice->seq();
	$self->query($slice);
	my $maskedslice =
	  $slice->get_repeatmasked_seq( $ANALYSIS_REPEAT_MASKING, $SOFT_MASKING )
	  or throw("Unable to fetch contig");
	my $alignAdaptor = $self->db->get_ProteinAlignFeatureAdaptor();
	my $stateinfocontainer = $self->db->get_StateInfoContainer;
	my $analysis_adaptor = $self->db->get_AnalysisAdaptor();

	foreach my $database ( @{$GB_SIMILARITY_DATABASES} ) {
		my $analysis_ln = $database->{'type'};
		my $threshold   = $database->{'threshold'};
		my $fps      = [];
		my $features =
		  $alignAdaptor->fetch_all_by_Slice_and_pid( $slice,
			$threshold,
			$analysis_ln );
		print STDERR
		  "Number of features matching threshold $database->{'threshold'} = "
		  . scalar(@$features) . "\n";
		foreach my $f (@$features) {
			if ( UNIVERSAL::isa( $f, "Bio::EnsEMBL::FeaturePair" )
				&& defined( $f->hseqname ) )
			{
				push( @$fps, $f );
			}
		}

		my $runnable =
		  Bio::EnsEMBL::Analysis::Runnable::Finished::HalfwiseHMM->new(
			'-query'    => $maskedslice,
			'-features' => $fps,
			'-pfamdb'   => $self->getPfamDB(),
			'-hmmdb'	=> $self->analysis->db_file,
			'-options'  => $self->analysis->parameters(),
			'-program'  => $self->analysis->program(),
			'-analysis' => $self->analysis
		  );
		$self->runnable($runnable);

		# set the db version searched which is a concatenation of Uniprot and Pfam db versions
		my $ana = $analysis_adaptor->fetch_by_logic_name($analysis_ln);
		my $uniprot_db_version = $stateinfocontainer->fetch_db_version($self->input_id,$ana);
		my $pfam_db_version    = $self->get_db_version();
		$self->db_version_searched(join('_',$uniprot_db_version,$pfam_db_version));
	}

}

sub make_hash_from_meta_value {

	my ( $self, $string ) = @_;
	if ($string) {
		my $hash = { eval $string };
		if($@) {
			die "error evaluating $string [$@]";
		} else {
			return $hash || {};
		}
	}

	return {};
}

sub getPfamDB {

	my ($self) = @_;
	unless ( $self->{'_pfam_db'} ) {
		my $pfam_meta = $self->db->get_MetaContainer();
		my $value     = $pfam_meta->list_value_by_key('pfam_db')
		  || throw("please enter pfam_db key - value into meta table\n");
		my $pfam_db_conn = $self->make_hash_from_meta_value( $value->[0] );

                $pfam_db_conn->{-reconnect_when_connection_lost} = 1;
                $pfam_db_conn->{-disconnect_when_inactive} = 1;

		# Use the Blast tracking system to set the correct Pfam DB name
		my $db_file_path    = $self->analysis->db_file;
		my $db_file_version = $self->get_db_version($db_file_path);
		if ( $db_file_version =~ /^(\d+).(\d+)$/ ) {
			$pfam_db_conn->{'-dbname'} = "pfam_$1_$2";
		}
		elsif ( $db_file_version =~ /^(\d+)$/ ) {
			$pfam_db_conn->{'-dbname'} = "pfam_$1_0";
		}
		$self->{'_pfam_db'} =
		  Bio::EnsEMBL::DBSQL::DBConnection->new(%$pfam_db_conn);
                # nb. there is not yet a db_handle (DBI connection)
	}
	return $self->{'_pfam_db'};
}

=head2 get_db_version

    Title   :  get_db_version
               [ distinguished from RunnableDB::*::db_version_searched() ]
    Useage  :  $self->get_db_version('/data/base/path')
               $obj->get_db_version()
    Function:  Set a blast database version from the supplied path
               Get a blast database version from previously supplied path
    Returns :  String
    Args    :  String (should be a full database path)
    Caller  :  $self::fetch_databases()
               RunnableDB::Finished_EST::db_version_searched()

=cut

sub get_db_version {
	my ( $self, $db ) = @_;
        my $ver = Bio::EnsEMBL::Analysis::Tools::BlastDBTracking::get_db_version_mixin(
            $self, '_pfam_db_version', $db,
            );
	return $ver;
}

sub db_version_searched {
	my ( $self, $arg ) = @_;
	$self->{'_db_version_searched'} = $arg if $arg;
	return $self->{'_db_version_searched'};
}

=head2  runnable

    Arg      : a Bio::EnsEMBL::Analysis::Runnable
    Function : Gets/sets the runnable
    Exception: throws if argument passed isn't a runnable
    Caller   :
    Example  :'

=cut

sub runnable {
	my ( $self, $arg ) = @_;
	if ( !defined( $self->{'_seqfetchers'} ) ) {
		$self->{'_seqfetchers'} = [];
	}
	if ( defined($arg) ) {
		throw("[$arg] is not a Bio::EnsEMBL::Analysis::Runnable")
		  unless $arg->isa("Bio::EnsEMBL::Analysis::Runnable");

		push( @{ $self->{_runnable} }, $arg );
	}

	return @{ $self->{_runnable} };
}

sub run {
	my ($self) = @_;
	foreach my $runnable ( $self->runnable ) {
		$runnable || throw("Can't run - no runnable object");
		print STDERR "using " . $runnable . "\n";
		eval { $runnable->run; };
		if ( my $err = $@ ) {
			chomp $err;
			$self->failing_job_status($1)
			  if $err =~ /^\"([A-Z_]{1,40})\"$/i; # only match '"ABC_DEFGH"'
			throw("$@");
		}
	}
	$self->_convert_output();
}

=head2  write_output


    Arg      : none
    Function : writes the converted output to the database as genes
    Exception: none
    Caller   :
    Example  :

=cut

sub write_output {

	my ($self) = @_;
	my @times = times;
	print STDERR "started writing @times \n";
	my @genes        = $self->output();
	my $db           = $self->db();
	my $gene_adaptor = $db->get_GeneAdaptor;
	my $sliceid      = $self->input_id;
    my $sa           = $self->db->get_SliceAdaptor();
    my $slice        = $sa->fetch_by_name($sliceid);
	my $dbh          = $db->dbc->db_handle;
	$dbh->begin_work;
	eval {
		# delete old genes first
		foreach my $gene (@{$slice->get_all_Genes($self->analysis->logic_name)})
        {
            $gene_adaptor->remove($gene);
        }
		# now save new genes 
		foreach my $gene (@genes)
		{
			$gene_adaptor->store($gene,0);
		}
		$dbh->commit;
	};
	if ($@) {
		$dbh->rollback;
		throw("UNABLE TO WRITE GENES IN DATABASE\n[$@]\n");
	}
	@times = times;
	print STDERR "finished writing @times \n";
	return 1;
}

=head2  _convert_output

    Arg      : none
    Function : takes the features from the halfwise runnable and runs _make_genes to convert them into Bio::EnsEMBL::Genes with appropriately attached exons and supporting evidence
    Exception: thows if there are no analysis types
    Caller   :
    Example  :

=cut

sub _convert_output {

	my ($self) = @_;
	my @genes;
	my $analysis = $self->analysis();
	my $genetype = 'Halfwise';

	# make an array of genes for each runnable
	my @out;
	foreach my $runnable ( $self->runnable ) {
		push( @out, $runnable->output );
		$self->pfam_lookup( $runnable->pfam_lookup )
		  if $runnable->can('pfam_lookup');
	}
	my @g = $self->_make_genes( $genetype, $analysis, \@out );
	push( @genes, @g );

	$self->output(@genes);
}

# get/set for lookup multi-valued hash { pfam_id => [pfam_acc, pfam_desc], ... }
# can append multiple to the lookup (see { %{$self->{'_pfam_lookup'}}, %{$hash_ref} })

sub pfam_lookup {
	my ( $self, $hash_ref ) = @_;
	if ( ref($hash_ref) eq 'HASH' ) {
		$self->{'_pfam_lookup'} ||= {};
		$self->{'_pfam_lookup'} =
		  { %{ $self->{'_pfam_lookup'} }, %{$hash_ref} };
	}
	return $self->{'_pfam_lookup'};
}

=head2  _make_genes

    Arg      : runnable being run and analysis object being used
    Function : converts the seqfeatures outputed by the runnable and actually converts them into Bio::EnsEMBL::Genes
    Exception: none
    Caller   :
    Example  :

=cut

=head2 _make_genes

  Title   :   make_genes
  Usage   :   $self->make_genes
  Function:   makes Bio::EnsEMBL::Genes out of the output from runnables
  Returns :   array of Bio::EnsEMBL::Gene
  Args    :   $genetype: string
              $analysis_obj: Bio::EnsEMBL::Analysis
              $results: reference to an array of FeaturePairs

=cut

sub _make_genes {

	my ( $self, $genetype, $analysis_obj, $results ) = @_;
	my $sliceid = $self->input_id;
	my $sa      = $self->db->get_SliceAdaptor();
	my $slice   = $sa->fetch_by_name($sliceid);
	my @genes;
	my $info_type	= 'MISC';
	my $info_text	= 'Relationship generated from genewisedb search';

	# fetch lookup multi-valued hash { pfam_id => [pfam_acc, pfam_desc], ... }

	my $pfam_lookup  = $self->pfam_lookup();
	my $pfam_release = $self->get_db_version();
	$self->_check_that_external_db_table_populated( $pfam_release, 'PFAM',
		'XREF' );

	foreach my $tmp_gene (@$results) {
		my $pfam_id     = $tmp_gene->seqname();
		my $acc_ver     = $pfam_lookup->{$pfam_id}->[0];
		my @pfamacc_ver = split /\./, $acc_ver;
		my $dbentry     = Bio::EnsEMBL::DBEntry->new(
			-primary_id  => $pfamacc_ver[0],
			-display_id  => $pfam_id,
			-version     => $pfamacc_ver[1],
			-release     => $pfam_release,
			-dbname      => "PFAM",
			-description => $pfam_lookup->{$pfam_id}->[1],
			-info_type	=>	$info_type,
			-info_text	=>	$info_text
		);
		$dbentry->status('XREF');

		my $gene       = Bio::EnsEMBL::Gene->new();
		my $transcript =
		  $self->_make_transcript( $tmp_gene, $slice, $genetype,
			$analysis_obj );
		$gene->biotype($genetype);
		$gene->analysis($analysis_obj);
		$gene->add_Transcript($transcript);
        
        # Add XRef to DBEntry list and as the display XREf to ensure
        # compatibility with both old and new (GFF based) fetching for otterlace
        $gene->add_DBEntry($dbentry);
        $gene->display_xref($dbentry);

		push( @genes, $gene );
	}
	return @genes;
}

sub _check_that_external_db_table_populated {
	my ( $self, $release, $name, $status ) = @_;
	$status ||= 'XREF';
	my $db             = $self->db();
	my $find_tuple_sql = qq(SELECT count(*) AS tuple_exists
			  FROM external_db
			  WHERE db_name = ?
			  && db_release = ?);
	my $sth = $db->prepare($find_tuple_sql);
	$sth->execute( $name, $release );
	my $tuple = $sth->fetchrow_hashref() || {};
	$sth->finish();

	# if there is one return do nothing and the job can get on as normal
	return if $tuple->{'tuple_exists'};

	# else lock table external_db write
	$sth = $db->prepare("LOCK TABLES external_db WRITE");
	$sth->execute();
	$sth->finish();

	# get the next external_db_id
	my $max_db_id_sql =
	  q`SELECT MAX(external_db_id) + 1 AS next_db_id from external_db`;
	$sth = $db->prepare($max_db_id_sql);
	$sth->execute();
	my $row       = $sth->fetchrow_hashref || {};
	my $max_db_id = $row->{'next_db_id'}   || warning "Error";

	# insert the row
	my $insert_sql =
q`INSERT INTO external_db (external_db_id, db_name, db_release, status, priority) VALUES(?, ?, ?, ?, 5)`;
	$sth = $db->prepare($insert_sql);
	$sth->execute( $max_db_id, $name, $release, $status );
	$sth->finish();

	# unlock tables;
	$sth = $db->prepare("UNLOCK TABLES");
	$sth->execute();
	return $max_db_id;
}

=head2 _make_transcript

 Title   : make_transcript
 Usage   : $self->make_transcript($gene, $slice, $genetype)
 Function: makes a Bio::EnsEMBL::Transcript from a SeqFeature representing a gene,
           with sub_SeqFeatures representing exons.
 Example :
 Returns : Bio::EnsEMBL::Transcript with Bio::EnsEMBL:Exons(with supporting feature
           data), and a Bio::EnsEMBL::translation
 Args    : $gene: Bio::EnsEMBL::SeqFeatureI, $contig: Bio::EnsEMBL::RawContig,
  $genetype: string, $analysis_obj: Bio::EnsEMBL::Analysis


=cut

sub _make_transcript {

	my ( $self, $gene, $slice, $genetype, $analysis_obj ) = @_;
	$genetype = 'unspecified' unless defined($genetype);

	unless ( $gene->isa("Bio::EnsEMBL::SeqFeatureI") ) {
		print "$gene must be Bio::EnsEMBL::SeqFeatureI\n";
	}

	my $transcript  = Bio::EnsEMBL::Transcript->new();
	my $translation = Bio::EnsEMBL::Translation->new();

	$transcript->translation($translation);
	$transcript->analysis($analysis_obj);

	my $excount = 1;
	my @exons;

	foreach my $exon_pred ( $gene->sub_SeqFeature ) {

		# make an exon
		my $exon = Bio::EnsEMBL::Exon->new();
		my $sa   = $slice->adaptor();
		$exon->display_id( $sa->get_seq_region_id($slice) );
		$exon->start( $exon_pred->start );
		$exon->end( $exon_pred->end );
		$exon->strand( $exon_pred->strand );
		$exon->phase( $exon_pred->phase || 0 );
		$exon->end_phase(0);
		$exon->slice($slice);

		# sort out supporting evidence for this exon prediction

		foreach my $subf ( $exon_pred->sub_SeqFeature ) {
			$subf->feature1->seqname( $slice->get_seq_region_id );
			$subf->feature1->score(100);
			$subf->feature1->analysis($analysis_obj);
			$subf->feature2->score(100);
			$subf->feature2->analysis($analysis_obj);
			$exon->add_Supporting_Feature($subf);
		}

		push( @exons, $exon );
		$excount++;
	}

	if ( @exons < 0 ) {

		# printSTDERR "Odd.  No exons found\n";
	}
	else {
		if ( $exons[0]->strand == -1 ) {
			@exons = sort { $b->start <=> $a->start } @exons;
		}
		else {
			@exons = sort { $a->start <=> $b->start } @exons;
		}

		foreach my $exon (@exons) {
			$transcript->add_Exon($exon);
		}

		$translation->start_Exon( $exons[0] );
		$translation->end_Exon( $exons[$#exons] );

		if ( $exons[0]->phase == 0 ) {
			$translation->start(1);
		}
		elsif ( $exons[0]->phase == 1 ) {
			$translation->start(3);
		}
		elsif ( $exons[0]->phase == 2 ) {
			$translation->start(2);
		}

		$translation->end( $exons[$#exons]->end - $exons[$#exons]->start + 1 );
	}

	return $transcript;
}

=head2 output

 Title   : output
 Usage   :
 Function: get/set for output array
 Example :
 Returns : array of Bio::EnsEMBL::Gene
 Args    :


=cut

sub output {
	my ( $self, @genes ) = @_;
	if ( !defined( $self->{'_output'} ) ) {
		$self->{'_output'} = [];
	}
	if ( scalar(@genes) ) {
		push( @{ $self->{'_output'} }, @genes );
	}
	return @{ $self->{'_output'} };
}

1;
