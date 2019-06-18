# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

# Bio::EnsEMBL::Analysis::Runnable::Finished::Est2Genome

=pod

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Finished::Est2Genome

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Analysis::Runnable::Finished::Est2Genome->new(
                                             -genomic => $genseq,
                                             -est     => $estseq
                                             );
    or
    my $obj = Bio::EnsEMBL::Analysis::Runnable::Finished::Est2Genome->new()

=head1 DESCRIPTION

Object to store the details of an est2genome run.
Stores the est2genome matches as an array of Bio::EnsEMBL::FeaturePair

=head2 Methods:

 new,
 genomic_sequence,
 est_sequence,
 run,
 output.

=head1 CONTACT

Modified by Sindhu K. Pillai B<email> sp1@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Analysis::Runnable::Finished::Est2Genome;

use vars qw($verbose);
use strict;
use warnings;

use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::SeqFeature;
use Bio::PrimarySeq;
use Bio::SeqIO;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

#compile time check for executable
use Bio::EnsEMBL::Analysis::Programs qw(est2genome);
use base 'Bio::EnsEMBL::Analysis::Runnable';

$verbose = 0;

=head2 new

    Title   :   new
    Usage   :   $self->new(-GENOMIC       => $genomicseq,
			   -EST           => $estseq,
                           -E2G           => $executable,
                           -ARGS          => $args);
    Function:   creates a
                Bio::EnsEMBL::Analysis::Runnable::Est2Genome object
    Returns :   Bio::EnsEMBL::Analysis::Runnable::Est2Genome object
    Args    :   -genomic:    Bio::PrimarySeqI object (genomic sequence)
                -est:        Bio::PrimarySeqI object (est sequence),
                -e2g:        Path to Est2Genome executable
                -args:       Arguments when running est2genome
=cut

sub new {

	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);
	$self->{'_fplist'}      = [];      # create key to an array of feature pairs
	$self->{'_clone'}       = undef;   # location of Bio::Seq object
	$self->{'_est_genome'}  = undef;   # location of est2genome
	$self->{'_workdir'}     = undef;   # location of temp directory
	$self->{'_filename'}    = undef;   # file to store Bio::Seq object
	$self->{'_estfilename'} = undef;   # file to store EST Bio::Seq object
	$self->{'_results'}     = undef;   # file to store results of analysis
	$self->{'_protected'}   = [];      # a list of files protected from deletion
	$self->{'_arguments'}   = undef;   # arguments for est2genome

	my ( $genomic, $est, $est_genome, $arguments ) =
	  rearrange( [qw(GENOMIC EST E2G ARGS)], @args );

	$self->genomic_sequence($genomic) if $genomic;
	$self->est_sequence($est)         if $est;

	if ($est_genome) {
		$self->est_genome($est_genome);
	}
	else {
		eval { $self->est_genome( $self->locate_executable('est2genome') ); };
		throw("Can't find executable") if $@;
	}
	if ($arguments) {
		$self->arguments($arguments);
	}
	else {
		$self->arguments(' -reverse ');
	}

	return $self;
}

#################
# get/set methods
#################

=head2 genomic_sequence

    Title   :   genomic_sequence
    Usage   :   $self->genomic_sequence($seq)
    Function:   Get/set method for genomic sequence
    Returns :   Bio::Seq object
    Args    :   Bio::Seq object

=cut

sub genomic_sequence {
	my ( $self, $value ) = @_;
	if ($value) {

		#need to check if passed sequence is Bio::Seq object
		$value->isa("Bio::PrimarySeqI")
		  || throw("Input isn't a Bio::PrimarySeqI");
		$self->{'_genomic_sequence'} = $value;
		$self->filename( $value->id . ".$$.seq" );
		$self->results( $self->filename . ".est_genome.out" );
	}
	return $self->{'_genomic_sequence'};
}

=head2 est_sequence

    Title   :   est_sequence
    Usage   :   $self->est_sequence($seq)
    Function:   Get/set method for est sequence
    Returns :   Bio::Seq object
    Args    :   Bio::Seq object

=cut

sub est_sequence {
	my ( $self, $value ) = @_;

	if ($value) {

		#need to check if passed sequence is Bio::Seq object
		$value->isa("Bio::PrimarySeqI")
		  || throw("Input isn't a Bio::PrimarySeqI");
		$self->{'_est_sequence'} = $value;
		$self->estfilename( $value->id . ".$$.est.seq" );
	}
	return $self->{'_est_sequence'};
}

sub estfilename {
	my ( $self, $estfilename ) = @_;
	$self->{'_estfilename'} = $estfilename if ($estfilename);
	return $self->{'_estfilename'};
}

sub filename {
	my ( $self, $filename ) = @_;
	$self->{'_filename'} = $filename if ($filename);
	return $self->{'_filename'};
}

sub results {
	my ( $self, $filename ) = @_;
	$self->{'_results'} = $filename if ($filename);
	return $self->{'_results'};
}

=head2 arguments

    Title   :   arguments
    Usage   :   $obj->arguments(' -reverse ');
    Function:   Get/set method for est_genome arguments
    Args    :   File path (optional)

=cut

sub arguments {
	my ( $self, $args ) = @_;
	if ($args) {
		$self->{'_arguments'} = $args;
	}
	return $self->{'_arguments'};
}

###########
# Analysis methods
##########

=head2 run

  Title   : run
  Usage   : $self->run()
            or
            $self->run("genomic.seq", "est.seq")
  Function: Runs est2genome and stores results as FeaturePairs
  Returns : TRUE on success, FALSE on failure.
  Args    : Temporary filenames for genomic and est sequences

=cut

sub run {

	my ( $self, @args ) = @_;

	# some constant strings
	$self->{'_source_tag'} = "est2genome";
	my $dirname = "/tmp/";

	#check inputs
	my $genomicseq = $self->genomic_sequence
	  || throw("Genomic sequence not provided");
	my $estseq = $self->est_sequence || throw("EST sequence not provided");

	#extract filenames from args and check/create files and directory
	my ( $genname, $estname ) = rearrange( [ 'genomic', 'est' ], @args );
	my ( $genfile, $estfile ) =
	  $self->_createfiles( $genname, $estname, $dirname );

	#use appropriate Bio::Seq method to write fasta format files
	{
		my $genOutput =
		  Bio::SeqIO->new( -file => ">$genfile", '-format' => 'Fasta' )
		  or throw("Can't create new Bio::SeqIO from $genfile : $!");
		my $estOutput =
		  Bio::SeqIO->new( -file => ">$estfile", '-format' => 'Fasta' )
		  or throw("Can't create new Bio::SeqIO from $estfile : $!");

		#fill inputs
		$genOutput->write_seq( $genomicseq );
		$estOutput->write_seq( $estseq );
	}

	my $est_genome_command = $BIN_DIR
	  . "/est2genome -reverse -genome $genfile -est $estfile -space 500000 -out stdout 2>&1 |";

	# use -align to get alignment
	print STDERR "\nRunning command $est_genome_command\n" if $verbose ;

	unless ( open( ESTGENOME, $est_genome_command ) ) {
		$self->_deletefiles( $genfile, $estfile );
		throw("Can't open pipe from '$est_genome_command' : $!");
	}
	my $firstline = <ESTGENOME>;
	print STDERR "firstline: \t$firstline" if $verbose;
	if ( $firstline =~ m/Align EST and genomic DNA sequences/ ) {

# Catch 'Align EST and genomic DNA sequences'. This comes from STDERR!! [ 2>&1 ]
		$firstline = <ESTGENOME>;
		print STDERR "\$firstline (secondline!): \t$firstline" if $verbose;
		if ( $firstline =~ /insufficient memory available/ ) {
			close(ESTGENOME) or warning("problem closing est_genome: $!\n");
			$self->_deletefiles( $genfile, $estfile );
			die qq{"OUT_OF_MEMORY"\n};
		}
	}

	if ( $firstline =~ m/reversed\sest/ ) {
		$self->est_strand(-1);
	}
	else {
		$self->est_strand(1);
	}

	if ( $firstline =~ m/forward\sgenome/ ) {
		$self->gen_strand(1);
	}
	else {
		$self->gen_strand(-1);
	}

	while (<ESTGENOME>) {
		print STDERR $_ if $verbose;
		if ( $_ =~ /^Segment|^Exon/ ) {

			# We only care about Segments in Exons
			# "gen" = genomic sequence
			my (
				$primary, $score,     $percent_id, $gen_start, $gen_end,
				$gen_id,  $est_start, $est_end,    $est_id
			  )
			  = (split)[ 0, 1, 2, 3, 4, 5, 6, 7, 8 ];
			$self->$primary(
				-score      => $score,
				-percent_id => $percent_id,
				-gen_start  => $gen_start,
				-gen_end    => $gen_end,
				-gen_id     => $gen_id,
				-est_start  => $est_start,
				-est_end    => $est_end,
				-est_id     => $est_id,
				-primary    => $primary
			);

		}
		elsif ( $_ =~ /Segmentation fault/ ) {
			close(ESTGENOME) or warning("problem closing est_genome: $!\n");
			$self->_deletefiles( $genfile, $estfile );
			throw("Segmentation fault from est_genome:\n<$_>\n");
		}
		elsif ( $_ =~ /(ERROR:.+)/ ) {
			close(ESTGENOME) or warning("problem closing est_genome: $!\n");
			$self->_deletefiles( $genfile, $estfile );
			throw("Error from est_genome:\n<$_>\n");
		}
	}
	foreach my $seg_array ( keys( %{ $self->{'_exons'} } ) ) {
		my $dnafp =
		  Bio::EnsEMBL::DnaDnaAlignFeature->new(
			-features => $self->{'_exons'}->{$seg_array} );
		$self->add_output($dnafp);
	}

	if ( !close(ESTGENOME) ) {
            my $err = ($!) ? "error $!" :
              sprintf("%s %d", $? & 127 ? (signal => $? & 127) : ('exit code', $? >> 8));
            $self->_deletefiles( $genfile, $estfile );
            throw("Problems running est_genome when closing pipe: $err\n");
	}
	$self->_deletefiles( $genfile, $estfile );

	return 1;

}

sub est_strand {

	my ( $self, $sign ) = @_;
	if ( defined($sign) && ( $sign eq '1' || $sign eq '-1' ) ) {
		$self->{'_est_strand'} = $sign;
	}
	return $self->{'_est_strand'};

}

sub gen_strand {

	my ( $self, $sign ) = @_;
	if ( defined($sign) && ( $sign eq '1' || $sign eq '-1' ) ) {
		$self->{'_gen_strand'} = $sign;
	}
	return $self->{'_gen_strand'};
}

sub Segment {    # named to match output from est2genome

	my ($self) = shift;
	my $p = {
		-score      => undef,
		-percent_id => undef,
		-gen_start  => undef,
		-gen_end    => undef,
		-gen_id     => undef,
		-est_start  => undef,
		-est_end    => undef,
		-est_id     => undef,
		-primary    => undef,
		@_
	};

	if ( $p->{-gen_end} < $p->{-gen_start} ) {
		( $p->{-gen_end}, $p->{-gen_start} ) =
		  ( $p->{-gen_start}, $p->{-gen_end} );
	}

	if ( $p->{-est_end} < $p->{-est_start} ) {
		( $p->{-est_end}, $p->{-est_start} ) =
		  ( $p->{-est_start}, $p->{-est_end} );
	}

	my $fp = $self->_create_FeaturePair(
		$p->{-score}, $p->{-percent_id},
		$p->{-gen_start}, $p->{-gen_end}, $p->{-gen_id},
		$p->{-est_start}, $p->{-est_end}, $p->{-est_id},
		$self->{'_source_tag'},
		$self->gen_strand,
		$self->est_strand,    # est_strand IS stored in the db now!
		$p->{-primary}
	);

	$self->_add_segment_to_exon( $fp, $p );

}

sub Exon {                    # named to match output from est2genome

	my ($self) = shift;
	my $p = {
		-score      => undef,
		-percent_id => undef,
		-gen_start  => undef,
		-gen_end    => undef,
		-gen_id     => undef,
		-est_start  => undef,
		-est_end    => undef,
		-est_id     => undef,
		@_
	};

#push( @{$self->{'_exon_lines'}} , [ $p->{-gen_start}, $p->{-gen_end}, $p->{-score}, $p->{-percent_id} ] );

	$self->_get_add_exon_lines(
		[ $p->{-gen_start}, $p->{-gen_end}, $p->{-score}, $p->{-percent_id} ] );

}

sub _add_segment_to_exon {

	my ( $self, $fp, $p ) = @_;
	my $exon_line_count = scalar( @{ $self->{'_exon_lines'} } );

	for ( my $i = 0 ; $i < $exon_line_count ; $i++ ) {

		# get current __exon data
		my ( $start, $end, $score, $pid ) =
		  ( @{ $self->{'_exon_lines'}->[$i] } );

		# check if feature pair belongs to current __exon
		next if ( $p->{-gen_start} < $start || $p->{-gen_start} > $end );

		# get last feature pair of current __exon for sanity check

		if ( my $prev_fp = $self->_get_last_fp_of_exon($i) ) {

	  # sanity check to cater for Bio::EnsEMBL::BaseAlignFeature->_parsefeatures

			if (   $fp->start eq ( $prev_fp->end + $fp->strand )
				|| $fp->end    eq ( $prev_fp->start + $fp->strand )
				|| $fp->hstart eq ( $prev_fp->hend + $fp->hstrand )
				|| $fp->hend   eq ( $prev_fp->hstart + $fp->hstrand ) )
			{
				$fp->score($score);
				$fp->percent_id($pid);
				$self->_add_fp_to_exon( $fp, $i );
			}
			else {
				$self->_change_exon_line( $i,
					[ $start, $prev_fp->end, $score, $pid ] );
				my $new_exon_lines =
				  $self->_get_add_exon_lines(
					[ $fp->start - 1, $end, $score, $pid ] );
				$fp->score($score);
				$fp->percent_id($pid);

				# this adds fp to new exon
				$self->_add_fp_to_exon( $fp,
					( scalar( @{$new_exon_lines} ) - 1 ) );
			}
		}
		else {
			$fp->score($score);
			$fp->percent_id($pid);
			$self->_add_fp_to_exon( $fp, $i );
		}
	}
}

sub _get_last_fp_of_exon {

	my $self = shift;
	my $exon = shift;
	my $last = 0;
	if ( defined($exon) ) {
		$self->{'_exons'}->{$exon} ||= [];
		$last = scalar( @{ $self->{'_exons'}->{$exon} } );
		$last = $last - 1 if $last;
	}

	return $self->{'_exons'}->{$exon}->[$last] || undef;
}

sub _add_fp_to_exon {

	my $self = shift;
	my $fp   = shift;
	my $exon = shift;

	if ( defined($exon) && defined($fp) ) {
		$self->{'_exons'}->{$exon} ||= [];
		push( @{ $self->{'_exons'}->{$exon} }, $fp );
	}

}

sub _get_add_exon_lines {

	my $self      = shift;
	my $exon_line = shift;

	if ( @{$exon_line} ) {
		push( @{ $self->{'_exon_lines'} }, $exon_line )
		  if scalar( @{$exon_line} ) == 4;
	}

	return $self->{'_exon_lines'};

}

sub _change_exon_line {

	my $self      = shift;
	my $index     = shift;
	my $exon_line = shift;
	if ( defined($index) && @{$exon_line} ) {
		$self->{'_exon_lines'}->[$index] = $exon_line;
	}
	return $self->{'_exon_lines'};
}

sub add_output {

	my ( $self, $fp ) = @_;
	push( @{ $self->{'_output'} }, $fp ) if defined($fp);
}

sub output {
	my ($self) = @_;
	return $self->{'_output'} || [];
}

sub convert_output {
	my ($self) = @_;
	my @genes;
	my @exons;
	my @supp_feat;

	# split the different features up
	foreach my $f ( @{ $self->{'_fplist'} } ) {
		print "Feature " . $f->gffstring . "\n";

		if ( $f->primary_tag eq 'Span' ) {
			push( @genes, $f );
		}
		elsif ( $f->primary_tag eq 'Exon' ) {
			push( @exons, $f );
		}
		elsif ( $f->primary_tag eq 'Segment' ) {
			push( @supp_feat, $f );
		}
	}

	# now reassemble them
	# add exons to genes
	foreach my $ex (@exons) {
		my $added = 0;

		foreach my $g (@genes) {
			if (   $ex->start >= $g->start
				&& $ex->end <= $g->end
				&& $ex->strand == $g->strand
				&& !$added )
			{
				$g->feature1->add_sub_SeqFeature( $ex, '' );
				$added = 1;
			}
		}
		warning("Exon $ex could not be added to a gene ...\n") unless $added;
	}

	# add supporting features to exons
	foreach my $sf (@supp_feat) {
		my $added = 0;
		foreach my $e (@exons) {
			if (   $sf->start >= $e->start
				&& $sf->end <= $e->end
				&& $sf->strand == $e->strand
				&& !$added )
			{
				$e->feature1->add_sub_SeqFeature( $sf, '' );
				$added = 1;
			}
		}
		warning("Feature $sf could not be added to an exon ...\n")
		  unless $added;
	}
	push( @{ $self->{'_output'} }, @genes );

}

sub _create_FeaturePair {

	my (
		$self,     $f1score,  $f1percent_id, $f1start, $f1end,
		$f1id,     $f2start,  $f2end,        $f2id,    $f1source,
		$f1strand, $f2strand, $f1primary
	  )
	  = @_;

	my $analysis_obj = new Bio::EnsEMBL::Analysis(
		-db              => "none",
		-db_version      => "none",
		-program         => "est_genome",
		-program_version => "none",
		-gff_source      => $f1source,
		-gff_feature     => $f1primary,
	);

	#create features
	my $feat1 = new Bio::EnsEMBL::SeqFeature(
		-start       => $f1start,
		-end         => $f1end,
		-seqname     => $f1id,
		-strand      => $f1strand,
		-score       => $f1score,
		-percent_id  => $f1percent_id,
		-source_tag  => $f1source,
		-primary_tag => $f1primary,
		-analysis    => $analysis_obj
	);
	my $feat2 = new Bio::EnsEMBL::SeqFeature(
		-start       => $f2start,
		-end         => $f2end,
		-seqname     => $f2id,
		-strand      => $f2strand,
		-score       => $f1score,
		-percent_id  => $f1percent_id,
		-source_tag  => $f1source,
		-primary_tag => $f1primary,
		-analysis    => $analysis_obj
	);

	#create featurepair
	my $fp = new Bio::EnsEMBL::FeaturePair(
		-feature1 => $feat1,
		-feature2 => $feat2
	);
	return $fp;
}

sub _createfeatures {
	my (
		$self,     $f1score,  $f1start,   $f1end,    $f1id,
		$f2start,  $f2end,    $f2id,      $f1source, $f2source,
		$f1strand, $f2strand, $f1primary, $f2primary
	  )
	  = @_;

	my $analysis_obj = new Bio::EnsEMBL::Analysis(
		-db              => "none",
		-db_version      => "none",
		-program         => "est_genome",
		-program_version => "none",
		-gff_source      => $f1source,
		-gff_feature     => $f1primary,
	);

	#create features
	my $feat1 = new Bio::EnsEMBL::SeqFeature(
		-start       => $f1start,
		-end         => $f1end,
		-seqname     => $f1id,
		-strand      => $f1strand,
		-score       => $f1score,
		-percent_id  => $f1score,
		-source_tag  => $f1source,
		-primary_tag => $f1primary,
		-analysis    => $analysis_obj
	);

	my $feat2 = new Bio::EnsEMBL::SeqFeature(
		-start       => $f2start,
		-end         => $f2end,
		-seqname     => $f2id,
		-strand      => $f2strand,
		-score       => $f1score,
		-percent_id  => $f1score,
		-source_tag  => $f2source,
		-primary_tag => $f2primary,
		-analysis    => $analysis_obj
	);

	#create featurepair
	my $fp = new Bio::EnsEMBL::FeaturePair(
		-feature1 => $feat1,
		-feature2 => $feat2
	);

	push( @{ $self->{'_fplist'} }, $fp );
}

sub _createfiles {
	my ( $self, $genfile, $estfile, $dirname ) = @_;

	#check for diskspace
	my $spacelimit = 0.01;       # 0.01Gb or about 10 MB
	                             #my $dir ="./";
	my $dir        = $dirname;
	unless ( $self->_diskspace( $dir, $spacelimit ) ) {
		throw("Not enough disk space ($spacelimit Gb required)");
	}

	#if names not provided create unique names based on process ID
	$genfile = $self->_getname("genfile") unless ($genfile);
	$estfile = $self->_getname("estfile") unless ($estfile);
	$genfile = $dirname . $genfile;
	$estfile = $dirname . $estfile;

	# Should check we can write to this directory
	throw("No directory $dirname") unless -e $dirname;

	#chdir ($dirname) or throw ("Cannot change to directory '$dirname' ($?)");
	return ( $genfile, $estfile );
}

sub _getname {
	my ( $self, $typename ) = @_;
	return $typename . "_" . $$ . ".fn";
}

sub _diskspace {
	my ( $self, $dir, $limit ) = @_;
	my $Gb   = 1024**3;

	open DF, "df -k $dir |" || throw("FAILED to open 'df' pipe ".
                                   "Runnable::diskspace : $!\n");
  	my $count = 0;
  	my $status = 1;
	while (<DF>) {
		if($count && $count > 0){
      		my @values = split;
      		my $space_in_Gb = $values[3] * 1024 / $Gb;
      		$status = 0 if ($space_in_Gb < $limit);
    	}
    	$count++;
	}
	close DF || throw("FAILED to close 'df' pipe ".
                    "Runnable::diskspace : $!\n");

	return $status;
}

sub _deletefiles {
	my ( $self, @files ) = @_;
	my $unlinked = unlink(@files);

	if ( $unlinked == @files ) {
		return 1;
	}
	else {
		my @fails = grep -e, @files;
		warning("Failed to remove @fails : $!\n");
	}
}

sub est_genome {
	my ( $self, $arg ) = @_;

	if ( defined($arg) ) {
		$self->{'_est_genome'} = $arg;
	}
	return $self->{'_est_genome'};
}

1;
