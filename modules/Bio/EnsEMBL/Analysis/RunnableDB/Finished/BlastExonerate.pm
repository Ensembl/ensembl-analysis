# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
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


package Bio::EnsEMBL::Analysis::RunnableDB::Finished::BlastExonerate;

use warnings ;
use strict;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::Root::RootI;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_Pfetch;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Getseqs;
use Bio::EnsEMBL::Analysis::Runnable::Finished::BlastExonerate;
use Bio::EnsEMBL::Pipeline::SeqFetcher::xdget;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Analysis::Config::Blast;
use Bio::PrimarySeq;
use base 'Bio::EnsEMBL::Analysis::RunnableDB::Finished';

sub new {
	my ( $new, @args ) = @_;
	my $self = $new->SUPER::new(@args);

	# dbobj, input_id, seqfetcher, and analysis objects are all set in
	# in superclass constructor (RunnableDB.pm)
	# Also set the BLAST PARAMETERS from the Blast Config file
	$self->read_and_check_config($BLAST_CONFIG_BY_LOGIC);
	return $self;
}

sub require_module {
	my ( $self, $module ) = @_;
	my $class;
	( $class = $module ) =~ s/::/\//g;
	eval { require "$class.pm"; };
	throw( "Couldn't require " . $class . " Blast:require_module $@" ) if ($@);
	return $module;
}

sub read_and_check_config {
	my ( $self, $var_hash ) = @_;

	#call RunnableDBs method to fill in values first
	$self->SUPER::read_and_check_config($var_hash);

	#now for type checking and other sanity checks

	#must have a parser object and to pass to blast
	throw(  "BLAST_PARSER must be defined either in the DEFAULT entry or in "
		  . "the hash keyed on "
		  . $self->analysis->logic_name
		  . " Blast::read_and_check_config" )
	  if ( !$self->BLAST_PARSER );
	$self->require_module( $self->BLAST_PARSER );

	#load the filter module if defined
	if ( $self->BLAST_FILTER ) {
		$self->require_module( $self->BLAST_FILTER );
	}

	#if any of the object params exist, all are optional they must be hash
	#refs
	throw(  "PARSER_PARAMS must be a hash ref not "
		  . $self->PARSER_PARAMS
		  . " Blast::read_and_check_config" )
	  if ( $self->PARSER_PARAMS && ref( $self->PARSER_PARAMS ) ne 'HASH' );
	throw(  "FILTER_PARAMS must be a hash ref not "
		  . $self->FILTER_PARAMS
		  . " Blast::read_and_check_config" )
	  if ( $self->FILTER_PARAMS && ref( $self->FILTER_PARAMS ) ne 'HASH' );
	throw(  "BLAST_PARAMS must be a hash ref not "
		  . $self->BLAST_PARAMS
		  . " Blast::read_and_check_config" )
	  if ( $self->BLAST_PARAMS && ref( $self->BLAST_PARAMS ) ne 'HASH' );

	my $blast_params;
	if ( $self->BLAST_PARAMS ) {
		$blast_params = $self->BLAST_PARAMS;
	}
	else {
		$blast_params = {};
	}
	my %parameters;
	if ( $self->parameters_hash ) {
		%parameters = %{ $self->parameters_hash };
	}
	foreach my $key (%parameters) {
		$blast_params->{$key} = $parameters{$key};
	}
	$self->BLAST_PARAMS($blast_params);
}

sub BLAST_PARSER {
	my ( $self, $value ) = @_;

	if ( defined $value ) {
		$self->{'_CONFIG_BLAST_PARSER'} = $value;
	}

	return $self->{'_CONFIG_BLAST_PARSER'};
}

sub PARSER_PARAMS {
	my ( $self, $value ) = @_;

	if ( defined $value ) {
		$self->{'_CONFIG_PARSER_PARAMS'} = $value;
	}

	return $self->{'_CONFIG_PARSER_PARAMS'};
}

sub BLAST_FILTER {
	my ( $self, $value ) = @_;

	if ( defined $value ) {
		$self->{'_CONFIG_BLAST_FILTER'} = $value;
	}

	return $self->{'_CONFIG_BLAST_FILTER'};
}

sub FILTER_PARAMS {
	my ( $self, $value ) = @_;

	if ( defined $value ) {
		$self->{'_CONFIG_FILTER_PARAMS'} = $value;
	}

	return $self->{'_CONFIG_FILTER_PARAMS'};
}

sub BLAST_PARAMS {
	my ( $self, $value ) = @_;

	if ( defined $value ) {
		$self->{'_CONFIG_BLAST_PARAMS'} = $value;
	}

	return $self->{'_CONFIG_BLAST_PARAMS'};
}

sub fetch_input {
	my ($self) = @_;
	my @fps;
	throw("No input id") unless defined( $self->input_id );
	my $sliceid = $self->input_id;
	my $sa      = $self->db->get_SliceAdaptor();
	my $slice   = $sa->fetch_by_name($sliceid);
	$slice->{'seq'} = $slice->seq();
	$self->query($slice);
	my $hard_masked_slice = $slice->get_repeatmasked_seq( $ANALYSIS_REPEAT_MASKING, 0 )
	  or throw("Unable to fetch contig");
	my $soft_masked_slice = $slice->get_repeatmasked_seq( $ANALYSIS_REPEAT_MASKING, 1 )
	  or throw("Unable to fetch contig");
	my $hardmaskedseq = $hard_masked_slice->seq();

	if ( scalar( $hardmaskedseq =~ s/([CATG])/$1/g ) > 3 ) {
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

	my $runnable = Bio::EnsEMBL::Analysis::Runnable::Finished::BlastExonerate->new(
		'-hard_masked'    => $hard_masked_slice,
		'-soft_masked'	  => $soft_masked_slice,
		'-analysis' => $self->analysis,
	);
	$self->runnable($runnable);
	return 1;
}

sub check_with_seg {
	my ( $self, $seqObj_to_test ) = @_;

	warn "need a Bio::Seq Obj" unless $seqObj_to_test;

	my ($filename) = $self->_createfiles( '/tmp', [qw(seg_checking)] );
	my $file = Bio::SeqIO->new(
		-file   => ">$filename",
		-format => 'Fasta'
	  )
	  or throw("Can't create Bio::SeqIO $filename $!");
	$file->write_seq($seqObj_to_test);

	my $seg_cmd = "nseg $filename -x";
	my $seg     = Bio::SeqIO->new(
		-file   => "$seg_cmd |",
		-format => 'Fasta'
	  )
	  or throw("Can't create Bio::SeqIO $seg_cmd $!");
	my $seq;
	eval { $seq = $seg->next_seq->seq; };
	unlink($filename);

	if ($@) {
		throw(
			"There was a problem with SEG masking.\nI tried to '$seg_cmd'");
	}
	if ( $seq =~ /[CATG]{3}/i ) {
		$self->input_is_void(0);
	}
	else {
		$self->input_is_void(1);
		warning("Need at least 3 nucleotides after SEG filtering");
	}
}

sub _createfiles {
	my ( $self, $dirname, $filenames ) = @_;
	my $unique = {};
	$unique = { map { $_, $unique->{$_}++ } @$filenames };
	my @files = ();

	$dirname ||= '/tmp';
	$dirname =~ s!(\S+)/$!$1!;

	foreach my $file (@$filenames) {
		if ( $unique->{$file} ) {

			#name not unique add random
			$file .= ".$$." . int( rand(200) );
			push( @files, "$dirname/$file" );
		}
		else {

			#name was unique just add it
			push( @files, "$dirname/$file.$$" );
		}
	}

	return @files;
}


sub run {
	my ($self) = @_;
	my $runnables = $self->runnable;
	print scalar(@$runnables)," runnables in BlastExonerate module\n";
	foreach my $runnable(@$runnables){
		$runnable || throw("Can't run - no runnable object");
		my $blast  = $self->BLAST_PARAMS;
		my $parser = $self->make_parser;
		my $filter;
		if ( $self->BLAST_FILTER ) {
			$filter = $self->make_filter;
		}
		print "parser $parser\n";
		print "filter $filter\n";
		print "blast $blast\n";
		eval { $runnable->run( $parser, $filter, $blast ); };
		if ( my $err = $@ ) {
			chomp $err;
			$self->failing_job_status($1)
			  if $err =~ /^\"([A-Z_]{1,40})\"$/i
			  ;    # only match '"ABC_DEFGH"' and not all possible throws
			throw("$@");
		}
		$self->output($runnable->output);
	}
	1;
}

=head2 db_version_searched

Est2genome

-options=>W=15 T=25 -hitdist=40 E=1e-20 M=1 N=-3 Q=3 R=3 B=100000 Z=500000000 cpus=1 -wordmask=seg,-exo_options=>-m e2g --forcescan q --softmasktarget yes  -M 1500 --dnahspthreshold 120 -s 2000 --geneseed 300,-query_type=>dna

-putenv "BLASTFILTER=/usr/local/ensembl/bin/"

Protein2genome

-options=>W=4  T=20 E=1e-4  B=100000 Z=500000000 cpus=1 -wordmask=seg,-exo_options=>-m p2g --forcescan q --softmasktarget yes  -M 1500 --dnahspthreshold 120 -s 2000 --geneseed 300,-query_type=>protein

 -putenv "BLASTMAT=/usr/local/ensembl/data/blastmat/"
 -putenv "BLASTFILTER=/usr/local/ensembl/bin/"


    Title   :  db_version_searched
               [ distinguished from Runnable::*::get_db_version() ]
    Useage  :  $self->db_version_searched('version string')
               $obj->db_version_searched()
    Function:  Get/Set a blast database version that was searched
               The actual look up is done in Runnable::Finished_Blast
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


sub make_parser {
	my ( $self, $hash ) = @_;
	if ( !$hash ) {
		$hash = $self->PARSER_PARAMS;
	}
	my %parser = %$hash;
	my $parser = $self->BLAST_PARSER->new(%parser);
	return $parser;
}

=head2 make_filter

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Blast
  Arg [2]   : hashref, parameters for filter constructor
  Function  : create a filter object
  Returntype: a filter object
  Exceptions:
  Example   :

=cut

sub make_filter {
	my ( $self, $hash ) = @_;
	if ( !$hash ) {
		$hash = $self->FILTER_PARAMS;
	}
	my %filter = %$hash;
	my $filter = $self->BLAST_FILTER->new(%filter);
	return $filter;
}

1;


=head1 NAME - Bio::EnsEMBL::Analysis::RunnableDB::Finished::BlastExonerate

=head1 AUTHOR

Mustapha Larbaoui B<email> ml6@sanger.ac.uk
