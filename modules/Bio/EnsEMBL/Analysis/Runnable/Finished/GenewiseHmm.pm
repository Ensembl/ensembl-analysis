
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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
Bio::EnsEMBL::Analysis::Runnable::Finished::GenewiseHmm 
- Hmms to genomic sequence

=head1 SYNOPSIS

    # Initialise the GenewiseHmm module  
    my $genewise = new  Bio::EnsEMBL::Analysis::Runnable::Finished::GenewiseHmm  ('query'  => $genomic,
								       'hmmfile'  => $hmmfile,
								       'memory'   => $memory);

   $genewise->run;
   my @genes = $blastwise->output

=head1 DESCRIPTION

This module is created and run by halfwisehmm.pm, it takes a pfam hmm and some genomic sequence and runs genewise using these then parses the output to produce
 feature pairs representing the genes with exons as subseq features and supporting evidence as subseqfeatures of the exons

=head1 CONTACT

refactored by Sindhu K. Pillai <sp1@sanger.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Analysis::Runnable::Finished::GenewiseHmm;
use warnings ;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use strict;
use vars qw(@ISA);

# Object preamble - inherits from Bio::EnsEMBL::Root

use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Analysis::Programs qw(/usr/local/ensembl/bin/genewisedb);
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::SeqIO;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

{
	my $options = {
		'Algorithm type:        ' => undef,
		'Search algorithm used: ' => undef,
		'Implementation:        ' => undef,
		'Search mode:           ' => undef,
		'Protein info from:     ' => undef,
		'Dna info from:         ' => undef,
		'Start/End (protein)    ' => undef,
		'Gene Paras:            ' => undef,
		'Codon Table:           ' => undef,
		'Subs error:            ' => undef,
		'Indel error:           ' => undef,
		'Model splice?          ' => undef,
		'Model codon bias?      ' => undef,
		'Model intron bias?     ' => undef,
		'Null model             ' => undef,
		'Alignment Alg          ' => undef
	};

	sub valid_options {
		return $options;
	}
}

=head2 new

  Args      : Bio::Seq, int, boolean, boolean, path, filename
  Function  : create a Finished_GenewiseHmm object
  Returntype: object created
  Exceptions: throws if give no query sequence or no hmmfile
  Caller    :
  Example   :

=cut

sub new {
	my ( $class, @args ) = @_;
	my $self = $class->SUPER::new(@args);
	$self->{'_output'}  = [];
	$self->{'_reverse'} = undef;
	my ( $query, $memory, $reverse, $endbias, $genewise, $hmmfile, $options ) =
	  rearrange( [qw(QUERY MEMORY REVERSE ENDBIAS GENEWISE HMMFILE OPTIONS)],
		@args );

	#  print $query . "\n";
	$genewise ||= 'genewise';
	$options  ||= '-ext 2 -genes';
	$self->genewise( $self->program($genewise) );
	$self->query($query) || throw("No query sequence entered for blastwise");
	$self->hmmfile($hmmfile) || throw("No Hmm file entered for Hmmgenewise");
	$self->is_reverse($reverse) if ( defined($reverse) );
	$self->endbias($endbias)    if ( defined($endbias) );
	$self->memory($memory)      if ( defined($memory) );
	$self->options($options);
	return $self;
}

####################
# accessor methods #
####################

=head2 genewise

  Arg [1]   : 
  Function  : 
  Returntype: genewise value
  Exceptions: 
  Caller    : 
  Example   : 

=cut

sub genewise {
	my ( $self, $arg ) = @_;
	$self->{_genewise} = $arg if $arg;
	return $self->{_genewise};
}

# These all set/get or initializing methods

sub is_reverse {
	my ( $self, $boolean ) = @_;
	$self->{'_reverse'} = ( $boolean == 1 ? $boolean : 0 ) if $boolean;
	return $self->{'_reverse'};
}

sub endbias {
	my ( $self, $boolean ) = @_;
	$self->{'_endbias'} = ( $boolean == 1 ? $boolean : 0 ) if $boolean;
	return $self->{'_endbias'};
}

sub memory {
	my ( $self, $arg ) = @_;
	$self->{'_memory'} = ' -kbyte ' . $arg if $arg;
	return $self->{'_memory'} || ' -kbyte 100000';
}

sub query {
	my ( $self, $arg ) = @_;

	if ( defined($arg) ) {
		throw(
"Genomic sequence input is not a Bio::SeqI or Bio::Seq or Bio::PrimarySeqI"
		  )
		  unless ( $arg->isa("Bio::SeqI")
			|| $arg->isa("Bio::Seq")
			|| $arg->isa("Bio::PrimarySeqI") );
		$self->{'_query'} = $arg;
	}
	return $self->{'_query'};
}

=head2 hmmfile

 Title   : hmmfile
 Usage   : $obj->hmmfile($newval)
 Function: 
 Returns : value of hmmfile
 Args    : newvalue (optional)


=cut

sub hmmfile {
	my ( $self, $file ) = @_;
	$self->{'hmmfile'} = $file if $file;
	return $self->{'hmmfile'};
}

######################
#running the analysis#
######################

=head2 set_environment

  Arg [1]   : none
  Function  : check that the WISECONFIGDIR is set
  Returntype: none
  Exceptions: thows if variable isn't a directory'
  Caller    : 
  Example   : 

=cut

sub set_environment {
	my ($self) = @_;
	if ( !-d $ENV{WISECONFIGDIR} ) {
		throw( "No WISECONFIGDIR [" . $ENV{WISECONFIGDIR} . "]" );
	}
}

=head2 run

  Arg [1]   : path to directory outputfiles are to be written too 
  Function  : calls methods which run genewise
  Returntype: none
  Exceptions: none
  Caller    : 
  Example   : 

=cut

sub run {
	my ( $self, $dir ) = @_;
	$self->set_environment;
	$self->_align_protein($dir);
}

sub parse_headers {
	my $self = shift;
	my @header = split( "\n", +shift );
	local $" = "\n";
	#print STDERR "HEADER: @header";
	my $information = $self->valid_options();
	foreach my $line (@header) {
		if ( $line =~ m!([a-zA-Z0-9 :?\(\)/]{23})(.+)! ) {
			if ( exists( $information->{$1} ) ) {
				$self->gw_options( $1, $2 );
			}
			else {

				#		print STDERR "NOT OPTION: $line\n";
			}
		}
		elsif ( $line =~ m!^Protein ! ) {
			my @high_scores = split( /\s+/, $line );

			#	    print join(" *** ", @high_scores) . "\n";
			$self->set_score(@high_scores);
		}
		else {

			#	    print STDERR "BLANK LINE: $line\n";
		}
	}
}

sub parse_alignments {
	my $self    = shift;
	my @results = split( "\n", +shift );

	#    local $" = "\n";
	#print STDERR "RESULTS: @results";

	#    local $" = " ";
	my $firstline = shift(@results);
	pop(@results) if $results[-1] eq '//';
	my ( $pdomain, $clone, $strand, $count, $gff_strand );

	if ( $firstline =~
		m!Results for (.+) vs (.+) \((forward|reverse)\) \[(\d)\]! )
	{
		$pdomain    = $1;
		$clone      = $2;
		$strand     = ( $3 eq 'forward' ? 1 : -1 );
		$gff_strand = ( $3 eq 'forward' ? '+' : '-' );
		$count      = $4;
	}
	else {
		warning("couldn't read firstline <$firstline>");
	}
	my $score = $self->get_score( $pdomain, "[$gff_strand]" );
	my @genes;
	my $allowed = { 'Gene' => 1, 'Exon' => 1, 'Supporting' => 1 };
	foreach (@results) {
		my @F = split;
		next unless ( defined $F[0] && $allowed->{ $F[0] } );
		my $check = 1;    # positive strand
		if ( $F[0] eq 'Gene' && $F[1] && $F[2] && !( $F[1] eq 'Paras:' ) ) {

			# swap reverse coords (make note this was done)
			( $F[1], $F[2], $check ) = ( $F[2], $F[1], -1 ) if $F[2] < $F[1];

			# sanity check strands and add gene
			if ( $check == $strand ) {
				my $gene = Bio::EnsEMBL::SeqFeature->new();
				$gene->id($clone);
				$gene->seqname($pdomain);
				push( @genes, $gene );
			}
			else { warning "strands don't match" }
		}
		elsif ( $F[0] eq 'Exon' ) {

			# swap reverse coords (make note this was done)
			( $F[1], $F[2], $check ) = ( $F[2], $F[1], -1 ) if $F[2] < $F[1];
			if ( $check == $strand ) {
				my $exon = Bio::EnsEMBL::SeqFeature->new();
				$exon->seqname($clone);
				$exon->id($pdomain);

				#$exon->phase  (0);
				$exon->start( $F[1] );
				$exon->end( $F[2] );
				$exon->strand($strand);
				$genes[-1]->add_sub_SeqFeature( $exon, 'EXPAND' );

				#[@F, $pdomain, $clone, $strand, $phase, $score]);
			}
			else { warning "strands don't match" }
		}
		elsif ( $F[0] eq 'Supporting' ) {

			# swap reverse coords (make note this was done)
			( $F[1], $F[2], $check ) = ( $F[2], $F[1], -1 ) if $F[2] < $F[1];
			if ( $check == $strand ) { }
		}
		else { }
	}

	return \@genes;
}

=head2 _align_protein

  Arg [1]   : path to directory output files are to be written too  
  Function  : runs genewise and parse output producing arrays of genes, exons and supporting features
  Returntype: none
  Exceptions: throws if problems opening or closing pipe to genewise 
  Caller    : 
  Example   : 

=cut

sub _align_protein {
	my ( $self, $dir ) = @_;
	my $memory = $self->memory;
	$self->workdir('/tmp') unless ( $self->workdir($dir) );
	$self->checkdir();

	#my $gwfile   = "gw." . $$ . ".out";
	my $genfile = "gw." . $$ . ".gen.fa";
	$self->write_fasta($genfile);
	my $hmm      = $self->hmmfile;
	my $genewise = $self->genewise;
	$self->files_to_delete($genfile);
	$self->files_to_delete($hmm);
	my $command  =
	  "$genewise $hmm $genfile " . $self->options() . $self->memory();
	$command .= " -init endbias -splice flat "
	  if ( $self->endbias && $self->endbias == 1 );
	$command .= " -trev " if ( $self->is_reverse && $self->is_reverse == 1 );
	print STDERR "Command is $command\n";    ##########
	    #my $outputfile = $self->hmmfile.".output";
	    #my $testing_output = "genewisedb.output";
	local $/ = "//\n";

# use one of these
#open(my $fh, "$command | tee $outputfile |") or throw("error piping to genewise: $!\n"); # test genewise keep output
#open(my $fh, "$testing_output") or throw('couldnt find the file'); # test a single file
	open( my $fh, "$command 2>&1 |" )
	  or throw("error piping to genewise: $!\n");    # production
	 #print "parseing output\n"; # making assumption of only 1 gene prediction ... this will change once we start using hmms

	while (<$fh>) {
		#    print STDERR "*" x 60 . "\n";
		if(/Out of memory/) {
			$self->delete_files;	
			die qq{"OUT_OF_MEMORY"\n};
		}
		my ( $header, $results ) = split( "\n>", $_ );
		$self->parse_headers($header) if $header;
		my $genes = $self->parse_alignments($results);
		$self->addGenes($genes);

		#    print STDERR "*" x 60 . "\n";
	}
	close($fh) or throw("Error running genewise:$!\n");
	print "there are ",
	  scalar( $self->output )
	  . " genes predicted using hmms in file "
	  . $hmm . "\n";
	$self->delete_files;
	#unlink $gwfile;
	#unlink $outputfile;
}

sub write_fasta {
	my ( $self, $genfile ) = @_;
	my $genio = Bio::SeqIO->new(
		'-file'   => ">$genfile",
		'-format' => 'fasta'
	);
	$genio->write_seq( $self->query );
	$genio = undef;
}

=head2 addGene

  Arg [1]   : Bio:EnsEMBL::SeqFeature
  Function  : adds the seqfeature to the output array
  Returntype: none
  Exceptions: throws if not passed a seqfeature
  Caller    : 
  Example   : 

=cut

sub addGene {
	my ( $self, $arg ) = @_;

	#print STDERR "adding ".$arg->seqname." to ouput\n";
	if ( $arg->isa("Bio::EnsEMBL::SeqFeature") ) {
		push( @{ $self->{'_output'} }, $arg );
	}
	else {
		throw("this, $arg, should be a seqfeature\n");
	}

	#foreach my $result(@{$self->{'_output'}}){
	#  print STDERR $result->seqname."\n";
	#}
	#print STDERR "there are ".scalar(@{$self->{'_output'}})." genes\n";
}

=head2 addGenes

  Arg [1]   : array ref of Bio:EnsEMBL::SeqFeature
  Function  : adds the seqfeatures to the output array like addGene does for singles
  Returntype: none
  Exceptions: throws if not passed a seqfeature
  Caller    : 
  Example   : 

=cut

sub addGenes {
	my ( $self, $args ) = @_;
	foreach my $arg (@$args) {
		if ( $arg->isa("Bio::EnsEMBL::SeqFeature") ) {
			push( @{ $self->{'_output'} }, $arg );
		}
		else {
			throw("this, $arg, should be a seqfeature\n");
		}
	}
}

################
#output methods#
################

=head2 output

  Arg [1]   : none
  Function  : returns the output array
  Returntype: 
  Exceptions: none
  Caller    : 
  Example   : 

=cut

sub output {
	my ($self) = @_;
	return @{ $self->{'_output'} };
}

sub gw_options {
	my ( $self, $option_name, $option_value ) = @_;
	unless ( $self->{'_gw_options'} ) {
		$self->{'_gw_options'} = $self->valid_options;
	}
	$self->{'_gw_options'}->{$option_name} = $option_value if $option_value;
	return $self->{'_gw_options'}->{$option_name};
}

sub set_score {
	my ( $self, @line ) = @_;
	$self->{'_scores'}->{ $line[1] }->{ $line[3] } = $line[5];
}

sub get_score {
	my ( $self, $p, $s ) = @_;
	return $self->{'_scores'}->{$p}->{$s};
}

1;
