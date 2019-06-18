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

=head1 NAME - Bio::EnsEMBL::Analysis::Runnable::Finished_EST

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

Modified by Sindhu K. Pillai B<email> sp1@sanger.ac.uk
=cut

package Bio::EnsEMBL::Analysis::Runnable::Finished::EST;

use strict;
use warnings;

#use vars qw(@ISA);

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Analysis::Runnable::Finished::MiniEst2Genome;
use Bio::EnsEMBL::Analysis::Runnable::Finished::Blast;
use Bio::EnsEMBL::Analysis::Config::Blast;
use Bio::EnsEMBL::Analysis::Tools::BPliteWrapper;
use Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_xdget;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Storable;
use base 'Bio::EnsEMBL::Analysis::Runnable';


#use vars qw(@ISA $verbose);
#@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

use vars qw($verbose);
$verbose = 0;
my $debug = 0;

sub new {
    my ( $class, @args ) = @_;
    my $self = bless {}, $class;
    my ( $query, $unmasked, $analysis, $seqfetcher ) = rearrange(
        [
        qw{
            QUERY
            UNMASKED
            ANALYSIS
            SEQFETCHER
        }
        ],
        @args
    );
    die "No QUERY (masked genomic sequence) given" unless $query;
    die "No UNMASKED (genomic sequence) given"     unless $unmasked;
    $self->query($query);
    $self->unmasked($unmasked);
    $self->analysis($analysis);

    # this should make a Finished_xdget SeqFetcher, but will default to Pfetch
    my $indices;
	$indices = [ split(",", $self->analysis->db) ] if($self->analysis->db);
	$seqfetcher = $self->_make_seqfetcher($indices);
	$self->seqfetcher($seqfetcher);

    return $self;
}

sub query {
    my ( $self, $seq ) = @_;
    if ($seq) {
        unless ( $seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::Seq") ) {
            throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
        }
        $self->{'_query'} = $seq;
    }
    return $self->{'_query'};
}

sub unmasked {
    my ( $self, $seq ) = @_;
    if ($seq) {
        unless ( $seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::Seq") ) {
            die Dumper($seq);
            throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
        }
        $self->{'_unmasked'} = $seq;
    }
    return $self->{'_unmasked'};
}

sub seqfetcher {
    my ( $self, $seqfetcher ) = @_;
    if ($seqfetcher) {
        $self->{'_seqfetcher'} = $seqfetcher;
    }
    return $self->{'_seqfetcher'};
}

sub seq_cache {
    my ( $self, $seq_cache ) = @_;
    if ($seq_cache) {
        $self->{'_seq_cache'} = $seq_cache;
    }
    return $self->{'_seq_cache'};
}

sub analysis {
    my ( $self, $analysis ) = @_;
    if ($analysis) {
        $self->{'_analysis'} = $analysis;
    }
    return $self->{'_analysis'};
}

sub run {

    my ($self) = shift;
    my $parser = shift;
    my $filter = shift;
    my $blastcon = shift;


    my $blast = Bio::EnsEMBL::Analysis::Runnable::Finished::Blast->new(
     -query => $self->query,
     -program => 'wublastn',
     -database => $self->analysis->db_file,
     -analysis => $self->analysis,
     -parser => $parser,
     -filter => $filter,
     %{$blastcon},
     );

    my $db_files = $blast->databases;

    # set the db_version_searched
    my $dbv = $blast->get_db_version();
    print "using $dbv\n";
    $self->get_db_version($dbv);

    DB_FILE: foreach my $db_file (sort @$db_files) {
    	print STDOUT "database file $db_file\n";
    	print STDOUT "DEBUG ".localtime().": Running blast against $db_file\n" if $debug;
    	$blast->clean_output();
    	$blast->clean_results_files();
    	$blast->clean_databases();
    	$blast->databases($db_file);
	    $blast->run();
	    my $features =  $blast->output ;
	    next DB_FILE unless @$features;
	    my %acc_hash = map{ $_->hseqname => 1 } @$features;
		my @accessions = keys %acc_hash;

		if (!$self->analysis->db)
		{
			# read the blast file and the wu-blast indexes to pre-warm disk cache
			system("cat $db_file ${db_file}.x* > /dev/null");
			my $seqfetcher = $self->_make_seqfetcher([$db_file]);
	    	$self->seqfetcher($seqfetcher);
		}

		# Fetch all the hit sequences from the blast db using xdget
		print STDOUT "DEBUG ".localtime().": Fetching ".scalar(@accessions)." sequences\n" if $debug;
		$self->seq_cache($self->seqfetcher->get_Seq_by_accs(\@accessions));
		print STDOUT "DEBUG ".localtime().": Storing ".scalar(@accessions)." sequences in /tmp/$$.desc\n" if $debug;
		# Store the sequence descriptions and lengths in a hash
		# then serialize it for subsequent use by Finished.pm
	  	# Advantage: Call xdget once => Keep Lustre FS load low => Everybody happy :-)
	  	my %descriptions = ();
	  	my $persistent = "/tmp/$$.desc";
	  	if(-e "$persistent") {
	  		my $hash = retrieve "$persistent";
	  		@descriptions{ keys %$hash } = values %$hash;
	  	}
		foreach(@accessions) {
			next unless $self->seq_cache->{$_};
			$descriptions{$_}{description} = $self->seq_cache->{$_}->description();
	    	$descriptions{$_}{length} = $self->seq_cache->{$_}->length();
		}
	  	store \%descriptions, "$persistent";

		print STDOUT "DEBUG ".localtime().": Running Est2Genome on Plus strand\n" if $debug;
	    $self->run_est_genome_on_strand( 1, $features );
		print STDOUT "DEBUG ".localtime().": Running Est2Genome on Minus strand\n" if $debug;
	    $self->run_est_genome_on_strand( -1, $features)
    }
}

sub run_est_genome_on_strand {

    my ( $self, $strand, $feat ) = @_;
    my $hit_features = {};
    for ( my $i = 0 ; $i < @$feat ; $i++ ) {

      my $f   = $feat->[$i];
      my $hid = $f->hseqname or throw("Missing hid");
      next unless $f->strand == $strand;
      $hit_features->{$hid} ||= [];
      push ( @{ $hit_features->{$hid} }, $f );
    }
    my ($is_linear);
    if ( $strand == 1 ) {
        $is_linear = sub {
            my ( $x, $y ) = @_;
            return $x->hstart <= $y->hstart;
        };
    }
    else {
        $is_linear = sub {
            my ( $x, $y ) = @_;
            return $x->hend >= $y->hend;
        };
    }


    while ( my ( $hid, $flist ) = each %$hit_features ) {

        print STDERR "$hid\n" if $verbose;
        print STDERR "PRESORTED\n" if $verbose; #
        foreach my $f (@$flist) {
            print STDERR "hit:  " . $f->gffstring . "\n" if $verbose;

        }
        # Sort features by start/end in genomic
        @$flist = sort { $a->start <=> $b->start or $a->end <=> $b->end } @$flist;
        print STDERR "SORTED\n" if $verbose;
        foreach my $f (@$flist) {
            print STDERR "hit:  " . $f->gffstring . "\n" if $verbose;

        }

        # Group into linear matches
        my @sets = ( [ $flist->[0] ] );
        my $curr = $sets[0];

        for ( my $i  = 1 ; $i < @$flist ; $i++ ) {
            my $prev = $flist->[ $i - 1 ];
            my $this = $flist->[$i];
            if ( &$is_linear( $prev, $this ) ) {
                push ( @$curr, $this );

            }
            else {
                $curr = [$this];
                push ( @sets, $curr );

            }

        }

	print STDERR "MADE ",scalar(@sets)," LINEAR MATCHES SET(S)\n" if $verbose;

	foreach my $lin (@sets) {
	  $self->do_mini_est_genome($lin);
        }

      }
}

sub do_mini_est_genome {

    my ( $self,$linear ) = @_;
    print STDERR "\nLinear Matches\n" if $verbose;
    foreach my $lin (@$linear) {
        print STDERR "before:  " . $lin->gffstring . "\n" if $verbose;
    }

    ### Is this merging features?  - May be cause of duplicate features if isn't
    ### Duplicate features found has been eliminated by setting the discard_overlaps to 1 , since the duplicates created were found to be due to many overlapping hits.

    my $e2g = new Bio::EnsEMBL::Analysis::Runnable::Finished::MiniEst2Genome(
        '-genomic'    => $self->unmasked,
        '-features'   => $linear,
        '-seqcache' => $self->seq_cache,
	'-analysis'   => $self->analysis,
    );

    $e2g->run;

    foreach my $fp ( @{$e2g->output} ) {
        print STDERR "after:  " . $fp->gffstring . "\n" if $verbose;
        # source_tag and primary_tag have to be set to
        # something, or validate method in FeaturePair
        # (callled by RunnableDB) throws!
        #        $fp->feature2->source_tag('I_am_valid');
        #        $fp->feature2->primary_tag('I_am_valid');
        $self->add_output($fp);
    }
}


sub add_output {
    my ( $self, $f ) = @_;
    my $ana = $self->analysis;
    $f->analysis($ana);
    $self->{'output'} ||= [];
    push ( @{ $self->{'output'} }, $f );

}

sub _make_seqfetcher {
    my ( $self, $indices ) = @_;
    my $seqfetcher;
    if ( ref($indices) eq "ARRAY"){
        $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::Finished_xdget->new(-db => $indices,
        															 -executable => '/software/farm/bin/xdget');
    }
    else{
        $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch->new();
    }
    return $seqfetcher;
}

=head2 get_db_version

    Title   :  get_db_version
               [ distinguished from RunnableDB::*::db_version_searched() ]
    Useage  :  $self->get_db_version('version string')
               $obj->get_db_version()
    Function:  Get/Set a blast database version that was searched
               The actual look up is done in Runnable::Finished_Blast
               This is just a holding place for the string in this
               module
    Returns :  String or undef
    Args    :  String
    Caller  :  $self::run()
               RunnableDB::Finished_EST::db_version_searched()

=cut


sub get_db_version{
    my ($self, $arg) = @_;
    $self->{'_db_version_searched'} = $arg if $arg;
    return $self->{'_db_version_searched'};
}
1;


