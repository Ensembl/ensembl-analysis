### Bio::EnsEMBL::Analysis::Runnable::Finished::EST


package Bio::EnsEMBL::Analysis::Runnable::Finished::EST;

use strict;
use warnings;

#use vars qw(@ISA);

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Analysis::Runnable::Finished::MiniEst2Genome;
use Bio::EnsEMBL::Analysis::Runnable::Finished::Blast;
use Bio::EnsEMBL::Analysis::Config::Blast;
use Bio::EnsEMBL::Analysis::Tools::BPliteWrapper;
use base 'Bio::EnsEMBL::Analysis::Runnable';


#use vars qw(@ISA $verbose);
#@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

use vars qw($verbose);
$verbose = 0;

sub new {
    my ( $class, @args ) = @_;
    my $self = bless {}, $class;
    my ( $query, $unmasked, $analysis, $seqfetcher ) = $self->_rearrange(
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

    # this should make a OBDAIndexSeqFetcher, but will default to Pfetch
    my $indices = [ split(",", $self->analysis->db) ];
    $seqfetcher = $self->_make_seqfetcher($indices);
    $self->seqfetcher($seqfetcher);
    return $self;
}

sub query {
    my ( $self, $seq ) = @_;
    if ($seq) {
        unless ( $seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::Seq") ) {
            $self->throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
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
            $self->throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
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
    $blast->run();
    # keep the features
    my $features =  $blast->output ;

    # set the db_version_searched
    my $dbv = $blast->get_db_version();
    print "using $dbv\n";
    $self->get_db_version($dbv);
    print "\nPlus strand est_genome\n" if $verbose;
    $self->run_est_genome_on_strand( 1, $features );

    print  "\nMinus strand est_genome\n" if $verbose;
    $self->run_est_genome_on_strand( -1, $features)
}

sub run_est_genome_on_strand {

    my ( $self, $strand, $feat ) = @_;
    my $hit_features = {};
    for ( my $i = 0 ; $i < @$feat ; $i++ ) {

      my $f   = $feat->[$i];
      my $hid = $f->hseqname or $self->throw("Missing hid");
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
        '-seqfetcher' => $self->seqfetcher,
	'-analysis'   => $self->analysis,
    );

    $e2g->run;

    foreach my $fp ( $e2g->output ) {
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
    $self->{'_output'} ||= [];
    push ( @{ $self->{'_output'} }, $f );

}

sub output {
    my ($self) = @_;
    if ( my $out = $self->{'_output'} ) {
      return @$out;
      #  return $out;
    }
    else {
        return;
    }
}

sub _make_seqfetcher {
    my ( $self, $indices ) = @_;
    my $seqfetcher;
    if ( ref($indices) eq "ARRAY"){
        $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher->new(-db => $indices);
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

__END__

=head1 NAME - Bio::EnsEMBL::Analysis::Runnable::Finished_EST

=head1 AUTHOR

James Gilbert B<email> jgrg@sanger.ac.uk

Modified by Sindhu K. Pillai B<email> sp1@sanger.ac.uk
