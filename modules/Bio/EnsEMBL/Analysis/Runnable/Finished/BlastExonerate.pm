# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Runnable::Finished::BlastExonerate;
use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Analysis::Runnable::Finished::Blast;
use Bio::EnsEMBL::Analysis::Runnable::Finished::Exonerate;
use Bio::EnsEMBL::Analysis::Config::Blast;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use base 'Bio::EnsEMBL::Analysis::Runnable';


use vars qw($verbose);
$verbose = 1;

sub new {
    my ( $class, @args ) = @_;
    my $self = bless {}, $class;
    my ( $soft_masked, $hard_masked, $analysis, $seqfetcher ) = rearrange(
        [
        qw{
            SOFT_MASKED
            HARD_MASKED
            ANALYSIS
            SEQFETCHER
        }
        ],
        @args
    );
    throw "No SOFT_MASKED query (soft masked genomic sequence) given" unless $soft_masked;
    throw "No HARD_MASKED query (hard masked genomic sequence) given" unless $hard_masked;
    $self->soft_masked($soft_masked);
    $self->hard_masked($hard_masked);
    $self->analysis($analysis);

    my $indices = [ split(",", $self->analysis->db) ];
    $seqfetcher = $self->_make_seqfetcher($indices);
    $self->seqfetcher($seqfetcher);
    return $self;
}


sub hard_masked {
    my ( $self, $seq ) = @_;
    if ($seq) {
        unless ( $seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::Seq") ) {
            die Dumper($seq);
            throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
        }
        $self->{'_hard_masked'} = $seq;
    }
    return $self->{'_hard_masked'};
}

sub soft_masked {
    my ( $self, $seq ) = @_;
    if ($seq) {
        unless ( $seq->isa("Bio::PrimarySeqI") || $seq->isa("Bio::Seq") ) {
            die Dumper($seq);
            throw("Input isn't a Bio::Seq or Bio::PrimarySeq");
        }
        $self->{'_soft_masked'} = $seq;
    }
    return $self->{'_soft_masked'};
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
     -query => $self->hard_masked,
     -program => $blastcon->{-program} || 'wublastn',
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

    foreach my $db_file (sort @$db_files) {
    	print STDOUT "database file $db_file\n";
    	$blast->clean_output();
    	$blast->clean_results_files();
    	$blast->clean_databases();
    	$blast->databases($db_file);
	    $blast->run();
	    my $features =  $blast->output ;

		if (!$self->analysis->db)
		{
			print STDOUT "indice file $db_file\n";
			my $seqfetcher = $self->_make_seqfetcher([$db_file]);
	    	$self->seqfetcher($seqfetcher);
		}

		$self->run_exonerate($features,$blastcon) if @$features;

    }
}

sub run_exonerate {
	my ( $self, $feat, $exon_options ) = @_;

	# create a hash of features
	my $hit_features = {};
    for ( my $i = 0 ; $i < @$feat ; $i++ ) {

      my $f   = $feat->[$i];
      my $hid = $f->hseqname or throw("Missing hid");
      $hit_features->{$hid} ||= [];
      push ( @{ $hit_features->{$hid} }, $f );
    }

	# get a sorted list of unique feature IDs
    my @features_ids = sort {$a cmp $b } keys %$hit_features;

    print STDOUT "Blast output: ".scalar(@features_ids)." features\n" if $verbose;

	my $chunk_size = 200;
	for (my $i = 0; $i < @features_ids; $i += $chunk_size) {
        my $j = $i + $chunk_size - 1;
        # Set second index to last element if we're off the end of the array
        $j = $#features_ids if $#features_ids < $j;
        # Take a slice from the array
        my $chunk = [@features_ids[$i..$j]];
		print STDOUT "Process chunk of ".scalar(@$chunk)." features\n" if $verbose;
        # fetch the sequences
    	my $seqs = [];
		foreach my $id (@$chunk) {
			push @$seqs, $self->get_Sequence($id);
		}

		# run exonerate
		my $exonerate =
		Bio::EnsEMBL::Analysis::Runnable::Finished::Exonerate->new(
																	-analysis => $self->analysis,
																	-program  => $self->analysis->program,
																	-target		=> $self->soft_masked ,
																	-query_seqs	=> $seqs,
																	-exo_options => $exon_options->{-exo_options},
																	-query_type => $exon_options->{-query_type},
																	);
		$exonerate->run();
		my $exon_output = $exonerate->output();
		print STDOUT "Exonerate output: ".scalar(@$exon_output)." features\n" if $verbose;
		$self->output($exon_output);
    }
}

sub get_Sequence {
    my ($self,$id) = @_;
    my $seqfetcher = $self->seqfetcher;

    my $seq;
    eval {
      $seq = $seqfetcher->get_Seq_by_acc($id);
    };
    if ((!defined($seq)) && $@) {
      throw("Couldn't find sequence for [$id]:\n $@");
    }

    return $seq;
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
        $seqfetcher = Bio::EnsEMBL::Pipeline::SeqFetcher::xdget->new(-db => $indices,
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
               RunnableDB::BlastExonerate::db_version_searched()

=cut


sub get_db_version{
    my ($self, $arg) = @_;
    $self->{'_db_version_searched'} = $arg if $arg;
    return $self->{'_db_version_searched'};
}
1;


