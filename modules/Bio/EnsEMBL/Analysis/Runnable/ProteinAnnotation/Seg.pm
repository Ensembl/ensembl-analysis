# Author: Marc Sohrmann (ms2@sanger.ac.uk)
# Copyright (c) Marc Sohrmann, 2001
# You may distribute this code under the same terms as perl itself
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod

=head1 NAME

  Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Seg

=head1 SYNOPSIS

  my $seqstream = Bio::SeqIO->new ( -file => $queryfile,
                                    -fmt => 'Fasta',
                                  );
  $seq = $seqstream->next_seq;

  my $seg = Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Seg->new ( -QUERY => $seq);
  $seg->workdir ($workdir);
  $seg->run;
  my @results = $seg->output;

=head1 DESCRIPTION

  Seg takes a Bio::Seq (or Bio::PrimarySeq) object
  and runs seg on it (detecting low complexity sequences).
  The resulting output file is parsed to produce a set of features.

=head1 CONTACT

  http://www.ensembl.org/Help/Contact

=head1 APPENDIX

  The rest of the documentation details each of the object methods.
  Internal methods are usually preceded with a _.

=cut

package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Seg;

use vars qw(@ISA);
use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation);

sub run_analysis {
  my ($self) = @_;

  throw( "Error running " . $self->program . " on " . $self->queryfile )
    unless ( ( system( $self->program . " " . $self->queryfile . " -l > " . $self->resultsfile ) ) == 0 );
}

sub parse_results {
  my ($self) = @_;

  my ($fh);

  my $resfile = $self->resultsfile;

  if ( -e $resfile ) {
    # it's a filename
    if ( -z $resfile ) {
      return;
    }
    else {
      open( $fh, "<$resfile" ) or throw("Error opening $resfile");
    }
  }
  else {
    # it'a a filehandle
    $fh = $resfile;
  }

  # parse
  my @pfs;
  while (<$fh>) {
    chomp;
    next if /^$/;
    if (/^\>/) {
      /^\>(\S+)?\((\d+)\-(\d+)\)\s*complexity=(\S+)/;
      my $tid   = $1;
      my $start = $2;
      my $end   = $3;
      my $score = $4;

      my $fp = $self->create_protein_feature( $start, $end, $score, $tid, 0, 0, 'Seg', $self->analysis, 0, 0 );
      push @pfs, $fp;
    }
  }
  close($fh);

  $self->output( \@pfs );
} ## end sub parse_results

=head2 get_low_complexity_length

 Title    : get_low_complexity_length
 Usage    : $len = $self->get_low_complexity_length;
 Function : returns *percentage* low complexity of protein
 Example  :
 Returns  : a percentage_id
 Args     :
 Throws   :
 Notes    : It only makes sense to call this method when the
    Runnable was created with a single Bio::Seq

=cut

sub get_low_complexity_length {
  my ($self) = @_;

  if ( $self->query->length > 0 ) {
    my $lc_length = 0;

    foreach my $feat ( @{ $self->output } ) {
      $lc_length += abs( $feat->end - $feat->start ) + 1;
    }

    my $low_complexity = ($lc_length)/( $self->query->length );

    $low_complexity *= 100;

    return $low_complexity;
  }
  else {
    return 0;
  }
}

1;
