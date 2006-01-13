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
  
  Marc Sohrmann: ms2@sanger.ac.uk

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


sub multiprotein{
  my ($self) = @_;
  return 1;
}


sub run_analysis {
  my ($self) = @_;

  throw ("Error running ".$self->program." on ".$self->queryfile) 
      unless ((system ($self->program." ".$self->queryfile." -l > ".
                       $self->resultsfile)) == 0); 
}


sub parse_results {
  my ($self) = @_;

  my ($fh);

  my $resfile = $self->resultsfile;
  
  if (-e $resfile) {
    # it's a filename
    if (-z $resfile) {        
      return;
    }else {
      open($fh, "<$resfile") or throw ("Error opening $resfile");
    }
  } else {
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
      my $tid = $1;
      my $start = $2;
      my $end = $3;
      my $score = $4;
      
      my $fp = $self->create_protein_feature($start, $end, $score, $tid, 
                                             0, 0, 'Seg', 
                                             $self->analysis, 0, 0);
      push @pfs, $fp;
    }
  }
  close($fh);

  $self->output(\@pfs);
}

1;
