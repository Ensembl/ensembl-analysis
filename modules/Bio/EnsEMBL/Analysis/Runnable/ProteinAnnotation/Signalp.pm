
=pod 

=head1 NAME

 Bio::EnsEMBL::Analysis::Runnable::Protein::Signalp

=head1 SYNOPSIS

 my $signalp = Bio::EnsEMBL::Analysis::Runnable::Protein::Signalp->new(
    -query => $seq,
    -analysis => $ana,
 );
 $signalp->run;
 my @results = @{$signalp->output};

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Signalp;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation;
use Bio::SeqIO;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation);



sub multiprotein{
  my ($self) = @_;
  return 0;
}


#
# we over-ride write_seq_file because Signalp works best
# if only the first 50 a.a. of the protein is supplied 

sub write_seq_file{
  my ($self, $seq, $filename) = @_;
  if(!$seq){
    $seq = $self->query;
  }
  if(!$filename){
    $filename = $self->queryfile;
  }
  my $seqout = Bio::SeqIO->new(
                               -file => ">".$filename,
                               -format => 'Fasta',
                              );
  eval{
    my $init_seq = substr($seq->seq, 0, 50);
    my $newseqobj = Bio::Seq->new(-id => $seq->display_id,
                                  -seq => $init_seq);
    $seqout->write_seq($newseqobj);
  };
  if($@){
    throw("FAILED to write $seq to $filename Runnable:write_seq_file");
  }
  return $filename;
}


sub run_analysis {
  my ($self) = @_;
  
  throw("Error running ".$self->program." on ".$self->queryfile) 
      unless ((system($self->program . " -t euk "
                      .$self->queryfile 
                      . " > ".$self->resultsfile)) == 0); 
}


sub parse_results {
  my ($self) = @_;

  my ($fh);

  my $resfile = $self->resultsfile;
  
  if (-e $resfile) {
    # it's a filename
    if (-z $resfile) {  
      # no hits
      return;
    }else {
      open ($fh, "<$resfile") or throw ("Error opening $resfile");
    }
  }else {
    # it'a a filehandle
    $fh = $resfile;
  }
  
  # parse
  my (@fps, $id, $fact1, $fact2, $end);
  while (<$fh>) {
    chomp;
    if (/^\>(\S+)/) {
      $id = $1;
    }
    elsif (/max\.\s+Y\s+(\S+)\s+\S+\s+\S+\s+(\S+)/) {
      $fact1 = $2;
    }
    elsif (/mean\s+S\s+(\S+)\s+\S+\s+\S+\s+(\S+)/) {
      $fact2 = $2;
      if ($fact1 eq "YES" && $fact2 eq "YES") {
        my $line = <$fh>;
        if ($line =~ /Most likely cleavage site between pos\.\s+(\d+)/) {
          $end = $1;
        }
        else {
          $self->throw ("parsing problem in ".$self->program);
        }
        my $fp = $self->create_protein_feature(1, $end, 0, $id, 0, 0,
                                               'Sigp', $self->analysis,
                                               0, 0);
        push @fps, $fp;
      }
    }
  }
  close($fh);

  $self->output(\@fps);
}



1;
