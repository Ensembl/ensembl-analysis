=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Signalp - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut


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
      unless ((system($self->program . " -t euk -trunc 200 "
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
	# signalp 2 output, expecting to parse the next line
	my $line = <$fh>;
        if ($line =~ /Most likely cleavage site between pos\.\s+(\d+)/) {
	    $end = $1;
        }
        else {
          # otherwise signalp3 output, expecting to parse the line after
	  my $line = <$fh>;
	  if ($line =~ /Most likely cleavage site between pos\.\s+(\d+)/) {
	      $end = $1;
	  }
	  else {
	      # otherwise throw an exception
	      $self->throw ("parsing problem in ".$self->program);
	  }
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
