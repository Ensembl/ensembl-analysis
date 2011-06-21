
=head1 NAME

PrositeProfile

=head1 SYNOPSIS


=head1 DESCRIPTION

=cut


package Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::PrositeProfile;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object


use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation);


sub multiprotein{
  my ($self) = @_;
  return 0;
}


sub run_analysis {
  my ($self) = @_;
  
  throw("Failed during Profile run $!\n") unless 
    (system ($self->program . ' -f ' . $self->queryfile. ' ' .
             $self->database . ' > ' .$self->resultsfile) == 0) ;
 
}


sub parse_results {
  my ($self,$seqid) = @_;
  
  my ($fh);
  my $resfile = $self->resultsfile;
  
  if (-e $resfile) {	
    if (-z $resfile) {  
      return; 
    } else {
      open ($fh, "<$resfile") or throw("Error opening ", $resfile,);
    }
  }
  
  my (@pfs);
  while (<$fh>) {
    if (/^\s*(\S+)\s+(\d+)\s*pos\.\s+(\d+)\s+\-\s+(\d+)\s+(\w+)\|/) {
      my ($sc, $rsc, $st, $en, $acc) = ($1, $2, $3, $4, $5);
      my $fp = $self->create_protein_feature($st,
                                             $en,
                                             $sc,
                                             $seqid,
                                             0, 0,
                                             $acc,
                                             $self->analysis,
                                             0, 0);
      push @pfs, $fp;
    }
  }
  
  $self->output(\@pfs);  
}


