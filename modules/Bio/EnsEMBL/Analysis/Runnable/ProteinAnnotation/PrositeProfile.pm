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

Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::PrositeProfile - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut


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


