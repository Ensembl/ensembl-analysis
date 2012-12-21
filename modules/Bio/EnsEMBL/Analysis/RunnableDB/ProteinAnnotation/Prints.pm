
=head1 NAME

Prints.pm - DESCRIPTION of Object

=head1 SYNOPSIS

 my $rsb = Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::Prints->new(
    -db       => $db
    -input_id    => $id
    -analysis    => $analysis);


=head1 DESCRIPTION


=cut


# Let the code begin...


# $Source: /tmp/ENSCOPY-ENSEMBL-ANALYSIS/modules/Bio/EnsEMBL/Analysis/RunnableDB/ProteinAnnotation/Prints.pm,v $
# $Version: $
package Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::Prints;
use warnings ;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation;
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Prints;


@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation);


sub fetch_input {
  my ($self, @args) = @_;

  $self->SUPER::fetch_input(@args);

  my $run =  Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Prints->new(-query     => $self->query,
                                                                              -analysis  => $self->analysis);
  $self->runnable($run);
}




1;










