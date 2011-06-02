
=head1 NAME

Prints.pm - DESCRIPTION of Object

=head1 SYNOPSIS

 my $rsb = Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::Prints_wormbase->new(
    -db       => $db
    -input_id    => $id
    -analysis    => $analysis);


=head1 DESCRIPTION


=cut


# Let the code begin...


package Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::Prints_wormbase;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation;
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Prints_wormbase;


@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation);


sub fetch_input {
  my ($self, @args) = @_;

  $self->SUPER::fetch_input(@args);

  my $run =  Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Prints_wormbase->new(-query     => $self->query,
                                                                              -analysis  => $self->analysis);
  $self->runnable($run);
}




1;










