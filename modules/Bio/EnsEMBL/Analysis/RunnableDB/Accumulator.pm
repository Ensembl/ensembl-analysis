# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::Accumulator
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Accumulator

=head1 SYNOPSIS

  my $accumulator = Bio::EnsEMBL::Analysis::RunnableDB::Accumulator->
  new(
      -input_id => 'ACCUMULATOR',
      -db => $db,
      -analysis => $analysis,
     );
  $accumulator->fetch_input;
  $accumulator->run;
  $accumulator->write_output;

=head1 DESCRIPTION

This is a simple place holder module to allow the accumulator wait for all
stages in the pipeline to work. It does nothing just

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut


package Bio::EnsEMBL::Analysis::RunnableDB::Accumulator;

use strict;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB);

=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Dummy method to comply to the interface
    Returns :   none
    Args    :   none

=cut

sub fetch_input {
    my( $self) = @_;
    
    throw("No input id") unless defined($self->input_id);

    return 1;

}

sub run {
    my ($self) = @_;
    print "Dummy RunnableDB - no runnable to run\n";

}

sub write_output {
    my ($self) = @_;

    print "Dummy RunnableDB - no output to write\n";

    return 1;
}

1;
