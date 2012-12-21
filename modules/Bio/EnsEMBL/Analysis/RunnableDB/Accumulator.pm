=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
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

Bio::EnsEMBL::Analysis::RunnableDB::Accumulator - 

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

=head1 METHODS

=cut


# $Source: /tmp/ENSCOPY-ENSEMBL-ANALYSIS/modules/Bio/EnsEMBL/Analysis/RunnableDB/Accumulator.pm,v $
# $Version: $
package Bio::EnsEMBL::Analysis::RunnableDB::Accumulator;

use warnings ;
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
