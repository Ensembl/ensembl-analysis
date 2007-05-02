# Ensembl module for Bio::EnsEMBL::Analysis::Runnable::Funcgen
#
# Copyright (c) 2007 Ensembl
#
=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Fungen

=head1 SYNOPSIS

=head1 DESCRIPTION

This module is the base class for the Fungen Runnables.

=head1 AUTHOR

This module was created by Stefan Graf. It is part 
of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::Runnable::Funcgen;

use strict;
use warnings;
use Data::Dumper;

use Bio::EnsEMBL::Analysis::Runnable;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

=head2 probe_features

  Arg [1]     : Bio::EnsEMBL::Analysis::RunnableDB::Chipotle
  Arg [2]     : arrayref of probe features
  Description : container for probe features
  Returntype  : arrayref
  Exceptions  : throws if no probe feature container is defined
  Example     : 

=cut

sub probe_features {

  my ($self, $features) = @_;

  if($features){
      $self->{'probe_features'} = $features;
  }

  throw("No probe features available in Runnable.") 
      if (!$self->{'probe_features'});

  return $self->{'probe_features'};

}

1;
