
=pod 

=head1 NAME

  Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::Coil

=head1 SYNOPSIS

  my $ncoils = Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::Coil->new ( -db      => $db,
    	  	                                                        -input_id   => $input_id,
                                                                        -analysis   => $analysis,
                                                                      );
  $ncoils->fetch_input;  # gets sequence from DB
  $ncoils->run;
  $ncoils->write_output; # writes features to to DB

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::Coil;

use warnings ;
use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation;
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Coil;

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation);


sub fetch_input {
  my ($self, @args) = @_;

  $self->SUPER::fetch_input(@args);

  my $run = Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Coil->new(-query     => $self->query,
                                                                           -analysis  => $self->analysis);
  
  $self->runnable($run);
}


1;
