=pod 

=head1 NAME

  Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::Signalp

=head1 SYNOPSIS

  my $signalp = Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::Signalp->new ( -db      => $db,
    	  	                                                            -input_id   => $input_id,
                                                                            -analysis   => $analysis,
                                                                          );
  $signalp->fetch_input;  # gets sequence from DB
  $signalp->run;
  $signalp->write_output; # writes features to to DB

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::Signalp;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation;
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Signalp;

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation);

sub fetch_input {
  my ($self, @args) = @_;

  $self->SUPER::fetch_input(@args);
  
  my $run = Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Signalp->new(-query     => $self->query,
                                                                              -analysis  => $self->analysis);
  $self->runnable($run);
}


1;
