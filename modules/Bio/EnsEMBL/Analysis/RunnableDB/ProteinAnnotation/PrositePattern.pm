
=pod 

=head1 NAME

  Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::PrositePattern

=head1 SYNOPSIS

  my $tmhmm = Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::PrositePattern->
    new ( 
    -db      => $db,
    -input_id   => $input_id,
    -analysis   => $analysis)
    );
  $tmhmm->fetch_input;  # gets sequence from DB
  $tmhmm->run;
  $tmhmm->write_output; # writes features to to DB

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::PrositePattern;

use warnings ;
use strict;
use vars qw(@ISA);


use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::PrositePattern;
use Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation;


@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation);


sub fetch_input {
  my ($self, @args) = @_;
  
  $self->SUPER::fetch_input(@args);

  my $run = Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::PrositePattern->
      new(-query     => $self->query,
          -analysis  => $self->analysis,
          %{$self->parameters_hash}
          );
  $self->runnable($run);
}

1;
