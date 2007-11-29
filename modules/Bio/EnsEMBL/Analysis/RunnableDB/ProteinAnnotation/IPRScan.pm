#
#
#
=pod 

=head1 NAME

  Bio::EnsEMBL::Pipeline::RunnableDB::ProteinAnnotation::IPRScan

=head1 SYNOPSIS

  my $seg = Bio::EnsEMBL::Pipeline::RunnableDB::ProteinAnnotation::IPRScan->
    new ( -db      => $db,
    -input_id   => $input_id,
    -analysis   => $analysis,
                                                                      );
  $seg->fetch_input;  # gets sequence from DB
  $seg->run;
  $seg->write_output; # writes features to to DB

=head1 DESCRIPTION

  This object wraps Bio::EnsEMBL::Pipeline::Runnable::Hmmpfam
  to add functionality to read and write to databases in 
  a IPRScan-specific way.

=head1 CONTACT

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::IPRScan;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation;
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::IPRScan;

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation);


sub fetch_input {
  my ($self) = @_;
 
  $self->SUPER::fetch_input;
  print "FETCHING INPUT\n";
  my $run = Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::IPRScan->
      new(
          -query     => $self->query,
          -analysis  => $self->analysis,
          -program => $self->analysis->program_file,
          %{$self->parameters_hash}
         );
  $self->runnable($run);
}
