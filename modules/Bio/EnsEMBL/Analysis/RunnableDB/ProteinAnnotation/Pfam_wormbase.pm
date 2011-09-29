#
#
#
=pod 

=head1 NAME

  Bio::EnsEMBL::Pipeline::RunnableDB::ProteinAnnotation::Pfam_wormbase

=head1 SYNOPSIS

  my $seg = Bio::EnsEMBL::Pipeline::RunnableDB::ProteinAnnotation::Pfam_wormbase->
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
  a Pfam-specific way.

=head1 CONTACT

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::Pfam_wormbase;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation;
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Pfam_wormbase;

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation);


sub fetch_input {
  my ($self) = @_;
 
  $self->SUPER::fetch_input;
 
  my $run = Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Pfam_wormbase->
      new(-query     => $self->query,
          -analysis  => $self->analysis,
          -database  => $self->analysis->db_file,
          %{$self->parameters_hash}
          );
  $self->runnable($run);
}



