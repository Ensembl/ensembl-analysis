#
#
#
=pod 

=head1 NAME

  Bio::EnsEMBL::Pipeline::RunnableDB::ProteinAnnotation::Superfamily_wormbase

=head1 SYNOPSIS

  my $seg = Bio::EnsEMBL::Pipeline::RunnableDB::ProteinAnnotation::Superfamily_wormase
    ->new ( -db      => $db,
    -input_id   => $input_id,
    -analysis   => $analysis,
    );
  $seg->fetch_input;  # gets sequence from DB
  $seg->run;
  $seg->output;
  $seg->write_output; # writes features to to DB

=head1 DESCRIPTION

  This object wraps Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Superfamily_wormbase
  to add functionality to read and write to databases in 


=head1 CONTACT

=cut

# $Source: /tmp/ENSCOPY-ENSEMBL-ANALYSIS/modules/Bio/EnsEMBL/Analysis/RunnableDB/ProteinAnnotation/Superfamily_wormbase.pm,v $
# $Revision: 1.3 $
package Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::Superfamily_wormbase;

use warnings ;
use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation;
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Superfamily_wormbase;

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation);


#
# overridden methods
#
sub fetch_input {
  my ($self) = @_;

  $self->SUPER::fetch_input;

  my $run = Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Superfamily_wormbase->
      new(-query     => $self->query,
          -analysis  => $self->analysis,
          -database  => $self->analysis->db_file,
          %{$self->parameters_hash}
          );
  $self->runnable($run);    
}


1;
