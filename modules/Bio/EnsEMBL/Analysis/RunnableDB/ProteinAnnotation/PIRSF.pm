#
#
#
=pod 

=head1 NAME

  Bio::EnsEMBL::Pipeline::RunnableDB::ProteinAnnotation::PIRSF

=head1 SYNOPSIS

  my $seg = Bio::EnsEMBL::Pipeline::RunnableDB::ProteinAnnotation::PIRSF->
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

# $Source: /tmp/ENSCOPY-ENSEMBL-ANALYSIS/modules/Bio/EnsEMBL/Analysis/RunnableDB/ProteinAnnotation/PIRSF.pm,v $
# $Version: $
package Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::PIRSF;

use warnings ;
use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation;
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::PIRSF;

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation);


sub fetch_input {
  my ($self) = @_;
 
  $self->SUPER::fetch_input;
 
  my $run = Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::PIRSF->
        new(-query     => $self->query,
            -analysis  => $self->analysis,
            -database  => $self->analysis->db_file,
            %{$self->parameters_hash}
            );
  $self->runnable($run);
}



