=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
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

Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::Seg - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut


=pod 

=head1 NAME

  Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::Seg

=head1 SYNOPSIS

  my $seg = Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::Seg->new ( -db      => $db,
	    	                                                    -input_id   => $input_id,
                                                                    -analysis   => $analysis,
                                                                  );
  $seg->fetch_input;  # gets sequence from DB
  $seg->run;
  $seg->output;
  $seg->write_output; # writes features to to DB

 NB: The input_id can either be a peptide id or the location for a protein file. 

=head1 DESCRIPTION

  This object wraps Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Seg
  to add functionality to read and write to databases.
  The query sequence is provided through the input_id.
  The appropriate Bio::EnsEMBL::Analysis object
  must be passed for extraction of parameters.

=head1 CONTACT

  Marc Sohrmann: ms2@sanger.ac.uk

=head1 APPENDIX

  The rest of the documentation details each of the object methods. 
  Internal methods are usually preceded with a _.

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::Seg;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation;
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Seg;

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation);

sub fetch_input {
  my ($self, @args) = @_;

  $self->SUPER::fetch_input(@args);

  my $run = Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Seg->new(-query     => $self->query,
                                                                          -analysis  => $self->analysis);
  $self->runnable($run);
}


1;
