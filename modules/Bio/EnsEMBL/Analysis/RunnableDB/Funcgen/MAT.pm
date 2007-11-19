# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::MAT
#
# Copyright (c) 2007 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::MAT

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::MAT->new
     (
         -db       => $db,
         -input_id => 'chromosome::20:1:100000:1',
         -analysis => $analysis,
     );
  $runnable->fetch_input;
  $runnable->run;
  $runnable->write_output;

=head1 DESCRIPTION

This module provides an interface between the ensembl functional genomics 
database and the Runnable MAT which wraps the program MAT (Model-based 
Analysis of Tiling-array, see http://chip.dfci.harvard.edu/~wli/MAT).

=head1 AUTHOR

Stefan Graf, Ensembl Functional Genomics - http://www.ensembl.org/

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::MAT;

use strict;
use warnings;
use Data::Dumper;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::RunnableDB::Funcgen;
use Bio::EnsEMBL::Analysis::Runnable::Funcgen::MAT;

use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Analysis::Config::Funcgen::MAT;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use vars qw(@ISA); 

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::Funcgen);

=head2 new

  Arg [1]     : 
  Arg [2]     : 
  Description : Instantiates new MAT runnabledb
  Returntype  : Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::MAT object
  Exceptions  : 
  Example     : 

=cut

sub new {

    print "Analysis::RunnableDB::Funcgen::MAT::new\n";

    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    $self->Bio::EnsEMBL::Analysis::RunnableDB::read_and_check_config($CONFIG);

	### need to do some variable/config checking here ###

	#print Dumper $self;

    return $self;
	
}

sub fetch_input {

    my ($self) = @_;
    print "Analysis::RunnableDB::Funcgen::MAT::fetch_input\n";

    my $runnable = 'Bio::EnsEMBL::Analysis::Runnable::Funcgen::'
        .$self->analysis->module;

    $runnable = $runnable->new
        (
         -query => $self->query,
         -program => $self->analysis->program_file,
         -analysis => $self->analysis,
         );
    
    $self->runnable($runnable);

    return 1;

}

1;
