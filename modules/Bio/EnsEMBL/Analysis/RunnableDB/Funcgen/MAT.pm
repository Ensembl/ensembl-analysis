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

    $self->read_and_check_config($CONFIG);

    return $self;
	
}

=head2 query

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : integer
  Function  : container for chip number 
  Returntype: integer
  Exceptions: throws if passed the incorrect value
  Example   : 

=cut


sub query{
  my ($self, $chip) = @_;
  if($chip){
    throw("Must pass RunnableDB::Funcgen::MAT::query an integer ".
          "specifying the chip to process not ".$chip) 
        unless($chip =~ m/^\d+$/);
    $self->{'chip'} = $chip;
  }
  return $self->{'chip'};
}

=head2 check_Sets

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Function  : fetch and set ResultSets of interest
  Returntype: 1
  Exceptions: none
  Example   : 

=cut


sub check_Sets
{

    warn("\nNEED TO IMPLEMENT SETTING OF RESULT/DATA/FEATURE SET!!!\n\n");

}

sub fetch_input {

    my ($self) = @_;
    print "Analysis::RunnableDB::Funcgen::MAT::fetch_input\n";

    my $runnable = 'Bio::EnsEMBL::Analysis::Runnable::Funcgen::'
        .$self->analysis->module;

    warn("No input to fetch since we are going to work directly on the CEL files.");

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
