# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::Chipotle
#
# Copyright (c) 2007 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::Chipotle

=head1 SYNOPSIS

  my $chipotle = Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::Chipotle->new
     (
        -input_id => 'chromosome::20:1:100000:1',
        
     );
  $chipotle->fetch_input;
  $chipotle->run;
  $chipotle->write_output;

=head1 DESCRIPTION

This module provides an interface between the ensembl database and
the Runnable Chipotle which wraps the program ChIPoTle

This module can fetch appropriate input from the database
pass it to the runnable then write the results back to the database
in the predicted_feature table with respect to the resut_set, data_set 
and feature_set table

=head1 AUTHOR

This module was created by Stefan Graf. It is part of the 
Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::Chipotle;

use strict;
use warnings;
use Data::Dumper;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::RunnableDB::Funcgen;
use Bio::EnsEMBL::Analysis::Runnable::Funcgen::Chipotle;

use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Analysis::Config::Funcgen::Chipotle;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use vars qw(@ISA); 

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::Funcgen);



=head2 new

  Arg [1]     : 
  Arg [2]     : 
  Description : Instantiates new Chipotle runnabledb
  Returntype  : Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::Chipotle object
  Exceptions  : 
  Example     : 

=cut

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  
  $self->read_and_check_config($CHIPOTLE_CONFIG);

  # Chipotle analysis defaults
  $self->analysis->description('Sliding window approach bases on Gaussian '.
                               'background model method for significance '.
                               'estimation.');
  $self->analysis->display_label('ChIPOTle');

  #print Dumper $self;

  return $self;
}

################################################################################
# Declare and set up config variables
################################################################################

sub read_and_check_config {
  my $self = shift;

  $self->SUPER::read_and_check_config($CHIPOTLE_CONFIG);

  ##########
  # CHECKS
  ##########

  # check that compulsory options have values
  foreach my $config_var 
      (
       qw(
          PROGRAM
          OPTIONS
          LOGIC_NAME
          EFG_EXPERIMENT
          )
       ){
          if ( not defined $self->$config_var ){
              throw("You must define $config_var in config.");
          }
          
          #print Dumper $self->$config_var;
          
      }
  
  # make sure EFG_EXPERIMENT exists, is a hash, and contains all 
  # compulsory options
  throw("EFG_EXPERIMENT for ".$self->experiment." is not defined.")
      if(! exists $self->EFG_EXPERIMENT->{$self->experiment});
  throw("EFG_EXPERIMENT must be a hash ref not ".$self->EFG_EXPERIMENT.
        " Chipotle::read_and_check_config")
      if(ref($self->EFG_EXPERIMENT) ne 'HASH');
  foreach my $config_var
      (
       qw(
          FT_NAME
          FT_CLASS
          FT_DESC
          CT_NAME
          CT_DESC
          )
       ){
          throw("Must define $config_var in EFG_EXPERIMENT config.")
          if (! defined $self->EFG_EXPERIMENT->{$self->experiment}->{$config_var});
      }

}

#############################################################
###     end of config
#############################################################

1;
