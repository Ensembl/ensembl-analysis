# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::Nessie
#
# Copyright (c) 2007 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::Nessie

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::Nessie->new
     (
         -db       => $db,
         -input_id => 'chromosome::20:1:100000:1',
         -analysis => $analysis,
     );
  $runnable->fetch_input;
  $runnable->run;
  $runnable->write_output;

=head1 DESCRIPTION

This module provides an interface between the ensembl database and
the Runnable Nessie which wraps the program Nessie

=head1 AUTHOR

Stefan Graf, Ensembl Functional Genomics - http://www.ensembl.org/

=head1 CONTACT

Post questions to the Ensembl development list: dev@ensembl.org

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::Nessie;

use strict;
use warnings;
use Data::Dumper;

use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Analysis::Config::Funcgen::Nessie;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::RunnableDB::Funcgen;
use Bio::EnsEMBL::Analysis::Runnable::Funcgen::Nessie;

use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);
use vars qw(@ISA); 

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::Funcgen);

=head2 new

  Arg [1]     : 
  Arg [2]     : 
  Description : Instantiates new Nessie runnabledb
  Returntype  : Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::Nessie object
  Exceptions  : 
  Example     : 

=cut

sub new {

    print "Analysis::RunnableDB::Funcgen::Nessie::new\n";
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);

    $self->read_and_check_config($CONFIG);

    # add some runnable/program special params to analysis
    $self->PARAMETERS(join('; ',
                           $self->PARAMETERS.$self->TRAIN_PARAMETERS,
                           $self->PEAK_PARAMETERS));
    #print Dumper $self->PARAMETERS;
    
    # make sure we have the correct analysis object
    $self->check_Analysis();

    # make sure we can store the correct feature_set, data_sets, and result_sets
    $self->check_Sets();

    return $self;

}

sub TRAIN_PARAMETERS {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{'_CONFIG_TRAIN_PARAMETERS'} = $value;
    }

    if ( exists( $self->{'_CONFIG_TRAIN_PARAMETERS'} ) ) {
        return $self->{'_CONFIG_TRAIN_PARAMETERS'};
    } else {
        return undef;
    }
}

sub PEAK_PARAMETERS {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{'_CONFIG_PEAK_PARAMETERS'} = $value;
    }

    if ( exists( $self->{'_CONFIG_PEAK_PARAMETERS'} ) ) {
        return $self->{'_CONFIG_PEAK_PARAMETERS'};
    } else {
        return undef;
    }
}




1;
