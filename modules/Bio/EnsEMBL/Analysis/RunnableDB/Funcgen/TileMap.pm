# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::TileMap
#
# Copyright (c) 2007 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::TileMap

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::TileMap->new
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
the Runnable TileMap which wraps the program TileMap

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Stefan Graf, Ensembl Functional Genomics - http://www.ensembl.org

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::TileMap;

use strict;
use warnings;
use Data::Dumper;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::RunnableDB::Funcgen;
use Bio::EnsEMBL::Analysis::Runnable::Funcgen::TileMap;

use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Analysis::Config::Funcgen::TileMap;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use vars qw(@ISA); 

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::Funcgen);



=head2 new

  Arg [1]     : 
  Arg [2]     : 
  Description : Instantiates new TileMap runnabledb
  Returntype  : Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::TileMap object
  Exceptions  : 
  Example     : 

=cut

sub new {

    print "Analysis::RunnableDB::Funcgen::TileMap::new\n";
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    $self->read_and_check_config($CONFIG);

    return $self;
	
}

1;
