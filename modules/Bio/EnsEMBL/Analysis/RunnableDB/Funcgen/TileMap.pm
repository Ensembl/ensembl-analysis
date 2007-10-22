# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::TileMap
#
# Copyright (c) 2007 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::TileMap

=head1 SYNOPSIS

  my $chipotle = Bio::EnsEMBL::Analysis::RunnableDB::Funcgen::TileMap->new
     (
        -input_id => 'chromosome::20:1:100000:1',
        
     );
  $chipotle->fetch_input;
  $chipotle->run;
  $chipotle->write_output;

=head1 DESCRIPTION

This module provides an interface between the ensembl database and
the Runnable TileMap which wraps the program ChIPoTle

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

    #print "TileMap::new\n";
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);

    $self->read_and_check_config($TILEMAP_CONFIG);

    return $self;
	
}

sub fetch_ResultSets
{
    my $self = shift;

    my $rsa = $self->db->get_ResultSetAdaptor();
    my $rsets = $rsa->fetch_all_by_Experiment_Analysis
        ($self->experiment, $self->norm_analysis);
    print Dumper ($self->experiment->name, $self->analysis->logic_name);
    print "No. of available ResultSets: ", scalar(@$rsets), "\n";

    my @rsets = ();
    my $regex = $self->RESULT_SET_REGEXP;
    foreach my $rset (@{$rsets}) {
        #print Dumper $rset->name();
        next if ($rset->name() !~ m/$regex$/);
        push(@rsets, $rset);
    }

    return \@rsets;;

}

1;
