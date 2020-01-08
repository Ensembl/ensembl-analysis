# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2020] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

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

    $self->PARAMETERS(join('; ',
                           'METHOD='.$self->METHOD,
                           'POSTPROB='.$self->POSTPROB,
                           'MAXGAP='.$self->MAXGAP,
                           'HYBLENGTH='.$self->HYBLENGTH,
                           'TEMPLATE_FILE='.$self->TEMPLATE_FILE));

    # make sure we have the correct analysis object
    $self->check_Analysis();

    # make sure we can store the correct feature_set, data_sets, and result_sets
    $self->check_Sets();

    return $self;
	
}

sub TEMPLATE_FILE {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{'_CONFIG_TEMPLATE_FILE'} = $value;
    }

    if ( exists( $self->{'_CONFIG_TEMPLATE_FILE'} ) ) {
        return $self->{'_CONFIG_TEMPLATE_FILE'};
    } else {
        return undef;
    }
}

sub POSTPROB {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{'_CONFIG_POSTPROB'} = $value;
    }

    if ( exists( $self->{'_CONFIG_POSTPROB'} ) ) {
        return $self->{'_CONFIG_POSTPROB'};
    } else {
        return undef;
    }
}

sub HYBLENGTH {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{'_CONFIG_HYBLENGTH'} = $value;
    }

    if ( exists( $self->{'_CONFIG_HYBLENGTH'} ) ) {
        return $self->{'_CONFIG_HYBLENGTH'};
    } else {
        return undef;
    }
}

sub MAXGAP {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{'_CONFIG_MAXGAP'} = $value;
    }

    if ( exists( $self->{'_CONFIG_MAXGAP'} ) ) {
        return $self->{'_CONFIG_MAXGAP'};
    } else {
        return undef;
    }
}

sub METHOD {
    my ( $self, $value ) = @_;

    if ( defined $value ) {
        $self->{'_CONFIG_METHOD'} = $value;
    }

    if ( exists( $self->{'_CONFIG_METHOD'} ) ) {
        return $self->{'_CONFIG_METHOD'};
    } else {
        return undef;
    }
}

1;
