# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::EPCR
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

=head1 SYNOPSIS

  my $runnabledb = Bio::EnsEMBL::Analysis::RunnableDB::EPCR->
  new(
      -input_id => 'contig::AL805961.22.1.166258:1:166258:1',
      -db => $db,
      -analysis => $analysis,
     );
  $runnabledb->fetch_input;
  $runnabledb->run;
  $runnabledb->write_output;


=head1 DESCRIPTION



=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::EPCR;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Runnable::EPCR;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB);


=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::EPCR
  Function  : fetch data out of database and create runnable
  Returntype: none
  Exceptions: throws if no markers in the database when no -STS_FILE
  is defined
  Example   : 

=cut

sub fetch_input{
  my ($self) = @_;
  my %parameters = %{$self->parameters_hash};
  if($self->analysis->db_file){
    $parameters{'-STS_FILE'} = $self->analysis->db_file 
      unless($parameters{'-STS_FILE'});
  }
  if(!$parameters{'-STS_FILE'}){
    my $sts = $self->db->get_MarkerAdaptor->fetch_all;
    throw("No markers in ".$self->db->dbname) unless(@$sts);
    $parameters{'-STS_FEATURES'} = $sts;
  }
  my $slice = $self->fetch_sequence;
  $self->query($slice);
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::EPCR->new
    (
     -query => $slice,
     -program => $self->analysis->program_file,
     -analysis => $self->analysis,
     %parameters
    );
  $self->runnable($runnable);
}


=head2 get_adaptor

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::EPCR
  Function  : get marker feature adaptor
  Returntype: Bio::EnsEMBL::Map::DBSQL::MarkerFeatureAdaptor
  Exceptions: none
  Example   : 

=cut

sub get_adaptor{
  my ($self) = @_;
  return $self->db->get_MarkerFeatureAdaptor;
}
