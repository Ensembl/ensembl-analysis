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

Bio::EnsEMBL::Analysis::RunnableDB::EPCR - 

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


=head1 METHODS

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
