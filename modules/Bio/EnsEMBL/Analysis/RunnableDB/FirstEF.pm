# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::FirstEF
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::FirstEF

=head1 SYNOPSIS

  my $runnabledb = Bio::EnsEMBL::Analysis::RunnableDB::FirstEF->
  new(
      -input_id => 'contig::AL805961.22.1.166258:1:166258:1',
      -db => $db,
      -analysis => $analysis,
     );
  $runnabledb->fetch_input;
  $runnabledb->run;
  $runnabledb->write_output;


=head1 DESCRIPTION

This module provides an interface between the ensembl database and
the Runnable FirstEF which wraps the program FirstEF

This module can fetch appropriate input from the database
pass it to the runnable then write the results back to the database
in the simple_feature table 

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::FirstEF;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Runnable::FirstEF;
use Bio::EnsEMBL::Analysis::Config::General qw(PARAMETERS_DIR PARSE_SCRIPT) ;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB);


=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::FirstEF
  Function  : fetch data out of database and create runnable
  Returntype: 1
  Exceptions: none
  Example   : 

=cut



sub fetch_input{
  my ($self) = @_;
  my $slice = $self->fetch_sequence($self->input_id, $self->db, ['']); 
  #hard coded array is used here as it always wants all repeats in 
  #the table masked at somepoint the general variables will be replaced by
  #analysis specific variables and this will move back into confi
  $self->query($slice);

  if(!$self->analysis->program_file){
    $self->analysis->program_file('firstef');
  }


  my %parameters;
  if($self->parameters_hash){
    %parameters = %{$self->parameters_hash};
  }

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::FirstEF->new
    (
     -query => $self->query,
     -program => $self->analysis->program_file,
     -analysis => $self->analysis,
     -param_dir => $PARAMETERS_DIR,
     -parse_script => $PARSE_SCRIPT,
     %parameters,
    );
  $self->runnable($runnable);
  return 1;
}


=head2 get_adaptor

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::FirstEF
  Function  : get simple feature adaptor
  Returntype: Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor
  Exceptions: none
  Example   : 

=cut


sub get_adaptor{
  my ($self) = @_;
  return $self->db->get_SimpleFeatureAdaptor;
}


1;
