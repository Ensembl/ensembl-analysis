package Bio::EnsEMBL::Analysis::RunnableDB::HiveGeneBuilder;

use warnings;
use strict;
use feature 'say';

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);

use parent ('Bio::EnsEMBL::Analysis::RunnableDB::HiveBaseRunnable');


sub fetch_input {
  my $self = shift;
#  $self->get_input_genes($self->param('iid'));
}



sub run {
  my $self = shift;

}



sub write_output {
  my $self = shift;

}



sub get_input_genes {
  my ($self,$slice_name) = @_;

  $self->db('input_db',$self->get_dba($self->param('reference_db')));
  my $sa = $self->db('input_db')->get_SliceAdaptor();
  my $slice = $sa->fetch_by_name($slice_name);

  $self->{'genes'} = $slice->get_all_Genes();
}


sub db {
  my ($self,$adaptor_name,$value) = @_;

  if($value){
    $self->{$adaptor_name} = $value;
  }
  return $self->{$adaptor_name};
}



sub get_dba {
   my ($self,$connection_info) = @_;
   my $dba;

   if (ref($connection_info)=~m/HASH/) {
      $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                                  %$connection_info,

                                               );
    }

  $dba->dbc->disconnect_when_inactive(1) ;
  return $dba;

}

1;
