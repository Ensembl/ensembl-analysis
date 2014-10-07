package Bio::EnsEMBL::Analysis::RunnableDB::ProjectionCandidates;

use warnings;
use strict;
use feature 'say';

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);

use parent ('Bio::EnsEMBL::Analysis::RunnableDB::HiveBaseRunnable');


sub fetch_input {
  my $self = shift;
  $self->get_input_genes($self->param('iid'));
}



sub run {
  my $self = shift;
  $self->build_single_transcript_genes();
  return 1;
}



sub write_output {
  my $self = shift;

  $self->db('output_db',$self->get_dba($self->param('proj_candidate_db')));
  my $output_gene_adaptor = $self->db('output_db')->get_GeneAdaptor();

  foreach my $gene (@{$self->{'output_genes'}}) {
    empty_Gene($gene);
    $output_gene_adaptor->store($gene);
  }

  return 1;
}



sub get_input_genes {
  my ($self,$slice_name) = @_;

  $self->db('input_db',$self->get_dba($self->param('reference_db')));
  my $sa = $self->db('input_db')->get_SliceAdaptor();
  my $slice = $sa->fetch_by_name($slice_name);

  $self->{'genes'} = $slice->get_all_Genes();
}



sub build_single_transcript_genes {
  my $self = shift;
  $self->{'output_genes'} = [];

  foreach my $gene (@{$self->{'genes'}}) {

    foreach my $transcript (@{$gene->get_all_Transcripts()}) {

      unless(scalar(@{$transcript->get_all_Attributes('gencode_basic')})) {
        next;
      }

      my $new_gene = Bio::EnsEMBL::Gene->new(
              -START  => $transcript->start(),
              -END    => $transcript->end(),
              -STRAND => $transcript->strand(),
              -SLICE  => $transcript->slice(),
              -ANALYSIS => $transcript->analysis(),
      );

      $new_gene->add_Transcript($transcript);
      push(@{$self->{'output_genes'}},$new_gene);
    }
  }
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
                                                  -dbname => $$connection_info{'-dbname'},
                                                  -user => $$connection_info{'-user'},
                                                  -species => $$connection_info{'-species'},
                                                  -host => $$connection_info{'-host'},
                                                  -port => $$connection_info{'-port'},
                                                  -pass => $$connection_info{'-pass'},
                                                );
    }

  $dba->dbc->disconnect_when_inactive(1) ;
  return $dba;

}

1;
