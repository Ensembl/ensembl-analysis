
package Bio::EnsEMBL::Analysis::Runnable::Spliced_elsewhere;


use strict;
use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Analysis::Config::Pseudogene;
use Bio::EnsEMBL::Analysis::Runnable::Blast;
use Bio::EnsEMBL::Analysis::Tools::BPliteWrapper;
use Bio::EnsEMBL::Analysis::Tools::FilterBPlite;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->{'_genes'} = [];	#array of genes to test;  

  my($genes) = $self->_rearrange([qw(
				     GENES
				    )], @args);
  if ($genes) {
    $self->genes($genes);
  }
  $self->db("$PS_MULTI_EXON_DB"."multi_exon_blastdb.fasta");
  return $self;
}



sub run {
  my ($self)=@_;
  my @genes = @{$self->genes};
  foreach my $gene (@genes){
    foreach my $trans(@{$gene->get_all_Transcripts}){
      if ($trans->translateable_seq){
	$self->run_blast($trans);
      }
    }
  }
  return 1;
}


sub  run_blast{
  my ($self,$trans)=@_;
  my $bplitewrapper = Bio::EnsEMBL::Analysis::Tools::BPliteWrapper-> new
    (
     -query_type => 'pep',
     -database_type => 'pep',
    );

  my $query = $trans->translate;
  my $blast =  Bio::EnsEMBL::Analysis::Runnable::Blast->new 
    ('-query'     => $query,
     '-program'   => 'blastp',
     '-database'  => $self->db,
     '-threshold' => 1e-6,
     '-parser'    => $bplitewrapper,
     '-options'   => 'V=10',
     '-analysis'  => $self->analysis,
    );

  $blast->run();
  my @array =($trans, $blast->output);
  # maybe sort output based on evalue and chuck out dodgy hits
  # or put them into the parameters for the blast?

  push @{$self->output},\@array;
  return 1;
}



########################################################
# Containers

=head2 genes

Arg [1]    : array ref
  Description: get/set gene set to run over
  Returntype : array ref to Bio::EnsEMBL::Gene objects
  Exceptions : throws if not a Bio::EnsEMBL::Gene
  Caller     : general

=cut

sub genes {
  my ($self, $genes) = @_;
  if ($genes) {
    foreach my $gene (@{$genes}) {
      unless  ($gene->isa("Bio::EnsEMBL::Gene")){
	$self->throw("Input isn't a Bio::EnsEMBL::Gene, it is a $gene\n$@");
      }
    }
    $self->{'_genes'} = $genes;
  }
  return $self->{'_genes'};
}

sub db {
  my ($self, $db) = @_;
  if ($db) {
    $self->find_file($db);
    $self->{'_db'} = $db;
  }
  return $self->{'_db'};
}

1;
