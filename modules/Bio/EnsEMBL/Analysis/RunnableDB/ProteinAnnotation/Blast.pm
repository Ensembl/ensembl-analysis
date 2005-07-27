#
#
#
=pod 

=head1 NAME

  Bio::EnsEMBL::Pipeline::RunnableDB::ProteinAnnotation::Blast

=head1 SYNOPSIS

  my $seg = Bio::EnsEMBL::Pipeline::RunnableDB::ProteinAnnotation::Blast
    ->new ( -db      => $db,
    -input_id   => $input_id,
    -analysis   => $analysis,
    );
  $seg->fetch_input;  # gets sequence from DB
  $seg->run;
  $seg->write_output; # writes features to to DB

=head1 DESCRIPTION

This object uses the Bio::Ensembl::Analysis::Blast runnable to identify
protein domasn using Blast (a la Prodom)

=head1 CONTACT

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::Blast;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation;

use Bio::EnsEMBL::Analysis::Runnable::Blast;
use Bio::EnsEMBL::Analysis::Tools::BPliteWrapper;
use Bio::EnsEMBL::Analysis::Tools::ProteinAnnotationFilter;


@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation);


#
# overridden methods
#
sub fetch_input {
  my ($self) = @_;

  $self->SUPER::fetch_input;

  my @dbfiles = split(/;/,$self->analysis->db_file);

  my $parser = Bio::EnsEMBL::Analysis::Tools::BPliteWrapper->
      new(-regex => '^[^\#]+\#([^\#]+)\#',
          -query_type => 'pep',
          -database_type => 'pep');
          
  my $filter = Bio::EnsEMBL::Analysis::Tools::ProteinAnnotationFilter->new();

  my @queries;
  if (ref($self->query) and $self->query->isa("Bio::PrimarySeqI")) {
    push @queries, $self->query;
  } else {
    # must be a file
    my $seqio = Bio::SeqIO->new(-format => 'fasta',
                                -file   => $self->query);
    while(my $seq = $seqio->next_seq) {
      push @queries, $seq;
    }
  }

  foreach my $q (@queries) {
    foreach my $f (@dbfiles) {
      my $run = Bio::EnsEMBL::Analysis::Runnable::Blast->
          new(-program   => $self->analysis->program_file,
              -query     => $q,
              -database  => $f,
              -analysis  => $self->analysis,
              -parser    => $parser,
              -filter    => $filter,
              %{$self->parameters_hash}
              );
      $self->runnable($run);    
    }
  }
}

sub write_output {
  my ($self) = @_;

  my @pfs;
  foreach my $f (@{$self->output}) {
    $f->analysis($self->analysis);
  }

  $self->output(\@pfs);
  $self->SUPER::write_output;
}

1;
