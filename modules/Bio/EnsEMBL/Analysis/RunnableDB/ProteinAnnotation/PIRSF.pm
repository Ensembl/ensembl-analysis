#
#
#
=pod 

=head1 NAME

  Bio::EnsEMBL::Pipeline::RunnableDB::ProteinAnnotation::PIRSF

=head1 SYNOPSIS

  my $seg = Bio::EnsEMBL::Pipeline::RunnableDB::ProteinAnnotation::PIRSF->
    new ( -db      => $db,
    -input_id   => $input_id,
    -analysis   => $analysis,
                                                                      );
  $seg->fetch_input;  # gets sequence from DB
  $seg->run;
  $seg->write_output; # writes features to to DB

=head1 DESCRIPTION

  This object wraps Bio::EnsEMBL::Pipeline::Runnable::Hmmpfam
  to add functionality to read and write to databases in 
  a Pfam-specific way.

=head1 CONTACT

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation::PIRSF;

use strict;
use vars qw(@ISA);

use Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation;
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::PIRSF;

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::ProteinAnnotation);


sub fetch_input {
  my ($self) = @_;
 
  $self->SUPER::fetch_input;
 
  # PIRSF needs sequence lengths, so we must create one Runnable
  # for each input sequence, giving a Bio::PrimarySeq as input

  my @query_seqs;
  if (ref($self->query) and $self->query->isa("Bio::PrimarySeq")) {
    push @query_seqs, $self->query;
  } elsif (-e $self->query) {
    my $seqio = Bio::SeqIO->new(-format => 'fasta',
                                -file   => $self->query);
    while(my $seq = $seqio->next_seq) {
      push @query_seqs, $seq;
    }
    $seqio->close;
  } else {
    throw($self->query . " is neither an object ref nor a file name\n");
  }

  foreach my $seq (@query_seqs) {
    my $run = Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::PIRSF->
        new(-query     => $seq,
            -analysis  => $self->analysis,
            -database  => $self->analysis->db_file,
            %{$self->parameters_hash}
            );
    $self->runnable($run);
  }
}



