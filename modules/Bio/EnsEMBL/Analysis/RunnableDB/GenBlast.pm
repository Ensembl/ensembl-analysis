# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::GenBlast
#
# Copyright (c) 2009 WormBase
#

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::GenBlast

=head1 SYNOPSIS

  my $blat = Bio::EnsEMBL::Analysis::RunnableDB::GenBlast->
  new(
      -input_id => 'file_name',
      -db => $db,
      -analysis => $analysis,
     );
  $blat->fetch_input;
  $blat->run;
  $blat->write_output;

=head1 DESCRIPTION

This module provides an interface between the ensembl database and
the Runnable GenBlast which wraps the program GenBlast

This module can fetch appropriate input from the database
pass it to the runnable then write the results back to the database
in the dna_align_feature  tables 

=head1 CONTACT

Post questions to the Ensembl development list: dev@ensembl.org

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::GenBlast;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Runnable::GenBlast;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB);


=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::GenBlast
  Function  : fetch data out of database and create runnable
  Returntype: 1
  Exceptions: none
  Example   : 

=cut



sub fetch_input {
  my ($self) = @_;
  my (%parameters, %genome_slices);

  if($self->parameters_hash){
    %parameters = %{$self->parameters_hash};
  }

  my $genome_file = $self->analysis->db_file;
  open(my $fh, $genome_file) or throw("Could not open $genome_file for reading");
  while(<$fh>) {
    /^\>(\S+)/ and do {
      my $seq_name = $1;
      my $slice = $self->db->get_SliceAdaptor->fetch_by_region('toplevel', $seq_name);
      if (not defined $slice) {
        throw("Could not extract slice for $seq_name from database");
      }
      $genome_slices{$seq_name} = $slice;
    }
  }

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::GenBlast->new
    (
     -query => $self->input_id,
     -program => $self->analysis->program_file,
     -analysis => $self->analysis,
     -database => $self->analysis->db_file,
     -refslices => \%genome_slices,
     %parameters,
    );
  $self->runnable($runnable);
  return 1;
}


=head2 write_output

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::GenBlast
  Function  : writes the prediction transcripts back to the database
  after validation
  Returntype: none
  Exceptions: 
  Example   : 

=cut



sub write_output{
  my ($self) = @_;
  my $adaptor = $self->db->get_PredictionTranscriptAdaptor;
  my @output = @{$self->output};
  my $ff = $self->feature_factory;
  foreach my $pt(@output){
    $pt->analysis($self->analysis);
    $pt->slice($self->query) if(!$pt->slice);
#    if ($self->test_translates()) {
#      print "The GenBlast transcript ",$pt->display_label," doesn't translate correctly\n";
#      next;
#    } # test to see if this transcript translates OK
    $adaptor->store($pt);
  }
}

=head2 test_translates

  Arg [1]   : Bio::EnsEMBL::PredictionTranscript
  Function  : tests whether a transcript translates correctly
  Returntype: int 1 for failure, 0 for OK
  Exceptions: 
  Example   : 

=cut

sub test_translates {
  my ($pt) = @_;
  my $result = 0;
  my $tseq;
  eval{
    $tseq = $pt->translate;
  };
  if (!$tseq || $tseq->seq =~ /\*/) {
    print "$tseq->seq\n" if $tseq;
    $result = 1;
  }
  return $result; 
}


1;
