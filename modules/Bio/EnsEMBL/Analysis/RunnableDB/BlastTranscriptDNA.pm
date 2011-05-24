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

Bio::EnsEMBL::Analysis::RunnableDB::BlastTranscriptDNA - 

=head1 SYNOPSIS

my $btdna     = Bio::EnsEMBL::Analysis::RunnableDB::BlastTranscriptDNA->new ( 
                                                    -db         => $db,
			                            -input_id   => $input_id
                                                    -analysis   => $analysis );

$btdna->fetch_input();
$btdna->run();
$btdna->output();

=head1 DESCRIPTION

This object runs Bio::EnsEMBL::Analysis::Runnable::Blast of 
the translation of a representative transcript if the genes in
the region, against a DNA database. The resulting blast hits are 
written back as DnaDnaAlignFeature's.


=head1 METHODS

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


package Bio::EnsEMBL::Analysis::RunnableDB::BlastTranscriptDNA;

use strict;

use Bio::EnsEMBL::Analysis::RunnableDB::Blast;
use Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptDNA;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::Blast);


=head2 fetch_input

  Args       : none
  Example    : $runnable->fetch_input
  Description: Fetches input data for BlastTranscriptDNA and makes runnable
  Returntype : none
  Exceptions : $self->input_id is not defined
  Caller     : run_RunnableDB, Bio::EnsEMBL::Pipeline::Job

=cut

sub fetch_input {
  my($self) = @_;
  
  $self->throw("No input id") unless defined($self->input_id);

  my $slice = $self->fetch_sequence;
  $self->query($slice);

  my %blast = %{$self->BLAST_PARAMS};
  my $parser = $self->make_parser;
  my $filter;
  if($self->BLAST_FILTER){
    $filter = $self->make_filter;
  }

  my (@representative_trans);
  my $genes = $self->db->get_GeneAdaptor->fetch_all_by_Slice($self->query);
  # get longest transcript; only sensible to do it for that
  foreach my $gene (@$genes) {
    next if $gene->start < 1;

    my ($longest_pep_tran, $longest_pep_tran_len);
    foreach my $tran (@{$gene->get_all_Transcripts}) {
      if ($tran->translation) {
        my $pep_string = $tran->translate->seq;
        
        if (not defined $longest_pep_tran or length($pep_string) > $longest_pep_tran_len) {
          $longest_pep_tran = $tran;
          $longest_pep_tran_len = length($pep_string);
        }
      }
    }
    push @representative_trans, $longest_pep_tran;
  }
  # $self->_transcripts(@representative_trans);


    
  foreach my $t (@representative_trans) {
    foreach my $db (split ',', ($self->analysis->db_file)) {
      $self->runnable(Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptDNA
                        ->new(
                              -transcript     => $t,
                              -query          => $self->query,                              
                              -program        => $self->analysis->program_file,
                              -parser         => $parser,
                              -filter         => $filter,
                              -database       => $db,
                              -analysis       => $self->analysis,
                              %blast,
                             ));
    }
  }

  return 1;

}


1;
