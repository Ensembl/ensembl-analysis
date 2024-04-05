=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.  

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::BlastTranscriptPep - 

=head1 SYNOPSIS

my $db          = Bio::EnsEMBL::DBAdaptor->new($locator);
my $btpep     = Bio::EnsEMBL::Analysis::RunnableDB::BlastTranscriptPep->new ( 
                                                    -dbobj      => $db,
			                            -input_id   => $input_id
                                                    -analysis   => $analysis );

$btpep->fetch_input();
$btpep->run();
$btpep->output();

=head1 DESCRIPTION

This object runs Bio::EnsEMBL::Analysis::Runnable::Blast on peptides
obtained by translating a representative transcript from each gene
in the region. The resulting blast hits are written back as
DnaPepAlignFeature's.

The appropriate Bio::EnsEMBL::Analysis object must be passed for
extraction of appropriate parameters.

=head1 METHODS

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _'

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::BlastTranscriptPep;

use warnings ;
use strict;

use Bio::EnsEMBL::Analysis::RunnableDB::Blast;
use Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::Blast);


=head2 fetch_input

  Args       : none
  Example    : $runnable->fetch_input
  Description: Fetches input data for BlastTranscriptPep and makes runnable
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
    
  foreach my $t (@representative_trans) {
    foreach my $db (split ',', ($self->analysis->db_file)) {
      $self->runnable(Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep
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
