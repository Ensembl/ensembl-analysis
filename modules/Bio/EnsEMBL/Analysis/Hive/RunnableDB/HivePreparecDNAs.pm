=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2024] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

Please email comments or questions to the public Ensembl
developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

Questions may also be sent to the Ensembl help desk at
<http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivePreparecDNAs

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HivePreparecDNAs;

use strict;
use warnings;

use Bio::SeqIO;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: Default parameters are:
               format => 'fasta',
               _suffix => 'cleaned', # Extension added to the new file created
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    format => 'fasta',
    _suffix => 'cleaned',
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Check few parameters and set the output file name to be <input filename>.<_suffix>
 Returntype : None
 Exceptions : Throws if 'data_type' is not set or not equal to 'ena' or 'refseq'
              Throws if 'description_filter' is not set when data_type = refseq
              Throws if 'filename' does not exist
              Throws if 'filename'.'_suffix' exists

=cut

sub fetch_input {
  my ($self) = @_;

  my $data_type = $self->param_required('data_type');
  if ($data_type eq 'refseq') {
    $self->param_required('description_filter');
  }
  elsif ($data_type ne 'ena') {
    $self->throw('Cannot process data from '.$data_type);
  }
  my $source_file = $self->param_required('filename');
  $self->throw('File does not exist') unless (-e $source_file);
  my $new_file = $source_file.'.'.$self->param('_suffix');
  $self->throw("File $new_file exists") if (-e $new_file);
  $self->param('target_file', $new_file);
}


=head2 run

 Arg [1]    : None
 Description: Parse the input file and write only the sequence which
              satisfy the filters. At the moment we only take NM_ and NR_
              from RefSeq for the species specified in 'description_filter'.
              We only take mRNA from the ENA.
              If no sequence is kept, no file is created and dataflow is stopped
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;

  my $method = '_is_'.$self->param('data_type').'_sequence_ok';
  my $source = Bio::SeqIO->new(-format => $self->param('format'), -file => $self->param('filename'));
  my $target = Bio::SeqIO->new(-format => 'fasta', -file => '>'.$self->param('target_file'));
  my $seq_count = 0;
  my $filter = $self->param_is_defined('description_filter') ? $self->param('description_filter') : undef;
  while(my $seq = $source->next_seq) {
    if ($self->$method($seq, $filter)) {
      $seq->id($seq->id.'.'.$seq->version) if ($seq->version);
      $target->write_seq($seq);
      ++$seq_count;
    }
  }
  if ($seq_count) {
    $self->output([$self->param('target_file')]);
  }
  else {
    unlink $self->param('target_file');
    $self->input_job->autoflow(0);
    $self->complete_early('No sequence to write to file');
  }
}


=head2 _is_refseq_sequence_ok

 Arg [1]    : Bio::Seq, the sequence to test
 Arg [2]    : String $regex, a regex to test on the description of the sequence
 Description: Return 1 if the sequence is a NM|NR and from the species specified in Arg[2].
              Otherwise return 0
 Returntype : Boolean
 Exceptions : None

=cut

sub _is_refseq_sequence_ok {
  my ($self, $seq, $filter) = @_;

  return ($seq->id =~ /N[MR]_/ and $seq->desc =~ /$filter/);
}


=head2 _is_ena_sequence_ok

 Arg [1]    : Bio::Seq, the sequence to test
 Description: Return 1 if the sequence is a mRNA, otherwise return 0
 Returntype : Boolean
 Exceptions : None

=cut

sub _is_ena_sequence_ok {
  my ($self, $seq) = @_;

  return $seq->molecule eq 'mRNA';
}


=head2 write_output

 Arg [1]    : None
 Description: Write into branch '_branch_to_flow_to' the name of the file containing
              cDNA sequences as 'filename'
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;

  $self->dataflow_output_id([{filename => $self->output->[0]}], $self->param('_branch_to_flow_to'));
}


1;
