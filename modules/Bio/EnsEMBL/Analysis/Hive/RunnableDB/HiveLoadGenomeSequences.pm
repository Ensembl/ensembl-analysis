=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadGenomeSequences

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadGenomeSequences;

use strict;
use warnings;
use feature 'say';
use File::Spec::Functions qw(catfile);
use Bio::EnsEMBL::Hive::Utils qw(destringify);
use Bio::EnsEMBL::IO::Parser::Fasta;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: Defaults parameters:
               replace_ambiguous_bases => 0,
 Returntype : Hashref
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    replace_ambiguous_bases => 0,
    load_mitochondrion => 0,
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: Check that the FASTA file containing the sequence exists and
              create the parser. Connect to the 'target_db' to fetch the
              coordinate system which will have the sequence and check how
              many regions should have DNA added.
              If 'filename' is set, it will use this file instead of the
              Genbank genomic file to load the sequences.
 Returntype : None
 Exceptions : Throws if 'assembly_name' is not set
              Throws if 'target_db' is not set
              Throws if it cannot find a coordinate system
              Throws if there is no region in the coordinate system
              Throws if 'output_path' or 'assembly_accession' is not set when
                'filename' is not defined

=cut

sub fetch_input {
  my ($self) = @_;

  my $fasta_file;
  if ($self->param_is_defined('filename')) {
# This is when we load the full assembly
  }
  elsif (-e catfile($self->param_required('output_path'), $self->param_required('assembly_accession').'_'.$self->param_required('assembly_name').'_genomic.fna')) {
    $fasta_file = catfile($self->param('output_path'), $self->param_required('assembly_accession').'_'.$self->param_required('assembly_name').'_genomic.fna');
  }
  else {
    $self->throw('Cannot load data from "filename" or the Genbank genomic file');
  }
  my $db = $self->get_database_by_name('target_db');
  $self->hrdb_set_con($db, 'target_db');
  my $coord_system = $db->get_CoordSystemAdaptor->fetch_sequence_level;
  $self->throw('Could not find any coordinate system with version '.$self->param('assembly_name'))
    unless ($coord_system);
  $self->param('coord_system', $coord_system);
  my $sequence_count = @{$db->get_SliceAdaptor->fetch_all($coord_system->name, $coord_system->version)};
  $self->throw('No sequence is linked to coordinate system '.$coord_system->name.' '.$coord_system->version)
    if ($sequence_count == 0);
  $self->param('sequence_count', $sequence_count);
  my $parser = Bio::EnsEMBL::IO::Parser::Fasta->open($fasta_file);
  $self->param('fasta_parser', $parser);
}


=head2 run

 Arg [1]    : None
 Description: Do nothing as processing the file here would make the module
              use as much memory than the genome size
 Returntype : None
 Exceptions : None

=cut

sub run {

# Just doing nothing maybe we could process the file and print when
# debug and no_write is used
}


=head2 write_output

 Arg [1]    : None
 Description: Load the DNA into the database
              If 'replace_ambiguous_bases' is true, replace the ambiguous
              bases, default is false.
              If the coordinate system is 'toplevel', there is no need to
              check the assembly. Otherwise we need to create a job for
              checking toplevel on branch '_branch_to_flow_to'.
 Returntype : None
 Exceptions : Throws if the number of expected sequences is different from
                the number of loaded sequences

=cut

sub write_output {
  my ($self) = @_;

  my $region_count = 0;
  my $seq_count = 0;
  my $parser = $self->param('fasta_parser');
  my $db = $self->hrdb_get_con('target_db');
  my $slice_adaptor = $db->get_SliceAdaptor;
  my $sequence_adaptor = $db->get_SequenceAdaptor;
  my $replace_ambiguous_bases = $self->param('replace_ambiguous_bases');
  my $coord_system = $self->param('coord_system');
  my $coord_system_name = $coord_system->name;
  my $has_mitochondrion = 0;
  while ($parser->next) {
    my ($first, $second) = $parser->getHeader =~ /^\w*\|([A-Za-z0-9.]+)|^([A-Za-z0-9.]+)/;
    $has_mitochondrion = 1 if ($parser->getHeader =~ /mitochondrion/);
    my $accession = $first || $second;
    my $slice = $slice_adaptor->fetch_by_region($coord_system_name, $accession);
    if ($slice) {
      eval {
        $sequence_adaptor->remove($slice->get_seq_region_id);
      };
      if ($replace_ambiguous_bases and $parser->getSequence =~ /[^ACGTN]+/i) {
        $self->warning($slice->name.' had at least one non-ATGCN (RYKMSWBDHV) base. Changed all to N.');
        my $seq_clean = $parser->getSequence;
        $seq_clean =~ tr/RYKMSWBDHV/N/;
        $sequence_adaptor->store($slice->get_seq_region_id, \$seq_clean);
      }
      else {
        $sequence_adaptor->store($slice->get_seq_region_id, $parser->getSequence);
      }
      ++$region_count;
    }
    else {
      $self->warning("Did not found $accession in database");
    }
    ++$seq_count;
  }
  if (!$self->param('load_mitochondrion') and $has_mitochondrion) {
    --$seq_count;
  }
  if ($region_count != $seq_count) {
    $self->warning("You have a different number of sequence in your file and in your database: $seq_count !! $region_count");
  }
  if ($seq_count != $self->param('sequence_count')) {
    $self->throw("Not all sequence have DNA as $seq_count is different from expected ".$self->param('sequence_count'));
  }
  $self->dataflow_output_id(destringify($self->input_job->input_id), $self->param('_branch_to_flow_to'))
    unless ($coord_system->rank == 1 and $coord_system->is_sequence_level);
}

1;
