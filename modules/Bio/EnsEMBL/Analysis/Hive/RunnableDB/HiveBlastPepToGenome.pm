=head1 LICENSE

# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#Copyright [2016-2020] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastPepToGenome

=head1 SYNOPSIS

  my $blast = Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastPepToGenome->new
     (
      -analysis => $analysis,
      -db       => $db,
      -input_id => 'contig::AL1347153.1.3517:1:3571:1'
     );
  $blast->fetch_input;
  $blast->run;
  my @output =@{$blast->output};

=head1 DESCRIPTION

  Blast protein sequences against a genome one seqeunce at a time and
  store the alignment in a database

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastPepToGenome;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable::Blast;

use parent('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlast');



=head2 fetch_input

  Arg [1]   : None
  Function  : Fetch multiple protein sequence from a database using 'iid', add tr or sp
              to the logic_name to provide information about the status of the protein
              and create the Runnable
  Returntype: None
  Exceptions: None

=cut

sub fetch_input{
  my ($self) = @_;

  $self->create_analysis();
  $self->analysis->program($self->param('blast_program')) if ($self->param_is_defined('blast_program'));
  $self->analysis->program_file($self->param('blast_exe_path')) if ($self->param_is_defined('blast_exe_path'));
  $self->analysis->parameters($self->param('commandline_params')) if ($self->param_is_defined('commandline_params'));
  $self->analysis->db_file($self->param('blast_db_path')) if ($self->param_is_defined('blast_db_path'));
  $self->analysis->db($self->param('blast_db_name')) if ($self->param_is_defined('blast_db_name'));

  my $output_dba = $self->hrdb_get_dba($self->param('output_db'));
  if ($self->param_is_defined('dna_db')) {
    my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
    if($dna_dba) {
      $output_dba->dnadb($dna_dba);
    }
  }

  $self->hrdb_set_con($output_dba,'output_db');

  my $blast_params = $self->param('BLAST_PARAMS');
  my $blast_filter = $self->param('BLAST_FILTER');
  my $blast_parser = $self->param('BLAST_PARSER');
  my $filter_params = $self->param('FILTER_PARAMS');
  my $parser_params = $self->param('PARSER_PARAMS');

  my %blast = %{$blast_params};

  my $parser;
  if($parser_params) {
    $parser = $self->make_parser($parser_params);
  }

  my $filter;
  if($blast_filter){
    $filter = $self->make_filter($filter_params);
  }


  # submit blast module to use via analysis_parameters column of analysis table
  my $options_string ;
  my %options = %{$parser_params};

  my $query_seqs = $self->get_query_seqs($self->input_id);
  my %analysis_hash;
  my $analysis = $self->analysis;
  foreach my $db ('sp', 'tr') {
    my %new_analysis = %$analysis;
    $analysis_hash{$db} = ref($analysis)->new_fast(\%new_analysis);
    $analysis_hash{$db}->logic_name($analysis_hash{$db}->logic_name."_$db");
  }

  foreach my $seq (@{$query_seqs}) {
    if ($seq->namespace =~ /_tr/) {
      $analysis = $analysis_hash{tr};
    }
    else {
      $analysis = $analysis_hash{sp};
    }
    $parser->analysis($analysis);
    my $runnable = Bio::EnsEMBL::Analysis::Runnable::Blast->new(
                     -query => $seq,
                     -program => $self->param('blast_exe_path'),
                     -database => $self->param('blast_db_path'),
                     -parser => $parser,
                     -filter => $filter,
                     -analysis => $analysis,
                     -options => $self->param('commandline_params'),
                     %blast,
                   );
    if ($self->param_is_defined('timer')) {
      $runnable->timer($self->param('timer'));
    }
    $self->runnable($runnable);
  }
}


=head2 write_output

 Arg [1]    : None
 Description: Write the alignment information into the protein_align_feature table
              unless on of the protein blast lasted more than 'timer'. In this case,
              dataflow to the RUNLIMIT (-2) branch
 Returntype : None
 Exceptions : None

=cut

sub write_output {
  my ($self) = @_;

  if ($self->batch_failed) {
    $self->dataflow_output_id({iid => $self->input_id}, $self->param('_branch_to_flow_to_on_fail'));
  }
  else {
    my $out_dba = $self->hrdb_get_con('output_db');
    my $paf_adaptor = $out_dba->get_ProteinAlignFeatureAdaptor;
    my $slice_adaptor = $out_dba->get_SliceAdaptor;
    foreach my $hit ( @{$self->output} ) {
      my $slice = $slice_adaptor->fetch_by_name($hit->seqname);
      $hit->slice($slice);
      $paf_adaptor->store($hit);
    }
  }
}



=head2 get_query_seqs

 Arg [1]    : Arrayref of String, accessions to be fetched
 Description: Fetch the protein sequences provided in Arg[1] from the
              Hive pipeline database and add the biotype as namespace
              to the Bio::Seq object
 Returntype : Arrayref of Bio::Seq
 Exceptions : None

=cut

sub get_query_seqs {
  my ($self,$accession_array) = @_;

  my $table_adaptor = $self->db->get_NakedTableAdaptor();
  $table_adaptor->table_name($self->param_required('sequence_table_name'));

  my @seqs;

  foreach my $accession (@{$accession_array}) {
    my $db_row = $table_adaptor->fetch_by_dbID($accession);
    unless($db_row) {
      $self->throw("Did not find an entry in the uniprot_sequences table matching the accession. Accession:\n".$accession);
    }

    push(@seqs, Bio::Seq->new(-id => $accession, -seq => $db_row->{seq}, -namespace => $db_row->{'biotype'}));
  }
  return \@seqs;
}

1;
