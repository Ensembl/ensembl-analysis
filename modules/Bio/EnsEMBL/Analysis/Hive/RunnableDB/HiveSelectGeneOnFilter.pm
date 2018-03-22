=head1 LICENSE

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSelectGeneOnFilter

=head1 SYNOPSIS


=head1 DESCRIPTION

Module to filter genes by modifying the biotypes based on a filter object

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveSelectGeneOnFilter;

use strict;
use warnings;

use Bio::EnsEMBL::Gene;

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 fetch_input

 Arg [1]    : None
 Description: Fetch transcripts on the given slice. Autoflow will stop if
              no transcripts are found
 Returntype : None
 Exceptions : None

=cut

sub fetch_input {
  my ($self) = @_;

  $self->create_analysis;
  my $dna_db = $self->get_database_by_name('dna_db');
  my $db = $self->get_database_by_name('source_db', $dna_db);
  $self->hrdb_set_con($db);
  my $slice = $self->fetch_sequence;
  my $transcripts = $slice->get_all_Transcripts;
  if (scalar(@$transcripts)) {
    $self->param('transcripts', $slice->get_all_Transcripts);
  }
  else {
    $self->input_job->autoflow(0);
    $self->complete_early('No transcripts to process');
  }
  if ($self->param_is_defined('target_db')) {
    $self->hrdb_set_con($self->get_database_by_name('target_db', $dna_db));
  }
}


=head2 run

 Arg [1]    : None
 Description: Filter the transcripts based on the filter given in 'FILTER->{OBJECT}'
              and update the biotype of the transcripts
 Returntype : None
 Exceptions : None

=cut

sub run {
  my ($self) = @_;

  my @output;
  foreach my $transcript (@{$self->filter->filter_results($self->param('transcripts'))}) {
    my $gene;
    if (exists $transcript->{_gb_flag}) {
      $gene = Bio::EnsEMBL::Gene->new();
      $gene->add_Transcript($transcript);
    }
    else {
      $gene = $transcript->get_Gene;
    }
    $gene->analysis($self->analysis);
    $gene->biotype($self->analysis->logic_name);
    push(@output, $gene);
  }
  $self->output(\@output);
}


=head2 write_output

 Arg [1]    : None
 Description: Write the genes in a new database if 'target_db' is provided or
              update the biotype in the 'source_db'
 Returntype : None
 Exceptions : Throws if it fails to store genes

=cut

sub write_output {
  my ($self) = @_;

  my $dba = $self->get_adaptor;
  $self->hrdb_get_con->get_AnalysisAdaptor->store($self->analysis);
  if ($self->param_is_defined('target_db')) {
    foreach my $gene (@{$self->output}) {
      empty_Gene($gene);
      $dba->store($gene);
    }
  }
  else {
    foreach my $gene (@{$self->output}) {
      $dba->update($gene);
    }
  }
}


=head2 get_adaptor

 Arg [1]    : None
 Description: Getter for the Gene adaptor.
 Returntype : Bio::EnsEMBL::DBSQL::GeneAdaptor
 Exceptions : None

=cut

sub get_adaptor {
  my ($self) = @_;

  if (!$self->param_is_defined('_gene_adaptor')) {
    $self->param('_gene_adaptor', $self->hrdb_get_con->get_GeneAdaptor);
  }
  return $self->param('_gene_adaptor');
}


=head2 filter

 Arg [1]    : String $val, perl module name of the filter
 Description: Getter/setter to retrieve a filter object or to create the
              object based on Arg[1]
 Returntype : Object
 Exceptions : Throws if the object cannot be create when giving a parameter

=cut

# The module should inherit from HiveExonerate2Genes_cdna to get the filter
# method or from another model having a similar method
sub filter {
  my ($self, $val) = @_;
  if ($val) {
    $self->param('_transcript_filter',$val);
  }
  elsif (!$self->param_is_defined('_transcript_filter')
    and $self->param_is_defined('FILTER')
    and exists $self->param('FILTER')->{OBJECT}) {
    my $module = $self->require_module($self->param('FILTER')->{OBJECT});
    $self->param('_transcript_filter', $module->new(%{$self->param('FILTER')->{PARAMETERS}}));
  }
  if ($self->param_is_defined('_transcript_filter')) {
    return $self->param('_transcript_filter');
  }
  else {
    return;
  }
}

1;
