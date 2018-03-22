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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadUpdatecDNAs

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveLoadUpdatecDNAs;

use strict;
use warnings;

use Bio::Seq;
use Bio::EnsEMBL::IO::Parser::Fasta;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(hrdb_get_dba);
use Bio::EnsEMBL::Analysis::Tools::PolyAClipping qw(clip_if_necessary);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(empty_Object);
use Bio::EnsEMBL::KillList::KillList;
use Bio::EnsEMBL::UnmappedObject;
use Bio::EnsEMBL::Analysis;

use parent ('Bio::EnsEMBL::Hive::RunnableDB::JobFactory');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults()},
    sequence_biotype => 'cdna',
    column_names => ['iid'],
    sequence_table_name => 'cdna_sequences',
    min_length => 60,
    logic_name => 'cdna_update',
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: If 'strategy' is "update", it first retrieve all aligned and non aligned sequence
              from the cdna database 'old_cdna_db' and load the new sequences into the hive
              database. Then it writes a file of gene_ids 'retire_gene_file' to be able to
              remove from the 'new_cdna_db' any models based on sequences not present in 'cdna_file'.
              If 'strategy' is "complete", it simply loads all the sequences from 'cdna_file'.
 Returntype : None
 Exceptions : Throws if 'strategy' is not set
              Throws if it cannot open or close 'retire_gene_file' when in update mode
              Throws if it cannot query the databases when in update mode

=cut

sub fetch_input {
  my $self = shift;

  my $filename = $self->param_required('filename');
  my $format = $self->param_is_defined('format') ? $self->param('format') : 'fasta';
  if (!$self->param_is_defined('format') and $filename =~ /\.dat$/) {
    $format = 'embl'
  }
  my $parser = Bio::SeqIO->new(-format => $format, -file => $filename);

  # create a hash that has the cdnas that have already been aligned previously as the key
  my %seen_cdna;

# This is needed for the kill list
  my $new_db = hrdb_get_dba($self->param_required('new_cdna_db'));
  $self->param('new_db', $new_db);

  # only need to get the hit names if we're only doing an update and not complete
  if ($self->param_required('strategy') eq 'update') {
    my $old_db = hrdb_get_dba($self->param_required('old_cdna_db'));
    my $sth = $old_db->dbc()->prepare('SELECT DISTINCT(hit_name) FROM dna_align_feature');
    $sth->execute();
    while ( my $cdna  = $sth->fetchrow_array ) {
      $seen_cdna{$cdna} = 1;
    }
    $sth = $old_db->dbc()->prepare('SELECT DISTINCT(identifier) FROM unmapped_object');
    $sth->execute();
    while ( my $cdna = $sth->fetchrow_array ) {
      $seen_cdna{$cdna} = 1;
    }
    $sth->finish();
    $self->param('old_db', $old_db);
    my $dna_db = hrdb_get_dba($self->param_required('dna_db'));
    $old_db->dnadb($dna_db);
    $new_db->dnadb($dna_db);
  }
  my $biotype = $self->param('sequence_biotype');

  my $table_adaptor = $self->db->get_NakedTableAdaptor();
  $table_adaptor->table_name($self->param('sequence_table_name'));

  my @iids;
  my @to_copy;
  my $kill_list_object = Bio::EnsEMBL::KillList::KillList->new( -TYPE => 'cDNA', -KILL_LIST_DB => $self->param('killlist_db'), -FILTER_PARAMS => { -for_analyses => ['cdna_update'], -only_mol_type => 'cDNA'});
  my $kill_list = $kill_list_object->get_kill_list();
  my $min_length = $self->param('min_length');
  my @killed;
  while(my $seq = $parser->next_seq) {
    my $accession = $seq->display_id;
    my ($id) = $accession =~ /^(.*)\.\d+$/;

    if (exists $seen_cdna{$accession}) {
      delete $seen_cdna{$accession};
      push(@to_copy, $accession);
    }
    else {
      if (exists $kill_list->{$id}) {
        push(@killed, $accession);
        next;
      }
      my ($clipped, undef, undef) = clip_if_necessary($seq);
      next unless (defined $clipped and $clipped->length > $min_length);
      my $db_row = [{
        'accession'  => $accession,
        'seq'        => $clipped->seq,
        'biotype'    => $biotype,
      }];
      $table_adaptor->store($db_row);

      push(@iids, $accession);
    }
  }
  $self->param('inputlist', \@iids);
  $self->param('genes_to_copy', \@to_copy);
  if (@killed) {
    $self->param('genes_killed', \@killed) if (@killed);
  }
}

sub write_output {
  my ($self) = @_;

  $self->SUPER::write_output;
  if (scalar(@{$self->param('genes_to_copy')})) {
    my $old_db = $self->param('old_db');
    my $old_gene_adaptor = $old_db->get_GeneAdaptor;
    my $old_unmapped_adaptor = $old_db->get_UnmappedObjectAdaptor;
    my $new_db = $self->param('new_db');
    my $new_gene_adaptor = $new_db->get_GeneAdaptor;
    my $new_unmapped_adaptor = $new_db->get_UnmappedObjectAdaptor;
    foreach my $accession (@{$self->param('genes_to_copy')}) {
      my $genes = $old_gene_adaptor->fetch_all_by_transcript_supporting_evidence($accession, 'dna_align_feature');
      if ($genes and scalar(@$genes)) {
        foreach my $gene (@$genes) {
          empty_Gene($gene);
          $new_gene_adaptor->store($gene);
        }
      }
      else {
        foreach my $unmapped_object (@{$old_unmapped_adaptor->fetch_by_identifier($accession)}) {
          empty_Object($unmapped_object);
          $new_unmapped_adaptor->store($unmapped_object);
        }
      }
    }
  }
  if ($self->param_is_defined('genes_killed')) {
    my $analysis = $self->param('new_db')->get_AnalysisAdaptor->fetch_by_logic_name($self->param('logic_name'));
    my $unmapped_adaptor = $self->param('new_db')->get_UnmappedObjectAdaptor;
    if (!$analysis) {
      $analysis = $self->param('old_db')->get_AnalysisAdaptor->fetch_by_name($self->param('logic_name'));
      if ($analysis) {
        empty_Object($analysis);
      }
      else {
        $analysis = Bio::EnsEMBL::Analysis->new(-logic_name => $self->param('logic_name'));
      }
    }
    foreach my $id (@{$self->param('genes_killed')}) {
      $unmapped_adaptor->store(Bio::EnsEMBL::UnmappedObject->new(
       -type => 'cDNA',
       -identifier => $id,
       -summary => 'See kill-list database',
       -full_desc => 'This sequence has been excluded from the analysis - see the kill-list database for further details',
       -analysis => $analysis,
      ));
    }
  }
}

1;
