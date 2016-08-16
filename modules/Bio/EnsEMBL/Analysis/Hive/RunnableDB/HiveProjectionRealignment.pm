#!/usr/bin/env perl

# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveProjectionRealignment;

use strict;
use warnings;
use feature 'say';
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

sub fetch_input {
  my $self = shift;
  my $projection_dba = $self->hrdb_get_dba($self->param('projection_db'));
  my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
  if($dna_dba) {
    $projection_dba->dnadb($dna_dba);
  }

  my $projection_source_dba = $self->hrdb_get_dba($self->param('projection_source_db'));

  my $output_dba = $self->hrdb_get_dba($self->param('projection_realign_db'));
  if($dna_dba) {
    $output_dba->dnadb($dna_dba);
  }

  $self->hrdb_set_con($output_dba,'output_db');

  my $projection_transcript_id = $self->param("iid");
  my $transcript_adaptor = $projection_dba->get_TranscriptAdaptor();
  say "D1";
  my $projection_transcript = $transcript_adaptor->fetch_by_dbID($projection_transcript_id);
  say "D2";
  $self->projection_transcript($projection_transcript);
  say "D3 " .$projection_transcript->stable_id();

  my $source_transcript_adaptor = $projection_source_dba->get_TranscriptAdaptor();
  say "D4";
  my $parent_transcript = $source_transcript_adaptor->fetch_by_stable_id($projection_transcript->stable_id);
  unless($parent_transcript) {
    $self->throw("Could not retrieve the parent transcript from the projection source db using the stable id of the projected transcript!\n".
                 "Stable id: ".$projection_transcript->stable_id."\n".
                 "Projection source db: ".$projection_source_dba->dbc->dbname."@".$projection_source_dba->dbc->host);
  }

  say "D5";
  my $parent_protein = $parent_transcript->translation->seq();
  say "D6";
  unless($parent_protein) {
    $self->throw("Issue fetching the translation sequence for parent transcript from projection source db!\n".
                 "Stable id: ".$parent_transcript->stable_id.
                 "Projection source db: ".$projection_source_dba->dbc->dbname."@".$projection_source_dba->dbc->host);
  }

  $self->parent_protein($parent_protein);

  my $projection_protein_table = $self->param('protein_table_name');
  $self->projection_protein_table($projection_protein_table);

  return 1;
}

sub run {
  my $self = shift;
  my $is_bad = $self->assess_projection_transcript();
  $self->do_realignment($is_bad);
  return 1;
}

sub write_output {
  my $self = shift;

  my $transcript = $self->projection_transcript();
  if($self->do_realignment) {
    my $parent_protein = $self->parent_protein();
    my $protein_table_name = $self->projection_protein_table();

    # Store the stable id and the protein in the pipeline protein table for realignment later
    my $table_adaptor = $self->db->get_NakedTableAdaptor();
    $table_adaptor->table_name($protein_table_name);
    my $db_row = [{
                    'accession'  => $transcript->stable_id,
                    'seq'        => $parent_protein,
                 }];
    $table_adaptor->store($db_row);

    # Flow this transcript db id and the stable id out to the downstream analysis
    # At the moment the stable id is techincally unneeded, but we might use this technique elsewhere in future
    my $stable_id = $transcript->stable_id;
    my $db_id = $transcript->dbID;
    my $output_hash = {};
    $output_hash->{'iid'} = [$db_id,$stable_id];
    $self->dataflow_output_id($output_hash,2);

  } else {
    my $output_dba = $self->hrdb_get_con('output_db');
    my $gene = $transcript->get_Gene();
    my $gene_adaptor = $output_dba->get_GeneAdaptor();
    empty_Gene($gene);
    $gene_adaptor->store($gene);
  }

  return 1;
}

sub projection_transcript {
  my ($self,$transcript) = @_;
  if($transcript) {
    $self->param('_projection_transcript',$transcript);
  }

  return($self->param('_projection_transcript'));
}

sub parent_protein {
  my ($self,$protein) = @_;
  if($protein) {
    $self->param('_parent_protein',$protein);
  }

  return($self->param('_parent_protein'));
}

sub projection_protein_table {
  my ($self,$projection_protein_table) = @_;
  if($projection_protein_table) {
    $self->param('_projection_protein_table',$projection_protein_table);
  }

  return($self->param('_projection_protein_table'));
}


sub do_realignment {
  my ($self,$is_bad) = @_;
  if($is_bad) {
    $self->param('_is_bad',$is_bad);
  }

  return($self->param('_is_bad'));
}

sub assess_projection_transcript {
  my ($self) = @_;

  my $transcript = $self->projection_transcript();
  my $has_non_canonical_splice_sites = 0;
  my $has_frameshift_introns = 0;

  my @introns = @{$transcript->get_all_Introns()};
  foreach my $intron (@introns) {
    if($intron->length < 10) {
      $has_frameshift_introns++;
    } elsif($intron->is_splice_canonical() == 0) {
      $has_non_canonical_splice_sites++;
    }
  }

  say "Transcript has ".$has_non_canonical_splice_sites." splice sites";
  say "Transcript has ".$has_frameshift_introns." frameshift introns";

  if($has_non_canonical_splice_sites || $has_frameshift_introns) {
    return(1);
  } else {
    return(0);
  }
}

1;
