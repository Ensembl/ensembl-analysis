# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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


# Ensembl module for Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast

=head1 DESCRIPTION

This module provides an interface between the ensembl database and
the Runnable GenBlast which wraps the program GenBlast

=head1 CONTACT

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast;

use strict;
use warnings;
use feature 'say';
use Data::Dumper;

use Bio::Seq;

use Bio::EnsEMBL::Utils::IO::FASTASerializer;
use Bio::EnsEMBL::Analysis::Runnable::GenBlastGene;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(empty_Transcript);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 param_defaults

 Arg [1]    : None
 Description: It allows the definition of default parameters for all inherting module.
              These are the default values:
               projection_padding => 50000,
 Returntype : Hashref, containing all default parameters
 Exceptions : None

=cut

sub param_defaults {
  my ($self) = @_;
  return {
    %{$self->SUPER::param_defaults},
    projection_padding => 50000,
    timer => '2h',
  }
}


=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::GenBlast
  Function  : fetch data out of database and create runnable
  Returntype: 1
  Exceptions: none
  Example   :

=cut



sub fetch_input {
  my ($self) = @_;

  my $iid_type = $self->param_required('iid_type');
  my $genome_file;
  my $query;

  my $dba = $self->hrdb_get_dba($self->param_required('target_db'));
  my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
  if($dna_dba) {
    $dba->dnadb($dna_dba);
  }

  $self->hrdb_set_con($dba,'target_db');

  $self->create_analysis;
  $self->analysis->program_file($self->param_required('genblast_path'));
  $self->analysis->parameters($self->param('commandline_params'));
  unless($iid_type eq 'projection_transcript_id') {
    $genome_file = $self->param_required('genblast_db_path');
    $self->analysis->db_file($genome_file);
  }

  my %parameters;
  if($self->parameters_hash){
    %parameters = %{$self->parameters_hash};
  }

  if($iid_type eq 'db_seq') {
    $query = $self->get_query_seqs($self->input_id);
  } elsif($iid_type eq 'chunk_file') {
    $query = $self->input_id;
    if($self->param_is_defined('query_seq_dir')) {
      $query = catfile($self->param('query_seq_dir'), $query);
    }
  } elsif($iid_type eq 'projection_transcript_id') {
    my @iid = @{$self->input_id};
    my $projection_dba = $self->hrdb_get_dba($self->param('projection_db'));
    if($dna_dba) {
      $projection_dba->dnadb($dna_dba);
    }
    my $projection_transcript_id = $iid[0];
    my $projection_protein_accession = $iid[1];

    $query = $self->get_query_seqs([$projection_protein_accession]);
    $genome_file = $self->output_transcript_slice($projection_dba, $projection_transcript_id);
  } else {
    $self->throw("You provided an input id type that was not recoginised via the 'iid_type' param. Type provided:\n".$iid_type);
  }

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::GenBlastGene->new
    (
     -program => $self->analysis->program_file,
     -analysis => $self->analysis,
     -database => $genome_file,
     -max_rank => $self->param('max_rank'),
     -genblast_pid => $self->param('genblast_pid'),
     -workdir => File::Temp->newdir(),
     -database_adaptor => $dba,
     %parameters,
    );
  if (ref($query) eq 'ARRAY') {
    $runnable->write_seq_file($query);
  }
  else {
    $runnable->queryfile($query);
  }
  $runnable->genblast_program($self->param('genblast_program')) if ($self->param_is_defined('genblast_program'));
  $runnable->timer($self->param('timer'));
  $self->runnable($runnable);

  return 1;
}


sub run {
  my ($self) = @_;
  my $runnable = shift(@{$self->runnable});
  $self->runnable_failed(0);
  eval {
    $runnable->run;
  };

  if($@) {
    my $except = $@;
    $self->runnable_failed(1);
    if($except =~ /still running after your timer/) {
      $self->warning("genblast took longer than the timer limit (".$self->param('timer')."), will dataflow input id on branch -2. Exception:\n".$except);
      $self->param('_branch_to_flow_to_on_fail',-2);
    } else {
      $self->warning("Issue with running genblast, will dataflow input id on branch -3. Exception:\n".$except);
      $self->param('_branch_to_flow_to_on_fail',-3);
    }
  } else {
    $self->output($runnable->output);
  }

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

  my $adaptor = $self->hrdb_get_con('target_db')->get_GeneAdaptor;
  my $slice_adaptor = $self->hrdb_get_con('target_db')->get_SliceAdaptor;

  if($self->runnable_failed == 1) {
    # Flow out on -2 or -3 based on how the failure happened
    my $failure_branch_code = $self->param('_branch_to_flow_to_on_fail');
    my $output_hash = {};
    $output_hash->{'iid'} = $self->param('iid');
    $self->dataflow_output_id($output_hash,$failure_branch_code);
  } else {
    # Store results, first flag models that are considerably worse than the top hit for each protein
    my $analysis = $self->analysis;
    my $not_best_analysis = new Bio::EnsEMBL::Analysis(-logic_name => $analysis->logic_name()."_not_best",
                                                       -module     => $analysis->module);
    foreach my $transcript ($self->output) {
      $transcript->analysis($analysis);
      my $accession = $transcript->{'accession'};
      $transcript->biotype($self->get_biotype->{$accession});

      if($transcript->{'rank'} > 1) {
        $transcript->analysis($not_best_analysis);
      }
      else {
        $transcript->analysis($analysis);
      }

      if($self->param('flag_small_introns') && $self->avg_intron_size($transcript)) {
        $transcript->biotype($transcript->biotype()."_int");
      }

      empty_Transcript($transcript);
      my $gene = Bio::EnsEMBL::Gene->new();
      $gene->analysis($transcript->analysis);

      $gene->biotype($gene->analysis->logic_name);
      $gene->add_Transcript($transcript);
      $adaptor->store($gene);
    }

    # Now label anything stored that is subpar in terms of score compared to the top ranked hit
    if($self->param('flag_subpar_models')) {
      $self->flag_subpar_models($self->output);
    }
  }

  return 1;
}


sub runnable_failed {
  my ($self,$runnable_failed) = @_;
  if (defined $runnable_failed) {
    $self->param('_runnable_failed',$runnable_failed);
  }
  return ($self->param('_runnable_failed'));
}


sub get_biotype {
  my ($self,$biotype_hash) = @_;
  if($biotype_hash) {
    $self->param('_biotype_hash',$biotype_hash);
  }
  return($self->param('_biotype_hash'));
}


sub get_query_seqs {
  my ($self,$accession_array) = @_;


  my $table_adaptor = $self->db->get_NakedTableAdaptor();
  $table_adaptor->table_name($self->param_required('sequence_table_name'));

  my $biotypes_hash = {};
  my @seqs;

  foreach my $accession (@{$accession_array}) {
    my $db_row = $table_adaptor->fetch_by_dbID($accession);
    unless($db_row) {
      $self->throw("Did not find an entry in the uniprot_sequences table matching the accession. Accession:\n".$accession);
    }

    push(@seqs, Bio::Seq->new(-id => $accession, -seq => $db_row->{seq}));
    $biotypes_hash->{$accession} = $db_row->{'biotype'};
  }

  $self->get_biotype($biotypes_hash);

  return \@seqs;
}


sub output_transcript_slice {
  my ($self,$dba,$transcript_id) = @_;

  my $padding = $self->param('projection_padding');
  my $transcript_adaptor = $dba->get_TranscriptAdaptor();
  my $transcript = $transcript_adaptor->fetch_by_dbID($transcript_id);
  my $slice = $transcript->slice();
  my $new_start = $transcript->seq_region_start() - $padding;
  if($new_start < 1) {
    $new_start = 1;
  }
  my $new_end = $transcript->seq_region_end() + $padding;
  if($new_end > $slice->seq_region_length) {
    $new_end = $slice->seq_region_length;
  }

  my $outfile_path = $self->create_target_file('genblast_slice', 'db');
  my $writer = Bio::EnsEMBL::Utils::IO::FASTASerializer->new($self->param('target_file'), sub {return shift->name});
  $writer->print_Seq($slice->sub_Slice($new_start, $new_end, 1));
  return($outfile_path);
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



sub flag_subpar_models {
  my ($self,$transcripts) = @_;

  my $transcript_adaptor = $self->hrdb_get_con('target_db')->get_TranscriptAdaptor;
  my $uniprot_hitname_hash = {};

  # Store all the combined coverage and pid scores under hitname, then transcript dbID
  # This gives a hash with all the transcript ids for a hit name and their combined scores
  foreach my $transcript (@{$transcripts}) {
    my $sf = shift(@{$transcript->get_all_supporting_features});
    my $hitname = $sf->hseqname;
    $uniprot_hitname_hash->{$hitname}->{$transcript->dbID} = ($sf->percent_id + $sf->hcoverage);
  }

  # Work out the top scoring hit for each hitname and then mark anything under that for relabelling
  my @ids_to_relabel = ();
  foreach my $uniprot_protein (keys(%{$uniprot_hitname_hash})) {
    my $transcript_ids = $uniprot_hitname_hash->{$uniprot_protein};
    my $top_score = 0;
    foreach my $transcript_id (keys(%{$transcript_ids})) {
      my $combined_score = $transcript_ids->{$transcript_id};
      if($combined_score > $top_score) {
        $top_score = $combined_score;
      }
    }

    foreach my $transcript_id (keys(%{$transcript_ids})) {
      my $combined_score = $transcript_ids->{$transcript_id};
      unless($combined_score >= ($top_score - (($top_score * 5) / 100))) {
        push(@ids_to_relabel,$transcript_id);
      }
    }
  }

  # Fetch the transcript and relabel
  foreach my $id_to_relabel (@ids_to_relabel) {
    my $transcript = $transcript_adaptor->fetch_by_dbID($id_to_relabel);
    $transcript->biotype($transcript->biotype.'_sub');
    $transcript_adaptor->update($transcript);
  }
}


sub avg_intron_size {
  my ($self,$transcript) = @_;

  my $short_length = 75;
  my $max_short = 2;
  my $introns = $transcript->get_all_Introns();
  my $total_introns = scalar(@$introns);
  my $total_length = 0;
  my $total_short = 0;
  foreach my $intron (@$introns) {
    my $length = $intron->length();
    $total_length += $length;
    if($length < $short_length) {
      $total_short++;
    }
  }

  if($total_introns) {
    my $avg_length = $total_length / $total_introns;
    if($avg_length < $short_length) {
      say "Transcript avg intron length too short: ".$transcript->dbID." ".$avg_length." ".$total_short;
      return 1;
    } elsif($total_short > $max_short) {
      say "Too many short introns: ".$transcript->dbID." ".$avg_length." ".$total_short;
      return 1;
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}


1;
