# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016] EMBL-European Bioinformatics Institute
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

use Bio::EnsEMBL::DnaPepAlignFeature;
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

  my $dba = $self->hrdb_get_dba($self->param('target_db'));
  my $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
  if($dna_dba) {
    $dba->dnadb($dna_dba);
  }

  $self->hrdb_set_con($dba,'target_db');

  $self->create_analysis;
  my $genome_file = $self->param('genblast_db_path');
  $self->analysis->program_file($self->param('genblast_path'));
  $self->analysis->db_file($genome_file);
  $self->analysis->parameters($self->param('commandline_params'));

  my %parameters;
  if($self->parameters_hash){
    %parameters = %{$self->parameters_hash};
  }


  my $query_file;

  my $iid_type = $self->param_required('iid_type');
  # If there is a sequence table name, store for use later
  $self->param_required('sequence_table_name');

  if($iid_type eq 'db_seq') {
    $query_file = $self->output_query_file($self->param("iid"));
  } elsif($iid_type eq 'chunk_file') {
    $query_file = $self->param('iid');
    if($self->param_is_defined('query_seq_dir')) {
      $query_file = $self->param('query_seq_dir')."/".$query_file;
    }
  } elsif($iid_type eq 'projection_transcript_id') {
    my @iid = @{$self->param("iid")};
    my $projection_dba = $self->hrdb_get_dba($self->param('projection_db'));
    if($dna_dba) {
      $projection_dba->dnadb($dna_dba);
    }
    my $projection_transcript_id = $iid[0];
    my $projection_protein_accession = $iid[1];
    my $padding = $self->param('projection_padding');

    $query_file = $self->output_query_file([$projection_protein_accession]);
    $genome_file = $self->output_transcript_slice($projection_dba,$projection_transcript_id,$padding);
  } else {
    $self->throw("You provided an input id type that was not recoginised via the 'iid_type' param. Type provided:\n".$iid_type);
  }

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::GenBlastGene->new
    (
     -query => $query_file,
     -program => $self->analysis->program_file,
     -analysis => $self->analysis,
     -database => $genome_file,
     -max_rank => $self->param('max_rank'),
     -genblast_pid => $self->param('genblast_pid'),
     -workdir => $self->param('query_seq_dir'),
     -database_adaptor => $dba,
     %parameters,
    );
  $runnable->genblast_program($self->param('genblast_program')) if ($self->param_is_defined('genblast_program'));
  $self->runnable($runnable);

  return 1;
}


sub run {
  my ($self) = @_;
  $self->runnable_failed(0);
  foreach my $runnable (@{$self->runnable}) {
    eval {
      $runnable->run;
    };

    if($@) {
      my $except = $@;
      $self->runnable_failed(1);
      $self->warning("Issue with running genblast, will dataflow input id on branch -3. Exception:\n".$except);
    } else {
      $self->output($runnable->output);
    }
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
  my $failure_branch_code = -3;

  if($self->runnable_failed == 1) {
    my $output_hash = {};
    $output_hash->{'iid'} = $self->param('iid');
    $self->dataflow_output_id($output_hash,$failure_branch_code);
  } else {
    my @output = @{$self->output};

    my $analysis = $self->analysis;
    my $not_best_analysis = new Bio::EnsEMBL::Analysis(-logic_name => $analysis->logic_name()."_not_best",
                                                       -module     => $analysis->module);
    foreach my $transcript (@output){
      $transcript->analysis($analysis);
      my $accession = $transcript->{'accession'};
      $transcript->biotype($self->get_biotype->{$accession});

      if($transcript->{'rank'} > 1) {
        $transcript->analysis($not_best_analysis);
      }
      else {
        $transcript->analysis($analysis);
      }

      empty_Transcript($transcript);
      my $gene = Bio::EnsEMBL::Gene->new();
      $gene->analysis($transcript->analysis);

      $gene->biotype($gene->analysis->logic_name);
      $gene->add_Transcript($transcript);
      $adaptor->store($gene);
    }

  }

  if($self->files_to_delete()) {
    my $files_to_delete = $self->files_to_delete();
    `rm -r $files_to_delete/genblast_*`;
    `rmdir $files_to_delete`;
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


sub output_query_file {
  my ($self,$accession_array) = @_;


  my $table_adaptor = $self->db->get_NakedTableAdaptor();
  $table_adaptor->table_name($self->param('sequence_table_name'));


  my $rand = int(rand(1000000));
  # Note as each accession will occur in only one file, there should be no problem using the first one
  my $outfile_name = "genblast_".$$."_".$rand.".fasta";
  my $output_dir = $self->param('query_seq_dir')."/genblast_".$$."_".$rand;
  my $outfile_path = $output_dir."/".$outfile_name;

  my $biotypes_hash = {};


  unless(-e $output_dir) {
    `mkdir -p $output_dir`;
  }

  if(-e $outfile_path) {
    $self->warning("Found the query file in the query dir already. Overwriting. File path\n:".$outfile_path);
  }

  open(QUERY_OUT,">".$outfile_path);
  foreach my $accession (@{$accession_array}) {
    my $db_row = $table_adaptor->fetch_by_dbID($accession);
    unless($db_row) {
      $self->throw("Did not find an entry in the uniprot_sequences table matching the accession. Accession:\n".$accession);
    }

    my $seq = $db_row->{'seq'};
    $biotypes_hash->{$accession} = $db_row->{'biotype'};

    my $record = ">".$accession."\n".$seq;
    say QUERY_OUT $record;
  }
  close QUERY_OUT;

  $self->files_to_delete($output_dir);
  $self->output_dir_path($output_dir);

  $self->get_biotype($biotypes_hash);

  return($outfile_path);
}


sub output_transcript_slice {
  my ($self,$dba,$transcript_id,$padding) = @_;

  my $transcript_adaptor = $dba->get_TranscriptAdaptor();
  my $transcript = $transcript_adaptor->fetch_by_dbID($transcript_id);
  my $slice = $transcript->slice();
  my $new_start = $transcript->seq_region_start() - $padding;
  if($new_start <= 0) {
    $new_start = 1;
  }
  my $new_end = $transcript->seq_region_end() + $padding;
  if($new_end > $slice->seq_region_length) {
    $new_end = $slice->seq_region_length;
  }

  my $transcript_slice  = Bio::EnsEMBL::Slice->new
      (-seq_region_name   => $slice->seq_region_name,
       -seq_region_length => $slice->seq_region_length,
       -coord_system      => $slice->coord_system,
       -start             => $new_start,
       -end               => $new_end,
       -strand            => 1,
       -adaptor           => $slice->adaptor());

  my $name = $transcript_slice->name;
  my $seq = $transcript_slice->seq;

  my $outfile_name = "genblast_slice.db";
  my $output_dir = $self->output_dir_path;
  my $outfile_path = $output_dir."/".$outfile_name;

  unless(-e $output_dir) {
    `mkdir -p $output_dir`;
  }

  if(-e $outfile_path) {
    $self->warning("Found the db file in the query dir already. Overwriting. File path\n:".$outfile_path);
  }

  open(SLICE_OUT,">".$outfile_path);
  my $record = ">".$name."\n".$seq;
  say SLICE_OUT $record;
  close SLICE_OUT;

#  $self->files_to_delete($output_dir);

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


sub output_dir_path {
  my ($self,$val) = @_;
  if($val) {
    $self->param('_output_dir_path',$val);
  }

  return($self->param('_output_dir_path'));
}


sub files_to_delete {
  my ($self,$val) = @_;
  if($val) {
    $self->param('_files_to_delete',$val);
  }

  return($self->param('_files_to_delete'));
}

1;
