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
use Data::Dumper;

use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::Analysis::Runnable::GenBlastGene;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene attach_Slice_to_Gene);
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


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

  my $analysis = Bio::EnsEMBL::Analysis->new(
                                              -logic_name => $self->param('logic_name'),
                                              -module => $self->param('module'),
                                              -program_file => $self->param('genblast_path'),
                                              -db_file => $self->param('genblast_db_path'),
                                              -parameters => $self->param('commandline_params'),
                                            );
  $self->analysis($analysis);

  my %parameters;
  if($self->parameters_hash){
    %parameters = %{$self->parameters_hash};
  }

  my $genome_file = $self->analysis->db_file;

  my $query_file;

  my $iid_type = $self->param('iid_type');
  unless($iid_type) {
    $self->throw("You haven't provided an input id type. Need to provide one via the 'iid_type' param");
  }

  if($iid_type eq 'db_seq') {
    $query_file = $self->output_query_file();
  } elsif($iid_type eq 'chunk_file') {
    $query_file = $self->param('iid');
    if($self->param('query_seq_dir')) {
      $query_file = $self->param('query_seq_dir')."/".$query_file;
    }
  } else {
    $self->throw("You provided an input id type that was not recoginised via the 'iid_type' param. Type provided:\n".$iid_type);
  }


  my $genblast_program = $self->param('genblast_program');
  my $biotypes_hash = $self->get_biotype();
  my $max_rank = $self->param('max_rank');
  my $genblast_pid = $self->param('genblast_pid');

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::GenBlastGene->new
    (
     -query => $query_file,
     -program => $self->analysis->program_file,
     -analysis => $self->analysis,
     -database => $self->analysis->db_file,
     -genblast_program => $genblast_program,
     -max_rank => $max_rank,
     -genblast_pid => $genblast_pid,
     -work_dir => $self->param('query_seq_dir'),
     -database_adaptor => $dba,
     %parameters,
    );
  $self->runnable($runnable);

  return 1;
}


sub run {
  my ($self) = @_;

  foreach my $runnable (@{$self->runnable}) {
    eval {
      $runnable->run;
    };

    if($@) {
      my $except = $@;
      if($except =~ /Error closing exonerate command/) {
        warn("Error closing exonerate command, this input id was not analysed successfully:\n".$self->input_id);
      } else {
        $self->throw($except);
      }
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

  my @output = @{$self->output};
  my $ff = $self->feature_factory;

  foreach my $transcript (@output){
    my $gene = Bio::EnsEMBL::Gene->new();
    $gene->analysis($self->analysis);

    $transcript->analysis($self->analysis);
    my $accession = $transcript->{'accession'};
    my $transcript_biotype = $self->get_biotype->{$accession};
    $transcript->biotype($transcript_biotype);

    if($transcript->{'rank'} > 1) {
      my $analysis = new Bio::EnsEMBL::Analysis(-logic_name => $transcript->analysis->logic_name()."_not_best",
                                                -module     => $transcript->analysis->module);
      $transcript->analysis($analysis);
      $gene->analysis($analysis);
    }

    $gene->biotype($gene->analysis->logic_name);

#    if ($self->test_translates()) {
#      print "The GenBlast transcript ",$pt->display_label," doesn't translate correctly\n";
#      next;
#    } # test to see if this transcript translates OK
    $gene->add_Transcript($transcript);
    $gene->slice($transcript->slice);
    empty_Gene($gene);
    $adaptor->store($gene);
  }

  if($self->param('store_blast_results')) {
    $self->store_protein_align_features();
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

sub store_protein_align_features {
  my ($self) = @_;

  # Need to build the path to theblast report file. It is the query seq path with the genome db appended
  # on to it (with a .wublast.report extension)
  my $query_file_path = $self->files_to_delete();
  my $genome_db_path = $self->param('genblast_db_path');
  unless($genome_db_path =~ /[^\/]+$/) {
   $self->throw("Could not parse the genome db file name off the end of the path. Path used:\n".$genome_db_path);
  }

  my $genome_db_name = $&;
  my $report_file_path = $query_file_path."_".$genome_db_name.".wublast.report";

  unless(-e $report_file_path) {
    $self->throw("Parsed the report file info, but the constructed path points to a non-existant file. Offending path:\n".$report_file_path);
  }

  my $dba = $self->hrdb_get_con('target_db');
  my $protein_align_feature_adaptor = $dba->get_ProteinAlignFeatureAdaptor();
  my $analysis = $self->analysis();
  if($self->param('protein_align_feature_logic_name')) {
    my $logic_name = $self->param('protein_align_feature_logic_name');
    unless($self->analysis->logic_name eq $logic_name) {
      $analysis = Bio::EnsEMBL::Analysis->new(-logic_name => $logic_name);
    }
  }

  say "Parsing hits from:\n".$report_file_path;
  open(BLAST_REPORT,$report_file_path);
  my $remove_info_line = <BLAST_REPORT>;
  my $current_hit_name;
  while(<BLAST_REPORT>) {
    my $line = $_;
    if($line =~ /^\>(\S+)/) {
      $current_hit_name = $1;
    } else {
      my @hit_array = split('\s',$line);
      unless(scalar(@hit_array) == 14) {
        $self->throw("Did not get the expected number of columns when parsing a hit for ".$current_hit_name.", offending line:\n".$line);
      }

      my $seq_region_start = $hit_array[7];
      my $seq_region_end = $hit_array[8];
      my $seq_region_strand = $hit_array[10];
      my $hit_start = $hit_array[11];
      my $hit_end = $hit_array[12];
      my $hit_score = $hit_array[2];
      my $percent_identity = $hit_array[5];
      my $e_value = $hit_array[4];
      my $cigar = $hit_array[13]."M";
      my $slice = $self->genome_slices->{$hit_array[6]};
      unless($slice) {
        self->throw("Could not find a matching slice for:\n".$hit_array[6]);
      }

      my $hit = Bio::EnsEMBL::DnaPepAlignFeature->new(
                                                      -start      => $seq_region_start,
                                                      -end        => $seq_region_end,
                                                      -strand     => $seq_region_strand,
                                                      -hseqname   => $current_hit_name,
                                                      -hstart     => $hit_start,
                                                      -hend       => $hit_end,
                                                      -score      => $hit_score,
                                                      -percent_id => $percent_identity,
                                                      -cigar_string => $cigar,
                                                      -p_value    => $e_value,
                                                      -slice      => $slice,
                                                      -analysis   => $analysis);

      $protein_align_feature_adaptor->store($hit);
    }

  }
  close BLAST_REPORT;
}


sub get_biotype {
  my ($self,$biotype_hash) = @_;
  if($biotype_hash) {
    $self->param('_biotype_hash',$biotype_hash);
  }
  return($self->param('_biotype_hash'));
}



sub output_query_file {
  my ($self) = @_;

  my $accession_array = $self->param('iid');

  my $table_adaptor = $self->db->get_NakedTableAdaptor();
  $table_adaptor->table_name('uniprot_sequences');


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
  $self->get_biotype($biotypes_hash);

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


sub files_to_delete {
  my ($self,$val) = @_;
  if($val) {
    $self->param('_files_to_delete',$val);
  }

  return($self->param('_files_to_delete'));
}

1;
