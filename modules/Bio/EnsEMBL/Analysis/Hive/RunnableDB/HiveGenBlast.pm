# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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


# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::GenBlast
#
# Copyright (c) 2009 WormBase
#

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::GenBlast

=head1 SYNOPSIS

  my $blat = Bio::EnsEMBL::Analysis::RunnableDB::GenBlast->
  new(
      -input_id => 'file_name',
      -db => $db,
      -analysis => $analysis,
     );
  $blat->fetch_input;
  $blat->run;
  $blat->write_output;

=head1 DESCRIPTION

This module provides an interface between the ensembl database and
the Runnable GenBlast which wraps the program GenBlast

This module can fetch appropriate input from the database
pass it to the runnable then write the results back to the database
in the dna_align_feature  tables

=head1 CONTACT

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveGenBlast;

use strict;
use warnings;
use feature 'say';
use Data::Dumper;

use Bio::EnsEMBL::Analysis::Runnable::GenBlastGene;
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
  my $genome_slices = $self->get_genome_slices;
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
  my $transcript_biotype = $self->transcript_biotype();

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::GenBlastGene->new
    (
     -query => $query_file,
     -program => $self->analysis->program_file,
     -analysis => $self->analysis,
     -database => $self->analysis->db_file,
     -refslices => $genome_slices,
     -genblast_program => $genblast_program,
     -transcript_biotype => $transcript_biotype,
     %parameters,
    );
  $self->runnable($runnable);

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
  my @output = @{$self->output};
  my $ff = $self->feature_factory;

  foreach my $transcript (@output){
    my $gene = Bio::EnsEMBL::Gene->new();
    $gene->analysis($self->analysis);
    $gene->biotype($self->analysis->logic_name);

    $transcript->analysis($self->analysis);
#    $transcript->slice($self->query) if(!$transcript->slice);
#    if ($self->test_translates()) {
#      print "The GenBlast transcript ",$pt->display_label," doesn't translate correctly\n";
#      next;
#    } # test to see if this transcript translates OK
    $gene->add_Transcript($transcript);
    $adaptor->store($gene);
  }

  if($self->files_to_delete()) {
    my $files_to_delete = $self->files_to_delete();
    `rm $files_to_delete`;
    `rm ${files_to_delete}_*`;
  }

  return 1;
}


sub get_genome_slices {
  my ($self) = @_;
  my @slice_array;
  my $genomic_slices;
  my $slice_adaptor = $self->hrdb_get_con('target_db')->get_SliceAdaptor;
  # This was taken from exonerate module as it is much faster for loading the slices
  # when there are large numbers of slices involved

#
# NOTE: should re-implement the code I commented out below
#

  #also fetching non-reference regions like DR52 for human by default.
  #specify in Exonerate2Genes config-file.
#  if(defined($self->NONREF_REGIONS)){
#    @slice_array = @{$slice_adaptor->fetch_all('toplevel', undef, 1)};
#  }
#  else{
  @slice_array = @{$slice_adaptor->fetch_all('toplevel')};
#  }

  foreach my $slice (@slice_array) {
    $genomic_slices->{$slice->name} = $slice;
  }

  return $genomic_slices;
}


sub output_query_file {
  my ($self) = @_;

  my $accession = $self->param('iid');

  my $table_adaptor = $self->db->get_NakedTableAdaptor();

  $table_adaptor->table_name('uniprot_sequences');

  my $db_row = $table_adaptor->fetch_by_dbID($accession);
  unless($db_row) {
    $self->throw("Did not find an entry int eh uniprot_sequences table matching the accession. Accession:\n".$accession);
  }

  my $seq = $db_row->{'seq'};
  my $transcript_biotype = $db_row->{'biotype'};

  my $record = ">".$accession."\n".$seq;
  my $output_dir = $self->param('query_seq_dir');
  unless(-e $output_dir) {
    `mkdir $output_dir`;
  }

  my $outfile_name = "genblast_".$accession.".fasta";
  my $outfile_path = $output_dir."/".$outfile_name;
  if(-e $outfile_path) {
    $self->warning("Found the query file in the query dir already. Overwriting. File path\n:".$outfile_path);
  }
  open(QUERY_OUT,">".$outfile_path);
  say QUERY_OUT $record;
  close QUERY_OUT;

  $self->files_to_delete($outfile_path);
  $self->transcript_biotype($transcript_biotype);

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


sub transcript_biotype {
  my ($self,$val) = @_;
  if($val) {
    $self->param('_transcript_biotype',$val);
  }

  return($self->param('_transcript_biotype'));
}


sub files_to_delete {
  my ($self,$val) = @_;
  if($val) {
    $self->param('_files_to_delete',$val);
  }

  return($self->param('_files_to_delete'));
}

1;
