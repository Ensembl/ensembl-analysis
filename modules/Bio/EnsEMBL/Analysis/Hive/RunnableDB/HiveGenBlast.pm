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

use Bio::EnsEMBL::Analysis::Runnable::GenBlast;
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
  my $time = time();
  say "FM2 time before open: ".$time;
#  open(my $fh, $genome_file) or $self->throw("Could not open $genome_file for reading");
#  while(<$fh>) {
#    /^\>(\S+)/ and do {
#      my $seq_name = $1;
#      say "FM2 DUMPER: ".Dumper($dba);
#       my $slice = $dba->get_SliceAdaptor->fetch_by_name($seq_name);
#        my $slice = $self->fetch_sequence($seq_name,$dba);
#       $self->query($slice);
#        $time = time();
#       say "FM2 slice time ".$slice->name.": ".$time;
#      my $slice = $dba->get_SliceAdaptor->fetch_by_region('toplevel', $seq_name);
#      say "FM2 SLICE: ".Dumper($slice);
#      if (not defined $slice) {
#        $self->throw("Could not extract slice for $seq_name from database");
#      }
#      $genome_slices{$seq_name} = $slice;
#    }
#  }
#  close $fh;

  $time = time();
  say "FM2 time after open: ".$time;
  my $query_file = $self->input_id;
  if($self->param('query_seq_dir')) {
    $query_file = $self->param('query_seq_dir')."/".$query_file;
  }

  say "FM2 pre runnable: ";
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::GenBlast->new
    (
     -query => $query_file,
     -program => $self->analysis->program_file,
     -analysis => $self->analysis,
     -database => $self->analysis->db_file,
     -refslices => $genome_slices,
     %parameters,
    );
  $self->runnable($runnable);

#  say "FM2 DUMP DB SLICES: ".Dumper($genome_slices);
#  foreach my $slice_names (keys(%{$genome_slices})) {
#    say "FM2 SLICE NAME: ".$slice_names;
#  }
#  exit;

  $time = time();
  say "FM2 leaving fetch: ".$time;
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
  say "FM2 entering WO: ";
  my $adaptor = $self->hrdb_get_con('target_db')->get_GeneAdaptor;
  my @output = @{$self->output};
  my $ff = $self->feature_factory;

  foreach my $transcript (@output){
    my $gene = Bio::EnsEMBL::Gene->new();
    $gene->analysis($self->analysis);
    $gene->biotype($self->analysis->logic_name);

    $transcript->analysis($self->analysis);
    $transcript->slice($self->query) if(!$transcript->slice);
#    if ($self->test_translates()) {
#      print "The GenBlast transcript ",$pt->display_label," doesn't translate correctly\n";
#      next;
#    } # test to see if this transcript translates OK
    $gene->add_Transcript($transcript);
    $adaptor->store($gene);
  }
  say "FM2 leaving WO: ";
  return 1;
}


sub get_genome_slices {
  my ($self) = @_;
  my @slice_array;
  my $genomic_slices;
  my $slice_adaptor = $self->hrdb_get_con('target_db')->get_SliceAdaptor;

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


1;
