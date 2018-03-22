#!/usr/bin/env perl

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

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveComparecDNAfiles;

use strict;
use warnings;
use feature 'say';


#use Bio::EnsEMBL::IO::Parser::Fasta;
#use Bio::Phylo::Parsers::Fasta;
use Bio::SeqIO;
use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');



sub fetch_input {
  my $self = shift;
  return 1;
}

sub run {
  my ($self) = shift;

  my $query_hash = $self->param('compare_files');
  my $output_path = $query_hash->{'dest_dir'};
  my $new_cdna_file = $query_hash->{'new_cdna_file'};
  my $old_cdna_file = $query_hash->{'old_cdna_file'};
  my $filtered_cdnas = $query_hash->{'filtered_file'};
  my $retired_cdnas = $query_hash->{'retired_file'};

  $self->compare_cdnas($new_cdna_file,$old_cdna_file,$output_path,$filtered_cdnas,$retired_cdnas);

  return 1;
}

sub write_output {
  my $self = shift;
  return 1;
}

sub compare_cdnas {
  my ($self,$new_file,$old_file,$output_path,$filtered_cdnas,$retired_cdnas) = @_;

  my $start_index = 0;

  #my $new_parser = Bio::EnsEMBL::IO::Parser::Fasta->open($new_file);
  my $new_parser = Bio::Phylo::Parsers::Fasta->open($new_file);
  my $new_header;
  my $new_seq;

  #my $old_parser = Bio::EnsEMBL::IO::Parser::Fasta->open($old_file);
  my $old_parser = Bio::Phylo::Parsers::Fasta->open($old_file);
  my $old_header;
  my $old_seq;

  my %new_cdnas;
  my %old_cdnas;

  open( FC, ">", $filtered_cdnas ) or die "can't create file $filtered_cdnas\n";
  open( RETIRED, ">", $retired_cdnas ) or die "can't create file $retired_cdnas\n";
  
  # read the new file
  while($new_parser->next()) {
    $new_header = $new_parser->getHeader();
    $new_seq = $new_parser->getSequence();

    $new_cdnas{$new_header} = $new_seq;
  }  

  # read old file
  while($old_parser->next()) {
    $old_header = $old_parser->getHeader();
    $old_seq = $old_parser->getSequence();

    $old_cdnas{$old_header} = $old_seq;
  }

  # compare the 2 hashes to see which should be copied over and which shouldn't
  # first go through the new cdnas and see if they exist in the previous one
  # if they do we don't need to do anything, if they don't we add the cdna to the 
  # table or file and align using exonerate
  for (keys %new_cdnas) {
    if (exists $old_cdnas{$_}) {
      # perhaps remove the entry from the old_cdnas hash so we don't need to go through it to find which ones have been retired
      delete $old_cdnas{$_};
      next;
    } else {
      # either print to a new file or just do the loading of the cdnas at this point
      print FC ">", $_, "\n", $new_cdnas{$_}, "\n";
      next;
    }
  }
  # go through remaining cdnas in the old hash and delete the entries from the database
  for (keys %old_cdnas) {
    # delete from gene where dna_align_feature hit = '$old_cdnas{$_}';
    print RETIRED $_, "\n";

    #select distinct(gene.gene_id) from gene left join transcript on transcript.gene_id = gene.gene_id left join exon_transcript on exon_transcript.transcript_id = transcript.transcript_id left join supporting_feature on supporting_feature.exon_id = exon_transcript.exon_id left join dna_align_feature on dna_align_feature.dna_align_feature_id = supporting_feature.feature_id where hit_name = '';
    # delete from unmapped_object where identifier = '$old_cdnas{$_}';
    next;
  } 

  #close(RF);
  close(FC);
  close(RETIRED);
}


1;

