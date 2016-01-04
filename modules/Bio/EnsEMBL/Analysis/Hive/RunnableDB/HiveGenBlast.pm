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

use Bio::EnsEMBL::DnaPepAlignFeature;
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
  $self->genome_slices($genome_slices);

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

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::GenBlastGene->new
    (
     -query => $query_file,
     -program => $self->analysis->program_file,
     -analysis => $self->analysis,
     -database => $self->analysis->db_file,
     -refslices => $genome_slices,
     -genblast_program => $genblast_program,
     -biotypes => $biotypes_hash,
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
    my $accession = $transcript->{'accession'};
    my $transcript_biotype = $self->get_biotype->{$accession};
    $transcript->biotype($transcript_biotype);

    $self->get_supporting_features($transcript);

#    $transcript->slice($self->query) if(!$transcript->slice);
#    if ($self->test_translates()) {
#      print "The GenBlast transcript ",$pt->display_label," doesn't translate correctly\n";
#      next;
#    } # test to see if this transcript translates OK
    $gene->add_Transcript($transcript);
    $adaptor->store($gene);
  }

  if($self->param('store_blast_results')) {
    $self->store_protein_align_features();
  }

  if($self->files_to_delete()) {
    my $files_to_delete = $self->files_to_delete();
    `rm $files_to_delete`;
    `rm ${files_to_delete}_*`;
  }

  return 1;
}

sub get_supporting_features {
  my ($self,$transcript) = @_;

  my $accession = $transcript->{'accession'};
  my $genblast_id = $transcript->{'genblast_id'};
  my $seq_region_start = $transcript->seq_region_start;
  my $seq_region_end = $transcript->seq_region_end;
  my $strand = $transcript->seq_region_strand;
  my $seq_region_name = $transcript->seq_region_name;
  my $genblast_report_path = $self->files_to_delete()."_1.1c_2.3_s1_0_16_1";

  my @report_array = ();
  my $query_seq;
  my $target_seq;
  my $gene_info;
  my $transcript_percent_id;
  my $transcript_coverage;

  my $found = 0;
  open(GENBLAST_REPORT,$genblast_report_path);
  while(<GENBLAST_REPORT>) {
    my $line = $_;
    chomp $line;

    push(@report_array,$line);
    if(scalar(@report_array) > 5) {
      shift(@report_array);
    }

    if($line =~ /^Gene\:ID\=$genblast_id\|/) {
      $found = 1;
      $query_seq = $report_array[0];
      $query_seq =~ s/^query\://;
      $target_seq = $report_array[2];
      $target_seq =~ s/^targt\://;
      $gene_info = $report_array[4];
      $gene_info =~ /\|PID\:([^\|]+)$/;
      $transcript_percent_id = $1;
      last;
    }
  }
  close(GENBLAST_REPORT);

  unless($found) {
     $self->throw("Issue with parsing the genblast report, transcript id not found in report file. Transcript id:\n".$genblast_id.
                  "\nFile path:\n".$genblast_report_path);
  }

  unless($transcript_percent_id >= 0) {
    $self->throw("Issue with parsing the percent identity of the transcript. Transcript id in genblast file:\n".$genblast_id.
                 "\nFile path:\n".$genblast_report_path."\nOffending info line:\n".$gene_info);
  }

  $transcript_coverage = $self->calculate_coverage($query_seq,$target_seq);

  unless($transcript_coverage >= 0) {
    $self->throw("Issue with parsing the query coverage of the transcript. Transcript id in genblast file:\n".$genblast_id.
                 "\nFile path:\n".$genblast_report_path."\nOffending info line:\n".$gene_info);
  }

  say "Transcript coverage: ".$transcript_coverage;
  say "Transcript percent identity: ".$transcript_percent_id;


  my $exon_supporting_features = $self->get_exon_supporting_features($transcript,$query_seq,$target_seq,$transcript_coverage,$transcript_percent_id);
  my $exons = $transcript->get_all_Exons;
  $transcript->flush_Exons;

  for(my $i=0; $i<scalar(@{$exons}); $i++) {
    if(${$exon_supporting_features}[$i]) {
      my $exon =  ${$exons}[$i];
      my $exon_supporting_features = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => [${$exon_supporting_features}[$i]]);
      $exon->add_supporting_features($exon_supporting_features);
      $transcript->add_Exon($exon);
    }
  }

  my $transcript_supporting_features = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => $exon_supporting_features);
  $transcript->add_supporting_features($transcript_supporting_features);

}

sub get_exon_supporting_features {
  my ($self,$transcript,$query_seq,$target_seq,$cov,$pid) = @_;

  # The API/schema handles protein align features slightly differently than expected
  # The hit start and end do not refer to the corresponding positions in query sequence (the aligned protein) covered by
  # the alignment, instead they refer to the positions in the ungapped translation sequence of the transcript the query
  # was aligned to. This is odd for several reasons. But the main point is that there is a check in the API that expects
  # the feature length (number of bases the feature spans on the exon) to be exactly three times the lenght of the hit
  # length (the number of residues in the translation of the transcript covered by the alignment). What this means is
  # that there is no way to repesent split codons (codons than span two exons), because they would throw this 3 to 1 ratio
  # off and the API will throw. The code below is written with this limitation in mind
  # There is also the issue of coverage and percent id this throws up. For now the most reasonable was to present coverage
  # and percent id is on the transcript level. I have put code below to allow for the calculation on the exon level, but it
  # is currently disabled and the transcript values are used instead for every exon

  # This is best explained by a diagram:

  #  1   2   3   4   5   6   7   8   9   10  11  12  13    <- alignment col
  # .M. .P. .Q. .R. .S. .M. .P. .P. .-. .-. .S. .R. .S.    <- protein seq
  # .M. .-. .-. .R. .S. .M. .P. .P. .Q. .Q. .S. .R. .S.    <- transcript translation seq
  #  1           2   3   4   5   6   7   8   9   10  11    <- translation ungapped index
  # ATG         AGA AGT ATG CCT CCT CAA CAA AGT AGA AGT    <- codons
  # ###         ### #\\ \\\ \\\ \\# ### ### ### ### ###    <- alternating exons defined by # and \
  # 1                8            18                  33   <- genomic coordinates

  # For PAF what are the start, end, hit start and hit end coords for the corresponding DnaPepAlignFeature?
  # Exon 1:
  # Codon 1 and 2 are complete, codon 3 is split so have to discard because of how the API and the protein_align_feature table work
  # Feature start = 1
  # Feature end   = 6  ---> 6 - 1 + 1 = 6
  # Hit start     = 1
  # Hit end       = 2  ---> 2 - 1 + 1 = 2
  # Ratio         = 6:2 = 3:1
  #
  # Exon 2:
  # Covers codons 3 to 6. Codons 3 and 6 are split codons, so the feature will cover 4 and 5
  # Feature start = 10
  # Feature end   = 15 ---> 15 -10 + 1 = 6
  # Hit start     = 4
  # Hit end       = 5  ---> 5 - 4 + 1 = 2
  # Ratio         = 6:2 = 3:1
  #
  # Exon 3:
  # Covers codons 6 to 11. Codon 6 is split, so feature will cover codon 7 to 11
  # Feature start = 19
  # Feature end   = 33 ---> 33 - 19 + 1 = 15
  # Hit start     = 7
  # Hit end       = 11 ---> 11 - 7 + 1 = 5
  # Ratio         = 15:5 = 3:1

  my $exon_supporting_features;
  my $query_sub_align;
  my $target_sub_align;

  my $exons = $transcript->get_all_Exons();
  my $pep_index = 0;
  my $i=0;
  for($i=0; $i<scalar(@{$exons}); $i++) {
    my $exon = ${$exons}[$i];

    # If there is a split codon at the start of the exon we have to skip it as the API currently can't handle the idea of split
    # codons in protein align features, the feature length must be exactly three times the length of the hit itself
    # If an exon is phase 1 it means there are two bases in the split codon, if it is phase 2 then there is one base
    # Also worth remembering that the end phase of an exon is always the same as the phase of the new one
    # If there is a split codon then the index of the translation residue needs to be incremented to skip the residue over the
    # split codon
    my $bases_to_skip = 0;
    if($exon->phase) {
      if($exon->phase == 1) {
        $bases_to_skip = 2;
      } else {
        $bases_to_skip = 1;
      }
      $pep_index++;
    }

    # Make the feature start index of the first base of the first complete codon
    my $feature_start = $exon->start() + $bases_to_skip;

    # To calculate the feature end, use the end phase to work out if there is a split codon. The end phase is equal to the number of
    # base pairs in the split codon at the end of the exon (if present)
    my $feature_end;
    my $bases_to_remove = $exon->end_phase;

    # If it's the last exon and the end phase is 0 and the exon ends in a stop codon then remove this when building the feature coords
    if($i == (scalar(@{$exons} - 1)) && $exon->end_phase == 0 && ($exon->seq->seq =~ /T[AG]A$/ || $exon->seq->seq =~ /TAG$/)) {
      $bases_to_remove = 3;
    }

    # Remove any incomplete or stop codons from the feature end coords
    $feature_end = $exon->end() - $bases_to_remove;

    # This is a check to basically make sure that the feature start and end now cover only complete codons. The reason it is here is
    # that the core API will throw later if this isn't true but the message is a bit cryptic
    if(($feature_end - $feature_start + 1) % 3) {
      $self->throw("\nThe feature coords after removal of split codons were not divisable by 3. They should represent complete codons.\n".
                   "Feature start: ".$feature_start."\nFeature end: ".$feature_end."\nExon start: ".$exon->start."\nExon end: ".$exon->end.
                   "\nExon phase: ".$exon->phase."\nExon end phase: ".$exon->end_phase."\n");
    }

    # Now need to calculate the offset that the feature covers in the ungapped translation seq. It is important to remember
    # that this is not identical to the equivalent column in the alignment, which can be gapped
    my $pep_offset = ($feature_end - $feature_start + 1) / 3;


    # These are all useful for debugging so I will keep them
    #say "PEP INDEX: ".$pep_index;
    #say "PEP OFFSET: ".$pep_offset;
    #say "EXON PHASE: ".$exon->phase;
    #say "EXON END PHASE: ".$exon->end_phase;
    #say "BASES TO SKIP: ".$bases_to_skip;
    #say "BASES TO REMOVE: ".$bases_to_remove;
    #say "EXON START: ".$exon->seq_region_start;
    #say "EXON END: ".$exon->seq_region_end;
    #say "EXON SEQ: ".$exon->seq->seq;
    #say "FEATURE START: ".$feature_start;
    #say "FEATURE END: ".$feature_end;


    ###########################################################################################################
    # - code that could be used to calculate coverage and pid averages from the alignment on a per exon basis
    ###########################################################################################################
    # char_count is a variable to map the index of a residue in the ungapped translation to the equivalent position in the alignmnet
    #my $char_count = 0;

    # align_start_index and align_end_index are the start and end columns in the alignment that the feature
    # covers. As the feature is ungapped, the alignment columns may span a region of the alignment that is
    # more than just the feature length / 3
    #my $align_start_index;
    #my $align_end_index;

    # This loop finds the columns for the alignment start and end using the position in the ungapped translation sequence
    # (held in pep_index) and the number of complete codons the exon spans (held in pep_offset)
    #for(my $j=0; $j<length($target_seq); $j++) {
    #  my $char = substr($target_seq,$j,1);

    #  if($char_count == $pep_index) {
    #    $align_start_index = $j;
    #  }

     # if($char_count == ($pep_index + $pep_offset) || ($j == length($target_seq) - 1)) {
     #   $align_end_index = $j;
     #   last;
     # }

      # As we're mapping via the ungapped translation, skip characters that are gaps
      #unless($char eq '-') {
      #  $char_count++;
      #}
    #}

    #say "ASI: ".$align_start_index;
    #say "ASE: ".$align_end_index;

    # Pull out the sections of the alignment the feature covers and calculate the percent id and cov
    #my $sub_align_length = $align_end_index - $align_start_index;
    #$query_sub_align = substr($query_seq,$align_start_index,$sub_align_length);
    #$target_sub_align = substr($target_seq,$align_start_index,$sub_align_length);
    #my $exon_coverage = $self->calculate_coverage($query_sub_align,$target_sub_align);
    #my $exon_percent_id = $self->calculate_percent_id($query_sub_align,$target_sub_align);
    #say "EXON COV: ".$exon_coverage;
    #say "EXON PID: ".$exon_percent_id;

    # Build the feature
    my $paf = Bio::EnsEMBL::FeaturePair->new(
                                              -start      => $feature_start,
                                              -end        => $feature_end,
                                              -strand     => $exon->strand,
                                              -hseqname   => $transcript->{'accession'},
                                              -hstart     => $pep_index,
                                              -hend       => ($pep_index + $pep_offset - 1),
                                              -hcoverage  => $cov,
                                              -percent_id => $pid,
                                              -slice      => $exon->slice,
                                              -analysis   => $transcript->analysis);

    push(@{$exon_supporting_features},$paf);
    $pep_index = $pep_index + $pep_offset;
  }

  return($exon_supporting_features);
}

sub calculate_percent_id {
  my ($self,$query_seq,$target_seq) = @_;

  my $match_count = 0;
  for(my $j=0; $j<length($target_seq); $j++) {
    my $char_query = substr($query_seq,$j,1);
    my $char_target = substr($target_seq,$j,1);
    if($char_query eq $char_target) {
      $match_count++;
    }
  }

  my $percent_id = ($match_count / length($target_seq)) * 100;

  return($percent_id);

}

sub calculate_coverage {
  my ($self,$query_seq,$target_seq) = @_;

  my $coverage;
  my $query_gap_count = $query_seq =~ s/\-//g;

  my $ungapped_target_seq = $target_seq;
  $ungapped_target_seq  =~ s/\-//g;

  if(length($ungapped_target_seq) == 0) {
    return 0;
  }

  $coverage = 100 - (($query_gap_count/length($ungapped_target_seq)) * 100);

  return $coverage;
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

sub genome_slices {
  my ($self,$val) = @_;
  if($val) {
    $self->param('_genome_slices',$val);
  }
  return($self->param('_genome_slices'));
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

  my $accession_array = $self->param('iid');

  my $table_adaptor = $self->db->get_NakedTableAdaptor();
  $table_adaptor->table_name('uniprot_sequences');

  my $output_dir = $self->param('query_seq_dir');

  # Note as each accession will occur in only one file, there should be no problem using the first one
  my $outfile_name = "genblast_".${$accession_array}[0].".fasta";
  my $outfile_path = $output_dir."/".$outfile_name;

  my $biotypes_hash = {};

  unless(-e $output_dir) {
    `mkdir $output_dir`;
  }

  if(-e $outfile_path) {
    $self->warning("Found the query file in the query dir already. Overwriting. File path\n:".$outfile_path);
  }

  open(QUERY_OUT,">".$outfile_path);
  foreach my $accession (@{$accession_array}) {
    my $db_row = $table_adaptor->fetch_by_dbID($accession);
    unless($db_row) {
      $self->throw("Did not find an entry int eh uniprot_sequences table matching the accession. Accession:\n".$accession);
    }

    my $seq = $db_row->{'seq'};
    $biotypes_hash->{$accession} = $db_row->{'biotype'};

    my $record = ">".$accession."\n".$seq;
    say QUERY_OUT $record;
  }
  close QUERY_OUT;

  $self->files_to_delete($outfile_path);
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
