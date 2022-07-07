=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerateForGenewise

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerateForGenewise;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(Transcript_info set_stop_codon set_start_codon attach_Slice_to_Transcript evidence_coverage list_evidence);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ConfigDependent::TranscriptUtils qw( low_complexity_less_than_maximum ) ;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(contains_internal_stops);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(Gene_info print_Gene);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(id coord_string lies_inside_of_slice);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name write_seqfile);

use parent qw(Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBlastMiniGenewise);


sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    peptide_file => [],
  }
}

sub create_bmg_runnables{
  my ($self, $hits_hash, $seqfetcher, $analysis) = @_;
  $self->throw("RunnableDB::BlastMiniGenewise Can't create BlastMiniGenewise ".
        "runnables without a hash of hits") unless($hits_hash);
  $seqfetcher = $self->seqfetcher if(!$seqfetcher);
  $analysis = $self->analysis if(!$analysis);
  my %exonerate_parameters = %{$self->EXONERATE_PARAMETERS};
  my %params_hash = %{$self->parameters_hash};
  if($self->LIMIT_TO_FEATURE_RANGE == 0){
    my $genomic_file = write_seqfile($self->query,
                                     create_file_name("genomic", "fa"));
    $self->genomic_file($genomic_file);
    my @ids = sort keys(%$hits_hash);
    my $seqs = $self->get_protein_sequences(\@ids);
    my $peptide_file = write_seqfile($seqs,
                                     create_file_name("peptide", "fa"));
    $self->peptide_file($peptide_file);
    my $r = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript->new
      (
       -analysis => $self->analysis,
       -target_file => $genomic_file,
       -query_type => 'protein',
       -query_file => $peptide_file,
       -program => $self->analysis->program_file,
       %exonerate_parameters,
       %params_hash,
      );
    $self->runnable($r);
  }else{
    foreach my $id(keys(%$hits_hash)){
      my $features = $hits_hash->{$id};
      foreach my $feature(@$features){
        #This should be made a non overlapping range
        my $start = $feature->seq_region_start - $self->FEATURE_RANGE_PADDING;
        my $end = $feature->seq_region_end + $self->FEATURE_RANGE_PADDING;
        $start = 1 if($start < 1);
        $end = $self->query->seq_region_length
          if($end > $self->query->seq_region_length);
        my $name = ($self->query->coord_system->name.":".
                    $self->query->coord_system->version.":".
                    $self->query->seq_region_name.":".
                    $start.":".$end.":".
                    $self->query->strand);
        my $db = $self->get_dbadaptor($self->GENE_SOURCE_DB);
        my $query = $self->fetch_sequence($name, $db, $self->REPEATMASKING);
#        logger_info("Creating ExonerateTranscript Runnable over a limited range with "
#                    .$id." and ".$seqfetcher." to run on ".$query->name);
        my $genomic_file = write_seqfile($query,
                                            create_file_name("genomic", "fa"));
        $self->genomic_file($genomic_file);
        my $seqs = $self->get_protein_sequences([$id]);
        my $peptide_file = write_seqfile($seqs,
                                            create_file_name("peptide", "fa"));
        $self->peptide_file($peptide_file);
        my $r = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript->new
          (
           -analysis => $self->analysis,
           -target_file => $genomic_file,
           -query_type => 'protein',
           -query_file => $peptide_file,
           -program => $self->analysis->program_file,
           %exonerate_parameters,
           %params_hash,
          );
        $self->runnable($r);
      }
    }
  }
}


sub process_transcripts{
  my ($self, $transcripts) = @_;

  my @complete_transcripts;
  foreach my $transcript(@$transcripts){
    my $query_name = $transcript->start_Exon->seqname;
    my $unmasked_slice = $self->fetch_sequence($query_name,
                                               $self->output_db);
    attach_Slice_to_Transcript($transcript, $unmasked_slice);
    my $start_t = set_start_codon($transcript);
    my $end_t = set_stop_codon($start_t);
    if(contains_internal_stops($transcript)){
      $self->warning(Transcript_info($transcript)." contains internal stops will be thrown ".
              "This will be thrown away");
    }
    push(@complete_transcripts, $end_t);
  }
  foreach my $file($self->genomic_file, @{$self->peptide_file}){
    unlink $file;
  }
  return \@complete_transcripts
}

=head2 filter_genes

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastMiniGenewise
  Arg [2]   : arrayref of Genes
  Function  : filters the gene sets using the defined filter
  Returntype: n/a, stores the genes in the output array
  Exceptions:
  Example   :

=cut


sub filter_genes{
  my ($self, $genes) = @_;

  my $output;
  my $rejected_set;
  my $filtered_set;
  if($self->filter_object){
    ($filtered_set, $rejected_set) =
      $self->filter_object->filter_genes($genes);
    #$self->rejected_set($rejected_set);
    foreach my $gene(@{$filtered_set}){
      my $transcripts = $gene->get_all_Transcripts;
      my $transcript = $transcripts->[0];
      my $evidence = $self->get_Transcript_supporting_evidence($transcript);
      my $coverage = evidence_coverage($transcript, $evidence);
      if(low_complexity_less_than_maximum($transcript, $coverage)){
        push(@$output, $gene);
      }else{
        $gene->description("low_complexity is greater than evidence coverage");
        push(@$rejected_set, $gene);
      }
    }
  }else{
    $self->output($genes);
  }
  $self->output($output);
  $self->rejected_set($rejected_set);
}

sub genomic_file{
  my ($self, $file) = @_;

  if($file){
    $self->param('genomic_file', $file);
  }
  return $self->param('genomic_file');
}

sub peptide_file{
  my ($self, $file) = @_;

  if($file){
    push(@{$self->param('peptide_file')}, $file);
  }
  return $self->param('peptide_file');
}


sub get_protein_sequences{
  my ($self, $ids) = @_;
  my $seqfetcher = $self->seqfetcher;
  my $seqs;
  my $seq;
  foreach my $id(@$ids){
    eval {
      $seq = $seqfetcher->get_Seq_by_acc($id);
    };
    if ($@) {
      $self->throw("Problem fetching sequence for [$id]: [$@]\n");
    }
    if(!$seq){
      $self->throw("Can't fetch the sequence $id with ".$seqfetcher." - it's not in the index\n");
    }
    push(@$seqs, $seq);
  }
  return $seqs;
}


sub get_Transcript_supporting_evidence{
  my ($self, $transcript, $seqfetcher) = @_;

  $seqfetcher = $self->seqfetcher if(!$seqfetcher);
  $self->throw("Can't get ".$transcript." transcrpts supporting features without ".
        "a seqfetcher object") if(!$seqfetcher);
  my $ids = list_evidence($transcript);
  my $sequence;
  foreach my $id(@$ids){
    $sequence = $seqfetcher->get_Seq_by_acc($id);
    last if($sequence);
  }
  return $sequence;
}


1;
