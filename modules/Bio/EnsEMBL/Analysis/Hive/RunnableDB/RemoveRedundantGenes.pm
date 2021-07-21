=head1 LICENSE

# Copyright [2019-2020] EMBL-European Bioinformatics Institute
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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Hive::RunnableDB::RemoveRedundantGenes - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::RemoveRedundantGenes;

use warnings ;
use strict;
use feature 'say';

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;
use Bio::EnsEMBL::Utils::Argument qw (rearrange);
use Bio::EnsEMBL::Variation::Utils::FastaSequence qw(setup_fasta);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


sub fetch_input {
  my ($self) = @_;

  if($self->param('skip_analysis')) {
    $self->complete_early('Skip check flag is enabled, so no check will be carried out');
  }

  $self->create_analysis;


  my $target_dba = $self->hrdb_get_dba($self->param_required('target_db'));
  my $dna_dba;
  if($self->param('use_genome_flatfile')) {
    unless($self->param_required('genome_file') && -e $self->param('genome_file')) {
      $self->throw("You selected to use a flatfile to fetch the genome seq, but did not find the flatfile. Path provided:\n".$self->param('genome_file'));
    }
    setup_fasta(
                 -FASTA => $self->param_required('genome_file'),
               );
  } else {
    $dna_dba = $self->hrdb_get_dba($self->param('dna_db'));
    $target_dba->dnadb($dna_dba);
  }

  $self->hrdb_set_con($target_dba,'target_db');

  # Fetch the genes
#  my $gene_adaptor = $target_dba->get_GeneAdaptor;
#  my $genes = $gene_adaptor->fetch_all();
  my $genes = [];
  my $slice_adaptor = $target_dba->get_SliceAdaptor;
  foreach my $slice_name (@{$self->param('iid')}) {
    my $slice = $slice_adaptor->fetch_by_name($slice_name);
    my $initial_genes = $slice->get_all_Genes();
    foreach my $gene (@$initial_genes) {
      # Note that the below means we can use standardised slice lengths and only process genes once. There will be very rare edge cases where
      # a cluster will only be split because of varying terminal exon lengths, but this will not cause any actual problems
      unless($gene->end > $slice->length) {
        push(@$genes,$gene);
      }
    }
  }

  unless(scalar(@$genes)) {
    $self->input_job->autoflow(0);
    $self->complete_early('No genes to process');
  }

  $self->param('input_genes',$genes);
  if($self->param_required('target_type') eq 'biotype_priority') {
    $self->set_biotype_priorities();
  }
}


sub run {
  my ($self) = @_;

  my $target_dba = $self->hrdb_get_con('target_db');
  my $gene_adaptor = $target_dba->get_GeneAdaptor();
  my $genes = $self->param('input_genes');
  my $transcript_strings = {};
  my $count = 0;
  my $total_genes = scalar(@$genes);
  say "Total gene count: ".$total_genes;
  while(my $gene = pop(@$genes)) {
    $count++;
    if($count % 100 == 0) {
      say "Completed: ".$count."/".$total_genes;
    }
    # Only works on a one transcript per gene model
    my $transcript = ${$gene->get_all_Transcripts}[0];
#    my $transcript_string = $transcript->start.":".$transcript->end.":".$transcript->strand.":".$transcript->seq_region_name.":";
    my $transcript_string;
    my $intron_string;
    if ($transcript) {
      $transcript_string = $transcript->seq_region_name.":".$transcript->strand.":";
      $intron_string = $self->generate_intron_string($transcript->get_all_Introns());
      if ($intron_string) {
        $transcript_string .= $intron_string;
      }

      # This will process single exon genes. It won't handle cases where the single exon genes have different start/ends, but this isn't
      # exactly enough of a problem to warrant all the extra code it would take
      unless ($intron_string) {
        $transcript_string .= $transcript->start.":".$transcript->end;
      }
    }

    # For generic data we mostly just consider which has the longest cds if the introns are the same. Then if the cds is the
    # Same we consider which transcript has the most UTR (not in terms of length, but whether both ends have UTR or just one end or none), if
    # this comparison is the same, then the longest transcript is selected (so UTR length at that point)
    if($self->param('target_type') eq 'generic') {
      unless($transcript_strings->{$transcript_string}) {
        $transcript_strings->{$transcript_string} = $gene;
      } else {
        my $existing_transcript = ${$transcript_strings->{$transcript_string}->get_all_Transcripts}[0];
        if($self->compare_transcripts($existing_transcript,$transcript)) {
          $gene_adaptor->remove($transcript_strings->{$transcript_string});
          $transcript_strings->{$transcript_string} = $gene;
        } else {
          $gene_adaptor->remove($gene);
        }
      }
    } # end if($self->param('target_type') eq 'transcriptomic')

    # This is basically the same as the above, but duplicated for clarity. The major difference is that we consider
    # priority as the deciding factor. This means we can't do the same check to see if the start/end are identical
    # as priority is the most important thing to intially consider
    elsif($self->param('target_type') eq 'biotype_priority') {
      # If there's no equivalent biotype in layering then we should just remove the gene
      # There is a key called 'unrecognised_biotype' that has the worst priority that could be used
      # if we ever want to do anything with these
      unless(exists $self->param('biotype_priorities')->{$transcript->biotype}) {
        $gene_adaptor->remove($gene);
        next;
      }

      unless($transcript_strings->{$transcript_string}) {
        $transcript_strings->{$transcript_string} = $gene;
      } else {
        my $existing_transcript = ${$transcript_strings->{$transcript_string}->get_all_Transcripts}[0];
        if($self->compare_biotype_priorities($existing_transcript,$transcript)) {
          $gene_adaptor->remove($transcript_strings->{$transcript_string});
          $transcript_strings->{$transcript_string} = $gene;
        } else {
          $gene_adaptor->remove($gene);
	}
      }
    } # end elsif($self->param('target_type') eq 'biotype_priority')
  }
}


sub write_output {
  my ($self) = @_;
  return;
}


sub compare_transcripts {
  my ($self,$transcript1,$transcript2) = @_;

  my $transcript1_cds_length = $transcript1->cdna_coding_end - $transcript1->cdna_coding_start + 1;
  my $transcript2_cds_length = $transcript2->cdna_coding_end - $transcript2->cdna_coding_start + 1;
  if($transcript1_cds_length > $transcript2_cds_length) {
    return(0);
  } elsif($transcript2_cds_length > $transcript1_cds_length) {
    return(1);
  }

  # At this point the cds is identical so we should look at the UTRs. Give higher priority to having
  # UTR on one or both ends, after that you essentially look at the length
  my $utr1 = 0;
  my $utr2 = 0;

  if($transcript1->five_prime_utr_Feature) {
    $utr1++;
  }

  if($transcript1->three_prime_utr_Feature) {
    $utr1++;
  }

  if($transcript2->five_prime_utr_Feature) {
    $utr2++;
  }

  if($transcript2->three_prime_utr_Feature) {
    $utr2++;
  }

  if($utr1 == $utr2) {
    if($transcript1->length >= $transcript2->length) {
      return(0);
    } else {
      return(1);
    }
  } elsif($utr2 > $utr1) {
    return(1);
  } else {
    return(0);
  }
}


sub compare_biotype_priorities {
  my ($self,$transcript1,$transcript2) = @_;

  my $biotype_priorities = $self->param('biotype_priorities');

  my $transcript1_biotype = $transcript1->biotype;
  my $transcript2_biotype = $transcript2->biotype;

  # Start off by checking if either transcript has a better priority
  if($biotype_priorities->{$transcript1_biotype} < $biotype_priorities->{$transcript2_biotype}) {
    return(0);
  } elsif($biotype_priorities->{$transcript1_biotype} > $biotype_priorities->{$transcript2_biotype}) {
    return(1);
  } else {
    # At this point the priorities are the same. First examine the cds lenght and that the longest
    my $transcript1_cds_length = $transcript1->cdna_coding_end - $transcript1->cdna_coding_start + 1;
    my $transcript2_cds_length = $transcript2->cdna_coding_end - $transcript2->cdna_coding_start + 1;
    if($transcript1_cds_length > $transcript2_cds_length) {
      return(0);
    } elsif($transcript2_cds_length > $transcript1_cds_length) {
      return(1);
    } else {
      # At this point the cds lengths are the same. Get the supporting feature info and pick the transcript
      # with the best combined score, or just transcript1 if the scores are identical
      my $hcoverage1 = ${$transcript1->get_all_supporting_features}[0]->hcoverage;
      my $perc_ident1 = ${$transcript1->get_all_supporting_features}[0]->percent_id;
      my $hcoverage2 = ${$transcript1->get_all_supporting_features}[0]->hcoverage;
      my $perc_ident2 = ${$transcript1->get_all_supporting_features}[0]->percent_id;
      my $combined_score1 = $hcoverage1 + $perc_ident1;
      my $combined_score2 = $hcoverage2 + $perc_ident2;
      if($combined_score1 >= $combined_score2) {
        return(0);
      } else {
        return(1);
      }
    } # end else
  } # end else
}


sub set_biotype_priorities {
  my ($self, $val) = @_;

  my $layers = $self->param_required('layers');
  my $biotype_priorities = {};
  my $max_priority = 0;
  foreach my $layer (@$layers) {
    my $id = $layer->{'ID'};
    unless($id =~ /^LAYER(\d+)/) {
      $self->throw("Issue with parsing layer id. Expected id format of LAYER1, LAYER2... id used: ".$id);
    }
    my $priority = $1;
    if($priority > $max_priority) {
      $max_priority = $priority;
    }
    my $layer_biotypes = $layer->{'BIOTYPES'};
    foreach my $biotype (@$layer_biotypes) {
      if(exists $biotype_priorities->{$biotype}) {
        $self->throw("Found a repeated biotype in the layer array. Biotype: ".$biotype);
      }
      $biotype_priorities->{$biotype} = $priority;
    }
  }

  # It is possible that the db will have stuff that is not actually in layering. I have put in a score that
  # will make it a lower priority than anything else in layering, but in fact we will probably just remove
  # them regardless of whether they're redundant or not
  $biotype_priorities->{'unrecognised_biotype'} = $max_priority + 1;
  $self->param('biotype_priorities',$biotype_priorities);
}


sub generate_intron_string {
  my ($self,$intron_array) = @_;

  my $intron_string = "";
  foreach my $intron (@{$intron_array}) {
    my $start = $intron->start();
    my $end = $intron->end();
    $intron_string .= $start."..".$end.":";
  }

  return($intron_string);
}


1;

