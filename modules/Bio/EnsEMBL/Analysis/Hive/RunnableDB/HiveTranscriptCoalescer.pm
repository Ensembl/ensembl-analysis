=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveTranscriptCoalescer

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveTranscriptCoalescer;

use strict;
use warnings;

use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils qw(cluster_Genes get_single_clusters);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene attach_Analysis_to_Gene);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(has_polyA_signal);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(compute_translation);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils qw(clone_Exon);

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');

=head2 param_defaults

 Description: It allows the definition of default parameters for all inherting module.
              These are the default values:
               max_length_to_merge => 200,
 Returntype : Hashref, containing all default parameters
 Exceptions : None

=cut


sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    source_logic_name => undef,
    biotype_noPolyA => 'bad',
    biotype_PolyA => 'good',
    biotype_retained => 'retained',
    biotype_NMD => 'nonsense_mediated_decay',
    biotype_single => 'single',
    biotype_unclustered => 'low',
    exon_length_multiplier => 7,
    last_junction_length_modifier => 0.66,
    max_length_to_merge => 200,
    max_overlength => 20, # Random value might be around 20
    max_intron_wobble => 15,
    max_exon_wobble => 5,
  }
}


=head2 fetch_input

 Arg [1]    : None
 Description: It will get all the genes in the database for the biotypes requested
              or for all the biotypes
 Returntype : None
 Exceptions : None

=cut

sub fetch_input {
  my ($self) = @_;

  $self->create_analysis;
  my $dna_db = $self->get_database_by_name('dna_db');
  $self->hrdb_set_con($self->hrdb_get_dba($self->param_required('target_db'), $dna_db), 'target_db');
  my $logic_name = $self->param('source_logic_name');
  my @genes;
  foreach my $input_db (@{$self->input_dbs}) {
    my $db = $self->hrdb_get_dba($input_db, $dna_db);
    my $slice = $self->fetch_sequence($self->input_id, $db);
    if ($self->get_biotypes and scalar(@{$self->get_biotypes})) {
      foreach my $biotype (@{$self->get_biotypes}) {
        foreach my $gene (@{$slice->get_all_Genes_by_type($biotype, $logic_name, 1)}) {
          foreach my $transcript (@{$gene->get_all_Transcripts}) {
            $transcript->load;
          }
          push(@genes, $gene);
        }
      }
    }
    else {
      foreach my $gene ($slice->get_all_Genes($logic_name, undef, 1)) {
        foreach my $transcript (@{$gene->get_all_Transcripts}) {
          $transcript->load;
        }
        push(@genes, $gene);
      }
    }
  }
  if (scalar(@genes)) {
    print scalar(@genes), " genes to process\n";
    $self->param('genes', \@genes);
  }
  else {
    $self->input_job->autoflow(0);
    $self->complete_early('No genes to process');
  }
}

sub process_genes {
  my ($self, $good_gene, $gene_to_process) = @_;

  my $max_overlength = $self->param('max_overlength');
  my $exon_length_multiplier = $self->param('exon_length_multiplier');
  my $good_transcript = $good_gene->get_all_Transcripts->[0]; # We know that we only have 1 transcript per gene
  my $good_name = $good_transcript->get_all_supporting_features->[0]->hseqname;
  my $good_hashkey = build_hashkey($good_transcript, 'intron');
  my $transcript_to_process = $gene_to_process->get_all_Transcripts->[0]; # We know that we only have 1 transcript per gene
  my $hash_key_to_process = build_hashkey($transcript_to_process, 'intron');
  my $name_tp = $transcript_to_process->get_all_supporting_features->[0]->hseqname;
  if ($hash_key_to_process) { # If it is a single exon gene the hashkey will be empty
    if ($good_hashkey eq $hash_key_to_process) {
      print STDERR ' T P ', $name_tp, ' ', $transcript_to_process->seq_region_start, ' ', $transcript_to_process->seq_region_end, ' ', abs($good_transcript->end_Exon->length-$transcript_to_process->end_Exon->length), ' ', $good_gene->{full_length}, ' ', $hash_key_to_process, "\n";
      return 2;
    }
    elsif ($good_hashkey =~ /$hash_key_to_process/) { # If the transcript is a fragment we want to add it to our model unless it might be a retained intron
      my $first_exon_to_process = $transcript_to_process->start_Exon;
      my $last_exon_to_process = $transcript_to_process->end_Exon;
      my $strand = $good_transcript->strand;
      my $first_good_exon = $good_transcript->start_Exon;
      foreach my $good_exon (@{$good_transcript->get_all_Exons}) {
        if (($strand == 1
            and $first_exon_to_process->end == $good_exon->end
            and $good_exon->start - $first_exon_to_process->start > $max_overlength
            and $good_exon != $good_transcript->start_Exon)
          or ($strand == -1
            and $first_exon_to_process->start == $good_exon->start
            and $good_exon->end - $first_exon_to_process->end > $max_overlength
            and $good_exon != $good_transcript->start_Exon)) {
          if ($first_good_exon->length < 21 and $first_good_exon == $good_transcript->start_Exon) {
            print STDERR ' T FRS ', $name_tp, ' ', $transcript_to_process->seq_region_start, ' ', $transcript_to_process->seq_region_end, ' ', $good_exon->seq_region_start, ' ', $good_exon->seq_region_end, ' ', $first_exon_to_process->seq_region_start, ' ', $first_exon_to_process->seq_region_end, ' ', $max_overlength, ' ', $hash_key_to_process, "\n";
          }
          elsif ($first_exon_to_process->length < $good_exon->length*$exon_length_multiplier) {
            print STDERR ' T LRS ', $name_tp, ' ', $transcript_to_process->seq_region_start, ' ', $transcript_to_process->seq_region_end, ' ', $good_exon->seq_region_start, ' ', $good_exon->seq_region_end, ' ', $first_exon_to_process->seq_region_start, ' ', $first_exon_to_process->seq_region_end, ' ', $max_overlength, ' ', $hash_key_to_process, "\n";
          }
          elsif ($good_exon->length < 21) {
            print STDERR ' T SRS ', $name_tp, ' ', $transcript_to_process->seq_region_start, ' ', $transcript_to_process->seq_region_end, ' ', $good_exon->seq_region_start, ' ', $good_exon->seq_region_end, ' ', $first_exon_to_process->seq_region_start, ' ', $first_exon_to_process->seq_region_end, ' ', $max_overlength, ' ', $hash_key_to_process, "\n";
          }
          else {
            $gene_to_process->{retained} = 1;
            print STDERR ' T RS ', $name_tp, ' ', $transcript_to_process->seq_region_start, ' ', $transcript_to_process->seq_region_end, ' ', $good_exon->seq_region_start, ' ', $good_exon->seq_region_end, ' ', $first_exon_to_process->seq_region_start, ' ', $first_exon_to_process->seq_region_end, ' ', $max_overlength, ' ', $hash_key_to_process, "\n";
            return 0;
          }
        }
      }
      my $last_good_exon = $good_transcript->end_Exon;
      if (($strand == 1 and $last_exon_to_process->end < $last_good_exon->start)
          or ($strand == -1 and $last_exon_to_process->start > $last_good_exon->end)) {
        $gene_to_process->{partial} = 1;
        print STDERR ' T PARTIAL ', $name_tp, ' ', $transcript_to_process->seq_region_start, ' ', $transcript_to_process->seq_region_end, ' ', $last_good_exon->seq_region_start, ' ', $last_good_exon->seq_region_end, ' ', $last_exon_to_process->seq_region_start, ' ', $last_exon_to_process->seq_region_end, ' ', $max_overlength, ' ', $hash_key_to_process, "\n";
      }
      return 1;
    }
  }
  else {
# Here we process the single exon models
    if (@{$good_transcript->get_all_Exons} == 1) {
      print STDERR ' T SOS ', $name_tp, ' ', $transcript_to_process->seq_region_start, ' ', $transcript_to_process->seq_region_end, ' ', $good_transcript->seq_region_start, ' ', $good_transcript->seq_region_end, "\n";
      return 1;
    }
    else {
      my $good_last_exon = $good_transcript->end_Exon;
      if ($good_last_exon->start < $transcript_to_process->end
          and $good_last_exon->end > $transcript_to_process->start) {
        if ($good_last_exon->strand == 1) {
          if ($good_last_exon->start-$transcript_to_process->start < $max_overlength) {
  print STDERR ' T S ', $name_tp, ' ', $transcript_to_process->seq_region_start, ' ', $transcript_to_process->seq_region_end, ' ', ($good_last_exon->start-$transcript_to_process->start), ' < ', $max_overlength, ' ', $good_transcript->seq_region_start, ' ', $good_transcript->seq_region_end, "\n";
            return 1;
          }
          else {
            $gene_to_process->{retained} = 1; # Mark the gene as retained, might be a bit harsh
  print STDERR ' T SR ', $name_tp, ' ', $transcript_to_process->seq_region_start, ' ', $transcript_to_process->seq_region_end, ' ', ($good_last_exon->start-$transcript_to_process->start), ' > ', $max_overlength, ' ', $good_transcript->seq_region_start, ' ', $good_transcript->seq_region_end, "\n";
            return 0;
          }
        }
        else {
          if ($transcript_to_process->end-$good_last_exon->end < $max_overlength) {
  print STDERR ' T S ', $name_tp, ' ', $transcript_to_process->seq_region_start, ' ', $transcript_to_process->seq_region_end, ' ', ($transcript_to_process->end-$good_last_exon->end), ' < ', $max_overlength, ' ', $good_transcript->seq_region_start, ' ', $good_transcript->seq_region_end, "\n";
            return 1;
          }
          else {
            $gene_to_process->{retained} = 1; # Mark the gene as retained, might be a bit harsh
  print STDERR ' T SR ', $name_tp, ' ', $transcript_to_process->seq_region_start, ' ', $transcript_to_process->seq_region_end, ' ', ($transcript_to_process->end-$good_last_exon->end), ' < ', $max_overlength, ' ', $good_transcript->seq_region_start, ' ', $good_transcript->seq_region_end, "\n";
            return 0;
          }
        }
      }
      else {
  print STDERR ' T SB ', $name_tp, ' ', $transcript_to_process->seq_region_start, ' ', $transcript_to_process->seq_region_end, ' ', $good_transcript->seq_region_start, ' ', $good_transcript->seq_region_end, "\n";
        return 0;
      }
    }
  }
}

=head2 run

 Arg [1]    : None
 Description: Merge PacBio reads depending on several factors:
              If transcript 2 begin less than 'max_length_to_merge' downstream of transcript 1,
              the transcripts are merged. The value is cumulative, if transcript 3 would have been
              merged in transcript 2, it is merge in transcript 1.
 Returntype : 
 Exceptions : 

=cut

sub run {
  my ($self) = @_;

  if ($self->param('disconnect_jobs')) {
    $self->dbc->disconnect_if_idle;
    $self->hrdb_get_con('target_db')->dnadb->dbc->disconnect_if_idle;
  }
  print STDERR 'Clustering genes';
  my ($clusters, $unclustered) = cluster_Genes($self->param('genes'), $self->get_hashtypes);
  print STDERR " Done\n";
  print STDERR 'Working on unclustered genes', "\n";
  my $single_biotype = $self->param('biotype_unclustered');
  foreach my $cluster (@$unclustered) {
    foreach my $gene (@{$cluster->get_Genes}) {
      foreach my $transcript (@{$gene->get_all_Transcripts}) {
        $gene->{polyA_signal} = rank_polyA_signal($transcript);
        $gene->{full_length} = 1;
        compute_translation($transcript);
        $transcript->biotype($single_biotype);
      }
      $gene->biotype($single_biotype);
      $self->output([$gene]);
    }
  }
  print STDERR " Done\n";
  my $i = 0;
  foreach my $cluster (@$clusters) {
    print STDERR ('-' x10), "\n", 'DEBUG CLUSTER ', join(' ', ++$i, $cluster->start, $cluster->end, $cluster->strand), "\n";
    my $genes = $cluster->get_Genes;
    my $num_genes = scalar(@$genes);
    print STDERR "Working on $num_genes\n";
    print STDERR ('-'x5), 'STEP COLLAPSING BY INTRON STRUCTURE', "\n";
# First we are looking at all the models to collapse the one which have the same intron structure
    my @to_process = sort {$b->length <=> $a->length} @$genes;
    my @same_intron_structure;
    my $max_overlength = $self->param('max_overlength');
    my $exon_length_multiplier = $self->param('exon_length_multiplier');
    for (my $index = 0; $index < @to_process; $index++) {
      my $good_gene = $to_process[$index];
      $good_gene->{full_length} = 1;
      next if (exists $good_gene->{processed});
      my $good_transcript = $good_gene->get_all_Transcripts->[0]; # We know that we only have 1 transcript per gene
      my $good_name = $good_transcript->get_all_supporting_features->[0]->hseqname;
      $good_gene->{polyA_signal} = rank_polyA_signal($good_transcript);
      if (exists $good_gene->{partial} and !$good_gene->{polyA_signal}) {
        compute_translation($good_transcript);
        $good_gene->biotype('partial');
        $self->output([$good_gene]);
        print STDERR ' PARTIAL GENE ', $good_name, ' ', $good_transcript->seq_region_start, ' ', $good_transcript->seq_region_end, ' ', $good_transcript->strand, ' ', $good_gene->{polyA_signal}, "\n";
      }
      else {
        my $good_hashkey = build_hashkey($good_transcript, 'intron');
        print STDERR 'GOOD T ', $good_name, ' ', $good_transcript->seq_region_start, ' ', $good_transcript->seq_region_end, ' ', $good_transcript->strand, ' ', $good_gene->{polyA_signal}, "\t", $good_hashkey, "\n";
        for (my $jndex = $index+1; $jndex < @to_process; $jndex++) {
          my $gene_to_process = $to_process[$jndex];
          my $transcript_to_process = $gene_to_process->get_all_Transcripts->[0]; # We know that we only have 1 transcript per gene
          my $name_tp = $transcript_to_process->get_all_supporting_features->[0]->hseqname;
          my $hash_key_to_process = build_hashkey($transcript_to_process, 'intron');
          my $result = $self->process_genes($good_gene, $gene_to_process);
          if ($result == 1) {
            $gene_to_process->{processed} = 1; # Mark the gene as processed
            push(@{$good_gene->{genes}}, $gene_to_process);
            print STDERR ' T G ', $name_tp, ' ', $transcript_to_process->seq_region_start, ' ', $transcript_to_process->seq_region_end, "\t", $hash_key_to_process, "\n";
          }
          elsif ($result == 2) {
            $gene_to_process->{processed} = 1; # Mark the gene as processed
            ++$good_gene->{full_length}; # We want to keep the number of full length support
            push(@{$good_gene->{full_genes}}, $gene_to_process);
            print STDERR ' T P ', $name_tp, ' ', $transcript_to_process->seq_region_start, ' ', $transcript_to_process->seq_region_end, ' ', abs($good_transcript->end_Exon->length-$transcript_to_process->end_Exon->length), ' ', $good_gene->{full_length}, ' ', $hash_key_to_process, "\n";
          }
        }
        push(@same_intron_structure, $good_gene);
      }
    }
    print STDERR "\n", ('-'x5), 'STEP REMOVING SPLICE SITE DIFFERENCE', "\n";
# Now we are looking at the models to see if the difference in intron structure was caused by the error rate of the PacBio.
    # $$$$$$$$$$---------$$$$$$$
    # ######---------###########
    # ######---------###########
    # ######---------###########
    # ######---------###########
    my @splice_site_difference;
    my $gene_index = 1;
    my %good_splice_sites;
    foreach my $gene (@same_intron_structure) {
      foreach my $intron (@{$gene->get_all_Transcripts->[0]->get_all_Introns}) {
        $good_splice_sites{$intron->start.':'.$intron->end} += $gene->{full_length};
      }
    }
    my $max_intron_wobble = $self->param('max_intron_wobble');
    my $max_exon_wobble = $self->param('max_exon_wobble');
    foreach my $good_gene (sort {$a->start <=> $b->start || $b->end <=> $a->end} @same_intron_structure) {
      my $good_transcript = $good_gene->get_all_Transcripts->[0];
      my $good_introns = $good_transcript->get_all_Introns;
      my $good_name = $good_transcript->get_all_supporting_features->[0]->hseqname;
      for (my $index_to_process = $gene_index; $index_to_process < @same_intron_structure; $index_to_process++) {
        my $transcript_to_process = $same_intron_structure[$index_to_process]->get_all_Transcripts->[0];
        my $introns_to_process = $transcript_to_process->get_all_Introns;
        my $name_to_process = $transcript_to_process->get_all_supporting_features->[0]->hseqname;
        print STDERR "COMPARING $good_name WITH $name_to_process\n";
        foreach my $good_intron (@$good_introns) {
          my $good_intron_key = $good_intron->start.':'.$good_intron->end;
          foreach my $intron_to_process (@$introns_to_process) {
            if ($good_intron->start == $intron_to_process->start
                and $good_intron->end == $intron_to_process->end) {
              last;
            }
            elsif ($good_intron->end > $intron_to_process->start and $good_intron->start < $intron_to_process->end) {
              print STDERR "  INTRON OVERLAP\n";
              my $key_to_process = $intron_to_process->start.':'.$intron_to_process->end;
              my $good_prev = $good_intron->prev_Exon;
              my $good_next = $good_intron->next_Exon;
              my $prev_to_process = $intron_to_process->prev_Exon;
              my $next_to_process = $intron_to_process->next_Exon;
              if ($good_intron->length == $intron_to_process->length) {
                print STDERR "DEBUG SAME LENGTH\n";
                if ($good_transcript->strand == 1 and abs($good_prev->end-$prev_to_process->end) < $max_intron_wobble) {
                  print STDERR "DEBUG MATCH INTRON WOBBLE $max_intron_wobble\n";
                  if ($good_splice_sites{$good_intron_key} > $good_splice_sites{$key_to_process}) {
                    print STDERR $good_name, ' P ', $name_to_process, ' ', $prev_to_process->seq_region_end, ' ', $prev_to_process->end, ' ', $good_prev->end, ' ', $good_splice_sites{$good_intron_key}, ' ', $good_splice_sites{$key_to_process}, "\n";
                    $prev_to_process->{_gb_end} = $good_prev->end;
                    print STDERR $good_name, ' N ', $name_to_process, ' ', $next_to_process->seq_region_end, ' ', $next_to_process->end, ' ', $good_next->end, ' ', $good_splice_sites{$good_intron_key}, ' ', $good_splice_sites{$key_to_process}, "\n";
                    $next_to_process->{_gb_start} = $good_next->start;
                    ++$same_intron_structure[$index_to_process]->{to_process};
                  }
                  elsif ($good_splice_sites{$good_intron_key} < $good_splice_sites{$key_to_process}) {
                    print STDERR $good_name, ' GP ', $name_to_process, ' ', $good_prev->seq_region_end, ' ', $good_prev->end, ' ', $prev_to_process->end, ' ', $good_splice_sites{$good_intron_key}, ' ', $good_splice_sites{$key_to_process}, "\n";
                    $good_prev->{_gb_end} = $prev_to_process->end;
                    print STDERR $good_name, ' GN ', $name_to_process, ' ', $good_next->seq_region_start, ' ', $good_next->start, ' ', $next_to_process->start, ' ', $good_splice_sites{$good_intron_key}, ' ', $good_splice_sites{$key_to_process}, "\n";
                    $good_next->{_gb_start} = $next_to_process->start;
                    ++$good_gene->{to_process};
                  }
                }
                elsif ($good_transcript->strand == -1 and abs($good_prev->start-$prev_to_process->start) < $max_intron_wobble) {
                  if ($good_splice_sites{$good_intron_key} > $good_splice_sites{$key_to_process}) {
                    print STDERR $good_name, ' P ', $name_to_process, ' ', $prev_to_process->seq_region_start, ' ', $prev_to_process->start, ' ', $good_prev->start, ' ', $good_splice_sites{$good_intron_key}, ' ', $good_splice_sites{$key_to_process}, "\n";
                    $prev_to_process->{_gb_start} = $good_prev->start;
                    print STDERR $good_name, ' N ', $name_to_process, ' ', $next_to_process->seq_region_end, ' ', $next_to_process->end, ' ', $good_next->end, ' ', $good_splice_sites{$good_intron_key}, ' ', $good_splice_sites{$key_to_process}, "\n";
                    $next_to_process->{_gb_end} = $good_next->end;
                    ++$same_intron_structure[$index_to_process]->{to_process};
                  }
                  elsif ($good_splice_sites{$good_intron_key} < $good_splice_sites{$key_to_process}) {
                    print STDERR $good_name, ' GP ', $name_to_process, ' ', $good_prev->seq_region_start, ' ', $good_prev->start, ' ', $prev_to_process->start, ' ', $good_splice_sites{$good_intron_key}, ' ', $good_splice_sites{$key_to_process}, "\n";
                    $good_prev->{_gb_start} = $prev_to_process->start;
                    print STDERR $good_name, ' GN ', $name_to_process, ' ', $good_next->seq_region_end, ' ', $good_next->end, ' ', $next_to_process->end, ' ', $good_splice_sites{$good_intron_key}, ' ', $good_splice_sites{$key_to_process}, "\n";
                    $good_next->{_gb_end} = $next_to_process->end;
                    ++$good_gene->{to_process};
                  }
                }
                last;
              }
              elsif ($good_transcript->strand == 1) {
                if ($good_prev->end != $prev_to_process->end and abs($good_prev->end-$prev_to_process->end) < $max_exon_wobble) {
                  print STDERR "DEBUG MATCH PREV EXON WOBBLE $max_exon_wobble\n";
                  if ($good_splice_sites{$good_intron_key} > $good_splice_sites{$key_to_process}) {
                    print STDERR $good_name, ' P ', $name_to_process, ' ', $prev_to_process->seq_region_end, ' ', $prev_to_process->end, ' ', $good_prev->end, ' ', $good_splice_sites{$good_intron_key}, ' ', $good_splice_sites{$key_to_process}, "\n";
                    $prev_to_process->{_gb_end} = $good_prev->end;
                    ++$same_intron_structure[$index_to_process]->{to_process};
                  }
                  elsif ($good_splice_sites{$good_intron_key} < $good_splice_sites{$key_to_process}) {
                    print STDERR $good_name, ' GP ', $name_to_process, ' ', $good_prev->seq_region_end, ' ', $good_prev->end, ' ', $prev_to_process->end, ' ', $good_splice_sites{$good_intron_key}, ' ', $good_splice_sites{$key_to_process}, "\n";
                    $good_prev->{_gb_end} = $prev_to_process->end;
                    ++$good_gene->{to_process};
                  }
                }
                elsif ($good_next->start != $next_to_process->start and abs($good_next->start-$next_to_process->start) < $max_exon_wobble) {
                  print STDERR "DEBUG MATCH NEXT EXON WOBBLE $max_exon_wobble\n";
                  if ($good_splice_sites{$good_intron_key} > $good_splice_sites{$key_to_process}) {
                    print STDERR $good_name, ' N ', $name_to_process, ' ', $next_to_process->seq_region_start, ' ', $next_to_process->start, ' ', $good_next->start, ' ', $good_splice_sites{$good_intron_key}, ' ', $good_splice_sites{$key_to_process}, "\n";
                    $next_to_process->{_gb_start} = $good_next->start;
                    ++$same_intron_structure[$index_to_process]->{to_process};
                  }
                  elsif ($good_splice_sites{$good_intron_key} < $good_splice_sites{$key_to_process}) {
                    print STDERR $good_name, ' GN ', $name_to_process, ' ', $good_next->seq_region_start, ' ', $good_next->start, ' ', $next_to_process->start, ' ', $good_splice_sites{$good_intron_key}, ' ', $good_splice_sites{$key_to_process}, "\n";
                    $good_next->{_gb_start} = $next_to_process->start;
                    ++$good_gene->{to_process};
                  }
                }
              }
              else {
                if ($good_prev->start != $prev_to_process->start and abs($good_prev->start-$prev_to_process->start) < $max_exon_wobble) {
                  print STDERR "DEBUG MATCH PREV EXON WOBBLE $max_exon_wobble\n";
                  if ($good_splice_sites{$good_intron_key} > $good_splice_sites{$key_to_process}) {
                    print STDERR $good_name, ' P ', $name_to_process, ' ', $prev_to_process->seq_region_start, ' ', $prev_to_process->start, ' ', $good_prev->start, ' ', $good_splice_sites{$good_intron_key}, ' ', $good_splice_sites{$key_to_process}, "\n";
                    $prev_to_process->{_gb_start} = $good_prev->start;
                    ++$same_intron_structure[$index_to_process]->{to_process};
                  }
                  elsif ($good_splice_sites{$good_intron_key} < $good_splice_sites{$key_to_process}) {
                    print STDERR $good_name, ' GP ', $name_to_process, ' ', $good_prev->seq_region_start, ' ', $good_prev->start, ' ', $prev_to_process->start, ' ', $good_splice_sites{$good_intron_key}, ' ', $good_splice_sites{$key_to_process}, "\n";
                    $good_prev->{_gb_start} = $prev_to_process->start;
                    ++$good_gene->{to_process};
                  }
                }
                elsif ($good_prev->end != $prev_to_process->end and abs($good_prev->end-$prev_to_process->end) < $max_exon_wobble) {
                  print STDERR "DEBUG MATCH NEXT EXON WOBBLE $max_exon_wobble\n";
                  if ($good_splice_sites{$good_intron_key} > $good_splice_sites{$key_to_process}) {
                    print STDERR $good_name, ' N ', $name_to_process, ' ', $prev_to_process->seq_region_end, ' ', $prev_to_process->end, ' ', $good_prev->end, ' ', $good_splice_sites{$good_intron_key}, ' ', $good_splice_sites{$key_to_process}, "\n";
                    $prev_to_process->{_gb_end} = $good_prev->end;
                    ++$same_intron_structure[$index_to_process]->{to_process};
                  }
                  elsif ($good_splice_sites{$good_intron_key} < $good_splice_sites{$key_to_process}) {
                    print STDERR $good_name, ' GN ', $name_to_process, ' ', $good_prev->seq_region_end, ' ', $good_prev->end, ' ', $prev_to_process->end, "\n";
                    $good_prev->{_gb_end} = $prev_to_process->end;
                    ++$good_gene->{to_process};
                  }
                }
                else {
                  print STDERR '  ', $good_name, ' AGAINST ', $name_to_process, ' ', $prev_to_process->seq_region_start, ' ', $prev_to_process->seq_region_start, ' ', $good_prev->seq_region_start, ' ', $good_prev->seq_region_end, ' ', $good_prev->seq_region_end, ' ', $prev_to_process->seq_region_end, ' ', $good_splice_sites{$good_intron_key}, ' ', $good_splice_sites{$key_to_process}, "\n";
                }
              }
              last;
            }
          }
        }
      }
    }
    foreach my $gene (@same_intron_structure) {
      if (exists $gene->{to_process}) {
        foreach my $transcript (@{$gene->get_all_Transcripts}) {
          print STDERR ' NEW GENE ', $transcript->get_all_supporting_features->[0], "\n";
          my @exons;
          foreach my $exon (@{$transcript->get_all_Exons}) {
            my $new_exon = clone_Exon($exon);
            $new_exon->start($exon->{_gb_start}) if (exists $exon->{_gb_start});
            $new_exon->end($exon->{_gb_end}) if (exists $exon->{_gb_end});
            push(@exons, $new_exon);
          }
          my $new_transcript = Bio::EnsEMBL::Transcript->new(-exons => \@exons);
          $new_transcript->add_supporting_features($transcript->get_all_supporting_features->[0]);
          my $new_gene = Bio::EnsEMBL::Gene->new();
          $new_gene->{polyA_signal} = $gene->{polyA_signal};
          $new_gene->{full_length} = $gene->{full_length};
          $new_gene->{processed} = $gene->{processed} if (exists $gene->{processed});
          $new_gene->add_Transcript($new_transcript);
          $new_gene->biotype($new_transcript->biotype);
          push(@splice_site_difference, $new_gene);
        }
      }
      else {
        push(@splice_site_difference, $gene);
      }
    }
    print STDERR "\n", ('-'x5), 'STEP COLLAPSE MODELS AGAIN', "\n";
    my @to_reprocess = sort {$b->length <=> $a->length} @splice_site_difference;
    my @after_collapse_redone;
    for (my $index = 0; $index < @to_reprocess; $index++) {
      my $good_gene = $to_reprocess[$index];
      next if (exists $good_gene->{processed} and $good_gene->{processed} == 2);
      my $good_transcript = $good_gene->get_all_Transcripts->[0]; # We know that we only have 1 transcript per gene
      my $good_name = $good_transcript->get_all_supporting_features->[0]->hseqname;
      if (exists $good_gene->{partial}) {
        print STDERR ' PARTIAL GENE ', $good_name, ' ', $good_transcript->seq_region_start, ' ', $good_transcript->seq_region_end, ' ', $good_transcript->strand, ' ', $good_gene->{polyA_signal}, "\n";
      }
      my $good_hashkey = build_hashkey($good_transcript, 'intron');
      print STDERR 'GOOD T ', $good_name, ' ', $good_transcript->seq_region_start, ' ', $good_transcript->seq_region_end, ' ', $good_transcript->strand, ' ', $good_gene->{polyA_signal}, "\t", $good_hashkey, "\n";
      for (my $jndex = $index+1; $jndex < @to_reprocess; $jndex++) {
        my $gene_to_process = $to_reprocess[$jndex];
        my $transcript_to_process = $gene_to_process->get_all_Transcripts->[0]; # We know that we only have 1 transcript per gene
        my $name_tp = $transcript_to_process->get_all_supporting_features->[0]->hseqname;
        my $hash_key_to_process = build_hashkey($transcript_to_process, 'intron');
        my $result = $self->process_genes($good_gene, $gene_to_process);
        if ($result == 1) {
          $gene_to_process->{processed} = 2; # Mark the gene as processed
          push(@{$good_gene->{genes}}, $gene_to_process);
          push(@{$good_gene->{genes}}, @{$gene_to_process->{full_genes}}) if (exists $gene_to_process->{full_genes});
          delete $gene_to_process->{full_genes};
          delete $gene_to_process->{genes};
          print STDERR ' T G ', $name_tp, ' ', $transcript_to_process->seq_region_start, ' ', $transcript_to_process->seq_region_end, "\t", $hash_key_to_process, "\n";
        }
        elsif ($result == 2) {
          $gene_to_process->{processed} = 2; # Mark the gene as processed
          $good_gene->{full_length} += $gene_to_process->{full_length}; # We want to keep the number of full length support
          push(@{$good_gene->{full_genes}}, $gene_to_process);
          push(@{$good_gene->{full_genes}}, @{$gene_to_process->{full_genes}}) if (exists $gene_to_process->{full_genes});
          delete $gene_to_process->{full_genes};
          delete $gene_to_process->{genes};
          print STDERR ' T P ', $name_tp, ' ', $transcript_to_process->seq_region_start, ' ', $transcript_to_process->seq_region_end, ' ', abs($good_transcript->end_Exon->length-$transcript_to_process->end_Exon->length), ' ', $good_gene->{full_length}, ' ', $hash_key_to_process, "\n";
        }
      }
      push(@after_collapse_redone, $good_gene);
    }
    print STDERR "\n", ('-'x5), 'STEP REMOVE RETAINED INTRON', "\n";
# Now I'm looking at retained intron
# So I have to look at all possible intron then check if an exon is overlapping.
# If I have a exon overlapping this intron and it has the same boundaries as the surrounding
# exons, I have a retained exon
    my @for_step4 = sort {$b->length <=> $a->length} @after_collapse_redone;
    my $last_junction_length_modifier = $self->param('last_junction_length_modifier');
    foreach my $good_gene (@for_step4) {
      my $good_transcript = $good_gene->get_all_Transcripts->[0];
      my $good_introns = $good_transcript->get_all_Introns;
      foreach my $intron (@$good_introns) {
        ++$good_gene->{nc_count} unless ($intron->is_splice_canonical);
      }
      my $good_name = $good_transcript->get_all_supporting_features->[0]->hseqname;
      GENE: foreach my $gene_to_process (@for_step4) {
        next if ($good_gene == $gene_to_process or exists $gene_to_process->{retained});
        my $transcript_to_process = $gene_to_process->get_all_Transcripts->[0];
        my $name_tp = $transcript_to_process->get_all_supporting_features->[0]->hseqname;
        my $exons_to_process = $transcript_to_process->get_all_Exons;
        my $first_exon = $transcript_to_process->start_Exon;
        my $last_exon = $transcript_to_process->end_Exon;
#        $last_exon = $transcript_to_process->start_Exon if ($transcript_to_process->strand == -1);
      print STDERR 'GOOD T ', $good_name, ' ', $good_transcript->seq_region_start, ' ', $good_transcript->seq_region_end, ' ', $good_transcript->strand, ' ', $first_exon->seq_region_start, ' ', $first_exon->seq_region_end, ' ', $last_exon->seq_region_start, ' ', $last_exon->seq_region_end, "\n";
        my $exon_index = 0;
        foreach my $good_intron (@$good_introns) {
          foreach my $exon (@$exons_to_process) {
            if ($exon == $last_exon) {
              print STDERR ' L ', $name_tp, ' ', $good_intron->seq_region_start, ' ', $good_intron->seq_region_end, ' ', $exon->seq_region_start, ' ', $exon->seq_region_end, ' ', $last_exon->seq_region_start, ' ', $last_exon->seq_region_end, "\n";
#              print STDERR "Not looking, last exon\n";
            }
            elsif ($exon->start < $good_intron->end and $exon->end > $good_intron->start
                and $exon->start < $good_intron->start and $exon->end > $good_intron->end) {
#              print STDERR $exon->length, ' ', ($good_intron->length/3), ' ', $good_intron->is_splice_canonical, "\n" if ($exon == $first_exon);
              if ($exon == $first_exon) {
                if (!$good_intron->is_splice_canonical) {
#                  print STDERR 'Go next as first exon and NC ', $good_name, ' for ', $name_tp, "\n";
              print STDERR ' F NC ', $name_tp, ' ', $good_intron->seq_region_start, ' ', $good_intron->seq_region_end, ' ', $exon->seq_region_start, ' ', $exon->seq_region_end, ' ', $first_exon->seq_region_start, ' ', $first_exon->seq_region_end, "\n";
                  next;
                }
                if ($exon->length < $good_intron->length/3) {
#                  print STDERR 'Go next as first exon NB ', $good_name, ' for ', $name_tp, "\n";
              print STDERR ' SIL ', $name_tp, ' ', $good_intron->seq_region_start, ' ', $good_intron->seq_region_end, ' ', $exon->seq_region_start, ' ', $exon->seq_region_end, ' ', $first_exon->seq_region_start, ' ', $first_exon->seq_region_end, ' ', $exon->length, ' ', ($good_intron->length/3),"\n";
                  next;
                }
              }
              elsif ($good_intron->next_Exon == $good_transcript->end_Exon and $exon->length < ($good_intron->prev_Exon->length+$good_intron->length+$good_intron->next_Exon->length)*$last_junction_length_modifier) {
                print STDERR ' SAL ', $name_tp, ' ', $good_intron->seq_region_start, ' ', $good_intron->seq_region_end, ' ', $exon->seq_region_start, ' ', $exon->seq_region_end, ' ', $first_exon->seq_region_start, ' ', $first_exon->seq_region_end, ' ', $exon->length, ' ', (($good_intron->prev_Exon->length+$good_intron->length+$good_intron->next_Exon->length)*$last_junction_length_modifier),"\n";
                next;
              }
              else {
                print STDERR ' FSAL ', $name_tp, ' ', $good_intron->seq_region_start, ' ', $good_intron->seq_region_end, ' ', $exon->seq_region_start, ' ', $exon->seq_region_end, ' ', $first_exon->seq_region_start, ' ', $first_exon->seq_region_end, ' ', $exon->length, ' ', (($good_intron->prev_Exon->length+$good_intron->length+$good_intron->next_Exon->length)*$last_junction_length_modifier), ' ', $good_intron->next_Exon, ' ', $good_transcript->end_Exon, ' ', $good_intron->next_Exon->start, ' ', $good_intron->next_Exon->end, ' ', $good_transcript->end_Exon->start, ' ', $good_transcript->end_Exon->end, "\n";
              }
              $gene_to_process->{retained} = 1;
              print STDERR ' R ', $name_tp, ' ', $good_intron->seq_region_start, ' ', $good_intron->seq_region_end, ' ', $exon->seq_region_start, ' ', $exon->seq_region_end, "\n";
#              print STDERR 'Retained against ', $good_gene->dbID, ' ', $good_transcript->dbID, ' for ', $gene_to_process->dbID, ' ', $transcript_to_process->dbID, "\n";
#              print STDERR 'Retained against ', $good_name, ' ', $good_gene->{full_length}, ' for ', $name_tp, ' ', $gene_to_process->{full_length}, "\n";
#            print STDERR sprintf("%d < %d and %d > %d and %d > %d and %d < %d\n", $exon->start, $good_intron->end, $exon->end, $good_intron->start, $exon->start, $good_intron->start, $exon->end, $good_intron->end);
              next GENE; # Things you have to do for money...
            }
          }
        }
      }
    }
    my @retained = grep {exists $_->{retained}} @for_step4;
    my $biotype_retained = $self->param('biotype_retained');
    foreach my $gene (@retained) {
      $gene->biotype($biotype_retained.(exists $gene->{nc_count} ? '_'.$gene->{nc_count} : ''));
      foreach my $t (@{$gene->get_all_Transcripts}) {
        compute_translation($t);
      }
    }
    $self->output(\@retained);
#    print STDERR "\n", ('-'x5), 'STEP X', "\n";
# Check for non canonical splice sites
    @for_step4 = grep {!exists $_->{retained}} @for_step4;
    print STDERR "\n", ('-'x5), 'FINAL STEP', "\n";
    my @final_step;
# Now I'm looking at my evidences to get the correct start and end of my transcript
    my $no_polyA_biotype = $self->param('biotype_noPolyA');
    my $polyA_biotype = $self->param('biotype_PolyA');
    my $NMD_biotype = $self->param('biotype_NMD');
    my $single_exon = $self->param('biotype_single');
    foreach my $gene ( sort {$b->{full_length} <=> $a->{full_length}} @for_step4) {
      if ($gene->{polyA_signal}) {
        $gene->biotype($polyA_biotype.(exists $gene->{nc_count} ? '_'.$gene->{nc_count} : ''));
      }
      else {
        $gene->biotype($no_polyA_biotype.(exists $gene->{nc_count} ? '_'.$gene->{nc_count} : ''));
      }
      foreach my $transcript (@{$gene->get_all_Transcripts}) {
        compute_translation($transcript);
        if (defined($transcript->translation)) {
          if ($transcript->end_Exon != $transcript->translation->end_Exon) {
            $gene->biotype($NMD_biotype.(exists $gene->{nc_count} ? '_'.$gene->{nc_count} : '')) if (($transcript->translation->end_Exon->length-$transcript->translation->end) > 50);
          }
        } else {
          if ( (scalar(@{$transcript->get_all_Exons}) <2) or ($transcript->length()<90) ) {
            print STDERR "single exon without translation or short (less than 90bp) multiexon \n";
          } else {
            $self->throw('Multiexon_without_translation'); 
          }
        }
      }
      push(@final_step, $gene);
    }
    $self->output(\@final_step);
  }
  print_Gene_list($self->output);
}

sub print_Gene_list {
  my ($genes) = @_;

  my $i = 0;
  my $count = 0;
  foreach my $gene (@$genes) {
    print 'Gene ', $i++, ' ', $gene->{full_length}, ' ', $gene->biotype, ' ', $gene->{polyA_signal}, "\n";
    foreach my $transcript (@{$gene->get_all_Transcripts}) {
      foreach my $tsf (@{$transcript->get_all_supporting_features}) {
        print "\t", join(' ', $tsf->hseqname, $tsf->seq_region_start, $tsf->seq_region_end, $tsf->seq_region_strand), "\n";
        ++$count;
      }
    }
    foreach my $cgene (@{$gene->{genes}}) {
      foreach my $transcript (@{$cgene->get_all_Transcripts}) {
        foreach my $tsf (@{$transcript->get_all_supporting_features}) {
          print "\t", join(' ', $tsf->hseqname, $tsf->seq_region_start, $tsf->seq_region_end, $tsf->seq_region_strand), "\n";
          ++$count;
        }
      }
    }
  }
  print STDERR "Worked on $count\n";
}


sub rank_polyA_signal {
  my ($transcript) = @_;

  if (has_polyA_signal($transcript)) {
    return 2;
  }
  else {
    return has_polyA_signal($transcript, 1);
  }
}


sub build_hashkey {
  my ($transcript, $type) = @_;

  if (!exists $transcript->{$type.'_hashkey'}) {
    my $arrayref;
    if ($type eq 'exon') {
      $arrayref = $transcript->get_all_Exons;
    }
    elsif ($type eq 'coding_exon') {
      $arrayref = $transcript->get_all_translateable_Exons;
    }
    elsif ($type eq 'intron') {
      $arrayref = $transcript->get_all_Introns;
    }
    my @hashkey;
    foreach my $element (@$arrayref) {
      push(@hashkey, join(':', $element->seq_region_start, $element->seq_region_end, $element->seq_region_strand));
    }
    $transcript->{$type.'_hashkey'} = join('%', @hashkey);
  }

  return $transcript->{$type.'_hashkey'};
}


sub write_output {
  my ($self) = @_;

  my $analysis = $self->analysis;
  my $ga = $self->hrdb_get_con('target_db')->get_GeneAdaptor;
  foreach my $gene (@{$self->output}) {
    empty_Gene($gene);
    attach_Analysis_to_Gene($gene, $analysis);
    $ga->store($gene);
  }
  return 1;
}

=head2 genes

 Arg [1]    : Arrayref Bio::EnsEMBL::Gene (optional)
 Description: Getter/setter for the genes to process
 Returntype : Arrayref Bio::EnsEMBL::Gene
 Exceptions : Throws if Arg[1] is not an array ref

=cut

sub genes {
  my ($self, $value) = @_;

  if ($value) {
    if (ref($value) eq 'ARRAY') {
      $self->param('genes', $value);
    }
    else {
      $self->throw('You should provide an arrayref');
    }
  }
  if ($self->param_is_defined('genes')) {
    return $self->param('genes');
  }
  else {
    return;
  }
}


=head2 input_dbs

 Arg [1]    : None
 Description: Getter for the 'input_dbs' to query for fetching the genes
 Returntype : Arrayref Hashref, databases connection details
 Exceptions : None

=cut

sub input_dbs {
  my ($self) = @_;

  if ($self->param_is_defined('source_dbs')) {
    my $dbs = $self->param('source_dbs');
    $dbs = [$dbs] unless (ref($dbs));
    return $dbs;
  }
  else {
    return;
  }
}


=head2 get_biotypes

 Arg [1]    : None
 Description: Getter for the biotypes to fetch. If it is empty, all genes will be fetched
 Returntype : Arrayref String
 Exceptions : None

=cut

sub get_biotypes {
  my ($self) = @_;

  if ($self->param_is_defined('biotypes')) {
    return $self->param('biotypes');
  }
  else {
    return;
  }
}


=head2 get_hashtypes

 Arg [1]    : None
 Description: Return the biotypes given in 'biotypes' as the 'good' set
 Returntype : Hashref of Array
 Exceptions : None

=cut

sub get_hashtypes {
  my ($self) = @_;

  return {good => $self->get_biotypes};
}

1;
