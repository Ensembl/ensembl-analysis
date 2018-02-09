=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils qw(cluster_Genes get_single_clusters);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(empty_Gene attach_Analysis_to_Gene);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(has_polyA_signal);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw(compute_translation);

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
    max_length_to_merge => 200,
    source_logic_name => undef,
    biotype_noPolyA => 'bad',
    biotype_PolyA => 'good',
    biotype_retained => 'retained',
    max_overlength => 10, # Random value might be around 20
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
    $db->dbc->disconnect_if_idle if ($self->param('disconnect_jobs'));
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
  print 'Clustering genes';
  my ($clusters, $unclustered) = cluster_Genes($self->param('genes'), $self->get_hashtypes);
  print " Done\n";
  print 'Working on unclustered genes';
  foreach my $cluster (@$unclustered) {
    foreach my $gene (@{$cluster->get_Genes}) {
      foreach my $transcript (@{$gene->get_all_Transcripts}) {
      $gene->{polyA_signal} = rank_polyA_signal($transcript);
      compute_translation($transcript);
      }
      $self->output([$gene]);
    }
  }
  print " Done\n";
  my $i = 0;
  my $t1;
  my $t2;
  my $max_length_merge;
  foreach my $cluster (@$clusters) {
    print STDERR ('-' x10), "\n", 'DEBUG CLUSTER ', join(' ', ++$i, $cluster->start, $cluster->end, $cluster->strand), "\n";
    my $genes = $cluster->get_Genes;
    my $num_genes = scalar(@$genes);
    print STDERR "Working on $num_genes\n";
    print STDERR ('-'x5), 'STEP 1', "\n";
# First we are looking at all the models to collapse the one which have the same intron structure
    my $index_gene = 0;
    my @to_process = sort {$b->length <=> $a->length} @$genes;
    my @for_step2;
    my $max_overlength = $self->param('max_overlength');
    for (my $index = 0; $index < @to_process; $index++) {
      my $good_gene = $to_process[$index];
      $good_gene->{full_length} = 1;
      next if (exists $good_gene->{processed});
      my $good_transcript = $good_gene->get_all_Transcripts->[0]; # We know that we only have 1 transcript per gene
      my $good_name = $good_transcript->get_all_supporting_features->[0]->hseqname;
      my $good_hashkey = build_hashkey($good_transcript, 'intron');
      $good_gene->{polyA_signal} = rank_polyA_signal($good_transcript);
      for (my $jndex = $index+1; $jndex < @to_process; $jndex++) {
        my $gene_to_process = $to_process[$jndex];
        my $transcript_to_process = $gene_to_process->get_all_Transcripts->[0]; # We know that we only have 1 transcript per gene
        my $hash_key_to_process = build_hashkey($transcript_to_process, 'intron');
        my $name_tp = $transcript_to_process->get_all_supporting_features->[0]->hseqname;
        if ($hash_key_to_process) { # If it is a single exon gene the hashkey will be empty
          if ($good_hashkey eq $hash_key_to_process) {
            if (abs($good_transcript->end_Exon->length-$transcript_to_process->end_Exon->length) < 15) {
              $gene_to_process->{processed} = 1; # Mark the gene as processed
              ++$good_gene->{full_length}; # We want to keep the number of full length support
              push(@{$good_gene->{genes}}, $gene_to_process);
            }
          }
          elsif ($good_hashkey =~ /$hash_key_to_process/) { # If the transcript is a fragment we want to add it to our model unless it might be a retained intron
            my $first_exon_to_process = $transcript_to_process->start_Exon;
            my $last_exon_to_process = $transcript_to_process->end_Exon;
            my $add_to_good_gene = 1;
            foreach my $good_exon (@{$good_transcript->get_all_Exons}) {
#              print STDERR join(' ', $good_name, $name_tp, $first_exon_to_process->end, '==', $good_exon->end,
#                  'and', $good_exon->start, '-', $first_exon_to_process->start, '>', $max_overlength,
#                 'or', $last_exon_to_process->start, '==', $good_exon->start,
#                  'and', $first_exon_to_process->end, '-', $good_exon->end, '>', $max_overlength), "\n";
              if ($first_exon_to_process->end == $good_exon->end
                  and $good_exon->start - $first_exon_to_process->start > $max_overlength
                  and $good_exon != $good_transcript->start_Exon) {
                $gene_to_process->{retained} = 1;
#                print STDERR 'Retained against ', $good_name, ' ', $good_gene->{full_length}, ' for ', $name_tp, "\n";
                $add_to_good_gene = 0;
                last;
              }
              if ($last_exon_to_process->start == $good_exon->start
                   and $first_exon_to_process->end - $good_exon->end > $max_overlength) {
                 $add_to_good_gene = 0;
                 last;
               }
            }
            if ($add_to_good_gene) {
              push(@{$good_gene->{genes}}, $gene_to_process);
              $gene_to_process->{processed} = 1; # Mark the gene as processed
            }
            else {
              push(@{$gene_to_process->{match}}, $good_gene); # Our models match except for the last exons, it can be used for scoring later
            }
          }
        }
        else {
# Here we process the single exon models
          my $good_last_exon = $good_transcript->end_Exon;
          if ($good_last_exon->start < $transcript_to_process->end
              and $good_last_exon->end > $transcript_to_process->start) {
            if ($good_last_exon->strand == 1) {
              if ($good_last_exon->start-$transcript_to_process->start < $max_overlength) {
                push(@{$good_gene->{genes}}, $gene_to_process);
                $gene_to_process->{processed} = 1; # Mark the gene as processed
              }
              else {
                $gene_to_process->{retained} = 1; # Mark the gene as retained, might be a bit harsh
              }
            }
            else {
              if ($transcript_to_process->end-$good_last_exon->end < $max_overlength) {
                push(@{$good_gene->{genes}}, $gene_to_process);
                $gene_to_process->{processed} = 1; # Mark the gene as processed
              }
              else {
                $gene_to_process->{retained} = 1; # Mark the gene as retained, might be a bit harsh
              }
            }
          }
        }
      }
      compute_translation($good_transcript);
      push(@for_step2, $good_gene);
    }
#    print_Gene_list(\@for_step2);
#    print STDERR "\n", ('-'x5), 'STEP 2', "\n";
# Now we are looking at the models to see if the difference in intron structure was caused by the error rate of the PacBio.
    # $$$$$$$$$$---------$$$$$$$
    # ######---------###########
    # ######---------###########
    # ######---------###########
    # ######---------###########
#    my @for_step3 = sort {$a->start <=> $b->start} @for_step2;
#    my @intron_to_change;
#    for (my $i = 0; $i < @for_step3-1; $i++) {
#      my $good_gene = $for_step3[$i];
#      next if (exists $good_gene->{processed} and $good_gene->{processed} == 2);
#      my $good_transcript = $good_gene->get_all_Transcripts->[0];
#      my $good_introns = $good_transcript->get_all_Introns;
#      for (my $j = $i+1; $j < @for_step3; $j++) {
#        my $gene_to_process = $for_step3[$j];
#        if ($gene_to_process->{full_length} == 1) {
#          my $transcript_to_process = $gene_to_process->get_all_Transcripts->[0];
#          my $introns_to_process = $transcript_to_process->get_all_Introns;
#          my $diff = 0;
#          my $same = 0;
#          my $good_intron_match = 0;
#          my $good_intron_missing = 0;
#          my $good_intron_mismatch = 0;
#          my $intron_to_process_match = 0;
#          my $intron_to_process_missing = 0;
#          my $intron_to_process_mismatch = 0;
#          my $corrected = 0;
#          for (my $k = 0; $k < @$good_introns; $k++) {
#            my $good_intron = $good_introns->[$k];
#            for (my $l = $k+$diff; $l < @$introns_to_process; $l++) {
#              my $intron_to_process = $introns_to_process->[$l];
#              if ($intron_to_process->start == $good_intron->start
#                  and $intron_to_process->end == $good_intron->end) {
#                ++$same;
#                ++$good_intron_match;
#                ++$intron_to_process_match;
#              }
#              elsif ($good_intron->end < $intron_to_process->start) {
#                ++$good_intron_missing;
#                ++$diff;
#              }
#              elsif ($intron_to_process->end < $good_intron->start) {
#                ++$good_intron_missing;
#                last;
#              }
#              elsif ($good_intron->length == $intron_to_process->length) {
#                if ($good_intron->prev_Exon->start == $intron_to_process->prev_Exon->start
#                    and $good_intron->next_Exon->end == $intron_to_process->next_Exon->end) {
#                  if ($good_intron->is_splice_canonical and $intron_to_process->is_splice_canonical) {
#                    # If both are canonical maybe it's just a different splice site, hard to guess
#                  }
#                  elsif ($good_intron->is_splice_canonical) {
## If it is canonical we can just add the evidence if the structure is the same in the rest of the transcript
#                    ++$corrected;
#                    push(@intron_to_change, [$good_gene, $good_intron, $gene_to_process, $intron_to_process]);
#                  }
#                  elsif ($intron_to_process->is_splice_canonical) {
## If it is canonical we can just add the evidence if the structure is the same in the rest of the transcript
## But it will be a little bit more complicated
#                    ++$corrected;
#                    push(@intron_to_change, [$gene_to_process, $intron_to_process, $good_gene, $good_intron]);
#                  }
#                  else {
#                    # I don't think this would happened or they would definitely be different
#                  }
#                }
#              }
#            }
#          }
#          if ($corrected and $corrected < 3) {
#
#          }
#        }
#      }
##      push(@for_step3, $good_gene);
#    }
#    print_Gene_list(\@for_step3);
    print STDERR "\n", ('-'x5), 'STEP 3', "\n";
# Now I'm looking at retained intron
# So I have to look at all possible intron then check if an exon is overlapping.
# If I have a exon overlapping this intron and it has the same boundaries as the surrounding
# exons, I have a retained exon
#    my @for_step4 = sort {$b->length <=> $a->length} @for_step3;
    my @for_step4 = sort {$b->length <=> $a->length} @for_step2;
    foreach my $good_gene (@for_step4) {
      my $good_transcript = $good_gene->get_all_Transcripts->[0];
      my $good_introns = $good_transcript->get_all_Introns;
      my $good_name = $good_transcript->get_all_supporting_features->[0]->hseqname;
      GENE: foreach my $gene_to_process (@for_step4) {
        next if ($good_gene == $gene_to_process or exists $gene_to_process->{retained});
        my $transcript_to_process = $gene_to_process->get_all_Transcripts->[0];
        my $name_tp = $transcript_to_process->get_all_supporting_features->[0]->hseqname;
        my $exons_to_process = $transcript_to_process->get_all_Exons;
        my $first_exon = $transcript_to_process->start_Exon;
        my $last_exon = $transcript_to_process->end_Exon;
        $last_exon = $transcript_to_process->start_Exon if ($transcript_to_process->strand == -1);
        my $exon_index = 0;
        foreach my $good_intron (@$good_introns) {
          foreach my $exon (@$exons_to_process) {
            if ($exon == $last_exon) {
#              print STDERR "Not looking, last exon\n";
            }
            elsif ($exon->end < $good_intron->start) {
              next;
            }
            elsif ($exon->start > $good_intron->end) {
              last;
            }
            elsif ($exon->start < $good_intron->end and $exon->end > $good_intron->start
                and $exon->start < $good_intron->start and $exon->end > $good_intron->end) {
#              print STDERR $exon->length, ' ', ($good_intron->length/3), ' ', $good_intron->is_splice_canonical, "\n" if ($exon == $first_exon);
              if ($exon == $first_exon) {
                if (!$good_intron->is_splice_canonical) {
                  print STDERR 'Go next as first exon and NC ', $good_name, ' for ', $name_tp, "\n";
                  next;
                }
                if ($exon->length < $good_intron->length/3) {
#                  print STDERR 'Go next as first exon NB ', $good_name, ' for ', $name_tp, "\n";
                  next;
                }
              }
              $gene_to_process->{retained} = 1;
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
      $gene->biotype($biotype_retained);
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
    foreach my $gene ( sort {$b->{full_length} <=> $a->{full_length}} @for_step4) {
      if ($gene->{polyA_signal}) {
        $gene->biotype($polyA_biotype);
      }
      else {
        $gene->biotype($no_polyA_biotype);
      }
      push(@final_step, $gene);
    }
    $self->output(\@final_step);
  }
#  print_Gene_list($self->output);
}

sub print_Gene_list {
  my ($genes) = @_;

  my $i = 0;
  my $count = 0;
  foreach my $gene (@$genes) {
    print 'Gene ', $i++, ' ', $gene->{full_length}, ' ', $gene->biotype, "\n";
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
