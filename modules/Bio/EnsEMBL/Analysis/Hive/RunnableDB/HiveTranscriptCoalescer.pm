=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016] EMBL-European Bioinformatics Institute

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

  my $dna_db = $self->get_database_by_name('dna_db');
  my @genes;
  foreach my $input_db (@{$self->input_dbs}) {
    my $db = $self->hrdb_get_dba($input_db, $dna_db);
    my $slice = $self->fetch_sequence($self->input_id, $db);
    if (scalar(@{$self->get_biotypes})) {
      foreach my $biotype (@{$self->get_biotypes}) {
        push(@genes, @{$slice->get_all_Genes_by_type($biotype)});
      }
    }
    else {
      push(@genes, @{$slice->get_all_Genes});
    }
  }
  print STDERR scalar(@genes), " genes to process\n";
  if (scalar(@genes)) {
    $self->genes(\@genes);
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

  my @output;
  my ($clusters, $unclustered) = cluster_Genes($self->genes, $self->get_hashtypes);
  foreach my $cluster (@$unclustered) {
    push(@output,  @{$cluster->get_Genes});
  }
  my $i = 0;
  my $t1;
  my $t2;
  my $max_length_merge;
  foreach my $cluster (@$clusters) {
    print STDERR ('-' x10), "\n", 'DEBUG CLUSTER ', join(' ', ++$i, $cluster->start, $cluster->end, $cluster->strand), "\n";
    my $genes = $cluster->get_Genes;
    my $num_genes = scalar(@$genes);
    my $index_gene = 0;
    my %processed;
    for ($index_gene = 0; $index_gene < $num_genes-1; $index_gene++) {
      print STDERR " ====\t", join(' ', $genes->[$index_gene]->seq_region_start, $genes->[$index_gene]->seq_region_end, $genes->[$index_gene]->seq_region_strand, $genes->[$index_gene]->get_all_Transcripts->[0]->get_all_supporting_features->[0]->hseqname), "\n";
      my @merge;
      my @possible_merge;
      my %can_be_merge;
      my $k = 1;
      my %retained;
      my $check_retain = 0;
      $max_length_merge = $self->param('max_length_to_merge');
      for ($k = $index_gene+1; $k < $num_genes; $k++) {
        if (exists $processed{$k}) {
          print STDERR '  ', '+' x10, "\n", '  ALREADY DONE ', join(' ', $k, $genes->[$k]->seq_region_start, $genes->[$k]->seq_region_end, $genes->[$k]->seq_region_strand, $genes->[$k]->get_all_Transcripts->[0]->get_all_supporting_features->[0]->hseqname), "\n";
        }
        else {
          print STDERR '  ', '*' x10, "\n", '  WORKING ON ', join(' ', $k, $genes->[$k]->seq_region_start, $genes->[$k]->seq_region_end, $genes->[$k]->seq_region_strand, $genes->[$k]->get_all_Transcripts->[0]->get_all_supporting_features->[0]->hseqname), "\n";
          $t1 = $genes->[$index_gene]->get_all_Transcripts->[0];
          $t2 = $genes->[$k]->get_all_Transcripts->[0];
          my $ei1 = 0;
          my $ei2 = 0;
          my $exons1 = $t1->get_all_Exons;
          my $exons2 = $t2->get_all_Exons;
          my $introns1 = $t1->get_all_Introns;
          my $introns2 = $t2->get_all_Introns;
          my $count_introns = 0;
          my $check_splice = 0;
          my $can_first_be_longer = 0;
          my $can_last_be_longer = 0;
          if ($exons1->[0]->seq_region_end >= $exons2->[0]->seq_region_start
            and $exons1->[0]->seq_region_start <= $exons2->[0]->seq_region_end) {
            if ($exons1->[0]->seq_region_start < $exons2->[0]->seq_region_start) {
              if (($exons2->[0]->seq_region_start-$exons1->[0]->seq_region_start) < $max_length_merge) {
                $can_first_be_longer = 1;
              }
              else {
                print STDERR '    DEBUG S length is ', ($exons2->[0]->seq_region_start-$exons1->[0]->seq_region_start), ' for ', $exons1->[0]->seq_region_start, ' and ', $exons2->[0]->seq_region_start, "\n";
                $can_first_be_longer = 2;
              }
              
            }
          }
          if ($exons1->[-1]->seq_region_end >= $exons2->[-1]->seq_region_start
            and $exons1->[-1]->seq_region_start <= $exons2->[-1]->seq_region_end) {
            if ($exons1->[-1]->seq_region_end > $exons2->[-1]->seq_region_end) {
              if (($exons1->[-1]->seq_region_end-$exons2->[-1]->seq_region_end) < $max_length_merge) {
                $can_last_be_longer = 1;
              }
              else {
                print STDERR '    DEBUG E length is ', ($exons1->[-1]->seq_region_end-$exons2->[-1]->seq_region_end), ' for ', $exons2->[-1]->seq_region_end, ' and ', $exons1->[-1]->seq_region_end, "\n";
                $can_last_be_longer = 2;
              }
            }
          }
          for ($ei1 = 0; $ei1 < scalar(@$exons1); $ei1++) {
            my $exon1 = $exons1->[$ei1];
            my $intron1 = $introns1->[$ei1];
            for ($ei2 = 0; $ei2 < scalar(@$exons2); $ei2++) {
              my $exon2 = $exons2->[$ei2];
              my $intron2 = $introns2->[$ei2];
              # First let's check the first and end exon to see if they are similar
              if ($exon1->seq_region_end < $exon2->seq_region_start) {
                last;
              }
              elsif ($exon1->seq_region_end >= $exon2->seq_region_start
                and $exon1->seq_region_start <= $exon2->seq_region_end) {
                if ($intron1) {
                  if ($intron2) {
                    my $iid1 = $intron1->seq_region_start.':'.$intron1->seq_region_end.':'.$intron1->seq_region_strand;
                    my $iid2 = $intron2->seq_region_start.':'.$intron2->seq_region_end.':'.$intron2->seq_region_strand;
                    if ($ei1 > 0 and $ei2 == 0) {
                      if (($exon2->seq_region_end - $exon1->seq_region_end) < 10) {
                        ++$check_retain;
                        $retained{start}{$k} = {iid => $intron2->seq_region_start.':'.$intron2->seq_region_end, ee1 => $exon1->seq_region_end, start => $intron2->seq_region_start, end => $intron2->seq_region_end};
                        print STDERR '    DEBUG RETAINED S ', join(' ', $exon2->seq_region_end, $introns1->[$ei1-1]->seq_region_start, $introns1->[$ei1-1]->seq_region_end) , "\n";
                      }
                      else {
                        print STDERR '    DEBUG NO RETAINED S ', join(' ', $exon2->seq_region_end, $introns1->[$ei1-1]->seq_region_start, $introns1->[$ei1-1]->seq_region_end) , "\n";
                      }
                    }
                    if ($intron1->seq_region_start == $intron2->seq_region_start
                      and $intron1->seq_region_end == $intron2->seq_region_end) {
                      ++$count_introns;
                      print STDERR '    DEBUG same Intron for ', $ei1, ' and ', $ei2, "\n";
                      $can_be_merge{$iid1}->{canonical} = $intron1->is_splice_canonical;
                      $can_be_merge{$iid1}->{score}++;
                      push(@{$can_be_merge{$iid2}->{elm}}, $k);
                    }
                    else {
                      print STDERR '    DEBUG Intron differ for ', $ei1, ' and ', $ei2, "\n";
                      if ($intron1->length == $intron2->length
                        and (($intron1->seq_region_start-$intron2->seq_region_start) == ($intron1->seq_region_end-$intron2->seq_region_end))) {
                        if (abs($intron1->seq_region_start-$intron2->seq_region_start) < 10) {
                          print STDERR '    DEBUG Intron differ but are close for ', $ei1, ' and ', $ei2, "\n";
                          ++$check_splice;
#                          $can_be_merge{$k} = { ei1 => $ei1, ei2 => $ei2, i1 => $intron1->seq_region_start.':'.$intron1->seq_region_end.':'.$intron1->seq_region_strand, i2 => $intron2->seq_region_start.':'.$intron2->seq_region_end.':'.$intron2->seq_region_strand, e1can => $intron1->is_splice_canonical, e2can => $intron2->is_splice_canonical};
                          $can_be_merge{$iid1}->{canonical} = $intron1->is_splice_canonical;
                          $can_be_merge{$iid1}->{score}++;
                          push(@{$can_be_merge{$iid1}->{elm}}, $index_gene);
                          push(@{$can_be_merge{$iid1}->{other}}, $iid2);
                          push(@{$can_be_merge{$iid2}->{other}}, $iid1);
                          $can_be_merge{$iid2}->{canonical} = $intron2->is_splice_canonical;
                          $can_be_merge{$iid2}->{score}++;
                          push(@{$can_be_merge{$iid2}->{elm}}, $k);
                          ++$count_introns;
                        }
                        else {
                          print STDERR '    DEBUG diff is too long', abs($intron1->seq_region_start-$intron2->seq_region_start), "\n";
                        }
                      }
                    }
                  }
                }
                else {
                  if ($intron2 and $exon1->seq_region_end > $exon2->seq_region_start) {
                    $retained{end}{$k} = {iid => $intron2->seq_region_start.':'.$intron2->seq_region_end, ee1 => $exon1->seq_region_end, start => $intron2->seq_region_start, end => $intron2->seq_region_end};
                    ++$check_retain;
                    print STDERR '    DEBUG RETAINED E ', join(' ', $exon1->seq_region_start, $intron2->seq_region_start, $intron2->seq_region_end) , "\n";
                  }
                  else {
                    print STDERR '    DEBUG NO RETAINED E ', join(' ', $exon2->seq_region_start, $introns1->[$ei1-1]->seq_region_start, $introns1->[$ei1-1]->seq_region_end) , "\n";
                  }
                }
              }
            }
          }
          if ($count_introns == scalar(@$introns1) and $can_first_be_longer < 2 and $can_last_be_longer < 2) {
            if ($check_splice) {
              push(@possible_merge, $k);
              print STDERR "   SPLIT $index_gene and $k\n";
            }
            else {
              $processed{$k} = 1;
              push(@merge, $k);
              print STDERR "   MERGE $index_gene and $k\n";
            }
          }
          elsif ($check_retain and $count_introns == scalar(@$introns2)) {
            push(@possible_merge, $k);
            print STDERR "   REMER $index_gene and $k\n";
          }
          else {
            print STDERR "   LEAVE $index_gene and $k, $count_introns ".@$introns1.", $can_first_be_longer, $can_last_be_longer\n";
          }
        }
      }
      if ($check_retain) {
        my %score;
        my %good_to_merge;
        foreach my $key ('start', 'end') {
          my $num_keys = 0;
          foreach my $index (keys %{$retained{$key}}) {
            ++$score{$retained{$key}{$index}->{iid}};
            ++$num_keys;
            if (exists $good_to_merge{$index}) {
              ++$good_to_merge{$index};
            }
            else{
              $good_to_merge{$index} = 0;
            }
          }
          foreach my $scorekey (keys %score) {
            if ($score{$scorekey} == 1 or $score{$scorekey} < ($num_keys*0.1)) {
              print STDERR '    NOT RETAINED ', $scorekey, "\n";
            }
            else {
              print STDERR '    POSSIBLE RETAINED ', $scorekey, "\n";
            }
          }
        }
      }
        my @keys;
      my %tmp;
      foreach my $index (@merge) {
        foreach my $key (keys %can_be_merge) {
          my $best_key = $key;
          print STDERR '     DEBUG: checking ', $key, "\n";
          my $count_other = 0;
          foreach my $other (@{$can_be_merge{$key}->{other}}) {
            print STDERR '     DEBUG: checking other ', $other, "\n";
            ++$other;
            if ($can_be_merge{$best_key}->{score} > $can_be_merge{$other}->{score}) {
              if (!$can_be_merge{$best_key}->{canonical}) {
                if (($can_be_merge{$best_key}->{score}/2) < $can_be_merge{$other}->{score}) {
                  $best_key = $other;
          print STDERR '     DEBUG: best is now ', $best_key, "\n";
                }
              }
            }
            else {
              if ($can_be_merge{$other}->{canonical}) {
                $best_key = $other unless ($can_be_merge{$best_key}->{canonical});
          print STDERR '     DEBUG: best is now ', $best_key, "\n";
              }
            }
          }
#          push(@keys, @{$can_be_merge{$best_key}->{elm}});
          push(@keys, $index);
        }
      }
      foreach my $index (@possible_merge) {

      }
      %tmp = map { $_ => $_ } @keys;
      foreach my $k (values %tmp) {
        print STDERR '   MERGING: ', join(' ', $k, $genes->[$k]->seq_region_start, $genes->[$k]->seq_region_end, $genes->[$k]->seq_region_strand, $genes->[$k]->get_all_Transcripts->[0]->get_all_supporting_features->[0]->hseqname), "\n";
      }
    }
  }
}

sub write_output {
  my ($self) = @_;
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

  if ($self->param_is_defined('input_dbs')) {
    return $self->param('input_dbs');
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
