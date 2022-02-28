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

Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveIsoseqs2Genes

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveIsoseqs2Genes;

use strict;
use warnings;

use parent ('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveExonerate2Genes');

sub param_defaults {
  my ($self) = @_;

  return {
    %{$self->SUPER::param_defaults},
    max_padding => 20000,
  }
}

sub run {
  my ($self) = @_;

  $self->throw("Can't run - no runnable objects") unless ($self->runnable);
  my @results;

  my %genome_slices;
  my $slice_adaptor = $self->hrdb_get_con('target_db')->get_SliceAdaptor;
  my $max_padding = $self->param('max_padding');
  foreach my $runnable (@{$self->runnable}) {
    $runnable->run;
    foreach my $result (@{$runnable->output}) {
      if (@{$result->get_all_Exons} > 1) {
        my $slice_id = $result->start_Exon->seqname;
        if (not exists $genome_slices{$slice_id}) {
          # assumes genome seqs were named in the Ensembl API Slice naming
          # convention, i.e. coord_syst:version:seq_reg_id:start:end:strand
          $genome_slices{$slice_id} = $slice_adaptor->fetch_by_name($slice_id);
        }
        my $slice = $genome_slices{$slice_id};

        foreach my $exon (@{$result->get_all_Exons}){
          $exon->slice($slice);
          foreach my $evi (@{$exon->get_all_supporting_features}){
            $evi->slice($slice);
            $evi->analysis($self->analysis);
          }
        }
        foreach my $evi (@{$result->get_all_supporting_features}) {
          $evi->slice($slice);
          $evi->analysis($self->analysis);
        }
        $result->slice($slice);
        my @list_introns = sort {$a->length <=> $b->length} @{$result->get_all_Introns};
        my $biggest_intron = $list_introns[-1];
        if ($biggest_intron->length > $list_introns[int(scalar(@list_introns)/2)-1]->length*2) {
          if (($biggest_intron->prev_Exon->length < 40 and $biggest_intron->prev_Exon->rank($result) < 3) or
              ($biggest_intron->next_Exon->length < 40 and $biggest_intron->next_Exon->rank($result) > (@list_introns-3))) {
            my $start = $result->end;
            my $end = $result->start;
            foreach my $exon (@{$result->get_all_Exons}) {
              if ($exon->length > 40) {
                $start = $exon->start if ($start > $exon->start);
                $end = $exon->end if ($end < $exon->end);
              }
            }
            if ($biggest_intron->length > $max_padding) {
              $start -= 20000;
              $end += 20000;
            }
            else {
              $start -= $biggest_intron->length;
              $end += $biggest_intron->length;
            }
            if ($start < $end) {
              $start = 1 if ($start < 1);
              $end = $slice->end if ($end > $slice->end);
              my $target_slice = $slice->sub_Slice($start, $end);
              my %parameters = %{$self->parameters_hash};
              if (not exists($parameters{-options}) and
                  defined $self->OPTIONS) {
                $parameters{-options} = $self->OPTIONS;
              }
              if (not exists($parameters{-coverage_by_aligned}) and
                  defined $self->COVERAGE_BY_ALIGNED) {
                $parameters{-coverage_by_aligned} = $self->COVERAGE_BY_ALIGNED;
              }

              my $biotypes_hash = $self->get_biotype();
              my $rerun = Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript->new(
                        -program  => $runnable->program,
                        -analysis => $runnable->analysis,
                        -query_type     => $runnable->query_type,
                        -annotation_file => $runnable->annotation_file,
                        -biotypes => $biotypes_hash,
                        -calculate_coverage_and_pid => $self->param('calculate_coverage_and_pid'),
                        %parameters,
                        );
              if ($runnable->query_seqs) {
                my $seqname = $result->get_all_supporting_features->[0]->hseqname;
                foreach my $seq (@{$runnable->query_seqs}) {
                  if ($seq->id eq $seqname) {
                    $rerun->query_seqs([$seq]);
                    last;
                  }
                }
              }
              $self->say_with_header($target_slice->name);
              $rerun->target_seqs([$target_slice]);
              $rerun->_verbose($self->debug) if ($self->debug);
              $rerun->run;
              my $result2 = $rerun->output->[0];
              if ($result2) {
                foreach my $exon (@{$result2->get_all_Exons}){
                  $exon->slice($target_slice);
                  $self->say_with_header($exon->seq_region_start.' '.$exon->seq_region_end);
                  foreach my $evi (@{$exon->get_all_supporting_features}){
                    $evi->slice($target_slice);
                    $evi->analysis($self->analysis);
                  }
                }
                foreach my $evi (@{$result2->get_all_supporting_features}) {
                  $evi->slice($target_slice);
                  $evi->analysis($self->analysis);
                }
                $result2->slice($target_slice);
                $result = $result2;
              }
              else {
                $self->warning('No result for second run '.$slice_id.' '.$rerun->query_seqs->[0]->id);
              }
            }
            else {
              $self->warning('All small exons '.$slice_id.' '.$result->get_all_supporting_features->[0]->hseqname);
            }
          }
        }
      }
      push(@results, $result);
    }
  }
  if ($self->USE_KILL_LIST) {
    unlink $self->filtered_query_file;
  }
  $self->say_with_header(@results.' possible genes have been found');
  if ($self->filter) {
    my $filtered_transcripts = $self->filter->filter_results(\@results);
    @results = @$filtered_transcripts;
  }
  $self->say_with_header(@results.' possible genes have been found');

  my @genes = $self->make_genes(@results);
  $self->say_with_header(@genes.' possible genes have been found');
  $self->param('output_genes',\@genes);
}

1;
