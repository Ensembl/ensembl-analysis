=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Analysis::Runnable::SelenocysteineFinder

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Runnable::SelenocysteineFinder;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Tools::Utilities qw(locate_executable);
use Bio::EnsEMBL::Utils::Exception qw(throw);

use parent qw(Bio::EnsEMBL::Analysis::Runnable);

sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($genewise, $minimum_identity, $score_threshold, $coverage_threshold) =
      rearrange([qw(
                    GENEWISE
                    MINIMUM_IDENTITY
                    SCORE_THRESHOLD
                    COVERAGE_THRESHOLD
                    )
                 ], @args);
  $self->genewise($genewise || 'genewise');
  $self->score_threshold($score_threshold || 2);
  $self->coverage_threshold($coverage_threshold || 90);
  $self->minimum_identity($minimum_identity || 95);

  return $self;
}

sub run {
  my ($self) = @_;

    $self->runnable(Bio::EnsEMBL::Analysis::Runnable::ExonerateTranscript->new(
      -program => $exonerate,
      -analysis => $analysis,
      -query_seqs => $querys,
      -query_type => 'protein',
      -target_file => $target,
      -options => '--model protein2genome --bestn 3',
    ));
  my ($self) = @_;

  my $runnable = $self->runnable;
  $runnable->run;
  my $slice_adaptor = $self->hrdb_get_con('target_db')->get_SliceAdaptor;
  my %alignments;
  my $minimum_identity = $self->param('minimum_identity');
  my $coverage_threshold = $self->param('coverage_threshold')/100;
  my $score_threshold = $self->param('score_threshold')/100;
  foreach my $transcript (@{$runnable->clean_output}) {
    if (exists $alignments{$transcript->{accession}}) {
        my $accession = $transcript->{accession};
        my $tsf = $transcript->get_all_supporting_features->[0];
        if ($tsf->score > ($alignments{$accession}->{score}+$tsf->{score}*$score_threshold)
          or ($tsf->percent_id > $minimum_identity and              
              $tsf->percent_id > $alignments{$accession}->get_all_supporting_features->[0]->percent_id and
              (($tsf->q_end-$tsf->q_start) > $tsf->{query_length}*$coverage_threshold)) ) {
            $alignments{$accession} = {
                t_seqname => $t_seqname,
                t_start => $t_start,
                t_end => $t_end,
                t_strand => $t_strand,
                score => $score,
                identity => $perc_ident,
                };
        }
    }
    else {
        $alignments{$transcript->{accession}} = $transcript;
    }
  }
  foreach my $transcript (values %alignments) {
    my $transcripts = $self->selenocysteine_finder($alignments, $seq, $slice_adaptor);
    my $genewise_sequences = $self->mutate_genomic_sequence($transcripts, $slice_adaptor);

    my $selenocysteine_transcripts = $self->genewise_sequences($genewise_sequences, $seq, $alignments);
    $self->output($selenocysteine_transcripts);
  }
}

sub genewise {
  my ($self, $binary) = @_;

  if ($binary) {
    $self->{_genewise} = locate_executable($binary, $self->bindir);
  }
  return $self->{_genewise};
}


sub minimum_identity {
  my ($self, $value) = @_;

  if ($value) {
    $self->{_minimum_identity} = $value;
  }
  return $self->{_minimum_identity};
}


sub score_threshold {
  my ($self, $value) = @_;

  if ($value) {
    $self->{_score_threshold} = $value;
  }
  return $self->{_score_threshold};
}


sub coverage_threshold {
  my ($self, $value) = @_;

  if ($value) {
    $self->{_coverage_threshold} = $value;
  }
  return $self->{_coverage_threshold};
}

1;
