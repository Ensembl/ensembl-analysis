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

Bio::EnsEMBL::Analysis::Tools::PacBioTranscriptFilter

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Tools::PacBioTranscriptFilter;

use strict;
use warnings;

use parent ('Bio::EnsEMBL::Analysis::Tools::CdnaUpdateTranscriptFilter');

sub filter_results {
  my ($self, $transcripts) = @_;
  my @modified_transcripts;
  foreach my $transcript (@$transcripts ){
    my $real_strand = $self->_get_transcript_evidence_strand($transcript);
    if ($transcript->strand != $real_strand) {
      my $exons = $transcript->get_all_Exons;
      $transcript->flush_Exons();
      foreach my $exon (@$exons) {
        $exon->strand($real_strand);
        $transcript->add_Exon($exon);
      }
      $transcript->{_gb_flag} = 1;
    }
    push(@modified_transcripts, $transcript);
  }
  return $self->SUPER::filter_results(\@modified_transcripts);
}

sub _get_transcript_evidence_strand {
  my ($self,$tran) = @_;

  my ($sf) = @{$tran->get_all_supporting_features};

  if (!$sf) {
    ($sf) = @{$tran->get_all_Exons->[0]->get_all_supporting_features};
  }

  return $sf->strand*$sf->hstrand;
}

1;
