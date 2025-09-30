=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2024] EMBL-European Bioinformatics Institute

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

  Bio::EnsEMBL::Analysis::Tools::FilterBlastGenome

=head1 SYNOPSIS

  my $parser = Bio::EnsEMBL::Analysis::Tools::FilterBlastGenome->
  new(
      -regex => '^\w+\s+(\w+)'
      -query_type => 'pep',
      -database_type => 'dna',
      -threshold_type => 'PVALUE',
      -threshold => 0.01,
     );
 my @results = @{$parser->parse_results('blast.out')};

=head1 DESCRIPTION

This is a blast parser/filter that uses the protein sequence as the
query and the genome as the db.

=head1 CONTACT

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

=cut

package Bio::EnsEMBL::Analysis::Tools::FilterBlastGenome;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Analysis::Tools::FeatureFilterOnGenome;

use parent 'Bio::EnsEMBL::Analysis::Tools::FilterBPlite';


=head2 filter_hits

 Arg [1]    : Arrayref of Bio::EnsEMBL::Analysis::Tools::BPliteWrapper
 Description: Check the coverage of the full hits to remove low coverage
              alignments overlapping higher coverage alignement of the
              same sequence
 Returntype : Hashref of String, the name of the sequences to keep
 Exceptions : None

=cut

sub filter_hits {
  my ($self, $parsers) = @_;

  my %ids;
  my @features;
  foreach my $parser(@$parsers) {
    while(my $sbjct = $parser->nextSbjct) {
      my $name = $sbjct->name;
      while (my $hsp = $sbjct->nextHSP) {
        if($self->is_hsp_valid($hsp)) {
          my $qstart = $hsp->subject->start();
          my $hstart = $hsp->query->start();
          my $qend   = $hsp->subject->end();
          my $hend   = $hsp->query->end();
          my $qstrand = $hsp->subject->strand();
          my $hstrand = $hsp->query->strand();
          my $score  = $hsp->score;
          my $p_value = $hsp->P;
          my $percent = $hsp->percent;

          my $fp = $self->feature_factory->create_feature_pair
            ($qstart, $qend, $qstrand, $score, $hstart,
             $hend, $hstrand, $name, $percent, $p_value);


          push(@features,$fp);
        }
      }
    }
  }

  my $search = Bio::EnsEMBL::Analysis::Tools::FeatureFilterOnGenome->new(-coverage => $self->coverage);

  foreach my $f (@{$search->filter_results(\@features)}) {
    $ids{$f->hseqname} = 1;
  }
  return \%ids;
}



=head2 split_hsp

 Arg [1]    : Bio::EnsEMBL::Analysis::Tools::BPlite::HSP, representing the region of
              genome where Arg[2] aligns
 Arg [2]    : String, name of the protein or dna sequence
 Description: Process the BLAST hit of a protein or dna sequence onto the genome by
              retrieving the alignment information from Arg[1], create multiple Bio::EnsEMBL::FeaturePair
              when there are gaps in the alignment to finally create Bio::EnsEMBL::BaseAlignFeature to
              be able to store the alignment in a database. The query and the subject are swapped to store
              the genomic location correctly.
 Returntype : Bio::EnsEMBL::BaseAlignFeature, either Bio::EnsEMBL::DnaPepAlignFeature or Bio::EnsEMBL::DnaDnaAlignFeature
 Exceptions : Throws if database_type and query_type is not either 1 and 1 or 1 and 3
              Throws if start is greater than end

=cut

sub split_hsp {
  my ($self,$hsp,$name) = @_;

  my $qstrand = $hsp->query->strand;
  my $hstrand = $hsp->subject->strand;
  my ($qinc,   $hinc)    = $self->find_increments($qstrand,$hstrand);
  my @qchars = split(//,$hsp->querySeq);  # split alignment into array of chars
  my @hchars = split(//,$hsp->sbjctSeq);  # ditto for hit sequence
  my $qstart = $hsp->query->start(); # Start off the feature pair start
  my $hstart = $hsp->subject->start(); # ditto
  my $qend   = $hsp->query->start(); # Set the feature pair end also
  my $hend   = $hsp->subject->start(); # ditto
  if ($qstrand == -1) {
      $qstart = $hsp->query->end;
      $qend   = $hsp->query->end;
  }
  if ($hstrand == -1) {
    $hstart = $hsp->subject->end;
    $hend   = $hsp->subject->end;
  }

  my $count = 0; # counter for the bases in the alignment
  my $found = 0; # flag saying whether we have a feature pair

  my $query_seqname;
  if ($hsp->query->can('seq_id')) {
    $query_seqname = $hsp->query->seq_id;
  }
  else {
    $query_seqname = $hsp->query->seqname;
  }

  my @tmpf;
  while ($count <= $#qchars) {
    # We have hit an ungapped region.  Increase the query and hit
    #counters and flag that we have a feature pair.
    if ($qchars[$count] ne '-' &&
        $hchars[$count] ne '-') {
      $qend += $qinc;
      $hend += $hinc;

      $found = 1;
    }
    else {
    # We have hit a gapped region.  If the feature pair flag is set
    # ($found) then make a feature pair, store it and reset the start
    # and end variables.
      if ($found == 1) {
        # The genomic sequence is the hit, thus we need to use the hit data first
        my $fp = $self->convert_to_featurepair($hstart, $hend, $hstrand,
                                               $hinc, $qstart, $qend,
                                               $qstrand, $qinc, $query_seqname,
                                               $name,
                                               $hsp->score,
                                               $hsp->percent, $hsp->P,
                                               $hsp->positive,
                                               $hsp->match);
        push(@tmpf, $fp);
      }

      # We're in a gapped region.  We need to increment the sequence that
      # doesn't have the gap in it to keep the coordinates correct.
      # We also need to reset the current end coordinates.

      if ($qchars[$count] ne '-') {
        $qstart = $qend   + $qinc;
      }
      else {
        $qstart = $qend;
      }
      if ($hchars[$count] ne '-') {
        $hstart = $hend   + $hinc;
      }
      else {
        $hstart = $hend;
      }
      $qend = $qstart;
      $hend = $hstart;

      $found = 0;
    }
    $count++;
  }
  # Remember the last feature
  if ($found == 1) {
    # The genomic sequence is the hit, thus we need to use the hit data first
    my $fp = $self->convert_to_featurepair($hstart, $hend, $hstrand,
                                           $hinc, $qstart, $qend,
                                           $qstrand, $qinc, $query_seqname,
                                           $name,
                                           $hsp->score,
                                           $hsp->percent, $hsp->P,
                                           $hsp->positive,
                                           $hsp->match);
    push(@tmpf, $fp);
  }

  my $fp;
  $qinc = abs( $qinc );
  $hinc = abs( $hinc );

  if ($qinc == 1 && $hinc == 3) {
    $fp = Bio::EnsEMBL::DnaPepAlignFeature->new(-features => \@tmpf, -align_type => 'ensembl');
  }
  elsif ($qinc == 1 && $hinc == 1) {
    $fp = Bio::EnsEMBL::DnaDnaAlignFeature->new(-features => \@tmpf, -align_type => 'ensembl');
  }
  else {
    throw( "Hardcoded values wrong?? " );
  }

# for compara
  $fp->positive_matches($hsp->positive);
  $fp->identical_matches($hsp->match);
  if ($fp->hstart > $fp->hend) {
    throw("Failed start ".$fp->hstart." is greater than end ".$fp->hend." ".
        "for ".$fp->hseqname."\n");
  }
  return $fp;
}

1;
