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

Bio::EnsEMBL::Analysis::Tools::RepeatCoverage

=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Analysis::Tools::RepeatCoverage;

use strict;
use warnings;

use Exporter qw(import);
use Bio::EnsEMBL::Feature;

our @EXPORT_OK = qw(
  get_genes
  get_repeat_blocks
  get_masked_stats
);


=head2 get_genes

 Arg [1]    : Bio::EnsEMBL::Slice from a database with genes
 Arg [2]    : Arrayref of String representing the biotype of genes to fetch
 Description: Fetch all the genes from a region having the biotypes from Arg[2]
 Returntype : Arrayref of Bio::EnsEMBL::Gene
 Exceptions : None

=cut

sub get_genes {
  my ($geneslice, $genetypes) = @_;

  my @genes;
  print "Fetching genes\n";

  foreach my $genetype (@$genetypes) {
    my $genes_by_type = $geneslice->get_all_Genes_by_type($genetype);
    print "Got " . scalar(@$genes_by_type) . " $genetype genes\n";
    push(@genes,@$genes_by_type);
  }
  print "Done fetching genes (fetched " . scalar(@genes) .")\n";
  @genes = sort {$a->start <=> $b->start} @genes;
  return \@genes;
}


=head2 get_repeat_blocks

 Arg [1]    : Arrayref of Bio::EnsEMBL::Feature representing the repeats of each analyses
 Description: Collapse the features from Arg[1] to create non overlapping features
 Returntype : Arrayref of Bio::EnsEMBL::SeqFeature
 Exceptions : None

=cut

sub get_repeat_blocks {
  my ($repeats) = @_;

  my @repeat_blocks;
  my $curblock = undef;
  foreach my $repeat (sort {$a->start <=> $b->start} @$repeats) {
    if ($repeat->start <= 0) {
      $repeat->start(1);
    }
    if (defined($curblock) && $curblock->end >= $repeat->start) {
      if ($repeat->end > $curblock->end) {
        $curblock->end($repeat->end);
      }
    }
    else {
      $curblock = Bio::EnsEMBL::Feature->new(-START => $repeat->start, -END => $repeat->end);
      push(@repeat_blocks, $curblock);
    }
  }
  return \@repeat_blocks;
}


=head2 get_masked_stats

 Arg [1]    : Arrayref of Bio::EnsEMBL::Feature representing the repeats
 Arg [2]    : Arrayref of Bio::EnsEMBL::Gene
 Description: Calculate the number of erpeat masked bases.
              If Arg[2] is defined, it calculates the number of genes for which
              at least one exon overlaps a repeat feature and print the supporting
              evidences used for the gene.
 Returntype : Array of Int, number of bases masked, number of genes overlapped
 Exceptions : None

=cut

sub get_masked_stats {
  my ($repeat_blocks, $genes) = @_;

  my $novergene = 0;
  my $nmasked = 0;
  foreach my $block (@$repeat_blocks) {
    $nmasked += $block->length;
    if ($genes) {
      my $ngene = scalar(@$genes);
      GENE: for (my $ind=0; $ind < $ngene; $ind++) {
        my $gene = $genes->[$ind];
        if (!defined($gene)) {
          next GENE;
        }
        if ($block->end >= $gene->start &&
            $block->start <= $gene->end) {
          foreach my $exon (@{$gene->get_all_Exons}) {
            if ($exon->overlaps($block)) {
              print 'Overlap for gene '.$gene->display_id."\n";
              foreach my $trans (@{$gene->get_all_Transcripts}) {
                foreach my $tsf (@{$trans->get_all_supporting_features}) {
                  print ' Support '.$tsf->analysis->logic_name.' '.$tsf->hseqname."\n";
                }
              }
              ++$novergene;
              $genes->[$ind] = undef;
              next GENE;
            }
          }
        }
        elsif ($gene->start > $block->end) {
          last;
        }
        elsif ($gene->end < $block->start) {
          $genes->[$ind] = undef;
        }
      }
      print "Number genes overlapped = $novergene\n";
    }
  }
  print "Number bases masked = $nmasked\n";
  return $nmasked, $novergene;
}

1;
