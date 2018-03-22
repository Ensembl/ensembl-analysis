=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Runnable::Pseudogene2x - Adaptation of Pseudogene calling code for 2x genomes

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Runnable::Pseudogene2x;

use warnings ;
use strict;
use vars qw(@ISA);


use Bio::EnsEMBL::Analysis::Runnable::Pseudogene;
use Bio::EnsEMBL::Analysis::Config::Pseudogene;

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::Pseudogene);

# Object preamble - inherits from Bio::EnsEMBL::Root;


sub transcript_evidence{

  my ($self,$transcript,$gene) =@_;
  my $repeat_blocks = $self->get_repeats($gene);
  my $results;
  my  @exons =  @{$transcript->get_all_Exons};
  @exons = sort {$a->start <=> $b->start} @exons;
  my $prev_exon = undef;
  my $total_intron_len = 0;
  my $covered_intron_len = 0;
  my $total_exon_len = 0;
  my $covered_exon_len = 0;
  my $n_real_intron = 0;
  my $n_frameshift_intron = 0;
  my $covered_introns = 0 ;
  my $covered_exons = 0 ;
  my $num_introns = 0;

  # ignore small terminal exons; they are mostly
  # alignment artefacts rather than indicators of pseudogenes
  while(@exons and $exons[0]->length <= $self->PS_FRAMESHIFT_INTRON_LENGTH()) {
    shift @exons;
  }
  while(@exons and $exons[-1]->length <= $self->PS_FRAMESHIFT_INTRON_LENGTH()) {
    pop @exons;
  }

  #print STDERR "TRANSCRIPT ", $transcript->dbID, "\n";

  foreach my $exon (@exons) {
    # need to check whether the exon like on a gap or not
    my $exon_in_gap = 0;

    if ($self->seqlevel) {
      my $overlaps = 0;
      foreach my $c (@{$self->seqlevel}) {
        if ($exon->start <= $c->end and
            $exon->end >= $c->start) {
          $overlaps = 1;
          last;
        }
      }
      $exon_in_gap = 1 if not $overlaps;
    }

    if (defined($prev_exon)) {
      my $intron = Bio::EnsEMBL::Feature->new(
                                              -START => $prev_exon->end+1,
                                              -END => $exon->start-1,
                                              -STRAND => $exon->strand
                                              );
      # ignore the intron if it lies completely in a sequence gap;
      my $intron_in_gap = 0;
      if ($self->seqlevel) {
        my $overlaps = 0;
        foreach my $c (@{$self->seqlevel}) {
          if ($intron->start <= $c->end and
              $intron->end >= $c->start) {
            $overlaps = 1;
            last;
          }
        }
        $intron_in_gap = 1 if not $overlaps;
      }

      if (not $intron_in_gap) {
        if ($intron->length > $self->PS_FRAMESHIFT_INTRON_LENGTH()) {
          $n_real_intron++;
          # need to restrict ourselves to parts of intron that are
          # on real sequence
          my @seq_intron_bits = $self->_intersect_list($intron, @{$self->seqlevel});
          foreach my $bit (@seq_intron_bits) {          
            $total_intron_len+=$bit->length;
            $covered_intron_len+=$self->_len_covered($bit,$repeat_blocks);
          }
        } else {
          $n_frameshift_intron++;
        }

        $num_introns++;
        #print STDERR "Assessed REAL intron ", $intron->length, "\n";
      } else {
        #print STDERR "Skipping GAP intron ", $intron->start, " ", $intron->end, "\n";
      }
    }
    
    if (not $exon_in_gap) {
      my $seq_feature_exon = Bio::EnsEMBL::Feature->new(
                                                        -START => $exon->start,
                                                        -END => $exon->end,
                                                        -STRAND => $exon->strand
							);

      my @seq_exon_bits = $self->_intersect_list($seq_feature_exon, @{$self->seqlevel});
      foreach my $bit (@seq_exon_bits) {
        $total_exon_len+=$bit->length;
        $covered_exon_len+=$self->_len_covered($bit,$repeat_blocks);
      }
      
      #print STDERR "Assessed REAL exon\n";
    } else {
      #print STDERR "Skipping GAP exon\n";
    }

    $prev_exon = $exon;
  }

  #calculate percentage coverage
  
  if ($total_intron_len > 0) {
    $covered_introns = (($covered_intron_len/$total_intron_len)*100);
  }
  if ($total_exon_len >  0) {
    $covered_exons =   (($covered_exon_len/$total_exon_len)*100)
  }
  
  #print "COVERED INTRONS = $covered_introns\n";
  #print "COVERED EXONS = $covered_exons\n";

  $results = {
	      'num_introns' => $num_introns,
	      'total_exon_len' =>  $total_exon_len,
	      'covered_exons' => $covered_exons,
	      'total_intron_len' =>  $total_intron_len,
	      'covered_introns' => $covered_introns,
	      'real_introns' => $n_real_intron,
	      'frameshift_introns' => $n_frameshift_intron
	     };
  return $results;
}


sub _intersect_list {
  my ($self, $f1, @blocks) = @_;

  my @ret;

  foreach my $f (@blocks) {
    if ($f->end < $f1->start) {
      next;
    } elsif ($f->start > $f1->end) {
      next;
    } else {
      my $start = $f1->start;
      my $end   = $f1->end;
      $start = $f->start if $f->start > $start;
      $end   = $f->end if $f->end < $end;

      push @ret, Bio::EnsEMBL::Feature->new(-start => $start,
                                            -end   => $end,
                                            -strand => $f1->strand);
    }
  }
  
  @ret = sort { $a->start <=> $b->start } @ret;

  return @ret;
}

sub seqlevel {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_seqlevel} = $val;
  }

  return $self->{_seqlevel};
}


return 1;
