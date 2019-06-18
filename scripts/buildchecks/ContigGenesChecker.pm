# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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
#
# EnsEMBL module for ContigGenesChecker
#
# Cared for by Steve Searle <searle@sanger.ac.uk>
#

# POD documentation - main docs before the code

=pod

=head1 NAME

Bio::EnsEMBL::Test::ContigGenesChecker


=head1 SYNOPSIS
Module to check the validity of a transcript

=head1 DESCRIPTION
Performs various checks on the Genes in a Bio::EnsEMBL::Contig object. These
should only be checks which are not internal to a particular Gene but instead
dependant on the arrangement of genes on the Contig. These include:
  1. Genes have been clustered correctly
  2. Interlocking Genes
  3. Genes on both strands with overlapping exons.
  4. Transcripts on both strands with all overlapping exons

This class does not use stable_ids but instead dbIDs because stable_ids are
not set until after the gene build.

=head1 CONTACT

  Steve Searle <searle@sanger.ac.uk>

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package ContigGenesChecker;
use vars qw(@ISA);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use strict;

use Checker;
use TranscriptCluster;



@ISA = qw( Checker );


sub new {
  my ($class, @args) = @_;

  my $self = bless {},$class;

  my ( $genes, $ignorewarnings, $adaptor,
       $slice  ) = rearrange
	 ( [ qw ( GENES
                  IGNOREWARNINGS
                  ADAPTOR
                  SLICE
	      )], @args );


  if( !defined $genes ) {
    $self->throw("Genes must be set in new for ContigGenesChecker")
  }
  $self->genes($genes);

  if( defined $ignorewarnings) { $self->ignorewarnings($ignorewarnings); }
  if( defined $slice) { $self->slice($slice); }
  if( defined $adaptor ) { $self->adaptor( $adaptor )}

  $self->{_errors} = [];
  $self->{_warnings} = [];

  return $self;
}





sub slice {
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      if (!($value->isa('Bio::EnsEMBL::Slice'))) {
         $self->throw("slice passed a non Bio::EnsEMBL::Slice object\n");
      }
      $self->{_slice} = $value;
    }
    return $self->{_slice};

}

sub adaptor {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{_adaptor} = $arg;
  }
  return $self->{_adaptor};
}

sub genes {
  my ( $self, $arg ) = @_;
  if( defined $arg ) {
    $self->{_genes} = $arg;
  }
  return $self->{_genes};
}


sub output {
  my $self = shift;

  $self->SUPER::output;
}


sub check {
  my $self = shift;


  print "Checking clustering\n";
  my $corrected_clusters = $self->check_Clustering();

  print "Looking for strand overlaps\n";
  $self->find_StrandOverlaps($corrected_clusters);

  print "Looking for interlocks\n";
  $self->find_Interlocks($corrected_clusters);
}

sub check_Clustering {
  my $self = shift;

  my @transcripts_unsorted;
# reject non-translating transcripts
  my $genes = $self->{'_genes'};
  foreach my $gene (@$genes) {
    foreach my $tran (@{$gene->get_all_Transcripts}) {
      push(@transcripts_unsorted, $tran);
    }
  }

  my @transcripts = sort by_transcript_high @transcripts_unsorted;

  my $ga = $self->adaptor->get_adaptor("Gene");

  my @clusters;
# clusters transcripts by whether or not any exon overlaps with an exon in
# another transcript (came from prune in GeneBuilder)
  foreach my $tran (@transcripts) {
    my @matching_clusters;
    my ($trans_start, $trans_end) = ($tran->start,$tran->end);

    # print "transcript limits: $trans_start $trans_end \n";

    CLUSTER: foreach my $cluster (@clusters) {

      #print "Testing against cluster with limits " . $cluster->start .
      #      " to " . $cluster->end . "\n";

      if (!($trans_start > $cluster->end ||
            $trans_end < $cluster->start)) {
        # print "In range\n";
        foreach my $cluster_transcript (@{$cluster->get_Transcripts()}) {
          foreach my $exon1 (@{$tran->get_all_Exons}) {

            foreach my $cluster_exon (@{$cluster_transcript->get_all_Exons}) {
              if ($exon1->overlaps($cluster_exon) &&
                  $exon1->strand == $cluster_exon->strand) {
                push (@matching_clusters, $cluster);
                next CLUSTER;
              }
            }
          }
        }
      }
    }

    if (scalar(@matching_clusters) == 0) {
      # print STDERR "Found new cluster for " . $tran->dbID . "\n";
      my $newcluster = new TranscriptCluster;
      $newcluster->put_Transcripts($tran);
      push(@clusters,$newcluster);

    } elsif (scalar(@matching_clusters) == 1) {
       # print STDERR "Adding to cluster for " . $tran->dbID . "\n";
      $matching_clusters[0]->put_Transcripts($tran);

    } else {

# Merge the matching clusters into a single cluster
      # print STDERR "Merging clusters for " . $tran->dbID . "\n";
      my @new_clusters;
      my $merged_cluster = new TranscriptCluster;
      foreach my $clust (@matching_clusters) {
        $merged_cluster->put_Transcripts(@{$clust->get_Transcripts});
      }
      $merged_cluster->put_Transcripts($tran);
      push @new_clusters,$merged_cluster;

# Add back non matching clusters
      foreach my $clust (@clusters) {
        my $found = 0;
        MATCHING: foreach my $m_clust (@matching_clusters) {
          if ($clust == $m_clust) {
            $found = 1;
            last MATCHING;
          }
        }
        if (!$found) {
          push @new_clusters,$clust;
        }
      }
      @clusters = @new_clusters;
    }
  }

#Safety checks
  my $ntrans = 0;
  my %trans_check_hash;
  foreach my $cluster (@clusters) {

    my %gene_ids = ();
    $ntrans += scalar(@{$cluster->get_Transcripts});
    foreach my $trans (@{$cluster->get_Transcripts}) {
      if (defined($trans_check_hash{"$trans"})) {
        $self->throw("Transcript " . $trans->dbID . " added twice to clusters\n");
      }
      $trans_check_hash{"$trans"} = 1;

      #recording gene ids
      $gene_ids{$ga->fetch_by_transcript_id($trans->dbID)->dbID} = 1;

    }
    if (!scalar(@{$cluster->get_Transcripts})) {
      $self->throw("Empty cluster");
    }

    if(scalar keys %gene_ids > 1){
      print STDERR "Multiple gene cluster: ".(join(", ", keys %gene_ids))."\n";
    }
  }
  if ($ntrans != scalar(@transcripts)) {
    $self->throw("Not all transcripts have been added into clusters $ntrans and " . scalar(@transcripts). " \n");
  }
#end safety checks

  if (scalar(@clusters) < scalar(@$genes)) {
    $self->add_Error("Reclustering reduced number of genes from " .
                     scalar(@$genes) . " to " . scalar(@clusters). "\n");
  } elsif (scalar(@clusters) > scalar(@$genes)) {
    $self->add_Error("Reclustering increased number of genes from " .
                     scalar(@$genes) . " to " . scalar(@clusters). "\n");
  }
  return \@clusters;
}

sub by_transcript_high {
  my $alow;
  my $blow;
  my $ahigh;
  my $bhigh;

  $alow = $a->start;
  $ahigh = $a->end;

  $blow = $b->start;
  $bhigh = $b->end;

  if ($ahigh != $bhigh) {
    return $ahigh <=> $bhigh;
  } else {
    return $alow <=> $blow;
  }
}

sub by_transcript_low {
  my $alow;
  my $blow;
  my $ahigh;
  my $bhigh;

  $alow = $a->start;
  $ahigh = $a->end;

  $blow = $b->start;
  $bhigh = $b->end;

  if ($alow != $blow) {
    return $alow <=> $blow;
  } else {
    return $bhigh <=> $ahigh;
  }
}


sub find_Interlocks {
  my $self = shift;
  my $clusters = shift;

  my @forward_clusters;
  my @reverse_clusters;

  foreach my $cluster (@$clusters) {
    if ($cluster->strand == 1) {
      push (@forward_clusters,$cluster);
    } else {
      push (@reverse_clusters,$cluster);
    }
  }

# Now check whether any bounds overlap

  for (my $i = 0; $i < scalar(@forward_clusters)-1; $i++) {
    my $cluster1 = $forward_clusters[$i];
    for (my $j = $i+1; $j < scalar(@forward_clusters); $j++) {
      my $cluster2 = $forward_clusters[$j];
      $self->check_for_interlock($cluster1,$cluster2,"forward");
    }
  }
  for (my $i = 0; $i < scalar(@reverse_clusters)-1; $i++) {
    my $cluster1 = $reverse_clusters[$i];
    for (my $j = $i+1; $j < scalar(@reverse_clusters); $j++) {
      my $cluster2 = $reverse_clusters[$j];
      $self->check_for_interlock($cluster1,$cluster2,"reverse");
    }
  }
}

sub check_for_interlock {
  my $self = shift;
  my $cluster1 = shift;
  my $cluster2 = shift;
  my $strand_name = shift;

#Check for overlap
  if ($cluster1->overlaps($cluster2)) {

    if (($cluster1->start < $cluster2->start && $cluster1->end < $cluster2->end) ||
        ($cluster1->start > $cluster2->start && $cluster1->end > $cluster2->end)) {
      $self->add_Error("Interlocking genes on the $strand_name strand (bounds ".
                       $cluster1->start . "-" . $cluster1->end . " and ".
                       $cluster2->start . "-" . $cluster2->end . ")\n");
    } else {
# One of the genes is fully contained within the other, but they can still
# be interlocked if not all of the exons of one gene are contained within
# a single intron in the containing gene
      if ($cluster1->start > $cluster2->start) {
        my $tmp = $cluster1;
        $cluster1 = $cluster2;
        $cluster2 = $tmp;
      }

      # foreach my $trans ($cluster1->get_Transcripts) {
      #   print "Cluster1 transcript = " . $trans->dbID . "\n";
      # }
      # foreach my $trans ($cluster2->get_Transcripts) {
      #   print "Cluster2 transcript = " . $trans->dbID . "\n";
      # }
      my @containing_exons_unsort = $self->get_cluster_Exons($cluster1);
      my @contained_exons_unsort = $self->get_cluster_Exons($cluster2);
      my @containing_exons = sort {$a->start <=> $b->start}
                                @containing_exons_unsort;
      my @contained_exons = sort {$a->start <=> $b->start}
                                @contained_exons_unsort;

      my ($left_exon_first, $right_exon_first) =
        $self->find_enclosing_Exons(\@containing_exons, $contained_exons[0]);
      my ($left_exon_last, $right_exon_last) =
        $self->find_enclosing_Exons(\@containing_exons,
                                    $contained_exons[$#contained_exons]);
      if ($left_exon_first != $left_exon_last ||
          $right_exon_first != $right_exon_last) {
        $self->add_Error("Interlocking enclosed gene on the $strand_name strand (bounds ".
                         $cluster1->start . "-" . $cluster1->end . " and ".
                         $cluster2->start . "-" . $cluster2->end . ")\n");
      } else {
        my $warnstr= "Enclosed gene on the $strand_name strand (bounds ".
                         $cluster1->start . "-" . $cluster1->end . " and ".
                         $cluster2->start . "-" . $cluster2->end . ")\n";
        $warnstr .= "   first cluster contains :";
        foreach my $trans (@{$cluster1->get_Transcripts}) {
          $warnstr .= " " . get_transcript_id($trans);
        }
        $warnstr .= "\n";
        $warnstr .= "   second cluster contains:";
        foreach my $trans (@{$cluster2->get_Transcripts}) {
          $warnstr .= " " . get_transcript_id($trans);
        }
        $warnstr .= "\n";
        $self->add_Warning($warnstr);
      }
    }
  }
}

sub get_transcript_id {
  my ($trans) = @_;

  if (defined($trans->stable_id)) {
    return $trans->stable_id;
  } elsif (defined($trans->dbID)) {
    return $trans->dbID;
  } else {
    return "no id";
  }
}

sub get_gene_id {
  my ($gene) = @_;

  if (defined($gene->stable_id)) {
    return $gene->stable_id;
  } elsif (defined($gene->dbID)) {
    return $gene->dbID;
  } else {
    return "no id";
  }
}

sub find_enclosing_Exons {
  my $self = shift;
  my $containing_exons = shift;
  my $contained_exon = shift;
  my $left_exon;
  my $right_exon;
  my $left_exon_rank;
  my $right_exon_rank;

  # print "Looking for enclosing exons for exon (range " . $contained_exon->start
  #       . " to " . $contained_exon->end .")\n";
  my $rank = 0;
  CEXON: foreach my $exon (@$containing_exons) {
    # print "Comparing to exon (range " . $exon->start." to ". $exon->end .")\n";
    if ($exon->end < $contained_exon->start) {
      $left_exon = $exon;
      $left_exon_rank = $rank;
      # print "Found left\n";
    } elsif ($exon->start > $contained_exon->end) {
      $right_exon = $exon;
      $right_exon_rank = $rank;
      # print "Found right\n";
      last CEXON;
    }
    $rank++;
  }

  if (!defined($left_exon) || !defined($right_exon)) {
    $self->throw("Didn't find enclosing exons - should not be possible\n");
  }
  if ($right_exon_rank != $left_exon_rank+1) {
    $self->throw("Left and right ranks not consecutive - should not be possible\n");
  }
  return ($left_exon, $right_exon);
}

sub get_cluster_Exons {
  my $self = shift;
  my $cluster = shift;
  my %h;

  foreach my $trans (@{$cluster->get_Transcripts}) {
    foreach my $exon ( @{$trans->get_all_Exons} ) {
      $h{"$exon"} = $exon;
    }
  }

  return values %h;
}


sub find_StrandOverlaps {
  my ($self, $clusters) = @_;

#Get all the exons in all the genes
#Split into two arrays on strand
#Sort low to high in both
#Look for overlaps

  my @forward_exons;
  my @reverse_exons;
  my $genes = $self->{'_genes'};
  my %exon_gene_map;
  foreach my $gene (@$genes) {
    my @gene_exons = @{$gene->get_all_Exons};
    if (scalar(@gene_exons)) {
      my $strand = $gene_exons[0]->strand;
      foreach my $exon (@gene_exons) {
        if ($strand != $exon->strand) {
          $self->add_Error("Gene ". $gene->dbID .
                           " has exons on different strands\n");
        }
        if ($exon->strand == 1) {
          push (@forward_exons,$exon);
        } elsif ($exon->strand == -1) {
          push (@reverse_exons,$exon);
        } else {
          $self->add_Error("Gene ". $gene->dbID .
                           " has an unstranded exon (" . $exon->dbID . ")\n");
        }
        if (exists($exon_gene_map{$exon->dbID})) {
          $self->add_Error("Seen exon " . $exon->dbID .
                           " more than once (in gene " . $gene->dbID . " and gene " .
                           $exon_gene_map{$exon->dbID}->dbID . ")\n");
        }
        $exon_gene_map{$exon->dbID} = $gene;
      }
    } else {
      $self->add_Error("Gene ". $gene->dbID . " has no exons\n");
    }
  }

  my @sorted_forward_exons = map { $_->[1] } sort { $a->[0] <=> $b->[0] } map { [$_->start, $_] } @forward_exons;
  my @sorted_reverse_exons = map { $_->[1] } sort { $a->[0] <=> $b->[0] } map { [$_->start, $_] } @reverse_exons;


#This can be optimised further
  my $has_overlaps = 0;

  my $n_r_exon = scalar(@sorted_reverse_exons);

  FEXON: foreach my $f_exon (@sorted_forward_exons) {

    for (my $i=0; $i<$n_r_exon; $i++) {
      my $r_exon = $sorted_reverse_exons[$i];
      next if (!defined($r_exon));
      if ($r_exon->end < $f_exon->start) {
        $sorted_reverse_exons[$i] = undef;
      }

      if ($r_exon->overlaps($f_exon)) {
        $self->add_Error("Overlapping exons on two strands for exons ".
                          $f_exon->dbID . " (gene = " . get_gene_id($exon_gene_map{$f_exon->dbID}) .
                         ") and " . $r_exon->dbID . " (gene = " . get_gene_id($exon_gene_map{$r_exon->dbID}) . ")\n");
        $has_overlaps = 1;
      }
      if ($r_exon->start > $f_exon->end) {
        next FEXON;
      }
    }
  }

#Only do the detailed check if there are some overlapping exons
  if ($has_overlaps) {
    #This contig has overlaps - lets see if there are any really bad ones
    print "  doing detail check\n";
    $self->find_StrandOverlap_Transcripts($clusters);
  }
}

sub find_StrandOverlap_Transcripts {
  my ($self, $clusters) = @_;

  my @forward_clusters;
  my @reverse_clusters;

  foreach my $cluster (@$clusters) {
    if ($cluster->strand == 1) {
      push (@forward_clusters,$cluster);
    } else {
      push (@reverse_clusters,$cluster);
    }
  }

  foreach my $f_cluster (@forward_clusters) {
    foreach my $r_cluster (@reverse_clusters) {
      if ($f_cluster->overlaps($r_cluster,'ignore')) {
        # print "Checking clusters\n";
        foreach my $f_trans (@{$f_cluster->get_Transcripts}) {
          foreach my $r_trans (@{$r_cluster->get_Transcripts}) {

            my @f_exons = @{$f_trans->get_all_Exons};
            my @r_exons = reverse @{$r_trans->get_all_Exons};

            if ($self->_all_exons_overlap(\@f_exons,\@r_exons,"forward")) {
              $self->add_Error("Transcript with all (". scalar(@f_exons) .
                               ") overlapping exons on other strand for " .
                               get_transcript_id($f_trans) . " and " . get_transcript_id($r_trans) .
                               " (forward)\n");
            }
            if ($self->_all_exons_overlap(\@r_exons,\@f_exons,"reverse")) {
              $self->add_Error("Transcript with all (". scalar(@r_exons) .
                               ") overlapping exons on other strand for " .
                               get_transcript_id($f_trans) . " and " . get_transcript_id($r_trans) .
                               " (reverse)\n");
            }
          }
        }
      }
    }
  }
}

sub _all_exons_overlap {
  my ($self,$exons1, $exons2, $direction) = @_;

  my $nexon2 = scalar(@$exons2);
  my $exon2_ind = 0;

  #print "Check for overlaps $direction begins\n";
  EXON1: foreach my $exon1 (@$exons1) {
    my $found = 0;
    for (; $exon2_ind < $nexon2; $exon2_ind++) {
      if ($exon1->overlaps($exons2->[$exon2_ind],'ignore')) {
        #print "NOTE NOTE Found overlap for exon " . $exon1->dbID . "\n";
        $found = 1;
        next EXON1;
      }
    }
    if (!$found) {
      #print "Failed finding exon " . $exon1->dbID . "\n";
      return 0;
    }
  }
  return 1;
}
