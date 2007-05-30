# Cared for by Ensembl
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
=head1 NAME 

Bio::EnsEMBL::Analysis::Runnable::BestTargetted

=head1 SYNOPSIS

my $runnable = Bio::EnsEMBL::Analysis::Runnable::BestTargetted->new(
      -query => $slice,
     );
  $runnable->run;
  my @genes = @{$runnable->output};

=head1 DESCRIPTION

BestTargetted combines gene-structures from any number of different targetted runs 
generating a new set containing the best genes from the sets. 

NOTES
- Each gene has only a single transcript
- First, Genes are clustered into overlapping clusters. This does not mean that 
  every Transcript in the cluster will overlap with every other Transcript in 
  the cluster. It is more like TranscriptCoalescer where each Transcript in the
  cluster overlaps with at least one other Transcript.
- Unclustered Transcripts (those that do not overlap with any other Transcript)
  are automatically kept
- If Genes in a cluster consist of those from only one analysis (logic_name), 
  then all the Genes in this cluster will be stored.
- The tricky bit. Now we have a cluster that contains Genes from more than one
  analysis. The Transcripts in the this cluster are hashed by:
      Protein -> Block -> list of genes in block
  where Protein is the protein from which the Transcript is made and Block is a
  location where all the Transcripts in this location overlap with every other 
  Transcript, and all Transcripts are made from the same protein.
- Foreach block, we want to find the best protein. This is achieved by ordering 
  the Genes in the block by logic_name, with Genes from the preferred analysis
  being first in the list. The Transcript translations are compared against the
  actual source protein. The first Transcript to match exactly is stored. If no
  Transcripts match exactly, but all Transcripts have the same Tranlasion, then
  store the first Transcript (Gene). If no Transcripts match exactly and they 
  have different translation sequences, then use exonerate to find the best Gene.
  If there is a choice between 2 short genes from one analysis vs 1 long gene 
  from a second analysis, then choose the 2 short genes (provided they have a
  good enough score/percent_id according to Exonerate).
- also, whene there is only one gene in a protein cluster, we might get suspicious.
  In these cases, exonerate the gene to check it's score and percent_id to the
  original protein. (Often a bad alignment has been produced.) Flag if the 
  exonerate results indicate that we may have a dodgy alignmennt.


=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut



package Bio::EnsEMBL::Analysis::Runnable::BestTargetted;

use strict;
use warnings; 
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(clone_Transcript identical_Transcripts);

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning info);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);




sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

    
  my( $logic_names, $seqfetcher, $verbose, $genes ) = 
      rearrange([qw(
                     LOGIC_NAMES
                     SEQFETCHER
                     VERBOSE
                     GENES
                    )], @args);

  $self->seqfetcher($seqfetcher) if defined $seqfetcher;
  $self->logic_names($logic_names) if defined $logic_names;
  $self->verbose($verbose) if defined $verbose;
  $self->all_genes($genes) if defined $genes;

  return $self ; 
}


=head2 run

   Arg : none
   Function  : Runs BestTargetted
               - compares overlapping transcripts from both sets, constructed from the same evidence 
               and stores the best model
   Returnval : none 

=cut  


sub run { 
  my ($self) = @_ ; 
  
  my %outgenes;
   
  # get all genes using both logic names on slice and cluster them
  my %logic_hash;
  foreach my $logic (@{$self->logic_names}) {
    #print STDERR "got $logic\n";
    $logic_hash{$logic} = [$logic];
  }
 
 
  # get all genes
  my @allgenes = @{$self->all_genes};

 
  # do clustering
  my ($clusters, $non_clusters) = cluster_Genes(\@allgenes, \%logic_hash ) ; 
  if ($self->verbose){
    print @$non_clusters ." non_clusters and ".@$clusters ." clusters\n";
  }


  # # #
  # loop thru non-clusters (single genes)
  # # #
  foreach my $non_cluster ( @$non_clusters ) {
    my @genes = $non_cluster->get_Genes();
    foreach my $gene (@genes){
      if ($self->verbose){
        print "storing non_cluster_gene ".$gene->dbID."\n";
      }
      $outgenes{$gene} = make_outgenes($self->analysis, $gene);
    }
  }


  # # #  
  # Now look at the clusters 
  # # #
  my @twoways;  
  foreach my $cluster (@$clusters) {
    my @inc_sets = @{$cluster->get_sets_included};
    # does the cluster contain genes from more than one set?
    if (scalar @inc_sets >= 2) {
      push @twoways, $cluster;
    } elsif (scalar @inc_sets == 1) { 
      foreach my $gene ($cluster->get_Genes_by_Set($inc_sets[0])) {
        if ($self->verbose){
          print "storing single_set_cluster_gene ".$gene->dbID." ".$inc_sets[0]."\n";
        }
        $outgenes{$gene} = make_outgenes($self->analysis, $gene);
      } # end of foreach my $gene ($cluster->get_Genes_by_Set($inc_sets[0])) {
    } # end of if (scalar @inc_sets >= 2) { 
  } #end of foreach my $cluster (@$clusters) {
  
  print "Got " . scalar(@twoways) . " twoways\n";


  # # #
  # Now look at the clusters that contain genes of more than one biotype / logic_name
  # # #
  foreach my $cluster (@twoways) {
    # check that genes have only one transcript
    # ignore gene if it does not translate
    # make a hash with keys = protein
    # and values = an array of all genes whose transcripts are made from this protein
    my %protein_clusters = %{cluster_by_protein($cluster)};    
    foreach my $protein (keys %protein_clusters) {
      my $pep = $self->seqfetcher->get_Seq_by_acc($protein)->seq;
      print "Comparing transcripts made from protein $protein with seq $pep\n";

      my @doublesets = check_for_pairs($protein_clusters{$protein});
#      # do clustering. Cluster genes made form the same protein.
#      my ($prot_clusters, $prot_non_clusters) = cluster_Genes($protein_clusters{$protein}, \%logic_hash ) ;
#      if ($self->verbose){
#        print @$prot_non_clusters ." non_clusters and ".@$prot_clusters ." clusters\n";
#      }
#      # see how many genes from each logic_name
#      # we don't expect any non_clusters
#      foreach my $non_cluster ( @$prot_non_clusters ) {
#        foreach my $gene ($non_cluster->get_Genes()) {
#          warning("Got non_cluster gene ".$gene->dbID);
#        }
#      } 
#      foreach my $cluster (@$prot_clusters) {
#        my @doublesets;
#        foreach my $set (@{$cluster->get_sets_included}) {
#          print "  Got set $set with ".$cluster->get_Genes_by_Set($set)." genes\n";
#          if ($cluster->get_Genes_by_Set($set) > 1) {
#            push @doublesets, $set;
#          } # if ($cluster->get_Genes_by_Set($set) > 1) {
#        } # foreach my $set (@{$cluster->get_sets_included}) {
#        my @unordered = $cluster->get_Genes();

      my @unordered = @{$protein_clusters{$protein}};
      my $ordered_genes = $self->get_ordered_genes(\@unordered);
      my @best_genes = @{$self->get_best_gene($protein, $ordered_genes, \@doublesets)};
      foreach my $best_gene (@best_genes) {
        if (! exists $outgenes{$best_gene}){
          $outgenes{$best_gene} = make_outgenes($self->analysis, $best_gene);
        }
      } # foreach my $best_gene (@best_genes) {
    } #foreach my $prot (keys %protein_clusters) {
  } # oreach my $cluster (@twoways) {
  my @outgenes;
  foreach my $gene ( keys %outgenes ) {
    push @outgenes, $outgenes{$gene};
  }
  $self->output(\@outgenes);

  return ;
}


sub get_best_gene {
  my ($self, $protein, $genes, $doubles) = @_;

  # genes ordered by biotype, from most preferred to least preferred biotype
  if (scalar(@$doubles) > 0) {
    print "*** Have doubles so need to be careful ***\n";
  }
  my $pep = $self->seqfetcher->get_Seq_by_acc($protein)->seq;
  my @best;
  my $exonerate;
  print "Finding best gene for protein $protein...\n";


  # Go thru each gene
  # See whether its tln matches the peptide
  my @match_pep;
  my %tln_groups;
  my $tln_seq;
  foreach my $gene (@$genes) {
    foreach my $transcript (@{$gene->get_all_Transcripts}) {
      $tln_seq = $transcript->translate->seq;
      if ($tln_seq eq $pep) {
        push @match_pep, $gene;
      }
      push @{$tln_groups{$tln_seq}}, $gene; 
    }
  }
  my $num_tln_groups = scalar(keys %tln_groups);
  print "Have $num_tln_groups different possible translations represented.\n"; 
  foreach my $s (keys %tln_groups) {
    print "option: $s\n";
  }


# Now go thru the results and see whether we have any 
  # exact matches or whether we need to exonerate
  if (scalar(@match_pep == 1)) {
    # we don't care if we have doubles or not
    print "Only one exact match to protein. \nKeeping only gene ".$match_pep[0]->dbID." (".$match_pep[0]->biotype.")\n";
    @best = ( $match_pep[0] );

  } elsif (scalar(@match_pep > 1)) {
    print "More than one exact match to protein.\n"; 
    if (scalar(@$doubles) > 0) {
      print "Have doubles and more than one match to the protein.\n";
      $exonerate = 1;
    } else {
      print "Keeping first gene ".$match_pep[0]->dbID.
            " from array as we prefer it's analysis (".$match_pep[0]->biotype.")\n";
      @best = ( $match_pep[0] );
    }

  } else {
    print "No exact match to protein. ";
    # we may have the translated sequences being all the same, 
    # or we may have the translated sequences being different.
    # take the first gene if they're all the same
    if ((scalar(@$genes) > 1) && (scalar(keys %tln_groups) == 1) && (scalar(@$doubles < 1))) {
      print "All transcripts have the same translated sequence. \nKeeping first gene ".$tln_groups{$tln_seq}->[0]->dbID.
            " from array as we prefer it's analysis (".$tln_groups{$tln_seq}->[0]->biotype.")\n";
      print "*** Maybe we should exonerate or look at exon structure to decide what to choose?\n";
      @best = ( $tln_groups{$tln_seq}->[0] );
    } else {
      print "We have ".scalar(@$genes)." genes, ".scalar(keys %tln_groups)." different sequences and ".scalar(@$doubles)." doubles. Let's dig deeper.\n";
      $exonerate = 1;
    }
  }
 

  # More clarification is needed.
  # So we do exonerates
  # we keep track of all genes' scores and percent_ids for later use
  if (defined $exonerate) {
    print "Exonerate...\n";
    my $max_gene_score = 0;
    my $max_perc_id = 0;
    my $best;
    my %scores;
    my $target = Bio::Seq->new( -display_id => $protein, -seq => $pep);
    foreach my $gene (@$genes) {
      if (!defined $best) {
        # sometimes, when there is only one gene and it is very short, it has a score = 0
        print "First gene is ".$gene->dbID."\n";
        $best = $gene;
      }
      foreach my $transcript (@{$gene->get_all_Transcripts}) {
        my $query = Bio::Seq->new( -display_id => $transcript->dbID, -seq => $transcript->translate->seq);
        my ($gene_score, $perc_id) = @{$self->run_exonerate($query, $target)}; 
       # my $gene_score = shift @{$self->run_exonerate($query, $target)};

        print "Gene ".$gene->dbID." (".$gene->biotype.") has score $gene_score and percent_id $perc_id\n";
        $scores{$gene->biotype}{$gene->dbID}{'score'} = $gene_score;
        $scores{$gene->biotype}{$gene->dbID}{'percent_id'} = $perc_id; 
        if ($gene_score > $max_gene_score) {
          $max_gene_score = $gene_score;
          $max_perc_id = $perc_id;
          $best = $gene;
          print "Gene with max score $max_gene_score is now ".$gene->dbID." (with percent_id $perc_id)\n";
        } # $gene_score > $max_gene_score
      } # reach my $transcript (@{$gene->get_all_Transcripts}) {
    } # foreach my $gene (@$genes) {

    # ok, we have all info that we need. Now 
    # see what we're dealing with. Are there doubles 
    # or can we just take the gene with the highest score?
    my %allowed_doubles;
    if (scalar(@$doubles) == 0) {
      print "No doubles.\nKeeping gene ".$best->dbID." (".$best->biotype.
            ") with highest score of $max_gene_score and percent_id $max_perc_id\n";
      @best = ($best);
    } elsif (scalar(@$doubles) >= 1) {
      print "We have ".scalar(@$doubles)." sets of doubles\n";
      # check that all genes have score >= 90% of the max
      # and then choose the biotype with the highest collective score
      foreach my $dbl_logic (@$doubles) {
        my $count = 0;
        $allowed_doubles{$dbl_logic} = 0;
        foreach my $gid (keys %{$scores{$dbl_logic}}) {
          if ($scores{$dbl_logic}{$gid}{'score'} >= 0.9*$max_gene_score) {
            # gene has sufficiently high score that we want to leep it
            $count++;
            $allowed_doubles{$dbl_logic} += $scores{$dbl_logic}{$gid}{'score'};
          } else {
            print "  Not acceptable: gene $gid has score ".$scores{$dbl_logic}{$gid}{'score'}." which is too low\n";
          } #if
        } # gid
        if ($count != scalar(keys %{$scores{$dbl_logic}})) {
          # we don't want to use this biotype as one or more of its gene has a low score
          # NOTE: what happens if there are 3 genes for the one biotype, with two non-overlapping
          # genes and the third gene having a long intron that joins the other two?
          # Do we want to consider all 3 genes or only the non-overlapping ones?
          print "  Not allowed to use biotype $dbl_logic. (Have ".scalar(keys %{$scores{$dbl_logic}}).
                " genes but only $count pass the threshold 90%)\n";
          $allowed_doubles{$dbl_logic} = 0; 
        }
      } # foreach my $dbl_logic (@$doubles) { 
    
      # choose the one with the max score
      # note that if 2 biotypes have the same max score, then there is no
      # distinguishing which one is betetr
      my $best_biotype;
      my $max_dbl_score = 0;
      
      #sort allowed doubles by order in config
      my @sorted_allowed;
      foreach my $logic (@{$self->logic_names}) {
	if (exists $allowed_doubles{$logic}){
	  push @sorted_allowed, $logic;
	}        
      }
      foreach my $dbl_logic (@sorted_allowed) {
	
	if (defined $dbl_logic && ($allowed_doubles{$dbl_logic} > $max_dbl_score)) {
          $max_dbl_score = $allowed_doubles{$dbl_logic};
          $best_biotype = $dbl_logic;
        } 
      }
      if ($max_dbl_score > 0) {
        # we take the doubles
        foreach my $gene (@$genes) {
          if ($gene->biotype eq $best_biotype) {
            print "Keeping gene ".$gene->dbID." biotype ".$gene->biotype."\n";
            push @best, $gene;
          }
        }
      } else {
        @best = ($best);
        print "Keeping gene ".$best->dbID." biotype ".$best->biotype."\n";
      }
    } # if (scalar(@$doubles) > 0) {
    foreach my $b (@best) {
      if ($scores{$b->biotype}{$b->dbID}{'percent_id'} < 90) {
        print "FLAG: best_gene ".$b->dbID." has percent_id ".$scores{$b->biotype}{$b->dbID}{'percent_id'}."<90 (start ".$b->start." end ".$b->end.")\n";
      }
    }
  } #if (defined $exonerate) {
  return \@best;
}

 
sub check_for_pairs {
  my ($genes) = @_;
  my %logics;

  print "we have ". scalar @$genes ." genes to compare\n";

  
  # all the genes passed in are made from the same protein.
  # see if there are any NON-OVERLAPPING genes of the same biotype  
  foreach my $outer (@$genes) {
    foreach my $inner (@$genes) {
      if (($outer->biotype eq $inner->biotype) && ($outer->start > $inner->end || $outer->end < $inner->start)) {
        #they have the same biotype
        #they do not overlap
	$logics{$outer->biotype} = 1;
      }
    }
  }
  return keys %logics;
}


sub get_ordered_genes {
  my ($self, $genes) = @_;

  my @ordered;
  my %seen;  
  
  print "block has ".scalar(@$genes)." genes\n";
  foreach my $logic (@{$self->logic_names}) {
    foreach my $gene (@$genes) {  
      if (!exists $seen{$gene} && ($logic eq $gene->biotype)) {
        print "  Got in cluster: gene ".$gene->dbID." of biotype ".$gene->biotype."\n";
        push @ordered, $gene;
        $seen{$gene} = 1;
      }
    }
  }
  return \@ordered;
}

sub cluster_by_protein {
  my ($cluster) = @_;
  my %protein_hash;

  # get all genes made from 1 protein
  
  GENE: foreach my $gene ($cluster->get_Genes()) {
    my ($gene_ok, $protein) = check_gene($gene);
    
    if (!$gene_ok) {
      throw("Gene check failed for gene ".$gene->dbID."");
    }
    
    my $duplicate; 
    if (exists $protein_hash{$protein}) {
      GOT: foreach my $got_gene (@{$protein_hash{$protein}}){
        if ($gene->biotype ne $got_gene->biotype) {
          next GOT;
        }
	# we assume that each gene has only one Transcript
	my @transcripts1 = @{$gene->get_all_Transcripts()};
	my @transcripts2 = @{$got_gene->get_all_Transcripts()}; 

	if (identical_Transcripts($transcripts1[0], $transcripts2[0])) {
          $duplicate = 1; 
	}  
      }
      if (!$duplicate) {
        push @{$protein_hash{$protein}}, $gene;
      }    
    } else {
      push @{$protein_hash{$protein}}, $gene;
    }
    

    # If genes 1, 2 and 3 are all made from the same protein,
    # then they will be grouped together even if gene 2
    # and gene3 do not overlap.
    #                      gene 1
    #         <--------------------------------------->
    #   gene 2                                     gene 3
    # <---------->                              <-------------->

  }
  return \%protein_hash;
}
=head2 get_alternate_transcripts

  Arg [0]   : runnable
  Arg [1]   : ref to a hash
  Function  :
  Returntype: 2 hashrefs

=cut
sub get_alternate_transcripts {
  my ($cluster) = @_;

  my %alternates;
  my %seen;

  #print "\n\nStarting new cluster with ".scalar($cluster->get_Genes())." genes:\n";
  
  GENE: foreach my $gene ($cluster->get_Genes()) {
    #print "Doing GENE ".$gene->dbID."...";
    if (exists $seen{$gene->dbID}) {
      #print "Already seen. Skipping...\n";
      next GENE;
    }
    # CHECK THAT:
    # each gene has only 1 transcript
    # each transcript translates
    # each transcript has one tsf
    my ($gene_ok, $protein) = check_gene($gene);
    if (!$gene_ok) {
      throw("Gene check failed for gene ".$gene->dbID."");
    }
    
    # add to the hash
    push @{$alternates{$protein}{$gene->dbID}}, $gene;
    $seen{$gene->dbID} = 1;
    my $exon_string = make_exon_string(@{$gene->get_all_Transcripts}[0]);

    INNER: foreach my $inner_gene ($cluster->get_Genes()) {
      #print "Doing INNER gene ".$inner_gene->dbID."...";
      if (exists $seen{$inner_gene->dbID}) {
        #print "Already seen. Skipping...\n";
        next INNER;
      }
      # DO CHECKS as above
      my ($inner_gene_ok, $inner_protein) = check_gene($inner_gene);
      if (!$inner_gene_ok) {
        throw("Gene check failed for gene ".$inner_gene->dbID."");
      }
      if ($protein ne $inner_protein) {
        #print "Proteins from GENE $protein and INNER $inner_protein do not match. Skipping...\n";
        next INNER;
      }
      my $inner_exon_string = make_exon_string(@{$inner_gene->get_all_Transcripts}[0]);

      # Now look for unique genes made form the same protein
      if ($exon_string eq $inner_exon_string) {
        # they are identical
        #print "Genes ".$gene->dbID." and ".$inner_gene->dbID." have same exon string\n";
        $seen{$inner_gene->dbID} = 1; 
        next INNER;
      } elsif ($gene->start <= $inner_gene->end && $gene->end >= $inner_gene->start) {
        # they overlap and should be investigated
        # but does it overlap with everything in the array?
        #print "outer GENE ".$gene->dbID." overlaps with INNER gene ".$inner_gene->dbID."\n";
        if (overlaps_all($alternates{$protein}{$gene->dbID}, $inner_gene) ) {
          #print "  INNER gene ".$inner_gene->dbID." overlaps all in cluster so far\n";
          # not sure how this affects a "joiner" - a transcript that joins 2 clusters
          push @{$alternates{$protein}{$gene->dbID}}, $inner_gene;
          $seen{$inner_gene->dbID} = 1;
        }
        # what happens when we are dealing with a case like that of gene 1?
        # if gene 1 can only be used once, the we might miss out on something important
        # how often will this occur though, in reality?
        #                      gene 1
        #         <--------------------------------------->
        #   gene 2                                     gene 3
        # <---------->                              <-------------->

      } else {
        # they do not overlap so ignore
        #print "  INNER gene ".$inner_gene->dbID." does not overlap all in cluster so far\n";
        print "Genes ".$gene->dbID." (exon string ".$exon_string.") and ".
              $inner_gene->dbID." (exon string $inner_exon_string) do not overlap.\n"; 
      } #if ($exon_string eq $inner_exon_string) {
    } #inner
  } # genes
  return \%alternates;
}

sub overlaps_all {
  my ($gene_array, $gene)  = @_;
  my $num_overlapped;
  my $overlaps_all;

  foreach my $g (@$gene_array) {
    if ($g->start <= $gene->end && $g->end >= $gene->start) {
      $num_overlapped++;
    }
  }
  if ($num_overlapped == scalar(@$gene_array)) {
    $overlaps_all = 1;
    #print "overlaps all $num_overlapped genes\n";
  } else {
    print "does not overlap all genes\n";
  }
  return $overlaps_all;
}

sub make_exon_string {
  my ($transcript) = @_;
  my $exon_string = "E";
  
  my @exons = sort {$a->start <=> $b->start} @{$transcript->get_all_translateable_Exons};
  foreach my $exon (@exons) {
    $exon_string .= ":".$exon->start.":".$exon->end;
  }
  return $exon_string;
}


sub check_gene {
  my ($gene) = @_;

  my $pass_checks = 1;
  my $protein;
  my @transcripts = @{$gene->get_all_Transcripts};

  if (scalar(@transcripts) != 1) {
    print STDERR "Gene ".$gene->stable_id." with dbID ".$gene->dbID." has ".scalar(@transcripts)." transcripts\n";
    $pass_checks = 0;
  }

  foreach my $transcript (@transcripts) {
    my @tsfs = @{$transcript->get_all_supporting_features};
    if (scalar(@tsfs) == 1) {
      $protein = $tsfs[0]->hseqname;
    } else {
      print STDERR "WARNING: transcript ".$transcript->stable_id." (dbID ".$transcript->dbID.") has ".scalar(@tsfs)." tsfs. ".
                   "Should only have one transcript_supporting_feature\n";
      $pass_checks = 0;
    }
    # check for a translation
    if (!defined($transcript->translation)) {
      print STDERR "Transcript ".$transcript->stable_id." does not translate.\n";
      $pass_checks = 0;
    }
  }
  return $pass_checks, $protein;
}
 

=head2 make_outgenes

  Arg [0]   : ref to a Bio::EnsEMBL::Analysis 
  Arg [1]   : ref to a Bio::EnsEMBL::Gene 
  Arg [2]   : ref to a Bio::EnsEMBL::Transcript (optional)
  Function  : Clones the gene object and assigns the new analysis to it
  Returntype: new Bio::EnsEMBL::Gene 
  Example   : 

=cut




sub make_outgenes {
  my ($analysis, $gene, $trans) = @_;
  
  my $outgene = copy_empty_gene($analysis, $gene);

  if (not defined $trans){
    my @trans = @{$gene->get_all_Transcripts};
  
    if (scalar @trans > 1){
      print STDERR "gene ".$$gene->dbID." has >1 transcript\n";
    }else{
      $trans = $trans[0];
    }
  }

  my $outtrans = clone_Transcript($trans); 
  $outtrans->analysis($analysis);
  $outtrans->dbID(undef);
  
 
  #update the transcript_supporting_feature analysis_ids to match the transcript
  foreach my $sf(@{$outtrans->get_all_supporting_features}){
    $sf->dbID(undef);
    $sf->analysis($analysis);
  }
  
  #update the supporting_feature analysis_ids to match the transcript
  foreach my $exon (@{$outtrans->get_all_Exons}){
    $exon->dbID(undef);
    foreach my $sf (@{$exon->get_all_supporting_features}){
      $sf->dbID(undef);
      $sf->analysis($analysis);
    }
  }
  
  $outgene->add_Transcript($outtrans);
 
  return $outgene;
}

=head2 copy_empty_gene

  Arg [0]   : ref to a Bio::EnsEMBL::Analysis
  Arg [1]   : ref to a Bio::EnsEMBL::Gene 
  Function  : makes a gene and sets biotype and analysis
  Returntype: new Bio::EnsEMBL::Gene

=cut



sub copy_empty_gene {
  my ($analysis, $gene) = @_;
  
  my $newgene = new Bio::EnsEMBL::Gene;
  $newgene->biotype( $gene->biotype);
  $newgene->analysis($analysis);
  
  return $newgene;
}

=head2 run_exonerate

  Arg [0]   : query  Bio::Seq object
  Arg [1]   : target Bio::Seq object
  Function  : calls exonerate->run
  Returntype: score - integer

=cut
      


sub run_exonerate {
  my ($self, $query, $target) = @_;
  
  my $exonerate = Bio::EnsEMBL::Analysis::Runnable::BestTargettedExonerate->new
      (
       -program     => $self->program,
       -analysis    => $self->analysis,
       -query_seqs  => [$query],
       -target_seqs => [$target],
       -options     => "--model affine:local",
       );
       
  
  $exonerate->run;
  my ($score, $perc_id) = @{$exonerate->output} ;
  
  return [$score, $perc_id];
}

########################
sub seqfetcher {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_bt_seqfetcher} = $val;
  }

  return $self->{_bt_seqfetcher};
}

sub logic_names {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_bt_logic_names} = $val;
  }

  return $self->{_bt_logic_names};
}

sub verbose {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_bt_verbose} = $val;
  }

  return $self->{_bt_verbose};
}

sub all_genes {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_all_genes} = $val;
  }

  return $self->{_all_genes};
}




##########################################
### local class
##########################################

package Bio::EnsEMBL::Analysis::Runnable::BestTargettedExonerate;

use Bio::EnsEMBL::Analysis::Runnable::BaseExonerate;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::BaseExonerate);


sub target_type {
  my ($self) = @_;

  return 'protein';
}


sub query_type {
  my ($self) = @_;
  
  return 'protein';
}



sub parse_results {
  my ($self, $fh) = @_;

  my $score = 0;
  my $perc_id = 0;
  
  # my $basic_options = "--showsugar false --showvulgar false --showalignment false --ryo \"RESULT: %S %pi %ql %tl %g %V\\n\" ";
  while (<$fh>){
    if ($_=~/^RESULT:/){
      my @tmp = split/\s+/, $_;
#     my (
#        $tag, $q_id, $q_start, $q_end, $q_strand,
#        $t_id, $t_start, $t_end, $t_strand, $score,
#        $perc_id, $q_length, $t_length, $gene_orientation,
#        @vulgar_blocks
#        ) = split;
#     my($match_type, $query_match_length, $target_match_length,@rest_of_vulgar) = @vulgar_blocks;
      if ($score < $tmp[9]){
        $score = $tmp[9]; #return the best score
        $perc_id = $tmp[10];
      }
    }
  }  

  $self->output([$score, $perc_id]);
}


1; 
