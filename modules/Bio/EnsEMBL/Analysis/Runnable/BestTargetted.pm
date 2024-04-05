=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Runnable::BestTargetted - 

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
- If Genes in a cluster consist of those from only one analysis (biotype), 
  then all the Genes in this cluster will be stored.
- The tricky bit. Now we have a cluster that contains Genes from more than one
  analysis. The Transcripts in the this cluster are hashed by:
      Protein -> Block -> list of genes in block
  where Protein is the protein from which the Transcript is made and Block is a
  location where all the Transcripts in this location overlap with every other 
  Transcript, and all Transcripts are made from the same protein.
- Foreach block, we want to find the best protein. This is achieved by ordering 
  the Genes in the block by biotype, with Genes from the preferred analysis
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


=cut



package Bio::EnsEMBL::Analysis::Runnable::BestTargetted;

use strict;
use warnings; 
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(clone_Transcript identical_Transcripts);
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(align_proteins);

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning info );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use parent ('Bio::EnsEMBL::Analysis::Runnable');




sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

    
  my( $biotypes, $seqfetcher, $verbose, $genes, $keep_single_analysis, $cluster_on_coding_exons, $min_identity, $min_coverage ) =
      rearrange([qw(
                     BIOTYPES 
                     SEQFETCHER
                     VERBOSE
                     GENES
                     KEEP_SINGLE_ANALYSIS
                     CLUSTER_ON_CODING_EXONS
                     MIN_IDENTITY
                     MIN_COVERAGE
                    )], @args);

  $self->seqfetcher($seqfetcher) if defined $seqfetcher;
  $self->biotypes($biotypes) if defined $biotypes;
  $self->verbose($verbose) if defined $verbose;
  $self->all_genes($genes) if defined $genes;
  $self->keep_single_analysis($keep_single_analysis) if defined $keep_single_analysis;
  $self->cluster_on_coding_exons($cluster_on_coding_exons) if defined $cluster_on_coding_exons;
  $self->min_coverage($min_coverage || 0);
  $self->min_identity($min_identity || 0);

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
  foreach my $logic (@{$self->biotypes}) {
    print STDERR "got $logic\n";
    $logic_hash{$logic} = [$logic];
  }
 
 
  # get all genes
  my @allgenes = @{$self->all_genes};

 
  # do clustering
  if ($self->cluster_on_coding_exons) {
    print "Clustering on coding exons only\n";
  } else {
    print "Clustering on coding AND non-coding exons\n";
  }
  my ($clusters, $non_clusters) = cluster_Genes(\@allgenes, \%logic_hash , $self->cluster_on_coding_exons) ; 
  if ($self->verbose){
    print @$non_clusters ." non_clusters and ".@$clusters ." clusters\n";
  }


  # # #
  # loop thru non-clusters (single genes)
  # # #
  foreach my $non_cluster ( @$non_clusters ) {
    my @genes = @{ $non_cluster->get_Genes() };
    GENE: foreach my $gene (@genes){
      # need to check that they have a tln
      my ($gene_ok, $protein) = check_gene($self, $gene);
      if (!$gene_ok) {
        warning("Gene check failed for unclustered gene ".$gene->dbID." - removed from analysis");
        next GENE;
      }

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
      GENE: foreach my $gene (@{ $cluster->get_Genes_by_Set($inc_sets[0]) }) {
        # need to check that they have a tln
        my ($gene_ok, $protein) = check_gene($self, $gene);
        if (!$gene_ok) {
          warning("Gene check failed for gene ".$gene->dbID." (in cluster of only one Set) - removed from analysis");
          next GENE;
        }
        if ($self->verbose){ 
           my $gene_id = $gene->dbID ? $gene->dbID : "" ; 
           print "storing single_set_cluster_gene ".$gene_id ." ".$inc_sets[0]."\n";
         }
        $outgenes{$gene} = make_outgenes($self->analysis, $gene);
      } # end of foreach my $gene ($cluster->get_Genes_by_Set($inc_sets[0]))
    } # end of if (scalar @inc_sets >= 2) 
  } #end of foreach my $cluster (@$clusters) 
  
  print "Got " . scalar(@twoways) . " twoways\n";


  # # #
  # Now look at the clusters that contain genes of more than one biotype / biotype
  # # #
  foreach my $cluster (@twoways) {
    # check that genes have only one transcript
    # ignore gene if it does not translate
    # make a hash with keys = protein
    # and values = an array of all genes whose transcripts are made from this protein
    my $protein_clusters = cluster_by_protein($self, $cluster);    
    my $filtered = filter_transcripts($protein_clusters);

    foreach my $protein (keys %{$protein_clusters}) {  

      my $fetched_protein;  

      #   
      # Seqfetcher does NOT throw just returns undef so we have to use eval block here to  
      # make the job die - otherwise job dies quietly which is BAD as it will be recorded successfully 
      #    
      
      eval { 
            $fetched_protein = $self->seqfetcher->get_Seq_by_acc($protein) ;     
           }; 

      if ($@) {    
         throw ( " Could not get sequence for protein $protein \nfrom index ". $self->seqfetcher->index_name . 
                 " \nIf you use OBDA Index : Maybe your keys are redundant ? \n\n$@ " ) ; 
      }   

      unless ( $fetched_protein ) { 
        throw("SeqFetcher cannot get_Seq_by_acc for protein '$protein' in indexfile : ".  $self->seqfetcher->index_name );
      }  

      my $pep = $fetched_protein->seq; 
      print "Comparing transcripts made from protein $protein with seq $pep\n";
      if (!$pep) {
        throw("Unable to fetch peptide sequence for protein $protein");
      }
  
      # filter the genes
      my @best_genes = @{$self->get_best_gene($protein, $self->biotypes, $filtered)};
      foreach my $best_gene (@best_genes) {
        if (! exists $outgenes{$best_gene}){
          $outgenes{$best_gene} = make_outgenes($self->analysis, $best_gene);
        }
      } 
    } #foreach my $prot (keys %protein_clusters) 
  } # oreach my $cluster (@twoways) 
  my @outgenes;
  foreach my $gene ( keys %outgenes ) {
    push @outgenes, $outgenes{$gene};
  }
  $self->output(\@outgenes);

  #return ;
}


sub get_best_gene {
  my ($self, $protein, $biotypes, $hsh) = @_;

  my %hash = %{$hsh};

  # we don't wanna make a gene here if only one analysis has built here
  my $biotypes_represented;
  my $token_gene;
  foreach my $biotype (keys %{$hash{$protein}}) {
    if (scalar(@{$hash{$protein}{$biotype}}) > 0) {
      $biotypes_represented++;
      $token_gene = $hash{$protein}{$biotype}->[0];
    }
  }
  if ($biotypes_represented == 1 && scalar(@$biotypes) > 2) {
    print "FLAG_NULL1: For $protein, only one biotype represented. See gene ".$token_gene->dbID." (".$token_gene->biotype.")\n";
    if (!defined $self->keep_single_analysis) {
      print "FLAG_NULL2: For $protein, am returning no genes.\n";
      return [];
    }
  }

  # which biotypes have more than 1 gene made from $protein?
  my @doubles;
  my $num_genes = 0;
  #foreach my $type (@$biotypes) {
    foreach my $biotype (keys %{$hash{$protein}}) {
        #my @g = $hash->{$protein}->{$biotype};
        #print "have ".scalar(@g)." genes, existing...\n"; 
        print "biotype $biotype has ".scalar(@{$hash{$protein}{$biotype}})." genes \n";
        $num_genes += scalar(@{$hash{$protein}{$biotype}});
        if (scalar(@{$hash{$protein}{$biotype}}) > 1) {
          print "*** Have ".$biotype." doubles so need to be careful ***\n";
          push @doubles, $biotype;

        }
    }
  #} 

  # fetch the protein
  # for cdna2genome, we might need to fix things here. Getting a protein seq this way
  # we need to modify the input config file to have peptide seq with a cdna id header
  my $pep = $self->seqfetcher->get_Seq_by_acc($protein)->seq;
  if (!$pep) {
    throw("Unable to fetch peptide sequence");
  }
  my @best;
  my $exonerate; # flag
  print "Finding best gene for protein $protein...\n";

  if (scalar(@doubles) > 0) {
    print "We have doubles so let's just use them\n";
    OUTER: foreach my $biotype (@$biotypes) {
      INNER: foreach my $type (@doubles) {
        if ($biotype eq $type) {
          print "The most favoured set of doubles is $biotype so we takes them\n";
          foreach my $gene (@{$hash{$protein}{$biotype}}) {
            print "Keeping gene ".$gene->dbID."\n";
          }
          push @best, @{$hash{$protein}{$biotype}};
          return \@best;
        }
      }
    }
  }

  # Go thru each gene in order of favourite to least fav biotype
  # See whether its tln matches the peptide
  my @match_pep;
  my %tln_groups;
  my $tln_seq;
  OUTER: foreach my $type (@$biotypes) {
    INNER: foreach my $biotype (keys %{$hash{$protein}}) {
      if ($biotype ne $type) {
        next INNER;
      }
      if (scalar(@{$hash{$protein}{$biotype}}) > 0) {
        foreach my $gene (@{$hash{$protein}{$biotype}}) {
          print "Got gene dbid ".$gene->dbID." ".$gene->biotype."\n";
          foreach my $transcript (@{$gene->get_all_Transcripts}) {
            $tln_seq = $transcript->translate->seq;
            if ($tln_seq eq $pep) {
              push @match_pep, $gene;
            }
            push @{$tln_groups{$tln_seq}}, $gene; 
          }
        }
      } else {
        print "No genes for biotype $biotype\n";
      }
    }
  }
  my $num_tln_groups = scalar(keys %tln_groups);
  print "Have $num_tln_groups different possible translations represented.\n"; 
  foreach my $s (keys %tln_groups) {
    print "option ".$tln_groups{$s}[0]->biotype.": $s\n";
  }


# Now go thru the results and see whether we have any 
  # exact matches or whether we need to exonerate
  if (scalar(@match_pep == 1)) {
    # we don't care if we have doubles or not
    print "Only one exact match to protein. \nKeeping only gene ".$match_pep[0]->dbID." (".$match_pep[0]->biotype.")\n";
    @best = ( $match_pep[0] );

  } elsif (scalar(@match_pep > 1)) {
    print "Have ". scalar(@match_pep) ." exact matches to protein.\n"; 
    if (scalar(@doubles) > 0) {
      print "Have doubles and more than one match to the protein.\n";
      $exonerate = 1;
    } else {
      print "Keeping first gene ".$match_pep[0]->dbID.
            " from array as we prefer its analysis (".$match_pep[0]->biotype.")\n";
      @best = ( $match_pep[0] );
    }

  } else {
    print "No exact match to protein. ";
    # we may have the translated sequences being all the same, 
    # or we may have the translated sequences being different.
    # take the first gene if they're all the same
    if (($num_genes > 1) && (scalar(keys %tln_groups) == 1) && (scalar(@doubles < 1))) {
      print "All transcripts have the same translated sequence. \nKeeping first gene ".$tln_groups{$tln_seq}->[0]->dbID.
            " from array as we prefer its analysis (".$tln_groups{$tln_seq}->[0]->biotype.")\n";
      print "*** Maybe we should exonerate or look at exon structure to decide what to choose?\n";
      @best = ( $tln_groups{$tln_seq}->[0] );
    } else {
      print "We have $num_genes genes, ".scalar(keys %tln_groups)." different sequences ".
            " and ".scalar(@doubles)." doubles. Let's dig deeper.\n";
      $exonerate = 1;
    }
  }
 

  # More clarification is needed.
  # So we do exonerates
  # we keep track of all genes' scores and percent_ids for later use
  if (defined $exonerate) {
    print "Exonerate...\n";
    my $gene_tln_comp;
    my $max_gene_score = 0;
    my $max_perc_id = 0;
    my $best;
    my $first_best;
    my %scores;
    # make sure that you're using the correct protein and pep for cdna2genome
    # Exonerate doesn't like ambiguity codes and $pep is not used anymore
    $pep =~ tr/BJZ/XXX/;
    my $target = Bio::Seq->new( -display_id => $protein, -seq => $pep);
    foreach my $biotype (@$biotypes) {
      if (defined $hash{$protein}{$biotype}) {
        foreach my $gene (@{$hash{$protein}{$biotype}}) {
          if (!defined $best) {
            # sometimes, when there is only one gene and it is very short, it has a score = 0
            print "First gene is ".$gene->dbID."\n";
            $best = $gene;
            $first_best = $gene;
          }
          foreach my $transcript (@{$gene->get_all_Transcripts}) {
            print "Running transcript ".$transcript->dbID."\n";
            my $query = Bio::Seq->new( -display_id => $transcript->dbID, -seq => $transcript->translate->seq);
            my ($gene_score, $perc_id) = @{$self->run_exonerate($query, $target)}; 
           # my $gene_score = shift @{$self->run_exonerate($query, $target)};
    
            print "Gene ".$gene->dbID." (".$gene->biotype.") has score $gene_score and percent_id $perc_id (ratio ".
                  ($transcript->translation->length)/($transcript->end - $transcript->start).")\n";
            $scores{$gene->biotype}{$gene->dbID}{'score'} = $gene_score;
            $scores{$gene->biotype}{$gene->dbID}{'percent_id'} = $perc_id; 
            $scores{$gene->biotype}{$gene->dbID}{'ratio'} = ($transcript->translation->length)/($transcript->end - $transcript->start);
            $scores{$gene->biotype}{$gene->dbID}{'gene'} = $gene;
            if ($gene_score > $max_gene_score) {
              $max_gene_score = $gene_score;
              $best = $gene;
              print "Gene with max score $max_gene_score is now ".$gene->dbID." (with percent_id $perc_id)\n";
            } # $gene_score > $max_gene_score
            if ($perc_id > $max_perc_id) {
              $max_perc_id = $perc_id;
            }
          } # reach my $transcript (@{$gene->get_all_Transcripts}) 
        } # foreach my $gene (@$genes) 
      }
    } # foreach my $biotype (@$biotypes) 

    # ok, we have all info that we need. Now 
    # see what we're dealing with. Are there doubles 
    # or can we just take the gene with the highest score? 
    my %allowed_doubles;
    if (scalar(@doubles) == 0) {
      print "No doubles, Checking the ratio\n";
      if ($scores{$best->biotype}{$best->dbID}{'ratio'} >= 0.8*($scores{$first_best->biotype}{$first_best->dbID}{'ratio'})) {
        @best = ($best);
        print "ratio OK. Keeping gene ".$best->dbID." biotype ".$best->biotype.
             " , ratio ".$scores{$best->biotype}{$best->dbID}{'ratio'}." vs ".
              $scores{$first_best->biotype}{$first_best->dbID}{'ratio'}."\n";
      } else {
        print "Ratio too low. Befault to 1st biotype. Keeping gene ".$first_best->dbID." biotype ".$first_best->biotype."\n";
        @best = ($first_best);
      }

    } elsif (scalar(@doubles) >= 1) {
      print "We have ".scalar(@doubles)." sets of doubles\n";
      # check that all genes have score >= 90% of the max
      # and then choose the biotype with the highest collective score
      foreach my $dbl_biotype (@doubles) {
        $allowed_doubles{$dbl_biotype} = 0;
        foreach my $gid (keys %{$scores{$dbl_biotype}}) {
          if ($scores{$dbl_biotype}{$gid}{'score'} >= 0.9*$max_gene_score) {
            # gene has sufficiently high score that we want to leep it
            $allowed_doubles{$dbl_biotype} += $scores{$dbl_biotype}{$gid}{'score'};
          } else {
            print "  Not acceptable: gene $gid has score ".$scores{$dbl_biotype}{$gid}{'score'}." which is too low\n";
          } #if
        } # gid
      } # foreach my $dbl_biotype (@$doubles)  
    
      # choose the one with the max score
      # note that if 2 biotypes have the same max score, then there is no
      # distinguishing which one is betetr
      my $best_biotype;
      my $max_dbl_score = 0;
      
      foreach my $dbl_biotype (@$biotypes) {
	
	if (defined $dbl_biotype && ($allowed_doubles{$dbl_biotype} > $max_dbl_score)) {
          $max_dbl_score = $allowed_doubles{$dbl_biotype};
          $best_biotype = $dbl_biotype;
        } 
      }
      if ($max_dbl_score > 0) {
        # we take the doubles
        foreach my $gene (@{$hash{$protein}{$best_biotype}}) {
          if (($gene->biotype eq $best_biotype) && ($scores{$gene->biotype}{$gene->dbID}{'score'} >= 0.9*$max_gene_score)) {
            print "Keeping gene ".$gene->dbID." biotype ".$gene->biotype."\n";
            push @best, $gene;
            if ($scores{$gene->biotype}{$gene->dbID}{'percent_id'} < 95) {
              print "WARNING_FLAG: percent_id is less than max\n";
            }
          } else {
            print "(Not keeping ".$gene->dbID." as score too low)\n";
          }
        }
      }
    }
    foreach my $b (@best) {
      if ($scores{$b->biotype}{$b->dbID}{'percent_id'} < 90) {
        print "FLAG: best_gene ".$b->dbID." has percent_id ".$scores{$b->biotype}{$b->dbID}{'percent_id'}.
              "<90 (start ".$b->start." end ".$b->end.")\n";
      }
    }
  } #if (defined $exonerate)  
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
	$logics{$outer->biotype} += 1;
      }
    }
  }
  foreach my $k (keys %logics){
    print "Double count $k ".($logics{$k} / 2)."\n"; #divided by 2 because each pair counted twice
  }

 #print if a biotype has >2 transcripts in that cluster
  foreach my $k (keys %logics){
    if ($logics{$k} > 2){
      print "Eek - more than 1 pair of doubles for a single biotype\n";
      last;
    }
  }

  return keys %logics;
}


sub get_ordered_genes {
  my ($self, $genes) = @_;

  my @ordered;
  my %seen;  
  
  print "block has ".scalar(@$genes)." genes\n";
  foreach my $logic (@{$self->biotypes}) {
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
  my ($self, $cluster) = @_;
  my %protein_hash;

  # get all genes made from 1 protein
  
  # for this method, we want to use the actual protein
  # so, for cdna2genome, we use the protein id (which
  # is the second accession in the header line)
  GENE: foreach my $gene (@{ $cluster->get_Genes()} ) {
    my ($gene_ok, $protein) = check_gene($self, $gene);
    if (!$gene_ok) {
      warning("Gene check failed for gene ".$gene->dbID." (in cluster with 2 or more Sets) - removed from analysis");
      next GENE;
    }
    
    my $duplicate; 
    if (!exists $protein_hash{$protein}{$gene->biotype}) {
      push @{$protein_hash{$protein}{$gene->biotype}}, $gene;
    } else {
      GOT: foreach my $got_gene (@{$protein_hash{$protein}{$gene->biotype}}){
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
        push @{$protein_hash{$protein}{$gene->biotype}}, $gene;
      }    
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
  my ($self, $cluster) = @_;

  my %alternates;
  my %seen;

  GENE: foreach my $gene (@{ $cluster->get_Genes() }) {
    #print "Doing GENE ".$gene->dbID."...";
    if (exists $seen{$gene->dbID}) {
      #print "Already seen. Skipping...\n";
      next GENE;
    }
    # CHECK THAT:
    # each gene has only 1 transcript
    # each transcript translates
    # each transcript has one tsf
    my ($gene_ok, $protein) = check_gene($self, $gene);
    if (!$gene_ok) {
      throw("Gene check failed for gene ".$gene->dbID."");
    }
    
    # add to the hash
    push @{$alternates{$protein}{$gene->dbID}}, $gene;
    $seen{$gene->dbID} = 1;
    my $exon_string = make_exon_string(@{$gene->get_all_Transcripts}[0]);

    INNER: foreach my $inner_gene (@{ $cluster->get_Genes() }) {
      #print "Doing INNER gene ".$inner_gene->dbID."...";
      if (exists $seen{$inner_gene->dbID}) {
        #print "Already seen. Skipping...\n";
        next INNER;
      }
      # DO CHECKS as above
      my ($inner_gene_ok, $inner_protein) = check_gene($self, $inner_gene);
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


sub filter_transcripts {
  my ($h) = @_;
  my %hash = %{$h};
  my %hash2;
  #my %filtered = %{filter_transcripts(\%protein_clusters)};
  # fix doubles and overlappings
  foreach my $protein (keys %hash) {
    print "Filtering protein $protein\n";
    foreach my $biotype (keys %{$hash{$protein}}) {
      print "Filtering biotypes in random order $biotype\n";
      my @genes = @{$hash{$protein}{$biotype}};
      print "  started off with ".scalar(@genes)." genes\n";
      foreach my $g (@genes) {
        print "... gene ".$g->dbID."\n";
      }
      if (scalar(@genes) > 1) {
        print "We have doubles or overlapping\n";
        my @genes2 = @{tidy_up($hash{$protein}{$biotype})};
        print "  now have ".scalar(@genes2)." genes left\n";
        push @{$hash2{$protein}{$biotype}}, @{tidy_up($hash{$protein}{$biotype})}; 
      } else {
        push @{$hash2{$protein}{$biotype}},  @genes;
      }
    }
  }
  return \%hash2;
}

sub tidy_up {
  my ($genes) = @_;

  print "In tidy up...have ".scalar(@$genes)." genes\n";
  my @ordered = sort {($a->end - $a->start) <=> ($b->end - $b->start)} @$genes;

  # make an array holding the shortest gene
  my @non_overlapping;
  push @non_overlapping, $ordered[0];

  # add in genes if they don't overlap
  my $overlaps = 0;
  OUTER: foreach my $outer (@ordered) {
    INNER: foreach my $inner (@non_overlapping) {
      if ($outer->dbID == $inner->dbID) {
        next OUTER;
      }
      if (($outer->start <= $inner->end && $outer->end >= $inner->start) 
       || ($inner->start <= $outer->end && $inner->end >= $outer->start)) {
        print "overlaps... checking if it extends the smaller gene... ";
        my $outer_transc = $outer->get_all_Transcripts->[0];
        my $inner_transc = $inner->get_all_Transcripts->[0];
        #my $outer_transc_seq = $outer->get_all_Transcripts->[0]->translateable_seq;
        #my $inner_transc_seq = $inner->get_all_Transcripts->[0]->translateable_seq;
        my $exact_match = 0;
        for (my $i=0; $i< scalar(@{$outer_transc->get_all_Exons}); $i++) {
          for (my $j=0; $j<scalar(@{$inner_transc->get_all_Exons}); $j++) {
            if ($outer_transc->get_all_Exons->[$i]->start == $inner_transc->get_all_Exons->[$j]->start &&
                $outer_transc->get_all_Exons->[$i]->end == $inner_transc->get_all_Exons->[$j]->end) {
              $exact_match++;
              #next;
            }
          }
        }
        if ($exact_match == scalar(@{$inner_transc->get_all_Exons})) { 
        #if ($outer_transc->translateable_seq =~ m/$inner_transc_seq/) {
          print "EXTENDS\n";
        } else {
          print "OVERLAP\n";
          $overlaps++;
          next OUTER;
        }
      }
    }
    if ($overlaps < 1) {
      print "ADDING A DOUBLE\n";
      push @non_overlapping, $outer;
    }
  }

  # now, we might be able 
  return \@non_overlapping;
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
  my ($self, $gene) = @_;

  my $pass_checks = 1;
  my $hit_name;
  my @transcripts = @{$gene->get_all_Transcripts};
  my $protein_acc;

  if (scalar(@transcripts) != 1) {
    print STDERR "Gene with dbID ".$gene->dbID." has ".scalar(@transcripts)." transcripts\n";
    $pass_checks = 0;
  }

  # This bit requires a bit of tweaking for cdna2genome:
  # the genes from cdna2genome have cdnas as transcript_supporting_evidence
  # so will be dna_align_features with a hseqname=cdna_id. 
  # For targetted run, tsfs are protein_align_features
  # with hseqname=protein_id. 
  foreach my $transcript (@transcripts) {
    my @tsfs = @{$transcript->get_all_supporting_features};
    if (scalar(@tsfs) == 1) {
      $hit_name = $tsfs[0]->hseqname; 
      my $entry_obj1 = $self->seqfetcher->get_entry_by_acc($hit_name);
      #print "PRIMARY primary $hit_name ".$entry_obj1."\n";

     # my @namesp = @{$self->seqfetcher->secondary_namespaces};
     # if (scalar(@namesp ) >0) {
     #   foreach my $sn (@namesp) {
     #     print "Namespace '$sn'\n";
     #   }
     # } 
     
      #my $secondary_header_id;
      if ($entry_obj1 =~m/^\>(\S+)\n/) { 
        $protein_acc = $1; 
      } elsif ($entry_obj1 =~m/^\>(\S+\.\d+)\s+\n/) { 
        $protein_acc = $1; 
      } elsif ($entry_obj1 =~m/^\>(\S+)\s+(\S+)\n/) {
      #  $secondary_header_id = $2;
        $protein_acc = $2;
      } else {
        print "Couldn't find protein acc for cDNA $hit_name. Input ID ".$self->query->name()." skipped.\n";
        throw("Entry has unusual header : ID '$entry_obj1' for hit_name '$hit_name'. Check why ID is an empty string.");
      }
      if ( length($protein_acc) > 50 ) {  
        warning("unusally long acc-number. are you sure this is an acc : $protein_acc ??? \n") ; 
      } 
      #if ($secondary_header_id) {
      #  print STDERR "fetching secondary seq with 'ID' and acc $secondary_header_id\n\n";
      #  my @seq_objs = $self->seqfetcher->get_Seq_by_secondary('ID', $secondary_header_id);
      #  foreach my $seq_obj2 (@seq_objs) {
      #    print "SECONDARY SEQ $hit_name $secondary_header_id\n ".$seq_obj2->seq."\n";
      #  }
      #}

      #if ($secondary_header_id) {
      #  my $sseq_obj = $self->seqfetcher->get_Seq_by_acc($secondary_header_id);
      #  print ">$secondary_header_id\n".$sseq_obj->seq."\n";
      #}

      # print "index name ".$self->seqfetcher->index_name."\n"; # prints /lustre/work1/ensembl/ba1/cow4/Seq/BestTargetted_proteome_final/proteome.fa
       

      # check for a translation
      if (!defined($transcript->translation)) {
        print STDERR "Transcript with dbID ".$transcript->dbID." does not translate.\n";
        $pass_checks = 0;
      }
      else {
        my ($original_seq) = $entry_obj1 =~ /^>\S+.*\n(\w+)/;
        my ($coverage, $identity) = align_proteins($original_seq, $transcript->translation->seq);
        warning("$coverage $identity");
        if ($coverage < $self->min_coverage and $identity < $self->min_identity) {
          print STDERR "Transcript with dbID ".$transcript->dbID." has low coverage $coverage and identity $identity\n";
          $pass_checks = 0;
        }
      }
    } else {
      print STDERR "Transcript with dbID ".$transcript->dbID.
           " has ".scalar(@tsfs)." tsfs. ".
           "Should only have one transcript_supporting_feature\n";
      $pass_checks = 0;
    }
  }
  return $pass_checks, $protein_acc;
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

sub biotypes {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_bt_biotypes} = $val;
  }

  return $self->{_bt_biotypes};
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

sub keep_single_analysis {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_bt_keep_single_analysis} = $val;
  }

  return $self->{_bt_keep_single_analysis};
}

sub cluster_on_coding_exons {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_bt_cluster_on_coding_exons} = $val;
  }

  return $self->{_bt_cluster_on_coding_exons};
}

sub min_coverage {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_bt_min_coverage} = $val;
  }

  return $self->{_bt_min_coverage};
}

sub min_identity {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_bt_min_identity} = $val;
  }

  return $self->{_bt_min_identity};
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
      print "exonerate results @tmp ";
      print "\n";
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
