package Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;

use Exporter;
use vars qw(@ISA @EXPORT);
use strict;
use Carp;

use Bio::EnsEMBL::Analysis::Tools::Algorithms::ExonCluster;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster ; 

use Bio::EnsEMBL::Utils::Exception qw (warning throw ) ; 
@ISA=qw(Exporter);

@EXPORT=qw(
           cluster_exons_in_transcript_cluster 
           genes_to_Transcript_Cluster
           );



 



sub cluster_exons_in_transcript_cluster { 
  my ($tc) = @_;
  my @clusters;
  
  foreach my $trans (@{$tc->get_Transcripts}) {
    my $tr_biotype = $trans->biotype; 

    foreach my $exon (@{$trans->get_all_Exons}) {

      my @matching_clusters;
       print "\nExon " . $exon->dbID . " limits: " . $exon->start . 
       " and " .  $exon->end . "\n";
 
      CLUSTER: foreach my $cluster (@clusters) {
         #print "Testing against cluster with limits " . 
         #$cluster->start. " to " . $cluster->end . " ".$cluster->strand ."\t";
        if (!($exon->start >= $cluster->end ||
              $exon->end <= $cluster->start)) {
           if ($cluster->strand eq $exon->strand ){ 
              push (@matching_clusters, $cluster);  
              print "cl. matches " .$cluster->strand ."\t" .$exon->strand . "\t" .$trans->strand ."\n" ;
           }
        }
        #print "\n";
      }
      if (scalar(@matching_clusters) == 0) {
        # print STDERR "Created new cluster for " . $exon->stable_id . " " . $exon->dbID . "\n";
        # print "\ncreating new cluster for Exon " . $exon->dbID . 
        # " limits: " . $exon->start . " and " .  $exon->end . "\n";
        
        my $newcluster = Bio::EnsEMBL::Analysis::Tools::Algorithms::ExonCluster->new() ; 
       
        $newcluster->add_exon($exon,$trans);
        push(@clusters,$newcluster);
 
      } elsif (scalar(@matching_clusters) == 1) {
         #print STDERR "Adding to cluster for " . $exon->stable_id . " " . $exon->dbID . "\n";
        $matching_clusters[0]->add_exon($exon,$trans);
      } else {
         # Merge the matching clusters into a single cluster
        print STDERR "Merging clusters for " . $exon->dbID ."\n";
        my @new_clusters;
        my $merged_cluster = Bio::EnsEMBL::Analysis::Tools::Algorithms::ExonCluster->new() ; 

        foreach my $clust (@matching_clusters) {
          $merged_cluster->merge($clust);
        }
        $merged_cluster->add_exon($exon,$trans);
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
  }
  print "UtilsFunc.pm : Finished Exon-clustering -- having " .scalar(@clusters) . " exon clusters \n" ; 

  # setting exon/cluster relationship 
  for my $c (@clusters) {
    for my $e(@{ $c->get_all_Exons_in_ExonCluster} ) { 
      $e->cluster($c) ; 
    }
  }
  return @clusters;
}    



sub genes_to_Transcript_Cluster {
  my ($genes_or_predTrans) = @_;
  
  my $tc = Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster->new() ; 
  print "building new TranscriptCluster\n" ; 
  foreach my $gene (@$genes_or_predTrans) {
    if( ref($gene)=~m/Gene/){ 
      # is a Bio::EnsEMBL::Gene 
      foreach my $trans (@{$gene->get_all_Transcripts}) {
        if ($gene->strand ne $trans->strand ) { 
          throw("Weird - gene is on other strand than transcript\n") ; 
        }
        for (@{ $trans->get_all_Exons} ) {
           if ($_->strand ne $trans->strand ) { 
             print $trans->biotype . " " . $trans->seq_region_start . " "  . $trans->seq_region_end . " " . $trans->seq_region_strand ."\n" ; 
             print $_->biotype . " " . $_->seq_region_start . " "  . $_->seq_region_end . " " .$_->seq_region_strand . "\n"; 
             throw("Weird - exon is on other strand than transcript\n") ; 
           } 
        }
        # assure that transcript has same biotype as gene 
        $trans->biotype($gene->biotype) ; 
        $trans->sort;
        #print "Adding transcript " . $trans->stable_id . "\n";
        $tc->put_Transcripts($trans);
        $tc->register_biotype($gene->biotype) ; 
      }
    } else {
      # is not a Bio::EnsEMBL::Gene 
      warning("Not having a Bio::EnsEMBL::Gene-object : clustering $gene\n") ; 
      $gene->sort ; 
      $tc->put_Transcripts($gene) ; 
    }
  }
  return $tc;
}




1;
