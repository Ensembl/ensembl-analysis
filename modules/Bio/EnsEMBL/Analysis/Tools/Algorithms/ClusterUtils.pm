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

  Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;


  This package contains methods to make the GeneClusster easier to use.  

=head1 DESCRIPTION  



 Here's some ASCII art to understand this module better : 
 
   We will refer to 2 types of clusters in this module :  

    'one_way clusters'    - a cluster of genes [= genes which have overlapping exons ]
                            which contain only genes of one set 
                            use get_single_clusters method

    'two_way clusters'    - a cluster of genes [= genes which have overlapping exons ]
                            which contains genes of both sets, 
                            use get_twoway_clusters method

  Maybe it would have been better to call them one-set clusters and two-set clusters?



  'unclustered'         - these clusters contain genes which don't overlap any other gene 
                          in the same set, or in any other set. These genes are also known as 
                          non-clustering genes 



Example 1 :  

  This is a two-way cluster ( two-set cluster ) 
    - genes of 2 different sets overlap 

  SET 1 :          AAAAAAAAAAAA----------------AAAAAA-AAAAAAAAAAAA       
  SET 1        BBBBBBBB------------------------BBBB 
  SET 2 :    ZZZZZZZZZZZ-----------------------------------------------ZZZZZZZZZZZZ



Example 2: 

  Below you see 2 one-way clustering clusters  (2-set-clusters)
    - both clusters are called 'one-way clustering clusters' - they only cluster 'one way'  - they are homogenous clusters...
    - no cluster contains a gene which belongs to a different set 
      method : #  get_oneway_clustering_genes_of_set($clustered,"transformed")

  SET 1 :          AAAAAAAAAAAA----------------AAAAAA-AAAAAAAAAAAA       
  SET 1 :       BBBBBBBB------------------------BBBB 
  SET 2 :                                                                  ZZZZZZZZZZZZZ----------------ZZZZZZZZZZZZZ
  SET 2 :                                                                      YYYYYYYYYYYY


Example 3 : 

  All the genes below form one two-way cluster because the exons overlap : 

  SET 1 :          AAAAAA--------------------------------AAAAAAAAAAAA       
  SET 1 :       BBBBB                                                                                    BBBBBBBBBBBB-----------BBBBBBB
  SET 2 :                                                   ZZZZZZZZZZZZZZZZ--------------ZZZZZZZZZZZZZZZZZZ 



Example 4 : 

  The genes below are one-way-clustering : ( one set cluster ) 
      method : #  get_oneway_clustering_genes_of_set($clustered,"transformed")

  SET 1 :          AAAAAAAAAAAA----------------AAAAAA-AAAAAAAAAAAA       
  SET 1 :       BBBBBBBB------------------------BBBB 




Some definitions : 


    one-way clustering means : 
        - genes only have exon overlap within the same set, they only cluster among each other ('homogenous clusters')
          ( see Example 4 and Example 2 ) 

    two-way clustering means (=two set cluster )
       - the cluster contains genes of set A and set B, and the genes have exon overlap 
       -  example 1 and 3 



=cut

package Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;


use warnings ;
use Exporter;
use vars qw(@ISA @EXPORT);
use strict;
use Carp;

use Bio::EnsEMBL::Analysis::Tools::Algorithms::ExonCluster;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster ; 
use Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster ; 
use Bio::EnsEMBL::Analysis::Tools::Algorithms::AlignFeatureCluster ; 

use Bio::EnsEMBL::Utils::Exception qw (warning throw ) ; 
@ISA=qw(Exporter);

@EXPORT = qw(
  simple_cluster_Genes
  cluster_Genes
  cluster_Genes_without_strand
  cluster_Genes_by_coding_exon_overlap
  get_single_clusters
  get_twoway_clusters
  get_twoway_clustering_genes_of_set
  get_oneway_clustering_genes_of_set
  make_types_hash
  make_types_hash_with_genes
  cluster_AlignFeatures
);



=head2 make_types_hash 

   Arg[1]    : Array ref. to gene set 1 
   Arg[2]    : Array ref. to gene set 2 
   Arg[3]    : Name of gene set 1 
   Arg[4]    : Name of gene set 2 

   Function  : Creates a type hash to contain the names of the two sets

   Returnval : Hashreference of types to be used
 

=cut


sub make_types_hash {
    my ( $gene_set1,$gene_set2,$gene_set1_name, $gene_set2_name  ) = @_;

    my (%types_hash,%types_1, %types_2 ) ;

    unless ( $gene_set1_name ) { 
      $gene_set1_name = "gene_biotypes_set1" ;
    }

    unless ( $gene_set2_name ) {
      $gene_set2_name = "gene_biotypes_set2" ;
    }

    @types_1{ map { $_->biotype } @{ $gene_set1 } } = 1 ;
    $types_hash{$gene_set1_name} = [ keys %types_1 ]  ;

    @types_2{ map { $_->biotype } @{ $gene_set2 } } = 1 ;
    $types_hash{$gene_set2_name} = [ keys %types_2 ]  ; 

    #my @tmp  = map { $_->biotype } @{ $gene_set2 };  
    #for ( @tmp ) { 
    #  print join("\n", $_ ) ; 
    #}

    my @intersection ;
    for ( keys %types_1 ) {
       if ( exists $types_2{$_} ) {
         push @intersection, $_ ;
       }
    }
    if ( @intersection > 0 ) {
      warn("There are biotypes you try to cluster which are in both gene sets. This is fine is it is expected\n")  ;
    }
    return \%types_hash;
} 




=head2 make_types_hash_with_genes 

   Arg[1]    : Array ref. to gene set 1 
   Arg[2]    : Array ref. to gene set 2 
   Arg[3]    : Name of gene set 1 
   Arg[4]    : Name of gene set 2 

   Function  : Creates a type hash to contain the names of the two sets, along with an arrayreference of the genes

   Returnval : Hashreference of types used and Arrayreference of gene objects
 

=cut



sub make_types_hash_with_genes { 
    my ( $gene_set1,$gene_set2,$gene_set1_name, $gene_set2_name  ) = @_;
  
   my $types_ref_hash = make_types_hash( $gene_set1,$gene_set2,$gene_set1_name, $gene_set2_name  ) ; 
   my @all_genes_combined = ( @$gene_set1, @$gene_set2);
   return [ $types_ref_hash, \@all_genes_combined ]; 
}


=head2 simple_cluster_Genes 

   Arg[1]    : Array ref. to gene set 1 
   Arg[2]    : Name of gene set 1 
   Arg[3]    : Array ref. to gene set 2 
   Arg[4]    : Name of gene set 2 

   Function  : Simply clusters two sets of genes 
               

   Returnval : Arrayreferenc of  Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster which can contain genes of type set1 , or set2, 
               or both. 
 

=cut


sub simple_cluster_Genes {
  my ( $gene_set1, $gene_set1_name, $gene_set2, $gene_set2_name ) = @_;

  my ($types_hash,$all_genes) = @{make_types_hash_with_genes($gene_set1, $gene_set2, $gene_set1_name, $gene_set2_name)};
  my ($clustered,$unclustered) = cluster_Genes($all_genes, $types_hash);

  return [ $clustered, $unclustered ];
}




=head2 get_single_clusters ( should be called get_oneway_clusters ) 

   Arg[1]    : Array reference to  Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster objects  

   Function  : Filters out all clusters of the array which are twoway-clusters [ which contain set1 and set2 ] 
               This means all clusters, which contain genes of set 1 and set 2, are not returned. 
               Only clusters are returned which contain only one type of set, either set1 or set2          

   Returnval : Arrayreferenc of  Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster which can contain genes of type set1,
               or set2, or both. 
 

=cut



sub get_single_clusters {  
  my ( $cluster_ref ) = @_ ; 
 
  my @single_type_cluster;   

   check_cluster_ref($cluster_ref) ; 

   for my $c ( @$cluster_ref ) { 
     unless ( $c->is_twoway_cluster ) {  
        push @single_type_cluster , $c ; 
     }    
   } 
   return \@single_type_cluster ; 
} 



=head2 get_twoway_clusters 

   Arg[1]    : Array reference to  Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster objects  

   Function  : Out of a given array of GeneCluster objects, only those ones are returned which 
               contain genes of type set1 AND set2. 

   Returnval : Arrayreference of  Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster which are twoway
 
=cut


sub get_twoway_clusters {  
   my ($cluster_ref ) = @_ ;   

   check_cluster_ref($cluster_ref) ; 
   my @tw;  


   for my $c ( @$cluster_ref ) { 
     if ( $c->is_twoway_cluster ) {  
        push @tw, $c ; 
     }    
   }
   return \@tw; 
}  


=head2 get_twoway_clustering_genes_of_set

   Arg[1]    : Array reference to  Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster objects and a set name 

   Function  : Out of a given array of GeneCluster objects, clusters which are twoways are selected.
               Then, only genes of the specified target_set_name are returned. 

   Returnval : Arrayreferenc of  Bio::EnsEMBL::Analysis::Gene objects which belong to target set name and are in twoway clusters
 
=cut




sub get_twoway_clustering_genes_of_set {
   my ($cluster_ref,$target_set_name ) = @_ ;

   check_cluster_ref($cluster_ref) ;
   my @twoway_clustering_genes ;

   for my $twoway_cluster (@{ get_twoway_clusters($cluster_ref)} ) {
       push @twoway_clustering_genes, @{ $twoway_cluster->get_Genes_by_Set($target_set_name)};
   }
   if ( scalar(@twoway_clustering_genes) == 0 ) {
     warning (" no gene of set-type \"$target_set_name\" found - i only know these sets : " . join (",", @{ get_sets_included($cluster_ref)} )) ;
   }
   return \@twoway_clustering_genes;
}




=head2 get_oneway_clustering_genes_of_set 

   Arg[1]    : Array reference to Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster objects  

   Function  : Out of the given array of Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster object, 
               clusters which are two-way-cluster ( they contain genes of set 1 AND set 2 ) are excluded. 
               Then, only genes of the specified target_set_name are returned. 

   Returnval : Arrayreference of Bio::EnsEMBL::Gene objects which belong to set $target_set_name and are in oneway_clusters. 


=cut


sub get_oneway_clustering_genes_of_set {  
   my ($cluster_ref,$target_set_name ) = @_ ;    

   check_cluster_ref($cluster_ref) ;  
   my @single_clustering_genes ; 
   my $cnt=0;
   for my $single_type_cluster (@{ get_single_clusters($cluster_ref)} ) {   
       push @single_clustering_genes, @{ $single_type_cluster->get_Genes_by_Set($target_set_name) }; 
      $cnt++; 
   }  
 
   if ( scalar(@single_clustering_genes) == 0 ) {  
     warning (" no gene of set-type \"$target_set_name\" found in $cnt out of " . @$cluster_ref . " clusters }- i only know these sets : " . join (",", @{ get_sets_included($cluster_ref)} )) ;
    my %all_set; 
    for my $cr ( @$cluster_ref ) {  
       @all_set{@{ $cr->get_sets_included($cr)}}=1;
    }  
   }   

   return \@single_clustering_genes; 
} 


sub get_sets_included {  
   my ($cluster_ref)  = @_ ;     
   my %tmp ; 
   for my $c ( @$cluster_ref) {  
      @tmp{ @{$c->get_sets_included} } = 1 ; 
   }
   return [ keys %tmp ] ; 
}



sub check_cluster_ref {   
  my ( $cluster_ref ) = @_ ;  
  unless ( ref($cluster_ref) =~m/ARRAY/ ) { 
     throw "Need to hand over array-ref" ; 
  }  
} 




=head2 cluster_Genes_by_coding_exon_overlap

   Arg[1]    : Aref to Bio::EnsEMBL::Gene objects 
   Arg[2]    : ref to hash which builds up an Evidence-Set to biotype-relation :

   Function  : clusters all genes in Arg[1] according to their coding exon extent and 
              sets the type according to their sets
              Uses the cluster_Genes method but clusters genes together only if there is coding exon overlap

   Returnval : Array of 2 Arrayrefs : First arrayref holds an array of all genes which 
               have been clustered together, second ref holds an array of genes 
               which were not clustered.   

=cut



sub cluster_Genes_by_coding_exon_overlap {  
  my ($genes, $types_hash ) = @_ ;
  return cluster_Genes($genes,$types_hash,1); 
}




=head2 cluster_Genes_without_strand

   Arg[1]    : Aref to Bio::EnsEMBL::Gene objects 
   Arg[2]    : ref to hash which builds up an Evidence-Set to biotype-relation :

   Function  : clusters all genes in Arg[1] according to their genomic extent and 
              sets the type according to their sets
              Uses the cluster_Genes method but genes cluster together even if there are on different strands

   Returnval : Array of 2 Arrayrefs : First arrayref holds an array of all genes which 
               have been clustered together, second ref holds an array of genes 
               which were not clustered.   

=cut



sub cluster_Genes_without_strand {  
  my ($genes, $types_hash ) = @_; 
  return cluster_Genes($genes,$types_hash,0,1); 
}



=head2 cluster_Genes 

   Arg[1]    : Aref to Bio::EnsEMBL::Gene objects 
   Arg[2]    : ref to hash which builds up an Evidence-Set to biotype-relation :
               $hash{ cdna } = ['cdna_kyoto', 'cdna_other' ] 
               $hash{ simg } = ['simgw_100','simgw_200'] 
   Arg[3]    : flag if you want to cluster by overlap of protein_coding exons only
               ( 1 = clustering my overlap of protein_coding exons only )  

   Arg[4]    : flag if you want to cluster without strand information 
               ( if you ignore the strand, overlapping genes on diff strand will end up in same cluster )  

   Arg[5]    : flag if you want to cluster without exon information 
               ( if you ignore the exon, it only take the start and end of the gene)  

   Function  : clusters all genes in Arg[1] according to their genomic extent and 
              sets the type according to their sets

   Returnval : Array of 2 Arrayrefs : First arrayref holds an array of all genes which 
               have been clustered together, second ref holds an array of genes 
               which were not clustered.  

=cut



sub cluster_Genes {
  my ($genes, $types_hash, $check_coding_overlap, $ignore_strand, $ignore_exon_overlap) = @_ ;


  #print "GOT " . scalar(@$genes ) . " GENES tocluster \n" ; sleep(2) ; 

  #
  # steves old cluster-routine clusters genes of two types : 'ncbi' and 'hinxton' 
  # ( see get_twoay_cluster.pl) 
  # he uses  two gene-sets : genes and compare_genes (each set may contain differnt biotypes)
  #   
  # he uses the sets to see if a cluster only consists of genes out of one set (ncbi) or hinxton 
  # and retrieves all sets of a cluster with "get_sets_included" 
  #
  # we do something 'nearly similar : we are clustering genes of diffrent sets (simgw, est, abinitio) 
  # and have methods to access these sets 
  # --> GeneCluster has methods get_Genes_of_Type / get_Genes_by_Type / get_Genes_by_Set
  # all genes on slice are handed over and a %types_hash which holds the setname and the  

  return ([],[]) if (!scalar(@$genes));

  # sorting of ALL genes
  my @sorted_genes = sort { $a->start <=> $b->start ? $a->start <=> $b->start  : $b->end <=> $a->end }  @$genes;

  #print STDERR "Clustering ".scalar( @sorted_genes )." genes on slice\n" ;  
  # select count(*) , g.biotype from gene g left join transcript t on g.gene_id = t.gene_id where isnull(t.gene_id ) group by g.biotype  ;
  for ( @sorted_genes ) {  
     my $tr = scalar ( @{$_->get_all_Transcripts} ) ;   
    if ( $tr == 0 ) {
      throw("data error - gene with gene_id " . $_->dbID ." biotype " . $_->biotype .  " does not have any associated transcripts \n") ;  
    }
  } 
  my $count = 0;

  my @active_clusters;
  my @inactive_clusters;
  TRANSCRIPT: foreach my $gene (@sorted_genes) {
    $count++;

    # Every 50 genes divide clusters into an active (ones which could be have this gene added to it) and inactive 
    # (ones which can not be altered by any of the genes to come)
    if (!($count%50)) { 
      my @still_active_clusters;
      my $gene_start = $gene->start;
      foreach my $cluster (@active_clusters) {
        if ($cluster->end < $gene_start) {
          push @inactive_clusters,$cluster;
          #print "Cluster inactive\n";
        } else {
          push @still_active_clusters,$cluster;
        }
      }
      @active_clusters = @still_active_clusters;
    }

    my @matching_clusters;

    ## 
    ## if there are Clusters (initialisation below) than check  
    ## if gene lies in the boundaries of the cluster and has at least 
    ## one exon which overlaps with an exon of a gene which already 
    ## belongs to the cluster
    ##

    CLUSTER: foreach my $cluster (@active_clusters) {

    # 
    # if gene lies in the boundaries of the cluster......
    #

      if ($gene->end  >= $cluster->start && $gene->start <= $cluster->end) {

        # search for a gene in the cluster which overlaps the new gene, 

        foreach my $cluster_gene (@{ $cluster->get_Genes} ){

        # check if clustered gene overlaps 

          if ($gene->end  >= $cluster_gene->start && $gene->start <= $cluster_gene->end) {

            #                             CASE 1: 
            #
            #         START----------------$cluster_gene---------------END 
            # START-------------$gene-----------------------END
            #
            #                             CASE 2 : 
            #                              
            #         START----------------$cluster_gene---------------END 
            #               START-------------$gene-----------END
            #
            #                             CASE 3 : 
            #
            #
            #         START----------------$cluster_gene----------END 
            #                                              START------$gene-------END
            #
            # add gene target-gene to cluster if it has at least
            # one gene which overlaps with an exon of the clustered gene 
            # and add to cluster  
            #

            if ($ignore_exon_overlap) {
              if (!$ignore_strand) {
                if ($gene->strand == $cluster_gene->strand) {
                  push (@matching_clusters, $cluster);
                  next CLUSTER;
                }
              } else {
                push (@matching_clusters, $cluster);
                next CLUSTER;
              }
            }
            elsif (_compare_Genes( $gene, $cluster_gene, $check_coding_overlap, $ignore_strand)) {
              push (@matching_clusters, $cluster);
              next CLUSTER;
            }
          }
        }
      }
    } # CLUSTER 

    ##
    ## Initialization of we have no matching cluster (above) 
    ###############################################################

    #
    # if above was found NO matching cluster
    # than make a new one 
    #  
    if (scalar(@matching_clusters) == 0) {
      my $newcluster = Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster->new($ignore_strand);
      foreach my $set_name (keys %$types_hash) { 
        $newcluster->gene_Types($set_name,$types_hash->{$set_name});
      }  
      $newcluster->put_Genes([$gene], $ignore_strand); # xx
      push(@active_clusters,$newcluster);

      #
      # if above was found ONE matching cluster
      #
    } elsif (scalar(@matching_clusters) == 1) {
      $matching_clusters[0]->put_Genes([$gene], $ignore_strand);

    } else {
      # Merge the matching clusters into a single cluster
      my @new_clusters;
      my $merged_cluster = Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster->new($ignore_strand);

      foreach my $set_name (keys %$types_hash) {
        $merged_cluster->gene_Types($set_name,$types_hash->{$set_name});
      }

      my %match_cluster_hash;
      foreach my $clust (@matching_clusters) {
        $merged_cluster->put_Genes ($clust->get_Genes, $ignore_strand);
        $match_cluster_hash{$clust} = $clust;
      }
      $merged_cluster->put_Genes([$gene], $ignore_strand);
      push @new_clusters,$merged_cluster;

      # Add back non matching clusters
      foreach my $clust (@active_clusters) {
        if (!exists($match_cluster_hash{$clust})) {
          push @new_clusters,$clust;
        }
      }
      @active_clusters =  @new_clusters;
    }
  }
  # Separate genes which are UNclustered (only one gene in cluster) and
  # from clusters which hold more than one gene 

  # print "Have " . scalar(@active_clusters) . " active clusters and " . scalar(@inactive_clusters) . " inactive clusters\n";
#  my @clusters = (@active_clusters,@inactive_clusters);

  my (@new_clusters, @unclustered);
  foreach my $cl (@active_clusters,@inactive_clusters){
    if ( $cl->get_Gene_Count == 1 ){
      push @unclustered, $cl;
    } else{
      push( @new_clusters, $cl );
    }
  }
 # print STDERR "All Genes clustered\nGot " . scalar(@new_clusters) . " new Clusters\n"  ;
  
  return (\@new_clusters, \@unclustered);
}



=head2  _compare_Genes


  Title  :  _compare_Genes
  Usage  :  this internal function compares the exons of two genes on overlap
  Source :  Bio::EnsEMBL::Pipeline::GeneComparison::GeneComparison 


=cut



sub _compare_Genes {
  my ($gene1,$gene2,$translate, $ignore_strand) = @_;
  
  # quit if genes do not have genomic overlap 
  #
  # start-------gene1------end   start--------gene2----------end
  #  
  
  
  
  # $overlaps = ( $exon1->end >= $exon2->start && $exon1->start <= $exon2-> end );  

  if ($ignore_strand or $gene1->strand == $gene2->strand) {
    if ($translate) {
      #print "clustering by overlap of coding exons only\n";
      # exon-overlap only on coding exons !
      my $exons1 = get_coding_exons_for_gene($gene1);
      my $exons2 = get_coding_exons_for_gene($gene2);
      foreach my $exon1 (@$exons1) {
        foreach my $exon2 (@$exons2) {
          if ($exon1->overlaps_local($exon2)){
            return 1;
          }
        }
      }
    } else {
      #
      # overlap check based on all (noncoding + coding) Exons
      #
      foreach my $exon1 (@{$gene1->get_all_Exons}){
        foreach my $exon2 (@{$gene2->get_all_Exons}){
          if ($exon1->overlaps_local($exon2)){
            #print "Passed exon overlap check (noncod. + cod. exons checked)  - returning 1\n";
            return 1;
          }
        }
      }
    }
  }
   #print "Failed overlap check (translate = $translate) - returning 0\n";
  return 0;
}


sub get_coding_exons_for_gene {
  my ($gene) = @_;

  my @coding;

  foreach my $trans (@{$gene->get_all_Transcripts}) {
    next if (!$trans->translation);
    foreach my $exon (@{$trans->get_all_translateable_Exons}) {
      push @coding, $exon;
    }
  }

  return \@coding; 
}

sub cluster_AlignFeatures {
  my ($features, $types_hash, $ignore_strand) = @_ ;

  print "GOT " . scalar(@$features ) . " FEATURES to cluster \n" ; sleep(2) ; 
  return ([],[]) if (!scalar(@$features));

  # sorting of ALL features
  my @sorted_features = sort { $a->start <=> $b->start ? $a->start <=> $b->start  : $b->end <=> $a->end }  @$features;

  # now start working
  my $count = 0;
  my @active_clusters;
  my @inactive_clusters;

  FEATURE: foreach my $feature (@sorted_features) {
    $count++;

    # Every 50 features divide clusters into an active (ones which could be have this feature added to it) and inactive 
    # (ones which can not be altered by any of the features to come)
    if (!($count%50)) { 
      my @still_active_clusters;
      my $feature_start = $feature->start;
      foreach my $cluster (@active_clusters) {
        if ($cluster->end < $feature_start) {
          push @inactive_clusters,$cluster;
        } else {
          push @still_active_clusters,$cluster;
        }
      }
      @active_clusters = @still_active_clusters;
    }

    my @matching_clusters;

    ## 
    ## if there are Clusters (initialisation below) then check  
    ## if feature lies in the boundaries of the cluster and has at least 
    ## one exon which overlaps with an exon of a feature which already 
    ## belongs to the cluster
    ##

    CLUSTER: foreach my $cluster (@active_clusters) {
    # if feature lies in the boundaries of the cluster......
      if ($feature->end  >= $cluster->start && $feature->start <= $cluster->end) {
        # search for a feature in the cluster which overlaps the new feature, 
        foreach my $cluster_feature (@{ $cluster->get_AlignFeatures} ){
          # check if clustered feature overlaps 
          if ($feature->end  >= $cluster_feature->start && $feature->start <= $cluster_feature->end) {

            #                             CASE 1: 
            #
            #         START----------------$cluster_feature---------------END 
            # START-------------$feature-----------------------END
            #
            #                             CASE 2 : 
            #                              
            #         START----------------$cluster_feature---------------END 
            #               START-------------$feature-----------END
            #
            #                             CASE 3 : 
            #
            #
            #         START----------------$cluster_feature----------END 
            #                                              START------$feature-------END
            #

            if (!$ignore_strand) {
              if ($feature->strand == $cluster_feature->strand) {
                push (@matching_clusters, $cluster);
                next CLUSTER;
              }
            } else {
              push (@matching_clusters, $cluster);
              next CLUSTER;
            }
          }
        }
      }
    } # CLUSTER 

    ##
    ## Initialization of we have no matching cluster (above) 
    ###############################################################

    #
    # if above was found NO matching cluster
    # than make a new one 
    #  
    if (scalar(@matching_clusters) == 0) {
      my $newcluster = Bio::EnsEMBL::Analysis::Tools::Algorithms::AlignFeatureCluster->new($ignore_strand);
      foreach my $set_name (keys %$types_hash) { 
        $newcluster->feature_Types($set_name,$types_hash->{$set_name});
      }  
      $newcluster->put_AlignFeatures([$feature], $ignore_strand); # xx
      push(@active_clusters,$newcluster);

      #
      # if above was found ONE matching cluster
      #
    } elsif (scalar(@matching_clusters) == 1) {
      $matching_clusters[0]->put_AlignFeatures([$feature], $ignore_strand);

    } else {
      # Merge the matching clusters into a single cluster
      my @new_clusters;
      my $merged_cluster = Bio::EnsEMBL::Analysis::Tools::Algorithms::AlignFeatureCluster->new($ignore_strand);

      foreach my $set_name (keys %$types_hash) {
        $merged_cluster->feature_Types($set_name,$types_hash->{$set_name});
      }

      my %match_cluster_hash;
      foreach my $clust (@matching_clusters) {
        $merged_cluster->put_AlignFeatures ($clust->get_AlignFeatures, $ignore_strand);
        $match_cluster_hash{$clust} = $clust;
      }
      $merged_cluster->put_AlignFeatures([$feature], $ignore_strand);
      push @new_clusters,$merged_cluster;

      # Add back non matching clusters
      foreach my $clust (@active_clusters) {
        if (!exists($match_cluster_hash{$clust})) {
          push @new_clusters,$clust;
        }
      }
      @active_clusters =  @new_clusters;
    }
  }
  # Separate features which are UNclustered (only one feature in cluster) and
  # from clusters which hold more than one feature 

  print "Have " . scalar(@active_clusters) . " active clusters and " . scalar(@inactive_clusters) . " inactive clusters\n";
  my @clusters = (@active_clusters,@inactive_clusters);

  my (@new_clusters, @unclustered);
  foreach my $cl (@clusters){
    if ( $cl->get_AlignFeature_Count == 1 ){
      push @unclustered, $cl;
    } else{
      push( @new_clusters, $cl );
    }
  }
 # print STDERR "All AlignFeatures clustered\nGot " . scalar(@new_clusters) . " new Clusters\n"  ;
  
  return (\@new_clusters, \@unclustered);
}

1;
