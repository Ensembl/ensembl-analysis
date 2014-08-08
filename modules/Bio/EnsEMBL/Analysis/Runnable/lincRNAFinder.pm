=head1 LICENSE

# Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Runnable::lincRNAFinder - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut
package Bio::EnsEMBL::Analysis::Runnable::lincRNAFinder;

use strict;  
use warnings;
use vars   qw(@ISA);

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils qw(transfer_supporting_evidence Exon_info);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(Transcript_info print_Transcript_and_Exons);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(compute_6frame_translations) ; 
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils; 
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(Gene_info);
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils; 
use Bio::EnsEMBL::Registry; 

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);



sub new{
  my ($class,@args) = @_;
  
  my $self = $class->SUPER::new(@args); 
  return $self;
} 


sub run{
  my ($self) = @_; 

  # sorting genes in sets 1 + 2 according to biotype + config file settings + exon count 
  my ($single_exon_cdna , $multi_exon_cdna_genes ) = $self->filter_cdna_genes_by_exon_count() ; 
   
  my %types_hash = %{ make_types_hash(\@{$multi_exon_cdna_genes},\@{$self->set_2_prot_genes}, 'SET_1_CDNA','SET_2_PROT')} ; 

  # cluster_Genes returns 2 types of clusters / 4 types of data :
  # 
  #  1) a genomic region with more than one gene [ a set ] which can contain : 
  #  - both sets of genes [ set1 + set 2 ] , a 'twoway-cluster' 
  #  - only one set of genes from [ set1 ] OR [ set2 ], a 'oneway-cluster'  
  # 
  #  2) unclustered genes : a cluster with only one gene, which does not cluster with anything 
  #    - not even with other genes of its own type. 
  #     unclustered 'clusters' only contain ONE gene.   
  
  ### NOTE ###

  ### There are four major stages in lincRNAFinder, the first three are all clustering steps.

  ### The models to be considered as lincRNAs usually are cDNAs.  They can also be models built from Illumina RNASeq data.
  ### For simplicity, we refer to the source models as "cDNAs".

  ### The first two stages aim at reducing the search space of both the cDNAs and the H3 features.
  ### The actual clustering of cDNAs and H3 features occurs in stage three.

  # Stage 1 --- separate cDNAs into those which overlap with protein_coding genes and those which don't
  # -------
  # Cluster cDNAs with protein_coding genes. For those cDNAs which overlap with protein_coding genes, 
  # they will be used in stage 2 to filter H3 features which are not relevant to the search of lincRNAs.  
  # For those which do not overlap with protein_coding genes (denoted "good" cDNAs), they will be taken 
  # to stage 3 where they will be clustered against filtered H3 features from stage 2.
 
  print "Stage 1a) cluster MULTI-exon cDNAs vs protein_coding genes at all exons ";
  print "on the same strand...\n" if (!$self->ignore_strand);
  print "on either strand (i.e. strandedness ignored!)...\n" if ($self->ignore_strand);

  my ($step1_clusters, $step1_unclustered) = cluster_Genes( [@{$multi_exon_cdna_genes}, @{$self->set_2_prot_genes}] , \%types_hash , 0 , $self->ignore_strand ) ;  

  # I select for cDNAs which cluster with protein_coding genes - not because they are lincRNA candidates, but because
  # they are useful for filtering H3 features in stage 2:
  my @cdna_gene_clusters_with_pc =  @{ get_twoway_clustering_genes_of_set($step1_clusters,"SET_1_CDNA") } ;
 
  my ($single_exon_clustered, $single_exon_unclustered);
  if ($self->find_single_exon_candidates) {
    print "Stage 1b) cluster SINGLE-exon cDNAs vs protein_coding genes (strandedness of models always IGNORED) ...\n"; 
    ($single_exon_clustered, $single_exon_unclustered) = cluster_Genes([@$single_exon_cdna, @{$self->set_2_prot_genes}], \%types_hash , 0 , 1 ) ;   # 4th arg = ignore_strand
    push @cdna_gene_clusters_with_pc ,@{ get_twoway_clustering_genes_of_set($single_exon_clustered,"SET_1_CDNA") } ;
  }

  $self->update_and_copy_cdna(\@cdna_gene_clusters_with_pc,"cdna_update_protein_coding" );   # cDNAA DEBUG


  # Stage 2
  # -------
  # Cluster all H3 features against cDNAs in cdna_update_protein_coding clusters (from stage 1). The aim is to
  # remove H3 features which overlap with cDNAs and/or protein_coding genes, i.e. remove features which are irrelevant
  # to our search of lincRNAs, reducing the search space for stage 3.

  print "\nStage 2) cluster all H3 features vs genes in cdna_update_protein_coding clusters to get H3 features which do NOT overlap with the cDNAs...\n"; 
  
  my @efg_sfg = @{ $self->efg_simple_feature_genes};   
  my ($typref,$genes ) = @{make_types_hash_with_genes(\@cdna_gene_clusters_with_pc,\@efg_sfg,'cdna_update_protein_coding','efg')} ;  

  my ($step2_clusters, $step2_unclustered) = cluster_Genes($genes, $typref) ; 

  my @unclustered_efg_genes  = ( 
                                @{get_oneway_clustering_genes_of_set($step2_clusters,"efg")},    # "stacks" of H3 features in a region without cdna_update_protein_coding clusters
                                @{get_oneway_clustering_genes_of_set($step2_unclustered,"efg")}  # "standalone" H3 features in a region without cdna_update_protein_coding clusters
                               );  
  

  # Stage 3
  # -------
  # a) Take the "good" cDNAs from stage 1 which did not overlap with protein_coding genes, and cluster them against the filtered set
  # of H3 features from stage 2 to identify the final set of potential lincRNA genes.
  #
  # b) After the clustering, check with the cDNA overlaps with both H3K4me3 and H3K36me3 features.

  print "\nStage 3) cluster 'good' cDNAs from stage 1 (which did not overlap with protein_coding genes) vs filtered H3 features (from stage 2)...\n\n";     

  # We first separate the H3 feature "genes" into 2 sets, according to their logic name ( H3K4 and H3K36 ), so 
  # we can cluster them separately.
  
  my %logicname_2_efgfeat = %{ separate_efg_features_by_logic_name ( \@unclustered_efg_genes ) } ; 

  # get multi-and single exon "good" cDNAs ( which did not cluster with protein_coding genes in stage 1) 

  my @unclustered_cdna_genes = ( 
                                 @{get_oneway_clustering_genes_of_set($step1_clusters,"SET_1_CDNA")},  
                                 @{get_oneway_clustering_genes_of_set($step1_unclustered,"SET_1_CDNA")},
                               );

  if ($self->find_single_exon_candidates) {   
    push(@unclustered_cdna_genes, @{get_oneway_clustering_genes_of_set($single_exon_clustered,"SET_1_CDNA")} );
    push(@unclustered_cdna_genes, @{get_oneway_clustering_genes_of_set($single_exon_unclustered,"SET_1_CDNA")} );
  }

  my @redun_cdnas_cluster_with_K4_and_or_K36;
  my %check_which_feature_cdna_clusters_with;

  foreach my $lg ( keys %logicname_2_efgfeat ) {
  
    print "Stage 3a) cluster 'good' cDNAs from stage 1 vs filtered $lg features from stage 2...\n" ;   
      
    my @unclustered_efg = @{ $logicname_2_efgfeat{$lg} } ;   

    print "  Have " . @unclustered_efg . " filtered $lg features which don't cluster with genes in cdna_update_protein_coding clusters\n" ; 
    print "  Have " . @unclustered_cdna_genes. " 'good' cDNAs which don't cluster with protein_coding genes \n" ;  
 
    # In the types_hash below:
    # "unclust" cdna_update refers to cDNAs "not in twoway cluster" in stage 1;
    # "unclust" efg refers to features "not in twoway cluster" in stage 2;

    my ($typref_step3,$genes_step3 ) = 
      @{make_types_hash_with_genes(\@unclustered_cdna_genes,\@unclustered_efg,'unclust_cdna_update','unclust_efg')}; 
  
    my ($step3_clusters, $step3_unclustered) = cluster_Genes($genes_step3, $typref_step3 ) ;   

    # Before we identify the "good" cDNAs which overlap with H3 features (and hence will become the
    # final set of potential lincRNAs), two debugging features have been built in here. (More to follow later!)
 
    # Debugging feature (1):
    # For H3 feature "genes"  which don't cluster with any cDNAs and/or protein_coding genes, we count how many there are,
    # assign them a new analysis logic name as defined by "DEBUG_LG_EFG_UNCLUSTERED" in the lincRNAFinder config file, and 
    # write them to the debug database (as specified by the DEBUG_OUTPUT_DB option in the config).

    my @unclust_efg = ( 
                         @{get_oneway_clustering_genes_of_set($step3_clusters,"unclust_efg")},  
                         @{get_oneway_clustering_genes_of_set($step3_unclustered,"unclust_efg")}
                          );    

    print "  Got " . @unclust_efg . " $lg  features which don't cluster with any cDNAs or protein_coding genes, not even cluster with 'good' cDNAs from stage 1.\n" ;       
  
    $self->update_efg_feature_genes ( \@unclust_efg, $self->unclustered_efg_analysis ) ;   # H3 FEAT DEBUG

    # Debugging feature (2):
    # For H3 feature "genes" which do cluster with "good" cDNAs, we may also want to flag them with a new analysis logic_name
    # as specified by "DEBUG_LG_EFG_CLUSTERING_WITH_CDNA" and write them to DEBUG_OUTPUT_DB.
  
    my $efg_clustering_with_cdna_but_not_prot_cod = get_twoway_clustering_genes_of_set($step3_clusters,"unclust_efg"); 
    $self->update_efg_feature_genes ( $efg_clustering_with_cdna_but_not_prot_cod ,$self->efg_clustering_with_cdna_analysis ) ;   # H3 FEAT DEBUG
  
    # Preliminary set of cDNAs which overlap with H3 features:
    my @cdnas_clustering_with_efg_only = @{get_twoway_clustering_genes_of_set($step3_clusters,"unclust_cdna_update")};
    push (@redun_cdnas_cluster_with_K4_and_or_K36,  @cdnas_clustering_with_efg_only);
      
    # Now keeping track of the cDNA which overlapped with a certain type of feature for filtering later if required:
    foreach my $cdna(@cdnas_clustering_with_efg_only) { 
      if ($lg =~ /H3K4/) {
        $check_which_feature_cdna_clusters_with{"H3K4"}{$cdna->dbID} = 1;
      } elsif ($lg =~ /H3K36/) {
        $check_which_feature_cdna_clusters_with{"H3K36"}{$cdna->dbID} = 1;
      }
    }

    if ( ($lg =~ /H3K36/) &&  ($self->check_cdna_overlap_with_multiple_K36) ) {
      my @twoway_clusters = @{get_twoway_clusters($step3_clusters)};
      foreach my $twoway_cluster (@twoway_clusters) {
        my @K36_genes = @{$twoway_cluster->get_Genes_by_Set('unclust_efg')};
        my @cdna_genes = @{$twoway_cluster->get_Genes_by_Set('unclust_cdna_update')};
        # if there are more than one H3K36me3 feature in the cluster, check (by running clustering) to find out whether 
        # the features all pile up in one location or are distributed across a genomic region. We want "multi K36s" which
        # are distributed.
        if (scalar(@K36_genes) > 1) { 
          print "\n  There are more than one K36 feature overlapping with cDNA(s).  Now checking if cDNA(s) overlap with multiple K36 features distributed across a genomic region...\n";
          my (@passed_cdnas_dbIDs) = $self->score_cdna_overlap_with_multi_K36(\@K36_genes, \@cdna_genes);
          foreach my $cdna_dbID(@passed_cdnas_dbIDs) {
            $check_which_feature_cdna_clusters_with{"H3K36_multi"}{$cdna_dbID} ++;
          }
        }
      }
    } # end if $lg =~ /H3K36/ && $self->check_cdna_overlap_with_multiple_K36
  }  # end foreach my $lg...
  
  print "\nStage 3b) count how many cDNAs overlap with only one type or with both types of H3 features...\n" ;

  my @cdnas_overlapping_with_both_K4_and_K36;
  my @cdnas_overlapping_with_only_K4;
  my @cdnas_overlapping_with_only_K36;
  my @cdnas_overlapping_with_multiple_K36;


  my %non_redun_cdnas_clustering_with_efg_only;
  foreach my $cdna (@redun_cdnas_cluster_with_K4_and_or_K36) {
    $non_redun_cdnas_clustering_with_efg_only{$cdna->dbID} = $cdna;
  }

  foreach my $cdna_dbID(keys%non_redun_cdnas_clustering_with_efg_only) {
    my $cdna_to_check =  $non_redun_cdnas_clustering_with_efg_only{$cdna_dbID};
    if ( ($check_which_feature_cdna_clusters_with{"H3K4"}{$cdna_to_check->dbID}) && 
         ($check_which_feature_cdna_clusters_with{"H3K36"}{$cdna_to_check->dbID}) ) {
      push(@cdnas_overlapping_with_both_K4_and_K36, $cdna_to_check);
    } elsif ( ($check_which_feature_cdna_clusters_with{"H3K4"}{$cdna_to_check->dbID}) &&
         (!$check_which_feature_cdna_clusters_with{"H3K36"}{$cdna_to_check->dbID}) ) {
      push(@cdnas_overlapping_with_only_K4, $cdna_to_check);
    } elsif  ( (!$check_which_feature_cdna_clusters_with{"H3K4"}{$cdna_to_check->dbID}) &&
        ($check_which_feature_cdna_clusters_with{"H3K36_multi"}{$cdna_to_check->dbID}) ) {
      push(@cdnas_overlapping_with_multiple_K36, $cdna_to_check);
    } elsif  ( (!$check_which_feature_cdna_clusters_with{"H3K4"}{$cdna_to_check->dbID}) &&
         ($check_which_feature_cdna_clusters_with{"H3K36"}{$cdna_to_check->dbID}) ) {
      push(@cdnas_overlapping_with_only_K36, $cdna_to_check);
    } else {
      throw("cDNA ". $cdna_to_check->dbID . " overlaps with neither H3K4 nor H3K36 features? Something is wrong!");
    }
  }

  print "\nThere are altogether ". scalar(keys%non_redun_cdnas_clustering_with_efg_only)." cDNAs clustering with H3K4me3 and/or H3K36me3.\n";
  print "  ".scalar(@cdnas_overlapping_with_only_K4)." cDNAs overlap with only H3K4me3 features.\n";
  print "  ".scalar(@cdnas_overlapping_with_only_K36)." cDNAs overlap with only H3K36me3 features.\n";
  print "  ".scalar(@cdnas_overlapping_with_multiple_K36)." cDNAs overlap with multiple H3K36me3 features but not H3K4me3 features.\n";
  print "  ".scalar(@cdnas_overlapping_with_both_K4_and_K36)." cDNAs overlap with both H3K4me3 and H3K36me3 features.\n";

  if ($self->check_cdna_overlap_with_both_K4_K36) {
    if ($self->check_cdna_overlap_with_multiple_K36) {
      print "\ncDNAs overlapping with multiple H3K36me3 (but no H3K4me3) OR overlapping with H3K4me3 + H3K36me3 will be considered as lincRNA candidates.\n\n";
      $self->update_and_copy_cdna(\@cdnas_overlapping_with_multiple_K36, 'cdna_multi_K36');  # cDNA DEBUG
      $self->update_and_copy_cdna(\@cdnas_overlapping_with_both_K4_and_K36, 'cdna_both_K4_K36');  # cDNA DEBUG
      my @tmp_result_set = (@cdnas_overlapping_with_multiple_K36, @cdnas_overlapping_with_both_K4_and_K36);
      $self->result_set(\@tmp_result_set);
    } else {
      print "\nOnly cDNAs overlapping with both H3K4me3 and H3K36me3 features will be considered as lincRNA candidates.\n\n";
      $self->update_and_copy_cdna(\@cdnas_overlapping_with_both_K4_and_K36, 'cdna_both_K4_K36');  # cDNA DEBUG
      $self->result_set(\@cdnas_overlapping_with_both_K4_and_K36);
    }
  } else {
    print "\ncDNAs overlapping with H3K4me3 and/or H3K36me3 features will be considered as lincRNA candidates.\n\n";
    my @cdnas_overlapping_some_H3 = values %non_redun_cdnas_clustering_with_efg_only;
    $self->update_and_copy_cdna(\@cdnas_overlapping_some_H3, 'cdna_overlap_H3_not_pc');  # cDNA DEBUG
    $self->result_set(\@cdnas_overlapping_some_H3);
  }                                                       
 
  print "Stage 4) check for 6-frame translations in lincRNA candidates...\n";

  my @genes_with_translations ;  
  RG: for my $rg( @{$self->result_set}  ) {   
    my $new_gene = compute_6frame_translations($rg);
    if (!defined $new_gene->get_all_Transcripts) {
      print "  Could not compute translation for cDNA: gene dbID ". $rg->dbID . " " . $rg->seq_region_name . " " . 
             $rg->seq_region_start . " " . $rg->seq_region_end ."\n";
      next RG;
    }     
    push @genes_with_translations, $new_gene ; 
    # print scalar(@{ $new_gene->get_all_Transcripts} ) ." translations found for gene \n";  
  }
  print "  ".scalar(@genes_with_translations) . " genes with 6-frame-translations found\n" ;     
  # cap the number of transcripts per gene according to config 
  my $capped_genes = $self->cap_number_of_translations_per_gene(\@genes_with_translations) ;  
  my ($short, $long ) = $self->filter_genes_with_long_translations($capped_genes); 
  
  print "\n  ".scalar(@$long)." genes have long translations, discarded from lincRNA pipeline.\n";
  print "  ".scalar(@$short)." genes with short translations found. They will be our lincRNAFinder final output.\n\n";
  
  $self->update_and_copy_cdna($long, 'has_long_translation');   # cDNA DEBUG
  $self->output($short) ; 
}    


sub cap_number_of_translations_per_gene { 
  my ( $self, $genes_with_6f_translations ) = @_ ;   

    my @capped_longest_genes;  

    GENES: for my $g ( @$genes_with_6f_translations ) {    
      my %longest_translations ; 
      for my $t ( @{$g->get_all_Transcripts } ) {  
        my $tl_length = $t->translate->length * 3 ; 
        push @{ $longest_translations{$tl_length}}, $t;
      } 

      my @tl_length = sort { $b <=> $a } keys %longest_translations ;   
 
      # We only take the first n longest translations of a gene.
      # "n" could have been defined by "MAXIMUM_TRANSLATION_LENGTH_RATIO"
      # in the config, or be set by the size of the tl_length array
      # (i.e. all translations are taken, no limits set!)

      my $max_translations_stored_per_gene; 
      if ( defined $self->max_translations_stored_per_gene ) {  
        $max_translations_stored_per_gene = $self->max_translations_stored_per_gene ; 
      } else { 
        $max_translations_stored_per_gene = scalar(@tl_length);   
      } 
      @tl_length = splice @tl_length,  0, $max_translations_stored_per_gene ;  
      my $mgt = new Bio::EnsEMBL::Gene(); 
      $mgt->biotype($g->biotype);  
      for my $length  ( @tl_length ) { 
         for my $lt ( @{ $longest_translations{$length} } ) {  
           $mgt->add_Transcript($lt);  
         }
      }
     push @capped_longest_genes , $mgt ;   
   }   
  return \@capped_longest_genes; 
} 


sub filter_genes_with_long_translations {  
  my ( $self, $capped_longest_genes ) = @_ ; 

  my $max_trans_length_ratio;
   if ( defined $self->maximum_translation_length_ratio ){  
     if ( $self->maximum_translation_length_ratio < 101 ) {  
       $max_trans_length_ratio =  $self->maximum_translation_length_ratio; 
     } else {  
       throw("translation-length-to-transcript length ratio > 100 does not make sense.\n"); 
     } 
   }else { 
     $max_trans_length_ratio = 100 ; 
   }   

  # if a gene has a transcript with a translation over a specified % threshold transcript
  # length, the gene will not be taken as lincRNA candidate.
  
  # if WRITE_DEBUG_OPTION is turned on, the cDNA "genes" with long translations will
  # be written to the DEBUG_OUTPUT_DB with new biotype "long translation".
  
  my (@long_genes, @short_genes );  

  GENES: for my $g ( @$capped_longest_genes ) { 
    for my $t ( @{$g->get_all_Transcripts } ) {   

       my $tl_length = $t->translate->length * 3 ; 
       my $tr_length = $t->length ;    
       my $ratio = sprintf('%.1f',$tl_length*100 /  $tr_length);

       if ( $ratio > $max_trans_length_ratio ) { 
         # The cDNA "genes" have no usable display IDs, hence not printing any in the next line:
         print "  LONG_TRANSLATION: cDNA has translation-length to transcript-length ratio ($ratio) higher than max. value allowed in config ($max_trans_length_ratio).\n";
         my $analysis =Bio::EnsEMBL::Analysis->new(
                                       -logic_name => 'long_translation_cDNA' ,
                                       );
         $g->analysis($analysis);
         $g->biotype('long_translation_cDNA');
         push @long_genes, $g; 
         next GENES; 
       }  
    }     
    push @short_genes , $g ;
  }
  return ( \@short_genes, \@long_genes );  
} 



sub separate_efg_features_by_logic_name {  
  my ( $efg_genes ) = @_ ; 

  my %logic_name_2_efg_features ;  
  for my $e ( @$efg_genes ) {
      $e->biotype($e->analysis->logic_name);     
      push @{$logic_name_2_efg_features{$e->analysis->logic_name}}, $e ;  
  }    
  return \%logic_name_2_efg_features ; 
}



=head2 result_set

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : arrayref of result_set
  Function  : push array passed into onto result_set array
  Returntype: arrayref
  Exceptions: throws if not passed correct type
  Example   : $self->result_set(\@output);

=cut


sub result_set{
  my ($self, $result_set) = @_;
  if(!$self->{'result_set'}){
    $self->{'result_set'} = [];
  }
  if($result_set){
    if(ref($result_set) ne 'ARRAY'){
      throw('Must pass RunnableDB:result_set an array ref not a '.$result_set);
    }
    push(@{$self->{'result_set'}}, @$result_set);
  }
  return $self->{'result_set'};
}

sub update_efg_feature_genes {  
  my ($self,  $aref , $analysis ) = @_;  

  throw(" no analysis " ) unless $analysis ; 
  my @updated_gene; 
  for my $gene ( @$aref ) {     
    $gene->dbID(undef); 
    $gene->analysis($analysis);   
    push @updated_gene, $gene ; 
  }  
  $self->updated_efg_genes(\@updated_gene) ; 
} 
 
sub updated_efg_genes { 
  my ( $self,$arg ) = @_ ;    

  unless (defined $self->{updated_efg_genes} ) {
    $self->{updated_efg_genes} = [] ; 
  }

  if ( defined  $arg ) { 
   push @{$self->{updated_efg_genes}} , @$arg;  
  }  
  return $self->{updated_efg_genes} ; 
} 

sub updated_cdnas {  
  my ( $self,$arg ) = @_ ;  
  
  if (! $self->{cdna_updated} ) {  
    $self->{cdna_updated} = []; 
  } 

  if (defined $arg ) { 
   push @{$self->{cdna_updated}} , $arg;  
  }  
  return $self->{cdna_updated} ; 
} 

sub update_and_copy_cdna {   # This runnable method is only called by the runnableDB when DEBUG option is turned on.
  my ($self, $array, $new_bt) = @_ ;  

  for my $g ( @$array) { 
    $g->biotype($new_bt);  
    $self->updated_cdnas($g);
  }  
} 


sub filter_cdna_genes_by_exon_count {   
  my ( $self ) = @_ ; 

  my (@single_exon, @multi_exon )  ;  

  for my $cg ( @{ $self->set_1_cdna_genes} ) {     
    for my $ct ( @{ $cg->get_all_Transcripts } ) {  
       my @ce = @{$ct->get_all_Exons };  

       if ( scalar(@ce == 1 )) {  
         push @single_exon, $cg ; 
       } else { 
         push @multi_exon, $cg ; 
       } 
    }
  }   
  return ( \@single_exon, \@multi_exon );
} 


=head2 score_cdna_overlap_with_multi_K36

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::lincRNAFinder
  Arg [2]   : Reference to an array of H3K36me3 gene objects
  Arg [3]   : Reference to an array of cDNA gene objects
  Function  : Selects for cDNA gene objects which overlap with multiple
              H3K36me3 features distributed across a genomic region.
              cDNA genes overlapping with multiple H3K36me3 features
              piling up in the same genomic region will not be 
              selected.
  Returntype: An array dbIDs (integers) of cDNAs which have been
              selected.
  Example   : $self->score_cdna_overlap_with_multi_K36(\@K36_genes, \@cdna_genes);

=cut

sub score_cdna_overlap_with_multi_K36 {
  my ($self, $K36_genes, $cdna_genes) = @_;
  my @dummy_set2_models;   # blank
  my ($clustered, $unclustered) = @{simple_cluster_Genes($K36_genes, "K36" , \@dummy_set2_models, "dummy")};
  my @K36_regions = (@$clustered, @$unclustered);
  my $K36_cluster_cnt = scalar(@K36_regions);
  print "    Got a total of $K36_cluster_cnt clusters after local clustering, ". scalar@$clustered . " of which are K36 clusters and " . scalar@$unclustered . " are K36 singletons.\n";

  my @cdna_overlapping_multi_K36;
  my %scoring;
  
  if ($K36_cluster_cnt > 1) { # i.e. the H3K36me3 feats don't fall into the same genomic location but distributed across exons
    foreach my $cdna(@$cdna_genes) {
      my $curr_sr_start = $self->query->start;
      my $curr_sr_end = $self->query->end;
      foreach my $K36_region(@K36_regions) {
        my $K36_clust_start = $curr_sr_start + $K36_region->start;
        my $K36_clust_end = $curr_sr_start + $K36_region->end;
        print "    NEW LOCAL CLUSTER REGION: " . $self->query->seq_region_name . ",  $K36_clust_start - $K36_clust_end.\n";
       if ($cdna->seq_region_end >= $K36_clust_start && $cdna->seq_region_start <= $K36_clust_end) {   # overlap
          print "    FOUND MATCH in this region for cDNA: " .$cdna->seq_region_start . " -  " . $cdna->seq_region_end . " (dbID " . $cdna->dbID.")\n";
          $scoring{$cdna->dbID} ++;
        }
      }
    }
  } else {
    print "    This means only one block of H3K36me3 features are overlapping with cDNA(s).\n"
  }
  return keys%scoring;
}

use vars '$AUTOLOAD';
sub AUTOLOAD {
 my ($self,$val) = @_;
 (my $routine_name=$AUTOLOAD)=~s/.*:://; #trim package name
 if ( defined $val ) { 
   $self->{$routine_name}=$val; 
 }
 return $self->{$routine_name} ;
}
sub DESTROY {} # required due to AUTOLOAD

1;  
