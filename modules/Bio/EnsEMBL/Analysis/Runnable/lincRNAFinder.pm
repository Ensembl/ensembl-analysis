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
use Bio::EnsEMBL::Analysis::Runnable::Blast; 
use Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Pfam; 
use Bio::EnsEMBL::Analysis::Tools::FilterBPlite;   
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry; 

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);



sub new{
  my ($class,@args) = @_;
  
  my $self = $class->SUPER::new(@args); 
  return $self;
} 


sub run{
  my ($self) = @_; 

  # sorting genes in sets 1 + 2 according to biotype + config fle + exon count 
  my ($single_exon_cdna , $multi_exon_cdna_genes ) = $self->filter_cdna_genes_by_exon_count() ; 
   
  my %types_hash = %{ make_types_hash(\@{$multi_exon_cdna_genes},\@{$self->set_2_prot_genes}, 'SET_1_CDNA','SET_2_PROT')} ; 

  # cluster_Genes returns 2 types of clusters / 4 types of data :
  # 
  #  1) a genonmic region with more than one gene [ a set ] which can contain : 
  #  - both sets of genes [ set1 + set 2 ] , a 'twoway-cluster' 
  #  - only one set of genes [ set1 ] OR set2, a 'oneway-cluster'  
  # 
  #  2) unclustered genes : a cluster with only one gene, which does not clustere with anything 
  #    - not even with other genees of it's own type. 
  #     unclustered 'clusters' only contain ONE gene.   
  
  print "1a) 1st clustering of cdna_update-multi-exon-genes  vs protein_coding\n"; 

  my ($step1_clusters, $step1_unclustered) = cluster_Genes( [@{$multi_exon_cdna_genes},@{$self->set_2_prot_genes}] , \%types_hash ) ;  


  # I am lookign for cdna's which cluster with protein_coding genes - I only want the cdna's  
  # find genes which cluster with protein_coding and cdna  
  
  my @cdna_gene_clusters_with_pc =  @{ get_twoway_clustering_genes_of_set($step1_clusters,"SET_1_CDNA") } ;
 
 # cluster single_exon + protein-coding WITHOUT strand information 
 
  print "1b) 1st clustering of cdna_update-single-exon-genes  vs protein_coding\n"; 

  my ($single_exon_clustered, $unclustered) = cluster_Genes([@$single_exon_cdna, @{$self->set_2_prot_genes}], \%types_hash ) ;   
  push @cdna_gene_clusters_with_pc ,@{ get_twoway_clustering_genes_of_set($single_exon_clustered,"SET_1_CDNA") } ;


  $self->update_and_copy_cdna(\@cdna_gene_clusters_with_pc,"cdna_update_protein_coding" );   

  #
  # now cluster the cdnas [ the ones which cluster with protein coding genes from 1st clustering ] with efg features 
  #  
  # cdnas_genes  = [ where  cDNA and a protein_coding gene form a cluster together .... ]  
  # 
  
  my @efg_sfg = @{ $self->efg_simple_feature_genes};   
  my ($typref,$genes ) = @{make_types_hash_with_genes(\@cdna_gene_clusters_with_pc,\@efg_sfg,'cdna_update_protein_coding','efg')} ;  


  print "2nd clustering: cluster all EFG with cdna_update_protein_coding cluster to get oneway-EFG which do not cluster\n"; 
  # clustering EFG with cdna_update genes [only the cdna_update genes which cluster with protein_coding genes ] 
  my ($step2_clusters, $step2_unclustered) = cluster_Genes($genes, $typref) ; 

  my @unclustered_efg_genes  = ( 
                                @{get_oneway_clustering_genes_of_set($step2_clusters,"efg")},  
                                @{get_oneway_clustering_genes_of_set($step2_unclustered,"efg")}
                               );  
  # @unclustered_efg_genes = ONEWAY efg_domains in regions where they don't overlap cdna ( and the cdna don't overlap protein_coding )  
  
  
  print "3rd clustering - cluster unclustered cdna from 1st clustering with unclustered EFG features from 2nd clustering.. \n";     

  # We now separate out the EFG simple feature genes into 2 sets, according to their logic name ( H3K4 and H3K36 )
  # this is to cluster them separately 
  
  my %logicname_2_efgfeat = %{ separate_efg_features_by_logic_name ( \@unclustered_efg_genes ) } ; 

  # get mult-and single exon cDNA unclustered ( does not cluster with protein_coding ) 
  my @unclustered_cdna_genes = ( 
                                 @{get_oneway_clustering_genes_of_set($step1_clusters,"SET_1_CDNA")},  
                                 @{get_oneway_clustering_genes_of_set($step1_unclustered,"SET_1_CDNA")},
                                 @{get_oneway_clustering_genes_of_set($single_exon_clustered,"SET_1_CDNA")},
                                 @{get_oneway_clustering_genes_of_set($unclustered,"SET_1_CDNA")}
                               );   

  my %step_4_sets;
 
  for my $lg ( keys %logicname_2_efgfeat ) {
      print "clustering unclustered efg [ $lg ]  with step_1 unclustered cdna [cdna which does not cluster with protein_coding....\n" ;   
      # unclustered efg = efg which don't cluster with cdna_update_protein_coding (cdna genes which cluster with protein_coding genes
      # unclustered cdna = cdna which does not cluster with protein_coding 
      
      my @unclustered_efg = @{ $logicname_2_efgfeat{$lg} } ;   

      print "have " . @unclustered_efg . " EFG FEAT type $lg which don't cluster with cdna_update_protein_coding \n" ; 
      print "have " . @unclustered_cdna_genes. " CDNA which dont cluster with protein_coding \n" ;  
 
      my ($typref_step3,$genes_step3 ) = 
       @{make_types_hash_with_genes(\@unclustered_cdna_genes,\@unclustered_efg,'unclust_cdna_update','unclust_efg')}; 
  
      my ($step3_clusters, $step3_unclustered) = cluster_Genes($genes_step3, $typref_step3 ) ;   
  
      # just for counting efg which don't cluster with any cdna or pc : 
      my @unclust_efg = ( 
                         @{get_oneway_clustering_genes_of_set($step3_clusters,"unclust_efg")},  
                         @{get_oneway_clustering_genes_of_set($step3_unclustered,"unclust_efg")}
                        );    
      print "got " . @unclust_efg . " unclustered efg for $lg ( efg which don't cluster with uncl. cdna ) \n" ;       
  
     # change oneway-efg-clusters ( oneway-only one set is included )  
      $self->update_efg_feature_genes ( \@unclust_efg, $self->unclustered_efg_analysis ) ; 
  
     # twoway-cluster : a cluster which contains genes of set A and B 
     # oneway-cluser  : a cluster which contains only one set either A OR B 
  
     # unclustered EFG which do not cluster with protein_coding but with cdna : 
  
     # efg-features which cluster two-way with cdna ( efg_cdna_cluster )  -  efg_cdna_cluster_NO_protein 
     my $efg_clustering_with_cdna_but_not_prot_cod = get_twoway_clustering_genes_of_set($step3_clusters,"unclust_efg"); 
     $self->update_efg_feature_genes ( $efg_clustering_with_cdna_but_not_prot_cod ,$self->efg_clustering_with_cdna_analysis ) ; 
  
     # finally the result :  a list of cDNAs which cluster two-way with efg only 
     my $cdna_only_clustering_with_efg_only = get_twoway_clustering_genes_of_set($step3_clusters,"unclust_cdna_update"); 
     $self->update_and_copy_cdna($cdna_only_clustering_with_efg_only, 'cdna_efg_no_pc'); 
   
     $self->result_set($cdna_only_clustering_with_efg_only);    
  }     

  #
  # Step 4 - identify the genes which cluster with H3K4 as well as with H3K36 domains.   
  #  
  #my %cdnas_clustering_with_H3_and_H4; 
  #for my $cdna ( @{$self->result_set} ) {   
  #   $cdnas_clustering_with_H3_and_H4{$cdna}{num}++; 
  #   $cdnas_clustering_with_H3_and_H4{$cdna}{cdna}=$cdna;
  #}   
  # This section needs review - to identify which cDNA clusters with H3K4 and H3K36 you have to 
  # make sure you don't count cDNAs which cluster with multiple H3K4; 
  #for my $c (keys %cdnas_clustering_with_H3_and_H4 ) {  
  #  if ( $cdnas_clustering_with_H3_and_H4{$c}{num} > 1 ) {    
  #     print "clusters with H3K4 and H3K36 domain\n" ;  
  #  }  
  #}  
 
 
  my @genes_with_translations ;  
  for my $rg( @{$self->result_set}  ) {   
    push @genes_with_translations, compute_6frame_translations($rg); 
  }
  print scalar(@genes_with_translations) . " genes with 6-frame-translations found\n" ;    

  $self->output(\@genes_with_translations); 
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

#sub set_1_cdna_genes{
#  my ($self, $arg) = @_;
#  if($arg){
#    for my $g  ( @$arg ) {  
#      $g->biotype($g->biotype."_set1") ;    
#    } 
#    $self->{'set_1_genes'} = $arg; 
#  }
#  return $self->{'set_1_genes'};
#}
#
#
#sub set_2_prot_genes{
#  my ($self, $arg) = @_;
#  if($arg){
#    for my $g  ( @$arg ) { 
#      $g->biotype($g->biotype."_set2") ;   
#    } 
#    $self->{'set_2_genes'} = $arg;
#  }
#  return $self->{'set_2_genes'};
#}  
# 

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

sub update_and_copy_cdna { 
  my ($self, $array, $new_bt  ) = @_ ;  

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
