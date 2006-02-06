package Bio::EnsEMBL::Analysis::Runnable::TranscriptCoalescer;

use strict;
use warnings; 

use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Bio::Tools::CodonTable;

#check the cluster-modules 
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneCluster;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonCluster;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::UtilFuncs; 

use Bio::EnsEMBL::Analysis::Config::GeneBuild::Databases;  
use Bio::EnsEMBL::Analysis::Config::GeneBuild::TranscriptCoalescer;  

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw (return_translation) ;

use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended; 
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended; 

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);




sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->{'_modified_genes'} =[] ;        # array ref to modified genes to write to new db
  $self->{'_discarded_transcripts'} = []; # array ref to discarded transcripts
  $self->{'_genes'} = [];                 #array of genes to test;  
    
  my( $all_genes, $new_biotype, $evidence_sets ) = 
      rearrange([qw(
                     ALL_GENES
                     NEW_BIOTYPE
                     EVIDENCE_SETS
                    )], @args);

  $self->{merged_tr}=[] ; 
  #
  # we add one more evidence set called 'est_merged' to be able 
  # to re-cluster the newly constructed genes with the input-genes
  # ( of the sets 'simgw' and 'abinitio' genes ) 
  #
  # my %types_hash{'name_of_the_set'} = [gene_bt_100, simgw_200, simgw300, ...] 
  # types_set{'genes_with_est_source'} = [est_100, est_200, est_300....] 
  #
  $self->{new_biotype} = $new_biotype ;  
  ${$evidence_sets}{'est_merged'}=[$self->{new_biotype}] ; 

  # href to $hash{'evi_set'} = \@{biotype1,obiotype2}  
  $self->{evidence_sets} = $evidence_sets;    

  $self->{write_filtered_transcripts} = $WRITE_FILTERED_TRANSCRIPTS ;  
  $self->{write_alternative_transcripts} = $WRITE_ALTERNATIVE_TRANSCRIPTS  ; 
  $self->{adjudicate_simgw_est} = $ADJUDICATE_SIMGW_EST ; 
  
  $self->{all_genes_href} = $all_genes ;      # hashref $hash{'biotype'} = \@genes 
  $self->{v} = $VERBOSE ; # verbose or not 
  return $self ; 
}


sub run { 
  my ($self) = @_ ; 
  print " all genes fetched\n" if $VERBOSE  ;   
  # get all genes of evidence set 'est' on slice and cluster them
  
  my @allgenes = @{ $self->get_genes_by_evidence_set('est') }  ;
 
  my ($clusters, $non_clusters) = cluster_Genes(\@allgenes, $self->get_all_evidence_sets ) ; 

  push @$clusters, @$non_clusters if (scalar(@$non_clusters) > 0 ) ;   

  my (@all_merged_est_transcripts, @tr_to_delete ) ; 

  #
  # STAGE_1 : 
  # begin the process of merging est-gene and est-gene by using the conserved introns as evidence 
  #
  print "non_clusters : " . scalar(@$non_clusters) . "\n" ; 

  GENE_CLUSTER: foreach my $gene_cluster (@$clusters) {
    if ($self->{v}) { 
      print "\n" ; print_cluster_info($self->query,$gene_cluster ) ; print "="x80; print "\n" ;  
    }

    # cluster transcripts and exons and start the recursions on transcripts 
    # and exon clusters / this is the main routine to extend / merge the est-genes  
    #
    
    print "Starting main_clustering_and_recursion\n" if $self->{v} ;  
   
    my $tr = $self->main_clustering_and_recursion($gene_cluster) ; 


    # remove redundant transcripts and edit last exons

    if ( $tr ) {
      if ($self->{v}){  
        print "Stage 1: Number Genes resulting out of pure est-merging           : ".scalar(@$tr)."\n"; 
        print_transcripts_and_exons ($tr) ; 
      }
      # remove red. trans, edit terminal exons and remove redundant trans again 
      my @non_red_trans = 
      @{ remove_redundant_transcripts( $self->edit_terminal_exons(remove_redundant_transcripts($tr)))}; 

      print "Stage 2: Number Trans. after term exon edit / removing redundancy : "
       . scalar(@non_red_trans)."\n" if $self->{v}; 

      # transcripts are splitted in ones which we're going to cluster again and ones which are 
      # removed by now (biotype of them has been changed) 

      my ($non_overlapping_trans , $removed_trans )=remove_overlapping_transcripts(\@non_red_trans)  ;  

      print "Stage 3: Number non-overlapping Transcripts                       : "
       . scalar(@$non_overlapping_trans) 
       . " (" . scalar(@$removed_trans)  . " removed )\n\n" if $self->{v} ; 

      push @tr_to_delete  , @$removed_trans if $self->{write_filtered_transcripts}  ; 
      #
      # @$non_overlapping_trans are transcript stuctures which have been produced by est-merging
      # these might also be 'pure' ests / if we have a tr-structure later on (by abintio merging)
      # which overlap this transcript, this structure will be removed from this array  
      #
      push @all_merged_est_transcripts ,   @$non_overlapping_trans  ;
      
      #print "non-overlapping transcripts which have been made by est-merging :\n" ; 
      #print_transcripts_and_exons($non_overlapping_trans) ; 
    }
  } # GENE_CLUSTER 


 ##
 ##
 ##                 RE-clustering of genes with other evidence 
 ##                 to extend genes by abintio / simgw 
 ##
 ##

  if ($VERBOSE) { 
    print "after the first round\n" ; 
    print "="x80 ; print "\n" ; 
    print_transcripts_and_exons(\@all_merged_est_transcripts) ; 
  }
  
  # 
  # building new gene set which contains genes of set abinitio and simgw and the 
  # freshly merged genes (biotype of these genes is set in config) 
  #

  my @new_gene_set ; 
  push @new_gene_set , @{$self->convert_to_genes( \@all_merged_est_transcripts ) } ; 
  push @new_gene_set , @{ $self->get_genes_by_evidence_set('abinitio') }  ;
  push @new_gene_set , @{ $self->get_genes_by_evidence_set('simgw') }  ; 
  
  #
  # re-cluster freshly merged / extended of type 'est_merged' genes of other sources 
  #
  my ($clusters2, $non_clusters2) = cluster_Genes(\@new_gene_set, $self->get_all_evidence_sets ) ; 

  GENE_CLUSTERS: for my $gene_cluster (@$clusters2) { 
    #
    # check if there are genes of type 'est_merged' in result 
    #

    #my @est_merged_genes = $gene_cluster->get_Genes_by_Set('est_merged') ; 

    my %sets_in_gene_cluster ; 
    @sets_in_gene_cluster { @{$gene_cluster->get_sets_included} }= ()  ;  

    # If there are est_merged genes in the cluster (est's which have been merged in the first step)
    # try to merge them with abintio or simgw

    if ( exists $sets_in_gene_cluster{est_merged} ) { 
      # try to do recursive merge with simgw / abinitio, begin with clustering (Transcripts & Exons) 
      my @tr_merged = 
       map {@{$_->get_all_Transcripts}} $gene_cluster->get_Genes_by_Set ( 'est_merged' ); 
      my $exon_clusters = get_exon_clustering_from_gene_cluster($gene_cluster) ;    

      # recursive approach (merge already merged genes by abinitio or simgw)
      if ($gene_cluster->strand == 1 ) { 
        @tr_merged = sort { $a->seq_region_start <=> $b->seq_region_start } @tr_merged ;  
      }else { 
        @tr_merged = sort { $b->seq_region_start <=> $a->seq_region_start } @tr_merged ;  
      }
      my $start_transcript = shift @tr_merged ;  
      my $t =  $self->merge_transcripts_recursion($start_transcript , \@tr_merged, $exon_clusters ); 
  
      if ( @{ $self->{merged_tr} }>0 ){
        #print "cluster: have merged_tr genes \n" ; 
        my @cloned_tr = @{ $self->{merged_tr} } ; 

        $self->{merged_tr}=[] ;   # reset counter
        #print "gene_comes_from_merging\n" ;   
        #print_transcripts_and_exons(\@cloned_tr) ; 

        # remove redundant genes which result out of abintio/simgw merged transcripts 
        @cloned_tr = @{ remove_redundant_transcripts( \@cloned_tr ) }  ; 

        # remove overlapping transcripts 
        my ($non_overlapping_trans, $removed_trans )= remove_overlapping_transcripts( \@cloned_tr ) ;
 
        # if we want to write the filtered / removed /overlapping genes as well to DB 
        # i.e. for debugging   (biotype 'del_') (see config)  
        push @tr_to_delete,  @$removed_trans if ($self->{write_filtered_transcripts}); 

        # prevent to add the start-transcript a second to the set of all_merged_est_transcripts 

        push @all_merged_est_transcripts , @$non_overlapping_trans ;  

        @all_merged_est_transcripts = 
         @{ remove_redundant_transcripts( \@all_merged_est_transcripts ) }; 

        ($non_overlapping_trans, $removed_trans ) 
          = remove_overlapping_transcripts( \@all_merged_est_transcripts ) ;     
        @all_merged_est_transcripts =  @$non_overlapping_trans ; # this is ok

        push @tr_to_delete , @$removed_trans if $self->{write_filtered_transcripts}; 
      } else {
        # recover transcripts cause they haven't been merged  
        push @all_merged_est_transcripts , @tr_merged ;
      } 
    } else { 
      #print "have no est_merged genes\n" ; 
    }
  } 

 #
 # use rule-set to decide if we want to take simgw or est gene  
 #
 #######################################################################

 if ($self->{adjudicate_simgw_est}) { 
   print "deciding weather to use simgw or est-combined gene in final set\n" ;  

   #
   # 3rd re-clustering of simgw against est genes        
   #  
   my @est_simgw = @all_merged_est_transcripts ;

   push @est_simgw , @{ $self->get_genes_by_evidence_set('est_merged') }  ;
   push @est_simgw , @{ $self->get_genes_by_evidence_set('simgw') }  ;
   
   print "3rd reclustering to choose if we want to have simgw or est or both\n" ;  
   
   ($clusters, $non_clusters)  = cluster_Genes(\@est_simgw, $self->get_all_evidence_sets ) ; 
   GENE_CLUSTER: foreach my $gene_cluster (@$clusters) {
     # cluster exons and re-set the exon-cluster relation
     my @exon_clusters  = @{ get_exon_clustering_from_gene_cluster( $gene_cluster ) };  
     for my $ec (@exon_clusters) { 
       for my $e ( @{$ec->get_all_Exons_in_ExonCluster} ) {
         $e->cluster($ec) ; 
       }
     } 
     
     # get simgw genes and remove overlapping ones  
     my @tr_simgw = map {@{$_->get_all_Transcripts}} $gene_cluster->get_Genes_by_Set('simgw'); 
     my ($non_ov, $ov) = remove_overlapping_transcripts(\@tr_simgw) ; 
     push @tr_to_delete , @$ov if $self->{write_filtered_transcripts}; 
     @tr_simgw = @$non_ov ; 

     my @tr_est = $gene_cluster->get_Genes_by_Set('est_merged') ;  
     my %sets_in_gene_cluster ;
     @sets_in_gene_cluster{@{ $gene_cluster->get_sets_included }}=() ;   

     # FALL UNTERSCHEIDUNGEN 

     if (exists $sets_in_gene_cluster{simgw} && exists $sets_in_gene_cluster{est_merged} ) { 

       # 
       # filter genes and decide if we want simgw , est_merged or whatever / UTR addtion   
       #  
      print "\n"x4; 
      print "We have SIMGW and EST_MERGED genes in cluster :\n" ; 
      my ($keep, $removed) = compare_simgw_and_est_merged(\@tr_simgw, \@tr_est,$self->{write_alternative_transcripts}) ; 



    }elsif ( exists $sets_in_gene_cluster{est_merged} && !exists $sets_in_gene_cluster{simgw} ) { 

      print "have no simgw and but est_merged genes \n" ; 
      
#        
#    }elsif ( exists $sets_in_gene_cluster{simgw} && !exists $sets_in_gene_cluster{est_merged} ) { 
#      print "have simgw but no est_merged\n" ; 
#      
#      # there are only simgw-genes (and some abinitio ) in the GENE_CLUSTER 
#       
#      my @simgw_trans  =  
#       map { @{ $_->get_all_Transcripts } } $gene_cluster->get_Genes_by_Set('simgw') ;
##      #my $exon_clusters = get_exon_clustering_from_gene_cluster($gene_cluster) ;    
# 
##      for ( @tr_simgw ) { 
##        $_->biotype($self->{new_biotype} . "_" . $_->biotype) ; 
##      } 
##      push @all_merged_est_transcripts, @tr_simgw ; 
    } 
  }  
} 



  print "Having " . scalar( @all_merged_est_transcripts ) . " genes so far for this slice \n\n" ; 
  print_transcripts_and_exons(\@all_merged_est_transcripts) if $self->{v}; 
  
  # adding translations to transcripts 
  my @trans_with_tl = @{add_translation_and_trans_supp_features_to_transcripts(
       \@all_merged_est_transcripts) } ; 
 
  $self->output ( $self->convert_to_genes ( \@trans_with_tl  )) ;   
  return ; 
}



sub add_translation_and_trans_supp_features_to_transcripts  { 
   my ($trans) = @_ ; 
   my @trans_with_tl ; 
   # try to add translation  TRANSLATION
     for my $tr (@$trans) {   
       for my $ex (@{$tr->get_all_Exons} ) { 
           # print join ("\n", @{ $ex->transcript->get_all_supporting_features}  )  ;  
           if ($ex->transcript) { 
             if ($ex->transcript->get_all_supporting_features) { 
               $tr->add_supporting_features ( @{ $ex->transcript->get_all_supporting_features } ) ;  
             }else {  
              throw( "transcript misses transcript_supporting_features\n")  ; 
             }
           } else { 
              throw( "Exon has no transcript assigned :\n")  ; 
           }
        }
        $tr->translation(return_translation($tr)) ; 
        
#        print "=========================\n" ;  
#        print_transcripts_and_exons($tr) ;  
#        print "start_ex_TL: " ; 
#        print_object (  $tl->start_Exon) ;
#        print "  end_ex_TL: " ; 
#        print_object (  $tl->end_Exon) ;
#        print "\n" ;
       
        push @trans_with_tl , $tr; 
     }
   return \@trans_with_tl ;  
}



sub compare_simgw_and_est_merged {
  my ($simgw_tr, $est_tr , $write_alt_trans  ) = @_ ; 
  my (@keep, @removed, @all) ; 
  print scalar(@$simgw_tr) . " simgw genes\n" ; 
  print scalar(@$est_tr) . " EST genes\n" ; 
  print "\n\nEST-transcripts:\n--------------------------------\n" ; 
  print_transcripts_and_exons($est_tr) ; 
  print "\n\nSIMGW-transcripts:\n--------------------------------\n" ; 
  print_transcripts_and_exons($simgw_tr) ; 
  print "\n\n\n\n" ;  


  push @all, @$simgw_tr , @$est_tr ; 
  @all = @{ sort_transcripts_by_length ( \@all) } ; # sort tr descending their genomic extent length
  # 
  # check how supported the simgw transcripts are by the est transcripts
  #
  for my $simgw (@$simgw_tr) { 
    print "simgw $simgw\n" ; 
    my @sim_exons = @{ $simgw->get_all_Exons } ; 
    for my $est (@$est_tr) { 
      print "est $est\n" ; 
      my @est_exons = @{ $est->get_all_Exons } ; 
      for my $se (@sim_exons){ 
         print "simgw-exon  $se\n" ; 
        for my $ee (@est_exons) {
         print "ee-exon  $ee\n" ; 
          
          if ($ee->overlaps($se)) { 
            # ExonExtended
            $se->exon_is_overlapped($ee) ;
           
            # TranscriptExtended
            $simgw->nr_exons_overlapped_by_est(1) ; 
            print "\n\noverlaps... " ; 
            $simgw->different_est_support($est) ; 
            print "........ok\n" ;
            # info of 5prim end of transcript            
            if ( $se->is_5prim_exon ) { 
              $simgw->has_5prim_support($est) ;

              if ($simgw->seq_region_strand eq '1') { 
                if ($simgw->seq_region_start > $est->seq_region_start) { 
                  # est can be used for utr addtion 5prim end 
                  if ($simgw->extend_5prim_end) { 
                    # there was already another est which could be used -> take the longer one
                    if ($simgw->extend_5prim_end->seq_region_start > $est->seq_region_start) { 
                      $simgw->extend_5prim_end($est) ; 
                    }
                  }else { 
                      #no other est known so init with this one 
                      $simgw->extend_5prim_end($est) ; 
                 }     
                }
              } else { # minus strand 
              print "\n\nminus\n\n" ; 
                if ($simgw->seq_region_start > $est->seq_region_start ) {  
                  print "\n\ncheck1\n\n" ; 
                  # est can be used to add UTR
                  if ($simgw->extend_5prim_end) { 
                    if ($simgw->extend_5prim_end->seq_region_start < $est->seq_region_start) { 
                      $simgw->extend_5prim_end($est) ; 
                    }
                  }else { 
                    print "\n\ncheck2\n\n" ; 
                    $simgw->extend_5prim_end($est) ; 
                    print "\n\ncheck2\n\n" ; 
                  }
                }
              } 
            } # end is 5prim exon
                    print "\n\ncheck3\n\n" ; 

            # info of 3prim end of Transcript
            # set the est which can be used to extend transcript 3prim end UTR
            if ( $se->is_3prim_exon ) { 
             print "have 3prim exon : " ;
             print_object($se) ; 
             print_transcripts_and_exons($simgw) ; 
             print "\n\n\n" ; 
              $simgw->has_3prim_support($est) ; 
              if ($simgw->seq_region_strand eq '1') { 
                if ($simgw->seq_region_end < $est->seq_region_end) { 
                  if ($simgw->extend_3prim_end) { 
                    my $before_est = $simgw->extend_3prim_end; 
                    if ($before_est->seq_region_end < $est->seq_region_end ) { 
                      $simgw->extend_3prim_end($est) ; 
                    }                      
                  } else{ 
                  $simgw->extend_3prim_end($est) ; 
                  }
                }
              }else { # minus strand
                if( $est->seq_region_start <= $simgw->seq_region_start) {
                  #  |estestest start |------|estestest|
                  #               |simgw start |------|simgw|
                  if ( $simgw->extend_3prim_end) {
                    if ( $simgw->extend_3prim_end->seq_region_start > $est->seq_region_start) { 
                      # re-set the extend_3prim_end if this one goes further  
                      $simgw->extend_3prim_end($est) ; 
                    }
                  }else { 
                   # no extend_3prim_end defiend till now, let's set it ! 
                   $simgw->extend_3prim_end($est) ; 
                  }
                }else{
                  # est can't be used to extend simgw  
                }
              # DEBUG 
              }
            } # end 3prim exon
             print "\n\ncheck4\n\n" ; 
          }
        }       
      }
    } 
  } 
  print "processing simgws again\n" ; 

  SIMGW: for my $simgw (@$simgw_tr) {  
    print "processing simgw-gene :\n-------------------------------------------------\n" ;  
    print_transcripts_and_exons($simgw) ; 
    print "\n\nSUMMARY:\n--------------------------\n" ; 
    print "simgw_test nr exons overlapped ny est : " . $simgw->nr_exons_overlapped_by_est ."\n" ; 
    print "nr of all simgw exons                 : " . scalar(@{$simgw->get_all_Exons}) ."\n" ;

   # support from EST at 5prim end 
    print "simgw_test has 5prim support          : " ; 
    if ($simgw->has_5prim_support) { 
       print "OK\n" ;   
       print "simgw_test can extended to 5prim end  : ". $simgw->extend_5prim_end . "\n"  ;
    } else {
      print "---\n" ;
    } 

    print "simgw_test has 3prim support          : " ;
    if ($simgw->has_3prim_support) {
      print "OK\n" ; 
      print "simgw_test can extended to 3prim end  : " . $simgw->extend_3prim_end . "\n" ; 
    }else { 
      print "---\n" ;
    }    
    # support by different est's 
    my %tmp = %{ $simgw->different_est_support} ;
    my $nr_of_ests_which_support_simgw = scalar( keys %{$simgw->different_est_support }) ; 
    print "simgw_test \#nr ests supporting simgw  : $nr_of_ests_which_support_simgw \n" ; 
    if ($nr_of_ests_which_support_simgw == 1) { 
      print "nr of all est exons                   : " . scalar(@{$simgw->get_all_Exons}) ."\n" ;
    }
    print "\n\n\n\n" ; 

    
    # WORK 
    if ( $write_alt_trans ) {  
       print "Constructing alternative transcripts out of simgw-genes and UTR from EST\n\n" ; 
   
       if ($nr_of_ests_which_support_simgw == 1 ) {  
         
       }      

  
       # push @keep, $transcript ;   
       
    }
  }
  



  # case could be that one gene is overlapped by mult.partial others
  #
  # cases : 
  # - simgw is flanked by est on 1 or two terminal ends and middle
  # - est overlaps whole simgw and can be used to add UTR to simgw as well as stand-alone
  # - simgw hits very partial simgw 
  
  # testing if there is more than one est overlapping the SimGW 
  
  my %simgw_transcripts ;

 
#  # testing simgw-transcripts against est-transcripts
#  for my $sim_tr (@$simgw) {
#    for my $se (@{$sim_tr->start_Exon->cluster->get_all_Exons_in_ExonCluster}){
#      print $se->ev_set . "\n" ;  
#    }
#    for my $se (@{$sim_tr->end_Exon->cluster->get_all_Exons_in_ExonCluster}){
#      print $se->biotype . "\n" ;  
#      print $se->ev_set . "\n" ;  
#    }
#  }
#   for my $tr (@all) {  
#    for my $e (@{$tr->get_all_Exons} ) { 
#      print "cluster : " . $e->cluster . "\n" ; 
#    }
#  }
#
  return [\@keep, \@removed ] ;  
}




sub sort_transcripts_by_length { 
  my ($t) = @_ ; 
  my @result = @$t ; 

  if ($result[0]->seq_region_start eq '1' ) { 

  @result = sort { ($b->seq_region_end - $b->seq_region_start) 
                   <=> ( $a->seq_region_end - $a->seq_region_start ) } @result ;  
  } else { 
  @result = sort { ($a->seq_region_start - $a->seq_region_end) 
                   <=> ( $b->seq_region_start- $b->seq_region_end ) } @result ;  
  }
 return \@result ; 
}

sub check_if_lg_spans_st { 
  my ($lg, $st ) = @_ ; 
  #
  # check if long transcripts overlaps short transcript 
  # by looking at it's genomic extend
  # 
  if ($lg->seq_region_start <= $st->seq_region_start 
  && $lg->seq_region_end >= $st->seq_region_end
  && $lg->seq_region_strand == $st->seq_region_strand 
  && (@{$lg->get_all_Exons} >= @{$st->get_all_Exons } )
  && $lg ne $st  ) {   
     return 1 ;  
  }
  if ($st->seq_region_start <= $lg->seq_region_start 
  && $st->seq_region_end >= $lg->seq_region_end
  && $st->seq_region_strand == $lg->seq_region_strand 
  && (@{$st->get_all_Exons} >= @{$lg->get_all_Exons } )
  && $lg ne $st  ) {   
     return 1 ;  
  }
 
 

#
#   print "main_comparision_rules failed for a vs b :" ;  
#   print " SAME OBJECTS " if $lg eq $st ; 
#   print "\n" ; 
#   print "(a) (".scalar(@{$lg->get_all_Exons}) ." exons) "   ; print_object($lg) ; 
#   print "(b) (".scalar(@{$st->get_all_Exons}) ." exons) "   ; print_object($st) ; 
#   print "\n\n" ; 
   return 0 ; 
} 






sub remove_overlapping_transcripts {
  my ($tref ) = @_ ;
  
  my @transcripts ; 

  # only remove transcripts which are not already deleted 
  for my $t (@$tref) { 
     push @transcripts, $t unless ($t->biotype=~m/^del_/) ; 
  } 
 
  # sort by number of exons if tr have more than one exon 
  
  @transcripts  = sort { scalar(@{$a->get_all_Exons}) <=> scalar(@{$b->get_all_Exons}) } @transcripts;

  # case of having single exon transcripts sort them by length 
  if ( @{$transcripts[-1]->get_all_Exons} == 1){   
    @transcripts  = sort { $a->length  <=> $b->length } @transcripts;
  }

  # since we call this function out of the gene-cluster routine, 
  # all genes should be in a gene cluster which makes 
  # preocessing much easier ! 
  #
  # idea :  sort transcripts accourding to their number of exons and 
  # than begin with the shortest transcript and check if it's exon is 
  # covered by one one the longer transcripts 
  # we could also implement some other filters here 
  # FILTERS 
  # sort transcripts acc. to their number of exons 

  my @ltr = reverse @transcripts ;  # tr with most exons first 

  my @remove_tr ; 

  # MAIN LOOPS : loop through longer transcripts 
  for (my $i=0 ; $i<@ltr; $i++) { 
    my $lg = $ltr[$i] ; 
    # loop through shorter transcripts 
    for (my $j=$i+1 ; $j < @ltr; $j++ ) { 

    #print  " checking $i $j \n" ; 

      my $st = $ltr[$j] ; 
      if  ($lg ne $st) { 
        # check genomic extend , nr_exons and if  $lg eq st 
        if ( check_if_lg_spans_st ($lg,$st) ) {  
#          print "\n\nLG spans ST\n" ; 
#          print "comparing transcript lg (". scalar( @{ $lg->get_all_Exons } ) . " exons) " ;
#          print_object($lg) ;  
#          print "comparing transcript st (". scalar( @{ $st->get_all_Exons } ) . " exons) " ; 
#          print_object($st) ;  
#          print "#"x80; print "\n\n" ; 
#
#          print "\n lg\n" ; print_transcripts_and_exons($lg) ;      
#          print "\n st\n" ; print_transcripts_and_exons($st) ;      
#          print "\n\n" ; 
#
          if ( check_if_all_exons_are_overlapped ($lg,$st) ) { 
            # transcript is incorprated in other transcript all exons have been checked 
            #print "\nTranscript_will be removed : " ; print_object($st) ;  print "\n\n" ; 
            push @remove_tr, $st ; 
          }
        } else { 
#         print "lg does not span st\n" ;   
        }
      } else{ 
#        print "same transcript\n" ; 
      }
    }
  }

  
  my %seen ; 
  @seen{ @$tref } =  @$tref  ;
  delete @seen{ @remove_tr } ; 

  my @diff = values %seen ; 

  #
  # change biotype of the transcripts which are removed to get overview what has 
  # been removed  / remove redundant transcript out of @remove_tr before 
  #
  my %tmp ; 
  @tmp{@remove_tr} = @remove_tr ; 
  @remove_tr = values %tmp ; 

  # labelling the exons / transcripts 

  for my $t ( @remove_tr ) { 
    my $bt = "del_" . $t->biotype ; 
    $t->biotype($bt) ; 
    for my $e (@{ $t->get_all_Exons} ) { 
    # $e->biotype($bt) ;  
    }
  }
#  print "xxx  these tr will be removed\n===========================================\n" ; 
#  for (@remove_tr) { 
#   print_transcripts_and_exons($_) ; 
#  }
#
#  print "xxx  these tr will be NOT removed\n===========================================\n" ; 
#  for (@diff) { 
#   print_transcripts_and_exons($_) ; 
#  }
#
  return (\@diff, \@remove_tr) ;
}



# work 
sub check_if_all_exons_are_overlapped {
  my ($lg,$st) = @_ ; 

  # check if each exon of $lg has a counterpart in $st  
  my $test_ex_ok = 0 ;
  my %mark_exon ; 
  my %overlapped_test_exons  ;  
  my $nr_same_exons_conserved = 0 ; 
  my $exon_3prim_mismatch = 0 ; 
  my $exon_5prim_mismatch = 0   ; 

  for my $lge (@{$lg->get_all_Exons }) { 
    for my $ste (@{$st->get_all_Exons }) { 
       #print "testing (1)" ; print_object($lge) ;
       #print "   vs.  (2)" ; print_object($ste) ;
     
      if ( $lge->seq_region_start <= $ste->seq_region_start
      && $lge->seq_region_end >= $ste->seq_region_end
      && $lge->seq_region_strand eq $ste->seq_region_strand )
      {
        #
        # THIS EXON $ste is completely overlapped 
        # |xxxxx|   OR |xxxxxxxx| OR |xxxxxxxxx| OR |xxxxxx|
        # |xxxxx|      |xxxxx|         |xxxx|         |xxxx|
        #
         #print "long_exon (1)--> overlaps test_exon(2)\n" ;
         #print "(1)" ; print_object($lge) ;
         #print "(2)" ; print_object($ste) ;
        $nr_same_exons_conserved++ ; 
        $test_ex_ok++ ;
        $mark_exon{$lge->hashkey} = $lge ; 
        $overlapped_test_exons{$ste->hashkey} = $ste ;  

      } else {
        # no exact match 
        ##print "no exact match between exons\n" ; 
        # exon is a termial exon of ST 
        
        # watch out - here we compare each time the terimal exon of ste against all exons of $lge 
        # so we have to make sure that we have somehow an overlap or a distance between the 
        # exons compared !!!!!!! 
        #
        # perhaps we should just check if ste and lge overlap in their genomic extend than 
        # using an offset of 50 bp for this ! 
        #
        #
        my $start_diff = abs( $ste->seq_region_start - $lge->seq_region_start ) ; 
        my $end_diff = abs( $ste->seq_region_end - $lge->seq_region_end ) ; 
       
        # this is an extra rule for terminal exons  
        if ($ste->is_terminal_exon) {
           #print "exon is terminal : " ;  
          if ($start_diff < 25 && $end_diff < 25 ) { 
           #   print "start_end_diff smaller 50\n" ; 
            #
            # exon is 'in range 'of 50 bp to lge 
            # check how conserved the exon is   
            #
            my $exon_conservation = get_percentage_exon_conversation_in_exon_cluster($ste) ;
           # print "checking percentage conservation exon_conservation: $exon_conservation \n\n" ; 

            # if ( $exon_conservation <  0.1 && $exon_conservation != 0  )    
            if ( $exon_conservation <  0.1 ) { 
              $mark_exon{$lge->hashkey} = $lge ; 
              $overlapped_test_exons{$ste}=$ste ; 
              $test_ex_ok++ ; 
           #    print "test_ex_ok $test_ex_ok\n" ; 
            }else{
           #    print " percentage consersation too high or zero$exon_conservation\n" ; 
            }
          }
        } else {
          #
          # Exon is not terminal / if there is a boundary mismatch in 
          # an internal exon 
          #
          my $conservation = get_percentage_exon_conversation_in_exon_cluster($ste) ;
          if ($conservation < 0.1) { 
          #  print "not really conserved boundary\n" ; 
          }    
        }
  
        # new rule to remove transcripts which match all exons except the terminal one's 
        $start_diff = abs( $ste->seq_region_start - $lge->seq_region_start ) ; 
        $end_diff = abs( $ste->seq_region_end - $lge->seq_region_end ) ; 
        
        if ($ste->is_5prim_exon  ) { 
          # print "exon_is_5prim_exon\n" ; 
          # print_object($ste) ; 
           $exon_5prim_mismatch = 1 ; #if ( $start_diff < 50 || $end_diff < 50 ) ; 
        } 
        if ($ste->is_3prim_exon) { 
           $exon_3prim_mismatch = 1 ; #if ( $start_diff < 50 || $end_diff < 50 ) ; 
        } 
      } 
    }
  }
#  print "\n\n" ; 
#  print "all_exons_tested : test_ex_ok $test_ex_ok \n" ; 
#  print "all_exons_tested : exon_3prim_mismatch  $exon_3prim_mismatch  \n" ; 
#  print "all_exons_tested : exon_5prim_mismatch  $exon_5prim_mismatch  \n" ; 
#  print "all_exons_tested : nr_same_exons_conserved $nr_same_exons_conserved  \n" ; 
#
if ($test_ex_ok == scalar(@{ $st->get_all_Exons } ) ) {   
    # print "all exons of test-gene are overlapped by anotherone \n" ; 
    #
    # all exons of $st are overlapped by $lg 
    # now check for exon-skipping 
    # xxxxxxxxxxxxx-------------xxxxxxxxxxxxxxxx------------xxxxxxxxxxxxxxx
    # xxxxxxxxxxxxx-----------------------------------------xxxxxxxxxxxxxxx
    #
    for my $lge (@{$lg->get_all_Exons }) { 
      unless (exists $mark_exon{$lge->hashkey}) {  
        # if there are exons which are on the long one but not in the short one 
        
        # now check if mark-exon is outside of boundaries of st-transcript
        if  ($st->seq_region_start <= $lge->seq_region_end 
        && $st->seq_region_end >=$lge->seq_region_end ) { 
          # overlap of skipped exon and transcript 
          # print "checkxxx not removed  " ; print_object($st) ; 
          return 0  ;
        }else {
         # the skipped exon is outside of the short transcript 
        }   
      }
    }  
    return 1 ; # remove tr 
  } elsif (  (scalar(@{ $st->get_all_Exons } )  - $nr_same_exons_conserved) < 2 ) { 
    if ($exon_3prim_mismatch || $exon_5prim_mismatch) {
      #print "removing_tr\n" ; 
      #print_transcripts_and_exons($st) ; 
      return 1 ;  # remove tr 
    }
  }
  # print "not all exons are ok\n" ; 
  return 0 ;   #  dont remove transcript / not all exons are overlapped 
} 



sub get_percentage_exon_conversation_in_exon_cluster {
  my ( $exon ) = @_ ;  

  my @ex_clust = @{ $exon->cluster->get_all_Exons_in_ExonCluster } ; 

  # uniquify exon acc. to their hashkey 
  my %uniq_exons ; 
  for (@ex_clust) { 
    unless ($_->is_terminal_exon) {
      $uniq_exons { $_->hashkey }++ ; 
    }
  }

  # get maximum value (most consered boundaries) out of hash
  my @tmp = reverse sort (values %uniq_exons ) ;  

  my $max_cons = -1; 
  if (scalar(@tmp) > 0) { 
    $max_cons = shift @tmp ; 
  } 

  my $percentage_exon_conservation = 0 ; 
  if  ($max_cons > 0  ) { 
    if ($uniq_exons{$exon->hashkey} ){
     $percentage_exon_conservation = $uniq_exons{$exon->hashkey} / $max_cons ;  
    }
  } else {
   $percentage_exon_conservation = 0 ; 
  }
  return $percentage_exon_conservation ;  
}



sub print_transcript_info_quest {
  my ( $all_tr_ref ) = @_ ;  
  my $cnt_all=1 ; 
  for my $t (@$all_tr_ref) { 
    my $cnt_region = 1 ; 
    $cnt_region++ ; 
    $cnt_all++ ; 
  } 
}



sub merge_transcripts_recursion {
  my ( $self, $start_transcript , $transcripts_ref , $exon_clusters ) = @_ ;
      #print "recusrion start\n" ; 

  my $tr1 = $start_transcript ; 
  my $tr2 ;  

  if (@$transcripts_ref > 0 ) { 
    $tr2 = shift @$transcripts_ref ; 
    my @new_transcripts ;  
    my $ex_3prim = $tr1->end_Exon(); 
    my $merged_transcript ; 

    unless ($ex_3prim->next_exon ) {
  
      if ( $tr2->start_Exon->overlaps($ex_3prim) && !$tr2->start_Exon->prev_exon ) {
          # 
          # exon is in same cluster and has terminal 5prim exon of est
          #
          # xxxxx------------xxxxxxxxxxx                             tr1 (no next_exon) 
          #                     yyyyyyyyyyy----------yyyyyyyyyyyyyy  tr2 (no prev_exon) 
     if ($self->{v}){  
          print "no next exon - look for overlapping 5prm \n" ; 
          print "\nRECURSION : These transcripts could be merged:\n" ; print "="x80;  
          print "\ntr1 : " . $tr1->seq_region_start .  " -- " .$tr1->seq_region_end  .
            " " . $tr1->biotype . "\n" ;  ; 
          print "tr2 : " . $tr2->seq_region_start .  " -- " .$tr2->seq_region_end  . 
           " " . $tr2->biotype . "\n" ;  ; 
        } 
        my $this_ec = $tr1->end_Exon->cluster  ; 
          #
          # get evidence to merge $tr1->end_Exon AND  $tr2->start_Exon  
          #

        my @candidate_exons =  @{ $this_ec->get_all_Exons_of_EvidenceSet('simgw')}  ;
        push @candidate_exons, @{ $this_ec->get_all_Exons_of_EvidenceSet('abinitio')}  ;
  
          #
          # get ab-initios or sim-gws in this cluster which match the boundaries 
          #
          # aaaaaaaaaaaa|--------------|aaaaaaaaaaaaaaaa
          #                            |$tr1->end_Exon
          #
          #                                     bbbbbbbbbbbbbbbb|-------|bbbbbbbbbbbbbb--
          #                                     $tr2->start_Exon|
          #
          # simsimsimsim|--------------|simsimsimsimsimsimsimsim|-------|simsimsismisms--  CONSERVED, OK
          # abinitioabin|--------------|itioabinitioabinitioabin|-------|initioabinitio--  CONSERVED, OK
          #
          #
      
        for my $ce (@candidate_exons) { 
          if ($self->{v}){  
            print "RECURSION checking candidate exon which spans region :\n" ; 
            print "start : " . $ce->seq_region_start . "\t". 
            $ce->seq_region_end .  "  (" . $ce->analysis->logic_name . ")\n\n" ; 
           }
  
              if ( $ce->prev_exon && $ce->next_exon ) {
  
                # check the boundarie of previous exon of 3prim term exon of EST-1-gene : 
  
                my ( $tr1_boundary , $tr2_boundary  ) ;
                my $ce_start = $ce->seq_region_end ; 
                my $ce_end = $ce->seq_region_start ; 
  
  
                if ( $tr1->seq_region_strand eq '1' ){
  
                   $tr1_boundary = $tr2->start_Exon->seq_region_end ;  
                   $tr2_boundary = $tr1->end_Exon->seq_region_start ;  
  
                } else { 
                   $tr1_boundary = $tr1->end_Exon->seq_region_end ;  
                   $tr2_boundary = $tr2->start_Exon->seq_region_start ;  
                } 
                if ($self->{v}){  
                  print "check_exon has previous and next exon\n" ; 
                  print "\n" ; print "tr1_boundary  : " . $tr1_boundary . "\n" ; 
                  print "ce_start      : " . $ce_start . "\n" ; 
                  print "\n" ; print "tr2_boundary  : " . $tr2_boundary . "\n" ; 
                  print "ce_end        : " . $ce_end . "\n" ; print "\n\n" ; print "\n\n" ; 
                } 
  
                if ($ce_start == $tr1_boundary && $ce_end == $tr2_boundary  ) { 
                  print "Boundaries match - merging genes \n" if $self->{v}; 
                    # less strict way in finding a compatible exon : 
                  # only the intron has to be conserved between the flanking ab-initios and the est's 
                  # (ONLYE THE EXON BOUNDARIES WHICH HAVE DOUBLE '||' ARE CHECKED !!!) : 
                  #
                  #                    TR_1                                       TR_2
                  # --------estestest||-------------||estestestest 
                  #                                             estestestest||---------------||estestestest|---------- 
                  # --|abinitioabiont||-------------||abintioabintioabinitio||---------------||abinitioabinitioabinitio|----
                  #                  ^^  CONS.INTR. ^^                      ^^  CONS. INTRON ^^
                  #                  ||             ||                      ||               || 
                  # NOMATCH        match           match                   match            match           NOMATCH 
                  #
                  #
                  # It is only checked if there are two conserved Introns at 3prim and 5rim end 
                  # of 'broken' exon which we want to bridge.  The more 'strict' approach is handeld in the 
                  # second case below (look if the last exon of est-gene nr1 and next exon of est-gene nr2 
                  # have excatly the same coordinates as the simgw / abinitio exon we use to merge
                  # 
                  #
                  # 
                  # merge the two transcripts by using the coordinates of simgw/abinitio bridge exon 
                  #
  
                  # bondaries of abinotio/simgw exon match the start/ends of the exon to bridge  
                  my $tr1_last_intron = Bio::EnsEMBL::Intron->new($tr1->end_Exon->prev_exon,$tr1->end_Exon) ; 
                  my $ce_last_intron = Bio::EnsEMBL::Intron->new( $ce->prev_exon,$ce ) ; 
    
                  my $tr2_next_intron = Bio::EnsEMBL::Intron->new($tr2->start_Exon,$tr2->start_Exon->next_exon) ; 
                  my $ce_next_intron = Bio::EnsEMBL::Intron->new($ce, $ce->next_exon ) ; 
    
                  if ( compare ($tr1_last_intron, $ce_last_intron ) && compare ($tr2_next_intron, $ce_next_intron)) { 

                     push @new_transcripts , $self->merge_transcripts ( $tr1, $tr2, $ce, 'no_strict') ; 
                     $merged_transcript = $self->merge_transcripts ( $tr1, $tr2, $ce, 'no_strict') ; 

                     $self->merge_transcripts_recursion ( new_tr($merged_transcript) , $transcripts_ref , $exon_clusters ); 
                     # $self->merge_transcripts_recursion ( $merged_transcript , $transcripts_ref , $exon_clusters ); 
                  }
    
                  if ( $tr1->end_Exon->prev_exon->seq_region_end == $ce->prev_exon->seq_region_end
                    && $tr1->end_Exon->prev_exon->seq_region_start == $ce->prev_exon->seq_region_start  
                    && $tr2->start_Exon->next_exon->seq_region_start == $ce->next_exon->seq_region_start  
                    && $tr2->start_Exon->next_exon->seq_region_end == $ce->next_exon->seq_region_end ) {
  
        
                  # STRICT WAY in finding a compatible exon : 
                  # we need two flanking consreved Introns and two exons which match exactly the boundaries !! 
                  # (all double || marked starts/ends are checked 
                  #
                  # ---||esteststesteste||-------------||estestestest 
                  #                                                estestestest||---------------||estestestest||---------- 
                  # ---||abinitioabiniot||-------------||abintioabintioabinitio||---------------||abinitioabin||----------
                  #    ^^               ^^  CONS.INTR. ^^                      ^^ CONS. INTRON  ^^            ^^
                  #    ||               ||             ||                      ||               ||            ||
                  #   match           match          match                   match            match         match 
                  #
                  #
                  # It's checked if there are two conserved Introns at 3prim and 5rim end around the   
                  # 'broken' exon which we want to bridge and if the start/ends of the flanking exdons match as well
                  #
                 
                  # 
                  # merge the two transcripts by using the coordinates of simgw/abinitio bridge exon 
                  #


                    push @new_transcripts , $self->merge_transcripts ( $tr1, $tr2, $ce, 'strict' ) ; 
                    $merged_transcript = $self->merge_transcripts ( $tr1, $tr2, $ce, 'strict') ;
                    $self->merge_transcripts_recursion ( new_tr($merged_transcript) , $transcripts_ref , $exon_clusters ); 
  
                  } 
                } else { 
                  print "Boundaries does not match\n" if $self->{v}; 
                }
              } 
            } # for candidate exons 
          } 
        }  
      #### 
     } else { # there is no more transcriptin $transcripts_ref
       push @{$self->{merged_tr} } , $start_transcript ; 
       return  ; 
     }
   my $tr_new = new_tr($tr2) ;  
   $self->merge_transcripts_recursion ( $tr_new  , $transcripts_ref , $exon_clusters ); 
}






sub merge_transcripts {
  my ($self, $base_trans, $trans_to_add , $merge_exon , $string) = @_ ; 
  
 
  my $new_biotype = $self->{new_biotype} . "_" . $merge_exon->biotype . "_merged_$string" ; 
  # my $new_biotype = $self->{new_biotype} . "_abinitio_merged" ;

  #
  # we need to clone exon otherwise we end up with multipe transcripts which share the same Exon object 
  #

  my $cloned_ex = new Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended ( 
                     -start  => $merge_exon->start , 
                     -end    => $merge_exon->end , 
                     -phase  => $merge_exon->phase ,  
                     -end_phase  => $merge_exon->end_phase, 
                     -strand => $merge_exon->seq_region_strand , 
                     -slice  => $merge_exon->slice , 
                     -analysis => $base_trans->analysis ) ;  

    $cloned_ex->biotype( $new_biotype ) ; 
    $cloned_ex->transcript($merge_exon->transcript) ; 
    $cloned_ex->add_supporting_features( @{$merge_exon->get_all_supporting_features}) ;                
  #
  # register new biotype as a member of evidence_set 'est_merged'  and make sure that the new 
  # biotype is only stored once in the array of biotypes , otherwise gene clustering won't work
  #
  my %tmp = %{ $self->get_all_evidence_sets }; 
  my %tmp2 ; 
   @tmp2{ @{ $tmp{est_merged} } } =(); 
   unless (exists $tmp2{ $new_biotype } ) {  
     push @{ $tmp{'est_merged'}}, $new_biotype ; 
     $self->get_all_evidence_sets( \%tmp ) ;  
   }



  my $old_merge_exon_biotype = $merge_exon->biotype ;  
  if ( length($new_biotype) > 40) { 
    warning("The stringlength of the biotype is too loong....mysql will shorten it\n") ;  
  } 
   
  my @base_exons = @{$base_trans->get_all_Exons()} ; 

  #
  # chop of 3prim exon cause here we use the abintio/ simgw exon for merging  
  #
  pop @base_exons ; 

  my @exons_to_add = @{ $trans_to_add->get_all_Exons}  ; 

  #
  # chop of 5prim exon cause here we use the abintio/ simgw exon for merging  
  #
  shift @exons_to_add ;
  
  
  my $tr_merged  = new Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended(
                     -BIOTYPE => $new_biotype , 
                     -ANALYSIS =>$base_trans->analysis , 
                   ); 
  



  #
  # creating transcript supporting feature 
  #
  my $hseq_name ; 
  
  if ($merge_exon->stable_id) {
    $hseq_name = $merge_exon->stable_id ; 
  }elsif ($merge_exon->dbID) { 
    $hseq_name = $merge_exon->dbID ; 
  }else {
    $hseq_name = $merge_exon->biotype ; 
  }

  my $feat_pair = Bio::EnsEMBL::FeaturePair->new (
                    -start => $merge_exon->seq_region_start , 
                    -end   => $merge_exon->seq_region_end , 
                    -strand => $merge_exon->seq_region_strand ,  

                    -slice => $merge_exon->slice , 
                    -analysis =>$merge_exon->analysis , 

                    -hseqname =>$hseq_name , 
                    -hstart  => 1  ,  
                    -hend  => ( $merge_exon->seq_region_end - $merge_exon->seq_region_start +1 )  ,  
                    -hstrand  => $merge_exon->seq_region_strand  ,  
                    -score  => 0  ,  
                    -percent_id  => 0  ,  
                    ) ; 

  my $dna_align_feat = Bio::EnsEMBL::DnaDnaAlignFeature->new (-features =>[$feat_pair] , 
                                                              -analysis => $merge_exon->analysis ) ; 
                    

  $cloned_ex->add_supporting_features ( $dna_align_feat) ;

  #
  # Doing first part of merge  
  # 

  for my $be (@base_exons ) { 
    $be->biotype( $new_biotype ) ;  
    $tr_merged->add_Exon( $be ) ; 
    #    print "add base_exons to tr: " ; print_object($be) ; 
  }
 
  #    
  # add exon of abinitio / simgw 
  #
  
  #$merge_exon->biotype( $old_merge_exon_biotype ) ; 
  $tr_merged->add_Exon ( $cloned_ex )  ; 
  # $tr_merged->add_Exon ( $merge_exon )  ; 
  # print "add merge_exon to tr: " ; print_object($merge_exon) ; 



  #
  # adding rest of following est
  #
  for my $e (@exons_to_add) {
   $e->biotype($new_biotype) ; 
   #print "adding rest of second est: " ; print_object($e) ; 
    $tr_merged->add_Exon ($e) ; 
  }
  
  #print "\n\nThis_is_merged_trans :\n" ; 
  #print_transcripts_and_exons([$tr_merged],"tr_merged") ; 

  return $tr_merged ; 
}





sub edit_terminal_exons {
  my ($self,$aref) = @_ ; 
  for my $t (@$aref ) {
    my @exons = @{ $t->get_all_Exons } ; 

    # process 5prim terminal exon 
    my $ex_5prim = $exons[0];

    my $longest_exon_5 = $self->get_longest_5prim_term_exon_in_exon_cluster_of_this_exon($ex_5prim) ; 
    
    $t = $t->exchange_exon($ex_5prim, $longest_exon_5) ; 
     
    my $ex_3prim = $exons[$#exons]; 
    my $longest_exon_3 = $self->get_longest_3prim_term_exon_in_exon_cluster_of_this_exon($ex_3prim) ; 
    $t = $t->exchange_exon($ex_3prim, $longest_exon_3) ; 

  }

  print_transcripts_and_exons($aref) if $self->{v} ; 
  return $aref ; 
}


sub get_longest_5prim_term_exon_in_exon_cluster_of_this_exon { 
  my ($self, $small_exon ) = @_ ; 
  my $ec = $small_exon->cluster; 

  my @all_exons_in_clust = @{ $ec->get_all_Exons_in_ExonCluster } ; 
  my $longest_term_exon_in_clust = $small_exon  ; 

  for my $exon (@all_exons_in_clust) {

    if (!$exon->prev_exon && !$exon->next_exon) {
      warning("Skipping this exon because this is a single exon\n") ; 
      next ; 
    }

    if (!$exon->prev_exon && ($exon->length>$longest_term_exon_in_clust->length) ) { 

      if ($exon->seq_region_strand eq "1" ) { 
        if ($exon->seq_region_end == $small_exon->seq_region_end) {

          # we only want to include the most 5prim exon, no overlapping exons 
          # (extending b with a is ok but not with c) 
          #      aaaaaaaaa---aaa-------xxxxxxxx
          #           bbbb---bbb-------xxxxxxxx
          # cccccccccccccccccccc--------xxxxxxx

          $longest_term_exon_in_clust = $exon ;
        } 

      }elsif ($exon->seq_region_strand eq "-1") {
        if ($exon->seq_region_start == $small_exon->seq_region_start ) {
          $longest_term_exon_in_clust = $exon ;
        }
      }else {
#        $self->print_exon($exon) ; 
        throw("exon has no seq_region_strand assigned. exiting\n" ) ; 
      }
    }
  }
  return $longest_term_exon_in_clust ; 
}



sub get_longest_3prim_term_exon_in_exon_cluster_of_this_exon { 
  my ($self, $small_exon ) = @_ ; 
  my $ec = $small_exon->cluster; 

  my @all_exons_in_clust = @{ $ec->get_all_Exons_in_ExonCluster } ; 
 
  my $longest_term_exon_in_clust = $small_exon  ; 

  for my $exon (@all_exons_in_clust) {

    if (!$exon->prev_exon && !$exon->next_exon) {
      warning("Skipping this exon because this is a single exon\n") ; 
      next ; 
    }

    if (!$exon->next_exon && ($exon->length>$longest_term_exon_in_clust->length) ) { 

      if ($exon->seq_region_strand eq "1" ) { 
        if ($exon->seq_region_start == $small_exon->seq_region_start) {
          $longest_term_exon_in_clust = $exon ;
        } 
      }elsif ($exon->seq_region_strand eq "-1") {
        if ($exon->seq_region_end == $small_exon->seq_region_end ) {
          $longest_term_exon_in_clust = $exon ;
        }
      }else {
#        $self->print_exon($exon) ; 
        throw("exon has no seq_region_strand assigned. exiting\n" ) ; 
      }
    }
  }
  return $longest_term_exon_in_clust ; 
}


=head

Name: $self->get_exon_clustering_from_gene_cluster( $gene_cluster ) 
Arg :  Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneCluster;
Function  : gets a set of clustered genes and clusters the Transcripts of these genes as well as the Exons
Returnval :  Aref of  Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonCluster objects

=cut  

sub get_exon_clustering_from_gene_cluster {
  my ($gene_cluster) = @_ ;

  my @clg  = sort {$a->start <=> $b->start} $gene_cluster->get_Genes ; 

  # building Transcript-Cluster 
  my $tc = genes_to_Transcript_Cluster(\@clg); # UtilsFunc.pm
#   print "genes_to_transcript_cluster done\n" ; 
 
  my @exon_clusters = cluster_exons_in_transcript_cluster($tc); #UtilsFunc.pm
#   print "genes_to_transcript_cluster done\n" ; 
  
  if ($tc->strand eq '1') {
    @exon_clusters = sort { $a->start <=> $b->start } @exon_clusters ; 
  } else {  
    @exon_clusters = sort { $b->start <=> $a->start } @exon_clusters ; 
  }
  return \@exon_clusters ;  
}




sub main_clustering_and_recursion {
  my ($self, $gene_cluster) = @_ ;

  my @exon_clusters = @{ get_exon_clustering_from_gene_cluster ( $gene_cluster ) } ; 

  my @ex_cluster_est_evidence = @exon_clusters ; 
  $self->{start_terminal_ec} = $ex_cluster_est_evidence[0];
  $self->{end_terminal_ec} = $ex_cluster_est_evidence[$#ex_cluster_est_evidence];
  
  # use all possible exons with EST-evidence in first ExonCluster 

  # this is wrong - it hand's over wrong exons !!!
  #

  my @est_start_exons = @{ $ex_cluster_est_evidence[0]->get_all_Exons_of_EvidenceSet('est') } ; 

  # make sure we don't have 3prim exons in the start-exon-set ! (only use exons with have next exon)
  
  @est_start_exons  = grep {$_->next_exon} @est_start_exons ;   

  my @all_assembled_tr;
  
  #
  # we could do this stuff as well with a recursion / 
  # idea would be : 
  # check first eand second exon to avoid situation where it's like : 
  # xxxxx----------------------------------------------------xxxxxxxxxx
  #            xxxxxx--------xxxxxxxxx--------xxxxx-xxxxxxx--------------x-xxxxxxxxxxxxxxx
  #          xxxxxx--------xxxxxxxxx--------xxxxx-xxxxxxx--------------x-xxxxxxxxxxxxxxx
  # 
  # 
  #  ---> WHERE NICE GENES START IN THE SECOND EXON CLUSTER 
  #  implement a restarting of the algorithm ! 
  #
  #
 
  # returns array of arrays of transcripts  [ [ tra tra tra ] [ tra ] 
  my $trans_aref   = $self->transcript_recursion(
                                                       \@ex_cluster_est_evidence,
                                                       \@est_start_exons, 
                                                       []
                                                       ) ; 
 
  #
  # prune transcript (change biotype and re-adjust start/end of transcript) 
  #

  my @pruned_transcripts = @{ re_adjust_transcripts( $trans_aref , $self->{new_biotype}) } ; 

  return \@pruned_transcripts ; 
}



sub re_adjust_transcripts{ 
  my ($tref, $biotype) = @_ ; 

   my @pruned_transcripts ;   
   for my $transcript ( @$tref ) { 
    my $ntr = new Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended( -BIOTYPE => $biotype ) ; 

    for my $exon ( @{ $transcript->get_all_Exons } ) { 
      $exon->biotype( $biotype ) ; 
      $ntr->add_Exon( $exon ) ; 
    }
    push @pruned_transcripts , $ntr ; 
  } 
  return \@pruned_transcripts ; 
} 


# does recursion in region with clustered transcripts and selects 
# every unvisited exon as possible start exon for restart of recursion 
#
# i have a seq-range with all transcripts and exon-clusters 
# i have to select the start-exon-clusters if possible 
# i start to process the first start exon clusters 
# if i am finished with this i check if all exons have been visited 
# if not, i sort the exons, use the one at 5prim end and call exon_recursion on this exon 
# i check again if all exons have been visited -> if not, i choose the one which has'nt been 
# visited and call exon_recursion on this......
#

sub start_est_linking { 
  my ( $self,  $start_exons , $ex_cluster_real_ev ) = @_ ; 
  my @all_assembled_tr ; 

  if ($self->{v}){ 
    print "starting est-linking with set of start-exons\n" ; 
    for (@$start_exons) { 
      print "START_exon: " ; print_object($_) ; 
    }
  }

  for my $act_start_exon (@$start_exons) { 
    if ($self->{v}){    
      print "\n\nSTART: Starting with exon :\n" ; 
      $self->print_exon($act_start_exon) ;  
      print "starting the exon_recursion\n" ; 
    }
    my $t  =  $self->exon_recursion( 
                      $ex_cluster_real_ev,        # all ec with exon ranked > 0 
                      $act_start_exon,            # the exon where the story begins 
                       undef,                     # last exon
                       [],                        # transcript  
                     #  [],                        # ref to store exons  
                     #  {},                        # reference to transcript_hash 
                       0 ,                        # counter 
                       {}                         # href 
                     ); 

    push @all_assembled_tr, $t if $t;  
    print "exon_recursion finished" if $self->{v} ; 
  }

  print_transcripts_and_exons(\@all_assembled_tr,"After exon_recursion") if $self->{v}; 
  return \@all_assembled_tr ;  
}




sub transcript_recursion { 
  my ($self , $ex_clusters_real_ev,  $exons ,$all_tr ) = @_ ;
  
  print "\nTranscriptRecursion:\n=============================================\n" if $self->{v}; 
  my @ex_clusters_real_evidence = @$ex_clusters_real_ev ; 
  my @start_exons ; 


  if ( scalar(@$all_tr)==0) { 

   print "initialisation - using start exons which have been handend over _tmp_\n" if $self->{v} ;
   @start_exons = @$exons ;  

  } else { 

    # We have already transcripts and we have to check out for next 
    # cluster with unprocessed / unvistited exons and call the recursion again on this cluster 
    # how to find next startpoint ?  
    @start_exons  = @{ get_unvisited_exons_in_next_ec( $ex_clusters_real_ev ) } ; 

    #
    # sort start-exons such that we don't end up using terminal 3-prim exons
    #
    my @tmp ; 
    if ( ( @start_exons) > 0 ) { 
      for my $se (@start_exons) { 
       if ($se->next_exon) { 
          if ($self->{v}) { 
            print "this exon_has_next_exon : " ; 
            print_object($se) ; 
            print "this is next exon : " ; 
            print_object($se->next_exon) ; 
          } 
          push @tmp, $se if $se->next_exon() ; 
       }
      } 
    }
    @start_exons = @tmp ; 
  }
 
  if ( all_exons_visited ( $ex_clusters_real_ev) ) { 
    print "\nall_exons_visited_returning_from_transcript_recursion\n" if $self->{v};  
    return ; 
  } 
  
  if (@start_exons > 0 ) {
    print "starting est_linking\n" if $self->{v};

    my $trans_aref = $self->start_est_linking(\@start_exons , $ex_clusters_real_ev  )  ;  

    print_transcripts_and_exons($trans_aref , "trans_aref_smart" ) if $self->{v}; 

    push @$all_tr, @$trans_aref ;  

    print "finished est_linking\n" if $self->{v};

    for my $ref (@$all_tr) { 
      if ( ref($ref)=~m/ARRAY/){
       throw("having_array_that is not good");    
      }
    } 
  
  }
  $self->transcript_recursion ($ex_clusters_real_ev, \@start_exons , $all_tr ) ; 

  if ($self->{v}){
    print "END_OF_TRANSCRIPT_RECURSION\ntrying print_transcripts_and_exons\n"  ; 
    print_transcripts_and_exons($all_tr,"debug_trans_rec") ; 
  }
  return $all_tr ; 
} 








sub print_transcripts_and_exons  {
  my ($aref,$string ) = @_ ;
  my @tr_array ; 

  # check if array or only trnascript has been handed over 
  if ( ref($aref) =~m/ARRAY/) { 
   @tr_array = @$aref ;  
  } else {
   push @tr_array, $aref ;  
  }

  $string = "" unless $string ; 
  for my $t ( @tr_array ) { 
    print "$string (Tran): " ;
    print_object($t) ;
    print "\n" ; 
    for my $e ( @{ $t->get_all_Exons} ) { 
      print "$string (Exon): " ; print_object($e);  
    }
    print "\n\n" ;
  }
}



sub exon_recursion { 
  my ($self, $all_ec, $exon, $last_exon, $transcript , $cnt, $href ) = @_ ; 
  #my ($self, $all_ec, $exon, $last_exon, $transcript , $exon_aref, $trans_href,$cnt, $href ) = @_ ; 

  if ($self->{v}) { 
    print "--> Starting exon recursion with this exon : \n\n" ; 
    print_object($exon) ;  
    print "\n\n" ; 
  }
  
  $exon->visited("1") ; 

  my $tmp_biotype = $exon->biotype ."_temp" ; 

  my $init = 0 ;  

  if ( $cnt == 0 ) { 
    print "Initialisation\n" if $self->{v} ;  
    $init = 1 ; 
    $transcript = new Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended(
      -BIOTYPE => $tmp_biotype , 
     ); 
    #$$trans_href{$exon} = $transcript ; 
    $$href{$exon}=$transcript ; 
  }  
  $cnt++ ; 

  # STOP if exon is in end terminal cluster and has NO next exon 
  # exon could be in end-terminal cluster but there could be another exon cause of overlap 
  # xxxxxxxxxxxxxx-----------------------------xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  # xxxxxxxxxxxxx-------------------------------xxxxxxxxxxxx--------------xxxxxxxxxxxxxx

  if ( $exon->cluster eq $self->{end_terminal_ec} && !$exon->next_exon) { 
    #  $transcript->add_Exon($exon) ; 
    print "End terminal cluster reached\n" if $self->{v} ;  
    return $transcript  ; 
  }


  if ( $exon->next_exon) {
    $$href{$exon}=() ; # mark exon as visited
    print "Exon has next exon\n" if $self->{v} ;  
      
    unless ( exon_coordinates_match ($exon, $last_exon) ) {
      #$$trans_href{$exon} = $transcript ; 
      my $tr_tmp ; 

      if ($init) {
        # Initialisation : don't clone transcript, add exon
        $transcript->add_Exon($exon) ; 
        $tr_tmp = $transcript ; 
      }else{
        # Initialisation was done, clone transcript and add exon
        $tr_tmp  = new_tr($transcript) ;  
        $transcript->add_Exon($exon) ;
        #$$trans_href{$exon} = $transcript ; 
      } 
      $self->exon_recursion( $all_ec, $exon->next_exon , $exon,   $tr_tmp  ,  $cnt, $href); 
      #$self->exon_recursion( $all_ec, $exon->next_exon , $exon,   $tr_tmp  ,  $exon_aref, $trans_href, $cnt, $href); 

      # Recursion is finished now  - we have to check if we have to add an end-exon
      my @t_exons = @{ $transcript->get_all_Exons } ;  
  
      # checking if last Exon is really at the very end of transcript
      # or if we can add one more exon ($exon->next_Exon) 
        
      my $exon_to_add ; 

      if ($exon->next_exon) { 
        if ($transcript->end_Exon ne $exon->next_exon) {
            
          # check if exon we want to add is more 3prim than trans->end_Exon
          if ($transcript->seq_region_strand eq 1 ) { 
            if ($transcript->end_Exon->seq_region_end < $exon->next_exon->seq_region_start ) { 
             # xxxxxxxxxxxxxxxxx-------xxxxxx-----------xxxxxxxxxxxxxxxxxxx
              # $t->seq_reg_start                        $t->seq_reg_end
              #                                                                  YYYYYYYYYYYYYYYYYY
              #                                                                  $e->seq_reg_start
              $exon_to_add = $exon->next_exon ; 
              } 
            } else {
              if ($transcript->end_Exon->seq_region_start > $exon->next_exon->seq_region_end ) { 
                $exon_to_add = $exon->next_exon ; 
              }
            }
          }
        }
        if ($exon_to_add ) { 
          $transcript->add_Exon($exon_to_add) ; 
        } 
        return $transcript ;
      } # exon_coords_match 
    } else {   # exon->next_exon()
      print "Exon has no next exon - trying to find conserved intron to extend Transcript\n\n"  
       if $self->{v} ; 
    }


    #
    # exon has no next exon / getting all TEST-exons in exon_cluster
    #
    ####################################################################
    my @all_ex_in_ec = @{ $exon->cluster->get_all_Exons_in_ExonCluster() } ; 

    $$href{$exon}=() ; 
 
    print "Now walking through all exons in Cluster.... ( try to find last matching Intron)\n\n" 
     if $self->{v};


    ALL_EXONS_IN_EC: for my $te ( @all_ex_in_ec ) {

       my $nr_all_exons_in_ec = scalar(@all_ex_in_ec); 
       my $nr_exons_visited = keys %{$href} ; 
      
      if ($te eq $exon) { 
        print "exon is the same ---> skipping / next\n"  if $self->{v};   # exon is the same 
        next ;
      } 
      
      if (exists $$href{$te} )  { 
        # exon has already been processed  
        print "exon has been visited ---> skipping / next\n"  if $self->{v};   # exon is the same 
        next ; 
      }
      if ($self->{v}){  
        print "Now testing this exon for conseverd intron:\n" ; 
        $self->print_exon($te,"  (this exon will be tested for cons. intron) ") ; 
      }

      # INTRONS  
      my ( $te_intron, $ex_intron  ) ; 
      if ($te->prev_exon && $exon->prev_exon) {
        print "\nExon and test-exon have prev. exon - building last intron\n"  if $self->{v}; 
        $te_intron = Bio::EnsEMBL::Intron->new($te->prev_exon,$te) ; 
        $ex_intron = Bio::EnsEMBL::Intron->new($exon->prev_exon,$exon) ; 
      } else { 
        #print "Can't build intron cause test-exon OR exon have no prev exon\n" ; 
      }
      print "Comparing introns if there are introns...\n" if $self->{v};   
      if ( compare($te_intron, $ex_intron ) ) {
        print "Introns match !!\n" if $self->{v}; 
 
        if (!exists($$href{$te})) { 
          print "Exon hasn't been visited\n" if $self->{v}; 
          my $new_tr_2 ;  

          if ($last_exon) {
            if ($self->{v}) {  
              print "Last_exon :  $last_exon \n" ;   
              $self->print_exon($last_exon) ; 
              $self->print_exon($te->prev_exon) ; 
              print "switching transcript bcs conserved intron to this exon:\n"; 
              print_exon($te->prev_exon) ; 
              print "\nCloning transcript -has the foolowing exons now :\n" ;
              print "-------------------------------------------------\n" ; 
            }
            my $new_tr_2 = new_tr($transcript) ; 
            my @all_exons_new_tr = @{$new_tr_2->get_all_Exons} ; 
            for my $e (@all_exons_new_tr) { 
              $self->print_exon($e) ; 
            }  
            print "\n\n" if $self->{v}; 
          } else { 
            $new_tr_2 = new_tr($transcript) ; 
            $last_exon = $exon ; 
          }
          $new_tr_2 = new_tr($transcript) ; 
          $$href{$te}=1 ;                   # mark exon as visited
          ## $$href{$te}=() ;                   # mark exon as visited

          # start new recursion 
          if ($self->{v})  {  
            print "new_starting_point:" ;  $self->print_exon($te) ;  
          }
          if ($te->next_exon) {  
            $self->exon_recursion ($all_ec,$te,$last_exon,$new_tr_2,$cnt,$href) ;
          }else {
            print "test exon has no next exon, skipping to start new recursiox\n" if $self->{v};  
          }
          #$self->exon_recursion ($all_ec,$te,$last_exon,$new_tr_2,$exon_aref,$trans_href,$cnt,$href) ;
          return ; 
        }else { 
          if ($self->{v})  {  
            print "exon has already been visited\n\n" ; 
            $self->print_exon($te,"  VISITED " ) ; 
          }
        }
      }elsif ( ($te->seq_region_end eq $exon->seq_region_end) && $te->next_exon)  { # end match 
        #print "Introns does not_match_1\n" ; 
        #$self->print_exon($te,"six",2); 
      }elsif ( ($te->seq_region_start eq $exon->seq_region_start) && $te->next_exon){  # start match 
        #print "Introns does not_match_2\n" ;
        #$self->print_exon($te,"seven",2); 
      }else {
        #print "Introns does not_match_3.... processing next exon in exon_cluster\n\n" ;
      } 
      $$href{$te}=1 ;  
 
    } #: ALL_EXONS_IN_EC 
  return  ; 
}






sub new_tr {
  my ($old_tr) = @_ ;
  my $new_tr = new Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended(
                 -EXONS => $old_tr->get_all_Exons(), 
                 -BIOTYPE =>$old_tr->biotype(),  
                 ) ; 
   $new_tr->ev_set($old_tr->ev_set) if $old_tr->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended")  ; 
  return $new_tr ; 
}


=head2 cluster_Genes 

Arg[1] : ref to an array of genes to cluster
Arg[2] : ref to hash which builds up an Evidence-Set to biotype-relation :
         $hash{'cdna'}=['cdna_kyoto','cdna_other' ] 
         $hash{'simgw'}=['simgw_100','simgw_200'] 

Function : clusters all genes in Arg[1] according to their genomic extent and 
   sets the type according to their sets

Returntype : Array of 2 Arrayreference : First arrayref holds an array of all genes which 
     have been clustered together, second ref holds an array of genes 
     which were not clustered.  

=cut 





sub cluster_Genes {
  my ($genes, $types_hash) = @_ ; 

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
  my @sorted_genes = 
          sort { $a->start <=> $b->start ? $a->start <=> $b->start  : $b->end <=> $a->end }  @$genes;
  
  print "Clustering ".scalar( @sorted_genes )." genes on slice\n" ;

  my @clusters;
  GENE: foreach my $gene (@sorted_genes) {
  
    my @matching_clusters;
    
    ## 
    ## if there are Clusters (initialisation below) than check  
    ## if gene lies in the boundaries of the cluster and has at least 
    ## one exon which overlaps with an exon of a gene which already 
    ## belongs to the cluster
    ##
            
    CLUSTER: foreach my $cluster (@clusters) {
    
    # 
    # if gene lies in the boundaries of the cluster......
    #
    
      if ($gene->end  >= $cluster->start && $gene->start <= $cluster->end) {
    
        # search for a gene in the cluster which overlaps the new gene, 
      
        foreach my $cluster_gene ($cluster->get_Genes){
      
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
            #         START----------------$cluster_gene----------END 
            #                                              START------$gene-------END
            #
            # add gene target-gene to cluster if it has at least
            # one gene wich overlaps with an exon of the clustered gene 
            # and add to cluster  
            #
                  
            if (_compare_Genes( $gene, $cluster_gene)) {
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
      my $newcluster = Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneCluster->new();
      foreach my $set_name (keys %$types_hash) {
        $newcluster->gene_Types($set_name,$types_hash->{$set_name});
      }
      $newcluster->put_Genes($gene);
      push(@clusters,$newcluster);
    
      #
      # if above was found ONE matching cluster
      #
    } elsif (scalar(@matching_clusters) == 1) {
      $matching_clusters[0]->put_Genes($gene);
  
    } else {
      # Merge the matching clusters into a single cluster
      my @new_clusters;
      my $merged_cluster = Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneCluster->new();
  
      foreach my $set_name (keys %$types_hash) {
        $merged_cluster->gene_Types($set_name,$types_hash->{$set_name});
      }
      
      my %match_cluster_hash;
      foreach my $clust (@matching_clusters) {
        $merged_cluster->put_Genes($clust->get_Genes);
        $match_cluster_hash{$clust} = $clust;
      }
      $merged_cluster->put_Genes($gene);
      push @new_clusters,$merged_cluster;
      
      # Add back non matching clusters
      foreach my $clust (@clusters) {
        if (!exists($match_cluster_hash{$clust})) {
          push @new_clusters,$clust;
        }
      }
      @clusters =  @new_clusters;
    }
  }
  
  # Seperate genes which are UNclustered (only one gene in cluster ) and
  # from clusteres which hold more than one gene 

  my (@new_clusters, @unclustered);
  foreach my $cl (@clusters){
    if ( $cl->get_Gene_Count == 1 ){
      push @unclustered, $cl;
    } else{
      push( @new_clusters, $cl );
    }
  }
  print "All Genes clustered\nGot " . scalar(@new_clusters) . " new Clusters\n"  ; 
  return (\@new_clusters, \@unclustered);
}



=head2 _compare_Genes()

Title: _compare_Genes
Usage: this internal function compares the exons of two genes on overlap
Source : Bio::EnsEMBL::Pipeline::GeneComparison::GeneComparison; 
=cut

sub _compare_Genes {         
  my ($gene1,$gene2,$translate) = @_;
  # quit if genes do not have genomic overlap 
  #
  # start-------gene1------end   start--------gene2----------end
  #  
  
  if ($gene1->end < $gene2->start || $gene1->start > $gene2->end) {
    print "Gene 1  " . $gene1->start . " " . $gene1->end . " \n";
    print "Gene 2  " . $gene1->start . " " . $gene1->end . " \n";
    print "Failed extents check - returning 0\n";
    return 0;
  }
  
  
  # $overlaps = ( $exon1->end >= $exon2->start && $exon1->start <= $exon2-> end );  
  
  if ($translate) {
    # exon-overlap only on coding exons !
    my $exons1 = get_coding_exons_for_gene($gene1);
    my $exons2 = get_coding_exons_for_gene($gene2);
    foreach my $exon1 (@$exons1) {
      foreach my $exon2 (@$exons2) {
        if ( ($exon1->overlaps($exon2)) && ($exon1->strand == $exon2->strand) ){
          #print "Passed CDS overlap check - returning 1\n";
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
        if ( ($exon1->overlaps($exon2)) && ($exon1->strand == $exon2->strand) ){
          #print "Passed exon overlap check (noncod. + cod. exons checked)  - returning 1\n";
          return 1;
        }
      }
    }
  } 
   #print "Failed overlap check (translate = $translate) - returning 0\n";
  return 0;
}      
#########################################################################



sub print_cluster_info{
  my ($slice,$cluster) = @_; 
  my $name = "GENE-CLUSTER : ";
  $name = "EXON-CLUSTER" if ref($cluster)=~m/Exon/ ; 
  
  my $offset = $slice->start -1 ;  #correction for apollo
  my $cl_start = $offset + $cluster->start ; 
  my $cl_end= $offset + $cluster->end ; 
  print "$name  start: $cl_start\tend: $cl_end\t" . $cluster->strand ." \n"  ; 
  return ;
}



sub newAnalysis{
  my ($self,$logic_name ) = @_ ; 
  if ($logic_name) { 
    my $obj = Bio::EnsEMBL::Analysis->new ( -logic_name => $logic_name );  
    $self->{_new_logicname} = $obj ; 
  } 
  return $self->{_new_logicname} ; 
}





sub print_exon {
  my ($self, $ex,$string,$nr_tab) = @_ ; 
    if ($self->{v}) { 
      $string ="" unless $string ;
      my $tr_dbID = "" ;  
      my $hseq_name = "" ; 
      $nr_tab = 0  unless $nr_tab;
      my $tr = $ex->transcript() ;
      print $ex->biotype if $ex->biotype ; 
  
      my @sup_feat = @{ $ex->transcript->get_all_supporting_features} ; 
      $tr_dbID = "\ttr-db-id " . $ex->transcript->dbID ; 
      $hseq_name = $sup_feat[0]->hseqname ; 
      print  " " . $hseq_name . 
             "\t". $ex->seq_region_start .  "---". $ex->seq_region_end . 
             "\t". $ex->seq_region_end.  "---". $ex->seq_region_start . 
             "\t". $ex->seq_region_strand . 
             "\t". $ex->biotype ."$string" . 
             "\n" ; 
    }
}

sub print_intron {
  my ($self, $ex, $string ) = @_ ; 
   if ($self->{_verbose}) { 
     $string = "" unless $string ; 
     print "INTRON  " .  
       "\t". $ex->seq_region_start . 
       "---". $ex->seq_region_end . 
       "\t". $ex->seq_region_strand . 
       "\t". "$string" . 
       "\n" ; 
   } 
}

sub print_object {
  my ( $ex, $string  ) = @_ ; 
   print "$string :\n"  if $string ; 
   print "srs -- sre:" .  
       "\t". $ex->seq_region_start . 
       "---". $ex->seq_region_end . 
       "\t". $ex->seq_region_strand . 
       "\t". $ex->biotype  . 
#       "\t" . $ex .  
       "\n" ; 
}



=head2 

Arg : Arrayref to array of array of transcripts 
Function : removes redundant transcripts by comparing their exon-hashkeys 
Returns : Arrayref to Array of Transcripts 

=cut 

sub remove_redundant_transcripts {
  my ($tr_ref) = @_ ;

  my @non_red_trans ; 
  my %tmp ; 

  # make hashkey out of exon-start-end-strand for each transcripts to see if 
  # they are the same 

  for my $t (@$tr_ref){ 
    my @ex= @{ $t->get_all_Exons } ; 
    my $string ; 
    for my $e (@ex) { 
      my $exon_hk = $e->hashkey ; 
      $string.=$exon_hk ; 
    }
    
    push @{$tmp{$string}},  $t ; 
  }

  # since hashkeys are unique, get only one transcript for each hashkey
    
  for my $k (keys %tmp){
    my @t = @{$tmp{$k}} ; 
    push @non_red_trans, $t[0]; 
  }
  return \@non_red_trans ; 
}



=head2

Name : get_all_evidence_sets
Arg[1] : none
Function : Returns hash of evidence_set_names($key) and biotypes  
Returnval : Hashref $hash{evidence_set_name} = @\biotypes 

=cut 


sub get_all_evidence_sets { 
  my ($self) = @_ ; 
  return $self->{evidence_sets}; 
}



=head2

Name : get_genes_by_evidence_set($evidence_set)
Arg[1] : String 
Function : Returns all genes of an evidence set 
   (an evidence set contains genes of different biotypes ). 
   If there were any PredictionTranscripts specified as evidence set, 
   they were converted to genes and returned, too.

=cut


 
sub get_genes_by_evidence_set { 
  my ($self,$ev_set) = @_ ; 
  my @ev_set_genes ; 

  for my $biotype ( @{ $self->get_biotypes_of_evidence_set( $ev_set )} ) { 
    if ($self->get_genes_by_biotype($biotype)){ 
     push @ev_set_genes, @{ $self->get_genes_by_biotype( $biotype ) } ;  
    }
 }
  return \@ev_set_genes ; 
}
 




=head2

Name : get_genes_by_biotype($arg)  
Arg[1] : String 
Function : Returns all genes of a specific biotype (like all simgw_100 or genscan genes) 
Returnval : Arrayref of Bio::EnsEMBL::Gene objects

=cut 

sub get_genes_by_biotype {
  my ($self, $biotype ) = @_ ;
  return ${ $self->{all_genes_href}}{$biotype};  
}


=head2

Name : get_biotypes_of_evidence_set($arg)  
Arg[1] : String 
Function : Returns all biotypes specific evidence set (like simgw_100 , genscan ) 
Returnval : Arrayref of Strings 

=cut 


sub get_biotypes_of_evidence_set {
  my ($self, $ev_set ) = @_ ;
  return ${ $self->{evidence_sets}}{$ev_set};  
}



=head2 convert_to_genes( $tref ) 

Arg : Arrayref to Bio::EnsEMBL::Transcript objects
Name : convert_to_genes 
Function : converts all transcripts to Genes and set's new biotype as specified in the 
           Config 
Return : returns arrayref of Bio::EnsEMBL::Gene objects  

=cut 

sub convert_to_genes { 
  my ($self, $tref) = @_; 
  my @genes ; 

  for my $t (@$tref) {
    my $g = Bio::EnsEMBL::Gene->new() ;
    $g->biotype( $t->biotype ) ;
    $g->add_Transcript($t);
    $g->analysis($self->analysis) ; 
    push @genes, $g ; 
  }
  for my $g(@genes) { 
    throw("there are no transcripts for this gene\n") if scalar(@{$g->get_all_Transcripts}) == 0 ;  
    for my $tr ( @{$g->get_all_Transcripts} ) { 
      throw("there are no exons  for this transcript \n") if scalar(@{$tr->get_all_Exons}) == 0 ;  
    }
  } 
  return \@genes ;  
}



#
#sub all_exons_unvisited {
#  my ( $exon_cluster_aref) = @_ ;
#  for my $ec (@$exon_cluster_aref) {
#    for my $e ( @{ $ec->get_all_Exons_in_ExonCluster } ) {
#      return 0 if ( $e->visited ) ;
#    }
#  }
#  return 1 ;
#}
#


sub all_exons_visited {
  my ( $exon_cluster_aref) = @_ ;

  for my $ec (@$exon_cluster_aref) {
    for my $e ( @{ $ec->get_all_Exons_in_ExonCluster } ) {
      return 0 unless ( $e->visited ) ;
    }
  }
  return 1 ;
}

sub get_unvisited_exons_in_next_ec {
  my ( $exon_cluster_aref) = @_ ;
  my @unprocessed ;
  for my $ec (@$exon_cluster_aref) {
    for my $e ( @{ $ec->get_all_Exons_in_ExonCluster } ) {
       push @unprocessed, $e  unless ( $e->visited ) ;
     }
     return \@unprocessed if ( @unprocessed > 0 ) ;
  }
  return [] ;
}
























1; 

