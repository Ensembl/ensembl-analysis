=head1 LICENSE

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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::lincRNAEvaluator - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut
package Bio::EnsEMBL::Analysis::Runnable::lincRNAEvaluator;

use strict;  
use warnings;
use vars   qw(@ISA);

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils qw(transfer_supporting_evidence Exon_info);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(Transcript_info print_Transcript_and_Exons);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils; 
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(Gene_info get_multi_Exon_Genes get_single_Exon_Genes compute_6frame_translations ); 
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils; 
use Bio::EnsEMBL::Analysis::Tools::FilterBPlite;   
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Registry; 

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);



sub new{
  my ($class,@args) = @_;
  
  my $self = $class->SUPER::new(@args); 

  my($linc_rna_genes, $ensembl_genes ) = 
    rearrange([qw(linc_rna_genes ensembl_genes )], @args); 

  $self->linc_rna_genes($linc_rna_genes); 
  $self->ensembl_genes($ensembl_genes); 
  
  print "Found linc_rna: " . scalar(@$linc_rna_genes) . " \n"; 
  print "Found ensembl : " . scalar(@$ensembl_genes) . " \n";
  
  return $self;
} 



sub run{
  my ($self) = @_; 

  # NOTES ON cluster_Genes method:
 
  # cluster_Genes returns 2 types of clusters / 4 types of data :
  # 
  #  1) a genonmic region with more than one gene [ a set ] which can contain : 
  #  - both sets of genes [ set1 + set 2 ] , a 'twoway-cluster' 
  #  - only one set of genes [ set1 ] OR set2, a 'oneway-cluster'  
  # 
  #  2) unclustered genes : a cluster with only one gene, which does not cluster with anything - 
  #     not even with other genes of its own type. Each unclustered 'cluster' only contains ONE gene,
  #     which can be a set1 gene or a set2 gene.   
 


  # STAGE 1: Sort lincRNA candidate genes into two groups: those with protein domains and those without.
  #          The candidate genes have already been broken down into single-transcript genes
  #          when fetched by the RunnableDB. 
  #          Next, cluster genes in the two groups against each other.  
  #          The idea is to pick lincRNA candidates which (i) do not contain protein features
  #          themselves, and (ii) do not overlap with other candidates which do contain
  #          protein features.

  print "\n";
  print "Stage 1) Sort lincRNA candidate genes into two groups: those with protein domains and those without,";
  print " then cluster the groups against each other.\n";

  my ( $prot_feat_genes, $other ) = $self->sort_linc_rna_genes;  

  print "  Found " . scalar(@$other) . " lincRNA candidate genes without protein domains\n" ;    
  print "  Found " . scalar(@$prot_feat_genes) . " lincRNA candidate genes with protein domains\n\n" ;     

  print "  Stage1_set1_input: ". scalar(@$prot_feat_genes) . " genes with domains.\n";
  print "  Stage1_set2_input: ". scalar(@$other) . " genes without domains.\n\n";

  my %types_hash = %{ make_types_hash($prot_feat_genes, $other ,  'GENES_W_PROT_DOM','NO_PROT_DOM')} ; 
  my ($step1_clusters, $step1_unclustered) = cluster_Genes( [@$prot_feat_genes, @$other] , \%types_hash ) ;  

  my @no_prot_domain_genes = (  # lincRNA candidate genes which DO NOT cluster with any genes with a protein domain 
                                @{get_oneway_clustering_genes_of_set($step1_clusters,"NO_PROT_DOM")},  
                                @{get_oneway_clustering_genes_of_set($step1_unclustered,"NO_PROT_DOM")}
                               );  

  print "  Found ". scalar(@no_prot_domain_genes ) . " lincRNA candidate genes which do not contain protein domains ".
        " and do not cluster with any candidates which have protein domains.\n" ;    

  
  my @genes_with_prot_domain = ( # rejected lincRNA candidate genes which contain protein domains, or cluster with
                                 # protein-domain containing candidates 
                                  @{get_oneway_clustering_genes_of_set($step1_clusters,"GENES_W_PROT_DOM")},
                                  @{get_oneway_clustering_genes_of_set($step1_unclustered,"GENES_W_PROT_DOM")},
                                  @{get_twoway_clustering_genes_of_set($step1_clusters,"GENES_W_PROT_DOM")},
                                  @{get_twoway_clustering_genes_of_set($step1_clusters,'NO_PROT_DOM')},        
                               );

  print "  Found ". scalar(@genes_with_prot_domain) . " rejected lincRNA candidate genes which are somehow associated with protein domains.\n"; ##  cDNA DEBUG


  # STAGE 2: Cluster no-protein-domain lincRNA candidates from stage 1 vs Ensembl genes (the validation set).
  #          The idea is to find candidates which do not overlap with Ensembl genes
  #          (except for those with special gene biotypes, see stage 3) 
  #  
  #          Filtering of single-exon lincRNAs and artefactual two-exon lincRNAs (containing frameshift introns)
  #          happens at this stage if required.
  
  print "DEBUG::proble1\n";
  my $multi_exon_genes  = get_multi_Exon_Genes(\@no_prot_domain_genes) ;     
  print "DEBUG::proble2\n";
  my $single_exon_genes  =  get_single_Exon_Genes(\@no_prot_domain_genes) ;     

  print "  Up to this stage, out of all lincRNA candidates, " . scalar(@$multi_exon_genes) . " is multi-exon and " . scalar(@$single_exon_genes) . " is single-exon.\n";
  
  my $artefactual_2ex_gene_cnt = 0;
  my @multi_exon_genes_passed;

  # Each of the multi_exon_genes below is a single-transcript gene

  MULTI_EX_GENE: foreach my $gene(@$multi_exon_genes) {
    my $transcript = ${$gene->get_all_Transcripts()}[0];
    my $exon_cnt = scalar(@{$transcript->get_all_Exons}); 
    if ($exon_cnt == 2) {
      my $intron = ${$transcript->get_all_Introns()}[0];
      if ($intron->length < ($self->max_frameshift_intron_len) ) {
      	# print "DEBUG::lincRNA eval:: check length " . $self->max_frameshift_intron_len . "\n";
      	
        $artefactual_2ex_gene_cnt ++;
        next MULTI_EX_GENE;
      } else {
        push (@multi_exon_genes_passed, $gene);
        next MULTI_EX_GENE;
      }
    } else {  # more than 2 exons, automatic pass
      push (@multi_exon_genes_passed, $gene);
    }
  }

  my @filtered_no_prot_domain_genes;

  if ($self->exclude_single_exon_lincRNAs && $self->exclude_artefact_two_exon_lincRNAs) {
    print "  *** Excluding all single-exon and artefactual two-exon lincRNA candidates (contain frameshift introns)...\n"; 
    print "     " . scalar(@$single_exon_genes) . " single-exon lincRNA candidates discarded.\n" ;  
    print "     $artefactual_2ex_gene_cnt two-exon lincRNA candidates discarded.\n";
    print "     Only " . scalar(@multi_exon_genes_passed) . " multi-exon lincRNA candidates will proceed to stage 2 clustering.\n" ;  
    @filtered_no_prot_domain_genes = @multi_exon_genes_passed;  
  } elsif ($self->exclude_single_exon_lincRNAs && !$self->exclude_artefact_two_exon_lincRNAs) {
    print "  *** Excluding all single-exon lincRNA candidates...\n";
    print "     " . scalar(@$single_exon_genes) . " single-exon lincRNA candidates discarded.\n" ;
    print "     " . scalar(@$multi_exon_genes) . " multi-exon lincRNA candidates will proceed to stage 2 clustering.\n";
    @filtered_no_prot_domain_genes = @$multi_exon_genes;
  } elsif (!$self->exclude_single_exon_lincRNAs && $self->exclude_artefact_two_exon_lincRNAs) {
    print "  *** Excluding artefactual two-exon lincRNA candidates (contain frameshift introns)...\n";
    print "     $artefactual_2ex_gene_cnt two-exon lincRNA candidates discarded.\n";
    print "     " . scalar(@{$single_exon_genes}) . " single-exon and ". scalar(@multi_exon_genes_passed) . " multi-exon lincRNA candidates will proceed to stage 2 clustering.\n";
    @filtered_no_prot_domain_genes = (@{$single_exon_genes}, @multi_exon_genes_passed);
  } elsif (!$self->exclude_single_exon_lincRNAs && !$self->exclude_artefact_two_exon_lincRNAs) {
    print "     " . scalar(@{$single_exon_genes}) . " single-exon and ". scalar(@{$multi_exon_genes}) . " multi-exon lincRNA candidates will proceed to stage 2 clustering.\n";
    @filtered_no_prot_domain_genes = @no_prot_domain_genes;  # no modifications
  }
   
  print "\n"; 
  print "\nStage 2) Cluster lincRNA candidates (not associated with protein domains) vs Ensembl core genes...\n";
  print "  Stage2_set1_input: " . scalar(@{$self->ensembl_genes}) . " Ensembl core 'validation' genes.\n" ;    
  print "  Stage2_set2_input: " . scalar(@filtered_no_prot_domain_genes). " genes with NO protein_domain.\n" ;      

  my %types_hash2 = %{ make_types_hash($self->ensembl_genes, \@filtered_no_prot_domain_genes, 'ENSEMBL_SET', 'NO_PROT_DOM')} ;   

  my ($step2_clusters, $step2_unclustered) = cluster_Genes( [@{$self->ensembl_genes}, @filtered_no_prot_domain_genes] , \%types_hash2 ) ;  

  my @lincrna_unclustered = ( 
                                @{get_oneway_clustering_genes_of_set($step2_clusters,"NO_PROT_DOM")},  
                                @{get_oneway_clustering_genes_of_set($step2_unclustered,"NO_PROT_DOM")}
                            ); 
  print "\n";

  print "  Found " . scalar(@lincrna_unclustered). " lincRNAs which don't contain protein domains AND don't cluster with any Ensembl genes.\n" ;    
                                    
  my @clustered_lincrna  =  @{ get_twoway_clustering_genes_of_set($step2_clusters,'NO_PROT_DOM')};
  my @ensembl  =  @{ get_twoway_clustering_genes_of_set($step2_clusters,'ENSEMBL_SET')};

  print "  Found " . scalar(@clustered_lincrna). " lincRNAs which don't contain protein domains but cluster with Ensembl genes.\n" ;    
  print "  Found " . scalar(@ensembl). " Ensembl genes which cluster with lincRNA genes.\n" ;    


  # STAGE 3: For no-protein-domain lincRNA candidates which cluster with Ensembl genes in stage 2,
  #          look for special cases where the Ensembl genes involved are of "processed_transcript" or 
  #          "lincRNA" biotypes (the latter are existing lincRNAs already in the Ensembl DB).
  #
  #          lincRNA candidates which overlap with anything other than "processed_transcript" genes
  #          or existing "lincRNA" genes will be rejected.

  print "\nStage 3) Check stage 2 two-way clusters for cases where no-protein-domain lincRNAs overlap with Ensembl processed_transcript or existing lincRNA genes...\n";

  my @ncrna_clusters_with_processed_transcript;
  my @ncrna_clusters_with_existing_lincRNAs;
  my @ncrna_clusters_with_single_biotype;
  my @ncrna_clusters_with_multiple_biotypes; 

  # the tmp hashes help to keep track of proc_tran or lincRNA Ens genes which have already been matched
  # because one Ens gene can overlap with multiple single-transcript lincRNA "genes" in a given cluster.

  my %tmp_proc_trans_hash;
  my @proc_tran_genes_to_update; 
  my %tmp_old_lincRNA_hash;
  my @old_lincRNA_genes_to_update;

  CLUSTER: for my $cl ( @{ get_twoway_clusters($step2_clusters)} ) {  
    my @ens_genes = @{ $cl->get_Genes_by_Set('ENSEMBL_SET' )};    
    my @ncrnas_in_cluster = @{ $cl->get_Genes_by_Set('NO_PROT_DOM')};
      print "  Looking at a two_way cluster which contains " .scalar(@ens_genes)  . " Ensembl genes and " .
             scalar(@ncrnas_in_cluster)  . " lincRNA candidates without any protein domain \n" ; 

    my @ens_biotypes = @{ get_biotypes_of_ensembl_genes(\@ens_genes) }; 

    if ( scalar(@ens_biotypes) == 1 ) {
      if ( $ens_biotypes[0]=~m/^processed_transcript$/ ) {
        for my $ncrna_gene ( @ncrnas_in_cluster ) {  
          E_GENE: for my $ens ( @ens_genes ) {
            if ($tmp_proc_trans_hash{$ens}) {
              next E_GENE;
            } else {
              push (@proc_tran_genes_to_update, $ens); 
              $tmp_proc_trans_hash{$ens} = 1;
            }
          }
          push @ncrna_clusters_with_processed_transcript, $ncrna_gene; 
        }  
      } elsif ( $ens_biotypes[0]=~m/^lincRNA$/i) {
        for my $ncrna_gene ( @ncrnas_in_cluster ) {  
          L_GENE: for my $ens ( @ens_genes ) {  
            if ($tmp_old_lincRNA_hash{$ens}) {
              next L_GENE;
            } else {
              push (@old_lincRNA_genes_to_update, $ens);
              $tmp_old_lincRNA_hash{$ens} = 1;
            }
          }
          push (@ncrna_clusters_with_existing_lincRNAs, $ncrna_gene);  
        }
      } else {
        push (@ncrna_clusters_with_single_biotype, @ncrnas_in_cluster);
      }
    } elsif ( scalar(@ens_biotypes) > 1 ) {
      push (@ncrna_clusters_with_multiple_biotypes, @ncrnas_in_cluster) ; 
    } else {    
      throw("Gene clusters with Ensembl gene of NO biotype? This should not happen."); 
    }  
  }  # End foreach $cl

  foreach my $pt_gene(@proc_tran_genes_to_update) {
    print "    INFO : This processed_transcript gene clusters with lincRNA: ".  $pt_gene->dbID . "\t" . $pt_gene->analysis->logic_name . " " . 
           $pt_gene->seq_region_start . "\t" . $pt_gene->seq_region_end . "\n" ;  #  This will be printed out regardless of the config value of "PERFORM_UPDATES_ON_SOURCE_PROTEIN_CODING_DB"
  }
  
  foreach my $ol_gene(@old_lincRNA_genes_to_update) {
    print "    INFO : This existing lincRNA gene clusters with newly-identified lincRNA: ".  $ol_gene->dbID . "\t" . $ol_gene->analysis->logic_name . " " .
           $ol_gene->seq_region_start . "\t" . $ol_gene->seq_region_end . "\n" ;  #  This will be printed out regardless of the config value of "PERFORM_UPDATES_ON_SOURCE_PROTEIN_CODING_DB"
  }

  # Rejected data
  $self->ncrna_clusters_with_protein_domain(\@genes_with_prot_domain);  # rejected at stage 1 
  $self->ncrna_clusters_with_single_ens_biotype(\@ncrna_clusters_with_single_biotype);  # rejected at stage 3
  $self->ncrna_clusters_with_multiple_ens_biotypes(\@ncrna_clusters_with_multiple_biotypes);  # rejected at stage 3

  print "FINAL TALLY OF REJECTED lincNRAs:\n";
  print "  " . scalar(@{$self->ncrna_clusters_with_protein_domain}) . " rejected as they cluster with genes with protein domains\n" ;     
  print "  " . scalar(@{$self->ncrna_clusters_with_single_ens_biotype}). " rejected as they cluster with Ensembl gene of single biotype (other than proc_trans and lincRNA).\n";
  print "  " . scalar(@{$self->ncrna_clusters_with_multiple_ens_biotypes}). " rejected as they cluster with Ensembl genes of multiple biotypes (other than proc_trans and lincRNA).\n" ;     


  # data which will be used for output 
  $self->unclustered_ncrnas(\@lincrna_unclustered); 
  $self->ncrna_clusters_with_processed_transcript(\@ncrna_clusters_with_processed_transcript);   
  $self->ncrna_clusters_with_existing_lincRNAs(\@ncrna_clusters_with_existing_lincRNAs);   
  $self->proc_tran_genes_to_update(\@proc_tran_genes_to_update); 
  $self->old_lincRNA_genes_to_update(\@old_lincRNA_genes_to_update);


  
  # I now want to identify the ensembl gene which clusters with this gene.  
  # genes which have the same exons as ensembl will get biotype 'same as ensembl' 
  # @use= @{ compare_lincrna_vs_ensembl(\@ncrna_clusters_with_processed_transcript, \@ensembl) } ; 
  #

  #for ( @lincrna_unclustered ) {  
  #  $_->biotype("UN_clustered");  # DONT CHANGE BIOTYPE - used in RunnableDB  
  #} 

  #for ( @genes_with_prot_domain ) { 
  #  $_->biotype("clusters_with_prot_dom") ;  # DONT CHANGE BIOTYPE - used in RunnableDB  
  #} 
  #$self->output([@use ,@lincrna_unclustered,@genes_with_prot_domain,@not_use]);
}    



#  
# this routine sorts the potential lincrna genes into 2 subsets : 
# - genes which have a translation and a protein_domain annotation
# - genes which have no translation, or a translation and NO protein domain   
#  
# REMEMBER : you have to run protein_annotation before you run the lincRNAEvaluator
# 

sub sort_linc_rna_genes { 
  my ($self) = @_;  

  my @genes_to_sort = @{$self->linc_rna_genes } ;  
  ### print "DEBUG::lincRNAEval :: sort_linc_rna_genes:: Sorting " . @genes_to_sort . " genes by protein domains\n";  
  my ( @genes_w_pf, @no_pf ) ;  

  for my $g ( @genes_to_sort ) {   
    my @tr = @{$g->get_all_Transcripts}; 
    if ( @tr > 1 ) { 
      throw("Can only work on single-transcript genes\n") ; 
    } else { 
      my $t = ${$g->get_all_Transcripts}[0];    
      my $translation = $t->translation;
      ### print "DEBUG::lincRNAEval :: sort_linc_rna_genes:: tranlation_id of: " . $translation->display_id() . "has " . scalar(@{ $t->translation->get_all_ProteinFeatures }) . " Protein Features \n";

      my @pf ; 
      if ( $t->translation &&  scalar(@{ $t->translation->get_all_ProteinFeatures }) > 0 ) {  
      	### print "DEBUG::lincRNAEval :: sort_linc_rna_genes:: Sorting get all protein domains of the transcript and translation. They are: " . scalar(@{ $t->translation->get_all_ProteinFeatures })  . "\n" ;    
           $g->biotype("has_pf");
           push @genes_w_pf, $g ; 
       } else { 
          $g->biotype("no_pf");
         push @no_pf, $g ; 
       } 
     }  
  }
  ### print "DEBUG::lincRNAEval finish the sorting... \n "; 
  return ( \@genes_w_pf, \@no_pf ) ; 
} 


sub get_biotypes_of_ensembl_genes { 
   my ( $gene_aref ) = @_ ;    

  my %h ;   
  for my $g ( @$gene_aref ) {  
    if ( defined $g && ref($g)=~m/Gene/ ) {  
      $h{$g->biotype}++; 
    }
  }   
  return [keys %h];
} 


sub get_tr_hashkeys {  
   my ( $gene_aref ) = @_ ;    

   my %tr_hash_keys ; 

  for my $g ( @$gene_aref) { 
    my $hk ;   
    my @tr = @{ $g->get_all_Transcripts} ;   
    for my $t ( @tr ) { 
      my @ex = @{ $t->get_all_Exons };  
      for my $e ( @ex ) {  
         $hk.=$e->hashkey  ;
      } 
      push @{$tr_hash_keys{$hk}{gene}},$g;
      push @{$tr_hash_keys{$hk}{trans}},$t;  
      #if ( $t->stable_id) {  
      #  #### print $t->stable_id . "\t$hk\n" ;  
      #} 
    }    
  }  
  return \%tr_hash_keys ; 
} 

sub compare_lincrna_vs_ensembl {  
  my ( $lincrna, $ensembl ) = @_ ;  

  my %lincrna_tr_hk = %{ get_tr_hashkeys($lincrna)} ; 
  my  %ensembl= %{ get_tr_hashkeys($ensembl)} ; 

   my %trans_to_change; 

  for my $thk ( keys %lincrna_tr_hk ) {   
    if ( exists $ensembl{$thk} ) {  

      my @l_tr = @{ $lincrna_tr_hk{$thk}{'trans'}};
      my @l_g  = @{ $lincrna_tr_hk{$thk}{'gene'}};   

      my @e_tr = @{ $ensembl{$thk}{'trans'}};
      my @e_g  = @{ $ensembl{$thk}{'gene'}};  

      # concat stable_ids of genes with same hk 
      my $sids;
      for ( @e_tr ) {   
        $sids.=$_->stable_id;
      }
      for my $ltr ( @l_tr ) {  
        $trans_to_change{$ltr}=$sids; 
      }  
    } 
  }   
  # now change lincrna biotypes  
  for my $mlg ( @$lincrna ) {   
    for my $mltr ( @{ $mlg->get_all_Transcripts } ) { 
      if ( exists($trans_to_change{$mltr})) {   
        #push @same_as_ensembl, $mlg;
        print "changing biotype to same_as_ensembl for $trans_to_change{$mltr}\n";   
        $mltr->stable_id($trans_to_change{$mltr});
        $mltr->version("0"); 
        $mltr->biotype("same_as_ensembl");  # Don't change biotype as it's used in RunnableDB 
        $mlg->biotype("same_as_ensembl");   # Don't change biotype as it's used in RunnableDB
      }
    }
  }  
  return $lincrna;
}




#
#sub check_syntenic_regions_in_close_species { 
#   my ($self, $genes) = @_ ;  
# 
#   my @genes = sort { $a->seq_region_start <=> $b->seq_region_start } @$genes ; 
#
#   #### print scalar(@genes) . " genes found \n" ;
#
#   my $cmp_core = $self->cmp_core_dba(); 
#   my $cmp_cdna = $self->cmp_cdna_dba();
#   my $comp_species = $self->cmp_core_dba->get_MetaContainer->get_Species->binomial; 
#   my $compara = $self->cmp_compara_dba();
#   my %ncrna_protein_coding_overlap ; 
#
#
#   Bio::EnsEMBL::Registry->add_DBAdaptor("compara","compara", $compara); 
#   Bio::EnsEMBL::Registry->add_DBAdaptor($comp_species, ,"core", $cmp_core);
#
#   my $align_slice_adaptor = $compara->get_AlignSliceAdaptor(); 
#   my $mlssa = $compara->get_MethodLinkSpeciesSetAdaptor() ;
#   my $method_link_species_set = $mlssa->fetch_all_by_method_link_type("BLASTZ_NET") ;
#
#   my @result_genes ; 
#
#   GENE: for my $gene ( @genes ) {    
#
#     print_gene_info($gene);  
#     my $old_coord = create_old_coord_string_for_simple_feature($gene); 
#
#     my $reference_slice = $self->get_reference_slice($gene); 
#
#
#     SET: for my $set ( @$method_link_species_set ) {  
#        my $target_set = $self->mlss_name ; # "H.sap-M.mus blastz-net"
#
#        if (  $set->name=~m/$target_set/ ){  
#
#          my $align_slice = $align_slice_adaptor->fetch_by_Slice_MethodLinkSpeciesSet( $reference_slice, $set, "expanded",1); 
#
#          my $all_mus_aligning_slices = $align_slice->get_all_Slices($comp_species);  
#
#          print scalar(@$all_mus_aligning_slices) . " slices found for $comp_species...\n" ;    
#
#          if ( scalar(@$all_mus_aligning_slices) == 0 ) { 
#            print "no align slice found for $gene ...\n" ;    
#            push @result_genes, $gene ; 
#            next GENE ; 
#          }
#          #my ( $slices_with_mus_alignment, $no_mus_alignment ) = @{ $self->filter_align_slices($all_aligning_slices)} ; 
#          #print "XXX " . scalar(@$no_mus_alignment) . " slices don't have mouse alignment \n" ; 
#
#          THIS_ALIGN_SLICE: foreach my $this_align_slice (@$all_mus_aligning_slices) { 
#            my $species_name = $this_align_slice->genome_db->name(); 
#            my $exclude_species = $self->species_name; 
#
#            if ($species_name =~m/$exclude_species/  ) {   
#              next THIS_ALIGN_SLICE;  
#            } 
#
#            my @underlying_slices = @{ $this_align_slice->get_all_underlying_Slices() } ;
#             
#            # 3 sets : protein_coding genes on this_align_slice, 
#            #          cdnas on the slice of the protein_coding genes 
#            #          cdnas on UNDERLYING slices   
#            
#            my @all_mapped_prot_coding = @{$this_align_slice->get_all_Genes_by_type("protein_coding")};   
#
#            print  scalar(@all_mapped_prot_coding ) . " protein_coding genes found in alignment region $species_name\n\n"; 
#
#            my @cdna_all_underlying_slices = @{ $self->get_cdnas_on_underlying_slices(\@underlying_slices) } ;  
#            #my @cdnas_underlying_protein_coding = @{$self->get_cdnas_on_slices_of_mapped_protein_coding_genes(\@all_mapped_prot_coding,$species_name)  } ;   
#
#            # checking if there's an overlap with a protein_coding gene in mouse 
#            # best would be to implement some clusetering here ...  
#            if ( scalar(@all_mapped_prot_coding) ==  0  )  {
#              print "\nzero overlapping protein_coding in core of comparison_speces found.  ncrna ? \n" ; 
#              $ncrna_protein_coding_overlap{$gene}{$this_align_slice}{'gene'} = $gene ; 
#              $ncrna_protein_coding_overlap{$gene}{$this_align_slice}{'pc_overlap'} = scalar(@all_mapped_prot_coding) ; 
#            } else {  
#              $ncrna_protein_coding_overlap{$gene}{$this_align_slice}{'gene'} = $gene ; 
#              $ncrna_protein_coding_overlap{$gene}{$this_align_slice}{'pc_overlap'} = scalar(@all_mapped_prot_coding) ;  
#              print "adding " .join(",",@all_mapped_prot_coding). " to   $gene\n";
#              $ncrna_protein_coding_overlap{$gene}{$this_align_slice}{'pc_genes'} = \@all_mapped_prot_coding ;
#            } 
#            print "Shame.  $gene overlaps protein_coding genes in aligment region : " . scalar(@all_mapped_prot_coding) . "\n" ; 
#          }  
#
# 
#          # create simple features for the underlying slices - store ?    
#          
##          my %all_sf ; 
##          for my $us ( @underlying_slices_wo_gaps ) {  @ underlying_slices_wo_gaps can be computed in subroutine get_cdnas_on_underlying_slices
##            my @simple_features = @{$self->create_simple_features(\@underlying_slices,$old_coord)} ;
##            for ( @simple_features ) {
##              my $hk = $_->display_label . "===> ".$_->seq_region_start . "-" . $_->seq_region_end."-".$_->strand;
##              $all_sf{$hk}=$_;
##            } 
##          }
##          print " storing " . scalar(values %all_sf ) . " simple features...\n" ;
##          if ( scalar values %all_sf > 0 ) { 
##            my $sfa  =  $self->efg_out_dba->get_SimpleFeatureAdaptor() ; 
##            $sfa->store(values %all_sf) ; 
##          } 
#
#        }
#      }
#    } 
#
#  print "now processing genes which have mouse-synteny-regions and potential protei-coding overlap :\n"; 
#
#  for my $gene_ref( keys %ncrna_protein_coding_overlap ) {  
#     #print "hash output gene ref $gene_ref \n" ; 
#    
#     my %genes_on_ali_slice = %{ $ncrna_protein_coding_overlap{$gene_ref}};  
#
#     for my $tsi ( keys %genes_on_ali_slice ) {   
#     
#        print "BUMMER : " . $genes_on_ali_slice{$tsi}{'gene'} ."\n" ;  
#        print "BUMMER : " . $genes_on_ali_slice{$tsi}{'pc_overlap'} ."\n" ;     
#
#        if ( $genes_on_ali_slice{$tsi}{'pc_overlap'} > 0 ) {  
#          print  $genes_on_ali_slice{$tsi}{'gene'} . " overlaps protein-coding genes in mouse : "; 
#          for ( @{$genes_on_ali_slice{$tsi}{'pc_genes'}} ) {  
#              print $_->stable_id . ",";
#          } 
#          print  "\n";
#        } else { 
#          print  $genes_on_ali_slice{$tsi}{'gene'} . " does NOT overlap protein-coding genes in mouse : ";
#        } 
#
#        my $potential_ncrna =  $genes_on_ali_slice{$tsi}{'gene'}; 
#
#
#        #print ">" . $pot_tr."\n" ; 
#        #(my $ncrna_translation = $pot_tr->translate->seq() ) =~s/(.{1,60})/$1\n/g;    
#        #print $ncrna_translation ."\n" ; 
#        
#        my $potential_trans = ${$potential_ncrna->get_all_Transcripts}[0];  
#        my $pot_tr = compute_translations($potential_trans ) ; 
#
#       # my $run_pfam = 0 ; 
#
#        # we run pfam on the potential ncRNA's in second , extegrated step ( we run protein_annotation between lincRNAFinder.pm and lincRNAEvaluator.pm ) 
#        
#       # if ( $run_pfam ) { 
#          #my @pfam_domains = @{run_pfam_on_ncrna($pot_tr)} ; 
#          #if ( scalar(@pfam_domains ) == 0 ) {  
#          #  $potential_ncrna->biotype("no_pfam_domain") ; 
#          #} else { 
#            # this was to align my ncrna translation vs translation of mus-genes in syteny  region
#            #my @pc_genes = @{ $genes_on_ali_slice{$tsi}{pc_genes} } ;  
#            #for my $pcg ( @pc_genes ) {  
#            #    print ">" . $pcg->stable_id ."\n" ; 
#            #   my $mouse_gene = $self->cmp_core_dba->get_GeneAdaptor()->fetch_by_stable_id($pcg->stable_id) ;  
#            #  (my $canonical_translation = $mouse_gene->canonical_transcript->translation->seq ) =~s/(.{1,60})/$1\n/g;     
#            #  print $canonical_translation . "\n" ;  
#            #} 
#            #my $nr_pc_overlap = $genes_on_ali_slice{$tsi}{'pc_overlap'};  
#          #  $potential_ncrna->biotype("pc_overlap_pfam_domain");
#       #   }
#       # } else { 
#        #    $potential_ncrna->biotype("pc_overlap");
#       # }   
#       
#      #  my $make_6frame_translations=1;
#      #  if ( $make_6frame_translations ) {  
#          push @result_genes, compute_6frame_translations($potential_ncrna); 
#      #  }else { 
#      #    push @result_genes, $potential_ncrna; 
#      #  }  
#     } 
#  } 
#  return [\%ncrna_protein_coding_overlap, \@result_genes ] ; 
#}  
#
#sub run_pfam_on_ncrna{ 
#   my ( $pot_tr) = @_ ; 
# 
#          print "running pfam domain search on longest ORF of ncRNA ...\n" ;  
#          my $pfam_ana = Bio::EnsEMBL::Analysis->new( 
#              -logic_name => "pfam" , 
#              -program=> "hmmpfam", 
#              -db_file =>'/data/blastdb/Ensembl/Pfam_ls;/data/blastdb/Ensembl/Pfam_fs', 
#              -parameters => " --cut_ga  "
#          ) ; 
#
#          my $hmm = Bio::EnsEMBL::Analysis::Runnable::ProteinAnnotation::Pfam->new(  
#                   -analysis=>$pfam_ana, 
#                   -query => $pot_tr->translate, 
#                   -program => 'hmmpfam',
#                   -database => '/data/blastdb/Ensembl/Pfam_ls;/data/blastdb/Ensembl/Pfam_fs'
#                  ); 
#
#          $hmm->run(); 
#          my @pfam_domains = @{ $hmm->output()} ;  
#
#          print scalar(@pfam_domains ) . "pfam output found \n" ;  
#
#          for ( @pfam_domains) {  
#            print $_->start."\t".$_->end."\t".$_->hseqname."\t".$_->hstart."\t".$_->hend."\n" ; 
#          }  
#    return \@pfam_domains ; 
#}
#  
#
#sub create_simple_features {
#   my ( $self, $slices,$name ) = @_ ; 
#
#   my $adaptor = $self->cmp_cdna_dba();  
#
#   my @sf ;
#
#   my $analysis = Bio::EnsEMBL::Analysis->new() ;
#   $analysis->logic_name("ncrna_synt_region_human") ;
#
#   my $counter = 0 ;
#
#   for my $sl ( @$slices) {
#     my ($coordsys,$asm,$chr,$start,$end,$strand) = split/\:/,$sl->name ;
#
#     #print "TEST " . $sl->name . "\n" ; 
#     next if $sl->name=~m/GAP/ig;
#
#     my $slice = $adaptor->get_SliceAdaptor->fetch_by_region($coordsys,$chr) ;
#     #print $name . "   " . $slice->name . " $start  $end  \n" ; 
#     #for my $strand ( qw ( 1 -1 ) ) { 
#       $counter++;
#       my   $feature = Bio::EnsEMBL::SimpleFeature->new(
#        -start         => $start ,
#        -end           => $end,
#        -strand        => $strand,
#        -slice         => $slice,
#        -analysis      => $analysis,
#        -score         => 0,
#        -display_label => $name."_".$counter ,
#       );
#     #  print $name . "    " . $feature->start  . "    " . $feature->end . "\n" ; 
#       push @sf, $feature;
#     #}
#   }
#   return \@sf ;
#}
#

 






##
##sub find_short_translations { 
##   my ($self, $genes) = @_ ;  
##  
##   my @g = @$genes; 
##   my %result ; 
##
##   for my $g (@g) {   
##      # add translations to transcripts in all 6 reading frames 
##      for my $t ( @{$g->get_all_Transcripts} )  {   
##        my @translation_data  = @{ run_translate($t) }   ;    
##
##        if (  @translation_data > 0 ) { 
##          my $data = pop @translation_data ;         
##          my ( $l,$s,$e ) = @$data ;    
##          my $max_aa =  sprintf('%.0f',$t->seq->length/3) ;  
##          my $ratio = sprintf('%.0f', ( $l / $max_aa ) * 100) ;  
##          if ( $ratio < 30 ) { 
##             #print "translation length is quite short. flagging gene for blastx\n" ; 
##             #print "CHECK :  max = $max_aa    found_lengh = $l     $l / $max_aa =  $ratio \% \n" ;     
##             $result{$g}=$g;
##          }  
##        } else {  
##           print "no sensible translation found for " . $t->seq->seq . "\n" ; 
##           $result{$g}=$g;
##        } 
##      }   
##   }    
##   print scalar ( keys %result ) . " transcripts found with short translations < 30 % " ; 
##   return [values %result ]; 
##}
##
#sub blastx_results { 
#  my ( $self, $short_genes  )= @_ ; 
#  # OK noe lets check the result    
#  for my $g ( @{ $short_genes } ) {  
#     for my $tr ( @{ $g->get_all_Transcripts } ) {   
#
#      my $ana = Bio::EnsEMBL::Analysis->new ( 
#                                              -logic_name => "blastx_ncrna",   
#      ) ; 
#
#    my $parser = Bio::EnsEMBL::Analysis::Tools::BPliteWrapper->
#    new(
#        -regex => '^\w+\s+(\w+)',
#        -query_type => 'dna',
#        -database_type => 'pep',
#      );
#
#
#
#      my $blastx = Bio::EnsEMBL::Analysis::Runnable::Blast->new (  
#        -query => $tr->seq,
#        -database => 'ensembl_pep_r56' , 
#        -parser => $parser , 
#        -analysis => $ana, 
#        -program => "/usr/local/ensembl/bin/blastx",
#        );
#      $blastx->run() ;  
#       my $out = $blastx->output ;  
#       print " got " . @$out . " output \n" ;  
#       for ( @$out ) {  
#         print $_->hseqname ."\t"; 
#         print $_->percent_id."\t" ; 
#         print $_->score."\n" ;  
#       } 
#       print "\n\n" ; 
#    } 
#  } 
#}
# 
#
#sub separate_efg_features_by_logic_name {  
#  my ( $efg_genes ) = @_ ; 
#
#  my %logic_name_2_efg_features ; 
#  for my $e ( @$efg_genes ) {    
#      push @{$logic_name_2_efg_features{$e->analysis->logic_name}}, $e ;  
#  }    
#  return \%logic_name_2_efg_features ; 
#}

#
# if you got 2 genes with similar biotyeps in differnt sets it's quite useful to add the setname to them  
#
#sub get_all_genes_of_set_1_2 { 
#  my ( $self ) = @_;   
#
#  my @allgenes;  
#
#  for my $g ( @{$self->set_1_cdna_genes}){ 
#     $g->biotype($g->biotype."_set1") ;   
#    push @allgenes , $g ; 
#  } 
#  for my $g ( @{$self->set_2_prot_genes}){ 
#     $g->biotype($g->biotype."_set2") ; 
#     push @allgenes , $g ; 
#  }  
#  return \@allgenes ; 
#}
#



sub _make_types_hash { 
    my ( $self) = @_;  

    my (%types_hash,%tmp) ; 
    my @all_bt_1 = map { $_->biotype } @{ $self->set_1_cdna_genes} ;  
    @tmp{@all_bt_1}=1; 
    $types_hash{"set1"} = [ keys %tmp ]  ;  
    %tmp=(); 
    my @all_bt_2 = map { $_->biotype } @{ $self->set_2_prot_genes} ;  
    @tmp{@all_bt_2}=1; 
    $types_hash{"set2"} = [ keys %tmp ] ;   
    return \%types_hash ; 
}



sub input_genes{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'input_genes'} = $arg;
  }
  return $self->{'input_genes'};
}


sub output_biotype{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'output_biotype'} = $arg;
  }
  return $self->{'output_biotype'};
} 

sub cmp_core_dba{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'cmp_core_dba'} = $arg;
  }
  return $self->{'cmp_core_dba'};
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



sub this_core_db {
  my ($self, $arg) = @_;
  if($arg){
    $self->{'this_core_db'} = $arg;
  }
  return $self->{'this_core_db'};
} 

sub cmp_cdna_dba{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'cmp_cdna_dba'} = $arg;
  }
  return $self->{'cmp_cdna_dba'};
} 

sub efg_out_dba {
  my ($self, $arg) = @_;
  if($arg){
    $self->{'efg_out_dba'} = $arg;
  }
  return $self->{'efg_out_dba'};
}


sub cmp_compara_dba{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'cmp_compara_dba'} = $arg;
  }
  return $self->{'cmp_compara_dba'};
}

sub set_1_cdna_genes{
  my ($self, $arg) = @_;
  if($arg){
    for my $g  ( @$arg ) { 
      $g->biotype($g->biotype."_set1") ;   
    } 
    $self->{'set_1_genes'} = $arg; 
  }
  return $self->{'set_1_genes'};
}


sub set_2_prot_genes{
  my ($self, $arg) = @_;
  if($arg){
    for my $g  ( @$arg ) { 
      $g->biotype($g->biotype."_set1") ;   
    } 
    $self->{'set_2_genes'} = $arg;
  }
  return $self->{'set_2_genes'};
}  

 
sub efg_simple_feature_genes{
  my ($self, $arg) = @_;
  # don't think we need this fun now
  print "BK::lincRNAEvaluator::efg_simple_feature_genes:: with arg: $arg \n";

  if($arg){
    $self->{'efg_simple_feature_genes'} = $arg;
  }
  return $self->{'efg_simple_feature_genes'};
}

sub extend_efg_features{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'extend_efg_features'} = $arg;
  }
  return $self->{'extend_efg_features'};
}
 


sub filter_cdna_genes_by_exon_count {   
  # don't understand why this is here, since it is used in finder
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

sub print_gene_info { 
   my ($gene) = @_ ; 

   my $db_id = ''; 
   if ( $gene->dbID ) { 
     $db_id = $gene->dbID ; 
   }   
   my $sid='';
   if ( $gene->stable_id) { 
    $sid= $gene->stable_id; 
   }   
   print "\n\nprocessing potential ncRNA : " . $db_id . " [ " . $gene->biotype . " ] ".$sid. "\t".@{ $gene->get_all_Transcripts} . "exons\n"   ;   
}  


sub create_old_coord_string_for_simple_feature { 
  my ($gene) = @_; 

  my $first_tr = ${$gene->get_all_Transcripts}[0]; 
  my @tsf = @{$first_tr->get_all_supporting_features} ; 
  
  my $sf_name ='no_name';
  if ( @tsf > 0 ) { 
    $sf_name = ${${$gene->get_all_Transcripts}[0]->get_all_supporting_features}[0]->hseqname;
  } 
  my $old_coord = $sf_name."_".$gene->seq_region_name."_".$gene->seq_region_start."_".$gene->seq_region_end."_".$gene->dbID ;
  return $old_coord; 
} 


sub get_reference_slice { 
     my ($self, $gene) = @_ ;  

     my $sa = $self->this_core_db->get_SliceAdaptor;   

     my $reference_slice = $sa->fetch_by_region(
                                                $gene->slice->coord_system->name,
                                                $gene->slice->seq_region_name,
                                                $gene->seq_region_start,
                                                $gene->seq_region_end,
                                                );
      unless ( $reference_slice ) {
        throw("could not fetch reference slice.\n" ) ;
      } 
     return $reference_slice; 
} 




sub get_cdnas_on_underlying_slices { 
   my ($self,$underlying_slices, ) = @_;  

   my %underlying_cdna ;
   my @underlying_slices_wo_gaps ; 
 
   for my $us ( @$underlying_slices ) {
      next if $us->coord_system_name=~m/GAP/i;
       push @underlying_slices_wo_gaps, $us ; 

       print "underlying aligning slice : " . $us->name . "\n" ; 
       my $ms = $self->cmp_core_dba->get_SliceAdaptor->fetch_by_name($us->name) ; 
       my @cdna = @{  $self->cmp_cdna_dba->get_GeneAdaptor->fetch_all_by_Slice($ms)} ;   

       for my $c ( @cdna ) {    
          my $hk;
          for my $tc ( @{ $c->get_all_Transcripts} ) { 
            for my $ec ( @{ $tc->get_all_Exons } ) {   
              $hk.=$ec->hashkey()
            }  
          } 
          #print "$hk " . ${${$c->get_all_Transcripts}[0]->get_all_supporting_features}[0]->hseqname."\n";
          $underlying_cdna{$hk} = $c; 
       }  
           #$underlying_cdna{${${$_->get_all_Transcripts}[0]->get_all_supporting_features}[0]->hseqname}=$_;
   }  

    print scalar(keys %underlying_cdna) . "  cDNAs found on the underlying slices in alignment region.\n";
    for ( values %underlying_cdna ) { 
       print "CDNA: " . $_->biotype . "  ".$_->seq_region_name . "   " . $_->seq_region_start . "  " . $_->seq_region_end . "  " ;  
       print ${${$_->get_all_Transcripts}[0]->get_all_supporting_features}[0]->hseqname  ."\n" ; 
    }    
   return [ values %underlying_cdna ];
}


sub get_cdnas_on_slices_of_mapped_protein_coding_genes { 
   my ( $self , $all_mapped_prot_coding,$species_name) = @_;  

   my $cmp_core = $self->cmp_core_dba(); 
   my $cmp_ga = $cmp_core->get_GeneAdaptor() ;  
   my $cmp_cdna = $self->cmp_cdna_dba();

   my %unique_cdnas_overlapping_protein_coding ; 
 
   for my $mg ( @$all_mapped_prot_coding ) { 
      # get absolute coordinates. ( $mg has relavtive start /end on slice ... )  
      
       my $mouse_gene = $cmp_ga->fetch_by_stable_id($mg->stable_id) ; 

       print "\n\nFOUND:" . $mg->stable_id."\t".$mouse_gene->seq_region_name."\t".$mouse_gene->seq_region_start ;   


       print $mouse_gene->seq_region_end."\t$species_name\t".$mg->biotype."\n" ;  
 
       my $mouse_trimmed_slice = $cmp_core->get_SliceAdaptor->fetch_by_region(
       $mouse_gene->slice->coord_system_name(), 
       $mouse_gene->slice->seq_region_name, 
       $mouse_gene->seq_region_start,
       $mouse_gene->seq_region_end
       ); 
 
        my @cdna = @{  $cmp_cdna->get_GeneAdaptor->fetch_all_by_Slice($mouse_trimmed_slice)} ; 
        print scalar(@cdna) . " cdnas in region on mouse_trimmed_slice\n" ;  

        # unique by location 
        for my $c ( @cdna ) {    
          my $hk;
          for my $tc ( @{ $c->get_all_Transcripts} ) { 
            for my $ec ( @{ $tc->get_all_Exons } ) {   
              $hk.=$ec->hashkey()
            }  
          } 
          #print "$hk " . ${${$c->get_all_Transcripts}[0]->get_all_supporting_features}[0]->hseqname."\n";
          $unique_cdnas_overlapping_protein_coding{$hk} = $c; 
       }  

        for ( @cdna ) {
           print "\n" . $_->biotype . "  ".$_->seq_region_name . "   " . $_->seq_region_start . "  " . $_->seq_region_end . " " ;  
           print ${${$_->get_all_Transcripts}[0]->get_all_supporting_features}[0]->hseqname   ; 
        }
       print "\n" ; 
   } 
   return [ values %unique_cdnas_overlapping_protein_coding ];
} 


sub filter_align_slices { 
  my ( $self, $all_aligning_slices) = @_ ; 

  my (@align_in_mouse, @no_align_in_mouse ) ; 
  my %tmp; 

  ALIGN_SLICE: foreach my $this_align_slice (@$all_aligning_slices) {  

     my $species_name = $this_align_slice->genome_db->name(); 
     my $exclude_species = $self->species_name;
     if ($species_name =~m/$exclude_species/  ) {   
        print "species are the same : $species_name - $exclude_species\n" ; 
        next ALIGN_SLICE;  
     }
     push @align_in_mouse,  $this_align_slice; 
     $tmp{$this_align_slice}=1;
   }  
   for my $as ( @$all_aligning_slices ) {  
     unless ( exists $tmp{$as} ) {  
       push @no_align_in_mouse, $as; 
     } 
   } 
   return [\@align_in_mouse, \@no_align_in_mouse] ; 
}


sub ncrna_clusters_with_multiple_ens_biotypes { 
  my ( $self, $genes ) = @_; 

  if ( defined $genes ) { 
    for my $g ( @$genes ) {   
      $g->biotype($g->biotype."_reject_mult");
      for my $t( @{$g->get_all_Transcripts } ) { 
        $t->biotype($g->biotype);
      } 
    } 
    $self->{ncrna_clusters_with_multiple_ens_biotypes}=$genes ;
  }
  return $self->{ncrna_clusters_with_multiple_ens_biotypes};  
}   


sub ncrna_clusters_with_single_ens_biotype {
  my ( $self, $genes ) = @_;

  if ( defined $genes ) {
    for my $g ( @$genes ) {
      $g->biotype($g->biotype."_reject_single");
      for my $t( @{$g->get_all_Transcripts } ) {
        $t->biotype($g->biotype);
      }
    }
    $self->{ncrna_clusters_with_single_ens_biotype}=$genes ;
  }
  return $self->{ncrna_clusters_with_single_ens_biotype};
}



sub ncrna_clusters_with_protein_domain { 
  my ( $self, $genes ) = @_;  

  if ( defined $genes ) { 
    for my $g ( @$genes ) {  
      $g->biotype($g->biotype."_reject_prot_dom");
      for my $t( @{$g->get_all_Transcripts } ) { 
        $t->biotype($g->biotype);
      } 
    } 
    $self->{ncrna_clusters_with_protein_domain}=$genes;
  }
  return $self->{ncrna_clusters_with_protein_domain};  
} 



use vars '$AUTOLOAD';
sub AUTOLOAD {
 my ($self,$val) = @_;
 (my $routine_name=$AUTOLOAD)=~s/.*:://; #trim package name
 $self->{$routine_name}=$val if $val ;
 return $self->{$routine_name} ;
}
sub DESTROY {} # required due to AUTOLOAD

1;  
