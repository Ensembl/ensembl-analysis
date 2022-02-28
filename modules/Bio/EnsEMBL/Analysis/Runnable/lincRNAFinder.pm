=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

Main goal: get RNAseq models and separate them into those which overlap with protein_coding genes and those which don't. Report those that they don't overlap

=head1 DESCRIPTION


=head1 METHODS

=cut
package Bio::EnsEMBL::Analysis::Runnable::lincRNAFinder;

use strict;  
use warnings;

use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils; 
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw(compute_6frame_translations Gene_info) ; 

use parent ('Bio::EnsEMBL::Analysis::Runnable');


sub run{
  my ($self) = @_; 

  # sorting genes in sets 1 + 2 according to biotype + config file settings + exon count 
  my ($single_exon_cdna , $multi_exon_cdna_genes ) = $self->filter_cdna_genes_by_exon_count() ; 
  my $types_hash = make_types_hash($multi_exon_cdna_genes, $self->set_2_prot_genes, 'SET_1_CDNA','SET_2_PROT');

  ### Stage 1 --- separate cDNAs into those which overlap with protein_coding genes and those which don't -------
  print "\nStage 1a) cluster MULTI-exon cDNAs vs protein_coding genes at all exons ";
  if ($self->ignore_strand) {
    print "on either strand (i.e. strandedness ignored!)...\n";
  }
  else {
    print "on the same strand...\n";
  }
  print "GOING TO COMPARE:: " . scalar(@{$multi_exon_cdna_genes}) . " (number of multi_exon_cdna)  VS.  " . scalar(@{$self->set_2_prot_genes}) ."  (number of proteins) \n";

  my ($step1_clusters, $step1_unclustered) = cluster_Genes([@{$multi_exon_cdna_genes}, @{$self->set_2_prot_genes}], $types_hash, 0, $self->ignore_strand);
  my @unclustered_cdna_genes = (
                                 @{get_oneway_clustering_genes_of_set($step1_clusters,"SET_1_CDNA")},
                                 @{get_oneway_clustering_genes_of_set($step1_unclustered,"SET_1_CDNA")},
                               );

  # $step1_clusters has the clustered, means no lincRNA candidates. 
  # $step1_unclustered has the lincRNA candidates 
 
  if ($self->find_single_exon_candidates == 1) {
    print "Stage 1b) cluster SINGLE-exon cDNAs vs protein_coding genes (strandedness of models always IGNORED) ...\n"; 
    my ($single_exon_clustered, $single_exon_unclustered) = cluster_Genes([@$single_exon_cdna, @{$self->set_2_prot_genes}], $types_hash , 1 , 1 ) ;   # 4th arg = ignore_strand
    push(@unclustered_cdna_genes, @{get_oneway_clustering_genes_of_set($single_exon_clustered,"SET_1_CDNA")} );
    push(@unclustered_cdna_genes, @{get_oneway_clustering_genes_of_set($single_exon_unclustered,"SET_1_CDNA")} );
    # push @cdna_gene_clusters_with_pc ,@{ get_twoway_clustering_genes_of_set($single_exon_clustered,"SET_1_CDNA") } ;
  } 

  # $self->update_and_copy_cdna(\@cdna_gene_clusters_with_pc,"cdna_update_protein_coding" );   # cDNAA DEBUG. Select cdnas that cluster with protein coding - Don't need them as lincRNAs! 

  # Stage 2
  # -------
  print scalar(@unclustered_cdna_genes) . " genes including single exons genes\n";

  print "\nStage 2) check for 6-frame translations in lincRNA candidates\n";
  print "result_set: " . scalar(@{$self->result_set}) . " and unclustered set: " . scalar(@unclustered_cdna_genes) . " for 3 step \n";

  my @genes_with_translations ;  
  my @genes_withOUT_translations ; 
  RG: for my $rg( @unclustered_cdna_genes ) {
 	# WARNING: when stand specific models this should be only against 3 frames of correct strand... 
    my $new_gene = compute_6frame_translations($rg);  # compute_translation() # 
    $new_gene->biotype("pre_finder_round1");  
    # print " translations found for gene " . Gene_info($rg) . "::" . $rg->display_id() . "\n"; 
    print scalar(@{ $rg->get_all_Transcripts} ) ." translations found for old gene " . Gene_info($rg) . "::" . $rg->display_id() . "\n";  # " seq_region: " . $rg->seq_region_name . " start: " . $rg->seq_region_start . " end: " . $rg->seq_region_end . " strand: " . $rg->seq_region_strand . " \n" ;

    if (!defined $new_gene->get_all_Transcripts) {
      $self->warning('  Could not compute translation for cDNA: gene dbID '. $rg->dbID . ' ' . $rg->seq_region_name . ' ' .
             $rg->seq_region_start . ' ' . $rg->seq_region_end.' '.$rg->length);
      next RG;
      $rg->biotype('gene_WITHOUT_translation');
      push @genes_withOUT_translations, $rg ; 
    }     
    print scalar(@{ $new_gene->get_all_Transcripts} ) ." translations found for new gene " . Gene_info($new_gene) . "::" . $new_gene->display_id() . "\n";  # " seq_region: " . $rg->seq_region_name . " start: " . $rg->seq_region_start . " end: " . $rg->seq_region_end . " strand: " . $rg->seq_region_strand . " \n";

    push @genes_with_translations, $new_gene ; 
  }
  print "--  ".scalar(@genes_with_translations) . " genes with 6-frame-translations found\n" ;     
  print "--  ".scalar(@genes_withOUT_translations) . " genes WITHOUT 6-frame-translations found\n" ;  # I migth need to report them too!!! 

  # cap the number of transcripts per gene according to config 
  my $capped_genes = $self->cap_number_of_translations_per_gene(\@genes_with_translations) ;  
  print "--  ".scalar(@$capped_genes) . " after genes with 6-frame-translations found\n" ;     
 
  my ($short, $long ) = $self->filter_genes_with_long_translations($capped_genes); 
  $capped_genes = ""; # WARNING: I am keeping too many @ and objects. This is a test!! Probably I can delete this... 
  print "--  ".scalar(@$long)." genes have long translations, discarded from lincRNA pipeline.\n";
  print "--  ".scalar(@$short)." genes with short translations found. They will be our lincRNAFinder final output.\n\n";

  # $self->update_and_copy_cdna($long, 'has_long_translation');   # cDNA DEBUG
  $self->output($short); 
}


sub print_output {
	my ( $self) = @_ ; 
}

sub cap_number_of_translations_per_gene { 
  my ( $self, $genes_with_6f_translations ) = @_ ;   
  my @capped_longest_genes;
  GENES: for my $g ( @$genes_with_6f_translations ) {
    my $transcripts = $g->get_all_Transcripts;
    if (defined $self->max_translations_stored_per_gene and @$transcripts > $self->max_translations_stored_per_gene) {
      $g->flush_Transcripts;
      my @sorted_t = sort {$a->translate->length <=> $a->translate->length} @$transcripts;
      for my $lt (splice @sorted_t, 0, $self->max_translations_stored_per_gene) {
        $g->add_Transcript($lt);
      }
    }
    push(@capped_longest_genes, $g);
  }   
  return \@capped_longest_genes; 
} 


sub filter_genes_with_long_translations {  
  my ( $self, $capped_longest_genes ) = @_ ; 

  my $max_trans_length_ratio;
   if ( defined $self->maximum_translation_length_ratio ){  
     if ( $self->maximum_translation_length_ratio < 101 ) {  
       $max_trans_length_ratio =  $self->maximum_translation_length_ratio; 
       print "-- your MAXIMUM_TRANSLATION_LENGTH_RATIO is: " . $max_trans_length_ratio . "\n";
     } else {  
       $self->throw("translation-length-to-transcript length ratio > 100 does not make sense.\n"); 
     } 
   }else { 
     $max_trans_length_ratio = 100 ; 
   }   

  # if a gene has a transcript with a translation over a specified % threshold transcript
  # length, the gene will not be taken as lincRNA candidate.
  
  # if WRITE_DEBUG_OPTION is turned on, the cDNA "genes" with long translations could 
  # be written to the DEBUG_OUTPUT_DB with new biotype "long translation".
  
  my (@long_genes, @short_genes );  

  GENES: for my $g ( @$capped_longest_genes ) { 
    for my $t ( @{$g->get_all_Transcripts } ) {   

       my $tl_length = $t->translate->length * 3 ; 
       my $tr_length = $t->length ;    
       my $ratio = sprintf('%.1f',$tl_length*100 /  $tr_length);
       if ( $ratio > $max_trans_length_ratio ) { 
         # The cDNA "genes" have no usable display IDs, hence not printing any in the next line:
         print "LONG_TRANSLATION: " .  Gene_info($g)   .  "  translation-length to transcript-length ratio ($ratio) higher than max. value allowed in config ($max_trans_length_ratio).\n";
         my $analysis = Bio::EnsEMBL::Analysis->new(
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
      $self->throw('Must pass RunnableDB:result_set an array ref not a '.$result_set);
    }
    push(@{$self->{'result_set'}}, @$result_set);
  }
  return $self->{'result_set'};
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

# it is not nice, but it is working 
use vars '$AUTOLOAD';
sub AUTOLOAD {
 my ($self,$val) = @_;
 (my $routine_name=$AUTOLOAD)=~s/.*:://; #trim package name
 if ( defined $val ) { 
   $self->{$routine_name}=$val;
   print "routine name for AUTOLOAD:: $routine_name \n"; 
}
return $self->{$routine_name} ;
}
sub DESTROY {} # required due to AUTOLOAD

1;  
