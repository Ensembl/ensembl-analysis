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

Bio::EnsEMBL::Analysis::Runnable::TranscriptConsensus -

=head1 SYNOPSIS

  my $runnable =
    Bio::EnsEMBL::Analysis::Runnable::TranscriptConsensus->new(
    -query     => $query,
    -analysis  => $analysis,
    -all_genes => \%biotypes_to_genes,
    # ref to $hash{biotype_of_gene} = @all_genes_of_this_biotype
    -evidence_sets => $evidence_sets,
    -dnadb         => $db, );

=head1 DESCRIPTION

TranscriptConsensus is an extension of TranscriptCoalescer that combines
protein coding and EST transcripts to identify the best supported
transcript models.  The initial gene sets are clustered, then collapsed
into a non-redundant set of exons and introns which are assigned
scores according to the amount of supporting evidence they have.
The similarity genes have UTR added using the est transcripts and
the resulting models are assigned a score by summing the individual
exon and intron scores.  The transcripts are sorted by score and the
highest scoring models are made into gene objets and written to the
TranscriptCoalescer database.

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Runnable::TranscriptConsensus;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(throw warning info);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils
  qw( count_non_canonical_splice_sites
      are_phases_consistent Transcript_info);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::EvidenceUtils
  qw( clone_Evidence );
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils
  qw( clone_Translation compute_translation);
use Bio::EnsEMBL::Analysis::Config::GeneBuild::TranscriptConsensus;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::CollapsedCluster;

use Bio::EnsEMBL::Analysis::Runnable;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);




=head2 new

  Function  : creates TranscriptConsensus-object
  Returnval : returns TranscripotConsensus-object

=cut




sub new {
  my ($class,@args) = @_ ;
  my $self = $class->SUPER::new(@args) ;

  my ( $all_genes, $evidence_sets, $filter_singletons, $filter_non_consensus, $filter_ests,
       $add_utr, $min_consensus, $utr_penalty, $end_exon_penalty, $est_overlap_penalty,
       $short_intron_penalty, $short_exon_penalty, $good_percent, $good_biotype, $small_biotype,
       $bad_biotype, $rnaseq_introns, $verbose) =
     rearrange([qw(
               ALL_GENES
               EVIDENCE_SETS
               FILTER_SINGLETONS
               FILTER_NON_CONSENSUS
               FILTER_ESTS
               ADD_UTR
               MIN_CONSENSUS
               UTR_PENALTY
               END_EXON_PENALTY
               EST_OVERLAP_PENALTY
               SHORT_INTRON_PENALTY
               SHORT_EXON_PENALTY
               GOOD_PERCENT
               GOOD_BIOTYPE
               SMALL_BIOTYPE
               BAD_BIOTYPE
               RNASEQ_INTRONS
               VERBOSE
              )], @args);

  $self->{evidence_sets} = $evidence_sets ;
  $self->{all_genes_href} = $all_genes ;
  $self->{filter_singletons} = $filter_singletons ;
  $self->{filter_non_consensus} = $filter_non_consensus;
  $self->{filter_ests} = $filter_ests ;
  $self->{add_utr} = $add_utr ;
  $self->{min_consensus} = $min_consensus ;
  $self->{utr_penalty} = $utr_penalty ;
  $self->{end_exon_penalty} = $end_exon_penalty ;
  $self->{est_overlap_penalty} = $est_overlap_penalty ;
  $self->{short_intron_penalty} = $short_intron_penalty ;
  $self->{short_exon_penalty} = $short_exon_penalty ;
  $self->{good_percent} = $good_percent ;
  $self->{good_biotype} = $good_biotype ;
  $self->{small_biotype} = $small_biotype ;
  $self->{bad_biotype} = $bad_biotype ;
  $self->{rnaseq_introns} = $rnaseq_introns ;
  $self->{verbose} = $verbose ; # verbose or not

  return $self ;
}



=head2 run

   Name       : run
   Args       : none
   Function   : Pre filters genes and clusters them, then collapses each cluster to make non-redundant
              : sets of exons and introns, builds transcripts and adds UTR and assigns them scores,
              : makes genes out of the highest scoring transcripts
   Exceptions : none
   Returnval  : none

=cut

sub run {
  my ($self) = @_;

  my $filter_singletons    = $self->{filter_singletons};
  my $filter_non_consensus = $self->{filter_non_consensus};

  my $verbose = $self->{verbose};
  my $count;

  # first we want to cluster on protein coding genes only then cluster
  # them with ESTs afterwards
  my @coding = @{ $self->get_genes_by_evidence_set('simgw') };
  my @non_coding =
    sort { $a->start() <=> $b->start() }
    @{ $self->get_genes_by_evidence_set('est') };

  # cluster the coding genes first
  my ( $coding_clusters, $coding_non_clusters ) =
    cluster_Genes( \@coding, $self->get_all_evidence_sets() );

  if ( scalar(@$coding_non_clusters) > 0 ) {
    push( @$coding_clusters, @$coding_non_clusters );
  }

  foreach my $coding_cluster ( @{$coding_clusters} ) {
    my @allgenes = @{ $coding_cluster->get_Genes() };

    # add overlapping non_coding genes
    if ($verbose) {
      print( "\n***\n# add overlapping non_coding genes\n" .
             "CLUSTER " . $coding_cluster->start() . " " .
             $coding_cluster->end() . "\n" );
    }

    foreach my $ncgene (@non_coding) {
      if ( $ncgene->end() < $coding_cluster->start() )      { next }
      if ( $ncgene->strand() != $coding_cluster->strand() ) { next }

      push( @allgenes, $ncgene );

      if ($verbose) {
        print( "Adding EST " . $ncgene->start() . " " . $ncgene->end() .
               " " . $ncgene->analysis->logic_name() . "\n" );
      }

      if ( $ncgene->start() > $coding_cluster->end() ) { last }
    }

    if ( $filter_singletons or $filter_non_consensus ) {
      my @tmp;
      my ( $clusters, $non_clusters ) =
        cluster_Genes( \@allgenes, $self->get_all_evidence_sets() );

      foreach my $cluster ( @$clusters, @$non_clusters ) {
        my $genes = $cluster->get_Genes();
        $genes = @{ $self->filter_genes(@$genes) };
        push( @tmp, $genes );
      }
      @allgenes = @tmp;
    }

    my ( $clusters, $non_clusters ) =
      cluster_Genes( \@allgenes, $self->get_all_evidence_sets() );

    # Create a hash of genes by strand to use when looking for
    # overlapping genes.  Useful to look at ALL genes not just those
    # in the cluster to prevent cluster joining.
    my $genes_by_strand;
    foreach my $gene (@allgenes) {
      push( @{ $genes_by_strand->{ $gene->strand() } }, $gene );
    }

    if ( scalar(@$non_clusters) > 0 ) {
      push( @$clusters, @$non_clusters );
    }

    $count = 0;

    my @leftover_genes;

    foreach my $cluster (@$clusters) {
      # cluster has to contain at least one similarity gene to be worth
      # continuing with
      if ( !@{ $cluster->get_Genes_by_Set('simgw') } ) { next }

      if ($verbose) {
        print("\nCluster $count\n");
        print( $cluster->start() . " " . $cluster->end() . " " .
               $cluster->strand() . "\n" );
      }

      $count++;

      # collapse the cluster down to make a non redundat set of introns
      # and exons with scores
      my $collapsed_cluster =
        $self->collapse_cluster( $cluster, $genes_by_strand );

      # add UTR to the similarity transcripts and make scores for the
      # transcripts
      my $transcripts =
        $self->make_transcripts( $cluster, $collapsed_cluster );

      # sort the scored transcripts by score and make gene objects out
      # of the best ones
      $self->make_genes( $transcripts, $cluster, $collapsed_cluster );

    } ## end foreach my $cluster (@$clusters)
  } ## end foreach my $coding_cluster ...

} ## end sub run

=head2 collapse_cluster

   Name       : collapse_cluster
   Arg[0]     : Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster
   Arg[1]     : Array ref of Bio::EnsEMBL::Gene objects
   Function   : Cretaes a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::CollapsedCluster object and
              : populates it with exons and introns from the GeneCluster
   Exceptions : none
   Returnval  : Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::CollapsedCluster

=cut

sub collapse_cluster{
  my ($self,$cluster,$genes_by_strand) = @_ ;
  my $genes = $cluster->get_Genes;
  my @exon_clusters = @{$cluster->get_exon_clustering_from_gene_cluster};
  my $collapsed_cluster = Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::CollapsedCluster->new();

  foreach my $gene (@$genes){
    foreach my $trans (@{$gene->get_all_Transcripts}){
      my @exons = sort {$a->start <=> $b->start } @{$trans->get_all_Exons};
      for (my $i =0 ; $i < scalar(@exons) ; $i++){
    my $exon_cluster;
    my $exon = $exons[$i];
    # figure out which exon cluster the exon belongs to
    # but only once for each non redundant exon
    unless ($collapsed_cluster->contains_exon($exon)){
    EC: foreach my $ec (@exon_clusters){
        if ($ec->contains_exon($exon)){
          $exon_cluster = $ec;
          last EC;
        }
      }
    }
    # pass all the exons to the collapsed cluster object
    $collapsed_cluster->add_exon($exon,$exon_cluster);
    # cannot make an intron unless we have 2 exons to look at
    next if $i == 0;
    # make an intron object using exon extended so I can store the evidence sets
    my $intron = $collapsed_cluster->make_intron_from_exons($exons[$i-1],$exons[$i]);
    $intron->ev_set($exon->ev_set);
    $intron->biotype('intron');
    $collapsed_cluster->add_intron($intron);
      }
    }
  }
  # add any rnaseq
  if ( $self->rnaseq_introns ) {
    foreach my $rnaseq_intron ( @{$self->rnaseq_introns} ) {
      # just ones within this gene
      next if $rnaseq_intron->end < $cluster->start;
      last if $rnaseq_intron->start > $cluster->end;
      my $intron = $collapsed_cluster->make_intron_from_rnaseq($rnaseq_intron);
      # treat them like ests
      $intron->ev_set('est');
      $intron->biotype('intron');
      $collapsed_cluster->add_intron($intron);
      # add the extra weight
      $collapsed_cluster->add_intron_weight($intron,$rnaseq_intron->score);
    }
  }

  # once we have added all of them into a nice non-redundant set
  # we need to assign them a score
  foreach my $exon (@{$collapsed_cluster->get_all_exons}){
    my $score = $self->score($collapsed_cluster,$exon,$genes,'exon');
    $collapsed_cluster->exon_score($exon,$score);
  }

  foreach my $intron (@{$collapsed_cluster->get_all_introns}){
    # use all the genes on the strand when looking for introns overlapping genes as it helps
    # stop cluster joining
    my $score = $self->score($collapsed_cluster,$intron,$genes_by_strand->{$intron->strand},'intron');
    $collapsed_cluster->intron_score($intron,$score);
  }

  return $collapsed_cluster;
}

=head2 score

   Name       : score
   Arg[0]     : Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::CollapsedCluster
   Arg[1]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Arg[2]     : Array ref of Bio::EnsEMBL::Gene objects
   Arg[3]     : String 'exon' or 'intron'
   Function   : Assigns scores to introns and exons essentially the number of identical
              : features divided by the number of overlapping non-identical features.
              : Additional weights are added for end exons and est exons.
   Exceptions : none
   Returnval  : Scalar score

=cut

# scores are all assigned in the same way


sub score {
  my ($self,$collapsed_cluster,$feature,$genes,$type) = @_ ;
  my $min_consensus = $self->{min_consensus} ;
  my $end_exon_penalty = $self->{end_exon_penalty} ;
  my $est_overlap_penalty = $self->{est_overlap_penalty} ;
  my $short_intron_penalty = $self->{short_intron_penalty} ;
  my $short_exon_penalty = $self->{short_exon_penalty} ;
  my $verbose = $self->{verbose} ;
  my $total_score;
  # scores are calculated as the number of exoct matches / the number of overlaps
  # it is complicated by end exons and evidence sets however.
  $total_score = $collapsed_cluster->get_intron_count($feature) if $type eq 'intron';
  $total_score = $collapsed_cluster->get_exon_count($feature) if $type eq 'exon';
  $collapsed_cluster->exon_info($feature) if $verbose;
  print " " if $verbose;
  my $score = 0 ;
  # this way we penalise anything with less than exons than the minumum number of overlapping exons
  my $overlap = $min_consensus;
  my $est_trans = $min_consensus;


  # the weighted score is:
  # 50% - proportion of exons in the stack that come from ests / number of overlapping est transcripts
  # 50% - the non est exons over the total number of overlapping transcripts of any type

  # count the number of overlapping exons / introns
  foreach my $gene (@$genes){
    foreach my $tran (@{$gene->get_all_Transcripts}){
      # if its an exon just add up how many genes of each type are in the cluster
      if ($type eq 'exon'){
    $est_trans++ if $tran->ev_set eq 'est';
    $overlap++;
      }
      next if scalar(@{$tran->get_all_Exons}) == 1;
      if ($type eq 'intron'){
    # for introns; count the number of overlapping trans from the entire slice
    # this creates a penatly for introns that span genes that
    # are not in the cluster because they dont share exons
    if ($feature->start < $tran->end() && $feature->end > $tran->start()){
      $est_trans++ if $tran->ev_set eq 'est';
      $overlap++;
    }
      }
    }
  }
  # number of est exons that overlap the feature
  my $est_exons = $collapsed_cluster->ev_count_by_feature($feature,'est');

  # for end exon scores we allow the external boundary to be non identical
  if ($type eq 'exon' && $feature->is_terminal_exon &&  $feature->number_exons > 1){
    print "End exon $est_exons " if $verbose;
    # count the total number of overlapping end exons
    $est_exons = $collapsed_cluster->ev_count_by_terminal_feature($feature,'est');
    $total_score = $collapsed_cluster->count_end_exon($feature);
    $self->throw("No evidence found for exon " . $collapsed_cluster->exon_info($feature) . "\n")
      unless defined($est_exons);
  }
  print "$type " if $verbose;
  print "TOTAL $total_score OVERLAP $overlap " if $verbose;

  if ($est_exons){
    print "WEIGHTED GENES $est_exons EST TRANS $est_trans " if $verbose;
    # number of exons of the weighted type / number of overlapping trans of the weighted type;
    my $weighted_score = 0;
    my $unweighted_score = 0;
    if ($est_trans){
      $weighted_score = $est_exons  / $est_trans ;
      print "weighted score $weighted_score " if $verbose;
    }
    # other componnt of score comes from all the exons
    $unweighted_score =  $total_score  / $overlap  ;
    print "unweighted score $unweighted_score \t" if $verbose;
    $score = ($weighted_score + $unweighted_score) / 2;
  } else {
    # there are no est exons overlaping this feature score is simple identical / non-identical overlaps
    $score += ($total_score / $overlap ) ;
    # if the feature overlaps an est transcript but not est exons, add a penalty
    if ($est_trans &&  $est_trans > 1.5 ){
      $score -= $est_overlap_penalty;
    }
  }
  # favor long end exons
  if ($feature->is_terminal_exon && $type eq 'exon' && $feature->number_exons > 1){
    # longest end exon
    my $longest_end_exon = $collapsed_cluster->get_end_exon($feature)->length;
    my $multiplyer = $feature->length / $longest_end_exon;
    print "End exon penalty ".$feature->length." / ". $longest_end_exon." -.5\t" if $verbose;
    if ($score >= 0){
      $score *= $multiplyer;
    } else {
      $score /= $multiplyer;
    }
    # end exons are penalised to prevent 'spindly' exons
    $score = $score - $end_exon_penalty;
    print " ENDEXON " if $verbose;
  } else {
    # penalty for  internal short exons / introns
      $score-= 0.5 if $feature->length < $short_intron_penalty && $type eq 'intron';
      $score-= 0.5 if $feature->length < $short_exon_penalty   && $type eq 'exon';
  }
  print "$score \n" if $verbose;
  return $score;
}

=head2 make_genes

   Name       : make_genes
   Arg[0]     : Array ref of Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended
   Arg[1]     : Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster
   Arg[2]     : Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::CollapsedCluster
   Function   : Takes an array ref of transcript objects and makes them into gene objects.
              : Deals with translations and orders the genes by score stores a subset
              : of the highest scoring genes on the output array.
   Exceptions : Throws if it cannot work out the translation start or end exons
   Returnval  : none

=cut

sub make_genes{
  my ($self,$transcripts,$cluster,$collapsed_cluster)= @_ ;
  my $add_utr = $self->{add_utr} ;
  my $min_consensus = $self->{min_consensus} ;
  my $good_percent = $self->{good_percent} ;
  my $good_biotype = $self->{good_biotype} ;
  my $small_biotype = $self->{small_biotype} ;
  my $bad_biotype = $self->{bad_biotype} ;
  my $verbose = $self->{verbose} ;
  my @genes;
  my @final_genes;

  my $store_bad = defined($bad_biotype);
  if ( !$store_bad ) {
    $bad_biotype = 'bad';
  }

  $self->throw("No Transcripts to use\n") unless defined($transcripts);
  SCORE:  foreach my $t (sort {$b->score <=> $a->score } @$transcripts){
    my $transcript = $self->clone_transcript_extended($t);
    my $cds_start_exon;
    my $cds_end_exon;
    my $similarity_tran = $transcript->similarity;
    # each smilarity gene only allowed once, just the highest scoring one

    my $translation = Bio::EnsEMBL::Translation->new();
    my $strand= 1;
    my $slice =  $similarity_tran->slice;
    my $weighted_score;
    if ($add_utr){
      # compute the translation and then check that the old translation 
      # is a subset of the new one     
      print "Add UTR: Transcript ".$transcript->start .":".$transcript->end .":".$transcript->strand." ".$transcript->biotype." " if $verbose;
      print scalar(@{$transcript->get_all_Exons})." Exons\n" if $verbose;

      # dont compute the translation just set it
      # need the similarity gene that it was built from
      # so I can set the cds start and end on the transcript
      # if the start is in the same exon in the new transcript
      # just transfer it, if it is over the end of the new exon
      # trim back so that it is in the exon if you see what I mean
      # print stuff
      my $last_exon;
      foreach my $exon (@{$transcript->get_all_Exons}){
        if ($last_exon){
          my $intron = $collapsed_cluster->make_intron_from_exons($last_exon,$exon);
        }
        $last_exon = $exon;
      }
      if ($transcript->strand == -1){
        my $invertedslice = $similarity_tran->slice->invert;
        $strand = -1;
        $similarity_tran = $similarity_tran->transfer($invertedslice);
        $transcript = $transcript->transfer($invertedslice);
      }
      my $cds_start = $similarity_tran->coding_region_start;
      my $cds_end = $similarity_tran->coding_region_end;
      my @exon_array = @{$transcript->get_all_Exons};
      my @similarity_exons = sort { $a->start <=> $b->start } @{$similarity_tran->get_all_Exons};
      my @coding_similarity_exons = sort { $a->start <=> $b->start } @{$similarity_tran->get_all_translateable_Exons};


      my $number = scalar(@exon_array);
      # translation start and end is upside down when on the opposite strand
      for ( my $i = 0 ; $i < scalar(@exon_array) ; $i ++ ){
        # check for overlaps 
        my $exon = $exon_array[$i];	  
        # set phases to non coding by default
        $exon->phase( -1 );
        $exon->end_phase( -1 );
        foreach my $se (@similarity_exons){
          # set coding exons to have same phase as similarity exons
          if ($exon->start == $se->start){
            $exon->phase($se->phase);
          }
          if ( $exon->end == $se->end ){
            $exon->end_phase($se->end_phase);
          }
        }
        if ($exon->start <= $cds_start && $exon->end >= $cds_start){
          # translation start is in this exon
          $translation->start_Exon($exon);
          $translation->start($cds_start - $exon->start+1);
          $exon->phase( -1 ) 	if $exon->start < $cds_start;
        }
        if ($exon->start <= $cds_end && $exon->end >= $cds_end){
          # translation stop is in this exon
          $translation->end_Exon($exon);

          $translation->end($cds_end - $exon->start+1);

          if ( $exon->end_phase > 0) {
            $translation->end($translation->end - $exon->end_phase);
          }

          if ( $exon->end > $cds_end ||  $exon->end_phase > 0) {
            print "# Setting end phase -1\n";
            $exon->end_phase( -1 );
          }
        }


        ## find the utr exon that shares the internal boundary with the first simgw exon
        # set the start and end coding exon, using the coding exon array as a pointer

        my $coding_start_index;
        my $coding_end_index;  

        for ( my $i = 0 ; $i < scalar(@similarity_exons) ; $i ++ ) {

          if ( $similarity_exons[$i]->start <= $coding_similarity_exons[0]->start &&
               $similarity_exons[$i]->end >= $coding_similarity_exons[0]->end ) { 
            $coding_start_index = $i;
            print "Coding start exon index set: ". $coding_start_index  ."\n" if $verbose; 
          }         
          print $similarity_exons[$i]->start ." <= ". $coding_similarity_exons[-1]->start ." &&\n". 
                $similarity_exons[$i]->end ." >= ". $coding_similarity_exons[-1]->end."\n";


          if ( $similarity_exons[$i]->start <= $coding_similarity_exons[-1]->start &&
                $similarity_exons[$i]->end >= $coding_similarity_exons[-1]->end ) { 
            $coding_end_index = $i; 
            print "Coding end exon index set: ". $coding_end_index  ."\n" if $verbose;
          }
        }
        
        # find the utr exon that shares the internal boundary with the first simgw exon
        if ($exon->end == $similarity_exons[$coding_start_index]->end ){
          $cds_start_exon = $exon;
        }
        # find the utr exon that shares the internal boundary with the last simgw exon
        if ($exon->start == $similarity_exons[$coding_end_index]->start){
          $cds_end_exon = $exon;
        }
      }
      print scalar(@{$transcript->get_all_Exons})." Exons\n" if $verbose;

     
      unless ( $translation->start_Exon && $translation->start_Exon eq $cds_start_exon ) {

        # the start is not inside an exon
        throw("Cannot find the coding sequence start exon\n") unless $cds_start_exon;
        # trim back to the exon that matches the internal boundary
        # stay in frame though
        $translation->start_Exon($cds_start_exon);
        my $trim = $cds_start_exon->start - $cds_start;
        
        if ( $trim%3 == 0 ){
          $translation->start(1);
        } else {	
          $translation->start(1+(3-$trim%3));
        }
        $cds_start_exon->phase( -1 ) 	if $cds_start_exon->start < $cds_start;
      }

      unless ($translation->end_Exon && $translation->end_Exon eq $cds_end_exon ){
        # the end is not inside an exon
        throw("Cannot find the coding sequence end exon\n") unless $cds_end_exon;
        # trim back to the exon that matches the internal boundary
        # stay in frame 
        $translation->end_Exon($cds_end_exon);
        my $trim = $cds_end - $cds_end_exon->end;
        if ( $trim%3 == 0 ){
          $translation->end($cds_end_exon->length);
        } else {
          $translation->end($cds_end_exon->length-(3-$trim%3));
        }
        $cds_end_exon->end_phase( -1 ) 	if $cds_end_exon->end > $cds_end;
      }

      $transcript->translation($translation);
      if ($strand == -1){
        # flip it back onto the - strand
        my $slice = $transcript->slice->invert;
        $transcript = $transcript->transfer($slice);
        $similarity_tran = $similarity_tran->transfer($slice);
      }	
      # set the CDS start end to be the same as the similarity
      # now is the new trimmed transcript compatible the longest possible orf?
      # make a copy of the transcript so as not to disturb the exon phases of the original
      my $test_transcript =  $self->clone_transcript_extended($transcript);
      # translate the copy 
      my $new_tran = compute_translation($test_transcript); 
      my $check_seq = $transcript->translate->seq;
      # check if the modified translation has worked properly
      my $fail = 0 ;
      if ( $check_seq =~ /\*/ ) {
        print "Modified translation contains stops replacing it with computed translation\n";
        $fail = 1;
      }
      unless ( are_phases_consistent($transcript) ) {
        $fail = 1;
        print "Modified translation is replaced by computed translation\n"; 
      }
      if ( $fail ) {
        $transcript = $new_tran;
      } else {
        # check old translation is contained within the new one
        my $new_seq = $new_tran->translation->seq;
        #	 print "NEW $new_seq\nOLD $check_seq\n";
        my $length_change = length($new_seq) - length($check_seq);
        print "Length change $length_change\n";
        my $extended = 0;
        if ( $length_change > 0 ) {
          if ( $new_seq =~ /$check_seq/ ) {
            print "MATCH \n";
            # go with the new translation
            $transcript = $new_tran;
            print "Extending translation by $length_change residues\n";
          }
        }
      }
      # check for the number of non coding exons 
      if (  scalar(@{$transcript->get_all_Exons}) - scalar(@{$transcript->get_all_translateable_Exons}) > 3 ){
        # just stick with the similarity model
        $transcript = $similarity_tran;
      }
    } else {
      # dont  need to mess with translations if I am not adding UTR...
      my $newtranslation = clone_Translation($similarity_tran,$transcript);
      $transcript->translation($newtranslation);
      # fix phases
      my @exon_array = @{$transcript->get_all_Exons};
      my @similarity_exons = sort { $a->start <=> $b->start } @{$similarity_tran->get_all_Exons};
      for ( my $i = 0 ; $i < scalar(@exon_array) ; $i ++ ){
        # check for overlaps 
        my $exon = $exon_array[$i];
        foreach my $se (@similarity_exons){
          # set coding exons to have same phase as similarity exons
          if ($exon->start == $se->start){
            $exon->phase($se->phase);
          }
          if ( $exon->end == $se->end ){
            $exon->end_phase($se->end_phase);
          }
        }
      }
    }

    # penalise transcripts where a utr exon overlaps another exon that was previously coding
    $weighted_score = $self->compare_translation($transcript,$collapsed_cluster);
    print $transcript->score."\t" if $verbose;
    print "\npeptide weighted score $weighted_score\n\n" if $verbose;
    $transcript->score( $weighted_score );
    $transcript->analysis($similarity_tran->analysis);
    $transcript->slice($similarity_tran->slice);
    my $gene = Bio::EnsEMBL::Gene->new;
    $gene->slice($transcript->slice);
    $gene->add_Transcript($transcript);
    $gene->analysis($transcript->analysis);
    push @genes,$gene;
  }


  # add the models of the single exon genes
  push @genes, @{$self->add_single_exon_genes($cluster,$collapsed_cluster)};

  # SMJS Test length scaling
  $self->weight_scores_by_cdna_length(\@genes);

  my $iteration = 0;

  while (@genes) {
    @genes = sort {
      $b->get_all_Transcripts->[0]->score()
      <=> $a->get_all_Transcripts->[0]->score()
    } @genes;

    my $top_score    = $genes[0]->get_all_Transcripts->[0]->score();
    my $bottom_score = $genes[-1]->get_all_Transcripts->[0]->score();

    # prevent the whole division by zero thing
    if ( $bottom_score == 0 ) {
      $bottom_score -= 0.0001;
    }

    my $score_range = ( $top_score - $bottom_score );
    print("SCORE RANGE $top_score - $bottom_score = $score_range\n");

    my @genes_to_consider;

    foreach my $gene (@genes) {
      foreach my $transcript ( @{ $gene->get_all_Transcripts() } ) {
        # Sort out exon phases - transfer them from the similarity gene.
        # Coding region start/stop is independent of strand so I will
        # swap them if the strand is -ve.

        my $biotype = $bad_biotype;

        if ( scalar(@genes) < $min_consensus ) {
          $biotype = $small_biotype;
        }
        else {
          # how about if we include the top 5% ($good_percent) of the
          # gene scores
          if ($good_percent) {
            my $top = ( $transcript->score() - $bottom_score );
            if ( $score_range &&
              ( ( $top/$score_range )*100 >= 100 - $good_percent ) )
            {
              $biotype = $good_biotype;
            }
          }
          elsif ( $transcript->score() == $top_score ) {
            # we want the single top scoring model
            $biotype = $good_biotype;
          }
        }

        $gene->biotype($biotype);

        push( @genes_to_consider, $gene );
      } ## end foreach my $transcript ( @{...})
    } ## end foreach my $gene (@genes)

    my ( $genes_to_store, $genes_to_recluster ) =
    $self->_find_genes_to_recluster( \@genes_to_consider );

    printf( "Iteration %d: Cluster of %d models " .
      "broken up into %d finished models " .
      "and %d models to reconsider\n",
      ++$iteration,
      scalar(@genes_to_consider),
      scalar( @{$genes_to_store} ),
      scalar( @{$genes_to_recluster} ) );

    if ( @{$genes_to_recluster} ) {
      printf( "First model to reconsider is located here: %s\n",
        $genes_to_recluster->[0]->feature_Slice()->name() );
    }

    foreach my $gene ( @{$genes_to_store} ) {
      if ( ( $store_bad && $gene->biotype() eq $bad_biotype ) ||
        ( $gene->biotype() eq $good_biotype )  ||
        ( $gene->biotype() eq $small_biotype ))
      {
        push( @final_genes, $gene );
      }
    }
    @genes = @{$genes_to_recluster};
  } ## end while (@genes)

  foreach my $final (@final_genes) {
    foreach my $transcript ( @{ $final->get_all_Transcripts() } ) {
      if ( !are_phases_consistent($transcript) ) {
        throw(
          Transcript_info($transcript) . " has inconsistent phases" );
      }

      # Modify score to be independent of length i total / number of
      # introns+exons.  Should always be less than 1.
      my $final_score =
      $transcript->score()/
      ( scalar( @{ $transcript->get_all_Exons() } )*2 - 1 );

      # Add score into trancript supporting feature score feild and
      # move coverage into the new covergae field.
      foreach
      my $tsf ( @{ $transcript->get_all_supporting_features() } )
      {
        $tsf->hcoverage( $tsf->score() );
        $tsf->score( sprintf( "%.4f", $final_score ) );
      }

      print("Final score $final_score\n");
    }
  } ## end foreach my $final (@final_genes)

  $self->output( \@final_genes );
}

sub _find_genes_to_recluster {
  # This method finds the gene models that were tagged with the "bad"
  # biotype by make_genes() but that does not have exon overlap
  # with any of the gene models tagged as "good".  If there are any
  # such gene models, it usually means that there are actually two
  # (or more) clusters of models that are joined by one or a few
  # connecting models.  These "bad" models may then be re-considered by
  # make_genes() in one or several extra iterations.

  my ( $self, $genes ) = @_;

  my $good_biotype = $self->{'good_biotype'};
  my $bad_biotype  = $self->{'bad_biotype'};

  if ( !defined($bad_biotype) ) {
    $bad_biotype = 'bad';
  }

  # Build hash of exon coordinates from the "good" models.

  my %exon_coords;

  foreach my $gene ( @{$genes} ) {
    if ( $gene->biotype() eq $good_biotype ) {
      foreach my $transcript ( @{ $gene->get_all_Transcripts() } ) {
        foreach my $exon ( @{ $transcript->get_all_Exons() } ) {
          my $slice = $exon->feature_Slice();
          my $key   = $slice->seq_region_name() . $slice->strand();
          push( @{ $exon_coords{$key} },
                [ $slice->start(), $slice->end() ] );
        }
      }
    }
  }

  # Compare the "bad" models against the exon coordinates of the "good"
  # ones.  If there is overlap, add them to @genes_to_store, otherwise
  # to @genes_to_recluster.

  my @genes_to_store;
  my @genes_to_recluster;

  foreach my $gene ( @{$genes} ) {
    if ( $gene->biotype() eq $bad_biotype ) {
      my $overlaps = 0;

    TRANSCRIPT_LOOP:
      foreach my $transcript ( @{ $gene->get_all_Transcripts() } ) {
        foreach my $exon ( @{ $transcript->get_all_Exons() } ) {
          my $slice = $exon->feature_Slice();
          my $key   = $slice->seq_region_name() . $slice->strand();

          if ( exists( $exon_coords{$key} ) ) {
            foreach my $range ( @{ $exon_coords{$key} } ) {
              if ( ( $slice->start() <= $range->[0] &&
                     $slice->end() >= $range->[0] ) ||
                   ( $slice->start() <= $range->[1] &&
                     $slice->end() >= $range->[1] ) ||
                   ( $slice->start() >= $range->[0] &&
                     $slice->end() <= $range->[1] ) )
              {
                $overlaps = 1;
                last TRANSCRIPT_LOOP;
              }
            }
          }
        }
      }

      if ($overlaps) {
        # Store "bad" models that overlap with any "good" model.
        push( @genes_to_store, $gene );
      }
      else {
        # Reconsider the "bad" models that do not overlap any "good"
        # models.
        push( @genes_to_recluster, $gene );
      }

    } ## end if ( $gene->biotype() ...)
    else {
      # Store "good" and "small" models.
      push( @genes_to_store, $gene );
    }
  } ## end foreach my $gene ( @{$genes...})

  return ( \@genes_to_store, \@genes_to_recluster );
} ## end sub _find_genes_to_recluster

sub weight_scores_by_cdna_length {
  my ($self,$genes) = @_ ;

  my @score_sorted_genes = sort { $b->get_all_Transcripts->[0]->score <=> $a->get_all_Transcripts->[0]->score } @$genes;
  my $top_score = $score_sorted_genes[0]->get_all_Transcripts->[0]->score;
  my $bottom_score = $score_sorted_genes[-1]->get_all_Transcripts->[0]->score;
  # prevent the whole division by zero thing

  my $score_range = ($top_score - $bottom_score);

  my @length_sorted_genes = sort { $a->get_all_Transcripts->[0]->length <=> $b->get_all_Transcripts->[0]->length } @$genes;
  my $minlen = $length_sorted_genes[0]->get_all_Transcripts->[0]->length;
  my $maxlen = $length_sorted_genes[-1]->get_all_Transcripts->[0]->length;

  my $len_range = $maxlen-$minlen;

  if (!$len_range) { return };

  #if ($len_range/$maxlen < 0.7 ) {
  #  print "!!!!!!! Returning on range limit\n";
  #  return;
  #}

  my $scale_factor = $score_range/1;

  foreach my $gene (@$genes) {
    my $trans = $gene->get_all_Transcripts->[0];
    my $length =  $trans->length;

    my $factor = ($length - $minlen) / $len_range;

     print "factor = $factor length = $length minlen = $minlen len_range = $len_range scale_factor = $scale_factor\n";
     print "Changed score for " . $gene->analysis->logic_name . " gene from " . $trans->score;
    $trans->score( $trans->score + $factor*$scale_factor);
     print " to " . $trans->score . "\n";
  }
}

=head2 add_single_exon_genes

   Name       : add_single_exon_genes
   Arg[0]     : Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster
   Arg[1]     : Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::CollapsedCluster
   Function   : Assigns scores to single exon genes
   Exceptions : Throws if the exon has not been scored
   Returnval  : Array ref of Bio::EnsEMBL::Gene

=cut


sub add_single_exon_genes {
  my ($self,$cluster,$collapsed_cluster) = @_ ;
  my @single_exon_genes;
  my $genes = $cluster->get_Genes;
  foreach my $gene (@$genes){
    next unless (scalar(@{$gene->get_all_Exons}) == 1);
    my $trans = $gene->get_all_Transcripts->[0];
    next unless $trans->ev_set eq 'simgw';
    my $exon = $trans->get_all_Exons->[0];
    my $score =  $collapsed_cluster->exon_score($exon);
    $self->throw("Score not found for exon ".$collapsed_cluster->exon_info($exon)."\n")
      unless defined($score);
    $trans->score($score);
    push @single_exon_genes, $self->clone_gene($gene);
  }
  return \@single_exon_genes;
}

=head2 make_transcripts

   Name       : make_transcripts
   Arg[0]     : Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster
   Arg[1]     : Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::CollapsedCluster
   Function   : Combines est and similarity trascripts to make coding genes with UTR.
              : Joins transcripts if they match the internal boundary of the first or last
              : exon of the similarity gene.
   Exceptions : none
   Returnval  : Array ref of Bio::EnsEMBL::Transcript

=cut


sub make_transcripts {
  my ($self,$cluster,$collapsed_cluster) = @_ ;
  my $add_utr = $self->{add_utr} ;
  my $genes = $cluster->get_Genes;
  my @est_genes;
  my @similarity_genes;
  my $longest_similarity = 0;
  my @all_transcripts;
  foreach my $gene (@$genes){
    push @est_genes, $gene if  $gene->get_all_Transcripts->[0]->ev_set eq 'est';
    push @similarity_genes, $gene if $gene->get_all_Transcripts->[0]->ev_set  eq 'simgw';
  }

  foreach my $similarity (@similarity_genes) {
    foreach my $trans (@{$similarity->get_all_Transcripts}) {
      my @partial_transcripts;
      my @possible_transcripts;
      my @possible_start_utr;
      my @possible_end_utr;
      my @transcript_exons = sort {$a->start <=> $b->start } @{$trans->get_all_Exons};
      # single exon genes are dealt with elsewhere
      next if (scalar(@transcript_exons) <= 1);
      my $exon_count = $#transcript_exons;
      my $start_splice;
      my $end_splice;
      my @similarity_exons;
      my $almost_complete;
      my @start_utr;
      my @end_utr;
      # using start and end rather than 3' and 5' to save strand problems
      # start UTR
      $start_splice = $collapsed_cluster->make_intron_from_exons($transcript_exons[0],$transcript_exons[1]);
      $self->throw("Intron was not previously defined in the cluster\n")
    unless $collapsed_cluster->contains_intron($start_splice);
      my @joining_exons = @{$collapsed_cluster->last_exons_by_intron($start_splice)};
      foreach my $exon (@joining_exons) {
    # get all the transcripts belonging to the non redundant potential joining exons
    my @transcripts = @{$collapsed_cluster->transcript_by_feature($exon)};
    foreach my $trans (@transcripts) {
      # only want est transcripts
      next unless $trans->ev_set eq 'est' ;
      push @start_utr, $trans if $add_utr;
    }
      }
      # end UTR
      $end_splice = $collapsed_cluster->make_intron_from_exons($transcript_exons[$exon_count-1],$transcript_exons[$exon_count]);
      $self->throw("Intron was not previously defined in the cluster\n")
    unless $collapsed_cluster->contains_intron($end_splice);
      @joining_exons = @{$collapsed_cluster->next_exons_by_intron($end_splice)};
      foreach my $exon (@joining_exons) {
    # get all the transcripts belonging to the non redundant potential joining exons
    my @transcripts = @{$collapsed_cluster->transcript_by_feature($exon)};
    foreach my $trans (@transcripts) {
      # only want est transcripts
      next unless $trans->ev_set eq 'est' ;
      push @end_utr, $trans if $add_utr;
    }
      }

      # if we have start utr
      if ( scalar(@start_utr) > 0 ) {
    # make and score transcripts that represent all the start UTRs that will fit this similarity gene
    foreach my $est (@start_utr) {
      my @exon_array;
      # truncate the ests to the right length
      foreach my $exon (@{$est->get_all_Exons}){
        if ($exon->end <= $start_splice->start){
          push @exon_array, $exon;
        }
      }
      my $new_transcript = $self->transcript_from_exons(\@exon_array,undef,$collapsed_cluster)
        if scalar(@exon_array > 0);
      if ($new_transcript){
        # keep the transcript supporting features
        $new_transcript->add_supporting_features(@{$est->get_all_supporting_features});
        push @possible_start_utr, $new_transcript;
      }
    }
    # only add the highest scoring possible utr for this similarity gene
    @possible_start_utr = sort { $a->score <=> $b->score } @possible_start_utr if @possible_start_utr;
    push @partial_transcripts, $possible_start_utr[-1] if @possible_start_utr;
      }
      # make a new transcript conatining just the first similarity exon
      my @exon_array = (  $transcript_exons[0]  );
      my $new_transcript = $self->transcript_from_exons(\@exon_array,undef,$collapsed_cluster);
      $self->throw("Unable to make transcript using first exon of similarity gene\n")
    unless $new_transcript;
      push @partial_transcripts, $new_transcript if $new_transcript;

      # now we are in the similarity gene proper put together all the internal exons
      @similarity_exons = @{$trans->get_all_Exons};
      pop @similarity_exons;
      shift @similarity_exons;

      #  add the similarity exons to all the transcripts made so far
      foreach my $transcript (@partial_transcripts) {
    my $score;
    $almost_complete = $self->transcript_from_exons(\@similarity_exons,$transcript,$collapsed_cluster);
    if ( scalar(@end_utr) > 0 ) {
      # add all the possible end utrs
      # also make one transcript without end utr
      foreach my $est (@end_utr) {
        my @exon_array;
        # make an exon array out of the est starting at the end_splice intron
        foreach my $exon (@{$est->get_all_Exons}){
          if ($exon->start >= $end_splice->end){
        push @exon_array, $exon;
          }
        }
        # to all the possible transcripts made so far
        my $new_transcript = $self->transcript_from_exons(\@exon_array,$almost_complete,$collapsed_cluster)
          if scalar(@exon_array > 0);
        if ($new_transcript){
          # keep the supporting features
          $new_transcript->add_supporting_features(@{$est->get_all_supporting_features});
          $new_transcript->similarity($trans);
          push @possible_end_utr, $new_transcript;
        }
      }
      # only add the highest scoring possible utr for this similarity gene
       @possible_end_utr = sort { $a->score <=> $b->score } @possible_end_utr;
       push @possible_transcripts, $possible_end_utr[-1];
     }
     #  add the last exon of the similarity gene
     my @exon_array = (  $transcript_exons[-1]  );
     my $new_transcript = $self->transcript_from_exons(\@exon_array,$almost_complete,$collapsed_cluster);
     $new_transcript->similarity($trans);
     $self->throw("Unable to make transcript using last exon of similarity gene\n")
       unless $new_transcript;
     push @possible_transcripts, $new_transcript if $new_transcript;
       }
      # only store the highest scoring transcript for each similarity gene
      @possible_transcripts = sort { $a->score <=> $b->score } @possible_transcripts if @possible_transcripts;
      my $chosen_transcript = $possible_transcripts[-1];
      # keep the transcript supporting feature of the similarity gene
      $chosen_transcript->add_supporting_features(@{$trans->get_all_supporting_features});
      $chosen_transcript->biotype($trans->biotype);
      push @all_transcripts, $chosen_transcript ;
    }
  }
  return \@all_transcripts;
}

=head2 transcript_from_exons

   Name       : transcript_from_exons
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended
   Arg[1]     : (optional) Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::TranscriptExtended
   Arg[2]     : Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::CollapsedCluster
   Function   : Makes a transcript object from an array of exons or extends an existing transcript
              : Assigns the new transcript a score by summing all the intron and exon scores
   Exceptions : Throws unless ARG[1] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended object
              : Throws if an exon or intron has not been previously assigned a score
   Returnval  : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::TranscriptExtended

=cut

sub transcript_from_exons {
  my ($self, $exons, $transcript, $collapsed_cluster) = @_ ;
  my $utr_penalty = $self->{utr_penalty} ;
  my $verbose = $self->{verbose} ;
  my $new_transcript = Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended->new();
  my $prev_exon;
  my $score;
  @$exons = sort { $a->start <=> $b->start } @$exons;

  # if we are extending an existing transcript the previous exon is the last in the transcript
  if (defined($transcript)){
    $self->throw("Transcript needs to be a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended not a $transcript\n")
      unless $transcript->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended");
    $score = $transcript->score;
    # clone the transcript array
    foreach my $exon (sort { $a->start <=> $b->start } @{$transcript->get_all_Exons}){
      $new_transcript->add_Exon($exon);
      $prev_exon = $exon;
    }
  }

  foreach my $exon (@$exons){

    my $collapsed_exon = $collapsed_cluster->get_exon($exon);
      unless ( $exon ){
    $collapsed_cluster->exon_info($collapsed_exon) if $verbose;
    $self->throw("No score assigned to this exon\n");
      }
    if ($prev_exon){
      my $intron = $collapsed_cluster->make_intron_from_exons($prev_exon,$exon);
      unless ( $intron ){
    $collapsed_cluster->exon_info($intron) if $verbose;
    $self->throw("No score assigned to this intron\n");
      }
      $score += $collapsed_cluster->intron_score($intron);
      # add a penalty to UTR exons
      $score -= $utr_penalty if $exon->ev_set eq 'est';
    }
    $score += $collapsed_cluster->exon_score($collapsed_exon);
    $score -= $utr_penalty if $exon->ev_set eq 'est';
    $prev_exon = $exon;
    $new_transcript->add_Exon($collapsed_exon);
  }
  $new_transcript->score($score);
  return $new_transcript;
}

=head2 filter_genes

   Name       : filter_genes
   Arg[0]     : Array ref of Bio::EnsEMBL::Gene objects
   Function   : Pre-filters the gene set on the basis of number of consensus splice sites or unsuported exons
   Exceptions : Throws unless ARG[0] is a Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object
   Returnval  : Array ref of Bio::EnsEMBL::Gene objects

=cut

# does some pre-processing on the gene set to remove dodgy looking genes

sub filter_genes {
  my ( $self, $genes) = @_ ;
  my $filter_singletons = $self->{filter_singletons} ;
  my $filter_non_consensus = $self->{filter_non_consensus} ;
  my $filter_ests = $self->{filter_ests} ;
  my $verbose = $self->{verbose} ;
  #print "FILTERING TRANSCRIPTS\n";
  #print "Have ".@$genes." genes to filter\n";
  my @genes_to_filter;
  my %singletons;
  my %non_con;
  my @filtered_genes;
  my %delete;
  my %exon_hash;
  my $gene_count = 0;

  for ( my $g = 0 ; $g < scalar(@$genes) ; $g++ ){
    # dont run on ests unless specifically asked to
    unless ($filter_ests){
      next if ($genes->[$g]->get_all_Transcripts->[0]->ev_set eq 'est');
    }
    $gene_count++;
    foreach my $trans (@{$genes->[$g]->get_all_Transcripts}){
      # count non consensus splices
      $non_con{$g} =  count_non_canonical_splice_sites($trans);

      my @exons = @{$trans->get_all_Exons};
      foreach (my $i =0 ; $i < scalar(@exons) ; $i++){
    my $exon = $exons[$i];
    # need to be able to backtrack, which gene exon comes from..
    my $unique_key = $exon->start.":".$exon->end.":".$exon->strand;
    push @{$exon_hash{$unique_key}} , $g;
      }
    }
  }
  # dont run on very small clusters
  if ($gene_count <= 2){
    return $genes;
  }
  # once we have looked at all the exons in the cluster find the singletons
  foreach my $key (sort keys %exon_hash){
    if (scalar(@{$exon_hash{$key}}) == 1){
      # singleton
      $singletons{$exon_hash{$key}[0]}++;
      $delete{$exon_hash{$key}[0]} = 1;
    }
  }
  # gotcha
  if ( scalar (keys %delete) >= $gene_count - 2 ){
    # dont delete them if you will end up deleting everything!
    print "All transcripts look dodgy, not filtering\n" if $verbose;
    return $genes;
  } else {
    for ( my $g = 0 ; $g < scalar(@$genes) ; $g++ ){
      if ( $filter_singletons && $singletons{$g} && $singletons{$g}>= $filter_singletons ){
    print  "Filtering out ".$genes->[$g]->dbID." ".$genes->[$g]->biotype." ".$genes->[$g]->start.":".$genes->[$g]->end.":".$genes->[$g]->strand." " if $verbose;
    print  $singletons{$g}." singleton exons\n" if $verbose;
    next;

      }
      if ( $filter_non_consensus && $non_con{$g} && $non_con{$g} >= $filter_non_consensus ){
    print  "Filtering out ".$genes->[$g]->dbID." ".$genes->[$g]->biotype." ".$genes->[$g]->start.":".$genes->[$g]->end.":".$genes->[$g]->strand." " if $verbose;
    print  $non_con{$g}." non - consensus\n" if $verbose;
    next;
      }
      push @filtered_genes , $genes->[$g];
    }
  }
  return \@filtered_genes;
}

=head2 compare_translation

   Name       : compare_translation
   Arg[0]     : Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::TranscriptExtended
   Arg[1]     : Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::CollapsedCluster
   Function   : Post processing to penalise transcripts made from fragments of similarity genes
              : that have been extended with UTR
   Exceptions : Throws if no ExonCluster is stored w.r.t the exon in question.
   Returnval  : Scalar score

=cut

sub compare_translation {
  my ($self,$transcript,$collapsed_cluster) = @_ ;
  my $verbose = $self->{verbose} ;
  my $score = $transcript->score;
  # just look at non coding exons
 EXON: foreach my $exon (@{$transcript->get_all_Exons}){
    # ignore exon that contains the start or stop
    if ($exon->end > $transcript->coding_region_start &&
    $exon->start < $transcript->coding_region_end ){
      next;
    }
    my $exon_cluster = $collapsed_cluster->get_exon_cluster($exon);
    $self->throw("Exon cluster not found for exon ".$collapsed_cluster->exon_info($exon)."\n")
      unless ($exon_cluster);
    my $coding_exons = 0;
    # get all the overlapping exons
    foreach my $overlapping_exon (@{$exon_cluster->get_all_Exons_of_EvidenceSet('simgw')}){
      $coding_exons++;
    }
    if ($coding_exons > 0){
      $collapsed_cluster->exon_info($exon) if $verbose;
      print " overlaps $coding_exons  coding exons penalty -.5\n" if $verbose;
      $score -= 0.5;
    }
  }
  return $score;
}


sub clone_gene {
  my ($self,$gene) =@_ ;
  my $new_gene = Bio::EnsEMBL::Gene->new();
  foreach my $transcript ( @{$gene->get_all_Transcripts} ){
    my $new_transcript = $self->clone_transcript_extended( $transcript );
    $new_gene->analysis( $new_transcript->analysis );
    $new_gene->slice( $new_transcript->slice );
    $new_gene->add_Transcript( $new_transcript );
  }
  $new_gene->biotype( $gene->biotype );
  return $new_gene;
}

sub clone_transcript_extended {
  my ($self,$transcript) =@_ ;
  $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::TranscriptExtended not a ".ref($transcript)."\n")
    unless $transcript->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended");
  my $newtranscript = Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptExtended->new();
  foreach my $exon ( @{$transcript->get_all_Exons}){
    my $new_exon = $self->clone_exon_extended($exon);
    $newtranscript->add_Exon($new_exon);
  }
  foreach my $sf(@{$transcript->get_all_supporting_features}){
    my $newsf = clone_Evidence($sf);
    $newsf->dbID(undef);
    $newtranscript->add_supporting_features($newsf);
  }
  my $attribs = $transcript->get_all_Attributes();
  if ( $transcript->translation ) {
    $newtranscript->translation( clone_Translation( $transcript,$newtranscript ));
  }
  $newtranscript->add_Attributes(@$attribs);
  $newtranscript->slice($transcript->slice);
  $newtranscript->biotype($transcript->biotype);
  $newtranscript->score($transcript->score);
  $newtranscript->ev_set($transcript->ev_set);
  $newtranscript->similarity($transcript->similarity);
  $newtranscript->analysis($transcript->analysis);
 # $newtranscript->translation(clone_Translation($transcript,$newtranscript));
  return $newtranscript;
}

sub clone_exon_extended {
  my ($self,$exon) = @_ ;
  $self->throw("Feature must be a Bio::EnsEMBL::Analysis:Tools::GeneBuildUtils::ExonExtended not a ".ref($exon)."\n")
    unless $exon->isa("Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended");
  my $newexon = Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended->new();
  foreach my $sf(@{$exon->get_all_supporting_features}){
    my $newsf = clone_Evidence($sf);
    $newsf->dbID(undef);
    $newexon->add_supporting_features($newsf);
  }
  $newexon->start      ($exon->start);
  $newexon->end        ($exon->end);
  $newexon->phase      ($exon->phase);
  $newexon->end_phase  ($exon->end_phase);
  $newexon->strand     ($exon->strand);
  $newexon->slice      ($exon->slice);
  $newexon->stable_id  ($exon->stable_id);
  $newexon->analysis   ($exon->analysis);
  $newexon->biotype    ($exon->biotype);
  $newexon->transcript ($exon->transcript);
  $newexon->ev_set     ($exon->ev_set);
  $newexon->cluster    ($exon->cluster);
  return $newexon;
}


=head2 get_all_evidence_sets

Name : get_all_evidence_sets
Arg[1] : none
Function : Returns hash of evidence_set_names($key) and biotypes
Returnval : Hashref $hash{evidence_set_name} = @\biotypes

=cut


sub get_all_evidence_sets {
  my ($self) = @_ ;
  return $self->{evidence_sets};
}




=head2 get_genes_by_evidence_set

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





=head2 get_genes_by_biotype

Name : get_genes_by_biotype($arg)
Arg[1] : String
Function : Returns all genes of a specific biotype (like all simgw_100 or genscan genes)
Returnval : Arrayref of Bio::EnsEMBL::Gene objects

=cut

sub get_genes_by_biotype {
  my ($self, $biotype ) = @_ ;
  return ${ $self->{all_genes_href}}{$biotype};
}


=head2 get_biotypes_of_evidence_set

Name : get_biotypes_of_evidence_set($arg)
Arg[1] : String
Function : Returns all biotypes specific evidence set (like simgw_100 , genscan )
Returnval : Arrayref of Strings

=cut


sub get_biotypes_of_evidence_set {
  my ($self, $ev_set ) = @_ ;
  return ${ $self->{evidence_sets}}{$ev_set};
}





=head2 verbose

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::TranscriptConsensus
  Arg [2]   : Varies, tends to be boolean, a string, a arrayref or a hashref
  Function  : Getter/Setter for config variables
  Returntype: again varies
  Exceptions:
  Example   :

=cut

#Note the function of these variables is better described in the
#config file itself Bio::EnsEMBL::Analysis::Config::GeneBuild::TranscriptConsensus


sub verbose {
    my ($self, $arg) = @_ ;
  if($arg) {
    $self->{'verbose'} = $arg ;
  }
  return $self->{'verbose'} ;
}

sub filter_singletons {
    my ($self, $arg) = @_ ;
  if($arg) {
    $self->{'filter_singletons'} = $arg ;
  }
  return $self->{'filter_singletons'} ;
}

sub filter_non_consensus {
    my ($self, $arg) = @_ ;
  if($arg) {
    $self->{'filter_non_consensus'} = $arg ;
  }
  return $self->{'filter_non_consensus'} ;
}

sub filter_ests {
    my ($self, $arg) = @_ ;
  if($arg) {
    $self->{'filter_ests'} = $arg ;
  }
  return $self->{'filter_ests'} ;
}

sub add_utr {
    my ($self, $arg) = @_ ;
  if($arg) {
    $self->{'add_utr'} = $arg ;
  }
  return $self->{'add_utr'} ;
}

sub min_consensus {
    my ($self, $arg) = @_ ;
  if($arg) {
    $self->{'min_consensus'} = $arg ;
  }
  return $self->{'min_consensus'} ;
}

sub utr_penalty {
    my ($self, $arg) = @_ ;
  if($arg) {
    $self->{'utr_penalty'} = $arg ;
  }
  return $self->{'utr_penalty'} ;
}

sub end_exon_penalty {
    my ($self, $arg) = @_ ;
  if($arg) {
    $self->{'end_exon_penalty'} = $arg ;
  }
  return $self->{'end_exon_penalty'} ;
}

sub est_overlap_penalty {
    my ($self, $arg) = @_ ;
  if($arg) {
    $self->{'est_overlap_penalty'} = $arg ;
  }
  return $self->{'est_overlap_penalty'} ;
}

sub short_intron_penalty {
    my ($self, $arg) = @_ ;
  if($arg) {
    $self->{'short_intron_penalty'} = $arg ;
  }
  return $self->{'short_intron_penalty'} ;
}

sub short_exon_penalty {
    my ($self, $arg) = @_ ;
  if($arg) {
    $self->{'short_exon_penalty'} = $arg ;
  }
  return $self->{'short_exon_penalty'} ;
}

sub good_percent {
    my ($self, $arg) = @_ ;
  if($arg) {
    $self->{'good_percent'} = $arg ;
  }
  return $self->{'good_percent'} ;
}

sub good_biotype {
    my ($self, $arg) = @_ ;
  if($arg) {
    $self->{'good_biotype'} = $arg ;
  }
  return $self->{'good_biotype'} ;
}

sub small_biotype {
    my ($self, $arg) = @_ ;
  if($arg) {
    $self->{'small_biotype'} = $arg ;
  }
  return $self->{'small_biotype'} ;
}

sub bad_biotype {
    my ($self, $arg) = @_ ;
  if($arg) {
    $self->{'bad_biotype'} = $arg ;
  }
  return $self->{'bad_biotype'} ;
}


sub rnaseq_introns  {
    my ($self, $arg) = @_ ;
  if($arg) {
    $self->{'rnaseq_introns'} = $arg ;
  }
  return $self->{'rnaseq_introns'} ;
}



1;
