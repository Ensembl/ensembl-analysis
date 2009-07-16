=pod

=head1 NAME Bio::EnsEMBL::Analysis::Runnable::TranscriptConsensus

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::TranscriptConsensus->new
     (
     -query            => $query,
     -analysis        => $analysis,
     -all_genes       => \%biotypes_to_genes , # ref to $hash{biotype_of_gene} = @all_genes_of_this_biotype
     -evidence_sets   => $evidence_sets ,
     -dnadb           => $db ,
     )

=head1 DESCRIPTION

TranscriptConsensus is an extension of TranscriptCoalescer that combines protein coding and est
transcripts to identify the best supported transcript models.
The initial gene sets are clustered, then collapsed into a non-redundant set of exons
and introns which are assigned scores according to the amount of supporting evidence they have.
The similarity genes have UTR added using the est transcripts and the resulting models are assigned
a score by summing the individual exon and intron scores.
The transcripts are sorted by score and the highest scoring models are made into
gene objets and written to the TranscriptCoalescer database.

=head1 CONTACT

Post questions to the EnsEMBL developer list: <ensembl-dev@ebi.ac.uk>

=cut


package Bio::EnsEMBL::Analysis::Runnable::TranscriptConsensus;

use strict;
use warnings; 

use Bio::EnsEMBL::Utils::Exception qw(throw warning info);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Runnable::TranscriptCoalescer;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ExonCluster;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw (count_non_canonical_splice_sites are_phases_consistent Transcript_info) ;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::EvidenceUtils qw ( clone_Evidence );
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils qw( clone_Translation );
use Bio::EnsEMBL::Analysis::Config::GeneBuild::TranscriptCoalescer;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::TranscriptConsensus;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::CollapsedCluster;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::TranscriptCoalescer);



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
  my ($self) = @_ ; 
  my $count =1;
  my $genes_by_strand;

  my @allgenes = @{ $self->get_genes_by_evidence_set('simgw') }  ;
  push @allgenes , @{ $self->get_genes_by_evidence_set('est') }  ;
  
  if ($FILTER_SINGLETONS or   $FILTER_NON_CONSENSUS){
    #@allgenes = @{$self->filter_genes(\@allgenes)};
    my @tmp;
    my ($clusters, $non_clusters) = cluster_Genes(\@allgenes, $self->get_all_evidence_sets ) ;
    foreach my $cluster(@$clusters, @$non_clusters){
      my @genes = $cluster->get_Genes;
      @genes = @{$self->filter_genes(\@genes)};
      push(@tmp, @genes);
    }
    @allgenes = @tmp;
  }

  my ($clusters, $non_clusters) = cluster_Genes(\@allgenes, $self->get_all_evidence_sets ) ;
  
  # create a hash of genes by strand to use when looking for overlapping genes
  # useful to look at ALL genes not just those in the cluster to prevent cluster joining
  foreach my $gene (@allgenes){
    push @{$genes_by_strand->{$gene->strand}}, $gene;
  }

  push @$clusters, @$non_clusters if (scalar(@$non_clusters) > 0 ) ;
  my $found;
  $count = 0;
  foreach my $cluster (@$clusters){
    my $simgw;
    # cluster has to contain at least one similarity gene to be worth continuing with
    my @genes = $cluster->get_Genes;
    #foreach my $gene(@genes){
     # my $new_gene = $gene->transform("toplevel");
     # print Gene_info($new_gene)."\n"; #if($gene->start == 6633986 && $gene->end == 6649904);
    #}
    next unless $cluster->get_Genes_by_Set('simgw');
    
    print "\nCluster $count\n" if $VERBOSE;
    print $cluster->start." ".$cluster->end." ".$cluster->strand."\n" if $VERBOSE;
    $count++;
    # collapse the cluster down to make a non redundat set of introns and exons with scores
    my $collapsed_cluster = $self->collapse_cluster($cluster,$genes_by_strand);

    # add UTR to the similarity transcripts and make scores for the transcripts
    my $transcripts = $self->make_transcripts($cluster,$collapsed_cluster);

    # sort the scored transcripts by score and make gene objects out of the best ones
    $self->make_genes($transcripts,$cluster,$collapsed_cluster);
  }
}

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
  my ($self,$cluster,$genes_by_strand) = @_;
  my @genes = $cluster->get_Genes;
  my @exon_clusters = @{$cluster->get_exon_clustering_from_gene_cluster};
  my $collapsed_cluster = Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::CollapsedCluster->new();

  foreach my $gene (@genes){
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

  # once we have added all of them into a nice non-redundant set
  # we need to assign them a score
  foreach my $exon (@{$collapsed_cluster->get_all_exons}){
    my $score = $self->score($collapsed_cluster,$exon,\@genes,'exon');
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
  my ($self,$collapsed_cluster,$feature,$genes,$type) = @_;
  my $total_score;
  my $solexa_score = 0;
  # scores are calculated as the number of exoct matches / the number of overlaps
  # it is complicated by end exons and evidence sets however.
  $total_score = $collapsed_cluster->get_intron_count($feature) if $type eq 'intron';
  $total_score = $collapsed_cluster->get_exon_count($feature) if $type eq 'exon';
  $collapsed_cluster->exon_info($feature) if $VERBOSE;
  print " " if $VERBOSE;
  my $score = 0 ;
  # this way we penalise anything with less than exons than the minumum number of overlapping exons
  my $overlap = $MIN_CONSENSUS ;
  my $est_trans = $MIN_CONSENSUS;
  # get solexa data if it's there
  if ($SOLEXA &&  defined($self->solexa_slice) && $type eq 'exon'){
    #Fetch the slice from the solexa db that corresponds to the exon
    #Get all solexa reads with score above the cutoff ( 150 reccomended ) regardless of strand
    my $sub_slice = $self->solexa_slice->adaptor->fetch_by_region
      (
       'toplevel',
       $self->solexa_slice->seq_region_name,
       $feature->start + $self->solexa_slice->start,
       $feature->end + $self->solexa_slice->start,
       1);
    $self->throw("Sub slice not found " . $feature->start . " " .  $feature->end . " " .  $feature->strand . "\n")
      unless $sub_slice;
    # Fetch the features
    $solexa_score = scalar(@{$sub_slice->get_all_DnaAlignFeatures
			       (
				$$TRANSCRIPT_CONSENSUS_DB_CONFIG{"SOLEXA_DB"}->{"BIOTYPE"}[0],
				$SOLEXA_SCORE_CUTOFF)}
			  );
    print  "Solexa reads = $solexa_score\n" if $VERBOSE;
  }
  
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
    print "End exon $est_exons " if $VERBOSE;
    # count the total number of overlapping end exons
    $est_exons = $collapsed_cluster->ev_count_by_terminal_feature($feature,'est');
    $total_score = $collapsed_cluster->count_end_exon($feature);
    $self->throw("No evidence found for exon " . $collapsed_cluster->exon_info($feature) . "\n")
      unless defined($est_exons);
  }
  print "$type " if $VERBOSE;
  print "TOTAL $total_score OVERLAP $overlap " if $VERBOSE;

  if ($est_exons){
    print "WEIGHTED GENES $est_exons EST TRANS $est_trans " if $VERBOSE;
    # number of exons of the weighted type / number of overlapping trans of the weighted type;
    my $weighted_score = 0;
    my $unweighted_score = 0;
    if ($est_trans){
      $weighted_score = ($est_exons + $solexa_score) / ($est_trans + $solexa_score) ;
      print "weighted score $weighted_score " if $VERBOSE;
    }
    # other componnt of score comes from all the exons
    $unweighted_score =  ($total_score + $solexa_score) / ($overlap + $solexa_score) ;
    print "unweighted score $unweighted_score \t" if $VERBOSE;
    $score = ($weighted_score + $unweighted_score) / 2;
  } else {
    # there are no est exons overlaping this feature score is simple identical / non-identical overlaps
    $score += (($total_score + $solexa_score) / ($overlap + $solexa_score)) ;
    # if the feature overlaps an est transcript but not est exons, add a penalty
    if ($est_trans &&  $est_trans > 1.5 ){
      $score -= $EST_OVERLAP_PENALTY;
    }
  }
  # favor long end exons
  if ($feature->is_terminal_exon && $type eq 'exon' && $feature->number_exons > 1){
    # longest end exon 
    my $longest_end_exon = $collapsed_cluster->get_end_exon($feature)->length;
    my $multiplyer = $feature->length / $longest_end_exon;
    print "End exon penalty ".$feature->length." / ". $longest_end_exon." -.5\t" if $VERBOSE;
    if ($score >= 0){
      $score *= $multiplyer;
    } else {
      $score /= $multiplyer;
    }
    # end exons are penalised to prevent 'spindly' exons
    $score = $score - $END_EXON_PENALTY;
    print " ENDEXON " if $VERBOSE;
  } else {
    # penalty for  internal short exons / introns
      $score-= 0.5 if $feature->length < $SHORT_INTRON_PENALTY && $type eq 'intron';
      $score-= 0.5 if $feature->length < $SHORT_EXON_PENALTY   && $type eq 'exon';
  }
  print "$score \n" if $VERBOSE;
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
  my ($self,$transcripts,$cluster,$collapsed_cluster)= @_;
  my @genes;
  my @final_genes;
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
    if ($ADD_UTR){
      # dont compute the translation just set it
      # need the similarity gene that it was built from
      # so I can set the cds start and end on the transcript
      # if the start is in the same exon in the new transcript
      # just transfer it, if it is over the end of the new exon
      # trim back so that it is in the exon if you see what I mean
      # print stuff
      print "Transcript ".$transcript->start .":".$transcript->end .":".$transcript->strand." ".$transcript->biotype." " if $VERBOSE;
      print scalar(@{$transcript->get_all_Exons})." Exons\n" if $VERBOSE;
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
	  $exon->end_phase( -1 ) 	if $exon->end > $cds_end;
	}
	# find the utr exon that shares the internal boundary with the first simgw exon
	if ($exon->end == $similarity_exons[0]->end){
	  $cds_start_exon = $exon;
	}
	# find the utr exon that shares the internal boundary with the last simgw exon
	if ($exon->start == $similarity_exons[-1]->start){
	  $cds_end_exon = $exon;
	}
      }
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
	# the start is not inside an exon
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
    } else {
      # dont need to mess with translations if I am not adding UTR...
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
    print $transcript->score."\t" if $VERBOSE;
    print "\npeptide weighted score $weighted_score\n\n" if $VERBOSE;
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
  @genes = sort { $b->get_all_Transcripts->[0]->score <=> $a->get_all_Transcripts->[0]->score } @genes;
  my $top_score = $genes[0]->get_all_Transcripts->[0]->score;
  my $bottom_score = $genes[-1]->get_all_Transcripts->[0]->score;
  # prevent the whole division by zero thing
  $bottom_score-= 0.0001 if $bottom_score == 0;
  foreach my $gene (@genes){
    foreach my $transcript ( @{$gene->get_all_Transcripts}){
      # sort out exon phases - transfer them from the similarity gene
      # coding region start / stop is independent of strand so I will swap them if the strand is -ve
      my $biotype = $BAD_BIOTYPE;
      # how about if we include the top 5% of the gene scores
      if ($GOOD_PERCENT){
        my $top = ($transcript->score - $bottom_score );
        my $bottom = ($top_score - $bottom_score);
	if ($bottom && (($top/$bottom) * 100 >= 100 - $GOOD_PERCENT)){
	  $biotype = $GOOD_BIOTYPE;
	}
      } else {
	# we want the single top scoring model
	if ($transcript->score == $top_score){
	  $biotype = $GOOD_BIOTYPE;
	}
      }
      $biotype = $SMALL_BIOTYPE if (scalar(@genes) < $MIN_CONSENSUS);
      $gene->biotype($biotype);
      push @final_genes, $gene if $biotype;
    }
  }
  foreach my $final(@final_genes){
    foreach my $transcript(@{$final->get_all_Transcripts}){
      throw(Transcript_info($transcript)." has inconsistent phases") 
        unless(are_phases_consistent($transcript));
      #modify score to be independent of length i total / number of introns+exons
      #Should always be less than 1
      my $final_score = $transcript->score / (scalar(@{$transcript->get_all_Exons})*2-1);    
      # add score into trancript supporting feature score feild and move
      # coverage into the new covergae field
      foreach my $tsf ( @{$transcript->get_all_supporting_features}) {
	$tsf->hcoverage($tsf->score);
	$tsf->score(sprintf("%.4f", $final_score));
      }	
      print "Transcript ". ($transcript->start + $self->solexa_slice->start -1) .":".
	($transcript->end  + $self->solexa_slice->start-1) . ":".
	  $transcript->strand." ".
	    $final->biotype." " if $VERBOSE && $SOLEXA;
      print "Final score $final_score\n";
    }
  }
  $self->output(\@final_genes);
  return;
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
  my ($self,$cluster,$collapsed_cluster) = @_;
  my @single_exon_genes;
  my @genes = $cluster->get_Genes;
  foreach my $gene (@genes){
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
  my ($self,$cluster,$collapsed_cluster) = @_;
  my @genes = $cluster->get_Genes;
  my @est_genes;
  my @similarity_genes;
  my $longest_similarity = 0;
  my @all_transcripts;
  foreach my $gene (@genes){ 
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
          if ($ADD_UTR ) {  
            unless ( $trans->translation) { 
	      push @start_utr, $trans ;  
            } else { 
              throw("You can't add UTR with non-EST data");   
	    } 
	  } 
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
          if ($ADD_UTR ) {  
            unless ( $trans->translation) { 
	      push @start_utr, $trans ;  
            } else {  
              throw("You can't add UTR with non-EST data"); 
            } 
	  }
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
  my ($self, $exons, $transcript, $collapsed_cluster) = @_;
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
	$collapsed_cluster->exon_info($collapsed_exon) if $VERBOSE;
	$self->throw("No score assigned to this exon\n");
      }
    if ($prev_exon){
      my $intron = $collapsed_cluster->make_intron_from_exons($prev_exon,$exon);
      unless ( $intron ){
	$collapsed_cluster->exon_info($intron) if $VERBOSE;
	$self->throw("No score assigned to this intron\n");
      }
      $score += $collapsed_cluster->intron_score($intron);
      # add a penalty to UTR exons
      $score -= $UTR_PENALTY if $exon->ev_set eq 'est';
    }
    $score += $collapsed_cluster->exon_score($collapsed_exon);
    $score -= $UTR_PENALTY if $exon->ev_set eq 'est';
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
  my ( $self, $genes) = @_;
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
    unless ($FILTER_ESTS){
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
    print "All transcripts look dodgy, not filtering\n" if $VERBOSE;
    return $genes;
  } else {
    for ( my $g = 0 ; $g < scalar(@$genes) ; $g++ ){
      if ( $FILTER_SINGLETONS && $singletons{$g} && $singletons{$g}>= $FILTER_SINGLETONS ){
	print  "Filtering out ".$genes->[$g]->dbID." ".$genes->[$g]->biotype." ".$genes->[$g]->start.":".$genes->[$g]->end.":".$genes->[$g]->strand." " if $VERBOSE;
	print  $singletons{$g}." singleton exons\n" if $VERBOSE;
	next;
	
      }
      if ( $FILTER_NON_CONSENSUS && $non_con{$g} && $non_con{$g} >= $FILTER_NON_CONSENSUS ){
	print  "Filtering out ".$genes->[$g]->dbID." ".$genes->[$g]->biotype." ".$genes->[$g]->start.":".$genes->[$g]->end.":".$genes->[$g]->strand." " if $VERBOSE;
	print  $non_con{$g}." non - consensus\n" if $VERBOSE;	
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
  my ($self,$transcript,$collapsed_cluster) = @_;
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
      $collapsed_cluster->exon_info($exon) if $VERBOSE;
      print " overlaps $coding_exons  coding exons penalty -.5\n" if $VERBOSE;
      $score -= 0.5;
    }
  }
  return $score;
}


sub clone_gene {
  my ($self,$gene) =@_;
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
  my ($self,$transcript) =@_;
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
  my ($self,$exon) = @_;
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

sub solexa_slice{
  my ($self,$slice) = @_;
  if  (defined($slice)) {
    $self->{_solexa} = $slice;
  }
  return $self->{_solexa};
}


1;
