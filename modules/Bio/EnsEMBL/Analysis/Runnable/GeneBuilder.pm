package Bio::EnsEMBL::Analysis::Runnable::GeneBuilder;

use strict;  
use warnings;
use vars   qw(@ISA);

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange ); 

use Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils qw(transfer_supporting_evidence Exon_info);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils qw(Transcript_info print_Transcript_and_Exons);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::GeneUtils qw (Gene_info);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);



sub new{
  my ($class,@args) = @_;
  #print "In GeneBuilder constructor with super class" . ref($class) . "\n";
  my $self = $class->SUPER::new(@args);
  my($genes, $blessed_biotypes, $max_transcript_number,
     $min_short_intron_len, $max_short_intron_len,$output_biotype) = 
    rearrange([qw(GENES BLESSED_BIOTYPES MAX_TRANSCRIPT_PER_CLUSTER
                  MIN_SHORT_INTRON_LEN MAX_SHORT_INTRON_LEN OUTPUT_BIOTYPE)], @args);
  print "HERE ARE THE ARGS: ", join("  ",@args),"\n";
  print "HERE I AM: ", $min_short_intron_len,"\t", $max_short_intron_len,"\n";

  ########################
  ###SETTING DEFAULTS#####
  $self->max_transcript_number(10);
  #$self->min_short_intron_len(7);
  #$self->max_short_intron_len(15);
  $self->output_biotype($self->analysis->logic_name);
  #########################

  $self->input_genes($genes);
  $self->blessed_biotypes($blessed_biotypes);
  $self->max_transcript_number($max_transcript_number);
  $self->min_short_intron_len($min_short_intron_len);
  $self->max_short_intron_len($max_short_intron_len);
  $self->output_biotype($output_biotype);

  #data sanity
  warning("Strange running the Genebuilder without any genes") 
    if(!$genes || scalar(@$genes) == 0);

  return $self;
}

sub input_genes{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'input_genes'} = $arg;
  }
  return $self->{'input_genes'};
}

sub blessed_biotypes{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'blessed_biotypes'} = $arg;
  }
  return $self->{'blessed_biotypes'};
}

sub max_transcript_number{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'max_transcript_number'} = $arg;
  }
  return $self->{'max_transcript_number'};
}

sub min_short_intron_len{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'min_short_intron'} = $arg;
  }
  return $self->{'min_short_intron'};
}

sub max_short_intron_len{
  my ($self, $arg) = @_;
  if(defined $arg){
    $self->{'max_short_intron'} = $arg;
  }
  return $self->{'max_short_intron'};
}

sub output_biotype{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'output_biotype'} = $arg;
  }
  return $self->{'output_biotype'};
}

sub run{
  my ($self) = @_;
  my $transcripts = $self->get_Transcripts;
  if(!$transcripts || @$transcripts == 0){
    print "GeneBuilder seems to have no transcripts to cluster\n";
    return;
  }
  #cluster transcripts
  #print "Have ".@$transcripts." transcripts to build from\n";
  my $transcript_clusters = $self->cluster_Transcripts($transcripts);
  #print "Have ".@$transcript_clusters." transcript clusters\n";
  #prune transcripts
  my $pruned_transcripts = $self->prune_Transcripts($transcript_clusters);
  #print "Have ".@$pruned_transcripts." pruned transcripts\n";
  #first gene cluster
  my $initial_genes = $self->cluster_into_Genes($pruned_transcripts);
  #print "Have ".@$initial_genes." initial gene structures\n";
  #prune redundant cds
  $self->prune_redundant_CDS($initial_genes);
  #prune redundant transcripts
  $self->prune_redundant_transcripts($initial_genes);
  #select best transcripts
  my $best_transcripts = $self->select_best_transcripts($initial_genes);
  #print "Have ".@$best_transcripts." best transcripts\n";
  #final cluster of genes
  my $final_genes = $self->cluster_into_Genes($best_transcripts);
  my $pruned_genes = $self->make_shared_exons_unique($final_genes);
  $self->output($pruned_genes);
  #print "Have ".@$final_genes." final gene clusters\n";
  #return output
}

sub get_Transcripts{
  my ($self, $transcripts) = @_;
  if($transcripts){
    $self->{'transcripts'} = $transcripts;
  }
  if(!$self->{'transcripts'}){
    foreach my $gene(@{$self->input_genes}){
      foreach my $transcript(@{$gene->get_all_Transcripts}){
        push(@{$self->{'transcripts'}}, $transcript);
      }
    }
  }
  return $self->{'transcripts'};
}

sub cluster_Transcripts{
  my ($self, $transcripts) = @_;
  my @forward;
  my @reverse;
  foreach my $transcript(@{$transcripts}){
    if($transcript->strand == 1){
      push(@forward, $transcript);
    }else{
      push(@reverse, $transcript);
    }
  }
  my $forward_clusters = $self->cluster_Transcripts_by_genomic_range(\@forward);
  my $reverse_clusters = $self->cluster_Transcripts_by_genomic_range(\@reverse);
  my @clusters;
  push(@clusters, @{$forward_clusters}) if($forward_clusters);
  push(@clusters, @{$reverse_clusters}) if($reverse_clusters);
  return \@clusters;
}

sub cluster_Transcripts_by_genomic_range{
  my ($self, $transcripts) = @_;
  if(!$transcripts || @$transcripts == 0){
    return undef;
  }
  my @transcripts = sort { $a->start <=> $b->start ? $a->start <=> $b->start : $b->end <=> $a->end } @$transcripts;

  my @clusters;
  my $cluster_count = 0;
  my @cluster_starts;
  my @cluster_ends;
  
  my $cluster = Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster->new();
  $cluster->put_Transcripts([$transcripts[0]], 0);

  $cluster_starts[$cluster_count] = $transcripts[0]->start;
  $cluster_ends[$cluster_count] = $transcripts[0]->end;

  push(@clusters, $cluster);
  for (my $c=1; $c<=$#transcripts; $c++){
    #if the current transcript end is not less than the current cluster start
    #or the current transcript start is not greater than the current cluster end
    #add the current transcript to the current cluster
    #readjust the coords
    if ( !( $transcripts[$c]->end < $cluster_starts[$cluster_count] ||
            $transcripts[$c]->start > $cluster_ends[$cluster_count] ) ){
      $cluster->put_Transcripts([$transcripts[$c]], 0);
      
      # re-adjust size of cluster
      if ($transcripts[$c]->start < $cluster_starts[$cluster_count]) {
        $cluster_starts[$cluster_count] = $transcripts[$c]->start;
      }
      if ( $transcripts[$c]->end > $cluster_ends[$cluster_count]) {
        $cluster_ends[$cluster_count] =  $transcripts[$c]->end;
      }
    }else{
      #here the transcripts don't overlap so a new cluster is created
      #this involves creating a new object, incrementing the cluster count
      #and setting up new cluster start and ends
      $cluster_count++;
      $cluster = Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster->new();
      $cluster->put_Transcripts([$transcripts[$c]], 0);
      $cluster_starts[$cluster_count] = $transcripts[$c]->start;
      $cluster_ends[$cluster_count]   = $transcripts[$c]->end;
      
      # store it in the list of clusters
      push(@clusters,$cluster);
    }
  }
  return \@clusters;
}

sub prune_Transcripts {
  my ($self, $transcript_clusters) = @_;
  my @newtran;

  my %blessed_genetypes = %{$self->blessed_biotypes};;

  my $cluster_count = 0;
  
 CLUSTER:
  foreach my $transcript_cluster ( @$transcript_clusters ){
    $cluster_count++;
    #print STDERR "Cluster $cluster_count\n";
    my $mytranscripts = $transcript_cluster->get_Transcripts;
   
    ########################
    #
    # sort the transcripts
    #
    ########################
    my @transcripts = @{$self->bin_sort_transcripts( $mytranscripts )};
  
    ##############################
    #
    # deal with single exon genes
    #
    ##############################

    # do we really just want to take the first transcript only? What about supporting 
    # evidence from other transcripts?
    # also, if there's a very long single exon gene we will lose any underlying 
    # multi-exon transcripts
    # this may increase problems with the loss of valid single exon genes as mentioned
    # below. it's a balance between keeping multi exon transcripts and losing 
    # single exon ones
   
    my $maxexon_number = 0;
    foreach my $t (@transcripts){
      if ( scalar(@{$t->get_all_Exons}) > $maxexon_number ){
        $maxexon_number = scalar(@{$t->get_all_Exons});
      }
    }
    if ($maxexon_number == 1){
      # take the longest:
      @transcripts = map { $_->[1] } sort { $b->[0]->length <=> $a->[0]->length } map{ [ $_->start_Exon, $_ ] } @transcripts;
      my $tran = shift( @transcripts );
      push (@newtran, $tran);
      my @es = @{$tran->get_all_Exons};
      my $e  = $es[0];
      foreach my $transcript (@transcripts){
        # make sure we keep it if it's blessed
        if(exists $blessed_genetypes{$transcript->type}){
          push(@newtran, $transcript);
        }
        else{
          foreach my $exon ( @{$transcript->get_all_Exons} ){
            transfer_supporting_evidence($exon, $e);
          }
        }
      }
      next CLUSTER;
    }
    # otherwise we need to deal with multi exon transcripts and reject duplicates.
    # links each exon in the transcripts of this cluster with a hash of other exons 
    #it is paired with
    my %pairhash;

    # allows retrieval of exon objects by exon->id - convenience
    my %exonhash;

    # keep track of single exon transcripts (automatically rejected if the "top" 
    # transcript is multi exon), 
    # so we can check their supporting evidence at the very end.
    my %single_exon_rejects;    

    ##############################
    #
    # prune redundant transcripts
    #
    ##############################
  TRANSCRIPT:
    foreach my $tran (@transcripts) {
      my @exons = @{$tran->get_all_Exons};

      my $i     = 0;
      my $found = 1;
      # if this transcript has already been seen, this
      # will be used to transfer supporting evidence
      my @evidence_pairs;

      # 10.1.2002 VAC we know there's a potential problem here - single exon 
      # transcripts which are in a  cluster where the longest transcriopt has > 1 
      # exon are not going to be considered in this loop, so they'll always be 
      # marked "transcript already seen" How to sort them out? If the single exon 
      # overlaps an exon in a multi exon transcript then by our rules it probably 
      # ought to be rejected the same way transcripts with shared exon-pairs are.
      # 27/5/2003 VAC But the supporting evidence should be transferred if we can! 
      # Tough one.
      # more of a problem is that if the transcript with the largest number of exons 
      # is really a single exon with frameshifts, it will get rejected here based on
      # intron size but in addition any valid non-frameshifted single exon transcripts
      # will get rejected - which is definitely not right.We need code to represent 
      # frameshifted exons more sensibly so the frameshifted one doesn't get through 
      #the check for single exon genes above.

    EXONS:
      for ($i = 0; $i < $#exons; $i++) {
        my $foundpair = 0;
        my $exon1 = $exons[$i];
        my $exon2 = $exons[$i+1];
                
        # Only count introns > 50 bp as real introns
        my $intron;
        if ($exon1->strand == 1) {
          $intron = abs($exon2->start - $exon1->end - 1);
        }
        else {
          $intron = abs($exon1->start - $exon2->end - 1);
        }
        
        if ($intron < $self->max_short_intron_len && 
            $intron > $self->min_short_intron_len ) {
          print STDERR "Intron too short: $intron bp. Transcript will be rejected\n";
          $foundpair = 1;       
          # this pair will not be compared with other transcripts
        } else {
          # go through the exon pairs already stored in %pairhash. 
          # If there is a pair whose exon1 overlaps this exon1, and 
          # whose exon2 overlaps this exon2, then these two transcripts are paired
          
          foreach my $first_exon_id (keys %pairhash) {
            my $first_exon = $exonhash{$first_exon_id};
            foreach my $second_exon_id (keys %{$pairhash{$first_exon}}) {
              my $second_exon = $exonhash{$second_exon_id};
              if ( $exon1->overlaps($first_exon) && $exon2->overlaps($second_exon) ) {
                $foundpair = 1;
                # eae: this method allows a transcript to be covered by exon pairs
                # from different transcripts, rejecting possible
                # splicing variants. Needs rethinking
                
                # we put first the exon from the transcript being tested:
                push( @evidence_pairs, [ $exon1 , $first_exon  ] );
                push( @evidence_pairs, [ $exon2 , $second_exon ] );
                
                transfer_supporting_evidence($exon1, $first_exon);
                transfer_supporting_evidence($first_exon, $exon1);
                transfer_supporting_evidence($exon2, $second_exon);
                transfer_supporting_evidence($second_exon, $exon2);
              }
            }
          }
        }
        if ($foundpair == 0) {  
          # ie this exon pair does not overlap with a pair yet found in another 
          # transcript
          $found = 0;           
          # ie currently this transcript is not paired with another
          # store the exons so they can be retrieved by id
          $exonhash{$exon1} = $exon1;
          $exonhash{$exon2} = $exon2;
          # store the pairing between these 2 exons
          $pairhash{$exon1}{$exon2} = 1;
        }
      }                         # end of EXONS

      # decide whether this is a new transcript or whether it has already been seen
      # if it's blessed, we keep it and there's nothing more to do
      if(exists $blessed_genetypes{$tran->type}){
        push (@newtran, $tran);
      }
      elsif ($found == 0) {
        push(@newtran,$tran);
        @evidence_pairs = ();
      } elsif ($found == 1 && $#exons == 0){
        # save the transcript and check though at the end to see if we can transfer 
        # supporting evidence; if we try it now we may not (yet) have any exons that 
        # overlap in %exonhash
        $single_exon_rejects{$tran} = $tran;
      } else {
        if ( $tran == $transcripts[0] ){
          print STDERR "Strange, this is the first transcript in the cluster!\n";
        }
                
        ## transfer supporting feature data. We transfer it to exons
        foreach my $pair ( @evidence_pairs ){
          my @pair = @$pair;
        
          # first in the pair is the 'already seen' exon
          my $source_exon = $pair[0];
          my $target_exon = $pair[1];
        
          transfer_supporting_evidence($source_exon, $target_exon)
        }
      }
    } # end of this transcript

    # check to see if we can transfer evidence from rejected single exon transcripts
    foreach my $reject(values %single_exon_rejects){
      my @exons = @{$reject->get_all_Exons};

      # is there any supporting evidence?
      my @sf = @{$exons[0]->get_all_supporting_features};

      foreach my $stored_exon(values %exonhash){
        if($exons[0]->overlaps($stored_exon)){
          # note that we could end up with bizarre situations of a single exon 
          # transcript overlapping two exons in a multi exon transcript, so the 
          # supporting evidence would be transferred in entirety to both exons.
          transfer_supporting_evidence($exons[0], $stored_exon);
        
        }
      }
    }
  } #end CLUSTER
  return \@newtran;
} 


sub cluster_into_Genes{
  my ($self, $transcripts_unsorted) = @_;
  my $num_trans = scalar(@$transcripts_unsorted);
  my @transcripts = sort { $a->start <=> $b->start ? $a->start <=> $b->start  : $b->end <=> $a->end } @$transcripts_unsorted;

  my @clusters;
  # clusters transcripts by whether or not any exon overlaps with an exon in 
  # another transcript (came from original prune in GeneBuilder)
  foreach my $tran (@transcripts) {
    my @matching_clusters;
  CLUSTER: 
    foreach my $cluster (@clusters) {
      foreach my $cluster_transcript (@$cluster) {
        if ($tran->end  >= $cluster_transcript->start &&
            $tran->start <= $cluster_transcript->end) {
          
          foreach my $exon1 (@{$tran->get_all_Exons}) {
            foreach my $cluster_exon (@{$cluster_transcript->get_all_Exons}) {
              if ($exon1->overlaps($cluster_exon) 
                  && $exon1->strand == $cluster_exon->strand) {
                push (@matching_clusters, $cluster);
                next CLUSTER;
              }
            }
          }
        }
      }
    }

    if (scalar(@matching_clusters) == 0) {
      my @newcluster;
      push(@newcluster,$tran);
      push(@clusters,\@newcluster);
    } elsif (scalar(@matching_clusters) == 1) {
      push @{$matching_clusters[0]}, $tran;
    } else {
      # Merge the matching clusters into a single cluster
      my @new_clusters;
      my @merged_cluster;
      foreach my $clust (@matching_clusters) {
        push @merged_cluster, @$clust;
      }
      push @merged_cluster, $tran;
      push @new_clusters,\@merged_cluster;
      # Add back non matching clusters
      foreach my $clust (@clusters) {
        my $found = 0;
      MATCHING: 
        foreach my $m_clust (@matching_clusters) {
          if ($clust == $m_clust) {
            $found = 1;
            last MATCHING;
          }
        }
        if (!$found) {
          push @new_clusters,$clust;
        }
      }
      @clusters =  @new_clusters;
    }
  }
  # safety and sanity checks
  $self->check_Clusters(scalar(@transcripts), \@clusters);
  
  # make and store genes
  my @genes;
  foreach my $cluster(@clusters){
    my $count = 0;
    my $gene = new Bio::EnsEMBL::Gene;
    $gene->biotype($self->output_biotype);
    foreach my $transcript (@$cluster){
      $gene->add_Transcript($transcript);
    }
    push( @genes, $gene );
  }
  return \@genes;
}


sub bin_sort_transcripts{

  my ($self,$transcripts) = @_;

  my %lengths;
  my $max_num_exons = 0;

  # sizehash holds transcript length - based on sum of exon lengths
  my %sizehash;

  # orfhash holds orf length - based on sum of translateable exon lengths
  my %orfhash;

  my %tran2orf;
  my %tran2length;

  # keeps track of to which transcript(s) this exon belong
  my %exon2transcript;

  foreach my $tran (@$transcripts) {

    # keep track of number of exons in multiexon transcripts
    my @exons = @{ $tran->get_all_Exons };
    if(scalar(@exons) > $max_num_exons){
      $max_num_exons = scalar(@exons);
    }

    # total exon length
    my $length = 0;
    foreach my $e ( @exons ){
      $length += $e->end - $e->start + 1;
      push ( @{ $exon2transcript{ $e } }, $tran );
    }
    $sizehash{$tran} = $length;
    $tran2length{ $tran } = $length;

    # now for ORF length
    $length = 0;
    foreach my $e(@{$tran->get_all_translateable_Exons}){
      $length += $e->end - $e->start + 1;
    }
    $tran2orf{ $tran } = $length;
    push(@{$orfhash{$length}}, $tran);
  }

  # VAC 15/02/2002 sort transcripts based on total exon length - this
  # introduces a problem - we can (and have) masked good transcripts
  # with long translations in favour of transcripts with shorter
  # translations and long UTRs that are overall slightly longer. This
  # is not good.
  # better way? hold both total exon length and length of translateable exons. 
  # Then sort: long translation + UTR > long translation no UTR > short translation 
  # + UTR > short translation no UTR
  ##########
  my @transcripts = ();
  ##########

  my @orflengths = sort {$b <=> $a} (keys %orfhash);

  # strict sort by translation length is just as wrong as strict sort by UTR length
  # bin translation lengths - 4 bins (based on 25% length diff)? 10 bins 
  # (based on 10%)?
  my %orflength_bin;
  my $numbins = 4;
  my $currbin = 1;

  ORF: foreach my $orflength(@orflengths){
    last if $currbin > $numbins;  

    if ( $orflengths[0] == 0 ) { 
        push(@{$orflength_bin{$currbin}}, @{$orfhash{$orflength}});
       # $orflengths[0] =1 ; 
       next ORF ; 
    }
    my $percid = ($orflength*100)/$orflengths[0]; 
    if ($percid > 100) { 
      $percid = 100; 
    }
    my $currthreshold = $currbin * (100/$numbins);
    $currthreshold    = 100 - $currthreshold;

    if($percid <$currthreshold) { 
      $currbin++; 
    } 
    my @tmp = @{$orfhash{$orflength}};
    push(@{$orflength_bin{$currbin}}, @{$orfhash{$orflength}});
  }

  # now, foreach bin in %orflengthbin, sort by exonlength
  $currbin = 1;
 EXONLENGTH_SORT:
  while( $currbin <= $numbins){
    if(!defined $orflength_bin{$currbin} ){
      $currbin++;
      next EXONLENGTH_SORT;
    }   
    my @sorted_transcripts = sort { $sizehash{$b} <=> $sizehash{$a} } 
      @{$orflength_bin{$currbin}};
    push(@transcripts, @sorted_transcripts);
    $currbin++;
  }

  return \@transcripts;
}

sub check_Clusters{
  my ($self, $num_transcripts, $clusters) = @_;
  #Safety checks
  my $ntrans = 0;

  my $cluster_num = 0;

  my %trans_check_hash;
  foreach my $cluster (@$clusters) {
    $ntrans += scalar(@$cluster);

    foreach my $trans (@$cluster) {

      if (defined($trans_check_hash{$trans})) {
        throw("Transcript " . $trans->dbID ." ".$trans->adaptor->dbc->dbname. 
              " added twice to clusters\n");
      }
      $trans_check_hash{$trans} = 1;
    }
    if (!scalar(@$cluster)) {
      throw("Empty cluster");
    }
  }
  if ($ntrans != $num_transcripts) {
    throw("Not all transcripts have been added into clusters $ntrans and " . $num_transcripts. " \n");
  }
  #end safety checks
  return;
}


sub prune_redundant_CDS {
  my ( $self, $genes ) = @_;

  my $nremoved = 0;
  my %blessed_genetypes = %{$self->blessed_biotypes};

  # For each gene
  foreach my $gene (@$genes) {
    my @trans_with_utrs;
    my @trans_without_utrs;
    
    # Separate into ones with UTRs and ones without
    foreach my $trans (@{$gene->get_all_Transcripts}) {
      $trans->sort;
      my @exons = @{$trans->get_all_Exons};
      if ($trans->translation) {
        if ($trans->translation->start_Exon == $exons[0] &&
            $trans->translation->start == 1 &&
            $trans->translation->end_Exon == $exons[$#exons] &&
            $trans->translation->end == $exons[$#exons]->length) {
          push @trans_without_utrs, $trans;
        } else {
          push @trans_with_utrs, $trans;
        }
      } else {
        warning("No translation for transcript");
      }
    }
    
    # Generate CDS exons just once (saves CPU time)
    my %cds_exon_hash;
    foreach my $trans (@{$gene->get_all_Transcripts}) {
      if ($trans->translation) {
        my @cds_exons = @{$trans->get_all_translateable_Exons};
        $cds_exon_hash{$trans} = \@cds_exons;
      }
    }
      
    # Compare CDSs of ones with UTRs to CDSs of ones without
    foreach my $utr_trans (@trans_with_utrs) {
      my $utr_trans_cds_exons = $cds_exon_hash{$utr_trans};
      
    CDS_TRANS: foreach my $cds_trans (@trans_without_utrs) {
        my $cds_trans_cds_exons = $cds_exon_hash{$cds_trans};
        if (scalar(@$cds_trans_cds_exons) != scalar(@$utr_trans_cds_exons)) {
          # print "Different numbers of exons\n";
          next CDS_TRANS;
        }
        my @cds_exons = @$cds_trans_cds_exons;
        foreach my $utr_trans_exon (@$utr_trans_cds_exons) {
          my $cds_trans_exon = shift @cds_exons;
          if ($cds_trans_exon->start     != $utr_trans_exon->start  ||
              $cds_trans_exon->end       != $utr_trans_exon->end    ||
              $cds_trans_exon->strand    != $utr_trans_exon->strand) {
            next CDS_TRANS;
          }
        }
        if (  exists $blessed_genetypes{$cds_trans->biotype} && 
              ! exists $blessed_genetypes{$utr_trans->biotype}) {
          # Hack to make sure transcript gets through if its like a blessed one
          $utr_trans->biotype($cds_trans->biotype);
        }
        $nremoved++;
        $self->remove_transcript_from_gene($gene,$cds_trans);
      }
    }
    #Get the number of transcript in the @trans_with_utrs 
    my $num_of_utr_trans = scalar(@trans_with_utrs);
    # Extra paranoid test to check if any similarity overlaps targetted models and 
    # get rid of them
  }
}


sub remove_transcript_from_gene {
  my ($self, $gene, $trans_to_del)  = @_;
  my @newtrans;

  foreach my $trans (@{$gene->get_all_Transcripts}) {
    if ($trans != $trans_to_del) {
      push @newtrans,$trans;
    }
  }
  $gene->{_transcript_array} = [];
  foreach my $trans (@newtrans) {
    $gene->add_Transcript($trans);
  }
  return scalar(@newtrans);
}

sub prune_redundant_transcripts {
  my ($self,$genes) = @_;

  my $nremoved = 0;
   my %blessed_genetypes = %{$self->blessed_biotypes};
  # For each gene
  foreach my $gene (@$genes) {
    # sort transcripts by total cdna length (apologies for the complexity of this line
    #  its done so as only to calculate cdna length once for each transcript

    my @sorted_transcripts = map { $_->[1] } sort { $b->[0] <=> $a->[0] } map { [$_->length, $_] } @{$gene->get_all_Transcripts};

    # so now longer UTR transcripts (with same introns will come first)
    # Now generate the sorted (low to high) exon arrays and put in a hash (for speed)
    # and cds start and end (genomic) hashes (for speed again)
    my %sorted_exon_hash;
    my %cds_start_hash;
    my %cds_end_hash;

    foreach my $trans (@sorted_transcripts) {
      my @exons = sort {$a->start <=> $b->start} @{$trans->get_all_Exons};
      $sorted_exon_hash{$trans} = \@exons;
      if ($trans->translation) {
        $cds_start_hash{$trans} = $trans->coding_region_start;
        $cds_end_hash{$trans}   = $trans->coding_region_end;
      }
    }

    # We're going to generate a list of transcripts to remove from the gene
    # and then remove them at the end
    my @trans_to_remove;

    # for transcript

    for (my $i=0; $i<scalar(@sorted_transcripts); $i++) {
      my $trans = $sorted_transcripts[$i];
      next if (!defined($trans));
      next if (!$trans->translation);
      my @exons = @{$sorted_exon_hash{$trans}};
      #   for every other transcript
    COMPTRANS:
      for (my $j=$i+1; $j<scalar(@sorted_transcripts); $j++) {
        my $comp_trans = $sorted_transcripts[$j];
        next if (!defined($comp_trans));
        next if (!$comp_trans->translation);
        next if ($trans == $comp_trans);
        next if ($cds_start_hash{$trans} != 
                 $cds_start_hash{$comp_trans});
        next if ($cds_end_hash{$trans} != $cds_end_hash{$comp_trans});
        my @comp_exons = @{$sorted_exon_hash{$comp_trans}};
        if (scalar(@comp_exons) != scalar(@exons)){
          next;
        }

        # special case for single exon
        if (scalar(@exons) == 1) {
          if (  exists $blessed_genetypes{$comp_trans->biotype} && 
                ! exists $blessed_genetypes{$trans->biotype}) {
            #hack to ensure blessed biotypes are maintained
            $trans->biotype($comp_trans->biotype);
          }

          $sorted_transcripts[$j] = undef;
          push @trans_to_remove, $comp_trans;
        } else {
          for (my $k=0; $k<scalar(@exons) - 1; $k++) {
            if ($exons[$k]->end != $comp_exons[$k]->end ||
                $exons[$k+1]->start != $comp_exons[$k+1]->start) {
              #means exons are unique
              next COMPTRANS;
            }
          }
          if (  exists $blessed_genetypes{$comp_trans->biotype} && 
                ! exists $blessed_genetypes{$trans->biotype}) {
            # Hack to make sure transcript gets through if its like a blessed one
            $trans->biotype($comp_trans->biotype);
          }
          $sorted_transcripts[$j] = undef;
          push @trans_to_remove, $comp_trans;
        }
      }
    }
    foreach my $trans (@trans_to_remove) {
      $nremoved++;
      $self->remove_transcript_from_gene($gene,$trans);
    }
  }
}

sub select_best_transcripts{
  my ( $self, $genes ) = @_;
  my @selected_transcripts;

  # make sure we don't discard blessed transcripts
  my %blessed_genetypes = %{$self->blessed_biotypes};
 GENE:
  foreach my $gene ( @$genes ){
    # sort the transcripts, get the longest CDS + UTR first (like in 
    # prune_Transcripts(); blessed transcripts come at the top )
    my $sorted_transcripts = 
      $self->bin_sort_transcripts( $gene->get_all_Transcripts );
    my $count = 0;
  TRAN:
    foreach my $transcript( @$sorted_transcripts ){
      $count++;
      unless (exists $blessed_genetypes{$transcript->biotype}){
        next TRAN if ($count > $self->max_transcript_number);
      }
      push ( @selected_transcripts, $transcript );
    }
  }
  return \@selected_transcripts;
}


sub make_shared_exons_unique{
  my ($self, $genes) = @_;
  my @pruned;
  foreach my $gene(@$genes){
    my $new_gene = $self->prune_Exons($gene);
    push(@pruned, $new_gene);
  }
  return \@pruned;
}

sub prune_Exons{
  my ($self, $gene) = @_;
  my @unique_Exons;
  foreach my $tran (@{$gene->get_all_Transcripts}) {
    my @newexons;
    foreach my $exon (@{$tran->get_all_Exons}) {
      my $found;
      #always empty
    UNI:foreach my $uni (@unique_Exons) {
        if ($uni->start  == $exon->start  &&
            $uni->end    == $exon->end    &&
            $uni->strand == $exon->strand &&
            $uni->phase  == $exon->phase  &&
            $uni->end_phase == $exon->end_phase
           ) {
          $found = $uni;
          last UNI;
        }
      }
      if (defined($found)) {
        push(@newexons,$found); 
        if ( $tran->translation ) { 
          if ($exon == $tran->translation->start_Exon){
            $tran->translation->start_Exon($found);
          }
          if ($exon == $tran->translation->end_Exon){
            $tran->translation->end_Exon($found);
          } 
        }
      } else {
        push(@newexons,$exon);
        push(@unique_Exons, $exon);
      }
    }          
    $tran->flush_Exons;
    foreach my $exon (@newexons) {
      $tran->add_Exon($exon);
    }
  }
  return $gene;
}

1;
