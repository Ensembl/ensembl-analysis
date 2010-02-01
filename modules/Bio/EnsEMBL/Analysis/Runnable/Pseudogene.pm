#
# You may distribute this module under the same terms as perl itself
#

=pod

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Pseudogene

=head1 SYNOPSIS

 my $runnable = Bio::EnsEMBL::Analysis::Runnable::Pseudogene->new
      ( 
       '-genes' => \@_genes,
       '-repeat_features' => \%repeat_blocks, 
        );
    $runnable->run;
    $runnable->output;

Where output returns an array of modified genes.
Repeat blocks is a hash of repeats covering each gene merged into blocks.


=head1 DESCRIPTION

Runnable for PseudogeneDB

Runs tests to identiy pseudogenes:
Calls it a pseudo gene if:

1. All of the introns are frameshifted.
2. Real introns are covered with repeats.
3. Real introns in a two exon gene are covered with a protein feature

Pseudogene takes a Bio::EnsEMBL::Slice object and assesses genes and transcripts for evidence of retrotransposition.
In the case of genes being identified as pseudogenes, the gene objects have their type set to pseudogene and all but the longest transcripts and translations are deleted.
If the gene has 1 or more pseudo transcripts but has other transcritps that are valid the dubious transcripts are removed The resulting gene objects are returned in an array.

=head1 CONTACT

Mail to B<ensembl-dev@ebi.ac.uk>

=cut

package Bio::EnsEMBL::Analysis::Runnable::Pseudogene;

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::SeqFeature;
use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception; 
use Bio::EnsEMBL::Transcript; 
use Bio::EnsEMBL::Gene ;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

# Object preamble - inherits from Bio::EnsEMBL::Root;



=head2 new

  Args       : various
  Description: Runnable constructor
  Returntype : Bio::EnsEMBL::Analysis::Runnable::Pseudogene
  Caller     : general

=cut


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  #SET UP ANY INSTANCE VARIABLES

  $self->{'_modified_genes'} =[] ; # array ref to modified genes to write to new db
  $self->{'_discarded_transcripts'} = []; # array ref to discarded transcripts
  $self->{'_genes'} = [];	#array of genes to test;  
  $self->{'_repeats'} = {};	# hash of repeat blocks corresponding to each gene;
  $self->{'_real'} = 0;		# scalar number of real genes identified
  $self->{'_pseudogenes'} = 0;	#scalar number of pseudogenes identified
  $self->{'_indeterminate_genes'} = [];	#array of indeterminategenes dbIDs identified
  $self->{'_single_exon_genes'} = [];	#array of single exon  gene dbIDs identified
  $self->{'_multi_exon_genes'} = [];	#array of multiple exon  gene to make into a blast database for spliced elsewhere 
 
  my( $genes,$repeat_features,$PS_REPEAT_TYPES,$PS_FRAMESHIFT_INTRON_LENGTH,$PS_MAX_INTRON_LENGTH,
      $PS_MAX_INTRON_COVERAGE,$PS_MAX_EXON_COVERAGE,
      $PS_NUM_FRAMESHIFT_INTRONS,$PS_NUM_REAL_INTRONS,$SINGLE_EXON,
      $INDETERMINATE,$PS_MIN_EXONS,$PS_MULTI_EXON_DIR,
      $BLESSED_BIOTYPES,$PS_PSEUDO_TYPE,$PS_BIOTYPE,$PS_REPEAT_TYPE,$DEBUG) = rearrange
	([qw(
	     GENES
	     REPEAT_FEATURES
	     PS_REPEAT_TYPES
	     PS_FRAMESHIFT_INTRON_LENGTH 
	     PS_MAX_INTRON_LENGTH
	     PS_MAX_INTRON_COVERAGE
	     PS_MAX_EXON_COVERAGE
	     PS_NUM_FRAMESHIFT_INTRONS
	     PS_NUM_REAL_INTRONS
	     SINGLE_EXON
	     INDETERMINATE
	     PS_MIN_EXONS
	     PS_MULTI_EXON_DIR
	     BLESSED_BIOTYPES
	     PS_PSEUDO_TYPE
	     PS_BIOTYPE
       PS_REPEAT_TYPE
       DEBUG )], @args);

  $self->genes($genes);
  $self->repeats($repeat_features);
  $self->PS_REPEAT_TYPES($PS_REPEAT_TYPES);
  $self->PS_FRAMESHIFT_INTRON_LENGTH($PS_FRAMESHIFT_INTRON_LENGTH);
  $self->PS_MAX_INTRON_LENGTH($PS_MAX_INTRON_LENGTH);
  $self->PS_MAX_INTRON_COVERAGE($PS_MAX_INTRON_COVERAGE);
  $self->PS_MAX_EXON_COVERAGE($PS_MAX_EXON_COVERAGE);
  $self->PS_NUM_FRAMESHIFT_INTRONS($PS_NUM_FRAMESHIFT_INTRONS);
  $self->PS_NUM_REAL_INTRONS($PS_NUM_REAL_INTRONS);
  $self->SINGLE_EXON($SINGLE_EXON);
  $self->INDETERMINATE($INDETERMINATE);
  $self->PS_MIN_EXONS($PS_MIN_EXONS);
  $self->PS_MULTI_EXON_DIR($PS_MULTI_EXON_DIR);
  $self->BLESSED_BIOTYPES($BLESSED_BIOTYPES);
  $self->PS_PSEUDO_TYPE($PS_PSEUDO_TYPE);
  $self->PS_BIOTYPE($PS_BIOTYPE);
	$self->PS_REPEAT_TYPE($PS_REPEAT_TYPE);
  $self->DEBUG($DEBUG);
  return $self;
}

=head2 run

Arg [none] :
  Description: runs the  runnable
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub run {
  my ($self) = @_;
  $self->test_genes;
  $self->summary;
  if ($self->SINGLE_EXON){
    # Write out multiple exon genes for making into blast db for Spliced elsewhere
    my $filename = $self->create_filename('multi_exon_seq','fasta',$self->PS_MULTI_EXON_DIR);
    $self->write_seq_array($self->multi_exon_genes,$filename);
  }
  return 0;
}

=head2 summary

Arg [none] :
  Description: prints out some data about the results
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub summary {
  my ($self) = @_;
  if ( $self->real){print STDERR   $self->real." real genes identified \n";}
  if ( $self->pseudogenes){print STDERR   $self->pseudogenes." pseudogenes identified \n";}
  if ( $self->repeatgenes){print STDERR   $self->repeatgenes." repeat genes identified \n";}
  if ( $self->discarded_transcripts){print STDERR   scalar(@{$self->discarded_transcripts})." pseudotranscripts to be chucked \n";}
  
  foreach my $transcript (@{$self->discarded_transcripts}) {
    print STDERR   $transcript->dbID."\n";
  }
  if ($self->overlooked_genes){
    print STDERR "Overlooked genes\n";
    foreach my $hash_ref (@{$self->overlooked_genes}) {
      my %hash = %{$hash_ref};
      foreach my $gene (keys %hash){
	print STDERR  "$gene ".scalar(@{$hash{$gene}})."\t transcript ".@{$hash{$gene}}[0]->dbID." \n";
      }
    }
  }
  if ($self->SINGLE_EXON){
    print STDERR scalar(@{$self->single_exon_genes})." single exon genes identified and held back for further study\n";
  }
  if ($self->INDETERMINATE){
    print STDERR scalar(@{$self->indeterminate_genes})." indeterminate genes identified and held back for further study\n";
  }
  return 1;
}

=head2 test_genes

Arg [none] :
  Description: Check genes one transcript at at a time pushes test result for each transcript onto an array, tests array to make final decision  genes are classed as pseudogenes if the following criteria are met:
1. At least 80% of the introns are covered with repeats and the total intron length is smaller than 5kb and the gene has a least 1 real and 1 frameshifted intron (default values).
2. All of the introns are short frameshifted introns
3. A two exon gene with intron covered with a protein feature
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub test_genes{
  my $self = shift;
  my @evidence;
  my $num=0;
  my $pseudo= 0;
  my $possible = 0;
  my @genes = @{$self->genes};
 
 GENE:  foreach my $gene (@genes) {
    my %trans_type = ();

  TRANS: foreach my $transcript (@{$gene->get_all_Transcripts}) {

      my $evidence = $self->transcript_evidence($transcript,$gene);

      # store the single exon gene for further analysis unless they are all covered by repeats
      if ($self->SINGLE_EXON){
	if (scalar(@{$transcript->get_all_Exons()})==1) {
	  unless ( $evidence->{'covered_exons'} && $evidence->{'covered_exons'} >= $self->PS_MAX_EXON_COVERAGE ) {
	    push @{$trans_type{'single_exon'}}, $transcript;
	    next TRANS;
	  }
	}
      }

      $num++;

      #transcript tests

      #CALL PSEUDOGENE IF AT LEAST 80% COVERAGE OF INTRONS BY REPEATS
      #AT LEAST 1 F/S EXON AND 1 REAL EXON (?)
      #TOTAL INTRON LENGTH < 5K
	
      if ($evidence->{'total_intron_len'} < $self->PS_MAX_INTRON_LENGTH &&
	  $evidence->{'frameshift_introns'} >= $self->PS_NUM_FRAMESHIFT_INTRONS &&
	  $evidence->{'real_introns'} >= $self->PS_NUM_REAL_INTRONS &&
	  $evidence->{'covered_introns'} >= $self->PS_MAX_INTRON_COVERAGE  ) {
	push @{$trans_type{'pseudo'}}, $transcript;
	print STDERR $gene->dbID." - repeats in introns in transcript ".$transcript->dbID."\n" if $self->DEBUG;
	print STDERR join (', ',%{$evidence}),"\n" if $self->DEBUG;
	next TRANS;
      }

      #ALL FRAMESHIFTED - it is a pseudogene
	
      if ($evidence->{'num_introns'} && 
	  $evidence->{'frameshift_introns'} == $evidence->{'num_introns'}) {
	push@{$trans_type{'pseudo'}},  $transcript;
	next TRANS;
      }

    # repeats covering exons

    if ($evidence->{'covered_exons'} && 
	$evidence->{'covered_exons'} >= $self->PS_MAX_EXON_COVERAGE ) {
      push@{$trans_type{'repeat'}},  $transcript;
      next TRANS;
    }


      if ($self->INDETERMINATE){
	# Tests for genes that look odd to run PSILC over
	# anything with any frameshifts above the cutoff or
	# filled introns but lacking the other criteria
	
	if ($evidence->{'frameshift_introns'} >= $self->PS_NUM_FRAMESHIFT_INTRONS or
	    $evidence->{'covered_introns'} >= $self->PS_MAX_INTRON_COVERAGE  ) {
	 push @{$trans_type{'indeterminate'}}, $transcript;
	  print STDERR $gene->dbID." - looks a little dodgy ".$transcript->dbID."\n" if $self->DEBUG;
	  print STDERR join (', ',%{$evidence}),"\n" if $self->DEBUG;
	  next TRANS;
	}
      }

      # Tests for situation where 2 exon gene has a protein feasture
      # covering intron, ie it has spliced around somthing bad...

      if ($evidence->{'real_introns'} == 1) {
	my $judgement = $self->protein_covered_intron($transcript,$gene);
	if ($judgement eq "dodgy"){
	  push @{$trans_type{'pseudo'}}, $transcript;
	  next TRANS;
	}	
      }

      # transcript passes all tests, it is real
      push @{$trans_type{'real'}}, $transcript;
      push @{$trans_type{'not_multi_exon'}}, $transcript 
	if scalar @{$transcript->get_all_Exons} < $self->PS_MIN_EXONS;
    }

    #########################################
    # gene tests

    #########################################
    # All transcripts have only one exon
    # put them to be stored for spliced elsewhere
    # if specified in config
    # Dont store IG segment genes, just protein codinf

    if ($self->SINGLE_EXON){
      if ($trans_type{'single_exon'}){
	if ( $gene->biotype eq $self->PS_BIOTYPE ){
	  unless ($trans_type{'pseudo'} or
		  $trans_type{'indeterminate'} or
		  $trans_type{'real'} or
		  $trans_type{'repeat'}) {
	    $self->single_exon_genes($gene->dbID);
	    next GENE;
	  }
	}
      }
    }

    #################################################################
    # gene is pseudogene, all transcripts are pseudo
    # set type to pseudogene chuck away all but the longest transcript

      if ($trans_type{'pseudo'}){
	unless ($trans_type{'single_exon'} or
		$trans_type{'indeterminate'} or
		$trans_type{'real'} or
		$trans_type{'repeat'}) {
	  $gene->biotype($self->PS_PSEUDO_TYPE);
	  my @pseudo_trans = @{$trans_type{'pseudo'}};
	  @pseudo_trans = sort {$a->length <=> $b->length} @pseudo_trans;
	  my $only_transcript_to_keep = pop  @pseudo_trans;  
	  $only_transcript_to_keep = $self->transcript_to_keep($only_transcript_to_keep);

	  foreach my $pseudo_transcript (@pseudo_trans) {
	    my $blessed = $self->_remove_transcript_from_gene($gene,$pseudo_transcript);
	    if ( $blessed ) {
	      print STDERR "Blessed transcript " . $pseudo_transcript->display_id . 
		" is a pseudo transcript \n";
	    }
	  } 
          my $new_gene = Bio::EnsEMBL::Gene->new(); 
          $new_gene->analysis($self->analysis) ; 
          $new_gene->biotype($self->PS_PSEUDO_TYPE); 
          $only_transcript_to_keep->biotype($self->PS_PSEUDO_TYPE); 
          $new_gene->add_Transcript($only_transcript_to_keep) ;
          $gene = $new_gene ;
	  $self->modified_genes($gene);
	  $self->pseudogenes(1);
	  next GENE;
	}
      }
    ##################################################
    # gene is indeterminate, if it has no real 
    # transcripts - all either pseudo or indeterminate
    # dont add gene to list to return instead get its
    # dbID for later

    if ($self->INDETERMINATE){
      if ($trans_type{'indeterminate'}){
	unless ($trans_type{'real'} or
		$trans_type{'single_exon'} or
		$trans_type{'repeat'} ){
	  $self->indeterminate_genes($gene->dbID);
	  next GENE;
	}
      }
    }

    ###############################################
    # gene is real but has some dodgy transcripts
    # delete the dodgy transcripts from the gene
    # do you want to throw away indeterminate transcripts? 
    # - not at the moment

    if ($trans_type{'pseudo'}){
      if ($trans_type{'real'} or 
	  $trans_type{'single_exon'} or 
	  $trans_type{'indeterminate'} or 
	  $trans_type{'repeat'}){
	foreach my $trans (@{$trans_type{'pseudo'}}) {
	  my $blessed = $self->_remove_transcript_from_gene($gene,$trans);
	  if ( $blessed ) {
	    print STDERR "Blessed transcript " . $trans->display_id . 
	      " is a pseudo transcript \n";
	  }
	}
	$self->modified_genes($gene);
	$self->multi_exon_genes($gene) unless $trans_type{'not_multi_exon'};
	$self->pseudogenes(1);
	next GENE;
      }
    }

    ####################################
    # gene and transcripts are real

    unless ($trans_type{'pseudo'} or $trans_type{'repeat'}) {
      $self->modified_genes($gene);
      $self->multi_exon_genes($gene) unless $trans_type{'not_multi_exon'};
      $self->real(1);
      next GENE;
    }

    if ($trans_type{'repeat'}) {
      # if a repeat type is specified label the gene as a repeat and store it
      # otherwise we can just leave it - it wont get written to the final db
      if ($self->PS_REPEAT_TYPE) {
	my @pseudo_trans = @{$gene->get_all_Transcripts};
	@pseudo_trans = sort {$a->length <=> $b->length} @pseudo_trans;
	my $only_transcript_to_keep = pop  @pseudo_trans; 

        $only_transcript_to_keep = $self->transcript_to_keep($only_transcript_to_keep);

#	$self->transcript_to_keep($only_transcript_to_keep);
	foreach my $pseudo_transcript (@pseudo_trans) {
	  my $blessed = $self->_remove_transcript_from_gene($gene,$pseudo_transcript);
	  if ( $blessed ) {
	    print STDERR "Blessed transcript " . $pseudo_transcript->display_id . 
	      " is covered with repeats \n";
	  }
	}
	$gene->biotype($self->PS_REPEAT_TYPE);
	$self->modified_genes($gene);
	$self->repeatgenes(1);
      } else {
	print STDERR "Gene " . $gene->display_id . " is covered with repeats - ignoring it\n";
	$self->repeatgenes(1);
      }
      next GENE;
    }

    #########################################
    # should be nothing left after this point
    # might be good to catch anything left and flag it as weird

    $self->overlooked_genes(\%trans_type);

  }
  return 1;
}

=head2 transcript_evidence

Arg [none] : Bio::EnsEMBL::Transcript
  Description: Test individual transcripts return a hash containing results evidence
  Returntype : hash
  Exceptions : none
  Caller     : general

=cut

##############################################################
#test individual transcripts return a hash containing evidence


sub transcript_evidence{

  my ($self,$transcript,$gene) =@_;
  my $repeat_blocks = $self->get_repeats($gene);
  my $results;
  my  @exons =  @{$transcript->get_all_Exons};
  @exons = sort {$a->start <=> $b->start} @exons;
  my $prev_exon = undef;
  my $total_intron_len = 0;
  my $covered_intron_len = 0;
  my $total_exon_len = 0;
  my $covered_exon_len = 0;
  my $n_real_intron = 0;
  my $n_frameshift_intron = 0;
  my $covered_introns = 0 ;
  my $covered_exons = 0 ;
  
  foreach my $exon (@exons) {
    my $seq_feature_exon = Bio::EnsEMBL::Feature->new(
                                                      -START => $exon->start,
                                                      -END => $exon->end,
                                                      -STRAND => $exon->strand
                                                      );
    # Do intron
    if (defined($prev_exon)) {
      my $intron = Bio::EnsEMBL::Feature->new(
                                              -START => $prev_exon->end+1,
                                              -END => $exon->start-1,
                                              -STRAND => $exon->strand
                                              );
      if ($intron->length > $self->PS_FRAMESHIFT_INTRON_LENGTH) {
	$n_real_intron++;
      } else {
	$n_frameshift_intron++;
      }
      $total_intron_len+=$intron->length;
      $covered_intron_len+=$self->_len_covered($intron,$repeat_blocks);
    }
 
    # Do exon
    $total_exon_len+=$exon->length;
    $covered_exon_len+=$self->_len_covered($seq_feature_exon,$repeat_blocks);
    $prev_exon = $exon;
  }

  #calculate percentage coverage
  
  if ($total_intron_len > 0) {
    $covered_introns = (($covered_intron_len/$total_intron_len)*100);
  }
  if ($total_exon_len >  0) {
    $covered_exons =   (($covered_exon_len/$total_exon_len)*100)
  }
  
  $results = {
	      'num_introns' => $#exons,
	      'total_exon_len' =>  $total_exon_len,
	      'covered_exons' => $covered_exons,
	      'total_intron_len' =>  $total_intron_len,
	      'covered_introns' => $covered_introns,
	      'real_introns' => $n_real_intron,
	      'frameshift_introns' => $n_frameshift_intron
	     };
  return $results;
}


=head2 len_covered

 Args       : Bio::Seq::Feature object, reference to an array of repeat blocks
 Description: measures how much of the seq feature (intron or exon) is covered by repeat blocks
 Returntype : scalar

=cut 

sub _len_covered {
  my ($self,$feat,$repeat_blocks_ref) = @_;
#   print STDERR  "FT " . $feat->start . " " . $feat->end . "\n";
  my $covered_len = 0;
 RBLOOP: foreach my $repeat_block (@$repeat_blocks_ref) {
#    print STDERR  "RB " . $repeat_block->start . " " . $repeat_block->end . "\n";
    if ($repeat_block->overlaps($feat, 'ignore')) {
      my $inter = $self->intersection($feat,$repeat_block);
      $covered_len += $inter->length;
    } elsif ($repeat_block->start > $feat->end) {
      last RBLOOP;
    }
  }
  return  $covered_len;
}

=head2 protein_covered_intron

  Args       : Bio::EnsEMBL::Transcript object, Bio::EnsEMBL::Gene object 
  Description: decides if 'real' intron in transcript is covered with a protein feature
  Returntype : scalar

=cut 


sub protein_covered_intron{
  my ($self,$transcript,$gene) =@_;
  my %seq_features;
  my $identified;
  my @all_exons  = @{$transcript->get_all_Exons};
  @all_exons = sort {$a->start <=> $b->start} @all_exons;  

  my @exons;
  # find real intron
 EXON: for (my $i = 1 ; $i <= $#all_exons ; $i++){
    my $intron_length  = $all_exons[$i]->start - $all_exons[$i-1]->end;
    if ($intron_length > 9) {
      # real intron
      push @exons , $all_exons[$i-1];
      push @exons , $all_exons[$i];
      last EXON;
    }
  }
  $self->throw("real intron not found for gene " . $gene->dbID . " exons : @all_exons\n")  unless (scalar(@exons) == 2); 
  my @exon_features = @{$exons[0]->get_all_supporting_features};
  push @exon_features,@{$exons[1]->get_all_supporting_features};

  
  ########################################
  # make a seq feature represening the intron
  
  my $intron = Bio::EnsEMBL::Feature->new(
					  -START => $exons[0]->end+1,
					  -END => $exons[1]->start-1,
					  -STRAND => $exons[0]->strand
					 );
  if (@exon_features) {
    ###########################################################
    # get featues off both exons, split them into ungapped
    # sections and test them against the intron
    # feature see if there is an overlap (80% by default)
    # need to group features by sequence name

  FEATURES:   foreach my $feat (@exon_features) {
      my @sub_features = $feat->ungapped_features;
      foreach my $subfeature (@sub_features) {
	my $seq_feature = Bio::EnsEMBL::Feature->new(
						     -START => $subfeature->start,
						     -END => $subfeature->end,
						     -STRAND => $subfeature->strand
						       );
	push @{$seq_features{$feat->hseqname}}, $seq_feature;
      }
    }
  }
  foreach my $key (keys %seq_features){
    if ($intron && scalar(@{$seq_features{$key}}) > 0) {
      my @features = @{$seq_features{$key}};
      my @feature_blocks;

      #########################################################
      # merge overlapping features together for protein features
      # grouped by name

      my $curblock = undef;
      @features = sort {$a->start <=> $b->start} @features;
      foreach my $feature (@features){
	if (defined($curblock) && $curblock->end >= $feature->start) {
	  if ($feature->end > $curblock->end) { 
	    $curblock->end($feature->end); 
	  }
	} else {
	  $curblock = Bio::EnsEMBL::Feature->new(
						 -START => $feature->start,
						 -END => $feature->end, 
						 -STRAND => $feature->strand
						);
	  push (@feature_blocks,$curblock);
	}
      }
      my $coverage =  $self->_len_covered($intron,\@feature_blocks)."\n";
      if ($coverage/$intron->length*100 > $self->PS_MAX_INTRON_COVERAGE) {
	$identified++;
	print STDERR $transcript->dbID." two exon with $key covering intron ".$coverage/$intron->length*100 . "%.\t";

	# need more than one peice of protein evidence to make the call

	if ($identified >1){
	  print STDERR "\ncalling ".$transcript->dbID." as having a protein covered intron\n";
	  return "dodgy";
	}
	else{
	  print "need another piece of evidence though....\n";
	  }
      }
    }
  }
  return 1;
}

=head2 remove_transcript_from_gene

  Args       : Bio::EnsEMBL::Gene object , Bio::EnsEMBL::Transcript object
  Description: steves method for removing unwanted transcripts from genes
  Returntype : scalar

=cut 


sub _remove_transcript_from_gene {
  my ($self, $gene, $trans_to_del)  = @_;
  # transcript is blessed dont delete it
  return 'BLESSED' if $self->BLESSED_BIOTYPES->{$trans_to_del->biotype};

  $self->discarded_transcripts($trans_to_del);
  $trans_to_del->translation(undef);
  my @newtrans;
  foreach my $trans (@{$gene->get_all_Transcripts}) {
    if ($trans != $trans_to_del) {
      push @newtrans,$trans;
    }
  }
  # The naughty bit!
  $gene->{_transcript_array} = [];
  
  foreach my $trans (@newtrans) {
    $gene->add_Transcript($trans);
  }
  
  return ;
}

=head2 transcript_to_keep

  Args       : Bio::EnsEMBL::Transcript object
  Description: removes the translation provided it is not a blessed transcript
  Returntype : scalar

=cut 


sub transcript_to_keep {
  my ($self, $trans_to_keep)  = @_;
   if  (  $self->BLESSED_BIOTYPES->{$trans_to_keep->biotype} ) {   
    return ;   
  }else { 
   my $tr= Bio::EnsEMBL::Transcript->new( -EXONS => $trans_to_keep->get_all_Exons ) ;  
   return $tr;
  }
}

=head2 intersection

Arg [none] :
  Description: returns length of intersecion of two features
  Returntype : scalar
  Exceptions : throws if not given two Bio::EsEMBL::Feature objects
  Caller     : len_covered

=cut

sub intersection {
  my ($self,$feat1,$feat2) = @_;
    unless ($feat1->isa("Bio::EnsEMBL::Feature") && $feat2->isa("Bio::EnsEMBL::Feature")){
    $self->throw("object is not a Bio::EnsEMBL::Feature they are feat1 $feat1, feat2 $feat2\n$@\n");
  }
  my @start = sort {$a<=>$b}
    ($feat1->start(), $feat2->start());
    my @end   = sort {$a<=>$b}
    ($feat1->end(),   $feat2->end());

    my $start = pop @start;
    my $end = shift @end;

    if($start > $end) {
	return undef;
    } else {
	return Bio::EnsEMBL::Feature->new('-start'  => $start,
			                  '-end'    => $end,
			                  '-strand' => $feat1->strand
                       			  );
    }
}

=head2 output

Arg [none] :
  Description: returns output of running Pseudogene
  Returntype : @{Bio::EnsEMBL::Gene}
  Exceptions : none
  Caller     : general

=cut

sub output {
  my ($self) = @_;
  for my $g ( @{ $self->modified_genes } ) {   
     my @t =@{ $g->get_all_Transcripts} ; 
     for my $t ( @t ) {  
       if ( $g->biotype eq $self->PS_PSEUDO_TYPE ) {  
         if ( defined ( $t->translation ) ) {  
            throw ("pseudogene cant have tranlation. game over \n") ; 
         } 
      } 
    }   
  }   
  return $self->modified_genes;
}


=head2 write_seq_array

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable
  Arg [2]   : Array of Bio::Ensembl::Gene objects
  Arg [3]   : filename
  Function  : This uses Bio::SeqIO to dump a transcript sequence to a fasta file
  Returntype: string, filename
  Exceptions: throw if failed to write sequence
  Example   : 

=cut


sub write_seq_array{
  my ($self, $genes, $filename) = @_;
  return 0 unless (scalar(@{$genes}>0));
  if(!$filename){
    $self->throw("FAILED to write genes - no filename");
  }
  my $seqout = Bio::SeqIO->new(
                               -file => ">".$filename,
                               -format => 'Fasta',
                              );
  eval{
    foreach my $gene (@{$genes}){
      foreach my $transcript (@{$gene->get_all_Transcripts}){
        next unless ( $transcript->translateable_seq );
        $seqout->write_seq($transcript->seq);
      }
    }
  };
  if($@){
    $self->throw("FAILED to write to $filename Runnable:write_seq_file\n\n$@");
  }
  return $filename;
}



#################################################################################
# Container methods

=head2 modified_genes

Arg [1]    : array ref
  Description: get/set modified gene set to return 
  Returntype : array ref to Bio::EnsEMBL::Gene objects
  Exceptions : throws if not a Bio::EnsEMBL::Gene
  Caller     : general

=cut

sub modified_genes {
  my ($self, $modified_genes) = @_;
  if ($modified_genes) {
    unless  ($modified_genes->isa("Bio::EnsEMBL::Gene")){
      $self->throw("Input isn't a Bio::EnsEMBL::Gene, it is a $modified_genes\n$@");
    }
    push @{$self->{'_modified_genes'}}, $modified_genes;
  }
  return $self->{'_modified_genes'};
}

=head2 multi_exon_genes

Arg [1]    : array ref
  Description: get/set multi exon gene set to return 
  Returntype : array ref to Bio::EnsEMBL::Gene objects
  Exceptions : throws if not a Bio::EnsEMBL::Gene
  Caller     : general

=cut

sub multi_exon_genes {
  my ($self, $multi_exon_genes) = @_;
  if ($multi_exon_genes) {
    unless  ($multi_exon_genes->isa("Bio::EnsEMBL::Gene")){
      $self->throw("Input isn't a Bio::EnsEMBL::Gene, it is a $multi_exon_genes\n$@");
    }
    push @{$self->{'_multi_exon_genes'}}, $multi_exon_genes;
  }
  return $self->{'_multi_exon_genes'};
}


=head2 discarded transcripts

Arg [1]    : array ref
  Description: get/set modified gene set to throw away
  Returntype : array ref to Bio::EnsEMBL::Gene objects
  Exceptions : throws if not a Bio::EnsEMBL::Gene
  Caller     : general

=cut

sub discarded_transcripts {
  my ($self, $discarded_transcripts) = @_;
  if ( $discarded_transcripts) {
    unless  ($discarded_transcripts->isa("Bio::EnsEMBL::Transcript")){
      $self->throw("Input isn't a Bio::EnsEMBL::Gene, it is a $discarded_transcripts\n$@");
    }
    push @{$self->{'_discarded_transcripts'}}, $discarded_transcripts;
  }
  return $self->{'_discarded_transcripts'};
}

=head2 genes

Arg [1]    : array ref
  Description: get/set gene set to run over
  Returntype : array ref to Bio::EnsEMBL::Gene objects
  Exceptions : throws if not a Bio::EnsEMBL::Gene
  Caller     : general

=cut

sub genes {
  my ($self, $genes) = @_;
  if ($genes) {
    foreach my $gene (@{$genes}) {
      unless  ($gene->isa("Bio::EnsEMBL::Gene")){
	$self->throw("Input isn't a Bio::EnsEMBL::Gene, it is a $gene\n$@");
      }
    }
    $self->{'_genes'} = $genes;
  }
  return $self->{'_genes'};
}

=head2 repeats

Arg [1]    : array ref
  Description: set repeat set to test genes against
  Returntype : none
  Exceptions : throws if not a Bio::EnsEMBL::SeqFeature
  Caller     : general

=cut

sub repeats {
  my ($self, $repeats) = @_;
  foreach my $repeat_array (values %{$repeats}) {
    foreach my $repeat (@{$repeat_array}) {
      unless ($repeat->isa("Bio::EnsEMBL::Feature")){
        $self->throw("Input is not a Bio::EnsEMBL::Feature, it is a $repeat");
      }
    }
  }
  $self->{'_repeats'} = $repeats;
  return  1;
}

=head2 get_repeats

Arg [1]    : array ref
 Description: get repeat array using a gene object as the key
  Returntype : array ref to Bio::EnsEMBL::SeqFeature objects
  Exceptions : warns if no values corresponding to key
  Caller     : general

=cut

sub get_repeats {
  my ($self, $gene) = @_;
  return  $self->{'_repeats'}->{$gene};
}


=head2 real

Arg [1]    : none
 Description: scalar number of real genes found
  Returntype : scalar integer
  Exceptions : none
  Caller     : general

=cut

sub real {
  my ($self,$num) = @_;
  if ($num){
    $self->{'_real'} += $num;
  }
  return $self->{'_real'};
}


=head2 pseudogenes

Arg [1]    : none
 Description: scalar number of pseudogenes found
  Returntype : scalar integer
  Exceptions : none
  Caller     : general

=cut

sub pseudogenes{
  my ($self,$num) = @_;
  if ($num){
    $self->{'_pseudogenes'}+= $num;
  }
  return $self->{'_pseudogenes'};
}

sub repeatgenes{
  my ($self,$num) = @_;
  if ($num){
    $self->{'_repeatgenes'}+= $num;
  }
  return $self->{'_repeatgenes'};
}


sub indeterminate_genes{
  my ($self,$input) = @_;
  if ($input){
    push @{$self->{'_indeterminate_genes'}}, $input;
  }
  return $self->{'_indeterminate_genes'};
}

sub single_exon_genes{
  my ($self,$input) = @_;
  if ($input){
    push @{$self->{'_single_exon_genes'}}, $input;
  }
  return $self->{'_single_exon_genes'};
}

sub overlooked_genes{
  my ($self,$input) = @_;
  if ($input){
    push @{$self->{'_overlooked_genes'}}, $input;
  }
  return $self->{'_overlooked_genes'};
}

sub PS_FRAMESHIFT_INTRON_LENGTH{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'PS_FRAMESHIFT_INTRON_LENGTH'} = $arg;
  }
  return $self->{'PS_FRAMESHIFT_INTRON_LENGTH'};
}

sub PS_MAX_INTRON_LENGTH{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'PS_MAX_INTRON_LENGTH'} = $arg;
  }
  return $self->{'PS_MAX_INTRON_LENGTH'};
}

sub PS_REPEAT_TYPES{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'PS_REPEAT_TYPES'} = $arg;
  }
  return $self->{'PS_REPEAT_TYPES'};
}

sub PS_MAX_INTRON_COVERAGE{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'PS_MAX_INTRON_COVERAGE'} = $arg;
  }
  return $self->{'PS_MAX_INTRON_COVERAGE'};
}


sub PS_MAX_EXON_COVERAGE{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'PS_MAX_EXON_COVERAGE'} = $arg;
  }
  return $self->{'PS_MAX_EXON_COVERAGE'};
}

sub PS_NUM_FRAMESHIFT_INTRONS {
  my ($self, $arg) = @_;
  if($arg){
    $self->{'PS_NUM_FRAMESHIFT_INTRONS'} = $arg;
  }
  return $self->{'PS_NUM_FRAMESHIFT_INTRONS'};
}

sub PS_NUM_REAL_INTRONS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'PS_NUM_REAL_INTRONS'} = $arg;
  }
  return $self->{'PS_NUM_REAL_INTRONS'};
}


sub SINGLE_EXON{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'SINGLE_EXON'} = $arg;
  }
  return $self->{'SINGLE_EXON'};
}


sub INDETERMINATE{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'INDETERMINATE'} = $arg;
  }
  return $self->{'INDETERMINATE'};
}

sub PS_MIN_EXONS{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'PS_MIN_EXONS'} = $arg;
  }
  return $self->{'PS_MIN_EXONS'};
}

sub PS_MULTI_EXON_DIR{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'PS_MULTI_EXON_DIR'} = $arg;
  }
  return $self->{'PS_MULTI_EXON_DIR'};
}

sub BLESSED_BIOTYPES{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'BLESSED_BIOTYPES'} = $arg;
  }
  return $self->{'BLESSED_BIOTYPES'};
}

sub PS_BIOTYPE{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'PS_BIOTYPE'} = $arg;
  }
  return $self->{'PS_BIOTYPE'};
}

sub PS_PSEUDO_TYPE{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'PS_PSEUDO_TYPE'} = $arg;
  }
  return $self->{'PS_PSEUDO_TYPE'};
}

sub PS_REPEAT_TYPE{
  my ($self, $arg) = @_;
	if($arg){
		$self->{'PS_REPEAT_TYPE'} = $arg;
	}
	  return $self->{'PS_REPEAT_TYPE'};
}					

sub DEBUG{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'DEBUG'} = $arg;
  }
  return $self->{'DEBUG'};
}

1;
