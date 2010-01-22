# Cared for by Ensembl
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::RefineSolexaGenes

=head1 SYNOPSIS

  my $db      = Bio::EnsEMBL::DBAdaptor->new($locator);
  my $refine_genes = Bio::EnsEMBL::Analysis::RunnableDB::RefineSolexaGenes->new (
                                                    -db      => $db,
                                                    -input_id   => $input_id
                                                    -analysis   => $analysis );
  $refine_genes->fetch_input();
  $refine_genes->run();
  $refine_genes->write_output(); #writes to DB



=head1 DESCRIPTION

This module takes intron spanning dna_align_features and combines them with 
rough transcript models to build refined genes with CDS. The module produces 
transcripts representing all possible combinations of introns and exons which
are then filtered according to the specifications in the config.
The databases containing the various features to combine is defined in
Bio::EnsEMBL::Analysis::Config::Databases and the configuration for the 
module is defined in Bio::EnsEMBL::Analysis::Config::GeneBuild::RefineSolexaGenes

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::RefineSolexaGenes;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;

use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils ;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils ;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::RefineSolexaGenes;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw (rearrange);


@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);

my $limit = 0;

sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->read_and_check_config($REFINESOLEXAGENES_CONFIG_BY_LOGIC);
  # Hard limit to the number of possible paths to explore
  $self->recursive_limit(10000);
  # initialise intron feature cash
  my %feature_hash;
  $self->feature_cash(\%feature_hash);
  return $self;
}

=head2 fetch_input

    Title        :   fetch_input
    Usage        :   $self->fetch_input
    Returns      :   nothing
    Args         :   none

=cut

sub fetch_input {
  my( $self) = @_;
  my $logic = $self->analysis->logic_name;
  # fetch adaptors and store them for use later
  $self->intron_slice_adaptor($self->get_dbadaptor($self->INTRON_DB)->get_SliceAdaptor);
  $self->repeat_feature_adaptor($self->db->get_RepeatFeatureAdaptor);
  $self->gene_slice_adaptor($self->get_dbadaptor($self->MODEL_DB)->get_SliceAdaptor);

  # want a slice and a full chromsome to keep everything on in the same coords
  my $slice = $self->fetch_sequence($self->input_id);
  my $gene_slice =  $self->gene_slice_adaptor->fetch_by_region
    ( 
     'toplevel',
     $slice->seq_region_name,
     $slice->start,
     $slice->end,
     1
    );
  my $repeat_slice = $self->intron_slice_adaptor->fetch_by_region
    ('toplevel',
     $slice->seq_region_name,
     $slice->start,
     $slice->end,
     1
    );
  
  my $chr_slice = $self->db->get_SliceAdaptor->fetch_by_region
    (
     'toplevel',
     $slice->seq_region_name
    );
  $self->chr_slice($chr_slice);

  my @repeats = sort { $a->start <=> $b->start } @{$self->repeat_feature_adaptor->fetch_all_by_Slice($repeat_slice)} ;
  # put on chromosome coords
  foreach my $repeat ( @repeats ) {
    $repeat = $repeat->transfer($chr_slice);
  }
  $self->repeats($self->make_repeat_blocks(\@repeats));
  # pre fetch all the intron features
  $self->dna_2_simple_features($slice->start,$slice->end);

  # we want to fetch and store all the rough transcripts at the beginning - speeds things up
  # also we want to take out tiny introns - merge them into longer structures
  my @prelim_genes;
  my @genes;
  if ( $self->MODEL_LN ) {
    @genes = @{$gene_slice->get_all_Genes_by_type( undef,$self->MODEL_LN )};
    print STDERR "Got " .  scalar(@genes) . " genes with logic name " . $self->MODEL_LN ."\n";
  } else {
    @genes = @{$gene_slice->get_all_Genes};
    print STDERR "Got " .  scalar(@genes) . "  genes  \n"; 
  }
  foreach my $gene ( @genes ) {
    # put them on the chromosome
    $gene = $gene->transfer($chr_slice);
    push @prelim_genes,$gene ;
  }
  # deterine strandedness ( including splitting merged genes )
  $self->prelim_genes(\@prelim_genes);
}

sub run {
  my ($self) = @_;
  $self->refine_genes;
}

=head2 refine_genes

    Title        :   refine_genes
    Usage        :   $self->refine_genes
    Returns      :   nothing
    Args         :   none
    Description  :   Combines exons with introns in all possible combinations to 
                 :   Make a series of transcript models

=cut

sub refine_genes {
  my ($self) = @_;
 GENE:  foreach my $gene ( @{$self->prelim_genes} ) {
    my @models;
    my $single_exon = 0;
    # first run on the fwd strand then on the reverse
  STRAND: for ( my $strand = -1 ; $strand <=1 ; $strand+= 2 ) {
      if ( $self->recursive_limit > 10000 ) {
	# mset recursion to 10000 in case it was raised for a tricky gene
	$self->recursive_limit(10000);
	warn("lowering recursive limit after complex gene\n");
      }
      print STDERR "Running on strand $strand \n";
      my %intron_count;
      my @exon_intron;
      my %intron_hash;
      my @exon_prev_intron;
      my %intron_exon;
      my $most_real_introns = 0;
      my $highest_score = 0;
      print STDERR $gene->stable_id. " : " .  $gene->start . " " . $gene->end . ":\n";
      my @exons =  @{$self->merge_exons($gene,$strand)};
      my $exon_count = $#exons;
      my @fake_introns;
      my %known_exons; 
    EXON:   for ( my $i = 0 ; $i <= $exon_count ; $i ++ ) {
	my $exon = clone_Exon($exons[$i]);
	my $retained_intron;
	my $left_introns = 0;
	my $right_introns = 0;
	$exon->{'left_mask'} = 0;
	$exon->{'right_mask'} = $exon->length;
	print STDERR "$i : " . $exon->start . " " . $exon->end . ":\n";
	# make intron features by collapsing the dna_align_features
	my @introns = @{$self->fetch_intron_features($exon->seq_region_start,$exon->seq_region_end)};
	my @left_c_introns;
	my @right_c_introns;
	my @left_nc_introns;
	my @right_nc_introns;
	my @filtered_introns;
	my $intron_overlap;
	my @retained_introns;
      INTRON: foreach my $intron ( @introns ){
	  next unless $intron->strand == $strand;
	  next unless $intron->length > $self->MIN_INTRON_SIZE;
	  # disguard introns that splice over our exon
	  if ( $intron->start< $exon->start && $intron->end > $exon->end ) {
	    $intron_overlap++;
	    next;
	  }
	  # check to see if this exon contains a retained intron
	  if (  $intron->start > $exon->start && $intron->end < $exon->end ) {
	    $retained_intron = 1;
	    $exon->{'retained'} =1;
	    push @retained_introns, $intron;
	  } else {
	    # we need to know how many consensus introns we have to the 
	    # left and  right in order to determine whether to put in 
	    # a non consensus intron
	    if ( $intron->end <= $exon->end ) {
	      if ( $intron->display_label =~ /non-consensus/ ) {
		push @left_nc_introns, $intron if $intron->score > 1;
	      } else {
		push @left_c_introns, $intron;
	      }
	    }
	    if ( $intron->start >= $exon->start ) {
	      if ( $intron->display_label =~ /non-consensus/ ) {
		push @right_nc_introns, $intron if $intron->score > 1;
	      } else {
		push @right_c_introns, $intron;
	      }
	    }
	  }
	}	
	
	# Restrict internal exons splice sites to most common
	# that way our alt splices will all share the same boundaries
	# but have different combinations of exons
	if ( $self->STRICT_INTERNAL_SPLICE_SITES && 
	     # either we apply it equeally to all exons
	     ( $self->STRICT_INTERNAL_END_EXON_SPLICE_SITES or
	       # only apply to internal exons, leave out end exons
	       ( !$self->STRICT_INTERNAL_END_EXON_SPLICE_SITES && 
		 ( scalar(@left_c_introns)  + scalar(@left_nc_introns) ) > 0 &&
		 ( scalar(@right_c_introns) + scalar(@right_nc_introns)) > 0 ))){
	  # pick best left splice
	  my $best_left_splice;
	  my $best_left_score = 0;
	  my @all_left_introns =  @left_c_introns;
	  push @all_left_introns, @left_nc_introns;
	  foreach my $i ( @all_left_introns ) {
	    if ( $best_left_score < $i->score ) {
	      $best_left_score = $i->score;
	      $best_left_splice = $i->end;
	    }
	  }
	  
	  # pick best right  splice
	  my $best_right_splice;
	  my $best_right_score = 0;
	  my @all_right_introns =  @right_c_introns;
	  push @all_right_introns, @right_nc_introns;
	  foreach my $i ( @all_right_introns ) {
	    if ( $best_right_score < $i->score ) {
	      $best_right_score = $i->score;
	      $best_right_splice = $i->start;
	    }
	  }
	  # filter out introns that pick other splice sites
	  foreach my $i ( @all_left_introns ) {
	    push @filtered_introns, $i if $i->end == $best_left_splice;
	  }

	  foreach my $i ( @all_right_introns ) {
	    push @filtered_introns, $i if $i->start == $best_right_splice;
	  }
	  
	} else {
	  
	  # add non consensus introns only where there are no consensus introns
	  push @filtered_introns, @left_c_introns;
	  push @filtered_introns, @right_c_introns;
	  push @filtered_introns, @left_nc_introns  if scalar(@left_c_introns)  == 0;
	  push @filtered_introns, @right_nc_introns if scalar(@right_c_introns) == 0; ;
	}

	if ( scalar(@left_c_introns)  == 0 && scalar(@left_nc_introns)  > 0) {
	  print STDERR "using " . scalar(@left_nc_introns) . " NC left \n";
	} 
	if ( scalar(@right_c_introns)  == 0 && scalar(@right_nc_introns)  > 0 ) {
	  print STDERR "using " . scalar(@right_nc_introns) . " NC right \n";
	}
	
	# single exon models are a special case
	if ( scalar(@exons) == 1 &&  scalar(@filtered_introns)  == 0 &&  scalar(@retained_introns == 0 )) {
	  # at least on this strand this model looks like a single exon
	  $single_exon += 1;
	}
	
	# we dont want to allow left and right introns to overlap - 
	# it leads to -ve length exons
	
	# we put all the retained introns in at the end we want to do all the 
	# entrances and exits to each exon before we look at whether its 
	# retained or not
	@retained_introns = sort { $b->start <=> $a->start } @retained_introns;
	# push @filtered_introns, @retained_introns;
      INTRON:  foreach my $intron ( @filtered_introns ) {
	  print STDERR "\t" . $intron->start . " " . $intron->end . " " . $intron->strand . " " . $intron->display_label . " " . $intron->score . "\n";
	  # becasue we make a new exons where we have a reatained intron to 
	  # stop circular references we need to allow the final 
	  # intron splicing out of the exon to be used more than once
	  # by each new exon in fact
	  $intron_count{$intron->display_label}++ unless $retained_intron;
	  $intron_hash{$intron->display_label} = $intron;
	  # only use each intron twice once at the end and once at the start of
	  # an exon
	  # exon_intron links exons to the intron on their right ignoring strand
	  push @{ $exon_intron[$i]}  , $intron if $intron->end > $exon->end;
	  # intron exon links introns to exons on their right ignoring strand
	  if ( $intron->start < $exon->start ) {
	    push @{$intron_exon{$intron->display_label}} , $i;
	    # exon_prev_intron links exons to introns on their left ignoring strand
	    push @{ $exon_prev_intron[$i]}  , $intron ;
	  }
	}
	if ( scalar( @retained_introns ) ) {
	  print STDERR "Dealing with " . scalar( @retained_introns ) . " retained introns \n";
	  my @new_exons;
	  push @new_exons,  $exon  ;
	  # sort first by start then by end where start is the same
	  @retained_introns =  sort {$a->start <=> $b->start } @retained_introns;
	  for ( my $i = 0; $i < $#retained_introns ; $i++ ) {
	    if ( $retained_introns[$i]->start ==  $retained_introns[$i+1]->start &&
		 $retained_introns[$i]->end >  $retained_introns[$i+1]->end ) {
	      # reverse the order
	      my $temp =  $retained_introns[$i];
	      $retained_introns[$i] = $retained_introns[$i+1];
	      $retained_introns[$i+1] = $temp;
	    }
	  }
	  # now lets deal with any retained introns we have
	RETAINED: foreach my $intron ( @retained_introns ) {
	    # we dont need to make all new exons for each alternate splice
	    # check the intron is still retained given the new exons
	    my $retained = 1;
	    foreach my $new_exon ( @new_exons ) {
	      if  (  $intron->start > $new_exon->start && $intron->end < $new_exon->end ) {
	      } else {
		$retained = 0;
	      }
	      next RETAINED unless $retained;
	    }
	    my $reject_score = 0;
	    # intron is within the exon - this is not a true exon but a retained intron
	    if (  $intron->start > $exon->start && $intron->end < $exon->end && $intron->length > $self->MIN_INTRON_SIZE ) {
	      # we are going to make a new exon and chop it up
	      # add intron penalty
	      print STDERR "RETAINED INTRON PENALTY for " . $intron->display_id ." before " . $intron->score . " ";
	      $reject_score = $intron->score - $self->RETAINED_INTRON_PENALTY;
	      # intron penalty is doubled for nc introns 
	      if ( $intron->display_label  =~ /non-consensus/ ) {
		$reject_score = $reject_score - $self->RETAINED_INTRON_PENALTY;
	      }
	      print STDERR " after " . $reject_score ."\n";
	      if ( $reject_score < 1 ) {
		# treat as single exon
		if ( scalar(@exons) == 1 ) {
		  # at least on this strand this model looks like a single exon
		  $single_exon += 1;
		}
		next;
	      }
	      print STDERR  "Exon " . $exon->start ."\t". $exon->end . " has retained intron:\n     " . $intron->start ."\t" .  $intron->end ." "."\n";
	      # dont have circular references to exons or the paths
	      # will be infinite so clone this exon instead
	      # I guess we also want to keep the original exon too?
	      my $new_exon1 = clone_Exon( $exon );
	      my $new_exon2 = clone_Exon( $exon );
	      # chop it up a bit so it no longer overlaps the other introns
	      print STDERR  "TRIMMING EXON \n";
	      my $length = $intron->end - $intron->start;
	      $new_exon1->end( $intron->start + int( $length / 2 ) - 2 );
	      $new_exon2->start( $intron->end - int( $length / 2 ) + 2 );
	      
	      push @new_exons,$new_exon1 unless $known_exons{$new_exon1->start."-".$new_exon1->end."-".$new_exon1->strand} ;
	      push @new_exons,$new_exon2 unless $known_exons{$new_exon2->start."-".$new_exon2->end."-".$new_exon2->strand};
	      
	      $known_exons{$new_exon1->start."-".$new_exon1->end."-".$new_exon1->strand} = 1;
	      $known_exons{$new_exon2->start."-".$new_exon2->end."-".$new_exon2->strand} = 1;
	    }
	  }
	  if ( scalar (@new_exons > 1 ) ) {
	    # we want to split the score equally across the new exons
	    foreach my $e ( @new_exons ) {
	      foreach my $d ( @{$e->get_all_supporting_features} ) {
		$d->score($d->score / scalar(@new_exons));
	      }
	    }
	    
	    splice( @exons,$i,1,@new_exons);
	    for ( my $i = 0 ; $i<= $#exons ; $i++ ) {
	      my $e = $exons[$i];
	    }
	    print "ADDED " . scalar( @new_exons) . " new exons\n";
	    $exon_count+= $#new_exons;
	  }
	}
      }
      
      next unless @exon_intron;
      # Loop around the path generation, 
      # if there are too many paths to process return undef
      # then re-run the path processing but with increasing strictness
      # where strictness = elimianating alternate low scoring introns
      my $paths;
      my $strict = 0;
      while ( !$paths ) {
	$paths = $self->process_paths( \@exons, \@exon_intron, \%intron_exon, $strict );
	next GENE if $paths && $paths eq 'Give up';
	$strict++;
      }  
      print STDERR "STRAND $strand BEFORE COLLAPSING  PATHS  = " . scalar( keys %$paths ) . "\n";
      # lets collapse redundant paths
      foreach my $path ( sort keys %$paths ) {
	#  print "PATHS $path\n";
	my @array = split ( /\./,$path);
	my ($start,$end,$middle);
	for ( my $j = 0 ; $j < scalar( @array ) ; $j++ )  {
	  $start .= $array[$j] ."."  unless $j < 2;
	  $middle .= $array[$j] ."." unless $j < 2  or $j >= $#array-1 ;
	  $end .= $array[$j] . "." unless $j >= $#array-1;
	}
	# remove redunancy from the array
	delete $paths->{$start} if $start && $paths->{$start};
	delete $paths->{$end} if $end && $paths->{$end};
	delete $paths->{$middle} if $middle && $paths->{$middle};
      }
      
      print STDERR "AFTER COLLAPSING  PATHS  = " . scalar( keys %$paths ) . "\n";
      push @models, @{$self->make_models($paths,$strand, \@exons,$gene,\%intron_hash )};
      print STDERR "Now have " . scalar ( @models ) ." models \n";
    }
    $self->filter_models(\@models);
    

    # process single exon models
    # if it has no introns on either strand
    if ( $self->SINGLE_EXON_MODEL && $single_exon == 2 ) {
      my $exon = $self->merge_exons($gene,1)->[0];
      my $single_exon_model;
      #print STDERR " Single exon = $single_exon\n";
      next unless $exon->length+40 >= $self->MIN_SINGLE_EXON;
      #print STDERR "Passed length filter " . $exon->length ."\n";
      # trim padding 
      $exon->start($exon->start + 20);
      $exon->end  ($exon->end   - 20);
      # get the cds
      my $fwd_exon =  clone_Exon($exon);
      $fwd_exon->strand(1);
      my $rev_exon = clone_Exon($exon);
      my $fwd_t =  new Bio::EnsEMBL::Transcript(-EXONS => [$fwd_exon]);
      my $fwd_tran = compute_translation(clone_Transcript($fwd_t));
      my $fwd_t_len = $fwd_tran->translation->genomic_end - $fwd_tran->translation->genomic_start;
      #print STDERR "FWD t length $fwd_t_len\n";
      my $rev_t =  new Bio::EnsEMBL::Transcript(-EXONS => [$rev_exon]);
      my $rev_tran = compute_translation(clone_Transcript($rev_t));
      my $rev_t_len = $rev_tran->translation->genomic_end - $rev_tran->translation->genomic_start;
      #print STDERR "REV t length $rev_t_len\n";
      if ( $fwd_tran->translateable_seq &&  
	   ( $fwd_t_len / $fwd_tran->length )* 100 >= $self->SINGLE_EXON_CDS &&
	   $fwd_t_len >  $rev_t_len ) {
	# keep this one
	$single_exon_model =  $fwd_tran;
      }
      if ( $rev_tran->translateable_seq &&  
	   ( $rev_t_len / $rev_tran->length )* 100 >= $self->SINGLE_EXON_CDS &&
	   $rev_t_len >  $fwd_t_len ) {
	# keep this one
	$single_exon_model = $rev_tran;
      }
      if ( $single_exon_model ) {
	$single_exon_model->analysis($self->analysis);
	$single_exon_model->version(1);
	my ( $new_gene ) = @{convert_to_genes(($single_exon_model),$gene->analysis)};
	$new_gene->biotype($self->SINGLE_EXON_MODEL); 
	# score comes from exon supporting feature;
	my $score =  $exon->get_all_supporting_features->[0]->score;
	$new_gene->stable_id($gene->stable_id . "-v1-" . int($score) );
	push @{$self->output} , $new_gene;
      }
    }
  }
}

=head2 filter_models

    Title        :   filter_models
    Usage        :   $self->filter_models(\@models);
    Returns      :   nothing
    Args         :   Array ref of array references of Bio::EnsEMBL::Transcript
    Description  :   Labels or removes models overlapping better scoring models on the 
                     opposite strand

=cut

sub filter_models {
  my ( $self, $models )= @_;
  my @clusters = @{$models};
  my @fwd;
  my @rev;
  my @models;
  foreach my $cluster ( @clusters ) {
    next unless $cluster->{'final_models'};
    push @fwd, $cluster if $cluster->{'strand'} == 1;
    push @rev, $cluster if $cluster->{'strand'} == -1;
    push @models, $cluster;
  }
  # overlaps
  foreach my $fc ( @fwd ) {
    foreach my $rc ( @rev ) {
      # one is within the other  or they are the same
      # they proably need to be rejected on the basis of coding overlap
      if ( ( $fc->{'start'} >= $rc->{'start'} && 
	     $fc->{'end'} <= $rc->{'end'} ) or 
	   ( $rc->{'start'} >= $fc->{'start'} && 
	     $rc->{'end'} <= $fc->{'end'}  ) )  {
	
	# do they have coding overlap?
	my @fg = @{$fc->{'final_models'}};
	my @rg = @{$rc->{'final_models'}};
	
	# do they have coding overlap?
      FG: foreach my $fg ( @fg ) {
	  my $ft = $fg->get_all_Transcripts->[0];
	  next unless $ft->translateable_seq;
	  if (  length($ft->translate->seq) <=  100 ) {
	    $fg->biotype('bad');
	    next FG;
	  }
	  foreach my $fe ( @{$ft->get_all_translateable_Exons} ) {
	  RG: foreach my $rg ( @rg ) {
	      my $rt = $rg->get_all_Transcripts->[0];
	      next unless $rt->translateable_seq;
	      if (  length($rt->translate->seq) <=  100 ) {
		$rg->biotype('bad');
		next RG;
	      }
	      foreach my $re ( @{$rt->get_all_translateable_Exons} ) {
		if ( $fe->{'start'} <= $re->{'end'} && 
		     $fe->{'end'}  >=  $re->{'start'}) {
		  # coding overlap	
		  if ( $ft->{'_score'} <  $rt->{'_score'} ) {
		    # get rid of / label the reverse genes 
		    $fg->biotype('bad');
		  } else {
		    $rg->biotype('bad');
		  }
		  next FG;
		}
	      }
	    }
	  }
	}
	
#	# or else they overlap and the antisense model is 2 exon with a cds < 100 AA
#	if (  $fc->{'start'} <= $rc->{'end'} && 
#	      $fc->{'end'} >=  $rc->{'start'} )  {
#	  if ( scalar($fg[0]->get_all_Exons) == 2 && 
#	       scalar($rg[0]->get_all_Exons) > 2 && 
#	       length($fg[0]->get_all_Transcripts->[0]->translate->seq) <=  100 ) {
#	    foreach my $gene ( @fg ) {
#	      $gene->biotype('bad');
#	    }
#	  }
#	  if ( scalar($rg[0]->get_all_Exons) == 2 && 
#	       scalar($fg[0]->get_all_Exons) > 2 && 
#	       length($rg[0]->get_all_Transcripts->[0]->translate->seq) <=  100 ) {
#	    foreach my $gene ( @rg ) {
#	      $gene->biotype('bad');
#	    }
#	  }
#	}
      }
    }
  }

  foreach my $cluster ( @models ) {
    my %exon_use_hash;
    my %exon_starts;
    my %exon_ends;
    my %exon_pattern;
    my $count = 0;
    my $translation_start = 1000000000000;
    my $translation_end = 0;
    foreach my $gene ( @{$cluster->{'final_models'}} ) {
      my $transcript =  $gene->get_all_Transcripts->[0];
      if ( $transcript->translateable_seq ) {
	if (  $transcript->coding_region_start < $translation_start ) {
	  $translation_start =  $transcript->coding_region_start;
	}
	if (  $transcript->coding_region_end > $translation_end ) {
	  $translation_end =  $transcript->coding_region_end;
	}
      }
      my $exon_use = $transcript->{'_exon_use'};
      # if ( $gene->biotype eq $self->OTHER_ISOFORMS ) {
      if ( $exon_use_hash{$exon_use}  ) {
	$gene->biotype('bad');
      }
      #  }
      $exon_use_hash{$exon_use} = 1 ;
      #  print "TRANSCRIPT " . $gene->get_all_Transcripts->[0]->{'_depth'} .
      #" Exon use $exon_use Biotype " . $gene->biotype ."\n";
      my $es = 0;
      my $ee = 0;
      my $pattern;
      foreach my $exon ( @{$gene->get_all_Exons} ) {
	$pattern .= $exon->start.":".$exon->end.":";
	$es++ if $exon_starts{$exon->start};
	$ee++ if $exon_ends{$exon->end};
	$exon_starts{$exon->start} = 1;
	$exon_ends{$exon->end} = 1;
      }
      if ( $ee == scalar(  @{$gene->get_all_Exons} ) &&
	   $es == scalar(  @{$gene->get_all_Exons} ) ) {
	# seen it before - or something very much like it
	$gene->biotype('bad') ;
	#	print "CALLING it bad\n";
      }
      if ( $exon_pattern{$pattern} ) {
	# seen it before - or something very much like it
	$gene->biotype('duplicate') ;
	#	print "CALLING it bad\n";
      }
      $exon_pattern{$pattern} = 1;
    }
    # promote "bad" models that have a cds as long as the best cds to 
    # alt isoforms
    my @final_models = @{$cluster->{'final_models'}};
    my $best_cds = 0;
    for (  my $g = 0; $g < scalar(@final_models) ; $g++ ) {
      my $gene = $final_models[$g];
        my $transcript =  $gene->get_all_Transcripts->[0];
      print "$g - " . $transcript->{'_score'} ." tran length " .
	( $transcript->cdna_coding_end - $transcript->cdna_coding_start ) ."\n";
      if ( $g == 0 ) {
	# best scoring model 
	if  ( $transcript->translateable_seq ) {
	  $best_cds =  $transcript->cdna_coding_end - $transcript->cdna_coding_start;
	}
      }
      my $transcript =  $gene->get_all_Transcripts->[0];
      if ( $transcript->translateable_seq ) {
	if ( $gene->biotype eq 'bad' && 
	     $transcript->coding_region_start == $translation_start &&
	     $transcript->coding_region_end == $translation_end ) {
	  $gene->biotype( $self->OTHER_ISOFORMS );
	}
	if ($gene->biotype eq 'bad' && 
	    $transcript->cdna_coding_end - $transcript->cdna_coding_start  > $best_cds ) {
	  $gene->biotype( $self->OTHER_ISOFORMS );
	}
      } 
      if ( $gene->biotype eq 'bad' ) {
	# change type to  a bad model if it is bad 
	# store it on output array if the bad type is defined
	if ( $self->BAD_MODELS ) {
	  $gene->biotype( $self->BAD_MODELS ) ;
	  push @{$self->output} , $gene if $count <= $self->OTHER_NUM ;
	}
      } else {
	unless ( $gene->biotype eq 'duplicate' ) {
	  push @{$self->output} , $gene if $count <= $self->OTHER_NUM ;
	}
      }
      $count++ if $gene->biotype eq $self->OTHER_ISOFORMS ;
      $count++ if $gene->biotype eq $self->BEST_SCORE ;
    }
  }
}


=head2 make_models

    Title        :   
    Usage        :   $self->make_models($paths, $strand ,$exons,$gene, $intron_hash);
    Returns      :   Array ref of array references of Bio::EnsEMBL::Transcript
    Description  :   Turns abstract paths into Bio::EnsEMBL::Gene models. Paths are
                     clustered and sorted by score - only the top X models for 
                     each cluster of paths get built ( X is defined in config )

=cut

sub make_models {
  my ( $self, $paths, $strand ,$exons,$gene, $intron_hash) = @_;
 
  # paths are stored as text - turn them into arrays of features "models"
  my @clusters;
  my @models;
  my @genes;

  foreach my $path ( keys %$paths ) {
    my $exon_use;
    my @model;
    my $exon_score = 0;
    my $intron_score = 0;
    foreach my $feature ( split(/\./,$path ) ) {
      if ( $feature =~ /^INTRON/ ) {
	push @model, $intron_hash->{$feature};
	$intron_score+= $intron_hash->{$feature}->score;
      } else {
	$exon_use.= "$feature,";
	push @model, $exons->[$feature];
	foreach my $daf ( @{$exons->[$feature]->get_all_supporting_features} ) {
	  $exon_score += $daf->score;
	}
      }
    }
    my $total_score = int($exon_score)/100 + $intron_score;
    # last elements are the strand and score
    push @model, $exon_use;
    push @model, $total_score;
    push @models,\@model;
  }
  # now lets cluster the models so that they are non overlapping
  # and return the clusters arranged by score
  my @model_clusters = @{$self->model_cluster(\@models,$strand)};
  
  # Now we cycle through all the models and turn them into genes
  # we start with the highest scoring modes and work backwards
  # until we have enough 
  my $cluster_count = 0;
  foreach my $cluster ( @model_clusters ) {
    $cluster_count++;
    my @trans;
    my $version = 0 ;
    my $strand = $cluster->{'strand'};
    # we want the array in the reverse order, highest scoring first
    my @models_by_score = @{$cluster->{'models'}};
    # all the models with a particualar score highest first
  SCORES:   for ( my $s = $#models_by_score ; $s >=0; $s-- ) {
      next unless $models_by_score[$s];
      my $models = $models_by_score[$s];
    MODEL: foreach my $model (  @$models ) {
	# list of the rough exons used in the model
	my $exon_use = pop(@{$model});
	my $intron_count = 0;
	my $intron_score = 0;
	my $exon_score = 0;
	my $non_con_introns = 0;
	my @new_exons;
	# make an array containing cloned exons
	for ( my $i = 0 ; $i < scalar(@$model) ; $i++ ) {
	  unless ( $model->[$i]->isa("Bio::EnsEMBL::SimpleFeature") ) {	
	    my $new_exon = clone_Exon($model->[$i]);
	    # add in exon coverage scores from supporting features
	    foreach my $daf ( @{$model->[$i]->get_all_supporting_features} ) {
	      $exon_score += $daf->score;
	    }
	    $new_exon->strand($strand);
	    push @new_exons,$new_exon;
	  } else {
	    push @new_exons, $model->[$i];
	  }
	}
	# trim the exons using the intron features to get the splice sites correct
	for ( my $i = 0 ; $i < scalar(@$model) ; $i++ ) {
	  if ( $model->[$i]->isa("Bio::EnsEMBL::SimpleFeature") ) {
	    my $intron = $model->[$i];
	    next unless $intron->strand == $strand;
	    next unless $new_exons[$i-1] && $new_exons[$i+1];
	    # its an intron trim the exons accordingly
	    $new_exons[$i-1]->end( $intron->start );
	    $new_exons[$i+1]->start( $intron->end );
	    if ( $new_exons[$i-1]->start >=  $new_exons[$i-1] ->end ) {
	   #   warn("Problem with intron $i\nINTRON:" . $intron->start ."-" . $intron->end ."\n");
	   #   warn("EXONS " . $new_exons[$i-1]->start ."-" . $new_exons[$i-1] ->end ." " .( $new_exons[$i-1] ->end -$new_exons[$i-1]->start ) ."bp\n" );
	   #   warn("EXONS " . $new_exons[$i+1]->start ."-" .  $new_exons[$i+1] ->end ."\n");
	   #   warn("This model is not going to work, moving on to the next one\n");
	      next MODEL;
	    }
	    $intron_count++;
	    $intron_score+= $intron->score;
	    #  print "Adding INTRON " . $intron->display_id ." total score now $intron_score \n";
	    # keep tabs on how many fake introns we have
	    $non_con_introns++ if $intron->display_label =~ /non-consensus/;
	  }
	}
	next MODEL unless $intron_count;
	
	# trim padding from the start and end exons
	$new_exons[0]->start($new_exons[0]->start + 20) if $new_exons[0]->length  > 20;
	$new_exons[-1]->end ($new_exons[-1]->end  - 20) if $new_exons[-1]->length > 20;
	
	# make it into a gene
	my @modified_exons;
	foreach my $exon ( @new_exons ) {
	  next if $exon->isa("Bio::EnsEMBL::SimpleFeature");
	  push @modified_exons, clone_Exon($exon);
	}
	if ( $strand == 1 ) {
	  @modified_exons = sort { $a->start <=> $b->start } @modified_exons;
	} else {
	  @modified_exons = sort { $b->start <=> $a->start } @modified_exons;
	}
	# make it into a gene
	my $t =  new Bio::EnsEMBL::Transcript(-EXONS => \@modified_exons);
	# add a translation 
	my $tran = compute_translation(clone_Transcript($t));
	# keep track of the scores for this transcript
	$tran->analysis($self->analysis);
	$tran->version(1);
	# favor longer cds by adding doubling the score for coding exons
	# only use exons that are completely coding otherwise you also 
	# end up adding in score which is really UTR for long terminal exons
	# that have a bit of coding in them
	my $coding_bonus = 0;
	my $coding_exons =0;
	if ( $tran->translateable_seq ) {
	  foreach my $ce ( @{$tran->get_all_Exons} ) {
	    unless ( $ce->phase == -1 or $ce->end_phase == -1 ) {
	      $coding_bonus += $ce->get_all_supporting_features->[0]->score;
	      $coding_exons++;
	    }
	  }
	}
#	print "Coding Bonus of $coding_bonus from $coding_exons completely coding exons \n";
	$tran->{'_score'} =  ( (int ( $intron_score + $exon_score ) / 10 ) + $coding_bonus  );
#	print "Final score = $intron_score + int( $exon_score / 100 ) + $coding_bonus = " . $tran->{'_score'} ;
	$tran->{'_depth'} =  ( $intron_score + $exon_score );
#	print " for tran " .$tran->{'_depth'} . "\n";
	$tran->{'_NC_introns'} =  $non_con_introns ;
	$tran->{'_exon_use'} = $exon_use;
	#print STDERR " EXON count $exon_count\n";
	$tran->{'_intron_count'} = $intron_count;
	push @trans, $tran;
      }
      # we only want one model if its best score only
      if ( $self->BEST_SCORE and not   $self->OTHER_ISOFORMS  ) {
	last SCORES if scalar(@trans)  > 0 ;
      }
      # we want X number of models
      if ( $self->BEST_SCORE &&  $self->OTHER_ISOFORMS && $self->MAX_NUM  ) {
	last SCORES if scalar(@trans)  >= ( $self->MAX_NUM +1 )  ;
      }
    }
    # re-sort the transcripts to take account of the revised scores
    @trans = sort { $b->{'_score'} <=> $a->{'_score'} } @trans;
    my $best;
    foreach my $tran ( @trans ) {
      my ( $new_gene ) = @{convert_to_genes(($tran),$gene->analysis)};
      $version++;
      $new_gene->biotype($self->OTHER_ISOFORMS);
      if ( $version == 1 ) {
	$new_gene->biotype($self->BEST_SCORE);
	$best = $tran;
      } else {
	if ( $tran->translateable_seq && $best->translateable_seq && length($tran->translate->seq) > length($best->translate->seq) ) {
	  # not just longer peptide but at least as long span
	  next unless ( $tran->translateable_seq && $best->translateable_seq &&
			$tran->translation->genomic_end >= $best->translation->genomic_end &&
			$tran->translation->genomic_start <= $best->translation->genomic_start ) ;
	  $new_gene->biotype($self->BEST_SCORE);
	}
      }
      $new_gene->stable_id($gene->stable_id . "-v$cluster_count.$version-" .
			   int($tran->{'_score'}) ."-" .
			   int($tran->{'_depth'}) ."-" .  
			   $tran->{'_intron_count'} ."-NC-" . 
			   $tran->{'_NC_introns'} . "-" . $tran->strand );
      push @{$cluster->{'final_models'}} , $new_gene;
    }
  }
  return \@model_clusters;
}


sub write_output{
  my ($self) = @_;

  my $outdb = $self->get_dbadaptor($self->OUTPUT_DB);
  my $gene_adaptor = $outdb->get_GeneAdaptor; 
  
  my @output = @{$self->output};
  
  my $fails = 0;
  my $total = 0;
  foreach my $gene (@output){
    $gene->analysis($self->analysis);
    $gene->source($self->analysis->logic_name);
    eval {
      $gene_adaptor->store($gene);
    };    
    if ($@){
      warning("Unable to store gene!!\n$@");
      $fails++;
    }
    $total++;
  }
  if ($fails > 0) {
    throw("Not all genes could be written successfully " .
          "($fails fails out of $total)");
  }
}

=head2 ProcessTree

    Title        :   ProcessTree
    Usage        :   $self->ProcessTree
    Returns      :   String containing paths through the gene
    Args         :   A hash reference contianing the possible intron exons
                 :   Integer key for the hashref
                 :   String containing keys used up to this point
                 :   String containing the paths under construction
    Description  :   Recursive method that creates paths that explore all possible
                 :   routes through a hashref, uses a configurable recursion limit 
                 :   to prevent it running out of memory if the paths cannot be solved
                 :   or are too large to be practical

=cut

sub ProcessTree {
  my ($self,$hashref,$index,$sofar,$paths) = @_;
  # dont let it go on for ever eating memory
  if ($limit > $self->recursive_limit){
    print STDERR  "Too many recursive possibilities\n";
    return "ERROR";
  }
  my @node =  keys %{$hashref->{$index}} ;
  $sofar.= "$index.";
  foreach my $child (@node){
    $limit++;
     my $result =  $self->ProcessTree($hashref,$child,$sofar,$paths);
    if ( $result eq "ERROR" ) {
      $limit = 0;
      return "ERROR"; 
    }
   # $result->{$sofar} = 1;
  }
  if ( scalar(@node == 0) ) {
   # print "$sofar\n";
    $paths->{$sofar} = 1;
  }
  return $paths;
}

=head2 ProcessTree

    Title        :   ProcessTree
    Usage        :   $self->process_paths
    Returns      :   String containing paths through the gene
    Args         :   A hash reference contianing the possible intron exons
                 :   Integer key for the hashref
                 :   String containing keys used up to this point
                 :   Integer flag indicating filtering should take place
    Description  :   Filters paths to remove the lowest scoring intron
                 :   for a given pair of exons where more than one intron
                 :   is possible. Filters progressivley if the paths cannot be
                 :   made for the model until the paths can be created or the 
                 :   model cannot be filtered any more, in this case the number
                 :   of recursions can be raised and the process repeated
                 :   untill the max_recursions limit is reached

=cut

sub process_paths{
  my ( $self, $exons, $exon_intron, $intron_exon, $strict ) = @_;
  my $variants;
  my $removed;
  # now lets make a hash of hashes holding which exons connect to which
  for ( my $i = 0 ; $i < scalar(@{$exons}) ; $i ++ ) {
    if ( $exon_intron->[$i] ) {
      if ( $strict ) {
	# group all the introns by exon pairs
	# they all join to left exon 
	# so lets group them by right exon
	my %intron_groups;
	foreach my $intron ( @{$exon_intron->[$i]} ) {
	  foreach my $right_exon ( @{$intron_exon->{$intron->display_label}} ) { 
	    push @{$intron_groups{$right_exon}}, $intron;
	  }	
	}
	# now lets sort these groups by score
	foreach my $group ( keys %intron_groups ) {
	  @{$intron_groups{$group}}  = sort {$b->score <=> $a->score} @{$intron_groups{$group}};
	}
	# now lets see what they look like
	print "EXON $i:". $exons->[$i]->start ." - ". 
	  $exons->[$i]->end ." - ". 
	    $exons->[$i]->strand ."\n";
	foreach my $group ( keys %intron_groups ) {
	  foreach my $intron ( @{$intron_groups{$group}} ) {
	    print "$group " . $intron->display_label . " " . $intron->score ."\n";
	  }
	  if ( scalar( @{$intron_groups{$group}} ) > 1 ) {
	    #  remove the lowest scoring one
	    my $intron = pop( @{$intron_groups{$group}} ) ;
	    print "Eliminating " . $intron->display_label . " " .
	      $intron->score . "\n";
	    unless (  $intron->display_label =~ /REMOVED/ ) {
	      $intron->display_label($intron->display_label."-REMOVED") ;
	      $removed++;
	    }
	  }	
	}
      }
      
      foreach my $intron ( @{$exon_intron->[$i]} ) {
	# only allow each intron to connect to 1 exon
	foreach my $exon ( @{$intron_exon->{$intron->display_label}} ) {
	  if ( $intron->end > $exons->[$exon]->end ) {
	    next ;
	  }
	  # store the possible paths as a hash (splice)variants
	  $variants->{$i}->{$intron->display_label} = 1;
	  $variants->{$intron->display_label}->{$exon} = 1;
	  # check 
	  if ( $intron->end > $exons->[$exon]->end ) {
	    throw(" exon $i start end " . $intron->end ." - " . $exons->[$exon]->start );
	  }
	}
      }
    }
  }
  
  if ($strict &! $removed ) {
    warn ( "Cannot simplify this gene any more ".
	   "EXON 0:". $exons->[0]->start ." - ". 
	   $exons->[0]->end ." - ". 
	   $exons->[0]->strand ."\n" );
    if ( $self->recursive_limit < $self->MAX_RECURSIONS ) {
      $self->recursive_limit ( $self->recursive_limit * 10 ) ; 
      warn('Upping recursive limit to ' . $self->recursive_limit . ' see if it helps');
    } else {
      warn("Giving up on " .
	   "EXON 0:".  $exons->[0]->seq_region_name ." - ". 
	   $exons->[0]->start ." - ". 
	   $exons->[0]->end ." - ". 
	   $exons->[0]->strand ."\n" );
      return "Give up";
    }
  }
  
  # work out all the possible paths given the features we have
  my $result;
  my %paths;
  for ( my $i = 0 ; $i < scalar(@$exons) ; $i ++ ) {
    $limit = 0;
    $result = $self->ProcessTree($variants,$i,undef,\%paths );
    if ($result eq "ERROR"){
      warn("Could not process  cluster $i trying again with simpler cluster\n");
      return undef;
    }
  }
  return $result;
}

=head2 model_cluster

    Title        :   model_cluster
    Usage        :   $self->model_cluster($models,$strand);
    Returns      :   2D array ref of exons and intron features
    Args         :   Array ref of  exons and intron features
                 :   Integer indicating strand
    Description  :   Clusters the initial models by start end
                 :   orders the models in each cluster by score

=cut

sub model_cluster {
  my ($self,$models,$strand) = @_;
  my @clusters;
  # sort them by the start of the fist exon ( the first array element )
  my @models = sort { $a->[0]->start <=> $b->[0]->start }  @$models ;
  # $model->[0] = 1st exon 
  # $model->[-2] = last exon once the score has been popped off
  # the clusters have an array of features where each feature 
  # is stored at an index = score, makes it easy to pull out the 
  # highest scoring models in each non overlapping cluster
  
  foreach my $model ( @models ) {
    # get the score ( last array element ) 
    my $score  = pop(@$model);
     my $clustered = 0;
    foreach my $cluster ( @clusters ) {
      # do they overlap?
      if ( $model->[0]->start <= $cluster->{'end'} 
	   &&  $model->[-2]->end >= $cluster->{'start'}) {
	# Expand the cluster
	$cluster->{'start'} = $model->[0]->start 
	  if $model->[0]->start < $cluster->{'start'};
	$cluster->{'end'} = $model->[-2]->end
	  if $model->[-2]->end   > $cluster->{'end'}; 
	push @{$cluster->{'models'}->[$score]}, $model;
	$clustered = 1;
      }
    }
    unless ($clustered) {
      my $cluster;
      push @{$cluster->{'models'}->[$score]}, $model;
      $cluster->{'start'} = $model->[0]->start;
      $cluster->{'end'}   = $model->[-2]->end;
      $cluster->{'strand'} = $strand;
      push @clusters, $cluster;
    }
  }
  return \@clusters;
}


=head2 merge_exons

    Title        :   merge_exons
    Usage        :   $self->merge_exons($gene)
    Returns      :   Array ref of Bio::EnsEMBL::Exon
    Args         :   Bio::EnsEMBL::Gene
    Description  :   Merges adjacent exons where the intron is covered by repeats or
                 :   is very small

=cut


# lets us merge exons with tiny  introns between them  unless they contain an intron
sub merge_exons {
  my ( $self, $gene, $strand) = @_;
  my @exons;
  next unless $gene->get_all_Transcripts->[0];
  foreach my $exon ( @{$gene->get_all_Transcripts->[0]->get_all_Exons} ) {
    push @exons, clone_Exon($exon);
  }
  @exons =  sort { $a->start <=> $b->start } @exons;

  for ( my $i = 1 ; $i <= $#exons ; $i ++ ) {
    my $exon = $exons[$i];
    my $prev_exon = $exons[$i-1];
    my $intron_count = 0;
    # is the intron tiny?    
    my @introns = @{$self->fetch_intron_features($prev_exon->end,$exon->start)};
    
     # we know it lies across the boundary does it lie within the 2 exons?
    foreach my $intron ( @introns ) {
      print "INTRON " . $intron->start . " " . $intron->end . " " , $intron->strand ." " , $intron->score ."\n";
      # ignore non consensus introns at this point
      next if $intron->display_label =~ /non-consensus/ ;
      if (   $intron->start > $prev_exon->start &&
	    $intron->end <  $exon->end && 
	 $intron->strand == $strand ){
	$intron_count++;
      }
    }
    # is it covered by repeats?
    # remove very small introns if there is no direct evidence for them
    if ( $exon->start - $prev_exon->end <= 20  && $intron_count == 0)   {
      $exon->start($prev_exon->start);
      splice(@exons,$i-1,1);
      $i--;
      next;
    }
    # dont merge long introns even if they are covered with repeats
    # next if $repeat_slice->length > 200 ;
    # dont merge introns if there is evidence for splicing
    next if $intron_count;
    # is the intron covered by a repeat?
    my @repeats = @{$self->repeats};


    my $repeat_coverage = 0;
    # so the repeats are now non-overlapping blocks ( if there is more than one )
    foreach my $repeat ( @repeats ) {
      next unless $repeat->start <= $exon->start && $repeat->end >= $prev_exon->end;
      last if $repeat->start > $exon->start;
      $repeat->start($prev_exon->end) if  $repeat->start < $prev_exon->end;
      $repeat->end($exon->start) if $repeat->end > $exon->start;
      $repeat_coverage+= $repeat->end - $repeat->start;
    }
    $repeat_coverage /= ($exon->start - $prev_exon->end ) ;
     print   " Intron " . $exon->start ." " .  $prev_exon->end . " coverage  $repeat_coverage \n";
    # splice the exons together where repeat coverage is > 99%
    if ($repeat_coverage >= 0.99   ) {
      print "MERGING EXONS Intron " . $exon->start ." " .  $prev_exon->end . " coverage  $repeat_coverage \n";
      $exon->start($prev_exon->start);
      # combine the exon scores when we merge the exons
      $exon->get_all_supporting_features->[0]->score
	($exon->get_all_supporting_features->[0]->score + $prev_exon->get_all_supporting_features->[0]->score) if 
	  $exon->get_all_supporting_features->[0] && $prev_exon->get_all_supporting_features->[0];
      splice(@exons,$i-1,1);
      $i--;
    }
  }
  
  return \@exons;
}

=head2 dna_2_simple_features
    Title        :   dna_2_simple_features
    Usage        :   $self->dna_2_simple_features($start,$end)
    Returns      :   None
    Args         :   Int start
                 :   Int end
    Description  :   Fetches all dna_align_features from the intron db that lie within
                 :   the range determined by start and end, collapses them down into a 
                 :   non redundant set and builds a Bio::EnsEMBL::SimpleFeature to 
                 :   represent it, then stores it in $self->intron_features

=cut

sub dna_2_simple_features {
  my ($self,$start,$end) = @_;
  my @sfs;
  my %id_list;
  my $intron_slice = $self->intron_slice_adaptor->fetch_by_region
    ('toplevel',
     $self->chr_slice->seq_region_name,
     $start,
     $end,
    );
  # featch all the dna_align_features for this slice by logic name
  my @reads;
  print STDERR  "Fetching reads with logic names: ";
  foreach my $logic_name ( @{$self->LOGICNAME} ) {
    print STDERR "$logic_name ";
    push @reads, @{$intron_slice->get_all_DnaAlignFeatures($logic_name)}; 
  }
  print STDERR "\n";
  # fetch them all if no logic name is Supplied
  if (scalar( @{$self->LOGICNAME} ) == 0 ) {
    @reads =  @{$intron_slice->get_all_DnaAlignFeatures()}; 
  }
  print STDERR "Got " . scalar(@reads) . " reads\n";
  while ( scalar @reads > 0 ) {
    my $read = pop(@reads);
    my $type = 'consensus';
    $type = 'non-consensus' if ( $read->hseqname =~ /\:NC$/ ) ;
    $read = $read->transfer($self->chr_slice);
    my @ugfs = $read->ungapped_features;
    @ugfs = sort { $a->start <=> $b->start } @ugfs;
    next if ( scalar(@ugfs) == 0 ) ;
    for ( my $i = 0 ; $i < scalar(@ugfs) - 1 ; $i++ ) {
      # one read can span several exons so make all the features 
      # cache them by internal boundaries
      # we use . to deliminate entries in our paths so dont allow them in the seq_region_name or it wont work
      my $name = $read->seq_region_name;
      $name =~ s/\./*/g;
      my $unique_id = $name . ":" . 
	$ugfs[$i]->end . ":" .
	  $ugfs[$i+1]->start . ":" . 
	    $read->strand .":$type";
      $id_list{$unique_id} ++;
    }
  }
  print STDERR "Got " . scalar( keys %id_list ) . " collapsed introns\n";
  # collapse them down and make them into simple features
  foreach my $key ( keys %id_list ) {
    my @data = split(/:/,$key) ;
    my $sf = Bio::EnsEMBL::SimpleFeature->new
      (-start => $data[1],
       -end => $data[2],
       -strand => $data[3],
       -slice => $self->chr_slice,
       -analysis => $self->analysis,
       -score =>  0,
       -display_label => "INTRON-" .$key,
      );
    $sf->score($id_list{$key});
    push @sfs , $sf;
  }
  # sort them
  @sfs = sort {$a->start <=> $b->start} @sfs;
  $self->intron_features(\@sfs);
  print STDERR "Got " . scalar(@sfs) . " simple features\n";
  return;
}

=head2 fetch_intron_features

    Title        :   fetch_intron_features
    Usage        :   $self->fetch_intron_features($start,$end)
    Returns      :   Array ref of Bio::EnsEMBL::SimpleFeature
    Args         :   Int start
                 :   Int end
    Description  :   Accesses the pre computed simple features representing introns
                 :   Filters out non consensus models that overlap consensus models

=cut

sub fetch_intron_features {
  my ($self,$start,$end) = @_;
  my @chosen_sf;
  my @filtered_introns;
  my @sfs =  @{$self->intron_features};
  # sfs is a sorted array
  foreach my $intron ( @sfs ) {
    next unless $intron->start <= $end && $intron->end >= $start;
    last if $intron->start > $end;
    push @chosen_sf, $intron;
  }
  INTRON: foreach my $intron ( @chosen_sf) {
    if ($intron->display_label =~ /non-consensus/ ) {
      # check it has no overlap with any consensus introns
      # unless it out scores a consensus intron
      foreach my $i ( @chosen_sf) {
	unless ($i->display_label =~ /non-consensus/ ) {
	  if ($intron->end > $i->start && $intron->start < $i->end && $intron->strand == $i->strand ) {
	    next INTRON if $intron->score <= $i->score;
	  }
	}
      }
      push @filtered_introns, $intron;
    } else {
      push @filtered_introns, $intron;
    }
  }
  return \@filtered_introns;
}


=head2 make_repeat_blocks
    Title        :   make_repeat_blocks
    Usage        :   $self->($repeats)
    Returns      :   Array ref of merged Bio::EnsEMBL::Repeat_Feature
    Args         :   Array ref of Bio::EnsEMBL::Repeat_Feature
    Description  :   Merges adjacent repeat features

=cut


sub make_repeat_blocks {
  my ($self,$repeats_ref) = @_;
  my @repeats = sort { $a->start <=> $b->start }@$repeats_ref;
  # merge repeat blocks
  for ( my $j = 1 ; $j <= $#repeats ; $j ++ ) { 
    if ( $repeats[$j]->start <= $repeats[$j-1]->end+1 ){
     # print "merging repeat $j " . $repeats[$j]->start . "-"  . $repeats[$j]->end. " " . $repeats[$j-1]->start ."-" . $repeats[$j-1]->end ."\n";
      $repeats[$j-1]->end($repeats[$j]->end) if  $repeats[$j]->end > $repeats[$j-1]->end ;
      splice(@repeats,$j,1);
      $j--;
    }
  }  
  print STDERR " got " . scalar(@repeats) . " repeat blocks after merging\n";
  return \@repeats;
}

##################################################################
# Containers

sub recursive_limit {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_recursive_limit} = $val;
  }

  return $self->{_recursive_limit};
}

sub repeat_feature_adaptor {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_rfa} = $val;
  }

  return $self->{_rfa};
}

sub gene_slice_adaptor {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_gsa} = $val;
  }

  return $self->{_gsa};
}

sub intron_slice_adaptor {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_isa} = $val;
  }

  return $self->{_isa};
}

sub feature_cash {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_feature_cash} = $val;
  }

  return $self->{_feature_cash};
}

sub prelim_genes {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_prelim_genes} = $val;
  }

  return $self->{_prelim_genes};
}

sub repeats {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_repeats} = $val;
  }

  return $self->{_repeats};
}


sub chr_slice {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_chr_slice} = $val;
  }

  return $self->{_chr_slice};
}

sub intron_features {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_introns} = $val;
  }

  return $self->{_introns};
}

####################################
# config variable holders
####################################

sub INTRON_DB {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_INTRON_DB'} = $value;
  }
  
  if (exists($self->{'_CONFIG_INTRON_DB'})) {
    return $self->{'_CONFIG_INTRON_DB'};
  } else {
    return undef;
  }
}

sub OUTPUT_DB {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_OUTPUT_DB'} = $value;
  }
  
  if (exists($self->{'_CONFIG_OUTPUT_DB'})) {
    return $self->{'_CONFIG_OUTPUT_DB'};
  } else {
    return undef;
  }
}


sub MODEL_DB {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MODEL_DB'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MODEL_DB'})) {
    return $self->{'_CONFIG_MODEL_DB'};
  } else {
    return undef;
  }
}


sub LOGICNAME {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_LOGICNAME'} = $value;
  }
  
  if (exists($self->{'_CONFIG_LOGICNAME'})) {
    return $self->{'_CONFIG_LOGICNAME'};
  } else {
    return undef;
  }
}

sub RETAINED_INTRON_PENALTY {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_RETAINED_INTRON_PENALTY'} = $value;
  }
  
  if (exists($self->{'_CONFIG_RETAINED_INTRON_PENALTY'})) {
    return $self->{'_CONFIG_RETAINED_INTRON_PENALTY'};
  } else {
    return undef;
  }
}


sub MIN_INTRON_SIZE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MIN_INTRON_SIZE'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MIN_INTRON_SIZE'})) {
    return $self->{'_CONFIG_MIN_INTRON_SIZE'};
  } else {
    return undef;
  }
}


sub BEST_SCORE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_BEST_SCORE'} = $value;
  }
  
  if (exists($self->{'_CONFIG_BEST_SCORE'})) {
    return $self->{'_CONFIG_BEST_SCORE'};
  } else {
    return undef;
  }
}


sub OTHER_NUM {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_OTHER_NUM'} = $value;
  }
  
  if (exists($self->{'_CONFIG_OTHER_NUM'})) {
    return $self->{'_CONFIG_OTHER_NUM'};
  } else {
    return undef;
  }
}

sub OTHER_ISOFORMS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_OTHER_ISOFORMS'} = $value;
  }
  
  if (exists($self->{'_CONFIG_OTHER_ISOFORMS'})) {
    return $self->{'_CONFIG_OTHER_ISOFORMS'};
  } else {
    return undef;
  }
}

sub MODEL_LN {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MODEL_LN'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MODEL_LN'})) {
    return $self->{'_CONFIG_MODEL_LN'};
  } else {
    return undef;
  }
}

sub BAD_MODELS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_BAD_MODELS'} = $value;
  }
  
  if (exists($self->{'_CONFIG_BAD_MODELS'})) {
    return $self->{'_CONFIG_BAD_MODELS'};
  } else {
    return undef;
  }
}

sub MAX_NUM {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MAX_NUM'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MAX_NUM'})) {
    return $self->{'_CONFIG_MAX_NUM'};
  } else {
    return undef;
  }
}

sub MAX_RECURSIONS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MAX_RECURSIONS'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MAX_RECURSIONS'})) {
    return $self->{'_CONFIG_MAX_RECURSIONS'};
  } else {
    return undef;
  }
}

sub MIN_SINGLE_EXON {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MIN_SINGLE_EXON'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MIN_SINGLE_EXON'})) {
    return $self->{'_CONFIG_MIN_SINGLE_EXON'};
  } else {
    return undef;
  }
}

sub SINGLE_EXON_CDS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_SINGLE_EXON_CDS'} = $value;
  }
  
  if (exists($self->{'_CONFIG_SINGLE_EXON_CDS'})) {
    return $self->{'_CONFIG_SINGLE_EXON_CDS'};
  } else {
    return undef;
  }
}

sub SINGLE_EXON_MODEL {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_SINGLE_EXON_MODEL'} = $value;
  }
  
  if (exists($self->{'_CONFIG_SINGLE_EXON_MODEL'})) {
    return $self->{'_CONFIG_SINGLE_EXON_MODEL'};
  } else {
    return undef;
  }
}

sub STRICT_INTERNAL_SPLICE_SITES{
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_STRICT_INTERNAL_SPLICE_SITES'} = $value;
  }
  
  if (exists($self->{'_CONFIG_STRICT_INTERNAL_SPLICE_SITES'})) {
    return $self->{'_CONFIG_STRICT_INTERNAL_SPLICE_SITES'};
  } else {
    return undef;
  }
}

sub STRICT_INTERNAL_END_EXON_SPLICE_SITES {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_STRICT_INTERNAL_END_EXON_SPLICE_SITES'} = $value;
  }
  
  if (exists($self->{'_CONFIG_STRICT_INTERNAL_END_EXON_SPLICE_SITES'})) {
    return $self->{'_CONFIG_STRICT_INTERNAL_END_EXON_SPLICE_SITES'};
  } else {
    return undef;
  }
}


1;
