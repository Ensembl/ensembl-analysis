=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

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

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::RefineSolexaGenes;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils ;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils ;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranslationUtils;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::RefineSolexaGenes;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw (rearrange);
use Bio::DB::Sam;


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
  # explicitly attach the ref db
  $self->db($self->get_dbadaptor("REFERENCE_DB"));
  $self->intron_slice_adaptor($self->get_dbadaptor($self->INTRON_DB)->get_SliceAdaptor)
    if $self->INTRON_DB;
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
  
  my $chr_slice = $self->db->get_SliceAdaptor->fetch_by_region
    (
     'toplevel',
     $slice->seq_region_name
    );
  $self->chr_slice($chr_slice);

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

  if ( $self->INTRON_BAM_FILE ) {

    my $sam = Bio::DB::Sam->new(   -bam => $self->INTRON_BAM_FILE,
				   -autoindex => 1,
				   -expand_flags => 1,
			       );
    $self->throw("Bam file " . $self->INTRON_BAM_FILE . "  not found \n") unless $sam; 
    my $count = 0;
    my $segment = $sam->segment($slice->seq_region_name,$slice->start,$slice->end);
    $self->throw("Bam file segment not found for slice " .  $slice->name . "\n")
      unless $segment;
    # need to seamlessly merge here with the dna2simplefeatures code
    $self->bam_2_intron_features($segment);
  } else {
    # pre fetch all the intron features
    $self->dna_2_intron_features($slice->start,$slice->end);
  }
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
    # hack taking out weeny models
    next GENE if $gene->get_all_Transcripts->[0]->length < 300;
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
      foreach my $exon ( @exons ) {
	print STDERR "EXTRAEXON: " . 
	  $exon->seq_region_name ." " .
	    ($exon->start +20) ." " .
	      ($exon->end -20)." " .
		( $exon->end - $exon->start -40)  ."\n"
		  if $exon->{"_extra"} ;
      }
#      if ( $self->extra_exons ) {
#	foreach my $exon ( @{$self->extra_exons} ) {
#	  print STDERR "OTHEREXON: " . 
#	    $exon->seq_region_name ." " .
#	      ($exon->start +20) ." " .
#		($exon->end -20)." " .
#		  ( $exon->end - $exon->start -40)  ."\n";
#	}
#      }
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
	  next unless $intron->length >  $self->MIN_INTRON_SIZE;
	  next unless $intron->length <= $self->MAX_INTRON_SIZE;
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
	      if ( $intron->hseqname =~ /non canonical/ ) {
		push @left_nc_introns, $intron if $intron->score > 1;
	      } else {
		push @left_c_introns, $intron;
	      }
	    }
	    if ( $intron->start >= $exon->start ) {
	      if ( $intron->hseqname =~ /non canonical/ ) {
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
	  print STDERR "\t" . $intron->start . " " . $intron->end . " " . $intron->strand . " " . $intron->hseqname . " " . $intron->score . "\n";
	  # becasue we make a new exons where we have a reatained intron to 
	  # stop circular references we need to allow the final 
	  # intron splicing out of the exon to be used more than once
	  # by each new exon in fact
	  $intron_count{$intron->hseqname}++ unless $retained_intron;
	  $intron_hash{$intron->hseqname} = $intron;
	  # only use each intron twice once at the end and once at the start of
	  # an exon
	  # exon_intron links exons to the intron on their right ignoring strand
	  push @{ $exon_intron[$i]}  , $intron if $intron->end > $exon->end;
	  # intron exon links introns to exons on their right ignoring strand
	  if ( $intron->start < $exon->start ) {
	    push @{$intron_exon{$intron->hseqname}} , $i;
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
	      if ( $intron->hseqname  =~ /non canonical/ ) {
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
      my $fwd_t_len;
      $fwd_t_len = $fwd_tran->translation->genomic_end - $fwd_tran->translation->genomic_start 
	if $fwd_tran->translateable_seq;
      #print STDERR "FWD t length $fwd_t_len\n";
      my $rev_t =  new Bio::EnsEMBL::Transcript(-EXONS => [$rev_exon]);
      my $rev_tran = compute_translation(clone_Transcript($rev_t));
      my $rev_t_len;
      $rev_t_len = $rev_tran->translation->genomic_end - $rev_tran->translation->genomic_start
	if $rev_tran->translateable_seq;;
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
      if ( $exon_use_hash{$exon_use}  ) {
	$gene->biotype('bad');
      }
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
      if ( $feature =~ /canonical$/ ) {
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
	my @introns;
	my $intron_count = 0;
	my $intron_score = 0;
	my $exon_score = 0;
	my $non_con_introns = 0;
	my @new_exons;
	# make an array containing cloned exons
	for ( my $i = 0 ; $i < scalar(@$model) ; $i++ ) {
	  unless ( $model->[$i]->isa("Bio::EnsEMBL::DnaDnaAlignFeature") ) {	
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
	  if ( $model->[$i]->isa("Bio::EnsEMBL::DnaDnaAlignFeature") ) {
	    my $intron = $model->[$i];
	    next unless $intron->strand == $strand;
	    next unless $new_exons[$i-1] && $new_exons[$i+1];
	    push @introns,$intron;
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
	    print "Adding INTRON " . $intron->hseqname ." total score now $intron_score \n";
	    # keep tabs on how many fake introns we have
	    $non_con_introns++ if $intron->hseqname =~ /non canonical/;
	  }
	}
	next MODEL unless $intron_count;
	
	# trim padding from the start and end exons
	$new_exons[0]->start($new_exons[0]->start + 20) if $new_exons[0]->length  > 20;
	$new_exons[-1]->end ($new_exons[-1]->end  - 20) if $new_exons[-1]->length > 20;
	
	# make it into a gene
	my @modified_exons;
	foreach my $exon ( @new_exons ) {
	  next if $exon->isa("Bio::EnsEMBL::DnaDnaAlignFeature");
	  push @modified_exons, clone_Exon($exon);
	}
	if ( $strand == 1 ) {
	  @modified_exons = sort { $a->start <=> $b->start } @modified_exons;
	} else {
	  @modified_exons = sort { $b->start <=> $a->start } @modified_exons;
	}
	# make it into a gene
	my $t =  new Bio::EnsEMBL::Transcript(-EXONS => \@modified_exons);
	# check for dna
	my $s = $t->seq->seq ;
	my $Ns =  $s =~  s/N//g;
	if( length($t->seq->seq) == $Ns ){
	  $self->throw("There does not appear to be ay DNA in the database, transcript seq is all N's\n");
	}

	# add a translation 
	my $initial_tran = compute_translation(clone_Transcript($t));

	
	# trim UTR
	my $tran = $self->prune_UTR($initial_tran,\@introns);
	
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

sub prune_UTR {
  my ($self,$transcript,$introns) = @_;
  unless ( $self->TRIM_UTR ) {
    return $transcript;
  }
  unless ( $transcript->translateable_seq ) {
    return $transcript;
  }
  # otherwise trim the UTR according to the values set out in the config
  my %intron_hash;
  print STDERR "Got " . scalar(@{$introns}) . " introns - hashing...";
  foreach my $intron ( @{$introns} ) {
    my $key = $intron->start  .":". $intron->end .":". $intron->strand;
    print "INTRON $key\n";
    $intron_hash{$key} = $intron;
  }
  print STDERR " done \n Got ";
  my @new_fivep;
  my @new_threep;
  my @new_exons;
  my @features;
  my @exons = sort {$a->start <=> $b->start }  @{$transcript->get_all_Exons};

  # put everything into the features array
  push @features, @exons;
  for ( my $i =0 ; $i < $#exons  ; $i++ ) {
    my $key = ($exons[$i]->end) .":". ($exons[$i+1]->start ) . ":" . $exons[$i]->strand;

    if ( my $intron = $intron_hash{$key}  ) {
      push @features, $intron;
    }
  }
  @features = sort { $a->start <=> $b->start } @features;
  # so now we should have an array of alternating introns and exons
  print "=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

Transcript " . $transcript->stable_id ." " . 
  $transcript->seq_region_name ." " . 
    $transcript->start ." " . 
      $transcript->end ." " . 
	$transcript->strand ." " . 
	  scalar(@{$transcript->get_all_Exons}) ."

";
  throw("Something is wrong we are missing " . scalar(@{$introns}) ." introns " . scalar(@exons) . "  exons " . scalar(@features) . " exons and introns\n")
    unless scalar(@features) == (scalar(@exons) * 2) -1 ;
  my $average_intron = 0;
  my $intron_count = 0;
  # leave single exon genes alone for now
  if ( scalar(@features) == 1 or scalar(@{$transcript->get_all_translateable_Exons}) == 1 )  {
    # lets strip the UTR
    return  $self->modify_transcript($transcript,$transcript->get_all_translateable_Exons);
  }
  # first calculate the average
  foreach my $f ( @features ) {
    if ( $f->isa("Bio::EnsEMBL::DnaDnaAlignFeature")) {
      $average_intron += $f->score;
      $intron_count++;
    }
  }
  $average_intron /= $intron_count;

  foreach my $f ( @features ) {
    print $f->start . "\t" . $f->end ."\t$f\t";
    if ( $f->isa("Bio::EnsEMBL::DnaDnaAlignFeature")) {
      # make the intron features not overlap with the ends of the exons
      $f->start($f->start+1);
      $f->end($f->end-1);
      print $f->score;
      if ( $average_intron && ($f->score /  $average_intron) * 100 <= $self->REJECT_INTRON_CUTOFF ) {
	print " Potentially bad";
      }
    }
    print "\n";
  }
  throw("Something is wrong we are missing introns " . scalar(@exons) . "  exons  and $intron_count introns\n")
    unless $intron_count == scalar(@exons) -1 ;
  print  "Average intron depth = $average_intron \n";
  
  
  my @fivep;
  my @threep;
  my $coding =0 ;
  # need to account for strand
  @features = sort { $b->start <=> $a->start } @features if $transcript->strand == -1;

  for ( my $i = 0 ; $i <=  $#features ; $i += 2  ) {
    my $e = $features[$i];
    throw("Got a DNA align feature where I should have an exon\n")
      unless $e->isa("Bio::EnsEMBL::Exon");
    if ( $e->coding_region_start($transcript) ) {
      # first coding exon
      for ( my $j = 0 ; $j <= $i ; $j++ ) {
	push @fivep,$features[$j];
      }
      last;
    }
  }
  for ( my $i = $#features ; $i >= 0 ;  $i -= 2  ) {
    my $e = $features[$i];
    throw("Got a DNA align feature where I should have an exon\n")
      unless $e->isa("Bio::EnsEMBL::Exon");
    if ( $e->coding_region_end($transcript) ) {
	  # first coding exon
      for ( my $j = $i ; $j <= $#features ; $j++ ) {
	push @threep,$features[$j];
      }
      last;
    }
  }
  
  # want to start at last coding exon and work outwards so....
  @fivep = reverse @fivep;
  # now we should be good
  print "FIVE P \n";
  @new_exons = @{$transcript->get_all_translateable_Exons};
  my $fivep_cds = shift(@new_exons);
  my $threep_cds = pop(@new_exons);
  my $fiveplen;
  my $threeplen;
  my $fivepc = 0 ;
  my $threepc = 0 ;
  my $nmd;
  
  # FIVE PRIME RULES
  
 FIVEP: for ( my $i = 0 ; $i <= $#fivep ; $i++ ) {
    my $f =  $fivep[$i];
    print $f->start . "\t" . $f->end ."\t$f\t";
    if ( $i == 0 ) {
      throw("First feature is not an exon \n")
	unless $f->isa("Bio::EnsEMBL::Exon");
      # UTR starts in this exon - how long is it?
      my $cds_start = $f->coding_region_start($transcript);
      $cds_start = $f->coding_region_end($transcript)  if $transcript->strand == -1;
      throw("First coding exon has no CDS \n") unless $cds_start;
      print "CDS START $cds_start\t";
      $fiveplen = $cds_start - $f->start +1 if $transcript->strand == 1;
      $fiveplen = $f->end - $cds_start   +1 if $transcript->strand == -1;
      # is the coding exon too long
      if ( $fiveplen > $self->MAX_5PRIME_LENGTH ) {
	# replace it with the cds
	@new_fivep = ($fivep_cds);
	print " 5p too long $fiveplen \n";
	last FIVEP;
      }
      push @new_fivep,$f;
      $fivepc++;
    } else {
      if (  $f->isa("Bio::EnsEMBL::Exon") ) {
	$fivepc++;
	$fiveplen+= $f->end - $f->start +1;
	# does it make the UTR too long?
	if ( $fiveplen > $self->MAX_5PRIME_LENGTH ) {
	  # dont add it
	  print " 5p too long $fiveplen \n";
	  last FIVEP;
	}
	# is it too many exons?
	if ( $fivepc > $self->MAX_5PRIME_EXONS ) {
	  # dont add it
	  print " too many 5p  $fivepc cut them all as we are not sure \n";
	  @new_fivep = ($fivep_cds);
	  last FIVEP;
	}
	push @new_fivep,$f;
      }
    }
    # Does the intron score well enough to include the exon
    # apply rules and add successful exons into the mix
    if ( $f->isa("Bio::EnsEMBL::DnaDnaAlignFeature")) {
      print $f->score;
      if ( $average_intron && ($f->score /  $average_intron) *100 <= $self->REJECT_INTRON_CUTOFF ) {
	print " Rejecting on score cutoff " . $f->score ." vs  $average_intron\n";
	# dont add any more 
	last FIVEP;
      }
    }
    print "\n";
  }
  
  # three P
  print "THREE P \n";
 THREEP:   for ( my $i = 0 ; $i <= $#threep ; $i++ ) {
    my $f = $threep[$i];
    print $f->start . "\t" . $f->end ."\t$f\t";
    if ( $i == 0 ) {
      throw("First feature is not an exon \n")
	unless $f->isa("Bio::EnsEMBL::Exon");
      # UTR starts in this exon - how long is it?
      my $cds_end = $f->coding_region_end($transcript);
      $cds_end = $f->coding_region_start($transcript)  if $transcript->strand == -1;
      throw("last coding exon has no CDS \n") unless $cds_end;
      print "CDS END $cds_end\t";
      $threeplen = $cds_end - $f->start +1 if $transcript->strand == -1;
      $threeplen = $f->end - $cds_end   +1 if $transcript->strand == 1;
      # is the coding exon too long
      if ( $threeplen > $self->MAX_3PRIME_LENGTH ) {
	# replace it with the cds
	@new_threep = ($threep_cds);
	print " 3p too long $threeplen \n";
	last THREEP;
      }
      push @new_threep,$f;
      $nmd = $threeplen ;
      $threepc++;
    } else {
      if (  $f->isa("Bio::EnsEMBL::Exon") ) {
	# does it break the NMD rule?
	if ( $nmd > 55 ) {
	  print " splice is after $nmd bp from stop codon - rejected on NMD rule of maximum 55 bp \n";
	  @new_threep = ($threep_cds);
	  last THREEP;
	}
	$threepc++;
	$threeplen+= $f->end - $f->start +1;
	# does it make the UTR too long?
	if ( $threeplen > $self->MAX_3PRIME_LENGTH ) {
	  # dont add it
	  print " 3p too long $threeplen \n";
	  last THREEP;
	}
	# is it too many exons?
	if ( $threepc > $self->MAX_3PRIME_EXONS ) {
	  # dont add it
	  print " too many 3p  $threepc cut them all as we are not sure \n";
	  @new_threep = ($threep_cds);
	  last THREEP;
	}
	push @new_threep,$f;
      }
    }
    # Does the intron score well enough to include the exon
    # apply rules and add successful exons into the mix
    if ( $f->isa("Bio::EnsEMBL::DnaDnaAlignFeature")) {
      print $f->score;
      if ( $average_intron && ($f->score /  $average_intron) * 100 <= $self->REJECT_INTRON_CUTOFF ) {
	print " Rejecting on score cutoff " . $f->score ." vs  $average_intron\n";
	# dont add any more 
	last THREEP;
      }
    }
    print "\n";
  }
  
  push @new_exons, @new_fivep;
  push @new_exons, @new_threep;
  print " New transript has " . scalar(@new_exons) , " exons\n";
  my @clones;
  foreach my $e ( @new_exons ) {
    throw("Not is not an exon " . $e->start ." " . $e->end . " $e\n")
      unless $e->isa("Bio::EnsEMBL::Exon");
    push @clones, clone_Exon($e);
  }
  @clones = sort { $a->start <=> $b->start } @clones;
  @clones =  reverse(@clones) if $transcript->strand == -1;
  return $self->modify_transcript($transcript,\@clones);
}


sub modify_transcript {
  my ($self,$tran,$exons) = @_;
  my $cds_start = $tran->coding_region_start;
  $cds_start = $tran->coding_region_end if $tran->strand == -1;
  my $cds_end = $tran->coding_region_end;
  $cds_end = $tran->coding_region_start if $tran->strand == -1;
  print "CDS START END $cds_start  $cds_end \n";
  print "PHASE " . $tran->translation->start . " " . $tran->translation->end ."\n";
  my $t =  new Bio::EnsEMBL::Transcript(-EXONS => $exons);
  my $se;
  my $ee;
  foreach my $e ( @{$t->get_all_Exons} ) {
    $ee = $e if $e->start <= $cds_end &&  $e->end >= $cds_end;
    $se = $e if $e->start <= $cds_start &&  $e->end >= $cds_start;
  }
  my $ts;
  my $te;
  if ( $tran->strand == -1 ) {
    $ts =  $se->end - $cds_start+ 1;
    $te =  $ee->end - $cds_end  + 1;
  } else {
    $ts =   $cds_start - $se->start+ 1;
    $te =   $cds_end - $ee->start  + 1;
  }
  #  my $start_phase = $se->phase;
  my $translation =  new Bio::EnsEMBL::Translation->new( 
							-START_EXON => $se,
							-END_EXON   => $ee,
							-SEQ_START  => $ts,
							-SEQ_END    => $te,
						       );
  print "S-E $ts $te \n";#START PHASE $start_phase\n";
  print "GS " . $translation->genomic_start ." " . $translation->genomic_end ."\n";
  $t->translation($translation);
  # calculate_exon_phases($t,$start_phase);
  unless (  $tran->translation->seq eq $t->translation->seq ) {
    $self->throw("Translations do not match: Before " . $tran->translation->seq ."\nAfter  " . 	  $t->translation->seq ."\n");
  }
  return $t;
}

sub write_output{
  my ($self) = @_;

  my $outdb = $self->get_dbadaptor($self->OUTPUT_DB);
  my $gene_adaptor = $outdb->get_GeneAdaptor;
  my $intron_adaptor = $outdb->get_DnaAlignFeatureAdaptor;
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
  # store the intron features too
  if ( $self->WRITE_INTRONS ) {
    my $introns = $self->intron_features;
    my $fails = 0;
    my $total = 0;
    print "Writing " . scalar( @$introns ) ." introns\n";
    foreach my $intron ( @$introns ) {
      # make the introns so they dont overlap the exons
      $intron->start($intron->start+1);
      $intron->end($intron->end-1);
      eval {
	$intron_adaptor->store($intron);
      };    
      if ($@){
	warning("Unable to store intron!!\n$@");
	$fails++;
      }
      $total++;
    }
    if ($fails > 0) {
      throw("Not all introns could be written successfully " .
	    "($fails fails out of $total)");
    }
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
	# Throw out exons that have retained introns for a start
	next if  $exons->[$i]->{'retained'};
      }
      if ( $strict > 1 ) {
	# group all the introns by exon pairs
	my %intron_groups;
	foreach my $intron ( @{$exon_intron->[$i]} ) {
	  next if  $intron->hseqname =~ /REMOVED/  ;
	  push @{$intron_groups{$i}}, $intron;
	}
	# now lets sort these groups by score
	foreach my $group ( keys %intron_groups ) {
	  @{$intron_groups{$group}}  = sort {$b->score <=> $a->score} @{$intron_groups{$group}};
	}
	# now lets see what they look like
	print STDERR "EXON $i:". $exons->[$i]->start ." - ". 
	  $exons->[$i]->end ." - ". 
	    $exons->[$i]->strand ."\n";
	foreach my $group ( keys %intron_groups ) {
	  foreach my $intron ( @{$intron_groups{$group}} ) {
	    print STDERR "$group " . $intron->hseqname . " " . $intron->score ."\n";
	  }
	  if ( scalar( @{$intron_groups{$group}} ) > 1 ) {
	    #  remove the lowest scoring one
	    my $intron = pop( @{$intron_groups{$group}} ) ;
	    print STDERR "Eliminating " . $intron->hseqname . " " .
	      $intron->score . "\n";
	    unless (  $intron->hseqname =~ /REMOVED/ ) {
	      $intron->hseqname($intron->hseqname."-REMOVED") ;
	      $removed++;
	    }
	  }	
	}
      }
      
      foreach my $intron ( @{$exon_intron->[$i]} ) {
	# only allow each intron to connect to 1 exon
	foreach my $exon ( @{$intron_exon->{$intron->hseqname}} ) {
	  if ( $intron->end > $exons->[$exon]->end ) {
	    next ;
	  }
	  # store the possible paths as a hash (splice)variants
	  $variants->{$i}->{$intron->hseqname} = 1;
	  $variants->{$intron->hseqname}->{$exon} = 1;
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
  my $ec = scalar(@exons) ;
  my $extra_exons = $self->extra_exons;
  # the extra exon is a list of start end coords of the spliced intron sections
  # ie: end:start:end:start where the 1st and last coords are anchors to tie it 
  # into our rough model both must match before we can try and add any potentialy 
  # novel exons in
  foreach my $key ( keys %$extra_exons ) {
    my @coords = split(/:/,$key);
    my $start_anchor = shift(@coords);
    my $end_anchor = pop(@coords);
 #   print "START AND END $start_anchor  $end_anchor \n";
    # do the anchors lie within the model?
    foreach my $exon ( @exons ) {
      if ( $start_anchor <= $exon->end && 
	   $start_anchor >= $exon->start ) {
	$start_anchor = -1;
      }
      if ( $end_anchor <= $exon->end && 
	   $end_anchor >= $exon->start ) {
	$end_anchor = -1;
      }
    }
    if ( $start_anchor == -1 && $end_anchor == -1 ) {
      # now to make new the exon(s)
      for ( my $i = 0 ; $i <= $#coords ; $i += 2 ) {
	my $left = $coords[$i];
	my $right = $coords[$i+1];
	my $extra = $self->make_exon(undef,$left,$right,$extra_exons->{$key},$key );
	$extra->{"_extra"} = 1;
	push @exons,$extra;
      }
    }
  }
  
  
  @exons =  sort { $a->start <=> $b->start } @exons;
  # want to get rid of any overlapping exons
  while  ( $ec != scalar(@exons) ) {
    $ec = scalar(@exons);
    for ( my $i = 1 ; $i < scalar(@exons) ; $i++ ) {
      my $left_exon = $exons[$i-1];
      my $right_exon = $exons[$i];
      # do they overlap 
      if ( $left_exon->start <= $right_exon->end && 
	   $left_exon->end >= $right_exon->start ) {
	# merge them 
	if (   $right_exon->end >= $left_exon->end &&
	       $right_exon->start <= $left_exon->start ){
	  $left_exon->{"_extra"} = 0;
	}
	$left_exon->start($right_exon->start) 
	  if $right_exon->start < $left_exon->start;
	$left_exon->end($right_exon->end) 
	  if $right_exon->end > $left_exon->end;
	# get rid of right exon
	splice(@exons,$i,1);
	$i-- ;
	@exons =  sort { $a->start <=> $b->start } @exons;
      }
    }
  }
  
  

  for ( my $i = 1 ; $i <= $#exons ; $i ++ ) {
    my $exon = $exons[$i];
    my $prev_exon = $exons[$i-1];
    my $intron_count = 0;
    # is the intron tiny?    
    my @introns = @{$self->fetch_intron_features($prev_exon->end,$exon->start)};
    
     # we know it lies across the boundary does it lie within the 2 exons?
    foreach my $intron ( @introns ) {
 #     print "INTRON " . $intron->start . " " . $intron->end . " " , $intron->strand ." " , $intron->score ."\n";
      # ignore non consensus introns at this point
      next if $intron->hseqname =~ /non canonical/ ;
      if (   $intron->start > $prev_exon->start &&
	    $intron->end <  $exon->end && 
	 $intron->strand == $strand ){
	$intron_count++;
      }
    }
    # remove very small introns if there is no direct evidence for them
    if ( $exon->start - $prev_exon->end <= 20  && $intron_count == 0)   {
      $exon->start($prev_exon->start);
      splice(@exons,$i-1,1);
      $i--;
      next;
    }
  }
  
  return \@exons;
}

=head2 bam_2_intron_features
    Title        :   bam_2_intron_features
    Usage        :   $self->bam_2_intron_features($segment)
    Returns      :   None
    Args         :   Bam file segement
    Description  :   Fetches all alignments from the bam file segment, collapses them down into a 
                 :   non redundant set and builds a Bio::EnsEMBL::DnaDnaAlignFeature to 
                 :   represent it, then stores it in $self->intron_features
                 :   analyses splice sites for consensus and non consensus splices as this data is 
                 :   not stored in the BAM.
                 :   Also checks for small exons defined by a single read splicing at least twice
                 :   stores any additional exons found this way in $self->extra_exons

=cut

sub bam_2_intron_features {
  my ($self,$segment) = @_;
  my $slice_adaptor = $self->gene_slice_adaptor;
  my @ifs;
  my $extra_exons;
  my %id_list;
  my %read_groups;
  if ( scalar(@{$self->GROUPNAME} > 0 } ) {
    my @groups = @{$self->GROUPNAME};
    print "Limiting to read groups ";
    foreach my $group ( @groups ) {
      print " $group";
      $read_groups{$group} = 1;
    }
    print "\n";
  }
  my $iterator = $segment->features(-iterator=>1);
 READ:  while (my $read = $iterator->next_seq) {
    # need to recreate the ungapped features code as the
    # auto splitting code does not seem to work with > 2 features
    if ( $self->GROUPNAME ) {
      next unless ($read_groups{$read->get_tag_values('RG')}) ;
    }
    my @mates = sort { $a->[2] <=> $b->[2] } @{$self->ungapped_features($read)};

    # if mates > 2 then we have a possibility of adding in some extra exons into our rough models
    # as the read has spliced into and out of an exon
    # lets make them unique
    if ( scalar(@mates) > 2 ) {
      my $string;
      for ( my $i = 0 ; $i <= $#mates  ; $i++ ) {

	my $start  = $mates[$i]->[2];
	my $end    = $mates[$i]->[3];
	my $hstrand = $read->strand;
	$string .= $start .":" if $i > 0 ;
	$string .= $end .":" if $i < $#mates ;
      }
      $extra_exons->{$string} ++;
    }
    my $strand = $read->target->strand;
   # print "\nREAD " . $read->cigar_str;
    my $offset;
    for ( my $i = 0 ; $i <= $#mates  ; $i++ ) {
   #   print "\n";
      # intron reads should be split according to the CIGAR line
      # the default split function seems to ad
      # we want the ungapped features to make our introns
      my $name   = $read->seq_id;
      # we dont allow . in the seq region name as we use them to delimit our paths
      $name =~ s/\./*/g;   
      my $start  = $mates[$i]->[2];
      my $end    = $mates[$i]->[3];
      my $cigar  = $mates[$i]->[4];
      my $hstrand = $read->strand;
   #   print "$name $start $end $strand $hstrand $cigar\t";
   #   print $read->get_tag_values('FIRST_MATE') ."\t";
      next if $i == $#mates;
      my $unique_id = $name . ":" . 
	$mates[$i]->[3] . ":" .
	  $mates[$i+1]->[2] . ":" . 
	    $strand ;
      $id_list{$unique_id} ++;
   #   print "$unique_id";
    }
   # print "\n";
  }

  # collapse them down and make them into simple features
  foreach my $key ( keys %id_list ) {
    my @data = split(/:/,$key) ;
    my $length =  $data[2] - $data[1] -1;
    next unless $length > 0 ;
    my $name = $data[0]. ":" . ($data[1]+1) . ":" . ($data[2] -1) .":" . $data[3] .":";

    my $if = Bio::EnsEMBL::DnaDnaAlignFeature->new
      (
       -start => $data[1],
       -end => $data[2],
       -strand => $data[3],
       -hstart => 1,
       -hend => $length,
       -hstrand => 1,
       -slice => $self->chr_slice,
       -analysis => $self->analysis,
       -score =>  $id_list{$key},
       -hseqname => "$name",
       -cigar_string => $length ."M",
      );
    my $canonical = 1;
    # figure out if its cannonical or not
    my $left_splice = $slice_adaptor->fetch_by_region('toplevel',
						      $if->seq_region_name,
						      $if->start+1,
						      $if->start+2,
						      $if->strand
						     );
    my $right_splice = $slice_adaptor->fetch_by_region('toplevel',
						       $if->seq_region_name,
						       $if->end-2,
						       $if->end-1,
						       $if->strand
						      );
    #  print "KEY $key " . $if->score ."\n";;
    #  print "LEFT  " . $left_splice->start ." " . $left_splice->end  ." " . $left_splice->strand ." " . $left_splice->seq . "\n";
    #  print "RIGHT " . $right_splice->start ." " . $right_splice->end  ." " . $right_splice->strand  ." " . $right_splice->seq ."\n\n";
    
    
    if ( $left_splice->seq eq 'NN' && $right_splice->seq eq 'NN' ) {
      warn("Cannot find dna sequence for " . $key .
	   " this is used in detecting non cannonical splices\n");
    } else {
      # is it cannonical
      if ( $if->strand  == 1 ) {
	#	print "Splice type " . $left_splice->seq ."-".  $right_splice->seq ." ";
	# is it GTAG?
	unless ( $left_splice->seq eq 'GT' && $right_splice->seq eq 'AG' ) {
	  $canonical = 0;
	}
      } else {
	#	print "Splice type " . $right_splice->seq ."-".  $left_splice->seq ." ";
	# is it GTAG?
	unless ( $right_splice->seq eq 'GT' && $left_splice->seq eq 'AG' ) {
	  $canonical = 0;
	}
      }
    }
    if ( $canonical ) {
       $if->hseqname($if->hseqname."canonical");
    } else {
      $if->hseqname($if->hseqname."non canonical");
    }
     push @ifs , $if;
  }
  # sort them
  @ifs = sort {$a->start <=> $b->start} @ifs;
  $self->intron_features(\@ifs);
  $self->extra_exons($extra_exons);
  print STDERR "Got " . scalar(@ifs)  . " unique introns  " ;
  print STDERR " and " . scalar(keys %$extra_exons) . " potential novel exons from " . $self->INTRON_BAM_FILE . "\n";
  return;
}

sub ungapped_features {
  my ($self,$read) = @_;
  my @ugfs;
  
  my $string = $read->cigar_str;
  my $start = $read->start;
  my $end = $read->end;
 # print "THINGS $start $end $string\n";
  my @pieces = ( $string =~ /(\d*[MDN])/g );
  foreach my $piece (@pieces) {
    my ($length) = ( $piece =~ /^(\d*)/ );
    if( $length eq "" ) { $length = 1 }
    if( $piece =~ /M$/ ) {
      #
      # MATCH
      #
      my ( $qstart, $qend);
      $qstart = $start;
      $qend = $start + $length - 1;
      $start = $qend + 1;
      
      my $ugf;
      $ugf->[0] = $read->query->name;
      $ugf->[1] = $read->seq_id;
      $ugf->[2] = $qstart;
      $ugf->[3] = $qend;
      $ugf->[4] = $length."M";
      push @ugfs, $ugf;
#      print "UNGAPPED " .$ugf->[2] .
#	" " . $ugf->[3] . " " . $ugf->[4] ."\n";
    } elsif( $piece =~ /N$/ ) {
      #
      # INSERT
      #
      $start += $length;
    } elsif( $piece =~ /D$/ ) {
      #
      # DELETION
      #
      $start += $length;
      
    } else {
      throw( "Illegal cigar line $string!" );
    }
  }
  return \@ugfs;
}

=head2 dna_2_extra_exons
    Title        :   dna_2_extra_exons
    Usage        :   $self->dna_2_extra_exons($start,$end)
    Returns      :   None
    Args         :   Int start
                 :   Int end
    Description  :   Fetches all dna_align_features from the small intron db that lie within
                 :   the range determined by start and end, collapses them down into a 
                 :   non redundant set and builds extra exons to account for short exons
                 :   missed by genomic alignemnt of reads

=cut

sub dna_2_extra_exons {
  my ($self,$start,$end) = @_;
  my $extra_introns = $self->intron_features;
  my @ifs;
  push @ifs, @$extra_introns if $extra_introns;
  my $rough_genes = $self->prelim_genes;
  my %exon_list;
  my %id_list;
  
  # fetch all the dna_align_features for this slice by logic name
  my @extra_reads;
  # look for extra introns from ESLA
  my $small_intron_slice = $self->small_intron_slice_adaptor->fetch_by_region
    ('toplevel',
     $self->chr_slice->seq_region_name,
     $start,
     $end,
    );  
  print STDERR  "Fetching ESLA reads with logic names: ";
  foreach my $logic_name ( @{$self->LOGICNAME} ) {
    print STDERR "$logic_name ";
    push @extra_reads, @{$small_intron_slice->get_all_DnaAlignFeatures($logic_name)}; 
  }
  print STDERR "\n";
  # fetch them all if no logic name is Suplied
  if (scalar( @{$self->LOGICNAME} ) == 0 ) {
    push @extra_reads ,  @{$small_intron_slice->get_all_DnaAlignFeatures()}; 
  }
  
  print STDERR "Got " . scalar(@extra_reads) . " extra  reads\n";
  # process extra reads and assign them to rough models
  while ( scalar @extra_reads > 0 ) {
    my $read = pop(@extra_reads);
    my $type = 'canonical';
    my $roughid;
    $type = 'non canonical' if ( $read->hseqname =~ /\:NC$/ ) ;
    $read = $read->transfer($self->chr_slice);
    # assign a rough model to the reads
    my @ugfs = $read->ungapped_features;
    @ugfs = sort { $a->start <=> $b->start } @ugfs;
    next if ( scalar(@ugfs) == 0 ) ;
    # first and or last ugf should overlap an exon
  ROUGH:  foreach my $rough ( @$rough_genes ) {
      foreach my $exon ( @{$rough->get_all_Exons} ) {
	if( ( $exon->start <= $ugfs[0]->end &&  $exon->end >= $ugfs[0]->start ) or 
	    ( $exon->start <= $ugfs[-1]->end &&  $exon->end >= $ugfs[-1]->start ) ) {
	  #ONLY USE THIS READ WITH THIS MODEL
	  $roughid = $rough->stable_id;
	  last ROUGH;
	}
      }
    }
    unless ( $roughid ) {
      print STDERR "Ignoring read " . $read->hseqname . " cannot pair it with a rough model\n";
    }
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
      if  ( $i > 0 && $i < scalar(@ugfs) - 1 ) {
	# if the read splices into and out of an exon, store that exon in case it is not 
	# found in the rough model but is an extra short exon determined by ExonerateSolexaLocalAlignment
	# check that it is introns on both sides and not inserts, need to be bigger than the min intron length
	next unless $ugfs[$i+1]->start - $ugfs[$i]->end > $self->MIN_INTRON_SIZE;
	next unless $ugfs[$i]->start - $ugfs[$i-1]->end > $self->MIN_INTRON_SIZE;
	my $name = $read->seq_region_name;
	$name =~ s/\./*/g;
	my $unique_id = $name . ":" . 
	  $ugfs[$i]->start . ":" .
	    $ugfs[$i]->end . ":" . 
	      $read->strand .":$type";
	$exon_list{$roughid}{$unique_id} = $ugfs[$i];	
      }
    }
  }
  
  # collapse down the intron feaures and make them into simple features
  foreach my $key ( keys %id_list ) {
   # print "KEY $key\n";
    my @data = split(/:/,$key) ;
    my $length =  $data[2] - $data[1] -1;
    next unless $length > 0 ;
    my $name = $self->chr_slice->seq_region_name . ":" . ($data[1]+1) . ":" . ($data[2] -1) .":" . $data[3] .":".$data[4];
    my $if = Bio::EnsEMBL::DnaDnaAlignFeature->new
      (
       -start => $data[1],
       -end => $data[2],
       -strand => $data[3],
       -hstart => 1,
       -hend => $length,
       -hstrand => 1,
       -slice => $self->chr_slice,
       -analysis => $self->analysis,
       -score =>  $id_list{$key},
       -hseqname => $name,
       -cigar_string => $length ."M",
      );
    push @ifs , $if;
  }
  # sort them
  @ifs = sort {$a->start <=> $b->start} @ifs;
  $self->intron_features(\@ifs);
  print STDERR "Got " . scalar(@ifs) . " intron features\n";
  # make the exon features
  my @extra_exons;
  foreach my $rough ( keys %exon_list ) {
    foreach my $key ( keys %{$exon_list{$rough}} ) {
      my $extra_exon = $self->make_exon($exon_list{$rough}{$key});
      $extra_exon->{"_extra"} = 1;
      $extra_exon->{"_model"} = $rough;
      push @extra_exons, $extra_exon;
    }
  }
  print STDERR "Got " . scalar(@extra_exons) . " extra exons\n";
  # make the exon features
  $self->extra_exons(\@extra_exons);
  return;
}


=head2 dna_2_intron_features
    Title        :   dna_2_intron_features
    Usage        :   $self->dna_2_intron_features($start,$end)
    Returns      :   None
    Args         :   Int start
                 :   Int end
    Description  :   Fetches all dna_align_features from the intron db that lie within
                 :   the range determined by start and end, collapses them down into a 
                 :   non redundant set and builds a Bio::EnsEMBL::DnaAlignFeature to 
                 :   represent it, then stores it in $self->intron_features
                 :   also checks for small exons defined by a single read splicing at least twice
                 :   stores any additional exons found this way in $self->extra_exons

=cut

sub dna_2_intron_features {
  my ($self,$start,$end) = @_;
  my @ifs;
  my %id_list;
  my $rough_genes = $self->prelim_genes;
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
    my $type = 'canonical';
    $type = 'non canonical' if ( $read->hseqname =~ /\:NC$/ ) ;
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
    my $length =  $data[2] - $data[1] -1;
    next unless $length > 0;
    my $name = $self->chr_slice->seq_region_name . ":" . ($data[1]+1) . ":" . ($data[2] -1) .":" . $data[3] .":".$data[4];
    my $if = Bio::EnsEMBL::DnaDnaAlignFeature->new
      (
       -start => $data[1],
       -end => $data[2],
       -strand => $data[3],
       -hstart => 1,
       -hend => $length,
       -hstrand => 1,
       -slice => $self->chr_slice,
       -analysis => $self->analysis,
       -score =>  $id_list{$key},
       -hseqname => $name,
       -cigar_string => $length ."M",
      );
    push @ifs , $if;
  }
  # sort them
  @ifs = sort {$a->start <=> $b->start} @ifs;
  $self->intron_features(\@ifs);
  print STDERR "Got " . scalar(@ifs) . " intron features\n";
  return;
}

=head2 fetch_intron_features

    Title        :   fetch_intron_features
    Usage        :   $self->fetch_intron_features($start,$end)
    Returns      :   Array ref of Bio::EnsEMBL::DnaAlignFeature
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
    if ($intron->hseqname =~ /non canonical/ ) {
      # check it has no overlap with any consensus introns
      # unless it out scores a consensus intron
      foreach my $i ( @chosen_sf) {
	unless ($i->hseqname =~ /non canonical/ ) {
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



=head2 make_exon
    Title        :   pad_exons
    Usage        :   $self->($ungapped_feature)
    Returns      :   Bio::EnsEMBL::Exon
    Args         :   Bio::EnsEMBL::FeaturePair 
    Description  :   Takes an ungapped feature, pads it and builds a 
                 :   Exon from it 
=cut

sub make_exon {
  my ($self,$ugf,$start,$end,$score,$display_id) = @_;
  if ( $ugf) {
    $start = $ugf->start;
    $end = $ugf->end;
    $display_id = $ugf->display_id;
    $score = $ugf->score;
  }
  my $padded_exon =  create_Exon
    (
     $start - 20,
     $end + 20 ,
     -1,
     -1,
     -1,
     $self->analysis,
     undef,
     undef,
     $self->chr_slice,
    );
  # dont let it fall of the slice because of padding
  $padded_exon->start(1) if $padded_exon->start <= 0;
  $padded_exon->end($self->chr_slice->length - 1) 
    if $padded_exon->end >= $self->chr_slice->length;
  
  my $feat = new Bio::EnsEMBL::DnaDnaAlignFeature
    (-slice    => $self->chr_slice,
     -start    => $padded_exon->start,
     -end      => $padded_exon->end,
     -strand   => -1,
     -hseqname => $display_id,
     -hstart   => 1,
     -hstrand  => 1,
     -hend     => $padded_exon->length,
     -analysis => $self->analysis,
     -score    => $score,
     -cigar_string => $padded_exon->length.'M');
  my @feats;
  push @feats,$feat;
  $padded_exon->add_supporting_features(@feats);
  return $padded_exon;
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

sub small_intron_slice_adaptor {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_sisa} = $val;
  }

  return $self->{_sisa};
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

sub extra_exons {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_extra_exons} = $val;
  }

  return $self->{_extra_exons};
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


sub MAX_INTRON_SIZE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MAX_INTRON_SIZE'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MAX_INTRON_SIZE'})) {
    return $self->{'_CONFIG_MAX_INTRON_SIZE'};
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

sub GROUPNAME {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_GROUPNAME'} = $value;
  }
  
  if (exists($self->{'_CONFIG_GROUPNAME'})) {
    return $self->{'_CONFIG_GROUPNAME'};
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

sub INTRON_BAM_FILE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_INTRON_BAM_FILE'} = $value;
  }
  
  if (exists($self->{'_CONFIG_INTRON_BAM_FILE'})) {
    return $self->{'_CONFIG_INTRON_BAM_FILE'};
  } else {
    return undef;
  }
}


sub WRITE_INTRONS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_WRITE_INTRONS'} = $value;
  }
  
  if (exists($self->{'_WRITE_INTRONS'})) {
    return $self->{'_WRITE_INTRONS'};
  } else {
    return undef;
  }
}

sub TRIM_UTR {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_TRIM_UTR'} = $value;
  }
  
  if (exists($self->{'_CONFIG_TRIM_UTR'})) {
    return $self->{'_CONFIG_TRIM_UTR'};
  } else {
    return undef;
  }
}


sub MAX_3PRIME_EXONS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MAX_3PRIME_EXONS'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MAX_3PRIME_EXONS'})) {
    return $self->{'_CONFIG_MAX_3PRIME_EXONS'};
  } else {
    return undef;
  }
}


sub MAX_3PRIME_LENGTH {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MAX_3PRIME_LENGTH'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MAX_3PRIME_LENGTH'})) {
    return $self->{'_CONFIG_MAX_3PRIME_LENGTH'};
  } else {
    return undef;
  }
}


sub MAX_5PRIME_EXONS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MAX_5PRIME_EXONS'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MAX_5PRIME_EXONS'})) {
    return $self->{'_CONFIG_MAX_5PRIME_EXONS'};
  } else {
    return undef;
  }
}


sub MAX_5PRIME_LENGTH {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MAX_5PRIME_LENGTH'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MAX_5PRIME_LENGTH'})) {
    return $self->{'_CONFIG_MAX_5PRIME_LENGTH'};
  } else {
    return undef;
  }
}


sub REJECT_INTRON_CUTOFF {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_REJECT_INTRON_CUTOFF'} = $value;
  }
  
  if (exists($self->{'_CONFIG_REJECT_INTRON_CUTOFF'})) {
    return $self->{'_CONFIG_REJECT_INTRON_CUTOFF'};
  } else {
    return undef;
  }
}





1;
