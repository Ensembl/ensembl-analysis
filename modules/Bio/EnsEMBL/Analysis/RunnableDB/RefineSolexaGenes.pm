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
  $self->recursive_limit(50000);
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
  $self->prediction_exon_adaptor($self->db->get_PredictionExonAdaptor);
  $self->repeat_feature_adaptor($self->db->get_RepeatFeatureAdaptor);
  $self->gene_slice_adaptor($self->get_dbadaptor($self->MODEL_DB)->get_SliceAdaptor);

  # fetch splice matricies for making ab-initio introns where we have no intron support
  $self->donor_matrix($self->load_matrix($self->DONOR_MATRIX));
  $self->acceptor_matrix($self->load_matrix($self->ACCEPTOR_MATRIX));

  # want a slice and a full chromsome to keep everything on in the same coords
  my $slice = $self->fetch_sequence($self->input_id);
  my $gene_slice =  $self->gene_slice_adaptor->fetch_by_region
    ( 
     'toplevel',
     $slice->seq_region_name,
     $slice->start,
     $slice->end,
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
  foreach my $gene ( @{$gene_slice->get_all_Genes} ) {
    # put them on the chromosome
    $gene = $gene->transfer($chr_slice);
    push @prelim_genes,$gene;
  }
  print STDERR scalar(@prelim_genes) . " genes found \n";
  # deterine strandedness ( including splitting merged genes )
  my @genes = @{$self->strand_separation(\@prelim_genes)};
  print STDERR scalar(@genes) . " genes once they have been separated by strand \n";
  $self->prelim_genes(\@genes);
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
    
    my %intron_count;
    my %intron_hash;
    my @exon_intron;
    my @exon_prev_intron;
    my %intron_exon;
    my $variants;
    my $strand = $gene->strand;
    my @trans;

    my $most_real_introns = 0;
    my $highest_score = 0;
    print STDERR $gene->stable_id. " : " .  $gene->start . " " . $gene->end . ":\n";
    # merge exons to remove little artifactual introns
    my @exons =  @{$self->merge_exons($gene)};
    my $exon_count =  $#exons;
    my @fake_introns;
    
 EXON:   for ( my $i = 0 ; $i <= $exon_count ; $i ++ ) {
      my $exon = $exons[$i];
      my $retained_intron;
      my $left_introns = 0;
      my $right_introns = 0;
      $exon->{'left_mask'} = 0;
      $exon->{'right_mask'} = $exon->length;
      print STDERR "$i : " . $exon->start . " " . $exon->end . ":\n";
      # make intron features by collapsing the dna_align_features
      my @introns = @{$self->dna_2_simple_features($exon->seq_region_start,$exon->seq_region_end)};
      my @filtered_introns;
      my $intron_overlap;
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
	}
	push @filtered_introns, $intron;
      }
      
    INTRON:  foreach my $intron ( @filtered_introns ) {
	print STDERR "\t" . $intron->start . " " . $intron->end . " " . $intron->strand . " " . $intron->display_label . "\n";
	# becasue we make a new exons where we have a reatained intron to 
	# stop circular references we need to allow the final 
	# intron splicing out of the exon to be used more than once
	# by each new exon in fact
	$intron_count{$intron->display_label}++ unless $retained_intron;

	$intron_hash{$intron->display_label} = $intron;
	# only use each intron twice once at the end and once at the start of
	# an exon
	$self->exon_mask($exon,$intron);
        next if $intron_count{$intron->display_label} &&  $intron_count{$intron->display_label} > 2;
	#print STDERR "USING intron " . $intron->display_label . "\n";
	#print STDERR "E-I\n"  if $intron->end > $exon->end;
	#print STDERR "I-E\n"  if $intron->start < $exon->start;
	# exon_intron links exons to the intron on their right ignoring strand
	push @{ $exon_intron[$i]}  , $intron if $intron->end > $exon->end;
	# intron exon links introns to exons on their right ignoring strand
	if ( $intron->start < $exon->start ) {
	  push @{$intron_exon{$intron->display_label}} , $i;
	  # exon_prev_intron links exons to introns on their left ignoring strand
	  push @{ $exon_prev_intron[$i]}  , $intron ;
	}

	# intron is within the exon - this is not a true exon but a retained intron
	if (  $intron->start > $exon->start && $intron->end < $exon->end && $intron->length > 50 ) {
	  #print STDERR  "retained " . $intron->start ." " . $exon->start ." ". $intron->end ." ". $exon->end ."\n";
	  # dont have circular references to exons or the paths
	  # will be infinite so clone this exon instead
	  my $new_exon = clone_Exon( $exon );
	  # chop it up a bit so it no longer overlaps the other introns
	  #print STDERR  "TRIMMING EXON \n";
	  $exon->end($intron->start + 10 );
	  $new_exon->start( $intron->end -10 );
	  # split the exon into 2 new exons and replace the old retained exon with the 2 new ones
	  my @replacements = ( $exon, $new_exon);
	  splice( @exons,$i,1,@replacements) ;
	  $exon_count++;
	  $i--;
	  next EXON;
	}
      }
    }

    # if we have no intron support we can guess at introns 
    if ( $self->ABINITIO_INTRONS ) {
    EXON:   for ( my $i = 1 ; $i <= $exon_count ; $i ++ ) {

	my $left_exon = $exons[$i-1];
	my $exon = $exons[$i];
	my $right_exon = $exon;
	my $intron_start = 0;
	my $intron_end = 0;
	
	# we may already know the splice sites for the left exon even if the intron 
	# that splices it doesnt join onto this gene...
	if ($exon_intron[$i-1]){
	  my $prev_intron = $exon_intron[$i-1][0];
	  $intron_start = $prev_intron->start;
	}

	# likewise right intron
	if ($exon_prev_intron[$i]){
	  # just pick the first one, lets not get carried away 
	  my $prev_intron = $exon_prev_intron[$i][0];
	  $intron_end = $prev_intron->end;
	}

	next if $intron_start && $intron_end;
	#print STDERR "FOUND intron start $intron_start and intron end $intron_end \n";
	# otherwise one is missing lets fake the intron
	# lets find the splice sites using a variety of cunning techniques
	if ($strand == 1) {
	  $intron_start =  $self->find_splice_sites( $left_exon,  $strand,'donor') unless $intron_start;
	  $intron_end   =  $self->find_splice_sites( $right_exon, $strand,'acceptor') unless $intron_end;
	} else {
	  $intron_start =  $self->find_splice_sites( $left_exon,  $strand,'acceptor') unless $intron_start;
	  $intron_end   =  $self->find_splice_sites( $right_exon, $strand,'donor') unless $intron_end;
	}
	# check for problems with introns 
	if ( $intron_end <= $intron_start ) {
	  warn("coudnt get splices for this intron\n");
	  $intron_start = $left_exon->end;
	  $intron_end = $right_exon->start;
	}
	# catch other problems
	if ( $intron_start <= $left_exon->start or 
	     $intron_end   >= $right_exon->end ) {
	  warn("This intron would result with an exon start > end $intron_start $intron_end \n");
	  $intron_start = $left_exon->end;
	  $intron_end = $right_exon->start;
	}
	
	# use existing values if you cannot find a splice site
	$intron_start = $left_exon->end unless  $intron_start;
	$intron_end = $right_exon->start unless $intron_end;
	# make the fake intron
	my $fake_intron = Bio::EnsEMBL::SimpleFeature->new 
	  (-start => $intron_start ,
	   -end => $intron_end  ,
	   -strand => $strand,
	   -slice => $self->chr_slice,
	   -analysis => $self->analysis,
	   -score => 0
	  ); 
	
	# add the fake intron into the hash structure
	my $key =  $fake_intron->start . ":" . $fake_intron->end . ":" . $fake_intron->strand;
	#print STDERR "INTRON $key\n";
	$fake_intron->display_label("INTRON-FAKE-$key");
	$intron_hash{$fake_intron->display_label} = $fake_intron;
	push @{ $exon_intron[$i-1]}  , $fake_intron if $fake_intron->end > $left_exon->end;
	push @{$intron_exon{$fake_intron->display_label}} , $i  if $fake_intron->start < $exon->start;
	$self->exon_mask($left_exon,$fake_intron);
	$self->exon_mask($right_exon,$fake_intron);
	
	if (  $exons[$i]->{'left_mask'} >=  $exons[$i]->{'right_mask'} ) {
	  print STDERR  "Problem gene " . $gene->stable_id ." Exon $i has overlapping introns skipping this gene\n";
	  next GENE;
	}
      }
    }

    next unless @exon_intron;

    # now lets make a hash of hashes holding which exons connect to which
    for ( my $i = 0 ; $i < scalar(@exons) ; $i ++ ) {
	if ( $exon_intron[$i] ) {
	foreach my $intron ( @{$exon_intron[$i]} ) {	
	  # only allow each intron to connect to 1 exon
	  foreach my $exon ( @{$intron_exon{$intron->display_label}} ) {
	    if ( $intron->end > $exons[$exon]->end ) {
	      next ;
	    }
	    # store the possible paths as a hash (splice)variants
	    $variants->{$i}->{$intron->display_label} = 1;
	    $variants->{$intron->display_label}->{$exon} = 1;
	    # check 
	    if ( $intron->end > $exons[$exon]->end ) {
	      throw(" exon $i start end " . $intron->end ." - " . $exons[$exon]->start );
	    }
	  }
	}
      }
    }
    
    # work out all the possible paths given the features we have
    my @paths;
  CLUSTER: for ( my $i = 0 ; $i < scalar(@exons) ; $i ++ ) {
      $limit = 0;
      my $result;
      $result = $self->ProcessTree($variants,$i);
      if ($result eq "ERROR"){
	print STDERR ("Could not process cluster\n");
	next GENE;
      }
      push( @paths, (keys %$result));
    }
    @paths = sort { length($a) <=> length($b) } @paths;
    
    # lets collapse redundant paths
    foreach ( my $i = scalar(@paths)-1 ; $i >= 0  ; $i-- ) {
      foreach my $pathb ( @paths ) {
	next if  $pathb eq $paths[$i]  ;
	next unless $paths[$i];
	if ( $pathb =~ /$paths[$i]/ ) {
	  splice(@paths,$i,1);
	}
      }
    }
    # paths are stored as text - turn them into arrays of features "models"
    my @models;
    foreach my $path ( @paths ) {
      my @model;
      foreach my $feature ( split(/\./,$path ) ) {
	if ( $feature =~ /INTRON/ ) {
	  push @model, $intron_hash{$feature};
	} else {
	  push @model, $exons[$feature];
	}
      }
      push @models,\@model;
    }

    # Now we cycle through all the models and turn them into genes

 MODEL:  foreach my $model ( @models ) {
      my $intron_count = 0;
      my $intron_score = 0;
      my $fake_introns = 0;
      my @new_exons;
      # make an array containing cloned exons
      for ( my $i = 0 ; $i < scalar(@$model) ; $i++ ) {
	unless ( $model->[$i]->isa("Bio::EnsEMBL::SimpleFeature") ) {	
	  my $new_exon = clone_Exon($model->[$i]);
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
	  $intron_count++;
	  $intron_score+= $intron->score;
	  # keep tabs on how many fake introns we have
	  $fake_introns++ if $intron->display_label =~ /FAKE/;
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
	# check it works
	if ( $exon->start > $exon->end ) {
	  warn("Exons start > exon end " . $exon->start . " " . $exon->end ."\n");
	  next MODEL;
	}
	
	push @modified_exons, create_Exon
	  (
	   $exon->start,
	   $exon->end,
	   -1,
	   -1,
	   $strand,
	   $exon->analysis,
	   undef,
	   undef,
	   $self->chr_slice,
	  );
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
      $tran->biotype('modified');
      $tran->{'_score'} =  $intron_score;
      $tran->{'_fake_introns'} =  $fake_introns ;
      unless ($self->ABINITIO_INTRONS) {
	# if we are not using fake exons compare the introns to 
	# how many there are in the (merged) rough model
	$tran->{'_fake_introns'} = $exon_count - $intron_count;
	$intron_count = $exon_count ;
      }
      print STDERR " EXON count $exon_count\n";
      $tran->{'_intron_count'} = $intron_count;
      $tran->{'_proportion_real_introns'} = int((($tran->{'_intron_count'} - $tran->{'_fake_introns'}) /$tran->{'_intron_count'}) *100);
      # make note of best scores and best introns
      if ( $tran->{'_proportion_real_introns'} > $most_real_introns ) {
	$most_real_introns = $tran->{'_proportion_real_introns'}
      }
      if ( $tran->{'_score'} > $highest_score ) {
	$highest_score = $tran->{'_score'}
      }
      push @trans, $tran;
    }
    # cluster on peptides to get the trans with the least exons where the 
    # peptides are identical - ie less UTR exons
    my @cluster = @{$self->peptide_cluster(\@trans)};
    next unless scalar(@cluster) > 0;
    # label by score / introns
    foreach my $tran ( @cluster ) {
      my ( $new_gene ) = @{convert_to_genes(($tran),$gene->analysis)};
      $new_gene->biotype('disguard');
      $new_gene->stable_id($gene->stable_id . "-" .
			   $tran->{'_score'} ."-" .
			   $tran->{'_fake_introns'} . "/" .
			   $tran->{'_intron_count'} . "-" .
			   $tran->{'_proportion_real_introns'});
      # add the biotypes depending on the scores
      $new_gene->biotype($self->OTHER_ISOFORMS) if $self->OTHER_ISOFORMS;
      
      if ( $self->BEST_INTRONS && $tran->{'_proportion_real_introns'} == $most_real_introns) {
	$new_gene->biotype($self->BEST_INTRONS);
      }
      if ( $self->BEST_SCORE && $tran->{'_score'} == $highest_score) {
	$new_gene->biotype($self->BEST_SCORE);
      }
      push@{$self->output} , $new_gene unless $new_gene->biotype eq 'disguard';
    }
  }
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
                 :   routes through a hashref, has a hard limit of 500,000 recursions
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
    $paths = $self->ProcessTree($hashref,$child,$sofar,$paths);
    if ( $paths eq "ERROR" ) {
      $limit = 0;
      return $paths; 
    }
  }
  push @{$paths->{$sofar}}, split(/\./,$sofar) if scalar(@node == 0);
  return $paths;
}

=head2 peptide_cluster
    Title        :   peptide_cluster
    Usage        :   $self->peptide_cluster(\@array);
    Returns      :   Array ref of Bio::EnsEMBL::Transcript
    Args         :   Array ref of Bio::EnsEMBL::Transcript
    Description  :   Collapses transcripts that have identical coding
                 :   sequences and exons, picks the one with the fewest UTR exons
=cut


sub peptide_cluster {
  my ($self,$array_ref) = @_;
  my @clusters;
  my @chosen;
  my %exon_hash;
  foreach my $tran ( @$array_ref ) {
    my $key = '';
    my @exons = @{$tran->get_all_translateable_Exons};
    unless ( scalar(@exons > 0 ))  {
      print ($tran->start." ". $tran->end ." has no translateable exons \n");
      return \@chosen;
    }
    for ( my $i = 1 ; $i <=$#exons ; $i++ ) {
      $key .= $exons[$i-1]->end . ":" . $exons[$i]->start . ":" ;
    }
    push @{$exon_hash{$key}} , $tran;
  }

  foreach my $key ( keys %exon_hash ) {
    my @trans =  sort {  scalar(@{$a->get_all_Exons}) <=> (@{$b->get_all_Exons}) } @{$exon_hash{$key}};
    foreach my $tran ( @{$exon_hash{$key}} ) {
    }
    push @chosen, @trans;
  }
  return \@chosen;
}

=head2 find_splice_sites
    Title        :   find_splice_sites
    Usage        :   $self->find_splice_sites($exon,$strand,$matrix_type);
    Returns      :   None
    Args         :   Bio::EnsEMBL::Exon
                 :   Integer ( strand of the gene )
                 :   String descibing the type of matrix to use (acceptor or donor)
    Description  :   Attempts to define splice sites for exons when no other data is
                 :   availibale. Uses a position specific weight matrix to score
                 :   potential splice sites, also looks for any overlapping genscan
                 :   exons that might have correct splice sites
=cut



sub find_splice_sites {
  my ( $self, $exon,  $strand, $matrix_type  ) = @_;

  my $pea = $self->prediction_exon_adaptor;
  my $matrix; 
  if ( $matrix_type eq 'donor' ) {
    $matrix =  $self->donor_matrix;
  } else {
    $matrix =  $self->acceptor_matrix;
  }
  
  # which end of the exon do we wish to trim - just do one

  print "GOT EXON " . $exon->seq_region_name . " " .  $exon->start . " "  . $exon->end ." $strand \n";
  print " looking for a splice - $matrix_type \n";

  my $genomic_slice;
  if ($strand == 1 ) {
    $genomic_slice = $self->db->get_SliceAdaptor->fetch_by_region('toplevel',
					   $exon->seq_region_name,
					   $exon->start,
					   $exon->end,
					   1);
  } else {
    $genomic_slice = $self->db->get_SliceAdaptor->fetch_by_region('toplevel',
					   $exon->seq_region_name,
					   $exon->start,
					   $exon->end,
					   -1);
  }

  my $highest_donor;
  my $highest_acceptor;
  # first work out the coverage and the max coverage

  # lets see if we have any genscan exons to help us with our splice sites
  my @genscan_exons = @{$pea->fetch_all_by_Slice_constraint($genomic_slice,"seq_region_strand = $strand") };
  my $genscan_exon;

  if ( scalar(@genscan_exons) == 1) {
    $genscan_exon = $genscan_exons[0];
  }
  # get the dna sequence so we can score it against the splice matrix
  my $dna = $genomic_slice->seq;
  my @dna = split(//,$dna);
  
  # donor
  my $scores = $self->score_pssm($exon,$matrix,\@dna,$matrix_type);
  
  my %best_scores;
  if ( $genscan_exon ) {
    if ( $matrix_type eq 'donor') {
    } else {
      $scores->[$genscan_exon->start - 2] += .9 if $genscan_exon->start > 0;
    }
  }

  for ( my $i = 0 ; $i <= $#dna ; $i++ ) {
    #store them
    # how about we weight scores by virtue of their closeness to the expected position ?
    # adding 3 to the ends of the masking to prevent <=0 bp exons 
    if ($strand == 1){
      $scores->[$i] = -1 if $i < ($exon->{'left_mask'}+1) or $i > ($exon->{'right_mask'}-1);
    } else {
      $scores->[$i] = -1 if $i < ($exon->length - ($exon->{'right_mask'}-1)-1) or $i > ( $exon->length - ($exon->{'left_mask'}+1) -1);
    }
 
    my $factor;
    if ( $matrix_type eq 'acceptor' ) {
      # thats .25 off the score for every base 4 bases = 1 poimt 
      $factor = abs($i -20) / 100 ;
    } else {
      $factor = abs(($#dna - $i) -20 ) / 100;
    }
    $factor = .99 if $factor > .99;

    $scores->[$i] -= $factor if $scores->[$i] && $scores->[$i] > -1;

    $best_scores{$i} = $scores->[$i] if $scores->[$i] && $scores->[$i] > -1;
  }

  my $count;
  my $pick;
  # pick the best scoring splice site
  foreach my $position ( keys %best_scores ) {
    if ( $pick ) {
      $pick = $position if $best_scores{$position} > $best_scores{$pick};
    } else {
      $pick = $position;
    }

  }

  unless ($pick ) {
    return;
  }

  # the rough models are always -ve strand rememer
  if ( $strand == -1 && $matrix_type eq 'donor' ) {
    return $exon->end - $pick +1;
  }
  if ( $strand == 1 && $matrix_type eq 'donor' ) {
    return $exon->start+ $pick-1 ;
  }
  if ( $strand == 1 && $matrix_type eq 'acceptor' ) {
    return $exon->start+$pick+1 ;
  }
  if ( $strand == -1 && $matrix_type eq 'acceptor' ) {
    return $exon->end-$pick-1 ;
  }
 
  return;
}

=head2 score_pssm
    Title        :   score_pssm
    Usage        :   $self->score_pssm($exon,$matrix,$dna,$type)
    Returns      :   Array ref of scores
    Args         :   Bio::EnsEMBL::Exon
                 :   Array ref containing scoring matrix
                 :   Array ref containing dna 
                 :   String describing the type of matrix ( acceptor or donor )
    Description  :   Scores a DNA sequence against a PSWM to identify possible 
                 :   splice sites. Returns an array of scores
=cut


sub score_pssm {
  my ( $self,$exon,$matrix,$dna,$type) = @_;
  my $scores;
  my $adjust;
  $adjust = 3 if $type eq 'donor';
  $adjust = 13 if $type eq 'acceptor';
  throw("Matrix not recognised\n") unless $matrix;
  # calculate the score for all the points in the exon 
  for ( my $i = 0 ; $i <= $exon->length - ( scalar(keys %{$matrix} ) ) ; $i++ ) {
    my $offset = 0;
    my $index = $i+$adjust;
    foreach my $key (sort { $a <=> $b }   keys %{$matrix} ) {
      $scores->[$index] += ($matrix->{$key}->{$dna->[$i+$offset]} / 100 );
      $offset++;
    }
    # adjust for length of matrix
    $scores->[$index] /=  scalar(keys %{$matrix} );
    # make sure it has a GT if its a donor or an AG if its an acceptor
    if ($type eq 'donor' ) {
      $scores->[$index] = 0 unless $dna->[$index] eq 'G';
      $scores->[$index] = 0 unless $dna->[$index+1] eq 'T';
    }
    if ($type eq 'acceptor' ) {
      $scores->[$index] = 0 unless $dna->[$index-1] eq 'A';
      $scores->[$index] = 0 unless $dna->[$index] eq 'G';
    }
  }
  return $scores;
}

=head2 load_matrix
    Title        :   load_matrix
    Usage        :   $self->load_matrix($file)
    Returns      :   Array ref - score matrix
    Args         :   File handle
    Description  :   Parses a tet file containing a scoring matrix and
                 :   loads it into a multidimensional  array ref
=cut


sub load_matrix {
  my ( $self,$file ) = @_;
  my $matrix;
  my @columns;
  open(FILE,$file) or $self->throw("cannot find matrix file $file \n");
  while (<FILE>) {
    chomp;
    next if $_ =~ /^#/;
    $_ =~ s/ //g;
    my @array = split(/\t/,$_);
    unless ( scalar(@columns) > 0 ) {
      @columns = @array;
    } else {
      for ( my $i = 1 ; $i <= $#columns ; $i++ ) {
	$matrix->{$array[0]}->{$columns[$i]} = $array[$i];
      }
    }
  }
  return $matrix;
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
  my ( $self, $gene) = @_;
  my @exons;
  foreach my $exon ( sort { $a->start <=> $b->start } @{$gene->get_all_Transcripts->[0]->get_all_Exons} ) {
    push @exons, clone_Exon($exon);
  }
  for ( my $i = 1 ; $i <= $#exons ; $i ++ ) {
    my $exon = $exons[$i];
    my $prev_exon = $exons[$i-1];
    my $intron_count = 0;
    # is the intron tiny?
    my @introns = @{$self->dna_2_simple_features($prev_exon->start,$exon->end)};
    
     # it must lie across the boundary and be within the two exons
    foreach my $intron ( @introns ) {
      if (  $intron->start > 0 && $intron->end < ( $exon->end - $prev_exon->start) &&
	    $intron->start < $prev_exon->length &&
	    $intron->end >  $prev_exon->length &&
	    $intron->strand == $gene->strand) {
	$intron_count++;
      }
    }
    # is it covered by repeats?
    my $repeat_slice = $self->intron_slice_adaptor->fetch_by_region
      ('toplevel',
       $gene->slice->seq_region_name,
       $prev_exon->end,
       $exon->start,
      );
    
    if ( $repeat_slice->length <= 20  && $intron_count == 0)   {
      $exon->start($prev_exon->start);
      splice(@exons,$i-1,1);
      $i--;
      next;
    }
    # dont merge long introns even if they are covered with repeats
    next if $repeat_slice->length > 200 ;
    # dont merge introns if there is evidence for splicing
    next if $intron_count;
    # is the intron covered by a repeat?
    my @repeats = sort { $a->start <=> $b->start } @{$self->repeat_feature_adaptor->fetch_all_by_Slice($repeat_slice)} ;

    # merge repeat blocks
    for ( my $j = 1 ; $j <= $#repeats ; $j ++ ) { 
      if ( $repeats[$j]->start <= $repeats[$j-1]->end+1 ){
	$repeats[$j-1]->end($repeats[$j]->end) if  $repeats[$j]->end > $repeats[$j-1]->end ;
	splice(@repeats,$j,1);
	$j--;
      }
    }

    my $repeat_coverage = 0;
    # so the repeats are now non-overlapping blocks ( if there is more than one )
    foreach my $repeat ( @repeats ) {
      $repeat->start(0) if  $repeat->start < 0;
      $repeat->end($repeat_slice->length) if $repeat->end > $repeat_slice->length;
      $repeat_coverage+= $repeat->end - $repeat->start;
    }
    
    $repeat_coverage /= $repeat_slice->length;
    # splice the exons together where repeat coverage is > 95%
    if ($repeat_coverage > 0.95   ) {
      $exon->start($prev_exon->start);
      splice(@exons,$i-1,1);
      $i--;
    }
  }
  
  return \@exons;
}
=head2 strand_separation
    Title        :   strand_separation
    Usage        :   $self->strand_separation(\@genes)
    Returns      :   Array ref of Bio::EnsEMBL::Gene
    Args         :   Array ref of Bio::EnsEMBL::Gene
    Description  :   Examines a gene model with repect to overlapping intron alignments.
                 :   Looks for places where the intron alignments switch strand
                 :   indicating that the model is in fact 2 genes on opposite
                 :   strands with overlapping UTR. It then splits the gene into 2
=cut


sub strand_separation {
  my ( $self,$genes ) = @_;
  my @separated_genes;
  my $strand;
 GENE: foreach ( my $i = 0 ; $i< scalar(@$genes) ; $i++ ) {
    my $total_strand;
    my $gene = $genes->[$i];
    my @intron_strand;
    # examine strand on an exon by exon basis
    my @exons =   @{$self->merge_exons($gene)};
    for ( my $e =1 ; $e <= $#exons ; $e ++ ) {
      my $exon = $exons[$e];
      $intron_strand[$e] = 0;
	my $prev_exon = $exons[$e-1];
	# do we have an intron between these exons;
	my @introns = @{$self->dna_2_simple_features($prev_exon->start,$exon->end)};
	# here is where we figure out how to split the joined genes...
	foreach my $intron ( @introns ) {
	  # intron has to actually join the 2 exons
	  next unless ( ( $intron->start > $prev_exon->start && 
			  $intron->start <= $prev_exon->end)  or
			( $intron->end >= $exon->start && 
			  $intron->end <  $exon->end));
	  # add the intron scores 
	  if ( $intron->score ) {
	    $intron_strand[$e] += $intron->strand * $intron->score;
	    $total_strand+= $intron->strand * $intron->score;
	  }
	}
    }
    # do we have *any* introns at all - only continue if we do
    unless ( $total_strand ) {
      print STDERR "No strand data to go on - skipping " . $gene->display_id ."\n";
      next GENE;
    }
    my $last_intron_score = 0;
    my $last_intron_position = 0;
    my %switch_strand;
    
    # ok do we have a clearly defined strandedness or do we have a split?
    for ( my $e = 1 ; $e <= $#exons ; $e ++ ) {
      if ($last_intron_score && ( $intron_strand[$e] > 1 && $last_intron_score < -1 ) 
	  or  ( $intron_strand[$e] < -1 && $last_intron_score > 1 ) ) {
	$switch_strand{'score'}++;
	$switch_strand{'start'} = $last_intron_position;
	$switch_strand{'end'} = $e;
	$switch_strand{'strand'} = 1 if $last_intron_score > 0;
	$switch_strand{'strand'} = -1 if $last_intron_score < 0;	  
      } 
      if ( $intron_strand[$e] ) {
	$last_intron_score = $intron_strand[$e];
	$last_intron_position = $e;
      }
    }
    
    
    if ($switch_strand{'score'} && $switch_strand{'score'} == 1 ) {
      print STDERR "We have a strand switch between " . $switch_strand{'start'} ." and " .   $switch_strand{'end'} ."\n";
      # lets make 2 copies of this transcript
      # one with introns 1 to switch_strand end -1
      # one with switch_strand start +1  to end
      # then delete the merged one
      my @modified_exons;
      my $strand = $switch_strand{'strand'} ;
      
      # left half
      for ( my $e = 0;  $e <= $switch_strand{'start'}  ; $e ++ ) {
	my $exon = $exons[$e];
	push @modified_exons,  create_Exon
	  (
	   $exon->start,
	   $exon->end,
	   -1,
	   -1,
	   -1,
	   $exon->analysis,
	   undef,
	   undef,
	   $gene->slice,
	  );
      }
      my $t =  new Bio::EnsEMBL::Transcript(-EXONS => \@modified_exons);
      my ( $new_gene ) = @{convert_to_genes(($t),$gene->analysis)};
      $new_gene->strand($strand);
      $new_gene->stable_id($gene->stable_id."A");
      push @separated_genes,$new_gene;
      
      my @modified_right_exons;
      # right half
      for ( my $e = $switch_strand{'end'} -1  ;  $e <= $#exons ; $e ++ ) {
	my $exon = $exons[$e];
	push @modified_right_exons,  create_Exon
	  (
	   $exon->start,
	     $exon->end,
	   -1,
	   -1,
	   -1,
	   $exon->analysis,
	   undef,
	   undef,
	   $gene->slice,
	  );
      }
      $t =  new Bio::EnsEMBL::Transcript(-EXONS => \@modified_right_exons);
      my ( $new_gene2 ) = @{convert_to_genes(($t),$gene->analysis)};
      $new_gene2->strand(-$strand);
      $new_gene2->stable_id($gene->stable_id."B");
      push @separated_genes,$new_gene2;
    }
    
    if ( $switch_strand{'score'} &&  $switch_strand{'score'} > 1 ) {	
      print STDERR "Strand is ambiguous, lets go with majority rule\n";
      $gene->strand(1) if $total_strand > 0;
      $gene->strand(-1) if $total_strand < 0;
      push @separated_genes,$gene;
    }
    
    unless ( $switch_strand{'score'}  ) {	
      $gene->strand(1) if $total_strand > 0;
      $gene->strand(-1) if $total_strand < 0;
      push @separated_genes,$gene;
    }
  }
  return \@separated_genes;
}

=head2 exon_mask
    Title        :   exon_mask
    Usage        :   $self->exon_mask($exon,$intron)
    Returns      :   None
    Args         :   Bio::EnsEMBL::Exon
                 :   Bio::EnsEMBL::SimpleFeature
    Description  :   records positions of the slice donor and acceptor to 
                 :   prevent overlapping introns resulting in an exon with
                 :   start < end
=cut


sub exon_mask {
  my ( $self,$exon,$intron ) = @_;
  if ( $intron->end >= $exon->start &&  $intron->start < $exon->start ) {
    $exon->{'left_mask'} = $intron->end - $exon->start if  $exon->{'left_mask'} < $intron->end - $exon->start;
  }
  if ( $intron->end > $exon->end &&  $intron->start <= $exon->end ) {
    $exon->{'right_mask'} = $intron->start - $exon->start if   $exon->{'right_mask'} > $intron->start - $exon->start;
  }
  return ;
}

=head2 dna_2_simple_features
    Title        :   dna_2_simple_features
    Usage        :   $self->dna_2_simple_features($start,$end)
    Returns      :   Array ref of Bio::EnsEMBL::SimpleFeature
    Args         :   Int start
                 :   Int end
    Description  :   Fetches all dna_align_features from the intron db that lie within
                 :   the range determined by start and end, collapses them down into a 
                 :   non redundant set and builds a Bio::EnsEMBL::SimpleFeature to 
                 :   represent it
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
  foreach my $read ( @{$intron_slice->get_all_DnaAlignFeatures($self->LOGICNAME)} ) {
    $read = $read->transfer($self->chr_slice);
    my @ugfs = $read->ungapped_features;
    next unless ( scalar(@ugfs) == 2 );
    @ugfs = sort { $a->start <=> $b->start } @ugfs;
    # cache them by internal boundaries
    my $unique_id = $read->seq_region_name . ":" . 
      $ugfs[1]->start . ":" .
	$ugfs[0]->end . ":" . 
	  $read->strand .":" . 
	    $self->LOGICNAME;
    $id_list{$unique_id} ++;
 }
  # collapse them down and make them into simple features
  foreach my $key ( keys %id_list ) {
    my $sf = $self->feature_cash->{$key};
    unless($sf) {
      my @data = split(/:/,$key) ;
      $sf = Bio::EnsEMBL::SimpleFeature->new
	(-start => $data[2],
	 -end => $data[1],
	 -strand => $data[3],
	 -slice => $self->chr_slice,
	 -analysis => $self->analysis,
	 -score =>  0,
	 -display_label => "INTRON-" .$key . "-". scalar( keys %{$self->feature_cash} ),
	);
      $self->feature_cash->{$key} = $sf;
    }
    $sf->score($id_list{$key});
    push @sfs , $sf;
  }
  return \@sfs;
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

sub donor_matrix {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_donor_matrix} = $val;
  }

  return $self->{_donor_matrix};
}

sub acceptor_matrix {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_acceptor_matrix} = $val;
  }

  return $self->{_acceptor_matrix};
}

sub prediction_exon_adaptor {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_pea} = $val;
  }
  return $self->{_pea};
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


sub chr_slice {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_chr_slice} = $val;
  }

  return $self->{_chr_slice};
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


sub DONOR_MATRIX {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_DONOR_MATRIX'} = $value;
  }
  
  if (exists($self->{'_CONFIG_DONOR_MATRIX'})) {
    return $self->{'_CONFIG_DONOR_MATRIX'};
  } else {
    return undef;
  }
}


sub ACCEPTOR_MATRIX {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_ACCEPTOR_MATRIX'} = $value;
  }
  
  if (exists($self->{'_CONFIG_ACCEPTOR_MATRIX'})) {
    return $self->{'_CONFIG_ACCEPTOR_MATRIX'};
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

sub ABINITIO_INTRONS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_ABINITIO_INTRONS'} = $value;
  }
  
  if (exists($self->{'_CONFIG_ABINITIO_INTRONS'})) {
    return $self->{'_CONFIG_ABINITIO_INTRONS'};
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


sub BEST_INTRONS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_BEST_INTRONS'} = $value;
  }
  
  if (exists($self->{'_CONFIG_BEST_INTRONS'})) {
    return $self->{'_CONFIG_BEST_INTRONS'};
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
1;
