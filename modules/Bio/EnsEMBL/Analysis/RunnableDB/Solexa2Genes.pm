# Cared for by Ensembl
#
# Copyright GRL & EBI
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Solexa2Genes

=head1 SYNOPSIS

  my $db      = Bio::EnsEMBL::DBAdaptor->new($locator);
  my $refine_genes = Bio::EnsEMBL::Analysis::RunnableDB::Solexa2Genes->new (
                                                    -db      => $db,
                                                    -input_id   => $input_id
                                                    -analysis   => $analysis );
  $refine_genes->fetch_input();
  $refine_genes->run();
  $refine_genes->write_output(); #writes to DB


=head1 DESCRIPTION


The databases containing the various features to combine is defined in
Bio::EnsEMBL::Analysis::Config::Databases and the configuration for the 
module is defined in Bio::EnsEMBL::Analysis::Config::GeneBuild::Solexa2Genes

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Solexa2Genes;

use strict;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;

use Bio::EnsEMBL::Analysis::Config::GeneBuild::Solexa2Genes;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils ;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils ;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);


sub new{
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);
  $self->read_and_check_config($SOLEXA2GENES_CONFIG_BY_LOGIC);
  return $self;
}

=head2 fetch_input

    Title        :   fetch_input
    Usage        :   $self->fetch_input
    Returns      :   nothing
    Args         :   none

=cut

sub fetch_input {
  my ($self) = @_;

  # store some adaptors
  $self->repeat_feature_adaptor($self->db->get_RepeatFeatureAdaptor);
  $self->feature_slice_adaptor($self->get_dbadaptor($self->ALIGNMENT_DB)->get_SliceAdaptor);
  $self->repeat_slice_adaptor($self->db->get_SliceAdaptor);
  

  my $id = $self->input_id;
  my $slice = $self->fetch_sequence($id);
  my $feature_slice =  $self->feature_slice_adaptor->fetch_by_region
    ( 
     'toplevel',
     $slice->seq_region_name,
     $slice->start,
     $slice->end,
    );
  my $chr_slice = $self->feature_slice_adaptor->fetch_by_region('toplevel',
								$slice->seq_region_name,
							       );
  $self->chr_slice($chr_slice);
  my @features = @{$feature_slice->get_all_DnaAlignFeatures};
  my %reads;
   
  while ( scalar(@features) > 0 )  {
    my $read = pop(@features);
    $reads{$read->hseqname} = $read->transfer($chr_slice);
  }
  # store the reads internally
  $self->reads(\%reads);
}



sub run { 
  my ($self) =@_;
  my %reads = %{$self->reads};
  my @genes;
  print STDERR "Found " . scalar(keys %reads) . " reads in total\nRunning initial clustering...";
  # Do an intial clustering to group together the read pairs
  # should roughly correspond to genes
  my $gene_clusters = $self->gene_cluster;
  print STDERR "Found " . scalar(@$gene_clusters) . " clusters\n";
  # take one cluster of reads at a time
  foreach my $gene_cluster ( @$gene_clusters ) {
    # break the gene clusters down into exon clusters 
    # linked by the pair information
    my ($exon_clusters)  = $self->exon_cluster($gene_cluster);
    next unless ( $exon_clusters );

    # Now we have collapsed our reads we need to make sure we keep the connections between
    # them so we can make our fake transcripts
    print STDERR "Found " . scalar(@$exon_clusters) . " transcripts \n";
    foreach my $exon_cluster ( @$exon_clusters ) {
      #print  STDERR scalar(@$exon_cluster) ." exon clusters\n";
      next unless scalar(@$exon_cluster) > 0;
      # make the dna align feature
      my $padded_exons = $self->pad_exons($exon_cluster);
      push @genes , $self->make_gene($padded_exons) if $padded_exons ;
    }
  }
  # merge genes separated by long repeats
  if ( $self->MERGE_GENES ) {
    my $merged_genes = $self->merge_genes(\@genes);
    $self->output($merged_genes);
  } else {
    $self->output(\@genes);
  }
}


sub write_output{
  my ($self) = @_;
  my @genes = @{$self->output};

  my $outdb = $self->get_dbadaptor($self->OUTPUT_DB);
  my $gene_adaptor = $outdb->get_GeneAdaptor;
 
 my $fails = 0;
  my $total = 0;
  
  foreach my $gene ( @genes ) {
    $gene->analysis($self->analysis);
    $gene->source($self->analysis->logic_name);
    $gene->biotype('rough');
    my $tran = $gene->get_all_Transcripts->[0];
    print "FILTERING " . $tran->start ." " , $tran->end ." ";
    # Filter models before writing them
    if ( scalar(@{$tran->get_all_Exons}) < $self->MIN_EXONS ) {
      print "Rejecting because of exon count " .  scalar(@{$tran->get_all_Exons}) ."\n";
      next;
    }
    if( ( $tran->end - $tran->start ) / $tran->length <= $self->MIN_SPAN ) {
      print "Rejecting because of span " . ( $tran->end - $tran->start ) / $tran->length ."\n";
      next;
    }
    if (  $tran->length > $self->MIN_LENGTH ){
      print "Rejecting because of length " . $tran->length ."\n";
    }
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
  print STDERR "$total genes written after filtering\n";
}

=head2 exon_cluster
    Title        :   exon_cluster
    Usage        :   $self->exon_cluster($ugfs)
    Returns      :   Array ref Bio::EnsEMBL::Exon
    Args         :   Array ref Bio::EnsEMBL::Exon
    Description  :   clusters individual reads into blocks representing exons
                 :   uses pair information to link blocks into transcripts
                 :   filters out poorly supported blocks
=cut

sub exon_cluster {
  my ($self,$gene_cluster) = @_;
  my @gapped_reads =  @{$gene_cluster->{'features'}};
  my @ugfs;
  # get the ungapped features ( corresponds to individual reads
  foreach my $read ( @gapped_reads ) {
    foreach my $ugf ($read->ungapped_features) {
      push @ugfs, $ugf;
    }
  }
  my $reads = $self->reads;
  my @exon_clusters;
  my $exon_cluster_num = 0;
  my $clean_exon_clusters;
  my @final_exon_clusters;
  my @transcripts;
  my @ungapped_features = sort { $a->start <=> $b->start } @ugfs;
  my $cluster_data;
  my $cluster_hash;
  my $pairs;

  # make exon clusters and store the names of the reads and associated cluster number
  foreach my $ugf ( @ungapped_features ) {
    my $clustered = 0;
    foreach (my $i = 0 ; $i < scalar(@exon_clusters) ; $i++ )  {
      my $exon_cluster = $exon_clusters[$i];
      # do they overlap?
      if ( $ugf->start <= $exon_cluster->end+1 &&  $ugf->end >= $exon_cluster->start-1 ) {
        # Expand the exon_cluster
        $exon_cluster->start($ugf->start) if $ugf->start < $exon_cluster->start;
        $exon_cluster->end($ugf->end)     if $ugf->end   > $exon_cluster->end;    
	$exon_cluster->score($exon_cluster->score + 1);
        $clustered = 1;
        $cluster_data->{$ugf->hseqname}->{$i} = 1;
      }
    }
    # start a new cluster if there is no overlap
    unless ( $clustered  ) {
      $cluster_data->{$ugf->hseqname}->{$exon_cluster_num} = 1 ;
      $ugf->score(1);
      push @exon_clusters, $ugf;
      $exon_cluster_num++;
    }
  }

  if ( scalar(@exon_clusters) == 1 ) {
    # just keep single exon clusters for now - might be useful later
    my $exon =  $exon_clusters[0];
    push @transcripts, [$exon];
    return \@transcripts;
  }

  
  print STDERR scalar(@exon_clusters) . " Exons\n";

  # make the exon pairings store them in cluster hash
  foreach my $read ( keys %{$cluster_data} ) {
    my @clusters = keys %{$cluster_data->{$read}};
    for (my $i = 0; $i < @clusters ; $i ++ ) {
      my $cluster = $clusters[$i];
      # adjust for compressed reads
      my $read_depth = 1;
      $read_depth = $reads->{$read}->p_value if $reads->{$read}->p_value;
      $cluster_hash->{$clusters[0]}->{$cluster} += $read_depth
	unless $clusters[0] == $cluster;
    }
  }

  # join the exon_clusters so that terminal exons get penalised
  # but put in all the other exons

  my $clean_cluster_hash;
  my $average_read_depth = 0;
  foreach my $key ( sort keys %{$cluster_hash} ) {
    # which exon_clusters is it connected to ?
    # is it an end exon? ie only connects to exon_clusters on one side?
    foreach my $key2 (  keys %{$cluster_hash->{$key}} ){
      # how many reads join these clusters to the right and the left?
      # initialise read count first - stop error messages
      my $read_count = $cluster_hash->{$key}->{$key2} ;
      $clean_cluster_hash->{$key}->{'right'}+= $read_count ;
      $clean_cluster_hash->{$key2}->{'left'}+= $read_count ;
      $average_read_depth += $read_count;
      print "READ DEPTH $key $key2 $read_count\n";
      # dont count single coverage pairs when
      # working out averages
      $pairs++ if $read_count >1;
    }
  }
  # whats the average coverage of the exons?
  $average_read_depth /=   $pairs if $pairs;
  print STDERR "Average READ DEPTH $average_read_depth \n";

  # Filter out pairings with very litte support
  foreach my $key ( sort keys %{$cluster_hash} ) {
    foreach my $key2 (  keys %{$cluster_hash->{$key}} ){
      if ($cluster_hash->{$key}->{$key2} < ($average_read_depth/100)*$self->EXCLUDE_EXONS ) {
	print "Removing pair  $key $key2 with a read_depth of ".$cluster_hash->{$key}->{$key2}  ."\n";
      	$cluster_hash->{$key}->{$key2} = 0;
      }
    }
  }

  # now lets  cycle through our exon_clusters removing end exons 
  # with little support  all other exons go forward to the next stage
  # add in single exon genes might be useful later
  
  foreach my $key ( sort keys %{$clean_cluster_hash} ) {
    # keep them if they are supported on both sides 
    if ( $clean_cluster_hash->{$key}->{'left'} && $clean_cluster_hash->{$key}->{'right'} ) {
      $clean_exon_clusters->{$key} = 1;
      next;
    }
    # is it a well supported left end exon?
    if ( $clean_cluster_hash->{$key}->{'left'}  && 
	 $clean_cluster_hash->{$key}->{'left'} >= $self->END_EXON_COVERAGE) {
      $clean_exon_clusters->{$key} =  1 ;
      next;
    }
    # is it a well supported right end exon?
    if ( $clean_cluster_hash->{$key}->{'right'} && 
	 $clean_cluster_hash->{$key}->{'right'} >= $self->END_EXON_COVERAGE ) {
      $clean_exon_clusters->{$key} = 1 ;
      next;
    }
    # othwise ignore it
  }

  # now need to find little clusters sitting in introns that are not connected to the transcript
  # do a reclustering based on which exons are connected to each other
  
  my @clean_clusters =  keys %$clean_exon_clusters;
  return unless ( scalar(@clean_clusters) > 0) ;

  # put one exon into the first cluster
  my @temp;
  push @temp,  pop(@clean_clusters);
  push @final_exon_clusters, \@temp;
  my $trans_count = 0;
 LOOP:  while ( scalar(@clean_clusters) > 0 ) {
    my $clustered;
    my $final_exon_cluster = $final_exon_clusters[$trans_count];
    foreach my $cluster_num ( @{$final_exon_cluster} ) {
      $clustered = 0;
     # do ANY of our exons join to this exon?
     for ( my $i =0  ; $i <= $#clean_clusters; $i++ )  {
	my $index = $clean_clusters[$i];
	# is the current exon attached to any exon in our cluster?
	if ( $cluster_hash->{$index}->{$cluster_num} or $cluster_hash->{$cluster_num}->{$index}) {
	  push @{$final_exon_cluster}, $index;
	  # chop it out
	  splice(@clean_clusters,$i,1);
	  $i--;
	  $clustered = 1;
	}
      }
    }
    unless ($clustered) {
      next unless scalar(@clean_clusters) > 0;
      my @temp;
      push @temp,  pop(@clean_clusters);
      # start another cluster
      push @final_exon_clusters, \@temp;
      $trans_count++;
    }
  }

  # So far we have just dealt with array indecies
  # now store the actual features
  foreach my $exon_cluster ( @final_exon_clusters ) {
    my @transcript;
    foreach my $exon (  @$exon_cluster ) {
      push @transcript, $exon_clusters[$exon];
    }
    @transcript =   sort { $a->start <=> $b->start} @transcript;
    push @transcripts, \@transcript;
  }
  return \@transcripts;
}

=head2 pad_exons
    Title        :   pad_exons
    Usage        :   $self->($exon_clusters)
    Returns      :   Array ref of Bio::EnsEMBL::Exon 
    Args         :   Array ref of Bio::EnsEMBL::Exon 
    Description  :   Takes an array of Exons, pads them and builds a 
                 :   DnaAlignFeature from it that represents a transcript
=cut

sub pad_exons {
  my ($self,$exon_cluster_ref) = @_;

  my @padded_exons;
  my @exon_clusters = sort { $a->start <=> $b->start } @$exon_cluster_ref;
  # make a padded exon array
  foreach my $exon ( @exon_clusters ){
    my $padded_exon =  create_Exon
      (
       $exon->start - 20,
       $exon->end + 20 ,
       -1,
       -1,
       -1,
       $exon->analysis,
       undef,
       undef,
       $self->chr_slice,
      );
    # dont let it fall of the slice because of padding
    $padded_exon->start(1) if $padded_exon->start <= 0;
    $padded_exon->end($self->chr_slice->length - 1) 
      if $padded_exon->end >= $self->chr_slice->length;

    my $feat = new Bio::EnsEMBL::DnaDnaAlignFeature
      (-slice    => $exon->slice,
       -start    => $padded_exon->start,
       -end      => $padded_exon->end,
       -strand   => -1,
       -hseqname => $exon->display_id,
       -hstart   => 1,
       -hstrand  => 1,
       -hend     => $padded_exon->length,
       -analysis => $exon->analysis,
       -score    => $exon->score,
       -cigar_string => $padded_exon->length.'M');
    my @feats;
    push @feats,$feat;
    $padded_exon->add_supporting_features(@feats);
    push @padded_exons, $padded_exon;
  }

  # dont let adjacent exons overlap
  for ( my $i = 1 ; $i <= $#padded_exons; $i++ ) {
  my $exon =  $padded_exons[$i];
  my $last_exon = $padded_exons[$i-1];
  if ( $last_exon->end >= $exon->start ) {
      # trim the exons so they dont overlap
      my $trim = int(($exon->start - $last_exon->end) /2);
      $last_exon->end(   $last_exon->end  + ($trim) -1 );
      $exon->start( $exon->start - ($trim) );
    }
  }
  return \@padded_exons
}

=head2 simple_cluster
    Title        :   simple_cluster
    Usage        :   $self->($reads)
    Returns      :   Hash ref of clustered reads
    Args         :   None
    Description  :   Does a simple clustering of read pairs to identify 
                 :   groups of pairs that correspond to genes
=cut

sub gene_cluster {
  my ($self) = @_;
  my $reads = $self->reads;
  my @clusters;
  my @features = sort { $a->start <=> $b->start }  values %$reads ;
  foreach my $f ( @features ) {
    my $clustered = 0;
    foreach my $cluster ( @clusters ) {
      # do they overlap?
      if ( $f->start <= $cluster->{'end'} 
	   &&  $f->end >= $cluster->{'start'} ) {
	# Expand the cluster
	$cluster->{'start'} = $f->start 
	  if $f->start < $cluster->{'start'};
	$cluster->{'end'} = $f->end
	  if $f->end   > $cluster->{'end'}; 
	push @{$cluster->{'features'}}, $f;
	$clustered = 1;
      }
    }
    unless ($clustered) {
      my $cluster;
      push @{$cluster->{'features'}}, $f;
      $cluster->{'start'} = $f->start;
      $cluster->{'end'}   = $f->end;
      push @clusters, $cluster;
    }
  }
  return \@clusters;
}

=head2 make_gene
    Title        :   make_gene
    Usage        :   $self->make_gene($exons)
    Returns      :   Array ref of Bio::EnsEMBL::Gene objects
    Args         :   Array ref of Bio::EnsEMBL::Exon objects
    Description  :   Builds gene models from an array of exons
=cut

sub make_gene {
  my ($self,$exon_ref) = @_;
  # we are making reverse strand genes so the exons need to be in the reverse order
  my @exons = sort { $b->start <=>  $a->start } @$exon_ref;
  my $tran =  new Bio::EnsEMBL::Transcript(-EXONS => \@exons);
  $tran->analysis($self->analysis);
 return @{convert_to_genes(($tran),$self->analysis)};
}

=head2 merge_genes
    Title        :   merge_genes
    Usage        :   $self->make_genes($genes)
    Returns      :   Array ref of Bio::EnsEMBL::Gene objects
    Args         :   Array ref of Bio::EnsEMBL::Gene objects
    Description  :   Merges gene models when they are separated
                 :   by a gap of <= MERGE_GENES and the gap is 100%
                 :   covered by repeats
=cut

sub merge_genes { 
  my ($self, $genes_ref) = @_;
  my @merged_genes;
  # join together neighbouring genes separated by a gap of X
  # where the gap is 100% filled by a repeat
  # get rid of nested genes first
  my @genes = sort {$a->start <=> $b->start} @$genes_ref;
  print " Got ". scalar(@genes) ."\n";
  for ( my $i = 1 ; $i <= $#genes ; $i++ ) {

    if ( $genes[$i-1]->start < $genes[$i]->start && 
	 $genes[$i-1]->end > $genes[$i]->end ) {
      push @merged_genes, splice (@genes,$i,1);
      $i--;
    }
  }

  # is the gap between the genes short enough?
 GENE:  for ( my $i = 1 ; $i <= $#genes ; $i++ ) {
    my $left_gene = $genes[$i-1];
    my $right_gene = $genes[$i];
    my $gap = $right_gene->start - $left_gene->end;
    if ($gap <= $self->MERGE_GENES && $left_gene->end < $right_gene->start) {
      # is it covered by repeats?
      my $repeat_slice = $self->repeat_slice_adaptor->fetch_by_region
	('toplevel',
	 $left_gene->slice->seq_region_name,
	 $left_gene->end,
	 $right_gene->start,
	);
      #print "looking at " . $repeat_slice->name ."\n";
      # is the gap covered by a repeat?
      my @repeats = sort { $a->start <=> $b->start } @{$self->repeat_feature_adaptor->fetch_all_by_Slice($repeat_slice)} ;
      # merge repeat blocks
      #print scalar(@repeats) . " repeats found \n";
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
      # merge the genes together where repeat coverage is 100%
      if ($repeat_coverage >= 0.95   ) {
	print "Do some merging\n";
	print "Repeat Coverage = $repeat_coverage \n";
	print $repeat_slice->name ."\n";
	# to do the merge do we just link the genes or do we join them with a long exon?
	# how about we keep it as 2 exons but bring them together so they abut each other
	# if there is evidene of splicing they will get separated in the refine genes code
	# otherwise they will stay merged, given that we are looking at areas covered by repeat
	# it seems likely that the exons should read right through.
	# so we want to extend the last exon of the left gene to finish at the end
	# of the first exon of the right gene, we also want to lose the first right
	# gene exon
	my @merged_exons;
	my @left_exons = sort { $a->start <=> $b->start } @{$left_gene->get_all_Exons};
	my @right_exons = sort { $a->start <=> $b->start } @{$right_gene->get_all_Exons};	
	my $merged_exon = pop(@left_exons);
	my $disgarded_exon = shift(@right_exons);
	$merged_exon->end($disgarded_exon->end);
	push @merged_exons,@left_exons;
	push @merged_exons,$merged_exon;
	push @merged_exons,@right_exons;
	my @new_genes = $self->make_gene(\@merged_exons);
	my $new_gene = $new_genes[0];
	print "NEW GENE " . $new_gene->start ." ". $new_gene->end ."\n";
	# take them out of the array and carry on
	splice (@genes,$i-1,2,$new_gene);
	$i-= 1;
	print "Merged 1 successfully\n";
      }
    }
  }
  # add whatever is left to the return array
  push @merged_genes, @genes;
  print "Returning " . scalar (@merged_genes) . "\n";
  return \@merged_genes;
}

###########################################
# Containers


sub reads {
  my ($self, $value) = @_;
  if ($value ) {
     $self->{'_reads'} = $value;
  }
  return $self->{'_reads'};
}

sub feature_slice_adaptor {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_fsa} = $val;
  }

  return $self->{_fsa};
}

sub repeat_slice_adaptor {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_rsa} = $val;
  }

  return $self->{_rsa};
}

sub repeat_feature_adaptor {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_rfa} = $val;
  }

  return $self->{_rfa};
}

sub chr_slice {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_chr_slice} = $val;
  }

  return $self->{_chr_slice};
}
##########################################
# Config variables 

sub ALIGNMENT_DB {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_ALIGNMENT_DB'} = $value;
  }
  
  if (exists($self->{'_CONFIG_ALIGNMENT_DB'})) {
    return $self->{'_CONFIG_ALIGNMENT_DB'};
  } else {
    return undef;
  }
}

sub  MIN_LENGTH{
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MIN_LENGTH'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MIN_LENGTH'})) {
    return $self->{'_CONFIG_MIN_LENGTH'};
  } else {
    return 0;
  }
}

sub MIN_EXONS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MIN_EXONS'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MIN_EXONS'})) {
    return $self->{'_CONFIG_MIN_EXONS'};
  } else {
    return 0;
  }
}

sub MIN_SPAN {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MIN_SPAN'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MIN_SPAN'})) {
    return $self->{'_CONFIG_MIN_SPAN'};
  } else {
    return 0;
  }
}

sub END_EXON_COVERAGE {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_END_EXON_COVERAGE'} = $value;
  }
  
  if (exists($self->{'_CONFIG_END_EXON_COVERAGE'})) {
    return $self->{'_CONFIG_END_EXON_COVERAGE'};
  } else {
    return 0;
  }
}


sub EXCLUDE_EXONS {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_EXCLUDE_EXONS'} = $value;
  }
  
  if (exists($self->{'_CONFIG_EXCLUDE_EXONS'})) {
    return $self->{'_CONFIG_EXCLUDE_EXONS'};
  } else {
    return 0;
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


sub MERGE_GENES {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MERGE_GENES'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MERGE_GENES'})) {
    return $self->{'_CONFIG_MERGE_GENES'};
  } else {
    return 0;
  }
}



1;
