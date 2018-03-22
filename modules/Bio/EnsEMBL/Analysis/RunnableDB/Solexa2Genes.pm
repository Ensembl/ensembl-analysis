=head1 LICENSE

# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::Solexa2Genes - 

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

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Solexa2Genes;

use warnings ;
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
  #$self->feature_slice_adaptor($self->get_dbadaptor($self->ALIGNMENT_DB)->get_SliceAdaptor);
  $self->repeat_slice_adaptor($self->db->get_SliceAdaptor);

  my $id = $self->input_id;
  my $slice = $self->fetch_sequence($id); 
  my $chr_slice;
  #hack to use more than 1 input db for alignments 
  my @features; 
  $chr_slice = $self->repeat_slice_adaptor->fetch_by_region('toplevel',
						       $slice->seq_region_name,
						      );
    
    $self->chr_slice($chr_slice);
  my $repeat_slice = $self->repeat_slice_adaptor->fetch_by_region
    ('toplevel',
     $slice->seq_region_name,
     $slice->start,
     $slice->end,
     1
    );
  my @repeats = sort { $a->start <=> $b->start } @{$self->repeat_feature_adaptor->fetch_all_by_Slice($repeat_slice,$self->REPEAT_LN)} ;
  # put on chromosome coords
  foreach my $repeat ( @repeats ) {
    $repeat = $repeat->transfer($chr_slice);
  }
  $self->repeats($self->make_repeat_blocks(\@repeats));
  my $feature_slice_adaptor = $self->get_dbadaptor($self->ALIGNMENT_DB)->get_SliceAdaptor;
 
  my $feature_slice =  $feature_slice_adaptor->fetch_by_region
    ( 
     'toplevel',
     $slice->seq_region_name,
     $slice->start,
     $slice->end,
    );
  $self->slice($feature_slice);
  $chr_slice = $feature_slice_adaptor->fetch_by_region('toplevel',
						       $slice->seq_region_name,
						      );
  
  $self->chr_slice($chr_slice);
  
  # do them in batches and return the hash of dafs
  if ( $self->BATCH_FETCH ) {
    $self->batch_num(0);
    $self->batch_fetch;
    return;
  }
  
  if (scalar(@{$self->LOGIC_NAMES}) < 1) {
    # no logic names have been specified in config file
    # so just load all DAFs
    push @features , @{$feature_slice->get_all_DnaAlignFeatures};
    print STDERR "Fetched ".scalar(@features)." \n";
  } else {
    # config file specifies at least one logic name
    foreach my $logic_name (@{$self->LOGIC_NAMES}){
      my @tmp_features = @{$feature_slice->get_all_DnaAlignFeatures($logic_name)};
      print STDERR "Fetched ".scalar(@tmp_features)." for $logic_name\n";
      push @features, @tmp_features;
    }
  }
  my %reads;
  foreach my $read ( @features ) {
    $read = $read->transfer($chr_slice);
  }

  # store the reads internally
  $self->reads(\@features);
}

sub batch_fetch {
  my ($self) = @_;
  my $slice = $self->slice;
  my $chr = $self->chr_slice;
  my $daf_adaptor = $self->get_dbadaptor($self->ALIGNMENT_DB)->get_DnaAlignFeatureAdaptor;
  # help getting too many open database connections?
  $self->get_dbadaptor($self->ALIGNMENT_DB)->dbc->disconnect_when_inactive(1);
  my $start = $self->batch_start;
  my $end  = $self->batch_end;
  my @dafs;
  unless ( $start ) {
    # first figure out the range of features we are talking about
    my $seq_id =  $self->sql( "SELECT seq_region_id from seq_region WHERE name = '" .  $slice->seq_region_name . "'",$self->db)->[0] ;
    # get the range of features to fetch 
    $end =  $self->sql( "  select max( dna_align_feature_id) from dna_align_feature where seq_region_id = $seq_id and seq_region_start >= " .$slice->start . " 
and seq_region_start < " . $slice->end  ,$self->get_dbadaptor($self->ALIGNMENT_DB))->[0] ;
     $start =  $self->sql( "  select min( dna_align_feature_id) from dna_align_feature where seq_region_id = $seq_id and seq_region_start >= " .$slice->start . " 
and seq_region_start < " . $slice->end  ,$self->get_dbadaptor($self->ALIGNMENT_DB))->[0] ;
    $self->batch_start($start);
    $self->batch_end($end);
  }
  my $batch_start = $self->batch_start;
  my $batch_end   = $batch_start + ( $self->BATCH_FETCH -1 ) ;
  # next time start at end +1
  $self->batch_start($batch_end+1);
  @dafs = @{$daf_adaptor->fetch_all_by_Slice_constraint($slice," dna_align_feature_id >= $batch_start AND dna_align_feature_id <= $batch_end ")};
  $self->read_count(scalar(@dafs));
  foreach my $read ( @dafs ) {
    $read = $read->transfer($chr);
  }


  $self->reads(\@dafs);
  $self->batch_num($self->batch_num+1);
  return;
}

sub sql {
  my ($self,$query,$db) = @_;
  print STDERR "query $query\n";
  my $sth = $db->dbc->prepare($query);
  $sth->execute();
  my @array = $sth->fetchrow_array;
  $sth->finish();
  return \@array;
}

sub run { 
  my ($self) =@_;
  return unless $self->reads;
  my @genes;
  my %exon_clusters;
  $self->cluster_count(0);
  if ( $self->BATCH_FETCH ) {
    my $gene_clusters;
    while ($self->read_count > 0 ) {
      $gene_clusters = $self->gene_cluster($gene_clusters);
      print STDERR "Found " . scalar(@$gene_clusters) . " clusters\n";
      # take one cluster of reads at a time
      foreach my $gene_cluster ( @$gene_clusters ) {
	next unless  scalar( @{$gene_cluster->{'features'}} );
	print $gene_cluster->{'start'} . " " .
	  $gene_cluster->{'end'} . " ";
	print scalar( @{$gene_cluster->{'features'}} ) ."\n";
	# break the gene clusters down into exon clusters 
	# linked by the pair information
	my %ec  = %{$self->exon_cluster($gene_cluster)};
	foreach my $key ( keys %ec ) {
	  $exon_clusters{$key} =  $ec{$key};
	}
      }
      # if we are batch fetching we can get overlapping clusters due to difficulty
      # in sorting reads when they are coming from ungapped features
      # you can only order them by the first feature, you dont know
      # what order the 2nd features will come in - anyway we need to merge the 
      # clusters to be sure of getting the same results as you would get if
      # you didnt batch fetch them
      
      my @tmp_clusters;
      my %list;
      foreach my $key (  keys %exon_clusters ) {
	# have to catch the merged exons or we merge them again
	next if $list{$exon_clusters{$key}->hseqname};
	$list{$exon_clusters{$key}->hseqname} =1;
	push @tmp_clusters,$exon_clusters{$key};
      }
      @tmp_clusters =  sort { $a->start <=> $b->start }  @tmp_clusters;
      
      #########################
      # exon clusters is acting as a map
      # why not use the map and a second 
      # that actually have them in the right order


      for ( my $i = 1 ; $i < scalar(@tmp_clusters ) ; $i++ ) {
	my $left_cluster = $tmp_clusters[$i-1];
	my $right_cluster = $tmp_clusters[$i];
	if ($left_cluster->end >= $right_cluster->start && 
	    $right_cluster->end >= $left_cluster->start ) {

	  # merge them
	  $left_cluster->start($right_cluster->start)
	    if $right_cluster->start < $left_cluster->start;
	  $left_cluster->end($right_cluster->end)
	    if $right_cluster->end > $left_cluster->end;
	  # merge scores
	  my $score =  $left_cluster->score + $right_cluster->score ;
	  $left_cluster->score( $score ) ;

	  # add any keys from right cluster into left cluster
	  foreach my $key ( @{$right_cluster->{'_keys'}} ) {
	    push @{$left_cluster->{'_keys'}},$key;
	  }
	  # need to update any other keys that may have been
	  # associated with this cluster in previous loops
	  foreach my $key ( @{$left_cluster->{'_keys'}} ) {
	    $exon_clusters{$key} = $left_cluster;
	  }
	  splice( @tmp_clusters,$i,1);
	  $i--;
	}
	
      }
      $self->batch_fetch;
    }
  } else {
    my @reads = @{$self->reads};
    # do it all in one go
    print STDERR "Found " . scalar(@reads) . " reads in total\nRunning initial clustering...";
    # Do an intial clustering to group together the read pairs
    # should roughly correspond to genes
    my $gene_clusters = $self->gene_cluster;
    print STDERR "Found " . scalar(@$gene_clusters) . " clusters\n";
    # take one cluster of reads at a time
    foreach my $gene_cluster ( @$gene_clusters ) {
      # break the gene clusters down into exon clusters 
      # linked by the pair information
      my %tmp_exon_clusters  = %{$self->exon_cluster($gene_cluster)};
      foreach my $key ( keys %tmp_exon_clusters ) {
	$exon_clusters{$key} =  $tmp_exon_clusters{$key};
      }
    }
  }
  my $transcripts = $self->process_exon_clusters(\%exon_clusters);
  if ( $transcripts ) {
    # Now we have collapsed our reads we need to make sure we keep the connections between
    # them so we can make our fake transcripts
    print STDERR "Found " . scalar(@$transcripts) . " transcripts \n";
    foreach my $transcript ( @$transcripts ) {
      #print  STDERR scalar(@$exon_cluster) ." exon clusters\n";
      next unless scalar(@$transcript) > 0;
      # make the dna align feature
      my $padded_exons = $self->pad_exons($transcript);
      push @genes , $self->make_gene($padded_exons) if $padded_exons ;
    }
  } else {
    print STDERR "No transcripts found for this slice\n";
  }

  $self->output(\@genes);
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
    if (  $tran->length < $self->MIN_LENGTH ){
      print "Rejecting because of length " . $tran->length ."\n";
      next;
    }
    if ( scalar(@{$gene->get_all_Exons}) == 1){
      if ( $tran->length <  $self->MIN_SINGLE_EXON_LENGTH ){
	print "Rejecting single exon transcript because of length " . $tran->length ."\n";
	next;
      }
    } else {
      # filter span on multiexon genes
      if( ( $tran->end - $tran->start +1 ) / $tran->length < ($self->MIN_SPAN) ) {
	if ( $tran->length <  $self->MIN_SINGLE_EXON_LENGTH ){
	  print "Rejecting because of span " . ( $tran->end - $tran->start +1 ) / $tran->length ."\n";
	  next;
	}
      }
    }
    eval {
      $gene_adaptor->store($gene);
    };    
    if ($@){
      warning("Unable to store gene!!\n$@");
      print STDERR $gene->start ." " . $gene->end ."\n";
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

sub process_exon_clusters {
  my ( $self, $exon_clusters ) = @_;
  my $cluster_data = $self->cluster_data;
  my $clean_exon_clusters;
  my @final_exon_clusters;
  my @transcripts;
  my $cluster_hash;
  my $pairs;
  # need to check that we have all the reads in before continuing
  # split into processing and clustering
  # dont use any reads in the processing only clusters and read names


  print STDERR "Processing Clusters\n";

  if ( scalar(keys %{$exon_clusters}) == 1 ) {
    # just keep single exon clusters for now - might be useful later
    foreach my $cluster ( values %{$exon_clusters} ) {
      my @transcript;
      push @transcript, $cluster;
      push @transcripts, \@transcript;
    }
    return \@transcripts;
  }

  # make the exon pairings store them in cluster hash
  foreach my $read ( keys %{$cluster_data} ) {
    my @clusters = keys %{$cluster_data->{$read}};
    my $left_cluster = $exon_clusters->{$clusters[0]}->hseqname;
    for (my $i = 1; $i < @clusters ; $i ++ ) {
      # exon clusters always points at the right cluster
      # even if the original cluster has been merge
      # ignore within cluster pairings which 
      # can happen if 2 clusters that started off as seperate
      # later got merged into one cluster
      my $right_cluster =  $exon_clusters->{$clusters[$i]}->hseqname;
      $cluster_hash->{$left_cluster}->{$right_cluster} ++
	unless $left_cluster eq $right_cluster;    }
  }

  

  my %tidy_clusters;
  # now lets clean up the hash so they key is the same
  # as the name of the cluster for merged clusters
  foreach my $cluster_name ( keys %{$exon_clusters} ) {
    my $ec = $exon_clusters->{$cluster_name};
    $tidy_clusters{$ec->hseqname} = $ec;
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
      print "$average_read_depth $pairs READ DEPTH $key $key2 $read_count\n";
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
    print "Ignoring $key\n";
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
    # get a non redundant set of exons
    foreach my $exon ( @$exon_cluster  ) {
      print "Adding exon $exon \n";
      push @transcript, $tidy_clusters{$exon};
    }
    @transcript =   sort { $a->start <=> $b->start} @transcript;
    push @transcripts, \@transcript;
  }
  return \@transcripts;
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
  my ($self,$gene_cluster ) = @_;
  my @ugfs;
  my $cluster_data = $self->cluster_data;
  # get the ungapped features ( corresponds to individual reads )
  return unless (scalar@{$gene_cluster->{'features'}});
  print STDERR "EXON CLUSTER \n";
  print STDERR scalar(@{$gene_cluster->{'features'}}) . " gapped reads\n";
  while ( scalar ( @{$gene_cluster->{'features'}} ) > 0 ){
    # drins the reads out of the clusters
    # batch fetch adds more reads to the clusters where appropriate
    my $read = pop(@{$gene_cluster->{'features'}});
    my @ungapped_features = sort { $a->start <=> $b->start } $read->ungapped_features;
    if ( scalar(@ungapped_features) > 2 ) {
      warn ( $read->hseqname ." is odd has more than 2 un gapped features\n");
      next;
    }
    if ( scalar(@ungapped_features) == 1 ) {
      my $single = $ungapped_features[0];
      push @ugfs, $single;
      next;
    }
    push @ugfs,  @ungapped_features ;
  }
  print STDERR scalar(@ugfs) . " ungapped reads\n";

  @ugfs = sort { $a->start <=> $b->start } @ugfs;
  my $feature_count = scalar(@ugfs);
  print STDERR "Clustering\n";
  my %exon_clusters;

  # make exon clusters and store the names of the reads and associated cluster number
  # need to find a way to merge the clusters and update
  # the cluster data accordingly

  READ: foreach my $ugf ( @ugfs ) {
    #print "UGF " . $ugf->start . " " , $ugf->end ."\n";
    my $clustered = 0;
    foreach my $exon_cluster ( values %exon_clusters ) {
      # do they overlap?
      if ( $ugf->start <= $exon_cluster->end+1 &&  $ugf->end >= $exon_cluster->start-1 ) {
	# Expand the exon_cluster
	$exon_cluster->start($ugf->start) if $ugf->start < $exon_cluster->start;
	$exon_cluster->end($ugf->end)     if $ugf->end   > $exon_cluster->end;    
	$exon_cluster->score($exon_cluster->score + 1);
	$clustered = 1;
	$cluster_data->{$ugf->hseqname}->{$exon_cluster->hseqname} = 1;
	# only allow it to be a part of a single cluster
	next READ;
      }
    }
    # start a new cluster if there is no overlap
    unless ( $clustered  ) {
      $ugf->score(1);
      $self->cluster_count($self->cluster_count +1 );
      # make a feature representing the cluster
      my $feat = Bio::EnsEMBL::FeaturePair->new
	(
	 -start      => $ugf->start,
	 -end        => $ugf->end,
	 -strand     => $ugf->strand,
	 -slice      => $ugf->slice,
	 -hstart     => $ugf->hstart,
	 -hend       => $ugf->hend,
	 -hstrand    => $ugf->hstrand,
	 -score      => $ugf->score,
	 -percent_id => $ugf->percent_id,
	 -hseqname   => "Cluster ". $self->cluster_count,
	 -analysis   => $ugf->analysis,
	);
      # store the clusters in a hash with a unique identifier
      $exon_clusters{$feat->hseqname} = $feat;
      # store the key within the feature
      push@{$feat->{'_keys'}},$feat->hseqname;
      $cluster_data->{$ugf->hseqname}->{$feat->hseqname} = 1;
    }
  }
  # store the relationships between the clusters
  $self->cluster_data($cluster_data);
  return  \%exon_clusters;
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
      $exon->start( $exon->start - ($trim)+1 );
    }
  }
  return \@padded_exons
}

=head2 gene_cluster
    Title        :   simple_cluster
    Usage        :   $self->($reads)
    Returns      :   Hash ref of clustered reads
    Args         :   None
    Description  :   Does a simple clustering of read pairs to identify 
                 :   groups of pairs that correspond to genes
=cut

sub gene_cluster {
  my ($self,$last_clusters) = @_;
  my $reads = $self->reads;
  my @clusters;
  # lets you continue if you are batch fetching reads
  @clusters = @$last_clusters if $last_clusters;
  my @features = sort { $a->start <=> $b->start }  @$reads ;
  print "FEATURES " .scalar(@features) ."\n";
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


sub make_repeat_blocks {
  my ($self,$repeats_ref) = @_;
  my %repeat_count;
  my @repeats = sort { $a->start <=> $b->start } @$repeats_ref;
  foreach my $rep ( @repeats ) {
    $repeat_count{$rep->analysis->logic_name}++;
  }
  foreach my $ln ( keys %repeat_count ) {
    print "Got " . $repeat_count{$ln} . " repeats of type $ln \n";
  }

  # merge repeat blocks
  for ( my $j = 1 ; $j <= $#repeats ; $j ++ ) { 
    if ( $repeats[$j]->start <= $repeats[$j-1]->end+1 ){
     # print "merging repeat $j " . $repeats[$j]->start . "-"  . $repeats[$j]->end. " " . $repeats[$j-1]->start ."-" . $repeats[$j-1]->end ."\n";
      $repeats[$j-1]->end($repeats[$j]->end) if  $repeats[$j]->end > $repeats[$j-1]->end ;
      splice(@repeats,$j,1);
      $j--;
    }
   # print "REPEAT $j " . $repeats[$j]->start . "-"  . $repeats[$j]->end. " " . $repeats[$j]->display_id ."\n";
  }  
  print STDERR " got " . scalar(@repeats) . " repeat blocks after merging\n";
  return \@repeats;
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

sub read_count {
  my ($self, $value) = @_;
  if (defined $value ) {
     $self->{'_read_count'} = $value;
  }
  return $self->{'_read_count'};
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

sub repeats {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_repeats} = $val;
  }

  return $self->{_repeats};
}

sub cluster_count {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_cluster_count} = $val;
  }

  return $self->{_cluster_count};
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

sub slice {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_slice} = $val;
  }

  return $self->{_slice};
}

sub batch_num {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_batch_num} = $val;
  }

  return $self->{_batch_num};
}

sub batch_start {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_batch_start} = $val;
  }

  return $self->{_batch_start};
}

sub batch_end {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_batch_end} = $val;
  }

  return $self->{_batch_end};
}

sub last_read {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_last_read} = $val;
  }

  return $self->{_last_read};
}

sub cluster_data {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_cluster_data} = $val;
  }

  return $self->{_cluster_data};
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

sub  MIN_SINGLE_EXON_LENGTH{
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_MIN_SINGLE_EXON_LENGTH'} = $value;
  }
  
  if (exists($self->{'_CONFIG_MIN_SINGLE_EXON_LENGTH'})) {
    return $self->{'_CONFIG_MIN_SINGLE_EXON_LENGTH'};
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

sub LOGIC_NAMES {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_LOGIC_NAMES'} = $value;
  }

  if (exists($self->{'_CONFIG_LOGIC_NAMES'})) {
    return $self->{'_CONFIG_LOGIC_NAMES'};
  } else {
    return [];
  }
}

sub BATCH_FETCH {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_BATCH_FETCH'} = $value;
  }
  
  if (exists($self->{'_CONFIG_BATCH_FETCH'})) {
    return $self->{'_CONFIG_BATCH_FETCH'};
  } else {
    return 0;
  }
}

sub REPEAT_LN {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_REPEAT_LN'} = $value;
  }
  
  if (exists($self->{'_CONFIG_REPEAT_LN'})) {
    return $self->{'_CONFIG_REPEAT_LN'};
  } else {
    return 0;
  }
}

1;
