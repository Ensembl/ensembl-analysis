=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::RunnableDB::Solexa2GenesLiteNew - 

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

package Bio::EnsEMBL::Analysis::RunnableDB::Solexa2GenesLiteNew;

use warnings ;
use strict;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;

use Bio::EnsEMBL::Analysis::Config::GeneBuild::Solexa2GenesLiteNew;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::TranscriptUtils ;
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils ;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Time::HiRes qw(gettimeofday);  
#use Time::TimeTick qw(timetick) ; 
use vars qw(@ISA);

@ISA = qw (Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);


sub new{
  my ($class,@args) = @_; 
  my $self = $class->SUPER::new(@args);    

  my ($ignore_config_file) = rearrange (['IGNORE_CONFIG_FILE'], @args);

  if ( $ignore_config_file != 1 ) {  
    $self->read_and_check_config($SOLEXA2GENES_CONFIG_BY_LOGIC); 
  } 
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
  #timetick("fetch_input start"); 

  # store some adaptors

  if ( ref($self->ALIGNMENT_DB =~m/ARRAY/ )) {  
    if ( scalar(@{$self->ALIGNMENT_DB})>1)   {
      throw(" multiple alignment db support not implemented yet\n") ;
    } 
  } 
    
  #print "using output db " . $self->OUTPUT_DB . "\n";

  #print "\n\nthis is a bit of a hacky script - we assume that you only have data for one species in your dna align feature table \n"; 
  #$self->feature_slice_adaptor($self->get_dbadaptor(${$self->ALIGNMENT_DB}[0])->get_SliceAdaptor);
  #$self->repeat_slice_adaptor($self->db->get_SliceAdaptor);

  my $feature_slice_adaptor = $self->db->get_SliceAdaptor();
  $self->db->disconnect_when_inactive(1);
  $self->feature_slice_adaptor($feature_slice_adaptor); 
#  my $feature_slice_adaptor;
#  my %dbnames_2_logicnames = %{ $self->DNA_ALIGN_FEATURE_DATA} ; 
#  foreach my $db_hash_key ( keys %dbnames_2_logicnames )  {
#     my $set_db = $self->get_dbadaptor($db_hash_key);
#     my $sa = $set_db->get_SliceAdaptor(); 
#     print "sa : $sa\n"; 
#     $feature_slice_adaptor = $sa; 
#   }
#  

  my $id = $self->input_id;
  my $slice = $self->fetch_sequence($id);
  my $feature_slice =  $self->feature_slice_adaptor->fetch_by_region
    ( 
     'toplevel',
     $slice->seq_region_name,
     $slice->start,
     $slice->end,
    );
  my $chr_slice = $self->feature_slice_adaptor->fetch_by_region('toplevel', $slice->seq_region_name );
  $self->chr_slice($chr_slice);  

  #print "fetching daf \n";
  #my $fetched_daf_features = $self->get_dna_align_features_from_databases($slice);   
# begin simons old code 
  my @reads;
  my $seq_id =  $self->sql( "SELECT seq_region_id from seq_region WHERE name = '" .  $slice->seq_region_name . "'",$self->db)->[0] ;   
   
  my @alignment_dbs; 
  if (! ref($self->ALIGNMENT_DB)=~m/ARRAY/){ 
     push @alignment_dbs, $self->ALIGNMENT_DB;
  } 
  foreach my $DB ( @alignment_dbs){ 
    #print "processing $DB\n";
    # assumption here we have only  one lane per database; 
  my $handle =  $self->sql_array( "SELECT dna_align_feature_id, seq_region_id, seq_region_start, seq_region_end, cigar_line from dna_align_feature where seq_region_id = $seq_id and seq_region_end >= " . $slice->start .  " and seq_region_start <= " .$slice->end . " ;" ,$self->get_dbadaptor($DB));

  while ( my @read = $handle->fetchrow_array ) {
    push @reads, \@read;
  }
 }
#  print "Got " . scalar(@reads) . " features \n";
  @reads = sort { $a->[2] <=> $b->[2] } @reads;
  # store the reads internally
  $self->reads(\@reads);


  #print "fetching repeat features \n" ; 
  my $rfa = $self->db->get_RepeatFeatureAdaptor();
  my @repeats = sort { $a->start <=> $b->start } @{$rfa->fetch_all_by_Slice($feature_slice)} ;
  foreach my $repeat ( @repeats ) {
    $repeat = $repeat->transfer($chr_slice);
  }
  $self->repeats($self->make_repeat_blocks(\@repeats));

  print scalar(@repeats) . " REPEATS fetched \n"; 
#
#  my $seq_id =  $self->sql( "SELECT seq_region_id from seq_region WHERE name = '" .  $slice->seq_region_name . "'",$self->db)->[0] ;
#
#  my @reads;
#  foreach my $DB ( @{$self->ALIGNMENT_DB} ) {  # change here ! 
#    print "processing $DB\n"; 
#  #my $handle =  $self->sql_array( "SELECT dna_align_feature_id, seq_region_id, seq_region_start, seq_region_end, cigar_line from dna_align_feature where seq_region_id = $seq_id and seq_region_end >= " . $slice->start .  " and seq_region_start <= " .$slice->end . " and analysis_id =  ;" ,$self->get_dbadaptor($self->ALIGNMENT_DB)); 
#  #
#  # assumption here that you only have one lane per database; 
#
# 
#  my $handle =  $self->sql_array( "SELECT dna_align_feature_id, seq_region_id, seq_region_start, seq_region_end, cigar_line from dna_align_feature where seq_region_id = $seq_id and seq_region_end >= " . $slice->start .  " and seq_region_start <= " .$slice->end . " ;" ,$self->get_dbadaptor($DB)); 
#
#  while ( my @read = $handle->fetchrow_array ) {
#    push @reads, \@read;
#  }
# }
#  print "Got " . scalar(@reads) . " features \n";
#  @reads = sort { $a->[2] <=> $b->[2] } @reads;
#  # store the reads internally
#  $self->reads(\@reads); 
}



sub run { 
  my ($self) =@_;
  my @genes;

  my $t1 = gettimeofday(); 

  # Do an intial clustering to group together the read pairs
  # should roughly correspond to genes
  my $gene_clusters = $self->gene_cluster;
  #print STDERR "Found " . scalar(@$gene_clusters) . " clusters\n";
  # take one cluster of reads at a time
  foreach my $gene_cluster ( @$gene_clusters ) {
    # break the gene clusters down into exon clusters 
    # linked by the pair information 
    #print STDERR "processing gene cluster\n"; 
    my ($exon_clusters)  = $self->exon_cluster($gene_cluster);
    next unless ( $exon_clusters );

    # Now we have collapsed our reads we need to make sure we keep the connections between
    # them so we can make our fake transcripts
    #print STDERR "Found " . scalar(@$exon_clusters) . " transcripts \n";
    foreach my $exon_cluster ( @$exon_clusters ) {
      #print  STDERR scalar(@$exon_cluster) ." exon clusters\n";
      next unless scalar(@$exon_cluster) > 0;
      # make the dna align feature
      my $padded_exons = $self->pad_exons($exon_cluster); 
      #print STDERR "exons padded\n"; 
      push @genes , $self->make_gene($padded_exons) if $padded_exons ;
      #print STDERR "gene made\n"; 
    }
  }
  # merge genes separated by long repeats
  if ( $self->MERGE_GENES ) {
    my $merged_genes = $self->merge_genes(\@genes);
    $self->output($merged_genes);
  } else {
   $self->output(\@genes);
  } 
  #print "TIME : " . (gettimeofday() - $t1 ) . "\n";
}


sub write_output{
  my ($self) = @_;
  my @genes = @{$self->output};

  my $outdb = $self->get_dbadaptor($self->OUTPUT_DB);
  my $gene_adaptor = $outdb->get_GeneAdaptor;
 
  my $fails = 0;
  my $total = 0; 

  if ( $self->USE_ANALYSIS_LOGIC_NAME_AS_DEFAULT_GENE_OUTPUT_BIOTYPE == 1 ) {
    #print "using logic_name : " . $self->analysis->logic_name . " as biotype for genes and transcripts\n";
  } else { 
    #print "using generic biotype \"rough\" for genes and transcripts\n"; 
  } 
  foreach my $gene ( @genes ) {
    $gene->analysis($self->analysis);
    $gene->source($self->analysis->logic_name); 

    if ( $self->USE_ANALYSIS_LOGIC_NAME_AS_DEFAULT_GENE_OUTPUT_BIOTYPE == 1 ) {
      $gene->biotype($self->analysis->logic_name);
    }else {
      $gene->biotype('rough');
    }

    my $tran = $gene->get_all_Transcripts->[0];
    #print "FILTERING " . $tran->start ." " , $tran->end ." ";
    $tran->biotype($self->analysis->logic_name); # this is important for step3 
    # Filter models before writing them
    if ( scalar(@{$tran->get_all_Exons}) < $self->MIN_EXONS ) {
      #print "Rejecting because of exon count " .  scalar(@{$tran->get_all_Exons}) ."\n";
      next;
    }
    if( ( $tran->end - $tran->start ) / $tran->length <= $self->MIN_SPAN ) {
      #print "Rejecting because of span " . ( $tran->end - $tran->start ) / $tran->length ."\n";
      next;
    }
    if (  $tran->length < $self->MIN_LENGTH ){
      #print "Rejecting because of length " . $tran->length ."\n";
      next;
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
  #print STDERR "$total genes written after filtering\n";
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
  #print "got " .scalar(@gapped_reads) ." gapped reads\n";
  # get the ungapped features ( corresponds to individual reads
  foreach my $read ( @gapped_reads ) {
    foreach my $ugf (@{$self->ungapped_features($read)}) {
      push @ugfs, $ugf;
    }
  }
  #print "got " .scalar(@ugfs) ." ugfs\n";

  my @exon_clusters;
  my $exon_cluster_num = 0;
  my $clean_exon_clusters;
  my @final_exon_clusters;
  my @transcripts;
  my @ungapped_features = sort { $a->[2] <=> $b->[2] } @ugfs;
  my $cluster_data;
  my $cluster_hash;
  my $pairs;

  # make exon clusters and store the names of the reads and associated cluster number 
  my $jhv_cnt = 0;   #jhv
  my $last_cluster_id = 0; #jhv
  foreach my $ugf ( @ungapped_features ) { 
    $jhv_cnt++;  #jhv
    if ( ($jhv_cnt % 500) == 0 ) {   #jhv
      #print "$jhv_cnt ungapped feat. processed - " . scalar(@exon_clusters) . " exon clusters\n"; #jhv
    } #jhv 
    my $clustered = 0; 
    
#org foreach (my $i = 0 ; $i < scalar(@exon_clusters) ; $i++ )   
    foreach (my $i = $last_cluster_id ; $i < scalar(@exon_clusters) ; $i++ )  {
      my $exon_cluster = $exon_clusters[$i];  

      #print "ec: $jhv_cnt / ".scalar(@ungapped_features) . " $i / ".scalar(@exon_clusters) ." ".  $exon_cluster->start ." - " . $exon_cluster->end . " VS ";  
      #print " ".$ugf->[2] . " -- " . $ugf->[3]. "\n";
      # do they overlap? 
      # $ugf->[2]  ungapped feature start ? 
      # $ugf->[3]  ungapped feature end   ?  
      if ( $ugf->[2] <= $exon_cluster->end+1 &&  $ugf->[3] >= $exon_cluster->start-1 ) {
        # Expand the exon_cluster 
       # print " overlap \n";
        $exon_cluster->start($ugf->[2]) if $ugf->[2] < $exon_cluster->start;
        $exon_cluster->end($ugf->[3])     if $ugf->[3]   > $exon_cluster->end;    
	$exon_cluster->score($exon_cluster->score + 1);
        $clustered = 1;
        $cluster_data->{$ugf->[0]}->{$i} = 1;
      }#else{
       # print "no overlap\n";
      #}
      
      $last_cluster_id = $i ; 
    }

    # start a new cluster if there is no overlap
    unless ( $clustered  ) {
      $cluster_data->{$ugf->[0]}->{$exon_cluster_num} = 1 ;
      # make a daf maybe?
      my $feat =  Bio::EnsEMBL::FeaturePair->new
        (-SLICE      => $self->chr_slice,
         -START      => $ugf->[2],
         -END        => $ugf->[3],
         -STRAND     => 1,
         -HSLICE      => $self->chr_slice,
         -HSEQNAME   => $ugf->[0],
         -HSTART     => 1,
         -HEND       => $ugf->[3] - $ugf->[2] +1,
         -HSTRAND    => 1,
         -SCORE      => 1,
         -PERCENT_ID => 100,
         -ANALYSIS   => $self->analysis,
         -P_VALUE    => 1,
         -HCOVERAGE   => 100);

      push @exon_clusters, $feat;
      $exon_cluster_num++;
    }
  }
  #print "processed ungapped features \n";  
  # check they dont overlap
  @exon_clusters = sort { $a->start <=> $b->start } @exon_clusters; 

  for (my $i = 1; $i <= $#exon_clusters ; $i ++ ) {
    my $exon = $exon_clusters[$i];
    my $prev_exon = $exon_clusters[$i-1];
    if ( $prev_exon->end >= $exon->start ) {
      #print "OVERLAP at clustering" .  $prev_exon->end . " " .  $exon->start ."\n";
    }
  }


  if ( scalar(@exon_clusters) == 1 ) {
    # just keep single exon clusters for now - might be useful later
    my $exon =  $exon_clusters[0];
    push @transcripts, [$exon];
    return \@transcripts;
  }

  
  print scalar(@exon_clusters) . " Exon clusters\n";

  # make the exon pairings store them in cluster hash
  foreach my $read ( keys %{$cluster_data} ) {
    my @clusters = keys %{$cluster_data->{$read}};
    for (my $i = 0; $i < @clusters ; $i ++ ) {
      my $cluster = $clusters[$i];
      $cluster_hash->{$clusters[0]}->{$cluster}++
	unless $clusters[0] == $cluster;
    }
  }

  # join the exon_clusters so that terminal exons get penalised
  # but put in all the other exons

  my $clean_cluster_hash;
  my $average_read_depth = 0;
  #timetick("cluster-start"); 
  foreach my $key ( sort { $a <=> $b }  keys %{$cluster_hash} ) {
    # which exon_clusters is it connected to ?
    # is it an end exon? ie only connects to exon_clusters on one side?
#      my $e = $exon_clusters[$key];  # change jhv
    foreach my $key2 (  keys %{$cluster_hash->{$key}} ){
      # how many reads join these clusters to the right and the left?
      # initialise read count first - stop error messages
      my $read_count = $cluster_hash->{$key}->{$key2} ;
      $clean_cluster_hash->{$key}->{'right'}+= $read_count ;
      $clean_cluster_hash->{$key2}->{'left'}+= $read_count ;
      $average_read_depth += $read_count; 

      #my $e = $exon_clusters[$key]; # sw4 - jan: i moved this out of this loop as it's quicker 
#      
#      my $e2 = $exon_clusters[$key2];
#      print "EXON  $key " . $e->start . " " . $e->end . " rpb " ;
#      my $rpb1 =  $e->score / $e->length;
#      print "$rpb1 JOINS $key2 " . $e2->start . " " . $e2->end . " rpb " ;
#      my $rpb2 =  $e2->score / $e2->length; 
#
#      my $rpb_ratio;
#      if (  $rpb1 > $rpb2 ){
# 	$rpb_ratio =  $rpb1 / $rpb2;
#      }else {
#	$rpb_ratio =  $rpb2 / $rpb1;
#      }
#      my $separation;
#      if ( $e->start < $e2->start ) {
#	$separation =  $e2->start - $e->end ;
#      } else {
#	$separation =  $e->start - $e2->end ;
#      }
#      print "$rpb2 separated by $separation joined by $read_count  ratio $rpb_ratio "; 
#      my $depth_connectivity1 =  $read_count /  $rpb1; 
#      my $depth_connectivity2 =  $read_count /  $rpb2; 
#      print " $depth_connectivity1 $depth_connectivity2 \n";
#      # dont count single coverage pairs when
#      # working out averages
#      $pairs++ if $read_count >1; 
#      my $cutoff = 0;
#      $cutoff += $separation / 3000;
      # what is the cut off to use and how to bias it by length????
#      if (  $read_count /  $rpb1 < $cutoff  or  $read_count /  $rpb2 < $cutoff ) {
#	print "killing this pair with cutoff $cutoff \n";
#	$cluster_hash->{$key}->{$key2} = 0;
#     } 
    }
  }
  #timetick("cluster-end"); 
  # filtering want to compare the numbers of reads connecting exons
  # with the number of reads in the exon
  # if your exon has 30,000 reads covering it and only 1 conecting it
  # to something else it is probably a mistake - we expect around a third
  # of our reads to be intron reads usually ??
  # we want the exon coverage to be reads - per base so as to treat all exons equally
  


  # whats the average coverage of the exons?
  $average_read_depth /=   $pairs if $pairs;
#  print STDERR "Average READ DEPTH $average_read_depth \n";

  # Filter out pairings with very litte support
#  foreach my $key ( sort keys %{$cluster_hash} ) {
#    foreach my $key2 (  keys %{$cluster_hash->{$key}} ){
#      if ($cluster_hash->{$key}->{$key2} < ($average_read_depth/100)*$self->EXCLUDE_EXONS ) {
#	print "Removing pair  $key $key2 with a read_depth of ".$cluster_hash->{$key}->{$key2}  ."\n";
#      	$cluster_hash->{$key}->{$key2} = 0;
#      }
#    }
#  }

  #print "\ncycling trough exon_clusters " . scalar(keys %{$clean_cluster_hash}) . " clusters \n"; 
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
  #print "\ndone cycling\n"; 
  # now lets  cycle through our exon_clusters removing end exons 

  # now need to find little clusters sitting in introns that are not connected to the transcript
  # do a reclustering based on which exons are connected to each other
  
  my @clean_clusters =  keys %$clean_exon_clusters;
  return unless ( scalar(@clean_clusters) > 0) ;

  # put one exon into the first cluster
  my @temp;
  push @temp,  pop(@clean_clusters);
  push @final_exon_clusters, \@temp;
  my $trans_count = 0; 


 my $tx = gettimeofday(); 
# ORIGINAL CODE 
# LOOP:  while ( scalar(@clean_clusters) > 0 )  
# ORIGINAL CODE 
#    my $clustered;
#    my $final_exon_cluster = $final_exon_clusters[$trans_count];  
#
#    foreach my $cluster_num ( @{$final_exon_cluster} ) {
#      $clustered = 0;
#     # do ANY of our exons join to this exon?
#      for ( my $i =0  ; $i <= $#clean_clusters; $i++ )  { 
# 	my $index = $clean_clusters[$i];
# 	# is the current exon attached to any exon in our cluster?
# 	if ( $cluster_hash->{$index}->{$cluster_num} or $cluster_hash->{$cluster_num}->{$index}) {
# 	  push @{$final_exon_cluster}, $index;
# 	  # chop it out 
# 	  splice(@clean_clusters,$i,1);
# 	  $i--;
# 	  $clustered = 1;
#       }
#     }
#   } 

  # now need to find little clusters sitting in introns that are not connected to the transcript
  # do a reclustering based on which exons are connected to each other

#we have a list of clean clusters
#we make the first of the clean clusters a final_exon_cluster 
#  we loop trough the list of all clean clusters cc  
#    each clean cluster has an index, which is $clean_clusters[i] 
#    we check if this index exists in cluster_hash   
#    if 
#      it exists in the cluster_hash, we remove it from the list of clean_clusters and we mark 'clustered'  
#    else 
#       continue to loop 
#    end
#
#  we loop rough shorter list of clean clusters
#    we check if this index exists in cluster_hash   
# cluster hash is a hash of hashes 
# - idea is to identify elements which are not in $cluster_hash->{$index}->{$cluster_num} $cluster_hash->{$cluster_num}->{$index} 
#
#timetick("while-start"); 

LOOP:  while ( scalar(@clean_clusters) > 0 ) {  
    # final-exon_cluster will later contain exon clustes which we use to make transcripts 
    my $clustered;
    my $final_exon_cluster = $final_exon_clusters[$trans_count];  

#print "final_exon_cluster  [$trans_count] : " . join(" ", @$final_exon_cluster) . " (".scalar(@final_exon_clusters) . " elements in array)\n"; 
#print "trans-count : $trans_count\n"; 

    my $c1jhv=0;
    #print "STARTING outer loop - will loop over " . scalar( @{$final_exon_cluster} )  . " elements\n";  

    OUTER: foreach my $cluster_num ( @{$final_exon_cluster} ) {
      $c1jhv++;
      $clustered = 0; 
     # do ANY of our exons join to this exon?
      INNER: for ( my $i =0  ; $i <= $#clean_clusters; $i++ )  {  
  	my $index = $clean_clusters[$i];
        
        #print "OUTER: jhv cluster_nr: $c1jhv / " . scalar(@{$final_exon_cluster}) . "   " . $i . " / ".scalar(@clean_clusters) . "\t";
        #print "cluster_num $cluster_num -  index: $index   - i  = $i\n";  

 	# is the current exon attached to any exon in our cluster? 
 	if ( $cluster_hash->{$index}->{$cluster_num} or $cluster_hash->{$cluster_num}->{$index}) { 
          # extending the outer loop  
 	  push @{$final_exon_cluster}, $index; 
          #print "extension of outer loop - will loop over ". scalar(@{$final_exon_cluster}) . " elements now\n";
 	  # chop it out 
 	  #print "removing element $i from clean_clusters array \n";
 	  splice(@clean_clusters,$i,1);
 	  $i--;
 	  $clustered = 1;
 	}
       } 
       #print "next OUTER loop\n";
    } 
    #print "OUTER loops finished;\n";   
    # - If both loops ran through clean and clustered is set = 0 then we take the last elemenet of the the clean_clusters, 
    #   -> and we add it to final_exon_clusters; 
    #      we increment the trans_count   
    # if it's clustered, the trans_count stays the same.
    

    unless ($clustered) {     # if clustered == 0 
      next unless scalar(@clean_clusters) > 0; 
      # start another cluster
      push @final_exon_clusters,  [ pop(@clean_clusters)] ;  #jhv
      $trans_count++;
    }
  }   # end while LOOP 

#  timetick("while-end "); 

 #print "WHILE : " .  (gettimeofday() - $tx) . "\n"; 
  #print "\ndone while looping\n"; 

  # So far we have just dealt with array indecies
  # now store the actual features
  foreach my $exon_cluster ( @final_exon_clusters ) {
    my @transcript;
    foreach my $exon (  @$exon_cluster ) { 
      # print "exon $exon_clusters[$exon]\n";  # isa Bio::EnsEMBL::FeaturePair 
      push @transcript, $exon_clusters[$exon];
    }
    @transcript =   sort { $a->start <=> $b->start} @transcript;
    push @transcripts, \@transcript;
  } 

  #print "\nreturning transcripts\n"; 
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
  # dont pad the 3' and 5' ends
  foreach ( my $i = 0 ; $i <=  $#exon_clusters ; $i++ ){
    my $exon = $exon_clusters[$i];
    my $start_padding = 20;
    my $end_padding = 20;
    $start_padding = 0 if $i == 0;
    $end_padding = 0 if $i == $#exon_clusters;
    my $padded_exon =  create_Exon
      (
       $exon->start - $start_padding,
       $exon->end + $end_padding ,
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
      $last_exon->end(   $last_exon->end  + $trim );
      $exon->start( $last_exon->end+1 );
    }
    # check they dont overlap
    #if ( $last_exon->end >= $exon->start ) {
      #print "OVERLAP at padding" .  $last_exon->end . " " .  $exon->start ."\n";
    #}
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
  my ($self) = @_;
  my @features = @{$self->reads};
  my @clusters;
  foreach my $f ( @features ) {
    my $clustered = 0;
    foreach my $cluster ( @clusters ) {
      # do they overlap?
      if ( $f->[2] <= $cluster->{'end'} 
	   &&  $f->[3] >= $cluster->{'start'} ) {
	# Expand the cluster
	$cluster->{'start'} = $f->[2] 
	  if $f->[2] < $cluster->{'start'};
	$cluster->{'end'} = $f->[3]
	  if $f->[3]   > $cluster->{'end'}; 
	push @{$cluster->{'features'}}, $f;
	$clustered = 1;
      }
    }
    unless ($clustered) {
      my $cluster;
      push @{$cluster->{'features'}}, $f;
      $cluster->{'start'} = $f->[2];
      $cluster->{'end'}   = $f->[3];
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
  #print " Got ". scalar(@genes) ."\n";
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

      #print "merging between " .$left_gene->end ." and ".  $right_gene->start  ."\n";
      # is the gap covered by a repeat?
      my @repeats = sort { $a->start <=> $b->start } @{$self->repeats} ;
      my $repeat_coverage = 0;
      # so the repeats are now non-overlapping blocks ( if there is more than one )
      foreach my $repeat ( @repeats ) {
      next unless $repeat->start <= $right_gene->start && $repeat->end >= $left_gene->end;
      last if $repeat->start > $right_gene->start;
      $repeat->start($left_gene->end) if  $repeat->start < $left_gene->end;
      $repeat->end($right_gene->start) if $repeat->end > $right_gene->start;
      $repeat_coverage+= $repeat->end - $repeat->start;
      }  
      $repeat_coverage /= $gap;
      # merge the genes together where repeat coverage is 100%
      if ($repeat_coverage >= 0.95   ) {
	#print "Do some merging\n";
	#print "Repeat Coverage = $repeat_coverage \n";
	#print $left_gene->slice->name ."\n";
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
	#print "NEW GENE " . $new_gene->start ." ". $new_gene->end ."\n";
	# take them out of the array and carry on
	splice (@genes,$i-1,2,$new_gene);
	$i-= 1;
	#print "Merged 1 successfully\n";
      }
    }
  }
  # add whatever is left to the return array
  push @merged_genes, @genes;
  #print "Returning " . scalar (@merged_genes) . "\n";
  return \@merged_genes;
}

sub sql {
  my ($self,$query,$db) = @_;
#  print STDERR "query $query\n";
  my $sth = $db->dbc->prepare($query);
  $sth->execute();
  my @array = $sth->fetchrow_array;
  return \@array;
}

sub sql_array {
  my ($self,$query,$db) = @_;
#  print STDERR "query $query\n";
  my $sth = $db->dbc->prepare($query);
  $sth->execute();
  return $sth;
}

sub ungapped_features {
  my ($self,$feat) = @_;
  my @ugfs;
  
  my $string = $feat->[4];
  my $start = $feat->[2];
  my $end = $feat->[3];
 # print "THINGS $start $end $string\n";
  my @pieces = ( $string =~ /(\d*[MDI])/g );
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
      $ugf->[0] = $feat->[0];
      $ugf->[1] = $feat->[1];
      $ugf->[2] = $qstart;
      $ugf->[3] = $qend;
      $ugf->[4] = $length."M";
      push @ugfs, $ugf;
#      print "UNGAPPED " .$ugf->[2] .
	" " . $ugf->[3] . " " . $ugf->[4] ."\n";
    } elsif( $piece =~ /I$/ ) {
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

sub make_repeat_blocks {
  my ($self,$repeats_ref) = @_;

  my @repeats = sort { $a->start <=> $b->start }@$repeats_ref;
 # print " got " . scalar(@repeats) . " repeat blocks initially\n";
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
  #print " got " . scalar(@repeats) . " repeat blocks after merging\n";
  return \@repeats;
}

###########################################
# Containers


sub reads {
  my ($self, $value) = @_;
  if (defined $value ) {
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

#sub repeat_slice_adaptor {
#  my ($self, $val) = @_;
#
#  if (defined $val) {
#    $self->{_rsa} = $val;
#  }
#
#  return $self->{_rsa};
#}

sub repeats {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_repeats} = $val;
  }

  return $self->{_repeats};
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

sub DNA_ALIGN_FEATURE_DATA {
  my ($self,$value) = @_;

  if (defined $value) {
    $self->{'_CONFIG_DNA_ALIGN_FEATURE_DATA'} = $value;
  }
  
  if (exists($self->{'_CONFIG_DNA_ALIGN_FEATURE_DATA'})) {
    return $self->{'_CONFIG_DNA_ALIGN_FEATURE_DATA'};
  } else {
    return undef;
  }
}


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


sub USE_ANALYSIS_LOGIC_NAME_AS_DEFAULT_GENE_OUTPUT_BIOTYPE {
  my ($self,$value) = @_;

  if ( defined $value ) {
     $self->{'_CONFIG_USE_ANALYSIS_LOGIC_NAME_AS_DEFAULT_GENE_OUTPUT_BIOTYPE'} = $value;
  }
  if (exists($self->{'_CONFIG_USE_ANALYSIS_LOGIC_NAME_AS_DEFAULT_GENE_OUTPUT_BIOTYPE'})) {
    return $self->{'_CONFIG_USE_ANALYSIS_LOGIC_NAME_AS_DEFAULT_GENE_OUTPUT_BIOTYPE'};
  } else {
    return 0;
  }
}


sub INTRON_PAIR {  
   warning("not implemented in lite version\n"); 
}
sub MIN_SINGLE_EXON_LENGTH {  
   warning("not implemented in lite version\n"); 
}  



sub get_dna_align_features_from_databases  { 
  my ($self,$slice) = @_; 

  my %dbnames_2_logicnames = %{ $self->DNA_ALIGN_FEATURE_DATA} ; 

  my @all_fetched_daf_reads;  

  foreach my $db_hash_key ( keys %dbnames_2_logicnames )  {
     #rint "querying $db_hash_key\n";
     # get standard  DBAdaptor ('core') with no dna-database attached 
     my $set_db = $self->get_dbadaptor($db_hash_key,undef,0);
     $set_db->disconnect_when_inactive(0);
     # get analysis ids by configured logic_name 
            
     # get seq_region_id of slice  
     my $new_slice = $set_db->get_SliceAdaptor->fetch_by_region(undef,$slice->seq_region_name);   
     my $seq_id = $new_slice->dbID ; 
     #my $seq_id =  $self->sql( "SELECT seq_region_id from seq_region WHERE name = '" .  $slice->seq_region_name . "'",$self->db)->[0] ;

     my $sql = "SELECT dna_align_feature_id, seq_region_id, seq_region_start, seq_region_end, cigar_line from dna_align_feature ".
               " where seq_region_id = $seq_id and seq_region_end >= " . $slice->start .  " and seq_region_start <= " .$slice->end ." ";  

     my @logic_names  = @{$dbnames_2_logicnames{$db_hash_key}};  

     # PRIVATE 
     
     sub _get_daf_analysis_ids_from_logic_names { 
        my ($self,$set_db,$lg ) = @_ ; 
        my @daf_analysis_ids; 
        my @logic_names  = @$lg;
        for my $lg ( @logic_names ) {  
          my $analysis = $set_db->get_AnalysisAdaptor->fetch_by_logic_name($lg);  
          push @daf_analysis_ids, $analysis->dbID ; 
       } 
       return \@daf_analysis_ids;
     }



     if ( scalar(@logic_names) == 0 ) { 
          throw("you have't defined any data to fetch out of the dna_align_feature table of $db_hash_key - use [\"*\"]  ". 
          " if you want to fetch all data or specify a logic_name of the dna_align_features you like to use for the analysis \n");
     }

     if (scalar(@logic_names) == 1 && $logic_names[0] eq "\*" ) { 
          warning("fetching all data from dna align feature table \n");
          $sql .= " ;" # fetch all data from daf table 
     } else { 
         my $daf_analysis_ids = $self->_get_daf_analysis_ids_from_logic_names($set_db,\@logic_names); 

         if (scalar(@$daf_analysis_ids) == 1 ){  
            $sql .= " and analysis_id = $$daf_analysis_ids[0] ;";
         } else { 
            # we have to fetch daf' from multiple aanalysis wit h diff analysis ids  
            $sql .= " and analysis_id in (".join(",",@$daf_analysis_ids). ");";  
         }
     }

     #print "\n\nSQL $sql\n\n"; 
     
    my $handle =  $self->sql_array( $sql, $set_db); 
    my @fetched_daf_reads; 

    while ( my @read = $handle->fetchrow_array ) {
      push @fetched_daf_reads, \@read;
    }

    if ( scalar(@fetched_daf_reads) == 0 ) {  
      throw("Configuration-error : NO DNA_ALIGN_FEATURES FOUND IN $db_hash_key - check your config !\n");
    }  

    push @all_fetched_daf_reads, @fetched_daf_reads;
  }
  @all_fetched_daf_reads = sort { $a->[2] <=> $b->[2] } @all_fetched_daf_reads;
  #print "Got " . scalar(@all_fetched_daf_reads) . " in total \n"; 
  return \@all_fetched_daf_reads;
}



1;
