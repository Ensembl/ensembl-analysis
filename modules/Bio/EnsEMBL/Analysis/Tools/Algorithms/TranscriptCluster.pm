=head1 NAME

TranscriptCluster

=head1 SYNOPSIS


=head1 DESCRIPTION

This object holds one or more transcripts which have been clustered according to 
comparison criteria external to this class (for instance, in the 
method _compare_Transcripts of the class GeneComparison).
Each TranscriptCluster object holds the IDs of the transcripts clustered and the beginning and end coordinates
of each one (taken from the start and end coordinates of the first and last exon in the correspondig
get_all_Exons array)

It inherits from Bio::EnsEMBL::Root and Bio::RangeI. A TranscriptCluster is a range in the sense that it convers
a defined extent of genomic sequence. It is also possible to check whether two clusters overlap (in range),
is included into another cluster, etc...

=cut

# Let the code begin ...

package Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster;

use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::ExonCluster;
use Bio::EnsEMBL::Analysis::Tools::Algorithms::IntronCluster;
use strict;

use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Root;
use Bio::RangeI;
use Bio::EnsEMBL::Utils::Exception qw (throw warning ) ; 
use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Root Bio::RangeI);

=head1 METHODS

=cut

#########################################################################


=head2 new()

new() initializes the attributes:
_transcript_array
_transcriptID_array
_start
_end

=cut

sub new {
  my ($class,$whatever)=@_;

  if (ref($class)){
    $class = ref($class);
  }
  my $self = {};
  bless($self,$class);
  
  if ($whatever){
    throw( "Can't pass an object to new() method. Use put_Genes() to include Bio::EnsEMBL::Gene in cluster");
  }
  $self->{_biotype}={} ; 
  return $self;
}
#########################################################################

=head1 Range-like methods

Methods start and end are typical for a range. We also implement the boolean
and geometrical methods for a range.

=head2 start()

  Title   : start
  Usage   : $start = $transcript_cluster->end();
  Function: get/set the start of the range covered by the cluster. This is re-calculated and set everytime
            a new transcript is added to the cluster
  Returns : a number
  Args    : optionaly allows the start to be set

=cut

sub start{
  my ($self,$start) = @_;
  if ($start){
    throw( "$start is not an integer") unless $start =~/^[-+]?\d+$/;
    $self->{'_start'} = $start;
  }
  return $self->{'_start'};
}

############################################################

=head2 end()

  Title   : end
  Usage   : $end = $transcript_cluster->end();
  Function: get/set the end of the range covered by the cluster. This is re-calculated and set everytime
            a new transcript is added to the cluster
  Returns : the end of this range
  Args    : optionaly allows the end to be set
          : using $range->end($end
=cut

sub end{
  my ($self,$end) = @_;
  if ($end){
    throw( "$end is not an integer") unless $end =~/^[-+]?\d+$/;
    $self->{'_end'} = $end;
  }
  return $self->{'_end'};
}

############################################################

=head2 length

  Title   : length
  Usage   : $length = $range->length();
  Function: get/set the length of this range
  Returns : the length of this range
  Args    : optionaly allows the length to be set
          : using $range->length($length)

=cut

sub length{
  my $self = shift @_;
  if (@_){
    $self->confess( ref($self)."->length() is read-only");
  }
  return ( $self->{'_end'} - $self->{'_start'} + 1 );
}

############################################################

=head2 strand

  Title   : strand
  Usage   : $strand = $transcript->strand();
  Function: get/set the strand of the transcripts in the cluster.
            The strand is set in put_Transcripts when the first transcript is added to the cluster, in that
            that method there is also a check for strand consistency everytime a new transcript is added
            Returns 0/undef when strand not set
  Returns : the strandidness (-1, 0, +1)
  Args    : optionaly allows the strand to be set
        
=cut

sub strand{
  my ($self,$strand) = @_;
  if ($strand){
    $self->{'_strand'} = $strand;
  }
  return $self->{'_strand'};
}


############################################################

=head1 Boolean Methods

These methods return true or false. They throw an error if start and end are
not defined. They are implemented in Bio::RangeI.

 $cluster->overlaps($other_cluster) && print "Clusters overlap\n";

=head2 overlaps

  Title   : overlaps
  Usage   : if($cluster1->overlaps($cluster)) { do stuff }
  Function: tests if $cluster2 overlaps $cluster1 overlaps in the sense of genomic-range overlap,
            it does NOT test for exon overlap.
  Args    : arg #1 = a range to compare this one to (mandatory)
            arg #2 = strand option ('strong', 'weak', 'ignore') 
  Returns : true if the clusters overlap, false otherwise

=cut

=head2 contains

  Title   : contains
  Usage   : if($cluster1->contains($cluster2) { do stuff }
  Function: tests whether $cluster1 totally contains $cluster2 
  Args    : arg #1 = a range to compare this one to (mandatory)
	             alternatively, integer scalar to test
            arg #2 = strand option ('strong', 'weak', 'ignore') (optional)
  Returns : true if the argument is totaly contained within this range

=cut

=head2 equals

  Title   : equals
  Usage   : if($cluster1->equals($cluster2))
  Function: test whether the range covered by $cluster1 has the same start, end, length as the range 
            covered by $cluster2
  Args    : a range to test for equality
  Returns : true if they are describing the same range

=cut

=head1 Geometrical methods

These methods do things to the geometry of ranges, and return
Bio::RangeI compliant objects or triplets (start, stop, strand) from
which new ranges could be built. They are implemented  here and not in Bio::RangeI, since we
want them to return a new TranscriptCluster object.

=head2 overlap_extent

 Title   : overlap_extent
 Usage   : ($a_unique,$common,$b_unique) = $a->overlap_extent($b)
 Function: Provides actual amount of overlap between two different
           ranges. Implemented already in RangeI
 Example :
 Returns : array of values for 
           - the amount unique to a
           - the amount common to both
           - the amount unique to b
 Args    : 

=cut

#########################################################################

=head2 intersection

  Title   : intersection
  Usage   : $intersection_cluster = $cluster1->intersection($cluster2)
  Function: gives a cluster with the transcripts which fall entirely within the intersecting range of
            $cluster1 and $cluster2
  Args    : arg #1 = a cluster to compare this one to (mandatory)
            arg #2 = strand option ('strong', 'weak', 'ignore') (not implemented here)
  Returns : a TranscriptCluster object, which is empty if the intersection does not contain
            any transcript
=cut

sub intersection{
  my ($self, $cluster, $ignore_strand) = @_;

  # if either is empty, return an empty cluster
  if ( scalar( @{$self->get_Transcripts}) == 0 ){
    warning( "cluster $self is empty, returning an empty TranscriptCluster");
    my $empty_cluster = Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster->new() ;
    return $empty_cluster;
  }

  if ( scalar( @{$cluster->get_Transcripts} ) == 0 ){
    warning( "cluster $cluster is empty, returning an empty TranscriptCluster");
    my $empty_cluster= Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster->new() ;
    return $empty_cluster;
  }

  my @transcripts = @{$self->get_Transcripts} ; 
  push( @transcripts, @{$cluster->get_Transcripts} );

  # make an unique list of transcripts, in case they are repeated
  my %list;
  foreach my $transcript (@transcripts){
    $list{$transcript} = $transcript;
  }
  @transcripts = values( %list );

  my ($inter_start,$inter_end);
  my ($start1,$end1) = ($self->start   ,   $self->end);
  my ($start2,$end2) = ($cluster->start,$cluster->end);

  # my $strand = $cluster->strand;

  # if clusters overlap, calculate the intersection
  if ( $self->overlaps( $cluster ) ){
    $inter_start = ($start2 > $start1) ? $start2 : $start1;
    $inter_end = ($end2 < $end1) ? $end2 : $end1;
  }
  else{
    warning( "clusters $self and $cluster do not intersect range-wise, returning an empty TranscriptCluster");
    my $empty_cluster =Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster->new(); 
    return $empty_cluster;
  }

  my $inter_cluster =Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster->new(); 
  my @inter_transcripts;

  # see whether any transcript falls within this intersection
  foreach my $transcript ( @transcripts ){
    my ($start,$end) = $self->_get_start_end($transcript);
    if ($start >= $inter_start && $end <= $inter_end ){
       $inter_cluster->put_Transcripts( [$transcript],  $ignore_strand );
    }
  }

  if ( scalar( @{$inter_cluster->get_Transcripts} ) == 0 ){
     warning( "cluster $inter_cluster is empty, returning an empty TranscriptCluster");
     return $inter_cluster;
  }
  else{
    return $inter_cluster;
  }
}

############################################################

=head2 union

  Title   : union
  Usage   : $union_cluster = $cluster1->union(@clusters);
  Function: returns the union of clusters 
  Args    : a TranscriptCluster or list of TranscriptClusters to find the union of
  Returns : the TranscriptCluster object containing all of the ranges

=cut

sub union{
my ($self,@clusters, $ignore_strand) = @_;

if ( ref($self) ){
 unshift @clusters, $self;
}

my $union_cluster = Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster->new() ;
my $union_strand;

foreach my $cluster (@clusters){
 unless ($union_strand){
  $union_strand = $cluster->strand;
 }
 unless ( $cluster->strand == $union_strand){
  warning("You're making the union of clusters in opposite strands");
 }
 $union_cluster->put_Transcripts( $cluster->get_Transcripts, $ignore_strand );
}

return $union_cluster;


}

############################################################

=head2 put_Transcripts()

  function to include one or more transcripts in the cluster.
  Useful when creating a cluster. It takes as argument an array of transcripts, it returns nothing.

=cut

sub put_Transcripts {
  my ($self, $trans_array_ref, $ignore_strand )= @_; 
 unless( ref($trans_array_ref)=~m/ARRAY/) {
    throw("Only take array ref $trans_array_ref!\n") ;
  }

  my @new_transcripts = @$trans_array_ref ; 

  throw("Can't add no transcripts to ".$self) if(!@new_transcripts || !$new_transcripts[0]);
  if ( !$new_transcripts[0]->isa('Bio::EnsEMBL::Transcript') ){
    throw( "Can't accept a [ $new_transcripts[0] ] instead of a Bio::EnsEMBL::Transcript");
  }

  #Get bounds of new transcripts
  my $min_start = undef;
  my $max_end   = undef;
  foreach my $transcript (@new_transcripts) {
    if ($transcript->stable_id) {
        #print "got transcript " . $transcript->stable_id . "\n" ;
    } else {
        #print "got transcript " . $transcript->display_id . "\n";
    }
    my ($start, $end) = $self->_get_start_end($transcript);
    if (!defined($min_start) || $start < $min_start) {
      $min_start = $start;
    }
    if (!defined($max_end) || $end > $max_end) {
      $max_end = $end;
    }
  }

  # check strand consistency among transcripts
  foreach my $transcript (@new_transcripts){
    if (!$ignore_strand) {
      my @exons = @{$transcript->get_all_Exons};
      unless ( $self->strand ){
        $self->strand( $exons[0]->strand );
      }
      if ( $self->strand != $exons[0]->strand ){
        warning( "You're trying to put $transcript in a cluster of opposite strand");
      }
    } else {
      # we can ignore the strand, and do nothing
    }
  }

  # if start is not defined, set it
  unless ( $self->start ){
    $self->start( $min_start );
  }

  # if end is not defined, set it
  unless ( $self->end ){
    $self->end( $max_end);
  }
  
  # extend start and end if necessary as we include more transcripts
  if ($min_start < $self->start ){
    $self->start( $min_start );
  }
  if ( $max_end > $self->end ){
    $self->end( $max_end );
  }

  push ( @{ $self->{'_transcript_array'} }, @new_transcripts );

}

#########################################################################

=head2 get_Transcripts()

  it returns the array of transcripts in the GeneCluster object

=cut

sub get_Transcripts {
  my $self = shift @_;
  return $self->{'_transcript_array'};
}

#########################################################################

=head2 to_String()

  it returns a string containing the information about the transcripts in the TranscriptCluster object

=cut

sub to_String {
  my $self = shift @_;
  my $data='';
  foreach my $tran ( @{ $self->{'_transcript_array'} } ){
    my @exons = @{$tran->get_all_Exons};
    my $id;
    if ( $tran->stable_id ){
      $id = $tran->stable_id;
    }
    else{
      $id = $tran->dbID;
    }
 
    $data .= sprintf "Id: %-16s"             , $id;
    $data .= sprintf "Contig: %-21s"         , $exons[0]->contig->id;
    $data .= sprintf "Exons: %-3d"           , scalar(@exons);
    my ($start, $end) = $self->_get_start_end($tran);
    $data .= sprintf "Start: %-9d"           , $start;
    $data .= sprintf "End: %-9d"             , $end;
    $data .= sprintf "Strand: %-3d"          , $exons[0]->strand;
    $data .= sprintf "Exon-density: %3.2f\n", $self->exon_Density($tran);
  }
  return $data;
}

#########################################################################

sub exon_Density{
  my ($self, $transcript) = @_;  
  my $density;
  my $exon_span;
  my @exons = @{$transcript->get_all_Exons};
  @exons = sort { $a->start <=> $b->start } @exons;
  my $transcript_length = $exons[$#exons]->end - $exons[0]->start;
  foreach my $exon ( @exons ){
    $exon_span += $exon->length;
  }
  $density = $exon_span/$transcript_length;
  return $density;
}

=head2 _get_start_end()

 function to get the start and end positions - written as one method
 for efficiency

=cut

sub _get_start_end {
  my ($self,$t) = @_;
  my @exons = sort { $a->start <=> $b->start } @{$t->get_all_Exons}; 
  my $start = $exons[0]->start;
  my $end   = $exons[-1]->end;
  return ($start, $end, $exons[0]->strand);
}


=head2 register_biotype

 function to set included biotypes in this cluster of transcripts 

=cut
 
sub register_biotype {
  my ($self,$biotype)  =@_ ; 
  ${$self->{_biotypes} }{$biotype}=() ; 
 
}

sub get_biotypes_included {
  my ($self) = @_ ; 
  my @bt = keys %{$self->{_biotypes} } ; 
  return \@bt ; 
}


sub get_ExonCluster {
  my ( $self , $ignore_strand) = @_;
  my @clusters;

  foreach my $trans (@{$self->get_Transcripts}) {
    my $tr_biotype = $trans->biotype;

    foreach my $exon (@{$trans->get_all_Exons}) {

      my @matching_clusters;
        #print "\nExon " . $exon->dbID . " limits: " . $exon->start . 
       # " and " .  $exon->end . "\n";

      CLUSTER: foreach my $cluster (@clusters) {
         #print "Testing against cluster with limits " . 
         #$cluster->start. " to " . $cluster->end . " ".$cluster->strand ."\t";
        if (!($exon->start >= $cluster->end ||
              $exon->end <= $cluster->start)) {
           if ($cluster->strand eq $exon->strand ){
              push (@matching_clusters, $cluster);
              #print "cl. matches " .$cluster->strand ."\t" .$exon->strand . "\t" .$trans->strand ."\n" ;
           }
        }
        #print "\n";
      }
      if (scalar(@matching_clusters) == 0) {
        # print STDERR "Created new cluster for " . $exon->stable_id . " " . $exon->dbID . "\n";
        # print "\ncreating new cluster for Exon " . $exon->dbID . 
        # " limits: " . $exon->start . " and " .  $exon->end . "\n";

        my $newcluster = Bio::EnsEMBL::Analysis::Tools::Algorithms::ExonCluster->new() ;

        $newcluster->add_exon($exon,$trans,$ignore_strand);
        push(@clusters,$newcluster);

      } elsif (scalar(@matching_clusters) == 1) {
         #print STDERR "Adding to cluster for " . $exon->stable_id . " " . $exon->dbID . "\n";
        $matching_clusters[0]->add_exon($exon,$trans, $ignore_strand);
      } else {
         # Merge the matching clusters into a single cluster
        print STDERR "Merging clusters for " . $exon->dbID ."\n";
        my @new_clusters;
        my $merged_cluster = Bio::EnsEMBL::Analysis::Tools::Algorithms::ExonCluster->new() ;

        foreach my $clust (@matching_clusters) {
          $merged_cluster->merge($clust, $ignore_strand);
        }
        $merged_cluster->add_exon($exon,$trans, $ignore_strand);
        push @new_clusters,$merged_cluster;

        # Add back non matching clusters
        foreach my $clust (@clusters) {
          my $found = 0;
          MATCHING: foreach my $m_clust (@matching_clusters) {
            if ($clust == $m_clust) {
              $found = 1;
              last MATCHING;
            }
          }
          if (!$found) {
            push @new_clusters,$clust;
          }
        }
        @clusters = @new_clusters;
      }
    }
  }

  # setting exon/cluster relationship 
  for my $c (@clusters) {
    for my $e(@{ $c->get_all_Exons_in_ExonCluster} ) {
      # exon has to be Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object 
      $e->cluster($c) ;
    }
  }
  return \@clusters;
}

sub get_ExonCluster_using_all_Exons{
  my ( $self, $ignore_strand) = @_;
  my @clusters;

  my @transcripts;
  if (ref($self) eq 'Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster') {
    @transcripts = @{$self->get_Transcripts};
  } elsif (ref($self) eq 'ARRAY' && $self->[0]->isa('Bio::EnsEMBL::Transcript')) {
    throw "get_ExonCluster_using_all_Exons called with not allowed argument !\n" ; 
  } else {
    throw("Not allowed");
  }


  foreach my $trans (@transcripts) {
    my $tr_biotype = $trans->biotype;

    foreach my $exon (@{$trans->get_all_Exons}) {

      my @matching_clusters;
       # print "\nStarting Exon " . $exon->dbID . " limits: " . $exon->start . 
       #" and " .  $exon->end . "\n";

      CLUSTER: foreach my $cluster (@clusters) {
        # print "Testing against cluster with limits " . 
         $cluster->start. " to " . $cluster->end . " ".$cluster->strand ."\t";
        if (!($exon->start >= $cluster->end ||
              $exon->end <= $cluster->start)) {
          if (!$ignore_strand) {
             if ($cluster->strand eq $exon->strand ){
                push (@matching_clusters, $cluster);
                #print "cl. matches " .$cluster->strand ."\t" .$exon->strand . "\t" .$trans->strand ."\n" ;
             }
          } else {
            # ignore strand; do nothing
            push (@matching_clusters, $cluster);
          }
        }
        #print "\n";
      }
      if (scalar(@matching_clusters) == 0) {
        # print STDERR "Created new cluster for " . $exon->stable_id . " " . $exon->dbID . "\n";
        # print "\ncreating new cluster for Exon " . $exon->dbID . 
        # " limits: " . $exon->start . " and " .  $exon->end . "\n";

        my $newcluster = Bio::EnsEMBL::Analysis::Tools::Algorithms::ExonCluster->new() ;

        $newcluster->add_exon($exon,$trans, $ignore_strand);
        push(@clusters,$newcluster);

      } elsif (scalar(@matching_clusters) == 1) {
        # print STDERR "Adding to cluster for " . $exon->stable_id . " " . $exon->dbID . "\n";
        #$matching_clusters[0]->add_exon($exon,$trans);
        $matching_clusters[0]->add_exon_if_not_present($exon,$trans, $ignore_strand);
      } else {
         # Merge the matching clusters into a single cluster
        print STDERR "Merging clusters for " . $exon->dbID ."\n";
        my @new_clusters;
        my $merged_cluster = Bio::EnsEMBL::Analysis::Tools::Algorithms::ExonCluster->new() ;

        foreach my $clust (@matching_clusters) {
          $merged_cluster->merge_new_exon($clust, $ignore_strand);
        }
        $merged_cluster->add_exon_if_not_present($exon,$trans, $ignore_strand);
        push @new_clusters,$merged_cluster;

        # Add back non matching clusters
        foreach my $clust (@clusters) {
          my $found = 0;
          MATCHING: foreach my $m_clust (@matching_clusters) {
            if ($clust == $m_clust) {
              $found = 1;
              last MATCHING;
            }
          }
          if (!$found) {
            push @new_clusters,$clust;
          }
        }
        @clusters = @new_clusters;
      }
    }
  }

  # setting exon/cluster relationship
  for my $c (@clusters) {
    for my $e(@{ $c->get_all_Exons_in_ExonCluster} ) {
      # exon has to be Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonExtended object 
      $e->cluster($c) ;
    }
  }

  return \@clusters;
}


sub get_coding_ExonCluster{
  my ( $self, $ignore_strand ) = @_;
  my @clusters;

  foreach my $trans (@{$self->get_Transcripts}) {
    my $tr_biotype = $trans->biotype;

    foreach my $exon (@{$trans->get_all_translateable_Exons}) {

      my @matching_clusters;
       # print "\nStarting Exon " . $exon->dbID . " limits: " . $exon->start . 
       #" and " .  $exon->end . "\n";

      CLUSTER: foreach my $cluster (@clusters) {
        # print "Testing against cluster with limits " . 
         # $cluster->start. " to " . $cluster->end . " ".$cluster->strand ."\t";
        if (!($exon->start >= $cluster->end ||
              $exon->end <= $cluster->start)) {
           if (!$ignore_strand) {
             if ($cluster->strand eq $exon->strand ){
                push (@matching_clusters, $cluster);
                #print "cl. matches " .$cluster->strand ."\t" .$exon->strand . "\t" .$trans->strand ."\n" ;
             }
          } else {
            # we don't care if exons are on different strands
            push (@matching_clusters, $cluster);
          }
        }
        #print "\n";
      }
      if (scalar(@matching_clusters) == 0) {
        # print STDERR "Created new cluster for " . $exon->stable_id . " " . $exon->dbID . "\n";
        # print "\ncreating new cluster for Exon " . $exon->dbID . 
        # " limits: " . $exon->start . " and " .  $exon->end . "\n";

        my $newcluster = Bio::EnsEMBL::Analysis::Tools::Algorithms::ExonCluster->new() ;

        $newcluster->add_exon($exon,$trans, $ignore_strand);
        push(@clusters,$newcluster);

      } elsif (scalar(@matching_clusters) == 1) {
        # print STDERR "Adding to cluster for " . $exon->stable_id . " " . $exon->dbID . "\n";
        #$matching_clusters[0]->add_exon($exon,$trans);
        $matching_clusters[0]->add_exon_if_not_present($exon,$trans, $ignore_strand);
      } else {
         # Merge the matching clusters into a single cluster
        print STDERR "Merging clusters for " . $exon->dbID ."\n";
        my @new_clusters;
        my $merged_cluster = Bio::EnsEMBL::Analysis::Tools::Algorithms::ExonCluster->new() ;

        foreach my $clust (@matching_clusters) {
          $merged_cluster->merge_new_exon($clust, $ignore_strand);
        }
        $merged_cluster->add_exon_if_not_present($exon,$trans, $ignore_strand);
        push @new_clusters,$merged_cluster;

        # Add back non matching clusters
        foreach my $clust (@clusters) {
          my $found = 0;
          MATCHING: foreach my $m_clust (@matching_clusters) {
            if ($clust == $m_clust) {
              $found = 1;
              last MATCHING;
            }
          }
          if (!$found) {
            push @new_clusters,$clust;
          }
        }
        @clusters = @new_clusters;
      }
    }
  }

  return \@clusters;
}

sub get_IntronClusters {
  my ( $self , $ignore_strand) = @_;
  my @clusters;

  foreach my $trans (@{$self->get_Transcripts}) {

    my @introns = @{get_all_Introns_between_translateable_Exons($trans)};
    # print "\n\n\nTranscript ".$trans->stable_id." ".$trans->biotype." with ".scalar(@introns)." introns\n";
    my $tr_biotype = $trans->biotype;


    foreach my $intron (@{get_all_Introns_between_translateable_Exons($trans)}) {
      #my ($intron_start, $intron_end, $intron_strand) = get_Intron_info($intron);

      my @matching_clusters;
       # print "\nStarting Intron with limits: " . $intron_start . 
       #" and " .  $intron_end . ", ".$intron->start." and ".$intron->end."\n";

      CLUSTER: foreach my $cluster (@clusters) {
        # print "Testing against cluster with limits " . 
        # $cluster->start. " to " . $cluster->end . " ".$cluster->strand ."\n";
        if ($intron->start == $cluster->start && $intron->end == $cluster->end) {
          if (!$ignore_strand) {
             if ($cluster->strand eq $intron->strand ){
                push (@matching_clusters, $cluster);
                #print "cl. matches " .$cluster->strand ."\t" .$intron_strand . "\t" .$trans->strand ."\n" ;
             }
          } else {
            push (@matching_clusters, $cluster);
          }
        }
        #print "\n";
      }
      if (scalar(@matching_clusters) == 0) {
        # print STDERR "Created new cluster for " . $intron->stable_id . " " . $intron->dbID . "\n";
        # print "\ncreating new cluster for Intron " . $intron->dbID . 
        # " limits: " . $intron_start . " and " .  $intron_end . "\n";

       # print "Making new intron cluster...\n";
        my $newcluster = Bio::EnsEMBL::Analysis::Tools::Algorithms::IntronCluster->new() ;
        $newcluster->{'_types_sets'} = $self->{'_types_sets'};

        $newcluster->put_Introns([$intron],$trans);
        push(@clusters,$newcluster);
        
      } elsif (scalar(@matching_clusters) == 1) {
        # print STDERR "Adding to cluster for " . $intron->stable_id . " " . $intron->dbID . "\n";
        #$matching_clusters[0]->put_Intron($intron,$trans);
        $matching_clusters[0]->put_Introns([$intron],$trans);
      } else {
        throw("Intron cannot match more than one cluster or it will not be unique");
      }
    }
  }

  return \@clusters;
}

sub get_all_Introns_between_translateable_Exons {
  my ($transcript) = @_;
  my @introns;

  my @exons = @{$transcript->get_all_translateable_Exons}; 

  for (my $i = 0; $i < scalar(@exons) - 1; $i++) {
    my $intron = new Bio::EnsEMBL::Intron($exons[$i],$exons[$i+1]);
    push @introns, $intron;
  }

  return \@introns;
}

#sub get_Intron_info {
#  my ($intron) = @_;
#
#  if (ref($intron) !~ m/Intron/) {
#    throw("Cannot get_Intron_info for object $intron");
#  }
#
#  my $intron_strand = $intron->prev_Exon->strand;
#  if ($intron->next_Exon->strand != $intron_strand) {
#    throw("Intron bounded by 2 exons on different strands");
#  }
#
#  my ($intron_start, $intron_end);
#  if ($intron_strand == 1) {
#    $intron_start = $intron->prev_Exon->end + 1;
#    $intron_end = $intron->next_Exon->start - 1;
#  } elsif ($intron_strand == -1) {
#    $intron_start = $intron->next_Exon->end + 1;
#    $intron_end = $intron->prev_Exon->start - 1;
#  } else {
#    throw("Strand not recognised");
#  }
#  return ($intron_start, $intron_end, $intron_strand);
#}
1;
