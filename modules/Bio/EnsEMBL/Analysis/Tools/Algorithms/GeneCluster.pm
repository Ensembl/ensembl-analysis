=head1 NAME

GeneCluster

=head1 SYNOPSIS


=head1 DESCRIPTION

This object holds one or more genes which has been clustered according to 
comparison criteria external to this class (for instance, in the 
methods compare and _compare_Genes methods of the class GeneComparison).
Each GeneCluster object holds the IDs of the genes clustered and the beginning and end coordinates
of each one (taken from the start and end coordinates of the first and last exon in the correspondig
get_all_Exons array)

=head1 CONTACT

eae@sanger.ac.uk

=cut

# Let the code begin ...

package Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Utils::Exception qw(throw warning );
#use Bio::EnsEMBL::Analysis::Tools::Algorithms::ClusterUtils qw( genes_to_Transcript_Cluster ) ;

@ISA = qw(Bio::EnsEMBL::Root);

=head1 METHODS

=cut

#########################################################################


=head2 new()

new() initializes the attributes:

$self->{'_benchmark_types'}
$self->{'_prediction_types'}
$self->{'_benchmark_genes'}
$self->{'_prediction_genes'}

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

  $self->{_cached_start}  = undef;
  $self->{_cached_end}    = undef;
  $self->{_cached_strand} = undef;
  return $self;
}

#########################################################################

=head2 put_Genes()

  function to include one or more genes in the cluster.
  Useful when creating a cluster. It takes as argument an array of genes, it returns nothing.

=cut

sub put_Genes {
  my ($self, @new_genes)= @_;
  if ( !defined( $self->{'_types_sets'} ) ){
    throw( "Cluster lacks references to gene-types, unable to put the gene");
  }


 GENE:
  foreach my $gene (@new_genes){
    throw("undef for gene") if (!$gene);

    foreach my $set_name ( keys %{$self->{'_types_sets'}}) {

      my $set = $self->{'_types_sets'}{$set_name};
      foreach my $type ( @{$set} ){

        if ($gene->type eq $type) {
          push ( @{ $self->{'_gene_sets'}{$set_name} }, $gene );
          next GENE; 
        } elsif (ref($gene)=~m/PredictionTranscript/)  { 

          # gene is prediction_transcript 
          push ( @{ $self->{'_gene_sets'}{$set_name} }, $gene );
          next GENE ;  

        }
      }
    }
    throw("Failed putting gene of type " . $gene->type . "\n");
  }
  $self->{_cached_start}  = undef;
  $self->{_cached_end}    = undef;
  $self->{_cached_strand} = undef;
}



sub get_sets_included {
  my $self = shift;
  my @included_sets;

  foreach my $set_name ( keys %{$self->{'_types_sets'}}) {
    if (defined( $self->{'_gene_sets'}{$set_name})) {
      push @included_sets,$set_name;
    }
  }
  return \@included_sets;
}


#########################################################################

=head2 get_Genes()

  it returns the array of genes in the GeneCluster object

=cut

sub get_Genes {
  my $self = shift @_;

  my @genes;
  if (!defined( $self->{'_gene_sets'} ) ) {
    $self->warning("The gene array you try to retrieve is empty");
    @genes = ();
  }

  foreach my $set_name (keys %{$self->{'_gene_sets'}}) {
    push( @genes, @{ $self->{'_gene_sets'}{$set_name} } );
  }

  return @genes;
}

sub get_all_Exons {
  my $self = shift @_;

  my @exons;
  if (!defined( $self->{'_gene_sets'} ) ) {
    $self->warning("The gene array you try to retrieve exons for is empty");
    @exons = ();
  }

  foreach my $set_name (keys %{$self->{'_gene_sets'}}) {
    foreach my $gene (@{ $self->{'_gene_sets'}{$set_name} }) {  
      push @exons, @{$gene->get_all_Exons};
    }
  }

  return \@exons;
}

############################################################

sub strand{
  my $self = shift;

  if (!defined($self->{_cached_strand})) {
    my @genes = $self->get_Genes;
    unless (@genes){
      $self->warning("cannot retrieve the strand in a cluster with no genes");
    }
    my $strand;
    foreach my $gene (@genes){
      if (ref($gene)=~m/Gene/) {
        foreach my $transcript (@{$gene->get_all_Transcripts}){
          unless (defined($strand)) {
            $strand = $transcript->start_Exon->strand;
            next;
          }
          if ( $transcript->start_Exon->strand != $strand ){
           throw("You have a cluster with genes on opposite strands");
          }
        }
      }else {
       # prediction transcript
        unless (defined($strand)) {
         $strand = $gene->start_Exon->strand;
        }
        if ( $gene->start_Exon->strand != $strand ){
           throw("You have a cluster with genes on opposite strands");
        }
      }
    }
    $self->{_cached_strand} = $strand;
  }
  return $self->{_cached_strand};
}

#########################################################################


=head2 get_Gene_Count()

  it returns the number of genes in the GeneCluster object

=cut

sub get_Gene_Count {
  my $self = shift @_;

  my @genes = $self->get_Genes;
  return scalar(@genes);
}

#########################################################################


sub gene_Types {
  my ($self, $set_name, $types) = @_;
  $self->{'_types_sets'}{$set_name} = $types;
  

  return $self->{'_types_sets'}{$set_name};
}

#########################################################################

=head2 get_Genes_of_Type()

  We can get the genes in each cluster of one type. 
  We pass one string identifying the genetype.
  The difference with get_Genes_by_Type is that that as an arrayref as argument.

=cut
  
sub get_Genes_of_Type() {
  my ($self,$type) = @_;

  unless ($type){
    throw( "must provide a type");
  }

  my @genes = $self->get_Genes;  # this should give them in order, but we check anyway
  my @selected_genes;
  push ( @selected_genes, grep { $_->type eq $type } @genes );
  return @selected_genes;
}


sub get_Genes_by_Set() {
  my ($self,$set) = @_;

  unless ($set){
    throw( "must provide a set");
  }

  my @selected_genes;
  #for (keys %{ $self->{_gene_sets} } ) { 
  #  print " i know the following sets : $_\n" ;  
  #}
  if (!defined($self->{'_gene_sets'}{$set})) {
    # throw("No genes of set name $set");
    warning("No genes of set name $set in cluster");
  }else{
    push @selected_genes, @{$self->{'_gene_sets'}{$set}};
  }
  return @selected_genes;
}
#########################################################################

=head2 get_Genes_by_Type()

  We can get the genes in each cluster of a given type. 
  We pass an arrayref containing the types we want to retrieve.

=cut
  
sub get_Genes_by_Type() {
  my ($self,$types) = @_;
  unless ($types){
    throw( "must provide a type");
  }
  my @genes = $self->get_Genes;  # this should give them in order, but we check anyway
  my @selected_genes;
  foreach my $type ( @{ $types } ){
    push ( @selected_genes, grep { $_->type eq $type } @genes );
  }
  return @selected_genes;
}

#########################################################################



=head2 to_String()

  it returns a string containing the information about the genes contained in the
  GeneCluster object

=cut

sub to_String {
  my $self = shift @_;
  my $data='';
  foreach my $gene ( $self->get_Genes ){
    my @exons = @{ $gene->get_all_Exons };
     
    $data .= sprintf "Id: %-16s"      , $gene->stable_id;
    $data .= sprintf "Contig: %-20s"  , $exons[0]->contig->id;
    $data .= sprintf "Exons: %-3d"    , scalar(@exons);
    $data .= sprintf "Start: %-9d"    , $self->_get_start($gene);
    $data .= sprintf "End: %-9d"      , $self->_get_end  ($gene);
    $data .= sprintf "Strand: %-2d\n" , $exons[0]->strand;
  }
  return $data;
}

#########################################################################

=head2 _get_start()

 function to get the start position of a gene - it reads the gene object and it returns
 the start position of the first exon

=cut

sub _get_start {
  my ($self,$gene) = @_;
  my @exons = @{ $gene->get_all_Exons };
  my $st;
  
  if ($exons[0]->strand == 1) {
    @exons = sort {$a->start <=> $b->start} @exons;
    $st = $exons[0]->start;
  } else {
    @exons = sort {$b->start <=> $a->start} @exons; # they're read in opposite direction (from right to left)
    $st = $exons[0]->end;                           # the start is the end coordinate of the right-most exon
  }                                                 # which is here the first of the list of sorted @exons
  return $st;
}

#########################################################################

=head2 _get_end()

 function to get the end position of a gene - it reads the gene object and it returns
 the end position of the last exon

=cut

sub _get_end {
  my ($self,$gene) = @_;
  my @exons = @{ $gene->get_all_Exons };
  my $end;
  
  if ($exons[0]->strand == 1) {
    @exons = sort {$a->start <=> $b->start} @exons;
    $end = $exons[$#exons]->end;
  } else {
    @exons = sort {$b->start <=> $a->start} @exons; # they're read in opposite direction (from right to left)
    $end = $exons[$#exons]->start;                  # the end is the start coordinate of the left-most exon  
  }                                                 # which is here the last of the list @exons
  return $end;
}


=head2 _translateable_exon_length()

 internal function that returns the length of the translateable exons 

=cut

sub _translateable_exon_length {
  my ($trans)= @_;
  my @exons = $trans->translateable_exons;
  my $length = 0;
  foreach my $ex (@exons) {
    $length += $ex->length;
  }
  return $length;
}

# method to get the start of the cluster, which we take to be the left_most exon_coordinate 
# i.e. the start coordinate of the first exon ordered as { $a->start <=> $b->start }, regardless of the strand

sub start{
  my ($self) = @_;

  if (!defined($self->{_cached_end})) {
    my $start;

    foreach my $gene ($self->get_Genes) {
      my $this_start = $gene->start;
      unless ( $start ){
        $start = $this_start;
      }
      if ( $this_start < $start ){
        $start = $this_start;
      }
    }
    $self->{_cached_start} = $start;
  }
  return $self->{_cached_start};
}
      
# method to get the end of the cluster, which we take to be the right_most exon_coordinate
# this being the end coordinate of the first exon ordered as { $b->end <=> $a->end }, regardless of the strand

sub end{
  my ($self) = @_;

  if (!defined($self->{_cached_end})) {
    my $end;

    foreach my $gene ($self->get_Genes) {
      my $this_end = $gene->end;
      unless ( $end ){
        $end = $this_end;
      }
      if ( $this_end > $end ){
        $end = $this_end;
      }
    }
    $self->{_cached_end} = $end;
  }
  return $self->{_cached_end};
}



=head2 get_exon_clustering_from_gene_cluster 

   Name      : $self->get_exon_clustering_from_gene_cluster()
   Arg[0]    : Bio::EnsEMBL::Analysis::Tools::Algorithms::GeneCluster;
   Function  : gets a GeneCluster and converts it by building a TranscriptCluster, than 
               clusters the exons of all Transcripts and returns an array-ref to 
               Bio::EnsEMBL::ExonCluster-objects 
   Returnval :  Aref of  Bio::EnsEMBL::Analysis::Tools::Algorithms::ExonCluster objects

=cut

sub get_exon_clustering_from_gene_cluster {
  my ($self) = @_ ;

  my @clg  = sort {$a->start <=> $b->start} $self->get_Genes ;

  # building Transcript-Cluster 
  #my $tc = genes_to_Transcript_Cluster(\@clg);
  my $tc = $self->get_TranscriptCluster ; 

  my @exon_clusters = $tc->get_ExonCluster() ; 

  if ($tc->strand eq '1') {
    @exon_clusters = sort { $a->start <=> $b->start } @exon_clusters ;
  } else {
    @exon_clusters = sort { $b->start <=> $a->start } @exon_clusters ;
  }
  return \@exon_clusters ;
}



=head2 get_TranscriptCluster 

   Name      : $self->get_TranscriptCluster() 
   Arg[0]    : Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster 
   Function  : gets a GeneCluster and converts it by building a TranscriptCluster, than 
               clusters the exons of all Transcripts and returns an array-ref to 
               Bio::EnsEMBL::ExonCluster-objects 
   Returnval :  Aref of  Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster objects

=cut

sub get_TranscriptCluster {
  my ($self) = @_;
  #my ($genes_or_predTrans) = @_;

  my $tc = Bio::EnsEMBL::Analysis::Tools::Algorithms::TranscriptCluster->new() ;
  print "building new TranscriptCluster\n" ;
    
  foreach my $gene ( $self->get_Genes ) { 
  #foreach my $gene (@$genes_or_predTrans) {

    if( ref($gene)=~m/Gene/){
      # is a Bio::EnsEMBL::Gene 
      foreach my $trans (@{$gene->get_all_Transcripts}) {
        if ($gene->strand ne $trans->strand ) {
          throw("Weird - gene is on other strand than transcript\n") ;
        }
        for (@{ $trans->get_all_Exons} ) {
           if ($_->strand ne $trans->strand ) {
             print $trans->biotype . " " . $trans->seq_region_start . " "
              . $trans->seq_region_end . " " . $trans->seq_region_strand ."\n" ;
             print $_->biotype . " " . $_->seq_region_start . " "  . $_->seq_region_end . " " .$_->seq_region_strand . "\n";
             throw("Weird - exon is on other strand than transcript\n") ;
           }
        }
        # assure that transcript has same biotype as gene 
        $trans->biotype($gene->biotype) ;
        $trans->sort;
        #print "Adding transcript " . $trans->stable_id . "\n";
        $tc->put_Transcripts($trans);
        $tc->register_biotype($gene->biotype) ;
      }
    } else {
      # is not a Bio::EnsEMBL::Gene  
      warning("Not having a Bio::EnsEMBL::Gene-object : clustering $gene\n") ;
      $gene->sort ;
      $tc->put_Transcripts($gene) ;
    }
  }
  return $tc;
}


1;
