#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#

=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::HaplotypeProjection

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Analysis::RunnableDB::HaplotypeProjection->new(
								    -db        => $db,
								    -input_id  => $id,
								    );
    $obj->fetch_input
    $obj->run

    my @newfeatures = $obj->output;


=head1 DESCRIPTION

This method is used to get project the genes annotated in the reference chromosome into the Haplotype regions.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Analysis::RunnableDB::HaplotypeProjection;

use vars qw(@ISA);
use strict;

# Object preamble
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Runnable::HaplotypeProjection;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);


#use Bio::EnsEMBL::Analysis::Config::HaplotypeProjection    qw (
                                                             
#                                                              );


@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);

############################################################

=head2 new

    Usage   :   $self->new(-DBOBJ       => $db,
                           -INPUT_ID    => $id,
                           -SEQFETCHER  => $sf,
                           -ANALYSIS    => $analysis,
                          );

                           
    Function:   creates a Bio::EnsEMBL::Analysis::RunnableDB::HaplotypeProjection object
    Returns :   A Bio::EnsEMBL::Analysis::RunnableDB::HaplotypeProjection object
    Args    :   -dbobj:      A Bio::EnsEMBL::DBSQL::DBAdaptor (required), 
                -input_id:   Hap_pair input id (required), 
                -seqfetcher: A Sequence Fetcher Object,
                -analysis:   A Bio::EnsEMBL::Analysis (optional) 
                -extend:     determines the extension of the virtual contig
                             note: not implemented yet!
                -golden_path: determines the name of the golden path to use
=cut

sub new {
    my ($class,@args) = @_;

    my $self = $class->SUPER::new(@args);    
           
    return $self;
}

############################################################

sub input_id {
    my ($self,$arg) = @_;
    
    if (defined($arg)) {
	$self->{_input_id} = $arg;
    }
    
    return $self->{_input_id};
}

############################################################

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   array of exons (with start and end)
    Args    :   none

=cut
    
    
sub write_output {
  my($self,@genes) = @_;
  
  my $db = $self->get_dbadaptor("REFERENCE_DB") ;
  # sort out analysis
  
  my $analysis = $self->analysis;
  unless ($analysis){
    $self->throw("an analysis logic name must be defined in the command line");
  }
  
  my $gene_adaptor = $db->get_GeneAdaptor;
  
    
  @genes = $self->output;
    
  foreach my $gene (@genes) { 
    my $before = scalar(@{$gene->get_all_Transcripts});
    unless (scalar(@{$gene->get_all_Transcripts})){
      warning("GENE DOES NOT HAVE TRANSCRIPTS:\n");
      next;
    }
    my %trans_types;
    $gene->analysis($analysis);
    $gene->{'stable_id'} = '';
    # poke the caches
    my %s_pfhash;
    foreach my $tran (@{$gene->get_all_Transcripts}) {
      $tran->{'stable_id'} = '';
      my @tsf = @{$tran->get_all_supporting_features};
      my @exons= @{$tran->get_all_Exons};
      my $tln = $tran->translation;
      if ($tln){
        $tln->{'stable_id'} = '';
      }
      
      foreach my $exon (@exons) {
        if ($tran->seq_region_name != $exon->seq_region_name){
          print "NO EXON WAS NOT TRANSFORMED BEFORE STORAGE\n";
          
        }else{
          print "TRANSFORMED transcript: ",$tran->seq_region_name , " exon: ", $exon->seq_region_name,"\n";
        }
        $exon->{'stable_id'} = '';
        my @esf = @{$exon->get_all_supporting_features};
      }
    }
    
    # store
   # eval {
    if ($before > scalar(@{$gene->get_all_Transcripts})){
      print "MISSING TRANSCRIPTS IN WRITTING\n";
    }
      $gene_adaptor->store($gene);
      #print STDERR "wrote gene " . $gene->dbID . " to database ".$gene->adaptor->db->dbname."\n";
   # }; 
   # if( $@ ) {
   #   warning("UNABLE TO WRITE GENE:\n$@");
   # }
  }
  
  my $genebuilders = $self->get_genebuilders;
  
  foreach my $target (keys %{ $genebuilders } ) {
    foreach my $query (keys %{$genebuilders->{$target}}){
      $genebuilders->{$target}->{$query}->clean_tables;
    }   
  } 
  
  return 1;   
}

############################################################

=head2 fetch_input

    Function:   It fetches the slice or contig according to the input_id, 
                and it defines the database where the
                previous annotations are stored and create a Bio::EnsEMBL::Pipeline::GeneBuilder
                object for that genomic, input_id and db
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
    my( $self) = @_;
    
    $self->throw("No input id") unless defined($self->input_id);

    
    my $discarded_db = $self->get_dbadaptor("DISCARDED_DB");

    print "DISCARDED GENE DB: ", $discarded_db->dbname,"\n";

    # database where the genebuild produced genes are
    
    my $ref_db = $self->get_dbadaptor("REFERENCE_DB");
    print "ENSEMBL DB : ",  $ref_db->dbname,"\n";

    print $self->input_id,"\n";

    my @input_id = split(/:/,$self->input_id);

    my $hap_slice = $ref_db->get_SliceAdaptor->fetch_by_region($input_id[0],$input_id[2],$input_id[3],$input_id[4],1,$input_id[2]);
    my $slice = $ref_db->get_SliceAdaptor->fetch_by_region($input_id[5],$input_id[7],$input_id[8],$input_id[9],1,$input_id[6]);

    #$self->fetch_sequence();

    print "HAP_slice: ",$hap_slice,"\n";
    print "REF_slice: ",$slice,"\n";
  
    $self->query($hap_slice);
    $self->target($slice);

    print "QUERY: ",$self->query->seq_region_name,"\n";
    print "TARGET: ",$self->target->seq_region_name,"\n";

    my $genebuilder = new Bio::EnsEMBL::Analysis::Runnable::HaplotypeProjection
      (
       '-hap_slice' => $self->query,
       '-slice'   => $self->target,
       '-input_id' => $self->input_id,
      );
    $genebuilder->discarded_db($discarded_db);
    $genebuilder->ensembl_db($ref_db);
     
    # store the object and the piece of genomic where it will run
    $self->addgenebuilder($genebuilder,$self->target,$self->query);
    
}

############################################################

sub addgenebuilder {
    my ($self,$arg,$target,$query) = @_;
    
    if (defined($arg) && defined($target) && defined($query)) {
	$self->{_genebuilder}{$target->id}{$query->id} = $arg;
    } 
    else {
	$self->throw("Wrong number of inputs [$arg,$target,$query]\n");
    }
}

############################################################

sub get_genebuilders {
    my ($self) = @_;
    
    return $self->{_genebuilder};
}

############################################################
	
sub run {
  my ($self) = @_;
  
  my @genes;

  # get a hash, with keys = contig/slice and value = genebuilder object
  my $genebuilders = $self->get_genebuilders;
  
  foreach my $target (keys %{ $genebuilders } ) {
    foreach my $query (keys %{$genebuilders->{$target}}){

      @genes = $genebuilders->{$target}->{$query}->project_genes;

    }
  }
  
  $self->output( @genes );
}

############################################################

# override the evil RunnableDB output method:

sub output{
    my ($self, @genes ) = @_;
    unless ( $self->{_output} ){
	$self->{_output} = [];
    }
    if (@genes){
	push( @{$self->{_output}}, @genes );
    }
    return @{$self->{_output}};
}

############################################################

sub target {
  my ($self,$slice) = @_;
  
  if (defined($slice)) {
    $self->{_target} = $slice;
  }
  return $self->{_target};
}


1;
