#
# Cared for by EnsEMBL  <ensembl-dev@ebi.ac.uk>
#
# written by Michele Clamp
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::HavanaAdder

=head1 SYNOPSIS

    my $obj = Bio::EnsEMBL::Analysis::RunnableDB::HavanaAdder->new(
								    -db        => $db,
								    -input_id  => $id,
								    );
    $obj->fetch_input
    $obj->run

    my @newfeatures = $obj->output;


=head1 DESCRIPTION

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Analysis::RunnableDB::HavanaAdder;

use vars qw(@ISA);
use strict;

# Object preamble
use Bio::EnsEMBL::Analysis::Tools::TranscriptUtils;
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Runnable::HavanaAdder;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;

use Bio::EnsEMBL::Analysis::Config::GeneBuild::General     qw (
							       GB_INPUTID_REGEX
							      );
#use Bio::EnsEMBL::Pipeline::Config::GeneBuild::GeneBuilder qw (
#							       GB_VCONTIG
#							      );
use Bio::EnsEMBL::Analysis::Config::HavanaAdder            qw (
                                                               GB_GENE_OUTPUT_BIOTYPE
                                                              );
use Bio::EnsEMBL::Analysis::Config::GeneBuild::Databases   qw (DATABASES
                                                              );


@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);

############################################################

=head2 new

    Usage   :   $self->new(-DBOBJ       => $db,
                           -INPUT_ID    => $id,
                           -SEQFETCHER  => $sf,
                           -ANALYSIS    => $analysis,
                           -VCONTIG     => 1,
                          );

                           
    Function:   creates a Bio::EnsEMBL::Analysis::RunnableDB::HavanaAdder object
    Returns :   A Bio::EnsEMBL::Analysis::RunnableDB::HavanaAdder object
    Args    :   -dbobj:      A Bio::EnsEMBL::DBSQL::DBAdaptor (required), 
                -input_id:   Contig input id (required), 
                -seqfetcher: A Sequence Fetcher Object,
                -analysis:   A Bio::EnsEMBL::Analysis (optional) 
                -vcontig:    determines whether it is running on virtual contigs
                             or RawContigs
                -extend:     determines the extension of the virtual contig
                             note: not implemented yet!
                -golden_path: determines the name of the golden path to use
=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);    
           
#    my( $use_vcontig) = $self->_rearrange([qw(VCONTIG)], @args);
#       
#    if (! defined $use_vcontig) {
#      $use_vcontig = $GB_VCONTIG;
#    }  
    
#    $self->use_vcontig($use_vcontig);

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
  
  # write genes out to a different database from the one we read genewise genes from.
  my $db = $self->get_dbadaptor("GENEBUILD_DB") ;
  # sort out analysis
  
  my $analysis = $self->analysis;
  unless ($analysis){
    $self->throw("an analysis logic name must be defined in the command line");
  }
  
  my %contighash;
  my $gene_adaptor = $db->get_GeneAdaptor;
  
  # this now assummes that we are building on a single VC.
  my $genebuilders = $self->get_genebuilders;
    
  foreach my $contig ( keys %$genebuilders ){
    my $vc = $genebuilders->{$contig}->query;
      
    @genes = $genebuilders->{$contig}->final_genes;
    
    return unless ($#genes >= 0);
    my @newgenes;
    
    foreach my $gene (@genes) { 
      $gene->analysis($analysis);
      $gene->type($GB_GENE_OUTPUT_BIOTYPE);
      # poke the caches
      my %s_pfhash;
      foreach my $tran (@{$gene->get_all_Transcripts}) {
        #$tran->stable_id(undef);
        my @tsf = @{$tran->get_all_supporting_features};
               
        my @exons= @{$tran->get_all_Exons};
        my $tln = $tran->translation;
        $tln->{'stable_id'} = undef;
        
        foreach my $exon (@exons) {
          my @esf = @{$exon->get_all_supporting_features};
          #$exon->{'stable_id'} = undef;
        }
      }  
      # store
      eval {
        $gene_adaptor->store($gene);
        #print STDERR "wrote gene " . $gene->dbID . " to database ".
        #   $gene->adaptor->db->dbname."\n";
      }; 
      if( $@ ) {
        $self->warn("UNABLE TO WRITE GENE:\n$@");
      }
    }   
  }    
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
    
    $self->fetch_sequence();
    # database where the genebuild produced genes are
    my $ensembl_db = $self->get_dbadaptor("PSEUDO_DB") ;

    my $havana_db = $self->get_dbadaptor("HAVANA_DB") ;
    
    #print STDERR "reading genewise and combined genes from $GB_COMB_DBNAME : $GB_COMB_DBHOST\n";
    
    my $genebuilder = new Bio::EnsEMBL::Analysis::Runnable::HavanaAdder
      (
       '-slice'   => $self->query,
       '-input_id' => $self->input_id,
      );
    $genebuilder->ensembl_db($ensembl_db);
    $genebuilder->havana_db($havana_db);
    
    # store the object and the piece of genomic where it will run
    $self->addgenebuilder($genebuilder,$self->query);
    
}

############################################################

#sub use_vcontig {
#    my ($self,$arg) = @_;
#    
#    if (defined($arg)) {
#	$self->{_vcontig} = $arg;
#    }
#
#    return $self->{_vcontig};
#}

############################################################

sub addgenebuilder {
    my ($self,$arg,$contig) = @_;
    
    if (defined($arg) && defined($contig)) {
	$self->{_genebuilder}{$contig->id} = $arg;
    } 
    else {
	$self->throw("Wrong number of inputs [$arg,$contig]\n");
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
    
    # get a hash, with keys = contig/slice and value = genebuilder object
    my $genebuilders = $self->get_genebuilders;
    
    my @genes;
    foreach my $contig (keys %{ $genebuilders } ) {
      my $query = $genebuilders->{$contig}->query;
      
      #print(STDERR "GeneBuilding for $contig\n");
      
      $genebuilders->{$contig}->build_Genes;
      
      @genes = $genebuilders->{$contig}->final_genes;
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



1;
