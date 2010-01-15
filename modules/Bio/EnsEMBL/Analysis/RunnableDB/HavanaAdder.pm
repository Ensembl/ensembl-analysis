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
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Runnable::HavanaAdder;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Analysis::Config::HavanaAdder            qw (
                                                               GB_ENSEMBL_INPUT_GENETYPE
                                                               HAVANA_GENE_OUTPUT_BIOTYPE
                                                               MERGED_GENE_OUTPUT_BIOTYPE
                                                               ENSEMBL_GENE_OUTPUT_BIOTYPE
                                                               MERGED_TRANSCRIPT_OUTPUT_TYPE
                                                              );


@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);

############################################################

=head2 new

    Usage   :   $self->new(-DBOBJ       => $db,
                           -INPUT_ID    => $id,
                           -SEQFETCHER  => $sf,
                           -ANALYSIS    => $analysis,
                          );

    Function:   creates a Bio::EnsEMBL::Analysis::RunnableDB::HavanaAdder object
    Returns :   A Bio::EnsEMBL::Analysis::RunnableDB::HavanaAdder object
    Args    :   -dbobj:      A Bio::EnsEMBL::DBSQL::DBAdaptor (required), 
                -input_id:   Contig input id (required), 
                -seqfetcher: A Sequence Fetcher Object,
                -analysis:   A Bio::EnsEMBL::Analysis (optional) 

=cut

sub new {
  my ( $class, @args ) = @_;

  my $self = $class->SUPER::new(@args);

  return $self;
}

############################################################

sub input_id {
  my ( $self, $arg ) = @_;

  if ( defined($arg) ) {
    $self->{_input_id} = $arg;
  }

  return $self->{_input_id};
}

############################################################

=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   Writes output data to db
    Returns :   nothing
    Args    :   none

=cut

sub write_output {

  my ( $self, @genes ) = @_;

  print "Starting the write output process now \n";

  # write genes out to a different database from the one we read
  # genewise genes from.
  my $db = $self->get_dbadaptor("GENEBUILD_DB");

  #print "WRITE DB IS:",%$db->dbc->dbname,"\n";

  # sort out analysis
  my $analysis = $self->analysis;
  unless ($analysis) {
    $self->throw(
                "an analysis logic name must be defined in the command line");
  }

  #  my %contighash;
  my $gene_adaptor = $db->get_GeneAdaptor;

  # this now assummes that we are building on a single VC.
  my $genebuilders = $self->get_genebuilders;

  foreach my $genesbuilt ( keys %$genebuilders ) {
    # my $vc = $genebuilders->{$contig}->query;

    @genes = @{ $genebuilders->{$genesbuilt}->final_genes } ; 
    print "I have ", scalar(@genes), " genes\n";

    return unless ( $#genes >= 0 );

    foreach my $gene (@genes) {
      $gene->analysis($analysis);
      # store
      eval {
        $gene_adaptor->store($gene);
        print STDERR "wrote gene " . $gene->dbID . " to database\n";
      };
      if ($@) {
        warning("UNABLE TO WRITE GENE:\n$@");
      }
    }
  }
  return 1;
} ## end sub write_output

############################################################

=head2 fetch_input

    Function:   It fetches the slice or contig according to the input_id, 
                and it defines the database where the
                previous annotations are stored and create a Bio::EnsEMBL::Analysis::Runnable::HavanaAdder
                object for that genomic, input_id and db
    Returns :   nothing
    Args    :   none

=cut

sub fetch_input {
  my ($self) = @_;

  $self->throw("No input id") unless defined( $self->input_id );

  $self->fetch_sequence();

  # Database with a number of genes which where identified as wrong
  # structures so if we build something similar we want to remove it
  my $discarded_db = $self->get_dbadaptor("DISCARDED_DB");

  print "DISCARDED GENE DB: ", $discarded_db->dbname, "\n";

  # database where the ensembl genebuild genes are located
  my $ensembl_db = $self->get_dbadaptor("PSEUDO_DB");

  print "ENSEMBL DB : ", $ensembl_db->dbname, "\n";

  # database with the Havana Vega genes to import
  my $havana_db = $self->get_dbadaptor("HAVANA_DB");

  print "HAVANA DB : ", $havana_db->dbname, "\n";

  # Database that contains the DNA sequence
  my $ref_db = $self->get_dbadaptor("REFERENCE_DB");

  print $self->input_id, "\n";

  my $slice = $ref_db->get_SliceAdaptor->fetch_by_name( $self->input_id );

  print $slice, "\n";

  $self->query($slice);

  print "QUERY: ", $self->query->seq_region_name, "\n";
  my $genebuilder =
    new Bio::EnsEMBL::Analysis::Runnable::HavanaAdder(
                                               '-slice'    => $self->query,
                                               '-input_id' => $self->input_id,
    );
  $genebuilder->discarded_db($discarded_db);
  $genebuilder->ensembl_db($ensembl_db);
  $genebuilder->havana_db($havana_db);

  # store the object and the piece of genomic where it will run
  $self->addgenebuilder( $genebuilder, $self->query );

  print "I finished fetching the database adaptors\n";

} ## end sub fetch_input

############################################################

sub addgenebuilder {
  my ( $self, $arg, $contig ) = @_;

  if ( defined($arg) && defined($contig) ) {
    $self->{_genebuilder}{ $contig->id } = $arg;
  } else {
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

  print "Now running the analysis\n";

  # get a hash, with keys = contig/slice and value = genebuilder object
  my $genebuilders = $self->get_genebuilders;

  print "Getting Gene adaptors again\n";
  my @genes;
  foreach my $region ( keys %{$genebuilders} ) {
    print "Starting to get some genes\n";
    my $query = $genebuilders->{$region}->query;

    print "GeneBuilding for $region\n";

    $genebuilders->{$region}->build_Genes;

    print "Genes build now getting the final set\n";
    @genes = @{ $genebuilders->{$region}->final_genes } ; 
  }
  print "OK now I have my genes, just need to write them\n";
  $self->output(@genes);
}

############################################################

# override the evil RunnableDB output method:

sub output {
  my ( $self, @genes ) = @_;
  unless ( $self->{_output} ) {
    $self->{_output} = [];
  }
  if (@genes) {
    push( @{ $self->{_output} }, @genes );
  }
  #return @{$self->{_output}};
  return $self->{_output};
}

############################################################



1;
