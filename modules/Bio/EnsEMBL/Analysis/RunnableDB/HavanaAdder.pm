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

  Questions may also be sent to the Ensembl help desk
  at <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

  Bio::EnsEMBL::Analysis::RunnableDB::HavanaAdder

=head1 SYNOPSIS

  my $obj =
    Bio::EnsEMBL::Analysis::RunnableDB::HavanaAdder->new( -db       => $db,
                                                          -input_id => $id,
    );
  $obj->fetch_input;
  $obj->run;

  my @newfeatures = $obj->output;

=cut

# Let the code begin...

package Bio::EnsEMBL::Analysis::RunnableDB::HavanaAdder;

use warnings ;
use vars qw(@ISA);
use strict;

# Object preamble
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Runnable::HavanaAdder;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);

use Bio::EnsEMBL::Analysis::Config::HavanaAdder qw (
  ENSEMBL_INPUT_CODING_TYPE
  HAVANA_GENE_OUTPUT_BIOTYPE
  MERGED_GENE_OUTPUT_BIOTYPE
  ENSEMBL_GENE_OUTPUT_BIOTYPE
  MERGED_TRANSCRIPT_OUTPUT_TYPE
);


@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild);

############################################################

=head2 new

  Usage      :   $self->new(-DBOBJ       => $db,
                          -INPUT_ID    => $id,
                          -SEQFETCHER  => $sf,
                          -ANALYSIS    => $analysis,
                        );
  Description: Creates a Bio::EnsEMBL::Analysis::RunnableDB::HavanaAdder object
  Returns    : A Bio::EnsEMBL::Analysis::RunnableDB::HavanaAdder object
  Args       : -dbobj:      A Bio::EnsEMBL::DBSQL::DBAdaptor (required),
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

  Usage      : $self->write_output
  Description: Writes output data to db
  Returns    : Nothing
  Args       : None

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

  Description: It fetches the slice or contig according
               to the input_id, and it defines the database where
               the previous annotations are stored and create a
               Bio::EnsEMBL::Analysis::Runnable::HavanaAdder object
               for that genomic, input_id and db
  Returns    : Nothing
  Args       : None

=cut

sub fetch_input {
  my ($self) = @_;

  $self->throw("No input id") unless defined( $self->input_id );

  $self->fetch_sequence();

  # database where the ensembl genebuild genes are located
  my $ensembl_db = $self->get_dbadaptor("PSEUDO_DB");

  print "ENSEMBL DB: ", $ensembl_db->dbname, "\n";

  # database with the Havana Vega genes to import
  my $havana_db = $self->get_dbadaptor("HAVANA_DB");

  print "HAVANA DB: ", $havana_db->dbname, "\n";

  # Database with the CCDS models
  my $ccds_db = $self->get_dbadaptor("CCDS_DB");

  print "CCDS DB: ";
  if ( defined($ccds_db) ) {
    printf( "%s\n", $ccds_db->dbname() );
  }
  else {
    print("(not defined)\n");
  }

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
  $genebuilder->ensembl_db($ensembl_db);
  $genebuilder->havana_db($havana_db);
  $genebuilder->ccds_db($ccds_db);

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
