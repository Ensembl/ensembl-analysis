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

Bio::EnsEMBL::Analysis::RunnableDB::HaplotypeMapper - 

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

=head1 METHODS


=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

# Let the code begin...

package Bio::EnsEMBL::Analysis::RunnableDB::HaplotypeMapper;

use warnings ;
use vars qw(@ISA);
use strict;

# Object preamble
use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Runnable::HaplotypeMapper;
use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);


use Bio::EnsEMBL::Analysis::Config::HaplotypeProjection    qw (
                                                             
                                                              );


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
  my($self,@mappings) = @_;
  
  #my $db = $self->get_dbadaptor("REFERENCE_DB") ;
  @mappings = $self->output;

  print "REACH TEST 2 position\n";
  print "now I have mappings: ",@mappings,"\n";
  my $genebuilders = $self->get_genebuilders;
  
  foreach my $target (keys %{ $genebuilders } ) {
    foreach my $query (keys %{$genebuilders->{$target}}){
      
      $genebuilders->{$target}->{$query}->load_tables(@mappings);
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

    my $hap_slice = $ref_db->get_SliceAdaptor->fetch_by_region($input_id[0],$input_id[2],$input_id[3],$input_id[4],1,$input_id[1]);
    my $slice = $ref_db->get_SliceAdaptor->fetch_by_region($input_id[5],$input_id[7],$input_id[8],$input_id[9],1,$input_id[6]);

    #$self->fetch_sequence();

    print "HAP_slice: ",$hap_slice,"\n";
    print "REF_slice: ",$slice,"\n";
  
    $self->query($hap_slice);
    $self->target($slice);

    print "QUERY: ",$self->query->seq_region_name,"\n";
    print "TARGET: ",$self->target->seq_region_name,"\n";

    my $genebuilder = new Bio::EnsEMBL::Analysis::Runnable::HaplotypeMapper
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
  
  my @mapping;

  # get a hash, with keys = contig/slice and value = genebuilder object
  my $genebuilders = $self->get_genebuilders;
  
  #my @genes;
  foreach my $target (keys %{ $genebuilders } ) {
    foreach my $query (keys %{$genebuilders->{$target}}){
        
      $genebuilders->{$target}->{$query}->create_alignment;

      $genebuilders->{$target}->{$query}->filter_alignment;

      @mapping = $genebuilders->{$target}->{$query}->make_map_regions;

      print "I got mappings: ",@mapping,"\n";
    }
  }
  
  $self->output( @mapping );
}

############################################################

# override the evil RunnableDB output method:

sub output{
    my ($self, @output ) = @_;
    unless ( $self->{_output} ){
	$self->{_output} = [];
    }
    if (@output){
	push( @{$self->{_output}}, @output );
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
