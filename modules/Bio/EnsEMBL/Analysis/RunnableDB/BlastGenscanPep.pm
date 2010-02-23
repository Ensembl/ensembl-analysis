# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::BlastGenscanPep
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

  Bio::EnsEMBL::Analysis::RunnableDB::Blast

=head1 SYNOPSIS

  my $blast = Bio::EnsEMBL::Analysis::RunnableDB::BlastGenscanPep->
  new(
      -analysis => $analysis,
      -db => $db,
      -input_id => 'contig::AL1347153.1.3517:1:3571:1'
     );
  $blast->fetch_input;
  $blast->run;
  my @output =@{$blast->output};

=head1 DESCRIPTION

 This module inherits from the Blast runnable and instantiates 
 BlastTranscriptPep passing in prediction transcript

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::BlastGenscanPep;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep;
use Bio::EnsEMBL::Analysis::RunnableDB::Blast;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Analysis::Config::Blast;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::Blast);

=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastGenscanPep
  Function  : fetch sequence and prediction transcripts of database, 
  read config files instantiate the filter, parser and finally the blast 
  runnables
  Returntype: none
  Exceptions: none
  Example   : 

=cut

sub fetch_input{
  my ($self) = @_;
  my $slice = $self->fetch_sequence($self->input_id, $self->db);
  $self->query($slice);
  my %blast = %{$self->BLAST_PARAMS};
  my $logic_name = $BLAST_AB_INITIO_LOGICNAME;
  $logic_name = 'Genscan' if(!$logic_name);
  my $pta = $self->db->get_PredictionTranscriptAdaptor;
  my $pts = $pta->fetch_all_by_Slice($self->query, $logic_name);
  my $parser = $self->make_parser;
  my $filter;
  if($self->BLAST_FILTER){
    $filter = $self->make_filter;
  }  

  # submit blast module to use via analysis_parameters column of analysis table 
  my $options_string ;  
  my %options = %{$self->PARSER_PARAMS};  

  if ( $options{-query_type}=~m/pep/ ) {  
    if ( $options{-database_type}=~m/pep/ ) {  
       $options_string = '-p blastp' ; 
    } elsif ( $options{-database_type}=~m/dna/ ) {  
       $options_string = '-p tblastn' ; 
    } 
  }   

  if ( $options{-query_type}=~m/dna/ ) {  
    if ( $options{-database_type}=~m/dna/ ) {  
       $options_string = '-p blastn' ; 
    }elsif ( $options{-database_type}=~m/pep/ ) {  
       $options_string = '-p blastx' ; 
    }
  }   


  foreach my $t(@$pts){
    my $runnable = Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptPep->
      new(
          -transcript => $t,
          -query => $self->query,
          -program => $self->analysis->program_file,
          -parser => $parser,
          -filter => $filter,
          -database => $self->analysis->db_file,
          -analysis => $self->analysis,   
          -options => $options_string, 
          %blast,
         );
    $self->runnable($runnable);
  }
}




