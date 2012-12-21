=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::BlastGenscanPep - 

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

=head1 METHODS

=cut

# $Source: /tmp/ENSCOPY-ENSEMBL-ANALYSIS/modules/Bio/EnsEMBL/Analysis/RunnableDB/BlastGenscanPep.pm,v $
# $Revision: 1.14 $
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
  my $pta = $self->db->get_PredictionTranscriptAdaptor;
  my $logic_names = $BLAST_AB_INITIO_LOGICNAME ;
  if ( !ref($logic_names) || scalar(@$logic_names) == 0 ) {
    $logic_names = ['Genscan'];
  }
  my @pts ;
  foreach my $logic_name (@$logic_names) {
    my $pt = $pta->fetch_all_by_Slice($self->query, $logic_name);
    push @pts, @$pt ;
  }
  my $parser = $self->make_parser;
  my $filter;
  if($self->BLAST_FILTER){
    $filter = $self->make_filter;
  }  

  # submit blast module to use via analysis_parameters column of analysis table 
  my $options_string ;  
  my %options = %{$self->PARSER_PARAMS};  

  if ( $blast{-type}=~m/ncbi/ ) { 
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
  }


  foreach my $t(@pts){
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




