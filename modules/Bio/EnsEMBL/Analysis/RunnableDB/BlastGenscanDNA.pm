# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::BlastGenscanDNA
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

  Bio::EnsEMBL::Analysis::RunnableDB::Blast

=head1 SYNOPSIS

  my $blast = Bio::EnsEMBL::Analysis::RunnableDB::BlastGenscanDNA->
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
 BlastTranscriptDNA passing in prediction transcript

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::BlastGenscanDNA;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptDNA;
use Bio::EnsEMBL::Analysis::RunnableDB::Blast;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Analysis::Config::Blast;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::Blast);


=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastGenscanDNA 
  Function  : fetch sequence and prediction transcripts of database, 
  read config files instantiate the filter, parser and finally the blast 
  runnables
  Returntype: none
  Exceptions: none
  Example   : 

=cut



sub fetch_input{
  my ($self) = @_;
  $self->setup_hashes;
  my $slice = $self->fetch_sequence($self->input_id, $self->db);
  $self->query($slice);
  my %blast = %{$self->blast_hash};
  my $logic_name = $BLAST_AB_INITIO_LOGICNAME;
  $logic_name = 'Genscan' if(!$logic_name);
  my $pta = $self->db->get_PredictionTranscriptAdaptor;
  my $pts = $pta->fetch_all_by_Slice($self->query, $logic_name);
  my $parser = $self->make_parser;
  my $filter;
  if($self->filter_object){
    $filter = $self->make_filter;
  }
  foreach my $t(@$pts){
    my $runnable = Bio::EnsEMBL::Analysis::Runnable::BlastTranscriptDNA->
      new(
          -transcript => $t,
          -query => $self->query,
          -program => $self->analysis->program_file,
          -parser => $parser,
          -filter => $filter,
          -database => $self->analysis->db_file,
          -analysis => $self->analysis,
          %blast,
         );
    $self->runnable($runnable);
  }
}


