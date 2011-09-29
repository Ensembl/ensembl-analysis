# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::BlastPep
#
# Copyright (c) 2004 Ensembl
# Copyright (c) 2007 Wormbase
#

=head1 NAME

  Bio::EnsEMBL::Analysis::RunnableDB::BlastPep

=head1 SYNOPSIS

  my $blast = Bio::EnsEMBL::Analysis::RunnableDB::BlastPep->
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
 BlastPep

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Analysis::RunnableDB::BlastPep;

use Bio::EnsEMBL::Analysis::RunnableDB::Blast;
use Bio::EnsEMBL::Analysis::Runnable::BlastPep;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Analysis::Config::Blast;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::Blast);

=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastPep
  Function  : fetch sequence and prediction transcripts of database, 
  read config files instantiate the filter, parser and finally the blast 
  runnables
  Returntype: none
  Exceptions: none
  Example   : 

=cut

sub fetch_input {
  my ($self) = @_;

  my $t=$self->fetch_translation();
  my %blast = %{$self->BLAST_PARAMS};
  my $parser = $self->make_parser;
  my $filter;
  if($self->BLAST_FILTER){
    $filter = $self->make_filter;
  }
 
  my $bio_pep= Bio::PrimarySeq->new(-seq => $t->seq,-id => $t->dbID);

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::BlastPep->new(
          -transcript => $t,
          -query => $bio_pep,
          -program => $self->analysis->program_file,
          -parser => $parser,
          -filter => $filter,
          -database => $self->analysis->db_file,
          -analysis => $self->analysis,
          %blast,
         );
   $self->runnable($runnable);
}

=head2 fetch_translation

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::BlastPep
  Function  : fetches translations from the database and puts them as query sequence
  read config files instantiate the filter, parser and finally the blast 
  runnables
  Returntype: Bio::EnsEMBL::Translation
  Exceptions: none
  Example   : 

=cut

sub fetch_translation {
	my $self= shift;
	my $id=$self->input_id;
	my $db=$self->db;
	my $ta=$db->get_TranslationAdaptor;
	$self->{'query'}=$ta->fetch_by_dbID($id);
	return $self->{'query'};
}

=head2 write_output

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::WormBlast
  Function  : get appropriate adaptors and write output to database
  after validating and attaching sequence and analysis objects
  Returntype: undef
  Exceptions: throws if the store fails / or segfaults
  Example   : 

=cut


sub write_output {
        
  my ($self) = @_;
  my $protein_fa = $self->db->get_ProteinFeatureAdaptor;
  foreach my $f(@{$self->output}){
    $f->analysis($self->analysis);
    if($f->isa('Bio::EnsEMBL::ProteinFeature')){
            eval{ $protein_fa->store($f)};
            throw("Blast:store failed failed to write ".$f." to the database $@") if($@);
    }
  }
  return ;
}

1;
