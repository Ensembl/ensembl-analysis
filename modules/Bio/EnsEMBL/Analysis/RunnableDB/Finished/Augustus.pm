=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
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

Bio::EnsEMBL::Analysis::RunnableDB::Finished::Augustus - 

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 METHODS

=cut

# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::Finished::Augustus
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

=head1 SYNOPSIS

  my $runnabledb = Bio::EnsEMBL::Analysis::RunnableDB::Finished::Augustus->
  new(
      -input_id => 'contig::AL805961.22.1.166258:1:166258:1',
      -db => $db,
      -analysis => $analysis,
     );
  $runnabledb->fetch_input;
  $runnabledb->run;
  $runnabledb->write_output;


=head1 DESCRIPTION

fetches sequence data from database an instantiates and runs the
augustus runnable, this inherits from the Genscan runnableDB an as such doesnt
implement much itself

=head1 CONTACT

Post questions to : anacode-people@sanger.ac.uk

=cut


package Bio::EnsEMBL::Analysis::RunnableDB::Finished::Augustus;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::RunnableDB::Genscan;
use Bio::EnsEMBL::Analysis::Runnable::Finished::Augustus;
use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::Genscan);



=head2 runnable_path

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Finished::Augustus
  Function  : return the runnable path
  Returntype: string
  Exceptions:
  Example   : my $runnable = $self->runnable_path->new
                               (
                                -query    => $self->query,
                                -program  => $self->analysis->program_file,
                                -analysis => $self->analysis,
                                %parameters,
                               );

=cut


sub runnable_path{
  my ($self);
  return "Bio::EnsEMBL::Analysis::Runnable::Finished::Augustus";
}

sub write_output{
  my ($self) = @_;
  my $adaptor = $self->db->get_PredictionTranscriptAdaptor;
  my $dbh = $self->db->dbc->db_handle;
  my @output = @{$self->output};
  my $ff = $self->feature_factory;

  $dbh->begin_work;

  eval {
	  foreach my $pt(@output){
	    $pt->analysis($self->analysis);
	    $pt->slice($self->query) if(!$pt->slice);
	    print STDERR "Validate transcript ".$pt->seqname."\n";
	    # dismiss transcript with invalid translation
	    eval {
	    	$ff->validate_prediction_transcript($pt, 1);
	    };
	    if($@){ warning($@); next; }
	    $adaptor->store($pt);
	  }
	  $dbh->commit;
  };
  if ($@) {
      $dbh->rollback;
	  throw("UNABLE TO WRITE PREDICTION TRANSCRIPTS IN DATABASE\n[$@]\n");
  }
}



1;
