# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::Genscan
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

=head1 SYNOPSIS

  my $runnabledb = Bio::EnsEMBL::Analysis::RunnableDB::Genscan->
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
genscan runnable


=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Genscan;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Runnable::Genscan;
use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info);
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB);


=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Genscan
  Function  : fetch data out of database and create runnable
  Returntype: 1
  Exceptions: none
  Example   : 

=cut



sub fetch_input{
  my ($self) = @_;
  my $slice = $self->fetch_sequence($self->input_id, $self->db, 
                                    $ANALYSIS_REPEAT_MASKING);
  $self->query($slice);
  my %args = %{$self->standard_args};
  my $runnable = $self->runnable_path->new
    (
     %args,
    );
  
  my $seq = $self->query->seq; 
  if ($seq =~ /[CATG]{3}/) {
     $self->input_is_void(0);
  } else {
     $self->input_is_void(1);
     warning("Need at least 3 nucleotides - maybe your sequence was fully repeatmasked ...");
  }
  
  $self->runnable($runnable);
  return 1;
}




=head2 write_output

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Genscan
  Function  : writes the prediction transcripts back to the database
  after validation
  Returntype: none
  Exceptions: 
  Example   : 

=cut



sub write_output{
  my ($self) = @_;
  my $adaptor = $self->db->get_PredictionTranscriptAdaptor;
  my @output = @{$self->output};
  my $ff = $self->feature_factory;
  foreach my $pt(@output){
    $pt->analysis($self->analysis);
    $pt->slice($self->query) if(!$pt->slice);
    $ff->validate_prediction_transcript($pt, 1);
    $adaptor->store($pt);
  }
}


sub runnable_path{
  my ($self);
  return "Bio::EnsEMBL::Analysis::Runnable::Genscan";
}



=head2 standard_args

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Genscan
  Function  : to return a hash which contains the 
  standard constructor args for the genscan runnable
  Returntype: hashref
  Exceptions: none
  Example   : 

=cut




sub standard_args{
  my ($self) = @_;
  my %parameters;
  if($self->parameters_hash){
    %parameters = %{$self->parameters_hash};
  }
  return {
          -query => $self->query,
          -program => $self->analysis->program_file,
          -analysis => $self->analysis,
          %parameters,
          -matrix => $self->analysis->db_file,
         };
}
1;
