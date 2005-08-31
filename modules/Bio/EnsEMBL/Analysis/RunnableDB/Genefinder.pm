# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::Genefinder
#
# Copyright (c) 2005 Ensembl
#

=head1 NAME

=head1 SYNOPSIS

  my $runnabledb = Bio::EnsEMBL::Analysis::RunnableDB::Genefinder->
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
genefinder runnable


=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::RunnableDB::Genefinder;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::RunnableDB::Genscan;
use Bio::EnsEMBL::Analysis::Runnable::Genefinder;
use Bio::EnsEMBL::Analysis::Config::General;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB::Genscan);



=head2 runnable_path

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Genefinder
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
  my ($self) = @_;
  return "Bio::EnsEMBL::Analysis::Runnable::Genefinder";
}




=head2 standard_args

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Genefinder
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
         };
}

1;
