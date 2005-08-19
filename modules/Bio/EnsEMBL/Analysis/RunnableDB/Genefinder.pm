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
