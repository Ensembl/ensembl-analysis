package Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  my ($matrix) = rearrange(['MATRIX'], @args);
  $self->matrix($matrix);

  return $self;
}




#containters



=head2 matrix

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio
  Arg [2]   : string, matrix file
  Function  : container for path to matrix file, when passed a filename it
  will use the find_file method to check for its existance
  Returntype: string
  Exceptions: 
  Example   : 

=cut


sub matrix{
  my ($self, $matrix) = @_;
  if($matrix){
    my $file = $self->find_file($matrix);
    $self->{'matrix'} = $file;
  }
  return $self->{'matrix'};
}


=head2 peptides

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio
  Arg [2]   : string, peptide sequence from genscan results
  Function  : container for the peptides sequences from genscan results
  Returntype: arrayref
  Exceptions: 
  Example   : 

=cut



sub peptides{
  my ($self, $peptides) = @_;
  if($peptides){
    $self->{'peptides'} = $peptides;
  }
  return $self->{'peptides'};
}


=head2 exon_groups

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio
  Arg [2]   : string, group name
  Arg [3]   : Bio::EnsEMBL::PredictionExon
  Function  : stores exons in arrayrefs in a hash keyed on group names
  Returntype: hashref
  Exceptions: throws if passed an exon which isnt a 
  Bio::EnsEMBL::PredictionExon
  Example   : 

=cut


sub exon_groups{
  my ($self, $group, $exon) = @_;
  if(!$self->{'exon_groups'}){
    $self->{'exon_groups'} = {};
  }
  if($group && $exon){   
    if(!$exon->isa('Bio::EnsEMBL::PredictionExon')){
      throw("must be passed an exon object ".
            "Bio::EnsEMBL::PredictionExon not an $exon ".
            "BaseAbInitio:exon_groups");
    }
    if(!$self->{'exon_groups'}->{$group}){
      $self->{'exon_groups'}->{$group} = [];
    }
    push(@{$self->{'exon_groups'}->{$group}}, $exon);
  }
  return $self->{'exon_groups'};
}


#utility methods


=head2 run_analysis

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio
  Arg [2]   : string, program name
  Function  : create and open a commandline for the program trf
  Returntype: none
  Exceptions: throws if the program in not executable or the system
  command fails to execute
  Example   : 

=cut


sub run_analysis{
  my ($self, $program) = @_;
  if(!$program){
    $program = $self->program;
  }
  throw($program." is not executable BaseAbInitio::run_analysis ") 
    unless($program && -x $program);
  
  my $command = $program." ".$self->matrix." ".$self->queryfile." > ".
    $self->resultsfile;
  print "Running analysis ".$command."\n";
  system($command) == 0 or throw("FAILED to run ".$command);
}

=head2 create_transcripts

  Arg [1]   : Bio::EnsEMBL::Analysis::Runnable::BaseAbInitio
  Function  : use the groups of exons to create prediction transcripts
  Returntype: none
  Exceptions: none
  Example   : 

=cut



sub create_transcripts{
  my ($self) = @_;
  my @transcripts;
  my $ff = $self->feature_factory;
  my %exon_groups = %{$self->exon_groups};
  foreach my $group(keys(%exon_groups)){
    my @exons = @{$exon_groups{$group}};
    my $transcript = $ff->create_prediction_transcript(\@exons, $self->query,
                                                       $self->analysis);
    $transcript->seqname($group);
    push(@transcripts, $transcript);
  }
  $self->output(\@transcripts);
}

