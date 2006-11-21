package Bio::EnsEMBL::Analysis::RunnableDB::Pmatch;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::Pmatch qw(PMATCH_BY_LOGIC);
use Bio::EnsEMBL::Analysis::Runnable::Pmatch;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::Databases;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw (rearrange);

@ISA = qw (
           Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild
           );


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->read_and_check_config($PMATCH_BY_LOGIC);

  return $self;
}



sub fetch_input{
  my ($self) = @_;
  my $slice = $self->fetch_sequence($self->input_id, $self->db, 
                                    $self->REPEAT_MASKING);
  $self->query($slice);
  my %parameters = %{$self->parameters_hash};
  my $program = $self->analysis->program_file;
  $program = $self->BINARY_LOCATION if(!$program);
  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Pmatch
    ->new(
          -query => $self->query,
          -program => $program,
          -analysis => $self->analysis,
          -protein_file => $self->PROTEIN_FILE,
          -max_intron_length => $self->MAX_INTRON_LENGTH,
          -min_coverage => $self->MIN_COVERAGE,
         );
  $self->runnable($runnable);
}


sub PROTEIN_FILE{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'PROTEIN_FILE'} = $arg;
  }
  return $self->{'PROTEIN_FILE'};
}


sub MIN_COVERAGE{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'MIN_COVERAGE'} = $arg;
  }
  return $self->{'MIN_COVERAGE'};
}

sub BINARY_LOCATION{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'BINARY_LOCATION'} = $arg;
  }
  return $self->{'BINARY_LOCATION'};
}


sub MAX_INTRON_LENGTH{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'MAX_INTRON_SIZE'} = $arg;
  }
  return $self->{'MAX_INTRON_SIZE'};
}


sub OUTPUT_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'OUTPUT_DB'} = $arg;
  }
  return $self->{'OUTPUT_DB'};
}

sub REPEAT_MASKING{
  my ($self, $arg) = @_;
  if($arg){
    throw("Runnable::Pmatch ".$arg." must be an array ref of logic names not .".$arg)
      unless(ref($arg) eq 'ARRAY');
    $self->{'REPEAT_MASKING'} = $arg;
  }
  return $self->{'REPEAT_MASKING'};
}



sub read_and_check_config {
  my $self = shift;

  $self->SUPER::read_and_check_config($PMATCH_BY_LOGIC);
  
  #######
  #CHECKS
  #######
  foreach my $config_var (qw(PROTEIN_FILE
                             OUTPUT_DB)){
    throw("You must define $config_var in config for logic '".
          $self->analysis->logic_name."'")
      if not defined $self->$config_var;
  }
};


sub get_adaptor{
  my ($self) = @_;
  my $output_db = $self->get_dbadaptor($self->OUTPUT_DB);
  return $output_db->get_ProteinAlignFeatureAdaptor;
}
