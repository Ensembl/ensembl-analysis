package Bio::EnsEMBL::Analysis::RunnableDB::BestPmatch;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild;
use Bio::EnsEMBL::Analysis::Config::GeneBuild::Pmatch qw(BESTPMATCH_BY_LOGIC);
use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw (rearrange);
use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils qw(empty_Object);
use Bio::EnsEMBL::Analysis::Runnable::BestPmatch;

@ISA = qw (
           Bio::EnsEMBL::Analysis::RunnableDB::BaseGeneBuild
           );


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->read_and_check_config($BESTPMATCH_BY_LOGIC);

  return $self;
}


sub fetch_input{
  my ($self) = @_;
  my $pafa = $self->get_dbadaptor($self->INPUT_DB)->get_ProteinAlignFeatureAdaptor;
  my $features = $pafa->fetch_all_by_logic_name($self->PMATCH_LOGIC_NAME);
  print "Have fetched ".@$features." with ".$self->PMATCH_LOGIC_NAME." from ".
    " ".$pafa->dbc->dbname."\n";
  $self->pmatch_features($features);
}

sub run{
  my ($self) = @_;
  my $pmf2 = Bio::EnsEMBL::Analysis::Runnable::BestPmatch->
    new( 
        '-protein_hits' => $self->pmatch_features,
        '-min_coverage' => $self->MIN_COVERAGE,
        '-analysis' => $self->analysis,
       );
  $pmf2->run;
  my @output = @{$pmf2->output};
  foreach my $output(@output){
    empty_Object($output);
  }
  #my @unique = $self->uniquify(@output);
  $self->output(\@output);
}


sub get_adaptor{
  my ($self) = @_;
  my $output_db = $self->get_dbadaptor($self->OUTPUT_DB);
  return $output_db->get_ProteinAlignFeatureAdaptor;
}

sub pmatch_features{
  my ($self, $features) = @_;
  if($features){
    if(ref($features) eq "ARRAY"){
      push(@{$self->{'pmatch_features'}}, @$features);
    }else{
      push(@{$self->{'pmatch_features'}}, $features);
    }
  }
  return $self->{'pmatch_features'};
}

sub PMATCH_LOGIC_NAME{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'PMATCH_LOGIC_NAME'} = $arg;
  }
  return $self->{'PMATCH_LOGIC_NAME'};
}

sub OUTPUT_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'OUTPUT_DB'} = $arg;
  }
  return $self->{'OUTPUT_DB'};
}



sub INPUT_DB{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'INPUT_DB'} = $arg;
  }
  return $self->{'INPUT_DB'};
}


sub MIN_COVERAGE{
  my ($self, $arg) = @_;
  if($arg){
    $self->{'MIN_COVERAGE'} = $arg;
  }
  return $self->{'MIN_COVERAGE'} = $arg;
}





sub read_and_check_config {
  my $self = shift;

  $self->SUPER::read_and_check_config($BESTPMATCH_BY_LOGIC);
  
  #######
  #CHECKS
  #######
  foreach my $config_var (qw(PMATCH_LOGIC_NAME
                             OUTPUT_DB
                             INPUT_DB)){
    throw("You must define $config_var in config for logic '".
          $self->analysis->logic_name."'")
      if not defined $self->$config_var;
  }
};


