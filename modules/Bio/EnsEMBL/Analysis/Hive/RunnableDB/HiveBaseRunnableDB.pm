package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB;

use strict;
use Carp;
use Bio::EnsEMBL::Hive::Utils ('stringify');
use feature 'say';

use parent ('Bio::EnsEMBL::Hive::Process','Bio::EnsEMBL::Analysis::RunnableDB');

sub runnable {
  my ($self, $runnable) = @_;
  if(!$self->param('runnable')){
    $self->param('runnable',[]);
  }
  if($runnable){
    throw("Must pass RunnableDB:runnable a ".
          "Bio::EnsEMBL::Analysis::Runnable not a ".$runnable)
      unless($runnable->isa('Bio::EnsEMBL::Analysis::Runnable'));
    push(@{$self->param('runnable')}, $runnable);
  }
  return $self->param('runnable');
}

sub analysis {
  my $self = shift;
  my $analysis = shift;
  if($analysis){
    throw("Must pass RunnableDB:analysis a Bio::EnsEMBL::Analysis".
          "not a ".$analysis) unless($analysis->isa
                                     ('Bio::EnsEMBL::Analysis'));
    $self->param('analysis',$analysis);
  }
  return $self->param('analysis');
}

sub parse_hive_input_id {
  my $self = shift;
  my $input_id_string = $self->Bio::EnsEMBL::Hive::Process::input_id;
  unless($input_id_string =~ /.+\=\>.+\"(.+)\"/) {
    throw("Could not parse the value from the input id. Input id string:\n".$input_id_string);
  }

  $input_id_string = $1;
  return($input_id_string);

}

sub _slurp {
  my ($self, $file_name) = @_;
  my $slurped;
  {
    local $/ = undef;
    open(my $fh, '<', $file_name) or $self->throw("Couldnt open file [$file_name]");
    $slurped = <$fh>;
    close($fh);
  }
  return $slurped;
}


1;
