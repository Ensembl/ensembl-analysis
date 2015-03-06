package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB;

use strict;
use Carp;
use Bio::EnsEMBL::Hive::Utils ('stringify');
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Analysis::Tools::FeatureFactory;
use feature 'say';

use parent ('Bio::EnsEMBL::Hive::Process');

sub run{
  my ($self) = @_;
  foreach my $runnable(@{$self->runnable}){
    $runnable->run;
    $self->output($runnable->output);
  }
  return $self->param('output');
}

sub output{
  my ($self, $output) = @_;
  if(!$self->param('output')){
    $self->param('output') = [];
  }
  if($output){
    if(ref($output) ne 'ARRAY'){
      throw('Must pass RunnableDB:output an array ref not a '.$output);
    }
    push(@{$self->param('output')}, @$output);
  }
  return $self->param('output');
}

sub runnable {
  my ($self, $runnable) = @_;
  if(!$self->param('runnable')){
    $self->param('runnable',[]);
  }

  if($runnable){
    unless($runnable->isa('Bio::EnsEMBL::Analysis::Runnable')) {
      throw("Must pass RunnableDB:runnable a Bio::EnsEMBL::Analysis::Runnable not a ".$runnable);
    }
    push(@{$self->param('runnable')}, $runnable);
  }
  return $self->param('runnable');
}

#sub db {
#  my $self = shift;
#  my $db = shift;
#  if($db) {
#    unless($db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
#      throw("Must pass RunnableDB:db a Bio::EnsEMBL::DBSQL::DBAdaptor not a ".$db);
#    }
#    $self->param('db',$db);
#  }
#  return $self->param('db');
#}

sub query {
  my $self = shift;
  my $slice = shift;
  if($slice) {
    unless($slice->isa('Bio::EnsEMBL::Slice')) {
      throw("Must pass RunnableDB:query a Bio::EnsEMBL::Slice not a ".$slice);
    }
    $self->param('slice',$slice);
  }
  return $self->param('slice');
}

sub analysis {
  my $self = shift;
  my $analysis = shift;
  if($analysis){
    unless($analysis->isa('Bio::EnsEMBL::Analysis')) {
      throw("Must pass RunnableDB:analysis a Bio::EnsEMBL::Analysis not a ".$analysis);
    }
    $self->param('analysis',$analysis);
  }
  return $self->param('analysis');
}


sub input_id {
  my $self = shift;
  my $value = shift;

  # Note this sub is special. It overrides Hive::Process::input_id and parses the hive input_id into
  # a normal genebuild style input_id. The issue here is that for the moment I'm going to make this
  # a getter and not a setter as the two functions are somewhat different in this context. It would
  # be wrong to have the get function look for and parse the hive input id but then allow the set
  # function to set a new input id. Also a point to note is that overriding input_id in Process
  # should be fine as it is not the input_id call that the hive itself uses
  my $input_id_string = $self->Bio::EnsEMBL::Hive::Process::input_id;
  unless($input_id_string =~ /.+\=\>.+\"(.+)\"/) {
    throw("Could not parse the value from the input id. Input id string:\n".$input_id_string);
  }

  $input_id_string = $1;
  return($input_id_string);
}

sub hrdb_con {
  my ($self,$dba,$dba_con_name) = @_;
  if($dba) {
    if($dba_con_name){
      $self->param('_'.$dba_con_name,$dba);
    } else {
      $self->param('_hrdbadaptor',$dba);
    }
  }

  if($dba_con_name) {
    return $self->param('_'.$dba_con_name);
  } else {
    return $self->param('_hrdbadaptor');
  }

}

sub hrdb_get_dba {
  my ($self,$connection_info) = @_;
  my $dba;

  if(ref($connection_info)=~ m/HASH/) {
    eval {
      $dba = new Bio::EnsEMBL::DBSQL::DBAdaptor(%$connection_info);
    };

    if(@!) {
      throw("Error while setting up database connection:\n".@!);
    }
  }

  $dba->dbc->disconnect_when_inactive(1);
  return $dba;
}

sub feature_factory{
  my ($self, $feature_factory) = @_;
  if($feature_factory) {
    $self->param('feature_factory',$feature_factory);
  }
  if(!$self->param('feature_factory')) {
    $self->param('feature_factory',Bio::EnsEMBL::Analysis::Tools::FeatureFactory->new());
  }
  return $self->param('feature_factory');
}


sub fetch_sequence{
  my ($self, $name, $db, $repeat_masking, $soft_masking) = @_;
  if(!$db){
    $db = $self->db;
  }
  if(!$name){
    $name = $self->parse_hive_input_id;
  }
  my $sa = $db->get_SliceAdaptor;
  my $slice = $sa->fetch_by_name($name);
  $repeat_masking = [] unless($repeat_masking);
  if(!$slice){
    throw("Failed to fetch slice ".$name);
  }
  if(@$repeat_masking){
    my $sequence = $slice->get_repeatmasked_seq($repeat_masking, $soft_masking);
    $slice = $sequence
  }
  return $slice;
}

sub parameters_hash{
  my ($self, $string) = @_;

  if(!$string){
    $string = $self->analysis->parameters;
  }
  my %parameters_hash;

  if ($string) {
    if($string =~  /,/ || $string =~ /=>/){
      my @pairs = split (/,/, $string);
      foreach my $pair(@pairs){
        my ($key, $value) = split (/=>/, $pair);
        if ($key && ($value || $value eq '0')) {
          $key   =~ s/^\s+//g;
          $key   =~ s/\s+$//g;
          $value =~ s/^\s+//g;
          $value =~ s/\s+$//g;
          $parameters_hash{$key} = $value;
        } else {
          $parameters_hash{$key} = 1;
        }
      }
    }else{
      $parameters_hash{'-options'} = $string;
    }
  }
  return \%parameters_hash;
}

sub require_module{
  my ($self, $module) = @_;
  my $class;
  ($class = $module) =~ s/::/\//g;
  eval{
    require "$class.pm";
  };
  throw("Couldn't require ".$class." Blast:require_module $@") if($@);
  return $module;
}

sub ignore_config_file {
  my $self = shift;
  my $value = shift;
  if($value) {
    $self->param('ignore_config',$value) = shift if(@_);
  }
  return $self->param('ignore_config');
}

sub no_config_exception{
  my $self = shift;
  my $value = shift;
  if($value) {
    $self->param('no_config_exception',$value);
  }
  return $self->param('no_config_exception');
}

sub input_is_void{
  my $self = shift;
  my $value = shift;
  if($value) {
    $self->param('input_is_void',$value);
  }
  return $self->param('input_is_void');
}


sub failing_job_status{
  my $self = shift;
  my $value = shift;
  if($value) {
    $self->param('failing_status',$value);
  }
  return $self->param('failing_status');
}

1;
