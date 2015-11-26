#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB;

use strict;
use Carp;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Analysis::Tools::FeatureFactory;
use Bio::EnsEMBL::Hive::Utils ('destringify');
use feature 'say';

use parent ('Bio::EnsEMBL::Hive::Process');

sub run {
  my ($self) = @_;
  foreach my $runnable(@{$self->runnable}){
    $runnable->run;
    $self->output($runnable->output);
  }
  return $self->param('output');
}

sub output {
  my ($self, $output) = @_;
  unless($self->param_is_defined('_output')){
    $self->param('_output',[]);
  }
  if($output){
    if(ref($output) ne 'ARRAY'){
      throw('Must pass RunnableDB:output an array ref not a '.$output);
    }
    push(@{$self->param('_output')}, @$output);
  }
  return $self->param('_output');
}


sub write_output {
  my ($self) = @_;

  my $adaptor  = $self->get_adaptor();
  my $analysis = $self->analysis();

  foreach my $feature ( @{ $self->output() } ) {
    $feature->analysis($analysis);

    if ( !defined( $feature->slice() ) ) {
      $feature->slice( $self->query() );
    }

    $self->feature_factory->validate($feature);

    eval { $adaptor->store($feature); };
    if ($@) {
      $self->throw("RunnableDB::write_output() failed: failed to store '".$feature."' into database '".
                   $self->hrdb_get_con('target_db')->dbname."': ".$@);
    }
  }

  return 1;
} ## end sub write_output

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


sub hrdb_set_con {
  my ($self,$dba,$dba_con_name) = @_;

  unless($dba->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
    throw("Expected a DBAdaptor object as input. If you want to retrieve a DBAdaptor then ".
          "use the getter sub instead (hrdb_get_con)");
  }

  if($dba_con_name){
      $self->param('_'.$dba_con_name,$dba);
  } else {
      $self->param('_hrdbadaptor',$dba);
  }

}


sub hrdb_get_con {
  my ($self,$dba_con_name) = @_;

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

    if($@) {
      throw("Error while setting up database connection:\n".$@);
    }
  } else {
    throw("DB connection info passed in was not a hash:\n".$connection_info);
  }

  $dba->dbc->disconnect_when_inactive(1);
  return $dba;
}


sub feature_factory {
  my ($self, $feature_factory) = @_;
  if($feature_factory) {
    $self->param('feature_factory',$feature_factory);
  }
  if(!$self->param('feature_factory')) {
    $self->param('feature_factory',Bio::EnsEMBL::Analysis::Tools::FeatureFactory->new());
  }
  return $self->param('feature_factory');
}


sub fetch_sequence {
  my ($self, $name, $dbcon, $repeat_masking, $soft_masking, $dbname) = @_;
  if(!$dbcon){
    $dbcon = $self->hrdb_get_con($dbname);
  }
  if(!$name){
    $name = $self->parse_hive_input_id;
  }
  my $sa = $dbcon->get_SliceAdaptor;
  my $slice = $sa->fetch_by_name($name);
  $repeat_masking = [] unless($repeat_masking);
  if(!$slice){
    $self->throw("Failed to fetch slice ".$name);
  }
  if(@$repeat_masking){
    my $sequence = $slice->get_repeatmasked_seq($repeat_masking, $soft_masking);
    $slice = $sequence;
  }
  return $slice;
}

sub parameters_hash {
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

sub require_module {
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

sub no_config_exception {
  my $self = shift;
  my $value = shift;
  if($value) {
    $self->param('no_config_exception',$value);
  }
  return $self->param('no_config_exception');
}

sub input_is_void {
  my $self = shift;
  my $value = shift;
  if($value) {
    $self->param('input_is_void',$value);
  }
  return $self->param('input_is_void');
}


sub failing_job_status {
  my $self = shift;
  my $value = shift;
  if($value) {
    $self->param('failing_status',$value);
  }
  return $self->param('failing_status');
}

sub create_analysis {
    my ($self, $add_module, $extra_params) = @_;

    my $logic_name = $self->param_is_defined('logic_name') ? $self->param('logic_name') : $self->input_job->analysis->logic_name;
    my $analysis = Bio::EnsEMBL::Analysis->new(-logic_name => $logic_name, %$extra_params);
    $analysis->module($self->input_job->analysis->module) if (defined $add_module);

    return $self->analysis($analysis);
}

sub get_database_by_name {
    my ($self, $name) = @_;

    return $self->hrdb_get_dba(destringify($self->param($name)));
}

sub is_slice_name {
    my ($self, $string) = @_;

    return $string =~ /^[^:]+:[^:]+:[^:]+:\d+:\d+:(1|-1)$/;
}

1;
