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
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Analysis::Tools::FeatureFactory;
use Bio::EnsEMBL::Hive::Utils ('destringify');
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(hrdb_get_dba);
use feature 'say';

use parent ('Bio::EnsEMBL::Hive::Process');

sub run {
  my ($self) = @_;
  foreach my $runnable(@{$self->runnable}){
    $runnable->run;
    $self->output($runnable->output);
  }
  return $self->output;
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
  if(!$self->param_is_defined('runnable')){
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

=head2 create_analysis

 Arg [1]    : Boolean $add_module, if set to 1, it will add the module used by the analysis in the eHive pipeline
 Arg [2]    : Hashref $extra_params, contains extra parameters like -program_file to use for populating the Bio::EnsEMBL::Analysis object
 Example    : $self->create_analysis;
 Description: Create an Bio::EnsEMBL::Analysis object. If your analysis has a logic_name parameter, it will be used. Otherwise
              it will use the logic_name from Hive. It wil also store the analysis created in $self->analysis
 Returntype : Bio::EnsEMBL::Analysis
 Exceptions : None


=cut

sub create_analysis {
    my ($self, $add_module, $extra_params) = @_;

    my $logic_name = $self->param_is_defined('logic_name') ? $self->param('logic_name') : $self->input_job->analysis->logic_name;
    my $analysis = Bio::EnsEMBL::Analysis->new(-logic_name => $logic_name, %$extra_params);
    $analysis->module($self->input_job->analysis->module) if (defined $add_module);

    return $self->analysis($analysis);
}

=head2 hrdb_get_dba

 Arg [1]    : String $name, name of a database as it stored in parameters
 Arg [2]    : Bio::EnsEMBL::DBSQL::DBAdaptor object, the database will have the dna (optional)
 Example    : $self->hrdb_get_dba($self->param('target_db'));
 Description: It's a wrapper for hrdb_get_dba from Bio::EnsEMBL::Analysis::Tools::Utilities
 Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
 Exceptions : Throws if it cannot connect to the database.
              Throws if $connection_info is not a hashref
              Throws if $dna_db is not a Bio::EnsEMBL::DBSQL::DBAdaptor object


=cut

sub  hrdb_get_dba {
    my ($self, $connection_info, $dna_db) = @_;

    return Bio::EnsEMBL::Analysis::Tools::Utilities::hrdb_get_dba($connection_info, $dna_db);
}
=head2 get_database_by_name

 Arg [1]    : String $name, name of a database as it stored in parameters
 Arg [2]    : Bio::EnsEMBL::DBSQL::DBAdaptor object, the database will have the dna (optional)
 Example    : $self->get_database_by_name('target_db');
 Description: It creates a object based on the information contained in $connection_info.
              If the hasref contains -dna_db or if the second argument is populated, it will
              try to attach the DNA database
 Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
 Exceptions : Throws if it cannot connect to the database.
              Throws if $connection_info is not a hashref
              Throws if $dna_db is not a Bio::EnsEMBL::DBSQL::DBAdaptor object

=cut

sub get_database_by_name {
    my ($self, $name, $dna_db) = @_;

    return $self->hrdb_get_dba(destringify($self->param($name)), $dna_db);
}

=head2 is_slice_name

 Arg [1]    : String, string to check
 Example    : $self->is_slice_name($input_id);;
 Description: Return 1 if the string given is an Ensembl slice name
 Returntype : Boolean
 Exceptions : None


=cut

sub is_slice_name {
    my ($self, $string) = @_;

    return $string =~ /^[^:]+:[^:]+:[^:]+:\d+:\d+:(1|-1)$/;
}

1;
