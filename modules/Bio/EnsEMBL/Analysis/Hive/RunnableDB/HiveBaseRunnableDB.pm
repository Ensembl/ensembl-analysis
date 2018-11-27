#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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
use warnings;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Analysis::Tools::FeatureFactory;
use Bio::EnsEMBL::Hive::Utils ('destringify');
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(create_file_name);
use feature 'say';

use parent ('Bio::EnsEMBL::Hive::Process');


=head2 param_defaults

 Description: It allows the definition of default parameters for all inherting module.
              These are the default values:
               _input_id_name => 'iid',
               disconnect_jobs => 0,
               _branch_to_flow_to => 2,
 Returntype : Hashref, containing all default parameters
 Exceptions : None

=cut

sub param_defaults {
    my $self = shift;

    return {
        %{$self->SUPER::param_defaults},
        _input_id_name => 'iid',
        disconnect_jobs => 0,
        _branch_to_flow_to => 2,
        _branch_to_flow_to_on_fail => -3,
        _output => [],
    }
}


=head2 run

 Arg [1]    : None
 Description: For each runnable, execute the run method of the runnable and store the output
              in output. If the parameter 'disconnect_jobs' is set to 1, it will disconnect
              the worker from the database. Only use this is your job is longer than 15-30 minutes
 Returntype : Arrayref of object to be stored
 Exceptions : None

=cut

sub run {
  my ($self) = @_;
  $self->dbc->disconnect_if_idle() if ($self->param('disconnect_jobs'));
  foreach my $runnable(@{$self->runnable}){
    $runnable->run;
    if ($self->can('filter_results')) {
      $self->output($self->filter_results($runnable->output));
    }
    else {
      $self->output($runnable->output);
    }
  }
  return $self->output;
}


=head2 output

 Arg [1]    : Arrayref of objects
 Description: Getter/setter for the data you want to use in the write_output method. You pass in
              an arrayref of data which will be added to the data already present
 Returntype : Arrayref of objects
 Exceptions : Throws if the parameter is not an arrayref

=cut

sub output {
  my ($self, $output) = @_;

  unless($self->param_is_defined('_output')) {
    $self->param('_output',[]);
  }

  if($output){
    if(ref($output) ne 'ARRAY'){
      $self->throw('Must pass RunnableDB:output an array ref not a '.$output);
    }
    push(@{$self->param('_output')}, @$output);
  }
  return $self->param('_output');
}


=head2 write_output

 Arg [1]    : None
 Description: Write the objects stored in output in the output database. It fetches the adaptor
              using the get_adaptor method
              It uses a Bio::EnsEMBL::Analysis::Tools::FeatureFactory to validate the objects
 Returntype : Integer 1
 Exceptions : Throws when storing a feature fails

=cut

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
                   $adaptor->dbc->dbname."': ".$@);
    }
  }

  return 1;
} ## end sub write_output


=head2 get_adaptor

 Arg [1]    : None
 Description: Getter for the adaptor to use for storing objects. It has to be implemented in the inheriting
              class
 Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
 Exceptions : Throws if the method has not been implemented in the inheriting class

=cut

sub get_adaptor {
    my ($self) = @_;

    $self->throw('You should implement the method to return a Bio::EnsEMBL::DBSQL::DBAdaptor for your output database');
}


=head2 runnable

 Arg [1]    : (optional) Bio::EnsEMBL::Analysis::Runnable
 Description: Getter/setter for your Runnable objects
 Returntype : Arrayref of Bio::EnsEMBL::Analysis::Runnable
 Exceptions : Throws if the object is not a Bio::EnsEMBL::Analysis::Runnable

=cut

sub runnable {
  my ($self, $runnable) = @_;
  if(!$self->param_is_defined('runnable')){
    $self->param('runnable',[]);
  }

  if($runnable){
    unless($runnable->isa('Bio::EnsEMBL::Analysis::Runnable')) {
      $self->throw("Must pass RunnableDB:runnable a Bio::EnsEMBL::Analysis::Runnable not a ".$runnable);
    }
    push(@{$self->param('runnable')}, $runnable);
  }
  return $self->param('runnable');
}


=head2 query

 Arg [1]    : (optional) Bio::EnsEMBL::Slice
 Description: Getter/setter for the slice you will work on
 Returntype : Bio::EnsEMBL::Slice
 Exceptions : Throws if the object is not a Bio::EnsEMBL::Slice

=cut

sub query {
  my $self = shift;
  my $slice = shift;
  if($slice) {
    unless($slice->isa('Bio::EnsEMBL::Slice')) {
      $self->throw("Must pass RunnableDB:query a Bio::EnsEMBL::Slice not a ".$slice);
    }
    $self->param('slice',$slice);
  }
  return $self->param('slice');
}


=head2 analysis

 Arg [1]    : (optional) Bio::EnsEMBL::Analysis
 Description: Getter/setter for the analysis you will work on
 Returntype : Bio::EnsEMBL::Analysis
 Exceptions : Throws if the object is not a Bio::EnsEMBL::Analysis

=cut

sub analysis {
  my $self = shift;
  my $analysis = shift;
  if($analysis){
    unless($analysis->isa('Bio::EnsEMBL::Analysis')) {
      $self->throw("Must pass RunnableDB:analysis a Bio::EnsEMBL::Analysis not a ".$analysis);
    }
    $self->param('analysis',$analysis);
  }
  return $self->param('analysis');
}


=head2 input_id

 Example    : my $input_id = $self->input_id;
 Description: It returns the input_id for this analysis. By default it will look for the parameter 'iid'
              You can change the name of the parameter by using a method param_defaults and specifying a
              value for _input_id_name:
              sub param_defaults {
                  my $super_param_defaults = $self->param_defaults;
                  $super_param_defaults->{_input_id_name} = 'filename';
                  return $super_param_defaults;
              }

 Returntype : String
 Exceptions : Throws if the input id is not set

=cut

sub input_id {
  my $self = shift;

  $self->throw("Could not fetch your input id ".$self->param('_input_id_name')) unless ($self->param_is_defined($self->param('_input_id_name')));
  return $self->param($self->param('_input_id_name'));
}


sub input_id_array {
  my ($self) = @_;

  if (ref($self->input_id) eq 'ARRAY') {
    return $self->input_id;
  }
  else {
    return [$self->input_id];
  }
}

=head2 hrdb_set_con

 Arg [1]    : Bio::EnsEMBL::DBSQL::DBAdaptor
 Arg [2]    : (optional) String
 Example    : $self->hrdb_set_con($self->hrdb_get_dba('input_db'), 'input_db');
 Description: Setter for a DBAdaptor you will use in run or write_output or any other method
 Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
 Exceptions : Throws if Arg[1] is not a Bio::EnsEMBL::DBSQL::DBAdaptor object

=cut

sub hrdb_set_con {
  my ($self,$dba,$dba_con_name) = @_;

  unless($dba->isa('Bio::EnsEMBL::DBSQL::DBAdaptor')) {
    $self->throw("Expected a DBAdaptor object as input. If you want to retrieve a DBAdaptor then ".
          "use the getter sub instead (hrdb_get_con)");
  }

  if($dba_con_name){
      $self->param('_'.$dba_con_name,$dba);
  } else {
      $self->param('_hrdbadaptor',$dba);
  }

}


=head2 hrdb_get_con

 Arg [1]    : (optional) String
 Example    : $self->hrdb_get_con('input_db');
 Description: Getter for a DBAdaptor you will use in run or write_output or any other method
 Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
 Exceptions : None

=cut

sub hrdb_get_con {
  my ($self,$dba_con_name) = @_;

  if($dba_con_name) {
    return $self->param('_'.$dba_con_name);
  } else {
    return $self->param('_hrdbadaptor');
  }
}


=head2 feature_factory

 Arg [1]    : (optional) Bio::EnsEMBL::Analysis::Tools::FeatureFactory
 Description: Getter/setter
 Returntype : Bio::EnsEMBL::Analysis::Tools::FeatureFactory
 Exceptions : None

=cut

sub feature_factory {
  my ($self, $feature_factory) = @_;
  if($feature_factory) {
    $self->param('feature_factory',$feature_factory);
  }
  if(!$self->param_is_defined('feature_factory')) {
    $self->param('feature_factory',Bio::EnsEMBL::Analysis::Tools::FeatureFactory->new());
  }
  return $self->param('feature_factory');
}


=head2 fetch_sequence

 Arg [1]    : (optional) String, Ensembl formatted sequence name, default is using the $self->input_id
 Arg [2]    : (optional) Bio::EnsEMBL::DBSQL::DBAdaptor if not present use Arg[5] to get the object
 Arg [3]    : (optional) Arrayref of String, the logic_names of repeat masking analysis
 Arg [4]    : (optional) Boolean 1 if you want to softmask (lowercase), 0 for hardmasking (Ns)
 Arg [5]    : (optional) String, database name which can be fetch with hrdb_get_con
 Description: Fetch the sequence specified by the input_id or using Ensembl formatted name
 Returntype : Bio::EnsEMBL::Slice
 Exceptions : Throws if it fails to fetch the slice

=cut

sub fetch_sequence {
  my ($self, $name, $dbcon, $repeat_masking, $soft_masking, $dbname) = @_;
  if(!$dbcon){
    $dbcon = $self->hrdb_get_con($dbname);
  }
  if(!$name){
    $name = $self->input_id;
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


=head2 parameters_hash

 Arg [1]    : (optional) String, default uses values from parameters stored in the analysis object
 Description: Break down the string in key values pairs using ',' and '=>' as separator and store
              them in a hash. If there is no separators, the key is '-options' and the value is the
              string
 Returntype : Hashref of String
 Exceptions : None

=cut

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


=head2 require_module

 Arg [1]    : String name of a module
 Example    : $self->require_module('Bio::EnsEMBL::Analysis::Tools::Utilities');
 Description: Load dynamically a module
 Returntype : String, name of the module
 Exceptions : Throws if it couldn't load the module

=cut

sub require_module {
  my ($self, $module) = @_;
  my $class;
  ($class = $module) =~ s/::/\//g;
  eval{
    require "$class.pm";
  };
  $self->throw("Couldn't require ".$class.' '.ref($self)."::require_module $@") if($@);
  return $module;
}


=head2 ignore_config_file

 Deprecated. This is from the old pipeline and probably useless now
 Arg [1]    : (optional) String
 Description: Getter/setter
 Returntype : String
 Exceptions : None

=cut

sub ignore_config_file {
  my $self = shift;

  $self->warning('DEPRECATED: ignore_config_file was used in the old pipeline, probably useless now, called from '.__PACKAGE__);
  my $value = shift;
  if($value) {
    $self->param('ignore_config',$value) = shift if(@_);
  }
  return $self->param('ignore_config');
}


=head2 no_config_exception

 Deprecated. This is from the old pipeline and probably useless now
 Arg [1]    : (optional) String
 Description: Getter/setter
 Returntype : String
 Exceptions : None

=cut

sub no_config_exception {
  my $self = shift;
  $self->warning('DEPRECATED: no_config_exception was used in the old pipeline, probably useless now, called from '.__PACKAGE__);
  my $value = shift;
  if($value) {
    $self->param('no_config_exception',$value);
  }
  return $self->param('no_config_exception');
}


=head2 input_is_void

 Deprecated. This is from the old pipeline and probably useless now
 Arg [1]    : (optional) String
 Description: Getter/setter
 Returntype : String
 Exceptions : None

=cut

sub input_is_void {
  my $self = shift;
  $self->warning('DEPRECATED: input_is_void was used in the old pipeline, probably useless now, called from '.__PACKAGE__);
  my $value = shift;
  if($value) {
    $self->param('input_is_void',$value);
  }
  return $self->param('input_is_void');
}


=head2 failing_job_status

 Deprecated. This is from the old pipeline and probably useless now
 Arg [1]    : (optional) String
 Description: Getter/setter
 Returntype : String
 Exceptions : None

=cut

sub failing_job_status {
  my $self = shift;
  $self->warning('DEPRECATED: failing_job_status was used in the old pipeline, probably useless now, called from '.__PACKAGE__);
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

    my $logic_name;
    if ($self->param_is_defined('logic_name')) {
        $logic_name = $self->param('logic_name');
    }
    else {
      $logic_name = $self->input_job->analysis->logic_name;
    }
    if ($self->param_is_defined('analysis_params')) {
      while( my ($key, $value) = each %{$self->param('analysis_params')}) {
        $extra_params->{$key} = $value;
      }
    }
    my $analysis = Bio::EnsEMBL::Analysis->new(-logic_name => $logic_name, %$extra_params);
    if ($add_module) {
      my $module;
      if ($self->param_is_defined('module')) {
        $module = $self->param('module');
      }
      else {
        $self->param('module', $self->input_job->analysis->module);
      }
      $analysis->module($self->param('module'));
    }
    $analysis->program_file($self->param('program_file')) if ($self->param_is_defined('program_file'));

    return $self->analysis($analysis);
}


=head2 hrdb_get_dba

 Arg [1]    : String $name, name of a database as it stored in parameters
 Arg [2]    : Bio::EnsEMBL::DBSQL::DBAdaptor object, the database will have the dna (optional)
 Arg [3]    : String $alternative_class (optional), Allowed class are Variation, Compara, Funcgen
 Example    : $self->hrdb_get_dba($self->param('target_db'));
 Description: It's a wrapper for hrdb_get_dba from Bio::EnsEMBL::Analysis::Tools::Utilities
 Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor
 Exceptions : Throws if it cannot connect to the database.
              Throws if $connection_info is not a hashref
              Throws if $dna_db is not a Bio::EnsEMBL::DBSQL::DBAdaptor object

=cut

sub  hrdb_get_dba {
    my ($self, $connection_info, $dna_db, $alternative_class) = @_;

    $self->throw($connection_info.' is not a HASHREF') unless (ref($connection_info) eq 'HASH');
    my $uniq_id = join(':', $connection_info->{-host},
                            $connection_info->{-dbname},
                            $connection_info->{-port},
                            $connection_info->{-user});
    if (exists $self->{_gb_cache}->{'_cache_lastlogicname'}
        and $self->{_gb_cache}->{'_cache_lastlogicname'} ne $self->input_job->analysis->logic_name) {
      delete $self->{_gb_cache};
    }
    if (!exists $self->{_gb_cache}->{'_cache_dba_'.$uniq_id}) {
      $self->{_gb_cache}->{'_cache_dba_'.$uniq_id} = Bio::EnsEMBL::Analysis::Tools::Utilities::hrdb_get_dba($connection_info, $dna_db, $alternative_class);
      $self->{_gb_cache}->{'_cache_lastlogicname'} = $self->input_job->analysis->logic_name;
    }
    return $self->{_gb_cache}->{'_cache_dba_'.$uniq_id};
}


=head2 get_database_by_name

 Arg [1]    : String $name, name of a database as it stored in parameters
 Arg [2]    : Bio::EnsEMBL::DBSQL::DBAdaptor object (optional), the database will have the dna
 Arg [3]    : String $alternative_class (optional), Allowed class are Variation, Compara, Funcgen
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
    my ($self, $name, $dna_db, $alternative_class) = @_;

    return $self->hrdb_get_dba(destringify($self->param($name)), $dna_db, $alternative_class);
}


=head2 is_slice_name

 Arg [1]    : String, string to check
 Example    : $self->is_slice_name($input_id);
 Description: Return 1 if the string given is an Ensembl slice name
 Returntype : Boolean
 Exceptions : None

=cut

sub is_slice_name {
    my ($self, $string) = @_;

    return Bio::EnsEMBL::Analysis::Tools::Utilities::is_slice_name($string);
}


=head2 create_target_file

 Arg [1]    : (optional) String, prefix for the file
 Arg [2]    : (optional) String, extension of the file
 Arg [3]    : (optional) String, directory to write to
 Description: Create a temporary file for the target sequence.
              It should be called without argument but mostly
              without Arg[3] as you want to use /tmp
              The object can be called with 'target_file'
 Returntype : String filename
 Exceptions : None

=cut

sub create_target_file {
  my ($self, $stem, $suffix, $dir) = @_;

  $stem = 'target' unless ($stem);
  return $self->_create_temporary_file('target_file', $stem, $suffix, $dir);
}


=head2 create_query_file

 Arg [1]    : (optional) String, prefix for the file
 Arg [2]    : (optional) String, extension of the file
 Arg [3]    : (optional) String, directory to write to
 Description: Create a temporary file for the query sequence.
              It should be called without argument but mostly
              without Arg[3] as you want to use /tmp
              The object can be called with 'query_file'
 Returntype : String filename
 Exceptions : None

=cut

sub create_query_file {
  my ($self, $stem, $suffix, $dir) = @_;

  $stem = 'query' unless ($stem);
  return $self->_create_temporary_file('query_file', $stem, $suffix, $dir);
}


=head2 _create_temporary_file

 Arg [1]    : String, name of parameter to store
 Arg [2]    : (optional) String, prefix for the file
 Arg [3]    : (optional) String, extension of the file
 Arg [4]    : (optional) String, directory to write to
 Description: Create a temporary file.
 Returntype : String filename
 Exceptions : None

=cut

sub _create_temporary_file {
  my ($self, $name, $stem, $suffix, $dir) = @_;

  $self->param($name, create_file_name($stem, $suffix, $dir));
  return $self->param($name)->filename;
}


=head2 post_cleanup

 Arg [1]    : (optional) int, secs to sleep
 Description: If a module use post_cleanup it needs to call this first
              otherwise the DBAdaptor cache will not be reset which might
              cause problems.
              Generic post_cleanup implementation. Note that this is designed to get
              runnable dbs to sleep for a specified amount of time (default is 5 secs)
              before finishing. We needed this because we have found issues during the
              with downstream jobs starting before the write from the upstream job is
              actually fully written to disk. This can have unusual effects in terms of
              creating situations where input is missing and the job silently fails
 Returntype : None
 Exceptions : None

=cut

sub post_cleanup {
  my ($self,$sleep_time) = @_;

  foreach my $key (keys %{$self->{_gb_cache}}) {
    next if ($key eq '_cache_lastlogicname');
    $self->{_gb_cache}->{$key}->clear_caches;
  }

  if(defined($sleep_time)) {
    unless($sleep_time =~ /^[0-9]+$/) {
      $self->throw("Value passed in for sleep time was not an positive integer or zero. Value: ".$sleep_time);
    }
    sleep($sleep_time);
  } else {
    sleep(5);
  }

}

1;
