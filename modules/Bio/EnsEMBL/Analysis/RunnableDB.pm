
package Bio::EnsEMBL::Analysis::RunnableDB;

# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB

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

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB

=head1 SYNOPSIS

  my $repeat_masker = Bio::EnsEMBL::Analysis::RunnableDB::RepeatMasker->
  new(
      -input_id => 'contig::AL805961.22.1.166258:1:166258:1',
      -db => $db,
      -analysis => $analysis,
     );
  $repeat_masker->fetch_input;
  $repeat_masker->run;
  $repeat_masker->write_output;

=head1 DESCRIPTION

This module acts as a base class for our RunnableDBs who act as an 
interface between the core database and our Runnables both fetching
input data and writing data back to the databases

The module provides some base functionality some of which must be
used

The constructor fits the model the pipeline which runs most of
our analyses expects. If a child runnabledb expects more arguments
to the constructor than this one it wont be directly runnable by the
pipeline

Most of the other methods provided are containers of some description
but there are some methods with specific functionality

parameters_hash is there to parse a string from the parameters varible
in the analysis object. This string should have the format
key => value, key => value where the key would be the Runnables
constructor argument and the value the variable. This is to allow
some flexibility in the arguments expected by and the way we run
Runnables.

fetch_sequence fetched a sequence using the fetch_by_name method of 
the slice adaptor from the given database. The name, database and an
array of logic_names to determine masking can be given. If no name
or database is provided the method defaults to input_id and db

validate, this is a method which does some basic validation of the
feature before storage. This checks if slice and analysis object are
defined and if start, end and strand are defined then that the start
is smaller than the end and both the start and end are > 0

All runnableDBs need to implement 3 methods to run within the pipeline

fetch_input, this always must be implemented by individual child
RunnableDBs as the input required for different analyses can vary so
widely that it is impossible to write a generic method

run, there is a run method implemented. To use this child runnabledbs
need to have added the runnables they want run to the runnable method
which holds an array of runnables which are each called and the output 
stored in this method

write_output, there is also a generic implementation of this. To use this
method the child runnabledb must implement a get_adaptor method which
returns the appropriate adaptor to be used in storage. 

=head1 CONTACT

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

=cut

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Analysis::Tools::FeatureFactory;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(parse_config parse_config_mini parse_config_value);
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_info logger_verbosity);
use Bio::EnsEMBL::Analysis::Config::General qw(CORE_VERBOSITY LOGGER_VERBOSITY);

use vars qw (@ISA);

@ISA = qw();


=head2 new

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
  Arg [3]   : Bio::EnsEMBL::Analysis
  Function  : create a Bio::EnsEMBL::Analysis::RunnableDB object
  Returntype: Bio::EnsEMBL::Analysis::RunnableDB
  Exceptions: throws if not passed either a dbadaptor, input id or
  an analysis object
  Example   : $rdb = $perl_path->new( -analysis => $self->analysis,
                                      -input_id => $self->input_id,
                                      -db => $self->adaptor->db );

=cut



sub new{
  my ($class,@args) = @_;
  my $self = bless {},$class;  

  my ($db, $input_id, $analysis,$ignore_config_file,$no_config_exception, $verbosity) = rearrange (['DB', 'INPUT_ID', 'ANALYSIS','IGNORE_CONFIG_FILE','NO_CONFIG_EXCEPTION', 'VERBOSITY'], @args);

  if(!$db || !$analysis || !$input_id){
    throw("Can't create a RunnableDB without a dbadaptor ".
          $db." an analysis object ".$analysis.
          " or an input_id ".$input_id);
  }
 

  #Clone analysis to prevent analysis reference problem when
  #using separate pipeline and output DBs
  #Do not use adaptor here as caching returns same reference
  my $cloned_analysis;
  %{$cloned_analysis} = %{$analysis};
  $analysis =  bless $cloned_analysis, ref ($analysis);


  $self->db($db);
  $self->analysis($analysis);
  $self->input_id($input_id); 
  $self->ignore_config_file($ignore_config_file) ;
  $self->no_config_exception($no_config_exception) ;

  verbose($CORE_VERBOSITY);
  logger_verbosity($LOGGER_VERBOSITY) unless ($verbosity);
  return $self;
}



=head2 db

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : Bio::EnsEMBL::Pipeline::DBAdaptor
  Function  : container for dbadaptor
  Returntype: Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor
  Exceptions: throws if not passed a Bio::EnsEMBL::DBSQL::DBConnection
  object
  Example   : 

=cut


sub db{
  my $self = shift;
  my $db = shift;
  if($db){
    throw("Must pass RunnableDB:db a Bio::EnsEMBL::DBSQL::DBAdaptor ".
          "not a ".$db) 
      unless($db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'));
    $self->{'db'} = $db;
  }
  return $self->{'db'};
}



=head2 analysis

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : Bio::EnsEMBL::Analysis
  Function  : container for analysis object
  Returntype: Bio::EnsEMBL::Analysis
  Exceptions: throws passed incorrect object type
  Example   : 

=cut



sub analysis{
  my $self = shift;
  my $analysis = shift;
  if($analysis){
    throw("Must pass RunnableDB:analysis a Bio::EnsEMBL::Analysis".
          "not a ".$analysis) unless($analysis->isa
                                     ('Bio::EnsEMBL::Analysis'));
    $self->{'analysis'} = $analysis;
  }
  return $self->{'analysis'};
}



=head2 query

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : Bio::EnsEMBL::Slice
  Function  : container for slice object
  Returntype: Bio::EnsEMBL::Slice
  Exceptions: throws if passed the incorrect object type
  Example   : 

=cut


sub query{
  my $self = shift;
  my $slice = shift;
  if($slice){
    throw("Must pass RunnableDB:query a Bio::EnsEMBL::Slice".
          "not a ".$slice) unless($slice->isa
                                     ('Bio::EnsEMBL::Slice'));
    $self->{'slice'} = $slice;
  }
  return $self->{'slice'};
}


=head2 runnable

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : Bio::EnsEMBL::Analysis::Runnable
  Function  : container for an array of runnables
  Returntype: arrayref
  Exceptions: throws if passed the wrong object type
  Example   : 

=cut



sub runnable{
  my ($self, $runnable) = @_;
  if(!$self->{'runnable'}){
    $self->{'runnable'} = [];
  }
  if($runnable){
    throw("Must pass RunnableDB:runnable a ".
          "Bio::EnsEMBL::Analysis::Runnable not a ".$runnable) 
      unless($runnable->isa('Bio::EnsEMBL::Analysis::Runnable'));
    push(@{$self->{'runnable'}}, $runnable);
  }
  return $self->{'runnable'};
}


=head2 input_id 

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : string/int
  Function  : container for specified variable. This pod refers to the
  three methods below, input_id, input_is_void, failing_job_status. This
  are simple containers which dont do more than hold and return an given
  value
  Returntype: string/int
  Exceptions: none
  Example   : 

=cut


sub input_id{
  my $self = shift;
  $self->{'input_id'} = shift if(@_);
  return $self->{'input_id'};
}


sub input_is_void{
  my $self = shift;
  $self->{'input_is_void'} = shift if(@_);
  return $self->{'input_is_void'};
}


sub failing_job_status{
  my $self = shift;
  $self->{'failing_status'} = shift if(@_);
  return $self->{'failing_status'};
}

=head2 output

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : arrayref of output
  Function  : push array passed into onto output array
  Returntype: arrayref
  Exceptions: throws if not passed correct type
  Example   : $self->output(\@output);

=cut


sub output{
  my ($self, $output) = @_;
  if(!$self->{'output'}){
    $self->{'output'} = [];
  }
  if($output){
    if(ref($output) ne 'ARRAY'){
      throw('Must pass RunnableDB:output an array ref not a '.$output);
    }
    push(@{$self->{'output'}}, @$output);
  }
  return $self->{'output'};
}


=head2 feature_factory

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : Bio::EnsEMBL::Analysis::Tools::FeatureFactory
  Function  : container for a feature factory object. If none is defined
  when one is requested a new one is created. 
  Returntype: Bio::EnsEMBL::Analysis::Tools::FeatureFactory
  Exceptions: none
  Example   : 

=cut



sub feature_factory{
  my ($self, $feature_factory) = @_;
  if($feature_factory){
    $self->{'feature_factory'} = $feature_factory;
  }
  if(!$self->{'feature_factory'}){
    $self->{'feature_factory'} = Bio::EnsEMBL::Analysis::Tools::FeatureFactory
      ->new();
  }
  return $self->{'feature_factory'};
}

#utility methods

=head2 fetch_sequence

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : string, name
  Arg [3]   : Bio::EnsEMBL::DBAdaptor
  Arg [4]   : arrayref of logic_name if sequence is to be masked
  Arg [5]   : Boolean for softmasking if sequence is to be softmasked
  Function  : gets sequence from specifed database
  Returntype: Bio::EnsEMBL::Slice
  Exceptions: none
  Example   : 

=cut


sub fetch_sequence{
  my ($self, $name, $db, $repeat_masking, $soft_masking) = @_;
  if(!$db){
    $db = $self->db;
  }
  if(!$name){
    $name = $self->input_id;
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


=head2 parameters_hash

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : string, parameters string
  Function  : parses the parameters string into a hash
  for the Runnables constructor. If neither of the delimiters
  are found in the string the string is given the key of options
  Returntype: hashref
  Exceptions: 
  Example   : 

=cut


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




=head2 run

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Function  : cycles through all the runnables, calls run and pushes
  their output into the RunnableDBs output array
  Returntype: array ref
  Exceptions: none
  Example   : 

=cut



sub run{
  my ($self) = @_;
  foreach my $runnable(@{$self->runnable}){
    $runnable->run;
    $self->output($runnable->output);
  }
  return $self->{'output'};
}



=head2 write_output

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Function  : set analysis and slice on each feature
  Returntype: 1
  Exceptions: none
  Example   : 

=cut


sub write_output {
  my ($self) = @_;

  my $adaptor  = $self->get_adaptor();
  my $analysis = $self->analysis();

  # Keep track of the analysis_id we have here, because the store()
  # method might change it if the analysis does not already exist in
  # the output database, which will make running the next job difficult
  # or impossible (because the analysis tables weren't in sync).  If we
  # find that the analysis ID changes, throw() so that the genebuilder
  # realizes that she has to sync thhe analysis tables.
  my $analysis_id = $analysis->dbID();

  foreach my $feature ( @{ $self->output() } ) {
    $feature->analysis($analysis);

    if ( !defined( $feature->slice() ) ) {
      $feature->slice( $self->query() );
    }

    $self->feature_factory->validate($feature);

    eval { $adaptor->store($feature); };
    if ($@) {
      throw( sprintf( "RunnableDB::write_output() failed: " .
                        "failed to store '%s' into database '%s': %s",
                      $feature, $adaptor->dbc()->dbname(), $@ ) );
    }
  }

  # Determine if the analysis ID changed, and throw() if it did.
  if ( $analysis->dbID() != $analysis_id ) {
    throw( sprintf( "The analysis ID for '%s' changed from %d to %d. " .
                      "Are the analysis tables in sync?\n",
                    $analysis->logic_name(), $analysis_id,
                    $analysis->dbID() ) );
  }

  return 1;
} ## end sub write_output


=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Function  : throw as it means child hasnt implement an essential method
  Returntype: none
  Exceptions: see function
  Example   : 

=cut



sub fetch_input{
  my ($self) = @_;
  throw("Must implement fetch input in ".$self." RunnableDB will ".
        "not provide this");
}


=head2 read_and_check_config

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : hashref, should be the hashref from which ever config file you are reading 
  Arg [3]   : label - the name of the config varible you're reading - useful for debugging
  Arg [4]   : flag (1/0 ) to throw or not throw if logic_name is missing from config. 
              ( useful for auto-setup in cloud ) 

  Function  : to on the basis of the entries of the hash in your specific
  config file set up instance variables first for the default values then for
  any values specific to you logic name
  Returntype: none
  Exceptions: none
  Example   : 

=cut


sub read_and_check_config{
  my ($self, $var_hash, $label ) = @_; 
 
  if ( defined $label ) { 
    print "READING CONFIG  : $label\n" ; 
  } 
  parse_config($self, $var_hash, $self->analysis->logic_name,$self->no_config_exception);
}


sub read_and_check_config_value{
  my ($self, $var_hash, $label, $values_to_get ) = @_; 
 
  if ( defined $label ) { 
    print "READING CONFIG  : $label\n" ; 
  } 
  parse_config_value($self, $var_hash, $self->analysis->logic_name, $values_to_get);
}

sub read_and_check_config_mini{
  my ($self, $var_hash, $label ) = @_; 
 
  if ( defined $label ) { 
    print "READING CONFIG  : $label\n" ; 
  } 
  parse_config_mini($self, $var_hash); 
}

=head2 require_module

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Blast
  Arg [2]   : string, module path
  Function  : uses perls require to use the past in module
  Returntype: returns module name with / replaced by ::
  Exceptions: throws if require fails
  Example   : my $parser = 
  $self->require('Bio/EnsEMBL/Analysis/Tools/BPliteWrapper');

=cut



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

=head2 ignore_config_file   

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : string  ( 1 or 0 ) 
  Function  : Getter/Setter for value if the configuration file for the module should 
              be ignored or not. Default is to not ignore ( =read / use ) the config file.
  Returntype: 1 or 0 

=cut

sub ignore_config_file {
  my $self = shift;
  $self->{'ignore_config'} = shift if(@_);
  return $self->{'ignore_config'};
} 

=head2 no_config_exception 

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]   : string  ( 1 or 0 ) 
  Function  : 
              
  Returntype: 1 or 0 

=cut

sub no_config_exception{
  my $self = shift;
  $self->{'no_config_exception'} = shift if(@_);
  return $self->{'no_config_exception'};
}



1;
