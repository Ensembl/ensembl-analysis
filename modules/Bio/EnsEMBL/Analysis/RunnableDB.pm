# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB
#
# Copyright (c) 2004 Ensembl
#

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

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut
package Bio::EnsEMBL::Analysis::RunnableDB;


use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Root;
use Bio::EnsEMBL::Analysis::Tools::FeatureFactory;

use vars qw (@ISA);

@ISA = qw(Bio::EnsEMBL::Root);


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
  my $self = $class->SUPER::new(@args);
  &verbose('WARNING');
  my ($db, $input_id, $analysis) = rearrange(['DB', 'INPUT_ID', 
                                              'ANALYSIS'], @args);
  if(!$db || !$analysis || !$input_id){
    throw("Can't create a RunnableDB without a dbadaptor ".$db." an ".
          "analysis object ".$analysis." or an input_id ".$input_id);
  }
  $self->db($db);
  $self->analysis($analysis);
  $self->input_id($input_id);
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
    throw("Must pass RunnableDB:db a Bio::EnsEMBL::DBSQL::DBConnection ".
          "not a ".$db) 
      unless($db->isa('Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor'));
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


=head2 containers

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
  Function  : gets sequence from specifed database
  Returntype: Bio::EnsEMBL::Slice
  Exceptions: none
  Example   : 

=cut


sub fetch_sequence{
  my ($self, $name, $db, $repeat_masking) = @_;
  if(!$db){
    $db = $self->db;
  }
  if(!$name){
    $name = $self->input_id;
  }
  my $sa = $db->get_SliceAdaptor;
  print "Fetching ".$name." from ".$db->dbname."\n";
  my $slice = $sa->fetch_by_name($name);
  $repeat_masking = [] unless($repeat_masking);
  if(@$repeat_masking){
    my $sequence = $slice->get_repeatmasked_seq($repeat_masking);
    $slice = $sequence
  }
  if(!$slice){
    throw("Failed to fetch slice ".$name);
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
  if(!$string){
    return;
  }
  if($string =~  /,/ || $string =~ /=>/){
     my @pairs = split (/,/, $string);
     foreach my $pair(@pairs){
       my ($key, $value) = split (/=>/, $pair); 
       if ($key && $value) {
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
    $parameters_hash{'options'} = $string;
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
    my $output =  $runnable->output;
    push(@{$self->{'output'}}, @{$runnable->output});
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


sub write_output{
  my ($self) = @_;
  my $adaptor = $self->get_adaptor;
  foreach my $feature(@{$self->output}){
    $feature->analysis($self->analysis);
    $feature->slice($self->query) if(!$feature->slice);
    $self->feature_factory->validate($feature);
    eval{
      $adaptor->store($feature);
    };
    if($@){
      throw("RunnableDB:store failed, failed to write ".$feature." to ".
            "the database $@");
    }
  }
  return 1;
}






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





1;
