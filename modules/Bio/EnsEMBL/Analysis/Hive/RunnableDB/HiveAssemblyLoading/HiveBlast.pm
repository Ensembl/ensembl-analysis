=head1 LICENSE

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

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Analysis::RunnableDB::Blast - 

=head1 SYNOPSIS

  my $blast = Bio::EnsEMBL::Analysis::RunnableDB::Blast->
  new(
      -analysis => $analysis,
      -db => $db,
      -input_id => 'contig::AL1347153.1.3517:1:3571:1'
     );
  $blast->fetch_input;
  $blast->run;
  my @output =@{$blast->output};

=head1 DESCRIPTION

  This module acts as an intermediate between the blast runnable and the
  core database. It reads configuration and uses information from the analysis
  object to setup the blast runnable and then write the results back to the 
  database

=head1 METHODS

=cut

package Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveAssemblyLoading::HiveBlast;

use strict;
use warnings;
use feature 'say';

use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Analysis::Runnable::Blast;

use parent('Bio::EnsEMBL::Analysis::Hive::RunnableDB::HiveBaseRunnableDB');


=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Blast
  Function  : fetch sequence out of database, instantiate the filter, 
  parser and finally the blast runnable
  Returntype: 1
  Exceptions: none
  Example   : 

=cut



sub fetch_input{
  my ($self) = @_;

  my $repeat_masking = $self->param('repeat_masking_logic_names');

  my $dba = $self->hrdb_get_dba($self->param('target_db'));
  $self->hrdb_set_con($dba,'target_db');

  $self->hive_set_config;

  my $input_id = $self->param('iid');
  my $slice = $self->fetch_sequence($input_id,$dba,$repeat_masking);
  $self->query($slice);

  my %blast = %{$self->BLAST_PARAMS};

  my $parser = $self->make_parser;
  my $filter;
  if($self->BLAST_FILTER){
    $filter = $self->make_filter;
  }

  unless($self->param('blast_db_path')) {
    $self->throw("You did not pass in the blast_db_path parameter. This is required to locate the blast db");
  }

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Blast->new
    (
     -query => $self->query,
     -program => $self->analysis->program_file,
     -parser => $parser,
     -filter => $filter,
     -database => $self->analysis()->db_file,
     -analysis => $self->analysis(),
     %blast,
    );
  $self->runnable($runnable);
  return 1;
}


=head2 run

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Blast
  Function  : go through the runnables, run each one and check for errors
  and push the output on to the output array
  Returntype: 1;
  Exceptions: throws if blast run fails
  Example   :

=cut

sub run {
  my ($self) = @_;
  my @runnables = @{$self->runnable};
  foreach my $runnable(@runnables){
    eval{
      $runnable->run;
    };
    if(my $err = $@){
      chomp $err;

      # only match '"ABC_DEFGH"' and not all possible throws
      if ($err =~ /^\"([A-Z_]{1,40})\"$/i) {
        my $code = $1;
        # treat VOID errors in a special way; they are really just
        # BLASTs way of saying "won't bother searching because
        # won't find anything"

        if ($code ne 'VOID') {
          $self->failing_job_status($1);
          $self->throw("Blast::run failed ".$@);
        }
      } elsif ($err) {
        $self->throw("Blast Runnable returned unrecognised error string: ".$err);
      }
    }
    $self->output($runnable->output);
  }

  return 1;
}



=head2 write_output

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Blast
  Function  : get appropriate adaptors and write output to database
  after validating and attaching sequence and analysis objects
  Returntype: undef
  Exceptions: throws if the store fails
  Example   :

=cut


sub write_output{
  my ($self) = @_;
  my $dna_a = $self->hrdb_get_con('target_db')->get_DnaAlignFeatureAdaptor;
  my $protein_a = $self->hrdb_get_con('target_db')->get_ProteinAlignFeatureAdaptor;
  my $ff = $self->feature_factory;
  foreach my $feature (@{$self->output}){
    $feature->analysis($self->analysis);

    unless($feature->slice) {
      $feature->slice($self->query);
    }

    $ff->validate($feature);

    # Put this in to stop any issues with the coordinates of the align feature. Was running into
    # this for some reason with the prediction transcript that had dbIDs as their input id
    my $slice = $feature->slice->seq_region_Slice;
    $feature->slice($slice);

    if($feature->isa('Bio::EnsEMBL::DnaDnaAlignFeature')){
      eval{
        $dna_a->store($feature);
      };
      if($@) {
        $self->throw("Blast:store failed failed to write ".$feature." to the database: ".$@);
      }
    } elsif($feature->isa('Bio::EnsEMBL::DnaPepAlignFeature')){
      eval{
        $protein_a->store($feature);
      };
      if($@) {
        $self->throw("Blast:store failed failed to write ".$feature." to the database: ".$@);
      }
    }
  }
  return 1;
}


=head2 make_parser

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Blast
  Arg [2]   : hashref, parameters for parser constructor
  Function  : create a parser object
  Returntype: a parser object
  Exceptions:
  Example   :

=cut


sub make_parser{
  my ($self, $hash) = @_;
  if(!$hash){
    $hash = $self->PARSER_PARAMS;
  }
  my %parser = %$hash;
  my $parser = $self->BLAST_PARSER->new(
                                        %parser
                                       );
  return $parser;
}

=head2 make_filter

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Blast
  Arg [2]   : hashref, parameters for filter constructor
  Function  : create a filter object
  Returntype: a filter object
  Exceptions:
  Example   :

=cut

sub make_filter{
  my ($self, $hash) = @_;
  if(!$hash){
    $hash = $self->FILTER_PARAMS;
  }
  my %filter = %$hash;
  my $filter = $self->BLAST_FILTER->new(
                                         %filter
                                        );
  return $filter;
}




#config methods

sub hive_set_config {
  my $self = shift;

  # Throw is these aren't present as they should both be defined
  unless($self->param_is_defined('logic_name') && $self->param_is_defined('module')) {
    $self->throw("You must define 'logic_name' and 'module' in the parameters hash of your analysis in the pipeline config file, ".
          "even if they are already defined in the analysis hash itself. This is because the hive will not allow the runnableDB ".
          "to read values of the analysis hash unless they are in the parameters hash. However we need to have a logic name to ".
          "write the genes to and this should also include the module name even if it isn't strictly necessary"
         );
  }

  # Make an analysis object and set it, this will allow the module to write to the output db
  my $analysis = new Bio::EnsEMBL::Analysis(
                                             -logic_name => $self->param('logic_name'),
                                             -module => $self->param('module'),
                                             -program_file => $self->param('blast_exe_path'),
                                             -db_file => $self->param('blast_db_path'),
                                             -parameters => $self->param('commandline_params'),
                                           );
  $self->analysis($analysis);

  # Now loop through all the keys in the parameters hash and set anything that can be set
  my $config_hash = $self->param('config_settings');
  foreach my $config_key (keys(%{$config_hash})) {
    if(defined &$config_key) {
      $self->$config_key($config_hash->{$config_key});
    } else {
      $self->throw("You have a key defined in the config_settings hash (in the analysis hash in the pipeline config) that does ".
            "not have a corresponding getter/setter subroutine. Either remove the key or add the getter/setter. Offending ".
            "key:\n".$config_key
           );
    }
  }

  unless($self->BLAST_PARSER) {
    #must have a parser object and to pass to blast
    $self->throw("BLAST_PARSER must be defined either in the DEFAULT entry or in the hash keyed on ".
                 $self->analysis->logic_name." Blast::read_and_check_config")
  }

  $self->require_module($self->BLAST_PARSER);

  #load the filter module if defined
  if($self->BLAST_FILTER){
    $self->require_module($self->BLAST_FILTER);
  }

  #if any of the object params exist, all are optional they must be hash refs
  if($self->PARSER_PARAMS && ref($self->PARSER_PARAMS) ne 'HASH') {
    $self->throw("PARSER_PARAMS must be a hash ref not ".$self->PARSER_PARAMS." Blast::set_hive_config");
  }

  if($self->FILTER_PARAMS && ref($self->FILTER_PARAMS) ne 'HASH') {
    $self->throw("FILTER_PARAMS must be a hash ref not ".$self->FILTER_PARAMS." Blast::set_hive_config");
  }

  if($self->BLAST_PARAMS && ref($self->BLAST_PARAMS) ne 'HASH') {
    $self->throw("BLAST_PARAMS must be a hash ref not ".$self->BLAST_PARAMS." Blast::set_hive_config");
  }

  my $blast_params;
  if($self->BLAST_PARAMS){
    $blast_params = $self->BLAST_PARAMS;
  }else{
    $blast_params = {};
  }

  my %parameters;
  if($self->parameters_hash){
    %parameters = %{$self->parameters_hash};
  }
  foreach my $key(%parameters){
    $blast_params->{$key} = $parameters{$key};
  }
  $self->BLAST_PARAMS($blast_params);

}


sub BLAST_PARSER{
  my ($self, $value) = @_;

  if(defined $value){
    $self->param('_CONFIG_BLAST_PARSER',$value);
  }

  if ($self->param_is_defined('_CONFIG_BLAST_PARSER')) {
        return $self->param('_CONFIG_BLAST_PARSER');
  }
  else {
      return;
  }
}


sub PARSER_PARAMS{
  my ($self, $value) = @_;

  if(defined $value){
    $self->param('_CONFIG_PARSER_PARAMS',$value);
  }

  if ($self->param_is_defined('_CONFIG_PARSER_PARAMS')) {
        return $self->param('_CONFIG_PARSER_PARAMS');
  }
  else {
      return;
  }
}



sub BLAST_FILTER{
  my ($self, $value) = @_;

  if(defined $value){
    $self->param('_CONFIG_BLAST_FILTER',$value);
  }

  if ($self->param_is_defined('_CONFIG_BLAST_FILTER')) {
        return $self->param('_CONFIG_BLAST_FILTER');
  }
  else {
      return;
  }
}



sub FILTER_PARAMS{
  my ($self, $value) = @_;

  if(defined $value){
    $self->param('_CONFIG_FILTER_PARAMS',$value);
  }

  if ($self->param_is_defined('_CONFIG_FILTER_PARAMS')) {
        return $self->param('_CONFIG_FILTER_PARAMS');
  }
  else {
      return;
  }
}


sub BLAST_PARAMS {
  my ($self, $value) = @_;

  if(defined $value){
    $self->param('_CONFIG_BLAST_PARAMS',$value);
  }

  if ($self->param_is_defined('_CONFIG_BLAST_PARAMS')) {
        return $self->param('_CONFIG_BLAST_PARAMS');
  }
  else {
      return;
  }
}

1;
