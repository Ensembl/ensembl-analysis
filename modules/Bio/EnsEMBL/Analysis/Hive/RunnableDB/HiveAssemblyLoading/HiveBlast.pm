=head1 LICENSE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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

  $self->create_analysis;

  my $input_id = $self->param('iid');
  my $slice = $self->fetch_sequence($input_id,$dba,$repeat_masking);
  $self->query($slice);

  my $parser = $self->make_parser;
  my $filter;
  if($self->BLAST_FILTER){
    $filter = $self->make_filter;
  }

  unless($self->param('blast_db_path')) {
    $self->throw("You did not pass in the blast_db_path parameter. This is required to locate the blast db");
  }

  $self->create_analysis;
  $self->analysis->program($self->param('blast_program')) if ($self->param_is_defined('blast_program'));
  $self->analysis->program_file($self->param('blast_exe_path')) if ($self->param_is_defined('blast_exe_path'));
  $self->analysis->parameters($self->param('commandline_params')) if ($self->param_is_defined('commandline_params'));
  $self->analysis->db_file($self->param('blast_db_path')) if ($self->param_is_defined('blast_db_path'));
  $self->analysis->db($self->param('blast_db_name')) if ($self->param_is_defined('blast_db_name'));

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Blast->new
    (
     -query => $self->query,
     -program => $self->analysis->program_file,
     -parser => $parser,
     -filter => $filter,
     -database => $self->analysis->db_file,
     -analysis => $self->analysis,
     %{$self->BLAST_PARAMS},
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

    if($@) {
      my $except = $@;
      chomp $except;

      if($except =~ /still running after your timer/) {
        $self->warning("BLAST took longer than the timer limit (".$self->param('timer')."), will dataflow input id on branch -2. Exception:\n".$except);
        $self->batch_failed(1);
        $self->param('_branch_to_flow_to_on_fail',-2);
        last;
      } elsif ($except =~ /^\"([A-Z_]{1,40})\"$/i) {
        # only match '"ABC_DEFGH"' and not all possible throws
        my $code = $1;
        # treat VOID errors in a special way; they are really just
        # BLASTs way of saying "won't bother searching because
        # won't find anything"

        if ($code ne 'VOID') {
          $self->failing_job_status($1);
          $self->throw("Blast::run failed ".$@);
        } elsif ($except) {
          $self->throw("Blast Runnable returned unrecognised error string: ".$except);
        }
      } # end elsif ($except =~ /^\"([A-Z_]{1,40})\"$/i)
    } # end if($@)
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


  if($self->batch_failed == 1) {
    # Flow out on -2 or -3 based on how the failure happened
    my $failure_branch_code = $self->param('_branch_to_flow_to_on_fail');
    my $output_hash = {};
    $output_hash->{'iid'} = $self->param('iid');
    $self->dataflow_output_id($output_hash,$failure_branch_code);
  } else {
    my $dna_a = $self->hrdb_get_con('target_db')->get_DnaAlignFeatureAdaptor;
    my $protein_a = $self->hrdb_get_con('target_db')->get_ProteinAlignFeatureAdaptor;
    my $ff = $self->feature_factory;
    foreach my $feature (@{$self->output}) {
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
  } # end else
  return 1;
}


=head2 make_parser

 Arg [1]    : (optional) Hashref
 Description: Load and create a parser object using the name in $self->BLAST_PARSER and using
              parameters in $self->PARSER_PARAMS if not Arg[1] given
 Returntype : Object
 Exceptions : None

=cut

sub make_parser{
  my ($self, $hash) = @_;
  if(!$hash){
    $hash = $self->PARSER_PARAMS;
  }
  my %parser = %$hash;
  $self->require_module($self->BLAST_PARSER);
  my $parser = $self->BLAST_PARSER->new(
                                        %parser
                                       );
  return $parser;
}


=head2 make_filter

 Arg [1]    : (optional) Hashref
 Description: Load and create a parser object using the name in $self->BLAST_FILTER and using
              parameters in $self->FILTER_PARAMS if not Arg[1] given
 Returntype : Object
 Exceptions : None

=cut

sub make_filter{
  my ($self, $hash) = @_;
  if(!$hash){
    $hash = $self->FILTER_PARAMS;
  }
  my %filter = %$hash;
  $self->require_module($self->BLAST_FILTER);
  my $filter = $self->BLAST_FILTER->new(
                                         %filter
                                        );
  return $filter;
}




#config methods

=head2 BLAST_PARSER

 Arg [1]    : (optional) String mmodule name
 Description: Getter/setter
 Returntype : String
 Exceptions : None

=cut

sub BLAST_PARSER{
  my ($self, $value) = @_;

  if(defined $value){
    $self->param('BLAST_PARSER',$value);
  }

  if ($self->param_is_defined('BLAST_PARSER')) {
        return $self->param('BLAST_PARSER');
  }
  else {
      return;
  }
}


=head2 PARSER_PARAMS

 Arg [1]    : (optional) Hashref
 Description: Getter/setter
 Returntype : Hashref
 Exceptions : Throws if the Arg[1] is not a hashref

=cut

sub PARSER_PARAMS{
  my ($self, $value) = @_;

  if(defined $value){
    $self->throw("PARSER_PARAMS must be a hash ref not ".$self->PARSER_PARAMS." Blast::set_hive_config")
      unless(ref($value) eq 'HASH');
    $self->param('PARSER_PARAMS',$value);
  }

  if ($self->param_is_defined('PARSER_PARAMS')) {
        return $self->param('PARSER_PARAMS');
  }
  else {
      return;
  }
}



=head2 BLAST_FILTER

 Arg [1]    : (optional) String module name
 Description: Getter/setter
 Returntype : String
 Exceptions : None

=cut

sub BLAST_FILTER{
  my ($self, $value) = @_;

  if(defined $value){
    $self->param('BLAST_FILTER',$value);
  }

  if ($self->param_is_defined('BLAST_FILTER')) {
        return $self->param('BLAST_FILTER');
  }
  else {
      return;
  }
}



=head2 FILTER_PARAMS

 Arg [1]    : (optional) Hashref
 Description: Getter/setter
 Returntype : Hashref
 Exceptions : Throws if Arg[1] is not a hashref

=cut

sub FILTER_PARAMS{
  my ($self, $value) = @_;

  if(defined $value){
    $self->throw("FILTER_PARAMS must be a hash ref not ".$self->FILTER_PARAMS." Blast::set_hive_config")
      unless (ref($value) eq 'HASH');
    $self->param('FILTER_PARAMS',$value);
  }

  if ($self->param_is_defined('FILTER_PARAMS')) {
        return $self->param('FILTER_PARAMS');
  }
  else {
      return;
  }
}


=head2 BLAST_PARAMS

 Arg [1]    : (optional) Hashref
 Description: Getter/setter, it always call $self->parameters_hash which get parameters
              from $self->analysis->parameters if defined
 Returntype : Hashref
 Exceptions : Throws if Arg[1] is not a hashref

=cut

sub BLAST_PARAMS {
  my ($self, $value) = @_;
  if(defined $value){
    $self->throw("BLAST_PARAMS must be a hash ref not ".$self->BLAST_PARAMS." Blast::set_hive_config")
      unless (ref($value) eq 'HASH');
    $self->param('BLAST_PARAMS',$value);
  }
  if($self->parameters_hash) {
    my %parameters = (%{$self->param('BLAST_PARAMS') || {}}, %{$self->parameters_hash});
    $self->param('BLAST_PARAMS', \%parameters);
  }

  if ($self->param_is_defined('BLAST_PARAMS')) {
        return $self->param('BLAST_PARAMS');
  }
  else {
      return;
  }
}

sub batch_failed {
  my ($self,$batch_failed) = @_;
  if($batch_failed) {
    $self->param('_batch_failed',$batch_failed);
  }
  return ($self->param('_batch_failed'));
}

1;
