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

package Bio::EnsEMBL::Analysis::RunnableDB::Blast;

use strict;
use warnings;

use Bio::EnsEMBL::Analysis::RunnableDB;
use Bio::EnsEMBL::Analysis::Runnable::Blast;
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Analysis::Config::General;
use Bio::EnsEMBL::Analysis::Config::Blast;
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::RunnableDB);


sub new {
  my ($class,@args) = @_;
  my $self = $class->SUPER::new(@args);

  $self->read_and_check_config($BLAST_CONFIG_BY_LOGIC);

  return $self;
}



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

  my $slice = $self->fetch_sequence($self->input_id, $self->db, 
                                    $ANALYSIS_REPEAT_MASKING);
  $self->query($slice);
  my %blast = %{$self->BLAST_PARAMS};

  my $parser = $self->make_parser;
  my $filter;
  if($self->BLAST_FILTER){
    $filter = $self->make_filter;
  }

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Blast->new
    (
     -query => $self->query,
     -program => $self->analysis->program_file,
     -parser => $parser,
     -filter => $filter,
     -database => $self->analysis->db_file,
     -analysis => $self->analysis,
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
          throw("Blast::run failed $@");
        }
      } elsif ($err) {
        throw("Blast Runnable returned unrecognised error string: $err");
      }
    }
    $self->output($runnable->output);
  }
  1;
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
  my $dna_a = $self->db->get_DnaAlignFeatureAdaptor;
  my $protein_a = $self->db->get_ProteinAlignFeatureAdaptor;
  my $ff = $self->feature_factory;
  foreach my $f(@{$self->output}){
    $f->analysis($self->analysis);
    $f->slice($self->query) if(!$f->slice);
    $ff->validate($f);
    if($f->isa('Bio::EnsEMBL::DnaDnaAlignFeature')){
      eval{
        $dna_a->store($f);
      };
      throw("Blast:store failed failed to write ".$f." to the database ".
            "$@") if($@);
    }elsif($f->isa('Bio::EnsEMBL::DnaPepAlignFeature')){
      eval{
        $protein_a->store($f);
      };
      throw("Blast:store failed failed to write ".$f." to the database ".
            "$@") if($@);
    }
  }
  return ;
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


sub read_and_check_config{
  my ($self, $var_hash) = @_;

  #call RunnableDBs method to fill in values first
  $self->SUPER::read_and_check_config($var_hash);

  #now for type checking and other sanity checks

  #must have a parser object and to pass to blast
  throw("BLAST_PARSER must be defined either in the DEFAULT entry or in ".
        "the hash keyed on ".$self->analysis->logic_name.
        " Blast::read_and_check_config") if(!$self->BLAST_PARSER);
  $self->require_module($self->BLAST_PARSER);

  #load the filter module if defined
  if($self->BLAST_FILTER){
    $self->require_module($self->BLAST_FILTER);
  }
  #if any of the object params exist, all are optional they must be hash 
  #refs
  throw("PARSER_PARAMS must be a hash ref not ".$self->PARSER_PARAMS.
        " Blast::read_and_check_config") 
    if($self->PARSER_PARAMS && ref($self->PARSER_PARAMS) ne 'HASH');
  throw("FILTER_PARAMS must be a hash ref not ".$self->FILTER_PARAMS.
        " Blast::read_and_check_config")
    if($self->FILTER_PARAMS && ref($self->FILTER_PARAMS) ne 'HASH');
  throw("BLAST_PARAMS must be a hash ref not ".$self->BLAST_PARAMS.
        " Blast::read_and_check_config")
    if($self->BLAST_PARAMS && ref($self->BLAST_PARAMS) ne 'HASH');

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
    $self->{'_CONFIG_BLAST_PARSER'} = $value;
  }

  return $self->{'_CONFIG_BLAST_PARSER'};
}


sub PARSER_PARAMS{
  my ($self, $value) = @_;

  if(defined $value){
    $self->{'_CONFIG_PARSER_PARAMS'} = $value;
  }

  return $self->{'_CONFIG_PARSER_PARAMS'};
}



sub BLAST_FILTER{
  my ($self, $value) = @_;

  if(defined $value){
    $self->{'_CONFIG_BLAST_FILTER'} = $value;
  }

  return $self->{'_CONFIG_BLAST_FILTER'};
}



sub FILTER_PARAMS{
  my ($self, $value) = @_;

  if(defined $value){
    $self->{'_CONFIG_FILTER_PARAMS'} = $value;
  }

  return $self->{'_CONFIG_FILTER_PARAMS'};
}


sub BLAST_PARAMS{
  my ($self, $value) = @_;

  if(defined $value){
    $self->{'_CONFIG_BLAST_PARAMS'} = $value;
  }

  return $self->{'_CONFIG_BLAST_PARAMS'};
}





1;
