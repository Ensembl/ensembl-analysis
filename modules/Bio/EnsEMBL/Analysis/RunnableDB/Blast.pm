# Ensembl module for Bio::EnsEMBL::Analysis::RunnableDB::Blast
#
# Copyright (c) 2004 Ensembl
#

=head1 NAME

  Bio::EnsEMBL::Analysis::RunnableDB::Blast

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

  This module is a wrapper for running blast. It knows how to construct
  the commandline and can call to other modules to run the parsing and 
  filtering. By default is constructs wublast commandlines but it can be
  told to construct ncbi command lines. It needs to be passed a Bio::Seq
  and a database name (this database should either have its full path 
  given or it should live in the location specified by the $BLASTDB 
  environment variable). It should also be given a parser object which has
  the method parse_file which takes a filename and returns an arrayref of
  results and optionally it can be given a filter object which has the 
  method filter_results which takes an arrayref of results and returns the
  filtered set of results as an arrayref. For examples of both parser
  objects and a filter object look in Bio::EnsEMBL::Analysis::Tools for
  BPliteWrapper, FilterBPlite and FeatureFilter

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

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



=head2 fetch_input

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Blast
  Function  : fetch sequence out of database, read config files
  instantiate the filter, parser and finally the blast runnable
  Returntype: 1
  Exceptions: none
  Example   : 

=cut



sub fetch_input{
  my ($self) = @_;

  $self->setup_hashes;
  my $slice = $self->fetch_sequence($self->input_id, $self->db, 
                                    $ANALYSIS_REPEAT_MASKING);
  $self->query($slice);
  my %blast = %{$self->blast_hash};

  my $parser = $self->make_parser;
  my $filter;
  if($self->filter_object){
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



=head2 containers

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Blast
  Arg [2]   : hashref/object
  Function  : container for hashrefs or object, this docs are for the
  5 methods, blast_hash, parser_hash, filter_hash, parser_object and
  filter object
  Returntype: hashref/object
  Exceptions: 
  Example   : 

=cut


sub blast_hash{
  my $self = shift;
  $self->{'blast_hash'} = shift if(@_);
  return $self->{'blast_hash'};
}

sub parser_hash{
  my $self = shift;
  $self->{'parser_hash'} = shift if(@_);
  return $self->{'parser_hash'};
}

sub filter_hash{
  my $self = shift;
  $self->{'filter_hash'} = shift if(@_);
  return $self->{'filter_hash'};
}

sub parser_object{
  my $self = shift;
  $self->{'parser_object'} = shift if(@_);
  return $self->{'parser_object'};
}


sub filter_object{
  my $self = shift;
  $self->{'filter_object'} = shift if(@_);
  return $self->{'filter_object'};
}


=head2 setup_hashes

  Arg [1]   : Bio::EnsEMBL::Analysis::RunnableDB::Blast
  Function  : parser the array of hashes from the config file to
  produce a parameter has for the blast, parser and filter object
  Returntype: none
  Exceptions: throws if do not find the analysis logic_name in the hash
  Example   : 

=cut



sub setup_hashes{
  my ($self) = @_;
  my %blast;
  my %parser;
  my %filter;
  my $found = 0;
  LOGIC_NAME:foreach my $hash(@$BLAST_CONFIG){
      if($hash->{logic_name} ne $self->analysis->logic_name){
        next LOGIC_NAME;
      }
      $found = 1;
      while(my ($k,$v) = each(%$hash)){
        if($k eq 'blast_parser'){
          $self->parser_object($v);
          $self->require_module($v);
        }elsif($k eq 'blast_filter'){
          $self->filter_object($v);
          $self->require_module($v);
        }elsif($k eq 'parser_params'){
          %parser = %$v;
        }elsif($k eq 'filter_params'){
          %filter = %$v;
        }elsif($k eq 'blast_params'){
          %blast = %$v;
        }else{
          throw("Don't recognise $k from ".
                "Bio::EnsEMBL::Analysis::Config::Blast::BLAST_CONFIG ")
            unless($k eq 'logic_name');
        }
      }
    }
  if($found == 0){
    throw("Failed to find info for ".$self->analysis->logic_name." in ".
          "Bio::EnsEMBL::Analysis::Config::Blast::BLAST_CONFIG ");
  }
  my %parameters;
  if($self->parameters_hash){
    %parameters = %{$self->parameters_hash};
  }
  foreach my $key(keys(%parameters)){
    $blast{$key} = $parameters{$key};
  }
  $self->blast_hash(\%blast);
  $self->parser_hash(\%parser);
  $self->filter_hash(\%filter);
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
      $self->failing_job_status($1) 
        if $err =~ /^\"([A-Z_]{1,40})\"$/i; 
      # only match '"ABC_DEFGH"' and not all possible throws
      throw("Blast::run failed $@");
    }
    my @output = @{$runnable->output};
    $self->output(\@output);
  }
  1;
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
    $hash = $self->parser_hash;
  }
  my %parser = %$hash;
  my $parser = $self->parser_object->new(
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
    $hash = $self->filter_hash;
  }
  my %filter = %$hash;
  my $filter = $self->filter_object->new(
                                         %filter
                                        );
  return $filter;
}
1;
