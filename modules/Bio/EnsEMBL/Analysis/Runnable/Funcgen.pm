# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2024] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Runnable::Fungen

=head1 SYNOPSIS

=head1 DESCRIPTION

This module is the base class for Fungen Runnables.

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Stefan Graf, Ensembl Functional Genomics - http://www.ensembl.org

=head1 CONTACT

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

=cut

package Bio::EnsEMBL::Analysis::Runnable::Funcgen;

use strict;
use warnings;
use Data::Dumper;

use Bio::EnsEMBL::Analysis::Runnable;

use Bio::EnsEMBL::Utils::Exception qw( throw warning stack_trace_dump );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (get_file_format);


use vars qw( @ISA );
@ISA = qw( Bio::EnsEMBL::Analysis::Runnable );

=head2 new

  Arg         : 
  Usage       : my $runnable = Bio::EnsEMBL::Analysis::Runnable::Funcgen->new()
  Description : Instantiates new Chipotle runnable
  Returns     : Bio::EnsEMBL::Analysis::Runnable::Funcgen::Nessie object
  Exceptions  : none

=cut

sub new {

    print "Analysis::Runnable::Funcgen::new\n";

    my ($class,@args) = @_;
    #print Dumper @args;
    my $self = $class->SUPER::new(@args);
    #warn $self;

    my ($result_features, $options, $workdir) = rearrange
        (['RESULT_FEATURES', 'OPTIONS', 'WORKDIR'], @args);

    $self->result_features($result_features);
    $self->workdir($workdir);
    $self->checkdir();

    #warn('RESULT_FEATURES: '.$self->result_features);

    #warn('OPTIONS:  '.$self->options);
    #warn('OPTIONS:  '.$options);
    #warn('A_PARAMS: '.$self->analysis->parameters);

    return $self;

}

=head2 run

    Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::Funcgen
    Arg [2]     : string, directory
    Description : Class specific run method. This checks the directory specifed
    to run it, write the data to file (tab-delimited), marks the query infile 
    and results file for deletion, runs the analysis, parses the results and 
    deletes any files
    Returntype  : 1
    Exceptions  : throws if no query sequence is specified
    Example     : 

=cut

sub run {
    my ($self, $dir) = @_;
     
    $self->write_infile();
    print "Infile:\t".$self->infile."\n";
    throw("Input file ".$self->infile." is empty.") if (-z $self->infile);
    $self->files_to_delete($self->infile);
	
    $self->run_analysis();
    $self->parse_results();
    return 1;
}

=head2 get_parameters

  Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::Funcgen
  Description : parse analysis parameters
  Returntype  : hash ref with parameter names as keys
  Exceptions  : none
  Example     : 

=cut


sub get_parameters {

    my ($self) = @_;
    my %parameters = ();
    
    throw ("Object needs to be a Bio::EnsEMBL::Analysis::Runnable::Funcgen") 
        if (! $self->isa("Bio::EnsEMBL::Analysis::Runnable::Funcgen"));

    my @parameters = split(/\s+/, $self->analysis->parameters());
    map { throw("Parameter $_ has not the correct format.") 
              if (! m/^(.+)=(.+)/);
          $parameters{$1} = $2; } @parameters;

    #print Dumper %parameters;

    return \%parameters;

}


=head2 infile

  Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::Funcgen
  Arg [2]     : filename (string)
  Description : will hold a given filename or if one is requested but none
  defined it will use the create_filename method to create a filename
  Returntype  : string, filename
  Exceptions  : none
  Example     : 

=cut


sub infile{

  my ($self, $filename) = @_;


  #This is a bit random!
  warn "Can we not name the files after the feature_sets?";

  if($filename){
    $self->{'infile'} = $filename;
	#Need to checkdir here?


  }
  if(!$self->{'infile'}){
    $self->{'infile'} = $self->create_filename
        ($self->analysis->module.'/'.$self->analysis->module, 'dat');
    (my $dir = $self->{'infile'}) =~ s,/[^/]+$,,;
    system("mkdir -p $dir") unless ( -d $dir );
  }

  return $self->{'infile'};

}

=head2 result_features

  Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::Funcgen
  Arg [2]     : arrayref of result features
  Description : container for result features
  Returntype  : arrayref
  Exceptions  : throws if no probe feature container is defined
  Example     : 

=cut

sub result_features {
    my ($self, $features) = @_;
    if($features){
        $self->{'result_features'} = $features;
    }
    #throw("No result features available in Runnable.") 
    #    if (!$self->{'result_features'});
    return $self->{'result_features'};
}

=head2 output_fields

  Arg [1]     : Bio::EnsEMBL::Analysis::RunnableDB
  Arg [2]     : arrayref of output fields
  Description : container for outout fields
  Returntype  : arrayref
  Exceptions  : throws if no output field container is defined
  Example     : 

=cut

sub output_fields {
    my ($self, $array) = @_;
    #print Dumper $array;
    if($array){
        $self->{'outfile_fields'} = $array;
    }
    throw("No outfile field string set in Runnable.") 
        if (!$self->{'outfile_fields'});
    return $self->{'outfile_fields'};
}


=head2 parse_results

  Arg [1]     : filename (string)
  Decription  : open and parse resultsfile
  Returntype  : none
  Exceptions  : throws 
  Example     : 

=cut

sub parse_results{
    my ($self, $resultsfile) = @_;
    
    if (! defined $resultsfile) {
        $resultsfile = $self->resultsfile;
    }
    throw ('No resultsfile defined in instance!')
        if(! defined $resultsfile);
    
    throw("parse_results: results file ".$resultsfile." does not exist.")
        if (! -e $resultsfile);
    
    warn('Resultsfile: '.$resultsfile);
    
    throw("parse_results: can't open file ".$resultsfile.": ". $!)
        unless (open(F, $resultsfile));
    
    my @output = ();
    
	#What format is this generic parse method expecting?
	#Getting failure for SWEmbl here as this expects seq_region_name, not slice name!
	#Which for SWEmbl is entirely dpendant on on the input
	#bed -> seq_region_name
	#sam -> slice name

	#Will this cause failure of other Runnables?
	
	#These regexs should be put in file parsers!
	my $regex;
	#warn "input_format to ".$self->{'input_format'}." $self";;
	#This is not the same instance of the Runnable which submitted the job
	#hence, we won't have set the input_format
	
	my $file_format = &get_file_format($self->query);
	warn "got file format $file_format from ".$self->query;


	#Changed seq_region match to handle all words here rather than just [0-9XYM]

	if($file_format eq 'sam'){
	  $regex = '^(\w+):[A-Za-z0-9]+:(\w+):[0-9]+:[0-9]+:[01]\s';
	}
	elsif($file_format eq 'bed'){
	  #defaulting to seq_region name
	  $regex = '^(chr|)(\w+)\s';
	}
	else{#Should have already caught this during submission
	  throw("$file_format file format not supported");
	}


    while (<F>) {
	  s/\"//;
	  next if($_ !~ /$regex/o);
        
	  chomp;
	  my @ft = split;
	  #This would be quicker doing string selwction directly from the regex match(as with seq_region_name)
	  push(@output, [ $2, @ft[@{$self->output_fields()}] ]);
        
    }
    
    $self->output(\@output);
    #print Dumper $self->output();
    
    throw("parse_results: can't close file ".$resultsfile.".")
        unless (close(F));
    
}

sub query {

    my ($self, $query) = @_;

    if ( $query ) {
        
        throw("Must pass RunnableDB:Funcgen:query a array ref not a ".
              ref($query)) unless (ref($query) eq 'ARRAY');

        map { 
            throw($_->name . " is not a Bio::EnsEMBL::Slice")
                unless ($_->isa("Bio::EnsEMBL::Slice"));
        } @$query;

        $self->{'query'} = $query;
    }

    return $self->{'query'};

}

1;
