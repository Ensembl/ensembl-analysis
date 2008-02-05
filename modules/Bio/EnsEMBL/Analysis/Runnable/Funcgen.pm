# Ensembl module for Bio::EnsEMBL::Analysis::Runnable::Funcgen
#
# Copyright (c) 2007 Ensembl
#
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

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::Runnable::Funcgen;

use strict;
use warnings;
use Data::Dumper;

use Bio::EnsEMBL::Analysis::Runnable;

use Bio::EnsEMBL::Utils::Exception qw( throw warning stack_trace_dump );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

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

    my ($class,@args) = @_;
    #print Dumper @args;
    my $self = $class->SUPER::new(@args);
    #print Dumper $self;

    my ($result_features) = rearrange(['RESULT_FEATURES'], @args);
    $self->result_features($result_features);
    #print Dumper($self->features);

    #print Dumper $self;
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
    print "Bio::EnsEMBL::Analysis::Runnable::Funcgen::run\n";
    my ($self, $dir) = @_;
  
    #print Dumper $self->options;

    #throw("Can't run ".$self." without result features") 
    #    unless($self->result_features);
    #print Dumper $self->probe_features;
   
    #$self->workdir($dir) if ($dir);
    ### DO NOT USE environment variable here!!! ###
    $self->workdir($ENV{ANALYSIS_WORK_DIR});
    $self->checkdir();
    print "work dir ".$self->workdir()." checked\n";

    $self->write_infile();
    
    throw("Input file ".$self->infile." is empty.") if (-z $self->infile);
    #$self->files_to_delete($self->infile);
	
    $self->run_analysis();

    #print "Parsing results ... ";
    $self->parse_results();
    #print "done!\n";

    $self->delete_files;
    return 1;
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

  Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::Funcgen
  Arg [2]     : filename (string)
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
    
    print Dumper $resultsfile;
    
    throw("parse_results: can't open file ".$resultsfile.": ". $!)
        unless (open(F, $resultsfile));
    
    my @output = ();
    
    while (<F>) {
        s/\"//;
        s/^chr//;
        next unless (/^[0-9XYM]+\s/);
        
        chomp;
        my @ft = split;
        push(@output, [ @ft[@{$self->output_fields()}] ]);
        
    }
    
    $self->output(\@output);
    #print Dumper $self->output();
    
    throw("parse_results: can't close file ".$resultsfile.".")
        unless (close(F));
    
}

1;
