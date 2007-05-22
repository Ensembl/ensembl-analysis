# Ensembl module for Bio::EnsEMBL::Analysis::Runnable::Funcgen
#
# Copyright (c) 2007 Ensembl
#
=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Fungen

=head1 SYNOPSIS

=head1 DESCRIPTION

This module is the base class for the Fungen Runnables.

=head1 AUTHOR

This module was created by Stefan Graf. It is part 
of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::Runnable::Funcgen;

use strict;
use warnings;
use Data::Dumper;

use Bio::EnsEMBL::Analysis::Runnable;

use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Analysis::Runnable);

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

    my ($probe_features) = rearrange(['PROBE_FEATURES'], @args);
    $self->probe_features($probe_features);
    #print Dumper($self->features);

    #print Dumper $self;
    return $self;

}

=head2 run

    Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::Funcgen::Chipotle
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
  
    #print Dumper $self->options;

    throw("Can't run ".$self." without probe features") 
        unless($self->probe_features);
    #print Dumper $self->probe_features;
   
    $self->workdir($dir) if ($dir);
    $self->checkdir();
    #print "work dir ".$self->workdir()." checked\n";
    
    my $infile = $self->write_infile();
    #print Dumper $infile;
    $self->files_to_delete($self->infile);

    $self->run_analysis();

    #print "Parsing results ... ";
    $self->parse_results();
    #print "done!\n";

    $self->delete_files;
    return 1;
}

=head2 probe_features

  Arg [1]     : Bio::EnsEMBL::Analysis::RunnableDB::Chipotle
  Arg [2]     : arrayref of probe features
  Description : container for probe features
  Returntype  : arrayref
  Exceptions  : throws if no probe feature container is defined
  Example     : 

=cut

sub probe_features {
    my ($self, $features) = @_;
    if($features){
        $self->{'probe_features'} = $features;
    }
    throw("No probe features available in Runnable.") 
        if (!$self->{'probe_features'});
    return $self->{'probe_features'};
}


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

  Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::Chipotle
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
  
  throw("parse_results: can't open file ".$resultsfile.": ". $!)
      unless (open(F, $resultsfile));

  my @output = ();

  while (<F>) {

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
