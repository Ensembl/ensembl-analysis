# Ensembl module for Bio::EnsEMBL::Analysis::Runnable::Funcgen::Nessie
#
# Copyright (c) 2007 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Funcgen::Nessie

=head1 SYNOPSIS

  my $nessie = Bio::EnsEMBL::Analysis::Runnable::Funcgen::Chipotle->new
  (
    -analysis => $analysis,
    -query => 'slice',
    -program => 'nessie',
  );
  $nessie->run;
  my @predicted_features = @{$nessie->output};

=head1 DESCRIPTION

Nessie expects to run the program Nessie and predicts features which 
can be stored in the predicted_feature table in the eFG database

=head1 AUTHOR

This module was created by Stefan Graf. It is part of the 
Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::Runnable::Funcgen::Nessie;

use strict;
use warnings;
use Data::Dumper;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Analysis::Runnable::Funcgen;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::Funcgen);


=head2 new

  Arg         : 
  Usage       : my $runnable = Bio::EnsEMBL::Analysis::Runnable::Nessie->new()
  Description : Instantiates new Chipotle runnable
  Returns     : Bio::EnsEMBL::Analysis::Runnable::Funcgen::Nessie object
  Exceptions  : none

=cut

sub new {

    my ($class,@args) = @_;
    #print Dumper @args;
    my $self = $class->SUPER::new(@args);
    #print Dumper $self;

    my ($features, $options) = 
        rearrange(['FEATURES', 'OPTIONS'], @args);
    $self->probe_features($features);
    #print Dumper($self->features);
    #print Dumper($options);
    $self->options($options);
    #print Dumper($self->options);

    ##################
    #SETTING DEFAULTS#
    ##################
    $self->program('nessie') if (!$self->program);
    ##################

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

    throw("Can't run ".$self." without a probe features") 
        unless($self->probe_features);
 
    #print Dumper $self->probe_features;
   
    $self->workdir($dir) if($dir);
    $self->checkdir();
    #print "work dir checked\n";

    my $infile = $self->infile();
    print Dumper $infile;
    
    foreach my $rset_name (sort keys %{$self->probe_features})
    {

        (my $datafile = $infile) =~ s,\.dat,.$rset_name.dat,;
        @{$self->datafiles($self->write_infile($datafile))};

    }
    print Dumper $self->datafiles;
    $self->write_filelist();

#    (my $resultsfile = $infile) =~ s/\.dat$/\_$alpha\_peaks\.tsv/;
#    $self->resultsfile($resultsfile);
#    #warn("infile:\t".$infile);
#    #warn("resultsfile:\t".$resultsfile);
#    
    $self->files_to_delete($self->infile);
    map {$self->files_to_delete($_)} @{$self->datafiles};
#    #$self->files_to_delete($self->resultsfile);
    $self->run_analysis();
#    #print "Parsing results ... ";
#    $self->parse_results;
#    #print "done!\n";
#    $self->delete_files;
    return 1;
}

=head2 run_analysis

  Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::Chipotle
  Arg [2]     : string, program name
  Usage       : 
  Description : 
  Returns     : 
  Exceptions  : 

=cut

sub run_analysis {
        
    my ($self, $program) = @_;
    
    if(!$program){
        $program = $self->program;
    }
    throw($program." is not executable Chipotle::run_analysis ") 
        unless($program && -x $program);
    
    my $command = $self->program.' --data="'.$self->infile().'" '.$self->options;
    
    warn("Running analysis " . $command . "\n");
    
#    eval { system($command) };
#    throw("FAILED to run $command: ", $@) if ($@);
#    
#    # Check if adjustments via --transform need to be applied
#    
#    (my $logfile = $self->resultsfile()) =~ s/_peaks.tsv/_STDOUT.txt/;
#    open(LOG, "$logfile")
#        or throw("Unable to open logfile $logfile.");
#    my $transform;
#    while (<LOG>) {
#        print;
#        next until (/^Re-run chipotle using \'(--transform -?\d+\.\d+)\'/);
#        $transform = $1;
#        #print "transform: ", $transform, "\n";
#    }
#    close LOG;
#
#    if (defined $transform) {
#        $command .= " $transform";
#        warn("Re-running analysis " . $command . "\n");
#        
#        eval { system($command) };
#        if ($@) {
#            throw("FAILED to run ".$command);
#        }
#    }
}

=head2 write_infile

    Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::Chipotle
    Arg [2]     : filename
    Description : 
    Returntype  : 
    Exceptions  : 
    Example     : 

=cut

sub write_infile {

    my ($self, $filename) = @_;

    if (! $filename) {
        $filename = $self->infile();
    }
    
    open(F, ">".$filename)
        or throw("Can't open file $filename.");

    foreach (keys %{$self->probe_features}) {
        next unless ($filename =~ m/$_/);
        foreach my $ft (@{$self->probe_features->{$_}}) {
            print F join("\t", @$ft), "\n";
        }
        last;
    }

    close F;

    return $filename;

}

=head2 infile

  Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::Chipotle
  Arg [2]     : filename (string)
  Description : will hold a given filename or if one is requested but none
  defined it will use the create_filename method to create a filename
  Returntype  : string, filename
  Exceptions  : none
  Example     : 

=cut


sub infile{

  my ($self, $filename) = @_;

  if($filename){
    $self->{'infile'} = $filename;
  }
  if(!$self->{'infile'}){
    $self->{'infile'} = $self->create_filename('nessie', 'dat');
  }

  return $self->{'infile'};

}

=head2 datafiles

  Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::Chipotle
  Arg [2]     : filename (string)
  Description : will hold/add a list of filename
  Returntype  : listref of filenames
  Exceptions  : none
  Example     : 

=cut


sub datafiles{

  my ($self, $filename) = @_;

  if(! $self->{'datafiles'}){
      $self->{'datafiles'} = [];
  }

  if($filename){
      push @{$self->{'datafiles'}}, $filename;
  }

  return \@{$self->{'datafiles'}};

}

=head2 write_filelist

    Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::Chipotle
    Description : writes the nessie specific wrapper file containing a 
                  list of replicate data file names
    Returntype  : none
    Exceptions  : 
    Example     : 

=cut

sub write_filelist {

    my ($self) = shift;

    throw("No infile found") if (! $self->infile());
    
    open(F, ">".$self->infile())
        or throw("Can't open file $self->infile.");

    foreach my $f (@{$self->datafiles}) {
        print F $f, "\n";
    }

    close F;

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

  if(!$resultsfile){
      $resultsfile = $self->resultsfile;
  }

  throw("parse_results: results file ".$resultsfile." does not exist.")
      if (! -e $resultsfile);
  
  throw("parse_results: can't open file ".$resultsfile.".")
      unless (open(F, $resultsfile));

  my $ff = $self->feature_factory;
  my @output = ();

  while (<F>) {

      chomp;
      my @ft = split;
      push(@output, \@ft);

  }

  $self->output(\@output);
  #print Dumper $self->output();

  throw("parse_results: can't close file ".$resultsfile.".")
      unless (close(F));


}

1;
