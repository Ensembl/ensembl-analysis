# Ensembl module for Bio::EnsEMBL::Analysis::Runnable::Funcgen::Chipotle
#
# Copyright (c) 2007 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Funcgen::Chipotle

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Funcgen::Chipotle->new
  (
    -analysis => $analysis,
    -query => 'slice',
    -program => 'program.pl',
  );
  $runnable->run;
  my @annotated_features = @{$runnable->output};

=head1 DESCRIPTION

Chipotle expects to run the program ChIPOTle (PMID: 16277752) and 
predicts features which can be stored in the predicted_feature table 
in the eFG database

=head1 AUTHOR

This module was created by Stefan Graf. It is part of the 
Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::Runnable::Funcgen::Chipotle;

use strict;
use warnings;
use Data::Dumper;

use Bio::EnsEMBL::Analysis::Runnable;
use Bio::EnsEMBL::Analysis::Runnable::Funcgen;

use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Analysis::Runnable::Funcgen);

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

=head2 run_analysis

  Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::Chipotle
  Arg [2]     : string, program name
  Usage       : 
  Description : 
  Returns     : 
  Exceptions  : 

=cut

sub run_analysis {
        
    print "Bio::EnsEMBL::Analysis::Runnable::Funcgen::Chipotle\n";
    my ($self, $program) = @_;
    
    if(!$program){
        $program = $self->program;
    }
    throw($program." is not executable Chipotle::run_analysis ") 
        unless($program && -x $program);


	#default parameters
	my %param = (
		'alpha' => 0.05,
		'stepSize' => 50,
		'windowSize' => 400,
		# Benjamini Hochberg
		'adjustPvalue' => 'BH',  
		# Bonferroni 
		#'adjustPvalue' => 'BON',
	);
    #print Dumper $self->analysis->parameters();

	foreach my $param (split(/,\s+/,$self->analysis->parameters())) {
		
		if ($param =~ m/-(.+)=>(.+)/) {
			
			#print Dumper $1;
			throw("Parameter $param is not valid.") if ! exists $param{$1};

			$param{$1} = $2;
			
		} else {

			throw("Parameter $param has not the correct format.");

		}

	}
	
    #if ($self->analysis->parameters =~ m/--alpha (\d+\.\d+)/) {
    #    $alpha = $1;
    #}
    #print $alpha;

    (my $resultsfile = $self->infile) =~ s/\.dat$/.peaks/;
    $self->resultsfile($resultsfile);
    $self->files_to_delete($resultsfile);

    my @fields = (0..3);
    $self->output_fields(\@fields);

    my $command = join(' ', $self->program, 'get_peaks', '0', '0', 
					   $self->infile(), $self->resultsfile(), 
					   $param{alpha}, $param{stepSize},
					   $param{windowSize}, $param{adjustPvalue});
    
    warn("Running analysis " . $command . "\n");
    
    eval { system($command) };
    throw("FAILED to run $command: ", $@) if ($@);

    # Check if adjustments via --transform need to be applied
	### Version 1 logfile processing disabled ###
    
    #(my $logfile = $self->resultsfile()) =~ s/_peaks.tsv/_STDOUT.txt/;
    #open(LOG, "$logfile")
    #    or throw("Unable to open logfile $logfile.");
    #my $transform;
    #while (<LOG>) {
    #    print;
    #    next until (/^Re-run chipotle using \'(--transform -?\d+\.\d+)\'/);
    #    $transform = $1;
    #    #print "transform: ", $transform, "\n";
    #}
    #close LOG;

    #if (defined $transform) {
    #    $command .= " $transform";
    #    warn("Re-running analysis " . $command . "\n");
    #    
    #    eval { system($command) };
    #    if ($@) {
    #        throw("FAILED to run ".$command);
    #    }
    #}
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

    #print Dumper $self->result_features;
    foreach (values %{$self->result_features}) {
        foreach my $ft (@$_) {
            #print join("\t", (@$ft)[4,0..3]), "\n";
            print F join("\t", (@$ft)[4,0..3]), "\n";
        }
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
    $self->{'infile'} = $self->create_filename('chipotle', 'dat');
  }

  return $self->{'infile'};

}

1;
