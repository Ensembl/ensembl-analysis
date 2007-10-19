# Ensembl module for Bio::EnsEMBL::Analysis::Runnable::Funcgen::TileMap
#
# Copyright (c) 2007 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Funcgen::TileMap

=head1 SYNOPSIS

  my $chipotle = Bio::EnsEMBL::Analysis::Runnable::Funcgen::TileMap->new
  (
    -analysis => $analysis,
    -query => 'slice',
    -program => 'chipotle.pl',
  );
  $chipotle->run;
  my @annotated_features = @{$chipotle->output};

=head1 DESCRIPTION

TileMap expects to run the program ChIPOTle (PMID: 16277752) and 
predicts features which can be stored in the predicted_feature table 
in the eFG database

=head1 AUTHOR

This module was created by Stefan Graf. It is part of the 
Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

package Bio::EnsEMBL::Analysis::Runnable::Funcgen::TileMap;

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

    Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::Funcgen::TileMap
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

  Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::TileMap
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
    throw($program." is not executable TileMap::run_analysis ") 
        unless($program && -x $program);

    my $command = $self->program . ' ' . $self->config_file();
    
    warn("Running analysis " . $command . "\n");
    
    eval { system($command) };
    throw("FAILED to run $command: ", $@) if ($@);

}

=head2 write_infile

    Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::TileMap
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
	print F join("\t", 'chromosome', 'position', 'array ids'), "\n";
    #print Dumper $self->result_features;
    foreach (values %{$self->result_features}) {
        foreach my $ft (@$_) {
            #print join("\t", (@$ft)[0..1,3]), "\n";
            print F join("\t", (@$ft)[0..1,3]), "\n";
        }
    }
    close F;

	# write config and cmpinfo files from templates
	(my $project = $filename) =~ s,.+/(.+)\.dat,$1,;

	my $workdir = $self->workdir;

	my $config_template = $self->analysis->parameters();
	#print Dumper $config_template;

	my $config = $workdir.'/'.$project.'_args.txt';
	#print Dumper $config;
	$self->config_file($config);

	open(IN, $config_template)
		or throw("Can't open config file $config_template");

	open(OUT, ">$config")
		or throw("Can't open config file $config");

	map { 
		s,^(O.1-.+=).+$,$1 $workdir,;
		s,^(O.2-.+=).+$,$1 $project,;
		s,^(I.2-.+=).+$,$1 $project.dat,;
		s,^(IV.5-.+=).+$,$1 ${project}_f_pb.sum,;
		print OUT;
	} <IN>;

	close IN;
	close OUT;

	(my $cmpinfo_template = $config_template) =~ s,_arg.txt,.cmpinfo,;
	#print Dumper $cmpinfo_template;
	my $cmpinfo = $self->workdir.'/'.$project.'.cmpinfo';
	system("cp $cmpinfo_template $cmpinfo");

	# UCSC *.bed file to report significant regions. 
	# Regions are sorted according to their genomic locations.
	my $results = $self->workdir.'/'.$project.'_hmm.bed';
	# tab-delimited file to report significant regions. Regions 
	# are ranked according to their significance levels.
	#my $results = $self->workdir.'/'.$project.'_hmm.reg';
	$self->resultsfile($results);

    # set columns (fields) for output
    my @fields = (0..2,4); # bed
    #my @fields = (0..2,6); # reg
    $self->output_fields(\@fields);

    return $filename;

}

=head2 infile

  Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::TileMap
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
    $self->{'infile'} = $self->create_filename($self->analysis->logic_name, 'dat');
  }

  return $self->{'infile'};

}


sub config_file {
    my $self = shift;
    $self->{'config_file'} = shift if(@_);
    return $self->{'config_file'};
}

1;
