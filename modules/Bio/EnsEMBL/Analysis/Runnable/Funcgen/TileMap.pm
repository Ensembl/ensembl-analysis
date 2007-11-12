# Ensembl module for Bio::EnsEMBL::Analysis::Runnable::Funcgen::TileMap
#
# Copyright (c) 2007 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Funcgen::TileMap

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Funcgen::TileMap->new
  (
    -analysis => $analysis,
    -query => 'slice',
    -program => 'program.pl',
  );
  $runnable->run;
  my @features = @{$runnable->output};

=head1 DESCRIPTION

TileMap expects to run the program TileMap (Ji and Wong (2005), PMID: 16046496) 
and predicts features which can be stored in the annotated_feature table 
in the eFG database

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Stefan Graf, Ensembl Functional Genomics - http://www.ensembl.org

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

	# determine both number of result sets (replicates) and features

	my $noa = scalar(keys %{$self->result_features});
	warn('no. of result sets: ', $noa);

	my $nof = scalar(@{(values %{$self->result_features})[0]});
	warn('no. of result fts: ', $nof);

	# dump features
	open(F, ">".$filename)
		or throw("Can't open file $filename.");
	print F join("\t", 'chromosome', 'position', 'array ids'), "\n";

	#print Dumper $self->result_features;
	
	for (my $i=0; $i<$nof; $i++) {
		my $coord=0;
		foreach my $rset (values %{$self->result_features}) {
			print F join("\t", ${$rset}[$i][0], ${$rset}[$i][1]) unless ($coord);
			$coord=1;    
			print F "\t".${$rset}[$i][3];
		}
		print F "\t0" if ($noa == 1);
		print F "\n";
	}
	close F;

	my $workdir = $self->workdir;

	# write config and cmpinfo files from templates
	(my $project = $filename) =~ s,.+/(.+)\.dat,$1,;

	my $config_template = $self->analysis->parameters();
	#print Dumper $config_template;

	my $config = $workdir.'/'.$project.'_args.txt';
	#print Dumper $config;
	$self->config_file($config);

	# program parameters
	my $method	= exists $ENV{METHOD}? $ENV{METHOD} : 0;
	my $postprob = exists $ENV{POSTPROB}? $ENV{POSTPROB} : 0.5;
	my $maxgap = exists $ENV{MAXGAP}? $ENV{MAXGAP} : 1000;
	my $hyblength = exists $ENV{HYBLENGTH}? $ENV{HYBLENGTH} : 28;

	open(IN, $config_template)
		or throw("Can't open config file $config_template");

	open(OUT, ">$config")
		or throw("Can't open config file $config");

	map { 
		s,^(O.1-.+=).+$,$1 $workdir,;     #[Working directory]
		s,^(O.2-.+=).+$,$1 $project,;     #[Project Title]
		s,^(I.2-.+=).+$,$1 $project.dat,; #[Raw data file]
		s,^(I.3-.+=).+$,$1 2,;            #[Range of test-statistics] (0: default; 1: [0,1], 2: (-inf, +inf))
		s,^(II.1-.+=).+$,$1 0,;           #[Apply local repeat filter?] (0:No; 1:Yes)
		s,^(II.2-.+=).+$,$1 NULL,;        #[*.refmask file]
		s,^(III.2-.+=).+$,$1 $method,;    #[Method to combine neighboring probes] (0:HMM, 1:MA)
		s,^(IV.1-.+=).+$,$1 $postprob,;   #[Posterior probability >]
		s,^(IV.2-.+=).+$,$1 $maxgap,;     #[Maximal gap allowed] (1000: default)
		s,^(IV.4-.+=).+$,$1 0,;           #[Provide your own selection statistics?] (0: No, use default; 1: Yes)
		s,^(IV.5-.+=).+$,$1 NULL,;        #[If Yes to IV.4, selection statistics file]
		s,^(IV.10-.+=).+$,$1 $hyblength,; #[Expected hybridization length]
		s,^(V.2-.+=).+$,$1 $maxgap,;      #[Maximal gap allowed] (500: default)
		print OUT;
	} <IN>;

	close IN;
	close OUT;
	
	(my $cmpinfo_template = $config_template) =~ s,_tilemap_arg.txt,.cmpinfo,;
	#print Dumper $cmpinfo_template;
	my $cmpinfo = $self->workdir.'/'.$project.'.cmpinfo';

	open(IN, $cmpinfo_template)
		or throw("Can't open config file $cmpinfo_template");

	open(OUT, ">$cmpinfo")
		or throw("Can't open cmpinfo file $cmpinfo");

	my $array_no = ($noa == 1)? '2' : $noa;
	my $group_no = ($noa == 1)? '2': $noa;
	my $groups = ($noa == 1)? '1 2': '1 'x $noa;

	map { 
		s,^(\[Array number\].+=).+$,$1 $array_no,;
		s,^(\[Group number\].+=).+$,$1 $group_no,;
		s,^1 1 1 2 2 2 3 3 3 1 1 1 2 2 2 3 3 3,$groups,;
		s,^\(1>2\) & \(1>3\),1>2,;
		s,^1 2 3,1 2,;
		print OUT;
	} <IN>;

	close IN;
	close OUT;
	
	#system("cp $cmpinfo_template $cmpinfo");

	my $ext;
	if ($method) {
		$ext = '_ma';
	} else {
		$ext = '_hmm';
	}
	# UCSC *.bed file to report significant regions. 
	# Regions are sorted according to their genomic locations.
	my $results = $self->workdir.'/'.$project.$ext.'.bed';
	# This .reg file is a tab-delimited file to report significant regions. 
	# Regions are ranked according to their significance levels.
	#my $results = $self->workdir.'/'.$project.$ext.'.reg';

	$self->resultsfile($results);

	# set columns (fields) for output
	my @fields = (0..2,4); # bed
	#my @fields = (0..2,6); # reg
	$self->output_fields(\@fields);

	return $filename;

}

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
