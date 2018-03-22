# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

=head1 AUTHOR

Stefan Graf, Ensembl Functional Genomics - http://www.ensembl.org

=head1 CONTACT

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

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

	# determine both number of result sets (replicates/arrays) and features

	my $noa = scalar(keys %{$self->result_features});
	warn("\tNo. of result sets: ". $noa);

	my $nof = scalar(@{(values %{$self->result_features})[0]});
	warn("\tNo. of result features: ". $nof);

	# dump features
	my $header = join("\t", 'chromosome', 'position', keys %{$self->result_features});
	$header .= "\tDUMMY" if ($noa == 1);
    
    open(F, ">".$filename)
		or throw("Can't open file $filename.");
    print F $header, "\n";
    
	#print Dumper $self->result_features;
	
	for (my $i=0; $i<$nof; $i++) {
		my $coord=0;
		foreach my $rset (values %{$self->result_features}) {
			print F join("\t", ${$rset}[$i][0], ${$rset}[$i][1]) unless ($coord);
			$coord=1;
            #why do we have to invert the score here???
			print F "\t".(-1*${$rset}[$i][3]);
		}
		print F "\t0" if ($noa == 1);
		print F "\n";
	}
	close F;
    warn("\tresult features written to ".$filename);

	my $workdir = $self->workdir.'/'.$self->analysis->module();
    #warn("\tworkdir: ".$workdir);

	# write config and cmpinfo files from templates
	(my $project = $filename) =~ s,.+/(.+)\.dat,$1,;


	# program parameters

    warn('ANALYSIS PARAMETERS: '.$self->analysis->parameters);

    my %parameters = ();
    map { 
        my ($key, $value) = split (/=/);
        $parameters{$key} = $value;
    } split(/;\s+/, $self->analysis->parameters);

    #print Dumper %parameters;

	my $config = $workdir.'/'.$project.'_args.txt';
	#print Dumper $config;
	$self->config_file($config);

    my $template_file = $parameters{TEMPLATE_FILE};
	open(IN, $template_file)
		or throw("Can't open config file $template_file");

	open(OUT, ">$config")
		or throw("Can't open config file $config");

	map { 
		s,^(O.1-.+=).+$,$1 $workdir,;     #[Working directory]
		s,^(O.2-.+=).+$,$1 $project,;     #[Project Title]
		s,^(I.2-.+=).+$,$1 $project.dat,; #[Raw data file]
		s,^(I.3-.+=).+$,$1 2,;            #[Range of test-statistics] (0: default; 1: [0,1], 2: (-inf, +inf))
		s,^(II.1-.+=).+$,$1 0,;           #[Apply local repeat filter?] (0:No; 1:Yes)
		s,^(II.2-.+=).+$,$1 NULL,;        #[*.refmask file]
		s,^(III.2-.+=).+$,$1 $parameters{METHOD},;    #[Method to combine neighboring probes] (0:HMM, 1:MA)
		s,^(IV.1-.+=).+$,$1 $parameters{POSTPROB},;   #[Posterior probability >]
		s,^(IV.2-.+=).+$,$1 $parameters{MAXGAP},;     #[Maximal gap allowed] (1000: default)
		s,^(IV.4-.+=).+$,$1 0,;           #[Provide your own selection statistics?] (0: No, use default; 1: Yes)
		s,^(IV.5-.+=).+$,$1 NULL,;        #[If Yes to IV.4, selection statistics file]
		s,^(IV.10-.+=).+$,$1 $parameters{HYBLENGTH},; #[Expected hybridization length]
		s,^(V.2-.+=).+$,$1 $parameters{MAXGAP},;      #[Maximal gap allowed] (500: default)
		print OUT;
	} <IN>;

	close IN;
	close OUT;
	
	my $cmpinfo = $workdir.'/'.$project.'.cmpinfo';
    warn("cmp info file: $cmpinfo");

	open(CMP, ">$cmpinfo")
		or throw("Can't open cmpinfo file $cmpinfo");

	my $array_no = ($noa == 1)? '2' : $noa;
	my $group_no = ($noa == 1)? '2' : 1;
	my $groups = ($noa == 1)? '1 2' : '1 'x $noa;
    warn("array_no: " . $array_no);
    warn("group_no: " . $group_no);
    warn("groups: " . $groups);


    print CMP <<EOCMP;
##############################
# TileMap Comparison Info    #
##############################
    
##############################
# Basic Info                 #
##############################
[Array number] = $array_no
[Group number] = $group_no
[Group ID]
$groups

##############################
# Patterns of Interest       #
##############################
[Comparisons]
1>2

##############################
# Preprocessing              #
##############################
[Truncation lower bound] = -1000000000000.0
[Take log2 before calculation?] (1:yes; 0:no) = 0

##############################
# Simulation Setup           #
##############################
[Monte Carlo draws for posterior prob.] = 0

##############################
# Common Variance Groups     #
##############################
[Common variance groups] = 1
$groups

##############################
# Permutation Setup          #
##############################
[Number of permutations] = 0
[Exchangeable groups] = 1
$groups
EOCMP



	close OUT;
	
	#system("cp $cmpinfo_template $cmpinfo");

	my $ext;
	if ($parameters{METHOD}) {
		$ext = '_ma';
	} else {
		$ext = '_hmm';
	}
	# UCSC *.bed file to report significant regions. 
	# Regions are sorted according to their genomic locations.
	my $results = $workdir.'/'.$project.$ext.'.bed';
	# This .reg file is a tab-delimited file to report significant regions. 
	# Regions are ranked according to their significance levels.
	#my $results = $workdir.'/'.$project.$ext.'.reg';

	$self->resultsfile($results);

	# set columns (fields) for output
	my @fields = (1..2,4); # bed
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
    
    system($command) == 0 
        or throw("FAILED to run $command: ".$?);

}

#=head2 infile
#
#  Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::TileMap
#  Arg [2]     : filename (string)
#  Description : will hold a given filename or if one is requested but none
#  defined it will use the create_filename method to create a filename
#  Returntype  : string, filename
#  Exceptions  : none
#  Example     : 
#
#=cut
#
#
#sub infile{
#
#  my ($self, $filename) = @_;
#
#  if($filename){
#    $self->{'infile'} = $filename;
#  }
#  if(!$self->{'infile'}){
#    $self->{'infile'} = $self->create_filename($self->analysis->logic_name, 'dat');
#  }
#
#  return $self->{'infile'};
#
#}


sub config_file {
    my $self = shift;
    $self->{'config_file'} = shift if(@_);
    return $self->{'config_file'};
}

1;
