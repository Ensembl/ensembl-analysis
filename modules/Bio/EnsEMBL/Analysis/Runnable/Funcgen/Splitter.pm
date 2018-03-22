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

=head1 NAME

Bio::EnsEMBL::Analysis::Runnable::Funcgen::Splitter

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Funcgen::Splitter->new
  (
    -analysis => $analysis,
    -query => 'slice',
    -program => 'program.pl',
  );
  $runnable->run;
  my @features = @{$runnable->output};

=head1 DESCRIPTION

Splitter expects to run the command line version of the Splitter program 
(http://zlab.bu.edu/splitter) and predicts features which can be stored in 
the annotated_feature table in the eFG database


=head1 AUTHOR

Stefan Graf, Ensembl Functional Genomics - http://www.ensembl.org

=head1 CONTACT

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

=cut

package Bio::EnsEMBL::Analysis::Runnable::Funcgen::Splitter;

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

    Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::Splitter
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

	# everything needs to happen in a particular subdirectory of the workdir,
	# filenames are set by the Splitter program
	
	(my $workdir = $filename) =~ s/\.dat$//;
	eval { system("mkdir -p $workdir") };
	throw ("Couldn't create directory $workdir: $@") if ($@);
	$self->workdir($workdir);

	$filename =~ s,.+/,,;
	
	# determine both number of result sets (replicates) and features

	my $noa = scalar(keys %{$self->result_features});
	warn('no. of result sets: ', $noa);

	my $nof = scalar(@{(values %{$self->result_features})[0]});
	warn('no. of result fts: ', $nof);

	# dump features
	open(F, "> $workdir/$filename")
		or throw("Can't open file $workdir/$filename.");

	#print Dumper $self->result_features;
	
	for (my $i=0; $i<$nof; $i++) {
		my $coord=0;
		foreach my $rset (values %{$self->result_features}) {
			printf F "chr%s\t%d\t%d\t", ${$rset}[$i][0], ${$rset}[$i][1], ${$rset}[$i][2] unless ($coord);
			$coord=1;    
			print F "\t".${$rset}[$i][3];
		}
		print F "\n";
	}
	close F;

	$self->infile($workdir.'/'.$filename);
	return "$workdir/$filename";
}

=head2 run_analysis

  Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::Splitter
  Arg [2]     : string, program name
  Usage       : 
  Description : 
  Returns     : 
  Exceptions  : 

=cut

### step 1
# ./splitter.cgi "id=SplitterJob.$(gdate -I)" 'submit=Step by step' 'signalfile=samplesplittersignal' 'sorted=yes' 'norm=qnorm' 'combine=median' 'splitterstep=10' 'splitterratio=2' 'cutoff=2.5sd' 'maxgap=100' 'minrun=5'
### step 2
# ./splitter.cgi "id=SplitterJob.$(gdate -I)" 'submit=Step 2' 'sorted=yes' 'norm=none' 'combine=median'
### step 3
# ./splitter.cgi "id=SplitterJob.$(gdate -I)" 'submit=Step 3' 'last=pipeline.median' 'splitterfrom=0.0149' 'splitterto=0.4884' 'splitterstep=10' 'splitterratio=2' 'cutoff=2.5sd' 'custom=0.4884' '2.5sd=0.4884' '1pct=0.395' '5pct=0.25' '2sd=0.3937' 'maxgap=100' 'minrun=5'

### One shot (all in one)
# ./splitter.cgi "id=SplitterJob.TEST" 'submit=One shot' 'signalfile=samplesplittersignal' 'sorted=yes' 'norm=none' 'combine=median' 'splitterstep=10' 'splitterratio=2' 'cutoff=2.5sd' 'maxgap=100' 'minrun=5'

#available Splitter parameters (for details see http://zlab.bu.edu/splitter)
# -Input already sorted by genomic positions of probes
#  'sorted=yes|no'
# -Perform normalization between replicates
#  'norm=none|qnorm|rnorm'
# -Combine replicates
#  'combine=mean|median'
# -Signal cutoff for the combined intensities>=
#  'cutoff=splitter|5pct|1pct|2sd|2.5sd'
#   'splitterfrom='
#   'splitterto='
#   'splitterstep=10'
#   'splitterratio=2'
# -Determine how the probes are clustered
#  Gap (Maxgap) <= n base pairs
#   'maxgap=100'
#  Clustering (Minrun) >= n probes
#   'minrun=4'

sub run_analysis {
        
    my ($self, $program) = @_;
    
    if(!$program){
        $program = $self->program;
    }
    throw($program." is not executable Splitter::run_analysis ") 
        unless($program && -x $program);

	(my $id = $self->workdir) =~ s,^.+/,,;

    my $command = $self->program." \'id=$id\' \'submit=One shot\' \'signalfile=".$self->infile."\' " .
		$self->analysis->parameters;
    
    warn("Running analysis " . $command . "\n");
    
    eval { system( $command ) };
    throw("FAILED to run $command: ", $@) if ($@);

	opendir(DIR, $self->workdir())
		or throw("Can't open workdir ".$self->workdir());
	my @resultfiles = grep { /\.scored\.txt$/ } readdir DIR;
	closedir(DIR);

	throw("More than one resultfile in ".$self->workdir()) 
		if (scalar(@resultfiles) > 1);

	# set resultsfile
	$self->resultsfile($self->workdir.'/'.$resultfiles[0]);

}

=head2 infile

  Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::Splitter
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

	  next unless (s/^chr//);
	  
      chomp;
      my @ft = split;
	  my $coord = shift(@ft);
	  my @coord = split(/[:-]/, $coord);
      push(@output, [ @coord, $ft[0] ]);

  }

  throw("No features to annotate  on slice ".$self->query->seq_region_name."!") if (scalar(@output) == 0);

  $self->output(\@output);
  #print Dumper $self->output();

  throw("parse_results: can't close file ".$resultsfile.".")
      unless (close(F));

}


1;
