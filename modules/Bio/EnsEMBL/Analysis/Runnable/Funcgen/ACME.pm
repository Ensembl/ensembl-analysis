# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

Bio::EnsEMBL::Analysis::Runnable::Funcgen::ACME

=head1 SYNOPSIS

  my $runnable = Bio::EnsEMBL::Analysis::Runnable::Funcgen::ACME->new
  (
    -analysis => $analysis,
    -query => 'slice',
    -program => 'program.pl',
  );
  $runnable->run;
  my @features = @{$runnable->output};

=head1 DESCRIPTION

ACME expects to run the R program ACME (Scacherie et al. (2006), PMID: 
16939795) and predicts features which can be stored in the annotated_feature 
table in the eFG database

=head1 LICENCE

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2022] EMBL-European Bioinformatics Institute
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

=head1 AUTHOR

Stefan Graf, Ensembl Functional Genomics - http://www.ensembl.org

=head1 CONTACT

Post questions to the Ensembl development list: http://lists.ensembl.org/mailman/listinfo/dev

=cut

package Bio::EnsEMBL::Analysis::Runnable::Funcgen::ACME;

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

    Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::ACME
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

    warn("filename: ", $filename);
    

	# determine both number of result sets (replicates) and features

	my $noa = scalar(keys %{$self->result_features});
	warn ("Processing $noa replicates. Resulting annotated features will be \n".
		  "linked all together to the same features set". $noa) if ($noa > 1) ;

	# dump features in GFF format
	
	my $workdir = $self->workdir.'/'.$self->analysis->module();
    warn($workdir);

	#print Dumper $self->result_features;

	my @fnames=();
	
	foreach my $rset_name (keys %{$self->result_features}) {
		(my $fname = $filename) =~ s,^.+/(.+)\.dat,$1_${rset_name}.dat,;
        warn("fname: ", $fname);
		push (@fnames, $fname);
		open(F, ">$workdir/$fname")
			or throw("Can't open file $filename.");
		foreach my $ft (@{${$self->result_features}{$rset_name}}) {
			printf F "%s\t%s\t%s\t%d\t%d\t%f\t.\t.\t.\n", 
			$ft->[0], $ENV{NORM_METHOD}, $rset_name, $ft->[1], $ft->[2], $ft->[3];
		}
		close F;
	}
    warn("fnames: ", @fnames);

	(my $Rfile = $filename) =~ s/\.dat/.R/;
	$self->infile($Rfile);

	# set resultsfile
	(my $resultsfile = $filename) =~ s/\.dat/.regions/;
	$self->resultsfile($resultsfile);
	
	(my $PDFfile = $filename) =~ s/\.dat/.pdf/;

	my $fnames = join('", "', @fnames);
    
	# parse parameters
	my $param = $self->get_parameters();

	my $window = $param->{WINDOW};
	my $threshold = $param->{THRESHOLD};

	my $pvalue = $param->{PVALUE};

	#my $chrom = $self->query()->seq_region_name;
	#my $start = $self->query()->start, 
	#my $end = $self->query()->end;
	
	open(R, ">".$Rfile)
		or throw("Can't open file $Rfile.");
	
	print R <<EOR;

	library("ACME")
	setwd("$workdir")
	fnames <- c("$fnames")
	a <- read.resultsGFF(fnames)
	colnames(a\@data) <- fnames
	#str(a)
	aGFFCalc <- do.aGFF.calc(a, window = $window, thresh = $threshold)
	regions <- findRegions(aGFFCalc,thresh=$pvalue)
	write.table(regions[regions\$TF,c('Chromosome', 'Start', 'End', 'Mean', 'Median')],
				file="$resultsfile", sep="\\t", 
				row.names = FALSE, col.names = FALSE)

EOR

    #f=file('${resultsfile}.bed','w')
    #writeLines('track name="ACME findRegions $pvalue" description="$pvalue"',con=f)
    #regions\$Chromosome <- paste('chr',regions\$Chromosome,sep="")
    #write.table(regions[regions\$TF==TRUE,c('Chromosome','Start','End')],f,sep="\t",col.names=F,row.names=F,quote=F)

	#pdf("$PDFfile", height = 10, width = 15);
	#plot(aGFFCalc, chrom = "$chrom", sample = 1)
    #dev.off()


	close R;

	# set columns (fields) for output
	my @fields = (1..3); # mean
	#my @fields = (0..2,4); # median
	$self->output_fields(\@fields);

}

=head2 run_analysis

  Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::ACME
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
    throw($program." is not executable ACME::run_analysis ") 
        unless($program && -x $program);

    my $command = $self->program . ' --no-save --slave < ' . $self->infile;
    
    warn("Running analysis " . $command . "\n");
    
    eval { system($command) };
    throw("FAILED to run $command: ", $@) if ($@);

}

#=head2 infile
#
#  Arg [1]     : Bio::EnsEMBL::Analysis::Runnable::ACME
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
